"""
MSseg CLI — 全切片 MCseg v2 分割指令列介面
==============================================

Usage:
    # 全切片 7-pass 分割（含 cpsam）
    uv run python -m backend.src.cli.segment \\
        --btf "K:/path/to/image.btf" \\
        --tp  "K:/path/to/tissue_positions.parquet" \\
        --h5  "K:/path/to/filtered_feature_bc_matrix.h5" \\
        --out "K:/path/to/output/" \\
        --tissue crc \\
        --cpsam

    # 快速 4-pass 分割（無 cpsam）
    uv run python -m backend.src.cli.segment \\
        --btf "K:/path/to/image.btf" \\
        --out "K:/path/to/output/"

    # 從已有 he_crop.tif 跳過裁切，直接分割
    uv run python -m backend.src.cli.segment \\
        --he-crop "K:/path/to/he_crop.tif" \\
        --out "K:/path/to/output/" \\
        --cpsam
"""

from __future__ import annotations

import argparse
import gc
import logging
import sys
import time
from pathlib import Path

import numpy as np

# ─── Logging ─────────────────────────────────────────────────────────────────
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger("msseg.cli")

# ─── Tissue presets ───────────────────────────────────────────────────────────
# 參數一律取自 `config/profiles/{tissue}.yaml`（CLAUDE.md §15：禁止硬編碼）。
#
# 這裡曾有一份硬編碼的 `TISSUE_PRESETS`，與 profile 各走一套而**靜默漂移**：
# 實測 CRC 的 `voronoi_distance` 差 1（8 vs 9）、LUAD 有 4 個參數不一致
# （voronoi 9/8、clahe 2.5/2.02、flow 0.4/0.415、max_size 5000/4000），
# 且 voronoi 在兩個組織是**反方向**差異 —— 看得出是某次複製貼上寫錯而非刻意調參。
# 後果是 CLI 與 Web UI 對「同一個 crc preset」會產出不同遮罩，而 README 卻聲稱
# 兩者「use exactly the same engine with consistent parameter semantics」。
#
# profile 的值才是有依據的那一組：LUAD 的 clahe 2.02 / flow 0.415 都標註為
# 「xenium_he_seg 最佳參數」，CLI 那組整數是最佳化**之前**的預設值。


def available_tissues() -> list[str]:
    """
    列出 `config/profiles/` 下可用的 tissue profile 名稱。

    必須濾掉 macOS 在 ExFAT/外接磁碟產生的 `._*` AppleDouble 檔（CLAUDE.md §4）
    —— 否則 `--tissue` 的選項會多出一個 `._default` 之類的幽靈組織。
    """
    from backend.src.utils.config import _PROFILES_DIR

    return sorted(
        p.stem for p in _PROFILES_DIR.glob("*.yaml")
        if not p.name.startswith("._")
    )


def load_tissue_preset(tissue: str) -> dict:
    """
    從 `config/profiles/{tissue}.yaml` 取 `segmentation.mcseg_v2` 參數。

    profile 缺少該區塊時視為設定錯誤而中止 —— 靜默跑一組空參數（等同 cellpose
    預設值）比直接失敗更糟，使用者會拿到一份看似正常但參數完全不對的遮罩。
    """
    from backend.src.utils.config import _load_profile

    profile = _load_profile(tissue)
    preset = (profile.get("segmentation") or {}).get("mcseg_v2") or {}
    if not preset:
        raise SystemExit(
            f"config/profiles/{tissue}.yaml 缺少 segmentation.mcseg_v2 區塊，"
            f"無法取得分割參數。可用 profile：{', '.join(available_tissues())}"
        )
    return dict(preset)


def _load_cli_config() -> dict:
    """
    取設定檔中與 CLI 相關的非路徑設定（目前為 `alignment:`）。

    CLI 的樣本一律由引數指定，故**不採用**設定檔的 `paths:` —— 設定檔與
    `state.json` 指向的可能是另一個樣本。讀不到設定檔時回空 dict，
    `resolve_bin_to_image_transform` 的預設值（自動偵測對位 JSON）仍然適用。
    """
    try:
        from backend.src.utils.config import load_config

        return {"alignment": load_config().get("alignment") or {}}
    except (OSError, ValueError, KeyError) as e:
        log.warning(f"讀取設定檔失敗（{e}），alignment 設定改用預設值")
        return {}


# ─── Step helpers ─────────────────────────────────────────────────────────────

def step_crop_btf(
    btf_path: Path,
    out_dir: Path,
    crop_y0: int,
    crop_y1: int,
    btf_col0: int,
    btf_col1: int,
) -> np.ndarray:
    """從 BTF 裁切 H&E 影像，若已存在則直接載入。"""
    import tifffile

    from backend.src.utils.slide_reader import open_slide

    crop_tif = out_dir / "he_crop.tif"
    if crop_tif.exists():
        log.info(f"[SKIP] 載入已存在的 H&E crop: {crop_tif.name}")
        img = tifffile.imread(str(crop_tif))
        if img.ndim == 3 and img.shape[-1] == 4:
            img = img[..., :3]
        log.info(f"  shape: {img.shape}")
        return img

    log.info(f"[1/4] 從 BTF 裁切 H&E (row {crop_y0}:{crop_y1}, col {btf_col0}:{btf_col1})")
    t0 = time.time()
    # 與 API `/run_full` 共用 SlideReader：只解壓被請求的 tile。
    # （原本走 `tifffile.aszarr()`，但它要求 zarr >= 3，而 anndata/scanpy 這條鏈
    #   還在 zarr 2.18.7 —— 在本環境會直接拋 ValueError，CLI 的 --btf 路徑全斷。）
    img = open_slide(btf_path).read_region(
        btf_col0, crop_y0, btf_col1 - btf_col0, crop_y1 - crop_y0
    )
    if img.ndim == 3 and img.shape[-1] == 4:
        img = img[..., :3]
    log.info(f"  shape: {img.shape}  ({time.time() - t0:.0f}s)")
    tifffile.imwrite(str(crop_tif), img, compression="zlib")
    log.info(f"  儲存: {crop_tif.name}")
    return img


def step_load_he_crop(he_crop_path: Path) -> np.ndarray:
    """直接從 he_crop.tif 載入影像。"""
    import tifffile
    log.info(f"[1/4] 載入 H&E crop: {he_crop_path.name}")
    img = tifffile.imread(str(he_crop_path))
    if img.ndim == 3 and img.shape[-1] == 4:
        img = img[..., :3]
    log.info(f"  shape: {img.shape}")
    return img


def step_segment(
    img: np.ndarray,
    cfg: dict,
    out_dir: Path,
    tile_size: int = 1024,
    overlap: int = 128,
    max_load_gb: float = 16.0,
) -> np.ndarray:
    """MCseg v2 分割，輸出 mcseg_mask.npy。"""
    mask_path = out_dir / "mcseg_mask.npy"
    if mask_path.exists():
        log.info(f"[SKIP] 載入已存在的遮罩: {mask_path.name}")
        mask = np.load(str(mask_path))
        log.info(f"  shape: {mask.shape}  cells: {int(mask.max()):,}")
        return mask

    from backend.src.fullslide.pipeline import (
        apply_full_seg_safety_clamp,
        check_full_seg_mask_memory,
    )
    from backend.src.segmentation.cellpose_runner import run_tiled_mcseg_v2

    h, w = img.shape[:2]
    check_full_seg_mask_memory(w, h, max_load_gb)

    # MPS batch_size≤2 安全鉗制與 API 全片路徑一致（docs/adr/0005：無條件套用）。
    # force_disable_cpsam=False：CLI 裁切窗格由使用者自行控制大小，不像 Web UI
    # 全片按鈕預設面對整張未知大小的切片，所以尊重既有的 --cpsam 旗標行為。
    safe_cfg = apply_full_seg_safety_clamp(cfg, force_disable_cpsam=False)
    if safe_cfg["batch_size"] != cfg.get("batch_size"):
        log.info(f"  batch_size {cfg.get('batch_size')} → {safe_cfg['batch_size']}"
                  "（MPS 安全鉗制，見 docs/adr/0005）")

    passes = 7 if safe_cfg.get("use_cpsam") else 4
    log.info(f"[2/4] MCseg v2 {passes}-pass 分割（GPU={safe_cfg.get('use_gpu', True)}）")

    def _progress(p: float, msg: str) -> None:
        bar = "█" * int(p * 30) + "░" * (30 - int(p * 30))
        log.info(f"  [{bar}] {p*100:.0f}%  {msg}")

    t0 = time.time()
    mask = run_tiled_mcseg_v2(
        img,
        safe_cfg,
        tile_size=tile_size,
        overlap=overlap,
        progress_callback=_progress,
    )
    elapsed = time.time() - t0
    log.info(f"  完成！{int(mask.max()):,} 個細胞  耗時 {elapsed/60:.1f} min")

    np.save(str(mask_path), mask)
    log.info(f"  儲存: {mask_path.name}")
    return mask


def step_bin_attribution(
    mask: np.ndarray,
    tp_path: Path,
    out_dir: Path,
    crop_y0: int,
    btf_col0: int,
    scale: tuple[float, float] = (1.0, 1.0),
    transform: "np.ndarray | None" = None,
    transform_source: str = "",
) -> "pd.DataFrame":  # noqa: F821
    """將 Visium HD 2µm bins 對齊到細胞遮罩（快取 + log 包裝）。"""
    import pandas as pd

    from backend.src.fullslide.pipeline import bin_attribution

    attr_path = out_dir / "bin_attribution.parquet"
    if attr_path.exists():
        log.info(f"[SKIP] 載入已存在的 bin attribution: {attr_path.name}")
        return pd.read_parquet(str(attr_path))

    log.info("[3/4] Bin attribution")
    log.info(f"  bins → 影像變換來源: {transform_source or '未指定'}")
    if transform is None:
        log.info(f"  近似縮放: x={scale[0]:.4f}, y={scale[1]:.4f}")
    attr = bin_attribution(
        mask, tp_path, crop_y0, btf_col0, out_path=attr_path,
        scale=scale, transform=transform,
    )
    log.info(f"  attributed bins: {len(attr):,}")
    log.info(f"  儲存: {attr_path.name}")
    return attr


def step_count_cells(
    mask: np.ndarray,
    attribution: "pd.DataFrame",  # noqa: F821
    h5_path: Path,
    out_dir: Path,
    pixel_size_um: float,
    origin_xy: tuple[int, int] = (0, 0),
) -> Path:
    """[4/6] 由 bin attribution 聚合 cells×genes，附加重心，輸出 cells.h5ad。"""
    from backend.src.fullslide.pipeline import add_centroids, aggregate_cells

    h5ad_path = out_dir / "cells.h5ad"
    if h5ad_path.exists():
        log.info(f"[SKIP] 載入已存在的 cells h5ad: {h5ad_path.name}")
        return h5ad_path

    log.info("[4/6] 聚合 cells×genes 矩陣")
    log.info(f"  讀取 h5 矩陣: {h5_path.name}")
    cells = aggregate_cells(attribution, h5_path)
    log.info(f"  unique cells with RNA: {cells.n_obs:,}")

    log.info("  計算細胞重心…")
    add_centroids(cells, mask, pixel_size_um, origin_xy=origin_xy)

    cells.write_h5ad(str(h5ad_path))
    log.info(f"  儲存: {h5ad_path.name}  ({cells.n_obs:,} cells × {cells.n_vars:,} genes)")
    return h5ad_path


def step_celltypist(
    cells_h5ad_path: Path,
    out_dir: Path,
    celltypist_model: str,
) -> "pd.DataFrame":  # noqa: F821
    """[5/6] CellTypist 標注；寫出 celltypist_labels.csv 並回寫標籤至 cells.h5ad。"""
    import pandas as pd
    import scanpy as sc

    csv_path = out_dir / "celltypist_labels.csv"
    if csv_path.exists():
        log.info(f"[SKIP] 載入已存在的 CellTypist 結果: {csv_path.name}")
        return pd.read_csv(str(csv_path))

    log.info(f"[5/6] CellTypist 標注（model={celltypist_model}）")
    cells = sc.read_h5ad(str(cells_h5ad_path))

    adata_norm = cells.copy()
    sc.pp.normalize_total(adata_norm, target_sum=1e4)
    sc.pp.log1p(adata_norm)

    import celltypist
    predictions = celltypist.annotate(
        adata_norm, model=celltypist_model, majority_voting=False,
    )
    ct_labels = predictions.predicted_labels["predicted_labels"].values

    df = pd.DataFrame(
        {"cell_id": cells.obs["cell_id"].values, "celltypist_label": ct_labels}
    )
    df.to_csv(str(csv_path), index=False)

    # 回寫標籤至 cells.h5ad，使其成為含註解的標準產物
    cells.obs["celltypist_label"] = ct_labels
    cells.write_h5ad(str(cells_h5ad_path))

    log.info(f"  儲存: {csv_path.name}（標籤亦回寫 {cells_h5ad_path.name}）")
    log.info(f"  細胞型態分佈:\n{df['celltypist_label'].value_counts().head(10).to_string()}")
    return df


def step_export_xenium(
    mask: np.ndarray,
    cells_h5ad_path: Path,
    out_dir: Path,
    pixel_size_um: float,
    he_image_path: Path | None = None,
) -> Path:
    """[6/6] 將整片遮罩 + cells.h5ad 匯出為 Xenium Explorer bundle。"""
    import json

    xen_dir = out_dir / "xenium_explorer"
    if (xen_dir / "experiment.xenium").exists():
        log.info(f"[SKIP] Xenium bundle 已存在: {xen_dir.name}")
        return xen_dir

    # 1. 細胞多邊形 GeoJSON（局部 µm，原點 = 裁切左上角；格式同 GUI 匯出）
    # 與 Web UI 共用同一套多邊形產生邏輯與 min_area_px 雜訊過濾
    # （backend/src/export/geometry.py），mask 已在記憶體中故直接吃陣列版本，
    # 不重新從硬碟 np.load 一次。
    geojson_path = out_dir / "cells_polygons.geojson"
    if geojson_path.exists():
        log.info(f"[SKIP] 載入已存在的多邊形: {geojson_path.name}")
    else:
        n_cells = int(mask.max())
        log.info(f"[6/6] 產生細胞多邊形 GeoJSON（{n_cells:,} cells，整片可能較久）…")
        from backend.src.export.geometry import _mask_array_to_geojson

        feature_collection = _mask_array_to_geojson(mask, pixel_size_um, min_area_px=20)
        with open(geojson_path, "w", encoding="utf-8") as f:
            json.dump(feature_collection, f)
        log.info(f"  多邊形數: {len(feature_collection['features']):,} → {geojson_path.name}")

    # 2. 組裝 Xenium Explorer bundle（多邊形 µm 座標與 cells.h5ad obs['cell_id'] 對齊）
    log.info("  匯出 Xenium Explorer bundle…")
    from backend.src.export.xenium_exporter import XeniumExporter

    exporter = XeniumExporter(
        poly_json_path=geojson_path,
        pixel_size_um=pixel_size_um,
        he_image_path=he_image_path if (he_image_path and he_image_path.exists()) else None,
    )
    exporter.export(cells_h5ad_path, xen_dir)
    log.info(f"  ✅ Xenium bundle: {xen_dir}")
    return xen_dir


# ─── CLI entry point ──────────────────────────────────────────────────────────

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="python -m backend.src.cli.segment",
        description=(
            "MSseg CLI — Visium HD BTF 全切片流程（分割 → 計數 → 細胞型注釈 → 選配 Xenium 匯出）\n"
            "輸出: he_crop.tif / mcseg_mask.npy / bin_attribution.parquet / "
            "cells.h5ad / celltypist_labels.csv / xenium_explorer/"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    # ── 輸入來源（二選一）
    src = p.add_argument_group("輸入影像（二選一）")
    src.add_argument("--btf",     type=Path, metavar="PATH", help="原始 BigTIFF (.btf) 路徑")
    src.add_argument("--he-crop", type=Path, metavar="PATH", help="已裁切的 he_crop.tif（略過 BTF 裁切步驟）")

    # ── BTF 裁切座標（--btf 時必填）
    crop = p.add_argument_group("BTF 裁切座標（--btf 時使用）")
    crop.add_argument("--crop-y0",    type=int, default=0,     metavar="PX", help="裁切起始 row（BTF 全圖座標，預設 0）")
    crop.add_argument("--crop-y1",    type=int, default=-1,    metavar="PX", help="裁切結束 row（-1 = 全圖）")
    crop.add_argument("--btf-col0",   type=int, default=0,     metavar="PX", help="裁切起始 col（BTF 全圖座標，預設 0）")
    crop.add_argument("--btf-col1",   type=int, default=-1,    metavar="PX", help="裁切結束 col（-1 = 全圖）")

    # ── RNA 計數（選填）
    rna = p.add_argument_group("RNA 計數（選填，需同時提供 --tp 與 --h5）")
    rna.add_argument("--tp",  type=Path, metavar="PATH", help="tissue_positions.parquet 路徑")
    rna.add_argument("--h5",  type=Path, metavar="PATH", help="filtered_feature_bc_matrix.h5 路徑")

    # ── 輸出
    p.add_argument("--out", type=Path, required=True, metavar="DIR",
                   help="輸出目錄（自動建立）")

    # ── 分割參數
    seg = p.add_argument_group("分割參數")
    _tissues = available_tissues()
    seg.add_argument("--tissue",  choices=_tissues, default="crc",
                     help=f"組織類型 profile（{' / '.join(_tissues)}），預設 crc。"
                          "參數取自 config/profiles/{tissue}.yaml")
    seg.add_argument("--cpsam",   action="store_true",
                     help="啟用 cpsam（7-pass，需更長時間）")
    seg.add_argument("--no-gpu",  action="store_true",
                     help="強制使用 CPU（預設自動偵測 GPU）")
    seg.add_argument("--batch-size", type=int, default=2, metavar="N",
                     help="Cellpose batch size（預設 2，VRAM 不足時調低）")
    seg.add_argument("--tile-size",  type=int, default=1024, metavar="PX",
                     help="Tile 大小（預設 1024）")
    seg.add_argument("--overlap",    type=int, default=128,  metavar="PX",
                     help="Tile 重疊寬度（預設 128）")
    seg.add_argument("--max-load-gb", type=float, default=16.0, metavar="GB",
                     help="輸出遮罩（int32）記憶體上限，超過則報錯而非默默 OOM（預設 16.0）")
    seg.add_argument("--dia-small",  type=float, metavar="PX",
                     help="小直徑 pass（覆寫 tissue preset）")
    seg.add_argument("--dia-mid",    type=float, metavar="PX",
                     help="主直徑 pass（覆寫 tissue preset）")
    seg.add_argument("--dia-large",  type=float, metavar="PX",
                     help="大直徑 pass（覆寫 tissue preset）")
    seg.add_argument("--voronoi-d",  type=int, metavar="PX",
                     help="Voronoi 擴張距離（覆寫 tissue preset）")
    seg.add_argument("--cellprob",   type=float, metavar="THRESH",
                     help="cellprob_threshold（覆寫 tissue preset）")

    # ── CellTypist
    ct = p.add_argument_group("CellTypist（需 --tp 與 --h5）")
    ct.add_argument("--celltypist-model", default="Human_Colorectal_Cancer.pkl",
                    metavar="MODEL",
                    help="CellTypist 模型名稱（預設 Human_Colorectal_Cancer.pkl）")
    ct.add_argument("--skip-celltypist", action="store_true",
                    help="跳過 CellTypist 標注")

    # ── Browser 匯出（選填）
    exp = p.add_argument_group("Browser 匯出（選填，需 --tp 與 --h5）")
    exp.add_argument("--export-xenium", action="store_true",
                     help="匯出 Xenium Explorer bundle（整片細胞數多時較耗時）")

    return p


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    # ── 驗證輸入
    if args.btf is None and args.he_crop is None:
        parser.error("請提供 --btf 或 --he-crop 其中之一")
    if args.btf and not args.btf.exists():
        parser.error(f"BTF 檔案不存在: {args.btf}")
    if args.he_crop and not args.he_crop.exists():
        parser.error(f"he_crop.tif 不存在: {args.he_crop}")

    out_dir: Path = args.out
    out_dir.mkdir(parents=True, exist_ok=True)
    log.info(f"MSseg CLI — 輸出目錄: {out_dir}")

    # ── 建立分割設定（基底來自 config/profiles/{tissue}.yaml）
    cfg = load_tissue_preset(args.tissue)
    cfg["use_gpu"]              = not args.no_gpu
    cfg["batch_size"]           = args.batch_size
    cfg["use_cpsam"]            = args.cpsam
    # use_hematoxylin 由 profile 決定（缺少時才預設啟用）—— 原本在此硬編碼 True，
    # 會讓刻意關閉它的 profile 被安靜忽略
    cfg.setdefault("use_hematoxylin", True)
    cfg["use_transcript_rescue"] = False
    if args.dia_small  is not None: cfg["dia_small"]           = args.dia_small
    if args.dia_mid    is not None: cfg["dia_mid"]             = args.dia_mid
    if args.dia_large  is not None: cfg["dia_large"]           = args.dia_large
    if args.voronoi_d  is not None: cfg["voronoi_distance"]    = args.voronoi_d
    if args.cellprob   is not None: cfg["cellprob_threshold"]  = args.cellprob

    passes = 7 if cfg["use_cpsam"] else 4
    log.info(
        f"設定摘要: tissue={args.tissue}  passes={passes}  "
        f"gpu={cfg['use_gpu']}  batch={cfg['batch_size']}  "
        f"dia={cfg['dia_small']}/{cfg['dia_mid']}/{cfg['dia_large']}px  "
        f"voronoi_d={cfg['voronoi_distance']}px"
    )

    # ── Step 1: 取得 H&E 影像
    # `image_shape` 是**來源影像**（未裁切）的 (h, w)，供座標變換推導使用；
    # 拿裁切後的遮罩尺寸當畫布會算出荒謬的縮放（見 fullslide.pipeline 註解）。
    if args.he_crop:
        img = step_load_he_crop(args.he_crop)
        crop_y0  = 0
        btf_col0 = 0
        # he_crop.tif 若本身是更大張影像的一塊，CLI 無從得知其原點與原尺寸
        image_shape = img.shape[:2]
    else:
        import tifffile

        from backend.src.fullslide.pipeline import _upper_is_open

        with tifffile.TiffFile(str(args.btf)) as tif:
            full_shape = tif.pages[0].shape
        image_shape = (int(full_shape[0]), int(full_shape[1]))
        # -1 代表「取到影像邊界」，與 API 全片路徑（FullSegParams）同一套哨兵語意
        crop_y1  = image_shape[0] if _upper_is_open(args.crop_y1)  else args.crop_y1
        btf_col1 = image_shape[1] if _upper_is_open(args.btf_col1) else args.btf_col1
        crop_y0  = args.crop_y0
        btf_col0 = args.btf_col0
        img = step_crop_btf(args.btf, out_dir, crop_y0, crop_y1, btf_col0, btf_col1)

    # ── Step 2: 分割
    mask = step_segment(img, cfg, out_dir, args.tile_size, args.overlap, args.max_load_gb)
    del img
    gc.collect()

    # ── Step 3: Bin attribution（有 tp & h5 才跑）
    # µm/px 取樣本實際值（scalefactors_json.json）優先，缺失才用預設常數。
    # --tp 指向 {binned_002}/spatial/tissue_positions.parquet，故其祖父層即 binned_002。
    from backend.src.fullslide.pipeline import (
        resolve_bin_to_image_transform,
        resolve_pixel_size,
    )
    _binned_002 = str(args.tp.parent.parent) if args.tp else ""
    # `alignment:` 等設定取自設定檔，但 paths 一律以 CLI 引數為準（CLI 可跑任意樣本，
    # 設定檔／state.json 指向的可能是完全不同的樣本）
    _cfg = _load_cli_config()
    _cfg.setdefault("paths", {})
    _cfg["paths"]["binned_002"] = _binned_002
    if args.btf:
        _cfg["paths"]["he_image"] = str(args.btf)
    pixel_size_um = resolve_pixel_size(_cfg)
    # SR fullres 與分割影像未必同座標系：優先用 Loupe/CytAssist 對位 JSON 組出的
    # homography（幾何正確），找不到才回退近似分軸縮放。以**來源影像**尺寸推導。
    bin_transform, bin_scale, transform_source = resolve_bin_to_image_transform(
        _cfg, image_shape
    )
    # transform 已含縮放；同時傳 scale 會讓 bin_attribution 每次正常執行都警告
    # 「同時指定 transform 與 scale」。與 resolve_full_count_inputs 的作法一致。
    if bin_transform is not None:
        bin_scale = (1.0, 1.0)
    log.info(f"  pixel_size_um = {pixel_size_um}")

    attribution = None
    if args.tp and args.h5:
        if not args.tp.exists():
            log.warning(f"tissue_positions 不存在，跳過 RNA 計數: {args.tp}")
        elif not args.h5.exists():
            log.warning(f"h5 矩陣不存在，跳過 RNA 計數: {args.h5}")
        else:
            attribution = step_bin_attribution(
                mask, args.tp, out_dir, crop_y0, btf_col0,
                scale=bin_scale, transform=bin_transform,
                transform_source=transform_source,
            )

    # ── Step 4: 聚合 cells×genes h5ad（有 attribution 才跑）
    cells_h5ad = None
    if attribution is not None:
        cells_h5ad = step_count_cells(
            mask, attribution, args.h5, out_dir, pixel_size_um,
            origin_xy=(btf_col0, crop_y0),
        )

        # ── Step 5: CellTypist
        if not args.skip_celltypist:
            step_celltypist(cells_h5ad, out_dir, args.celltypist_model)

    # ── Step 6: Xenium Explorer 匯出（選填）
    if args.export_xenium:
        if cells_h5ad is None:
            log.warning("--export-xenium 需要 --tp 與 --h5（產生 cells.h5ad）才能匯出，已跳過")
        else:
            he_for_img = args.he_crop if args.he_crop else (out_dir / "he_crop.tif")
            step_export_xenium(mask, cells_h5ad, out_dir, pixel_size_um, he_for_img)

    log.info("=" * 60)
    log.info(f"✅ MSseg CLI 完成！  結果目錄: {out_dir}")
    log.info(f"   細胞數量: {int(mask.max()):,}")
    log.info("=" * 60)
    return 0


if __name__ == "__main__":
    sys.exit(main())
