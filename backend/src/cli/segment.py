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
TISSUE_PRESETS: dict[str, dict] = {
    "crc": {
        "dia_small": 13.0,
        "dia_mid": 17.0,
        "dia_large": 22.0,
        "voronoi_distance": 8,
        "clahe_clip_limit": 3.0,
        "flow_threshold": 0.4,
        "cellprob_threshold": -2.0,
        "min_size": 20,
        "max_size": 6000,
    },
    "luad": {
        "dia_small": 10.0,
        "dia_mid": 14.0,
        "dia_large": 18.0,
        "voronoi_distance": 9,
        "clahe_clip_limit": 2.5,
        "flow_threshold": 0.4,
        "cellprob_threshold": -1.5,
        "min_size": 20,
        "max_size": 5000,
    },
    "default": {
        "dia_small": 13.0,
        "dia_mid": 17.0,
        "dia_large": 22.0,
        "voronoi_distance": 9,
        "clahe_clip_limit": 3.0,
        "flow_threshold": 0.4,
        "cellprob_threshold": -2.0,
        "min_size": 20,
        "max_size": 6000,
    },
}


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
    import zarr

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
    with tifffile.TiffFile(str(btf_path)) as tif:
        store = tif.aszarr()
    z = zarr.open(store, mode="r")
    arr = z[0] if z.ndim == 4 else z
    img = np.asarray(arr[crop_y0:crop_y1, btf_col0:btf_col1])
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


def step_segment(img: np.ndarray, cfg: dict, out_dir: Path) -> np.ndarray:
    """MCseg v2 分割，輸出 mcseg_mask.npy。"""
    mask_path = out_dir / "mcseg_mask.npy"
    if mask_path.exists():
        log.info(f"[SKIP] 載入已存在的遮罩: {mask_path.name}")
        mask = np.load(str(mask_path))
        log.info(f"  shape: {mask.shape}  cells: {int(mask.max()):,}")
        return mask

    from backend.src.segmentation.cellpose_runner import run_tiled_mcseg_v2

    passes = 7 if cfg.get("use_cpsam") else 4
    log.info(f"[2/4] MCseg v2 {passes}-pass 分割（GPU={cfg.get('use_gpu', True)}）")

    def _progress(p: float, msg: str) -> None:
        bar = "█" * int(p * 30) + "░" * (30 - int(p * 30))
        log.info(f"  [{bar}] {p*100:.0f}%  {msg}")

    t0 = time.time()
    mask = run_tiled_mcseg_v2(
        img,
        cfg,
        tile_size=1024,
        overlap=128,
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
) -> "pd.DataFrame":  # noqa: F821
    """將 Visium HD 2µm bins 對齊到細胞遮罩。"""
    import pandas as pd

    attr_path = out_dir / "bin_attribution.parquet"
    if attr_path.exists():
        log.info(f"[SKIP] 載入已存在的 bin attribution: {attr_path.name}")
        return pd.read_parquet(str(attr_path))

    log.info("[3/4] Bin attribution")
    tp = pd.read_parquet(str(tp_path), columns=[
        "barcode", "in_tissue", "pxl_row_in_fullres", "pxl_col_in_fullres"
    ])
    tp = tp[tp["in_tissue"] == 1]
    H, W = mask.shape
    row_local = (tp["pxl_row_in_fullres"].values - crop_y0).astype(np.int32).clip(0, H - 1)
    col_local = (tp["pxl_col_in_fullres"].values - btf_col0).astype(np.int32).clip(0, W - 1)
    tp = tp.copy()
    tp["cell_id"] = mask[row_local, col_local]
    attr = tp[tp["cell_id"] > 0][["barcode", "cell_id"]].reset_index(drop=True)
    log.info(f"  attributed bins: {len(attr):,} / {len(tp):,} ({len(attr)/len(tp):.1%})")
    attr.to_parquet(str(attr_path), index=False)
    log.info(f"  儲存: {attr_path.name}")
    return attr


def _aggregate_cells_raw(
    attribution: "pd.DataFrame",  # noqa: F821
    h5_path: Path,
) -> "ad.AnnData":  # noqa: F821
    """
    依 attribution（barcode→cell_id）把 2µm bins 聚合成 cells×genes 原始 counts。

    回傳 AnnData：X=原始 counts（稀疏），obs_names=cell_id 字串，
    obs['cell_id']（int）、obs['n_bins'］。不做 normalize（供下游自由運用）。
    """
    import gc as _gc

    import anndata as ad
    import numpy as np
    import scanpy as sc
    import scipy.sparse as sp

    log.info(f"  讀取 h5 矩陣: {h5_path.name}")
    adata_full = sc.read_10x_h5(str(h5_path))
    adata_full.var_names_make_unique()

    mask_obs = adata_full.obs_names.isin(attribution["barcode"].values)
    adata_crop = adata_full[mask_obs].copy()
    del adata_full
    _gc.collect()

    barcode_to_cell = attribution.set_index("barcode")["cell_id"]
    cell_ids = barcode_to_cell.reindex(adata_crop.obs_names).values.astype(np.int32)
    valid = cell_ids > 0
    adata_valid = adata_crop[valid]
    cell_ids_v = cell_ids[valid]
    unique_cells = np.unique(cell_ids_v)
    n_cells = len(unique_cells)
    log.info(f"  unique cells with RNA: {n_cells:,}")

    lut = np.zeros(int(unique_cells.max()) + 1, dtype=np.int64)
    lut[unique_cells] = np.arange(n_cells)
    rows = lut[cell_ids_v]
    cols = np.arange(len(cell_ids_v))
    A = sp.csr_matrix(
        (np.ones(len(cell_ids_v), dtype=np.float32), (rows, cols)),
        shape=(n_cells, adata_valid.n_obs),
    )
    X_agg = A @ adata_valid.X
    cells = ad.AnnData(
        X=X_agg.tocsr() if sp.issparse(X_agg) else sp.csr_matrix(X_agg),
        var=adata_valid.var.copy(),
    )
    cells.obs_names = [str(int(c)) for c in unique_cells]
    cells.obs["cell_id"] = unique_cells.astype(int)
    cells.obs["n_bins"] = np.asarray(A.sum(axis=1)).ravel().astype(int)
    del adata_crop, adata_valid, A
    _gc.collect()
    return cells


def step_count_cells(
    mask: np.ndarray,
    attribution: "pd.DataFrame",  # noqa: F821
    h5_path: Path,
    out_dir: Path,
    pixel_size_um: float,
) -> Path:
    """[4/6] 由 bin attribution 聚合 cells×genes，附加重心，輸出 cells.h5ad。"""
    import numpy as np
    from scipy.ndimage import center_of_mass

    h5ad_path = out_dir / "cells.h5ad"
    if h5ad_path.exists():
        log.info(f"[SKIP] 載入已存在的 cells h5ad: {h5ad_path.name}")
        return h5ad_path

    log.info("[4/6] 聚合 cells×genes 矩陣")
    cells = _aggregate_cells_raw(attribution, h5_path)
    unique_cells = cells.obs["cell_id"].values.astype(np.int64)

    # 重心（mask 局部 px，原點 = 裁切左上角）→ µm
    log.info("  計算細胞重心…")
    cen = np.asarray(
        center_of_mass(mask > 0, labels=mask, index=unique_cells.tolist()),
        dtype=float,
    )
    cy_px, cx_px = cen[:, 0], cen[:, 1]
    cells.obs["centroid_x_px"] = cx_px
    cells.obs["centroid_y_px"] = cy_px
    cells.obsm["spatial"] = np.stack(
        [cx_px * pixel_size_um, cy_px * pixel_size_um], axis=1
    )

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

    import numpy as np
    from skimage import measure

    xen_dir = out_dir / "xenium_explorer"
    if (xen_dir / "experiment.xenium").exists():
        log.info(f"[SKIP] Xenium bundle 已存在: {xen_dir.name}")
        return xen_dir

    # 1. 細胞多邊形 GeoJSON（局部 µm，原點 = 裁切左上角；格式同 GUI 匯出）
    geojson_path = out_dir / "cells_polygons.geojson"
    if geojson_path.exists():
        log.info(f"[SKIP] 載入已存在的多邊形: {geojson_path.name}")
    else:
        n_cells = int(mask.max())
        log.info(f"[6/6] 產生細胞多邊形 GeoJSON（{n_cells:,} cells，整片可能較久）…")
        features = []
        for prop in measure.regionprops(mask):
            cid = prop.label
            r0, c0, r1, c1 = prop.bbox
            cell_crop = (mask[r0:r1, c0:c1] == cid).astype(np.uint8)
            contours = measure.find_contours(np.pad(cell_crop, 1, mode="constant"), 0.5)
            if not contours:
                continue
            contour = max(contours, key=len)
            xy_um = np.column_stack([
                (contour[:, 1] - 1 + c0) * pixel_size_um,   # col → x
                (contour[:, 0] - 1 + r0) * pixel_size_um,   # row → y
            ])
            if not np.allclose(xy_um[0], xy_um[-1]):
                xy_um = np.vstack([xy_um, xy_um[0]])
            features.append({
                "type": "Feature",
                "geometry": {"type": "Polygon", "coordinates": [xy_um.tolist()]},
                "properties": {"full_id": str(int(cid)), "cell_id": int(cid)},
            })
        with open(geojson_path, "w", encoding="utf-8") as f:
            json.dump({"type": "FeatureCollection", "features": features}, f)
        log.info(f"  多邊形數: {len(features):,} → {geojson_path.name}")

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
    seg.add_argument("--tissue",  choices=list(TISSUE_PRESETS), default="crc",
                     help="組織類型 preset（crc / luad / default），預設 crc")
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
    seg.add_argument("--dia-small",  type=float, metavar="PX",
                     help="cyto3 小直徑（覆寫 tissue preset）")
    seg.add_argument("--dia-mid",    type=float, metavar="PX",
                     help="cyto3 中直徑（覆寫 tissue preset）")
    seg.add_argument("--dia-large",  type=float, metavar="PX",
                     help="cyto3 大直徑（覆寫 tissue preset）")
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

    # ── 建立分割設定
    cfg = dict(TISSUE_PRESETS[args.tissue])
    cfg["use_gpu"]              = not args.no_gpu
    cfg["batch_size"]           = args.batch_size
    cfg["use_cpsam"]            = args.cpsam
    cfg["use_hematoxylin"]      = True
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
    if args.he_crop:
        img = step_load_he_crop(args.he_crop)
        crop_y0  = 0
        btf_col0 = 0
    else:
        import tifffile
        if args.crop_y1 == -1 or args.btf_col1 == -1:
            with tifffile.TiffFile(str(args.btf)) as tif:
                full_shape = tif.pages[0].shape
            crop_y1  = full_shape[0] if args.crop_y1  == -1 else args.crop_y1
            btf_col1 = full_shape[1] if args.btf_col1 == -1 else args.btf_col1
        else:
            crop_y1  = args.crop_y1
            btf_col1 = args.btf_col1
        crop_y0  = args.crop_y0
        btf_col0 = args.btf_col0
        img = step_crop_btf(args.btf, out_dir, crop_y0, crop_y1, btf_col0, btf_col1)

    # ── Step 2: 分割
    mask = step_segment(img, cfg, out_dir)
    del img
    gc.collect()

    # ── Step 3: Bin attribution（有 tp & h5 才跑）
    from backend.src.utils.constants import VISIUM_UM_PX
    pixel_size_um = VISIUM_UM_PX

    attribution = None
    if args.tp and args.h5:
        if not args.tp.exists():
            log.warning(f"tissue_positions 不存在，跳過 RNA 計數: {args.tp}")
        elif not args.h5.exists():
            log.warning(f"h5 矩陣不存在，跳過 RNA 計數: {args.h5}")
        else:
            attribution = step_bin_attribution(
                mask, args.tp, out_dir, crop_y0, btf_col0
            )

    # ── Step 4: 聚合 cells×genes h5ad（有 attribution 才跑）
    cells_h5ad = None
    if attribution is not None:
        cells_h5ad = step_count_cells(
            mask, attribution, args.h5, out_dir, pixel_size_um
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
