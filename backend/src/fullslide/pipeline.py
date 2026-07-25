"""
全片流程核心函式（CLI 與 API 共用）
=====================================

本模組收斂「MCseg v2 遮罩 + Visium HD bins → cells×genes」這條配方，
使 `cli/segment.py`（指令列全片流程）與 `api/cellpose_count.py`（GUI 全圖計數）
共用同一份實作，避免邏輯分歧（CLAUDE.md §11 DRY）。

函式一律為純函式（不寫 log、不做快取），快取與進度回報留給呼叫端：

| 函式 | 職責 |
|------|------|
| `bin_attribution`   | 2µm bins → cell_id 對應表 |
| `aggregate_cells`   | 依對應表把 bins 聚合成 cells×genes 原始 counts |
| `add_centroids`     | 補細胞重心（裁切局部 px、全片 fullres px、µm） |
| `resolve_pixel_size`| 取樣本實際 µm/px（scalefactors 優先於預設常數） |
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np

logger = logging.getLogger("pipeline.fullslide")

if TYPE_CHECKING:  # pragma: no cover - 僅供型別檢查，避免匯入期拉進重量級套件
    import anndata as ad
    import pandas as pd


def resolve_bin_to_mask_scale(
    config: dict, mask_shape: tuple[int, int]
) -> tuple[float, float]:
    """
    推導 Space Ranger fullres px → 遮罩(raw TIFF) px 的**分軸**縮放。

    為何需要：Space Ranger 的 "fullres" 座標系是校準到**餵給它的那張影像**，
    未必等於 MSseg 分割所用的 raw TIFF。例如 dpcp01 樣本的 SR fullres 是
    0.5464 µm/px 而 TIFF 約 0.2737 → 相差近 2 倍。若不縮放，每個 bin 都會
    落在約一半的位置，**不報錯、不警告**（EP 2026-07-21 曾因此重跑計數）。

    推導不變量：`tissue_hires_image.png` 尺寸 ÷ `tissue_hires_scalef`
    ＝ SR fullres 畫布尺寸。此法算出的 (11266, 23552) 與 EP 已驗證的
    `segmentation_masks_fullslide_vfr_CORRECTED.npy` shape 完全一致。

    **必須分軸**：實測 dpcp01 為 col 1.9088、row 2.0000（差 4.6%）。強套等向
    2.0 會把最右側 bin 映射到 22484 px，超出 TIFF 寬度 21504 近千像素。

    TIFF 的 `XResolution` 標籤不可用作依據 —— 實測為 96 dpi 的通用預設值。

    Parameters
    ----------
    config : dict
        需含 `paths.binned_002`；`alignment.bin_to_mask_scale` 若指定 `[sx, sy]`
        則優先採用（供自動推導失準時人工覆寫）。
    mask_shape : tuple[int, int]
        遮罩的 `(height, width)`，即 raw TIFF 的像素尺寸。

    Returns
    -------
    tuple[float, float]
        `(scale_x, scale_y)`；無法推導時回退 `(1.0, 1.0)`。
    """
    import json

    override = (config.get("alignment") or {}).get("bin_to_mask_scale")
    if override:
        return float(override[0]), float(override[1])

    binned = config.get("paths", {}).get("binned_002", "")
    if not binned:
        return 1.0, 1.0

    spatial = Path(binned) / "spatial"
    sf_path = spatial / "scalefactors_json.json"
    hires_path = spatial / "tissue_hires_image.png"
    if not (sf_path.exists() and hires_path.exists()):
        return 1.0, 1.0

    try:
        scalef = float(json.loads(sf_path.read_text(encoding="utf-8"))["tissue_hires_scalef"])
        from PIL import Image

        Image.MAX_IMAGE_PIXELS = None   # hires 圖仍可能觸發 decompression bomb 警戒
        with Image.open(hires_path) as im:
            hires_w, hires_h = im.size
    except (ValueError, OSError, KeyError, TypeError):
        return 1.0, 1.0

    if scalef <= 0 or hires_w <= 0 or hires_h <= 0:
        return 1.0, 1.0

    sr_w, sr_h = hires_w / scalef, hires_h / scalef
    mask_h, mask_w = mask_shape
    return mask_w / sr_w, mask_h / sr_h


def bin_attribution(
    mask: np.ndarray,
    tp_path: str | Path,
    crop_y0: int,
    crop_x0: int,
    out_path: str | Path | None = None,
    scale: tuple[float, float] = (1.0, 1.0),
    transform: np.ndarray | None = None,
) -> "pd.DataFrame":
    """
    將 Visium HD 2µm bins 對應到 MCseg v2 細胞遮罩。

    座標換算順序：**先變換、後扣原點**（原點是遮罩空間的量）。

    ```text
    SR fullres px ──transform / ×scale──> 影像(raw TIFF) px ──−crop origin──> 遮罩局部 px
    ```

    Parameters
    ----------
    mask : np.ndarray
        (H, W) int 遮罩，像素值 = 細胞 ID，0 為背景。座標原點為裁切左上角。
    tp_path : str | Path
        `tissue_positions.parquet` 路徑。
    crop_y0, crop_x0 : int
        裁切左上角在 raw TIFF 座標系中的位置。
    out_path : str | Path | None
        給定時額外寫出 parquet 快取。
    scale : tuple[float, float]
        `(scale_x, scale_y)`，SR fullres → 遮罩 px 的分軸縮放。
        由 `resolve_bin_to_mask_scale` 取得；`(1.0, 1.0)` 表示兩者同座標系。
    transform : np.ndarray | None
        3×3 homography（作用於 `(x, y) = (col, row)` 齊次座標），由
        `registration.alignment.compose_bin_to_image` 取得。這是**幾何正確**
        的路徑；`scale` 僅為找不到對位 JSON 時的近似回退。兩者同時給定時
        `transform` 優先並記 warning。

    Returns
    -------
    pd.DataFrame
        欄位 `barcode`、`cell_id`；只含 `in_tissue == 1`、落在遮罩範圍內、
        且落在細胞內（`cell_id > 0`）的 bins。
        `.attrs["coverage"]` 記錄涵蓋率（見 `_coverage_attrs`）。

    Notes
    -----
    越界的 bin 會被**排除**而非夾到邊緣 —— 夾邊會把界外 bin 的 RNA 誤記到
    邊界細胞上，且不留痕跡。落在影像外是真實限制（dpcp01 正確變換下 SR 右緣
    約 980 px 寬確實不在高解析 TIFF 內），必須如實回報而非壓縮硬塞。
    """
    import pandas as pd

    tp = pd.read_parquet(
        str(tp_path),
        columns=["barcode", "in_tissue", "pxl_row_in_fullres", "pxl_col_in_fullres"],
    )
    tp = tp[tp["in_tissue"] == 1].copy()

    h, w = mask.shape
    src_row = tp["pxl_row_in_fullres"].values.astype(float)
    src_col = tp["pxl_col_in_fullres"].values.astype(float)

    if transform is not None:
        if scale != (1.0, 1.0):
            logger.warning(
                f"同時指定 transform 與 scale={scale}；採用 transform（幾何正確），忽略 scale。"
            )
        x, y = _apply_homography(np.asarray(transform, dtype=float), src_col, src_row)
        desc = "homography"
    else:
        x, y = src_col * scale[0], src_row * scale[1]
        desc = f"scale={scale}"

    row = np.rint(y - crop_y0).astype(np.int64)
    col = np.rint(x - crop_x0).astype(np.int64)

    in_bounds = (row >= 0) & (row < h) & (col >= 0) & (col < w)
    cell_ids = np.zeros(len(tp), dtype=mask.dtype)
    cell_ids[in_bounds] = mask[row[in_bounds], col[in_bounds]]
    tp["cell_id"] = cell_ids

    n_total = len(tp)
    n_oob = int((~in_bounds).sum())
    if n_oob and n_total and n_oob / n_total > 0.3:
        logger.warning(
            f"⚠️ {n_oob:,}/{n_total:,}（{n_oob/n_total:.1%}）個 bin 落在遮罩範圍外。"
            f"{desc}、crop=({crop_x0}, {crop_y0})、遮罩 {w}×{h}px —— "
            f"請確認 binned_outputs 是否為對應此影像的 CytAssist 註冊版本。"
        )

    attr = tp[tp["cell_id"] > 0][["barcode", "cell_id"]].reset_index(drop=True)
    attr.attrs["coverage"] = {
        "n_total": n_total,
        "n_in_bounds": int(in_bounds.sum()),
        "n_assigned": int(len(attr)),
        "frac_out_of_image": (n_oob / n_total) if n_total else 0.0,
    }

    if out_path is not None:
        attr.to_parquet(str(out_path), index=False)
    return attr


def _apply_homography(
    m: np.ndarray, x: np.ndarray, y: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """對 `(x, y)` 套用 3×3 homography（含齊次除法）。"""
    denom = m[2, 0] * x + m[2, 1] * y + m[2, 2]
    denom = np.where(denom == 0, np.nan, denom)   # 退化點 → NaN，後續判界時自然排除
    xt = (m[0, 0] * x + m[0, 1] * y + m[0, 2]) / denom
    yt = (m[1, 0] * x + m[1, 1] * y + m[1, 2]) / denom
    return np.nan_to_num(xt, nan=-1.0), np.nan_to_num(yt, nan=-1.0)


def aggregate_cells(attribution: "pd.DataFrame", h5_path: str | Path) -> "ad.AnnData":
    """
    依 attribution（barcode → cell_id）把 2µm bins 聚合成 cells×genes 原始 counts。

    不做 normalize —— 保留原始 counts 供下游自由運用。

    Parameters
    ----------
    attribution : pd.DataFrame
        `bin_attribution` 的輸出（欄位 `barcode`、`cell_id`）。
    h5_path : str | Path
        `filtered_feature_bc_matrix.h5` 路徑。

    Returns
    -------
    ad.AnnData
        X = 原始 counts（稀疏 CSR），`obs_names` 為 cell_id 字串，
        另含 `obs['cell_id']`（int）與 `obs['n_bins']`。
    """
    import gc

    import anndata as ad
    import scanpy as sc
    import scipy.sparse as sp

    adata_full = sc.read_10x_h5(str(h5_path))
    adata_full.var_names_make_unique()

    keep = adata_full.obs_names.isin(attribution["barcode"].values)
    adata_crop = adata_full[keep].copy()
    del adata_full
    gc.collect()

    barcode_to_cell = attribution.set_index("barcode")["cell_id"]
    cell_ids = barcode_to_cell.reindex(adata_crop.obs_names).values.astype(np.int32)
    valid = cell_ids > 0
    adata_valid = adata_crop[valid]
    cell_ids_v = cell_ids[valid]

    unique_cells = np.unique(cell_ids_v)
    n_cells = len(unique_cells)

    # 向量化 LUT：O(max_id) 建立、O(n) 查詢，比 dict 快 10-100x
    lut = np.zeros(int(unique_cells.max()) + 1, dtype=np.int64)
    lut[unique_cells] = np.arange(n_cells)
    rows = lut[cell_ids_v]
    cols = np.arange(len(cell_ids_v))

    A = sp.csr_matrix(
        (np.ones(len(cell_ids_v), dtype=np.float32), (rows, cols)),
        shape=(n_cells, adata_valid.n_obs),
    )
    x_agg = A @ adata_valid.X

    cells = ad.AnnData(
        X=x_agg.tocsr() if sp.issparse(x_agg) else sp.csr_matrix(x_agg),
        var=adata_valid.var.copy(),
    )
    cells.obs_names = [str(int(c)) for c in unique_cells]
    cells.obs["cell_id"] = unique_cells.astype(int)
    cells.obs["n_bins"] = np.asarray(A.sum(axis=1)).ravel().astype(int)

    del adata_crop, adata_valid, A
    gc.collect()
    return cells


def add_centroids(
    cells: "ad.AnnData",
    mask: np.ndarray,
    pixel_size_um: float,
    origin_xy: tuple[int, int] = (0, 0),
) -> None:
    """
    就地補上細胞重心（三種座標系並存）。

    | 欄位 | 座標系 | 用途 |
    |------|--------|------|
    | `obs['centroid_x_px']` / `_y_px`           | 裁切局部 px | 遮罩內定位（維持 CLI 既有語意） |
    | `obs['centroid_x_fullres']` / `_y_fullres` | 全片 fullres px | 跨 ROI／全片框選（Stage 3.5 區域選取） |
    | `obsm['spatial']`                          | µm（由局部換算） | 匯出與繪圖（維持既有語意） |

    fullres 欄位在此處落地，而非留給下游每次讀 metadata sidecar 自行補償 ——
    後者是重複且易錯的轉換（座標系錯配會讓框選偏移一整個裁切原點）。

    Parameters
    ----------
    cells : ad.AnnData
        `aggregate_cells` 的輸出，需含 `obs['cell_id']`。就地修改。
    mask : np.ndarray
        (H, W) 遮罩，座標原點為裁切左上角。
    pixel_size_um : float
        µm/px（建議由 `resolve_pixel_size` 取得樣本實際值）。
    origin_xy : tuple[int, int]
        裁切左上角於原始影像 fullres 座標系的位置 `(x0, y0)`。
    """
    from scipy.ndimage import center_of_mass

    unique_cells = cells.obs["cell_id"].values.astype(np.int64)
    cen = np.asarray(
        center_of_mass(mask > 0, labels=mask, index=unique_cells.tolist()),
        dtype=float,
    )
    cy_px, cx_px = cen[:, 0], cen[:, 1]
    x0, y0 = origin_xy

    cells.obs["centroid_x_px"] = cx_px
    cells.obs["centroid_y_px"] = cy_px
    cells.obs["centroid_x_fullres"] = cx_px + x0
    cells.obs["centroid_y_fullres"] = cy_px + y0
    cells.obsm["spatial"] = np.stack([cx_px * pixel_size_um, cy_px * pixel_size_um], axis=1)


# ── 全圖遮罩 metadata sidecar ────────────────────────────────────────────────

FULL_SEG_MASK_FILENAME = "full_image_segmentation_masks.npy"
FULL_SEG_META_FILENAME = "full_image_segmentation_meta.json"


def write_full_seg_meta(
    output_dir: str | Path,
    *,
    crop_x0: int,
    crop_y0: int,
    width: int,
    height: int,
    n_cells: int,
    pixel_size_um: float,
    passes: int,
) -> Path:
    """
    寫出全圖遮罩的 metadata sidecar。

    遮罩本身只有局部座標，下游（Stage 2 全圖計數、Stage 3.5 框選）必須靠
    `crop_x0/crop_y0` 才能還原回原始影像 fullres 座標系。
    """
    import json
    from datetime import datetime, timezone

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    meta = {
        "crop_x0": int(crop_x0),
        "crop_y0": int(crop_y0),
        "width": int(width),
        "height": int(height),
        "n_cells": int(n_cells),
        "pixel_size_um": float(pixel_size_um),
        "passes": int(passes),
        "created_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
    }
    path = output_dir / FULL_SEG_META_FILENAME
    path.write_text(json.dumps(meta, indent=2), encoding="utf-8")
    return path


def read_full_seg_meta(output_dir: str | Path) -> dict | None:
    """讀取 metadata sidecar；不存在或損壞時回傳 None（由呼叫端決定如何處理）。"""
    import json

    path = Path(output_dir) / FULL_SEG_META_FILENAME
    if not path.exists():
        return None
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (ValueError, OSError):
        return None


def resolve_full_count_inputs(config: dict) -> tuple[dict | None, str | None]:
    """
    解析全圖 RNA 計數所需的輸入，回傳 `(inputs, error)`。

    成功時 `inputs` 含 `mask_path`、`tp_path`、`h5_path`、`origin_xy`、
    `pixel_size_um`、`meta_missing`；失敗時 `inputs` 為 None 且 `error` 為
    可直接回傳給前端的訊息（不含絕對路徑）。

    `meta_missing=True` 代表找不到 sidecar，原點退為 (0, 0)。這是**警示而非
    錯誤** —— 舊版遮罩沒有 sidecar，且全圖模式的原點本來就是 (0, 0)。
    """
    from backend.src.utils.config import resolve_path

    paths = config.get("paths", {})
    output_dir = resolve_path(paths["output_dir"])
    mask_path = output_dir / FULL_SEG_MASK_FILENAME
    if not mask_path.exists():
        return None, "找不到全圖分割遮罩，請先完成全圖分割（Stage 1）"

    binned = paths.get("binned_002", "")
    tp_path = Path(binned) / "spatial" / "tissue_positions.parquet"
    if not tp_path.exists():
        return None, "找不到 tissue_positions.parquet，請確認 paths.binned_002 設定"
    h5_path = Path(binned) / "filtered_feature_bc_matrix.h5"
    if not h5_path.exists():
        return None, "找不到 filtered_feature_bc_matrix.h5，請確認 paths.binned_002 設定"

    meta = read_full_seg_meta(output_dir)
    origin_xy = (0, 0) if meta is None else (int(meta["crop_x0"]), int(meta["crop_y0"]))

    # 遮罩尺寸只讀 .npy header（mmap 不載入陣列本體，全片可達數 GB）
    mask_shape = np.load(str(mask_path), mmap_mode="r").shape

    return {
        "mask_path": mask_path,
        "tp_path": tp_path,
        "h5_path": h5_path,
        "origin_xy": origin_xy,
        "scale": resolve_bin_to_mask_scale(config, mask_shape),
        # 與 ROI 路徑（counter.py）採同一設定，否則全圖與 ROI 結果不可比
        "dilation_px": int((config.get("rna_counting") or {}).get("dilation_px", 0)),
        "pixel_size_um": resolve_pixel_size(config),
        "meta_missing": meta is None,
        "output_dir": output_dir,
    }, None


def resolve_pixel_size(config: dict) -> float:
    """
    取樣本實際的 µm/px：`scalefactors_json.json` 的 `microns_per_pixel` 優先，
    缺失／損壞時回退 `constants.VISIUM_UM_PX`。

    注意：此值只影響 µm 換算與匯出比例尺，**不影響 bin attribution**
    （後者為純像素運算）。
    """
    import json

    from backend.src.utils.constants import VISIUM_UM_PX

    binned = config.get("paths", {}).get("binned_002", "")
    if binned:
        sf = Path(binned) / "spatial" / "scalefactors_json.json"
        if sf.exists():
            try:
                mpp = json.loads(sf.read_text(encoding="utf-8")).get("microns_per_pixel")
                if mpp:
                    return float(mpp)
            except (ValueError, OSError, TypeError):
                pass   # 損壞不可中斷流程，回退預設常數
    return float(VISIUM_UM_PX)
