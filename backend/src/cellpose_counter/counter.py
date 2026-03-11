"""
Stage 2：Cellpose RNA 計數器

直接將 Visium HD 2µm bins 分配給 Cellpose 分割細胞（跳過 Proseg）。

核心邏輯：
  1. 載入 adata_002um.h5ad（bins × genes，obsm['spatial'] = 全域 fullres px 座標）
  2. 載入 segmentation_masks.npy（H×W，pixel 值 = 細胞 ID）
  3. 從 pipeline.yaml ROI 設定取得 ROI 裁切偏移（x_px, y_px）
  4. 對每個 bin，將座標換算到 ROI 像素空間，查詢對應細胞 ID
  5. 以稀疏矩陣乘法匯總每個細胞的 gene counts
  6. 計算細胞重心（scipy.ndimage.center_of_mass）
  7. 儲存為 cellpose_cells.h5ad（cells × genes）

輸出欄位：
  obs: cell_id, n_bins, centroid_x_px, centroid_y_px,
       centroid_x_um, centroid_y_um, cell_area_px, cell_area_um2
  obsm['spatial']: [x_um, y_um]（ROI 局部座標，原點 = ROI 左上角）
  var: 同 adata_002um
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

import numpy as np

from backend.src.utils.constants import VISIUM_UM_PX

logger = logging.getLogger("pipeline.cellpose_counter")

# ──────────────────────────────────────────────────────────────────────────────
# 核心計數函式
# ──────────────────────────────────────────────────────────────────────────────

def count_rna_per_cell(
    adata_path: Path,
    mask_path: Path,
    roi_x_px: int,
    roi_y_px: int,
    pixel_size_um: float = VISIUM_UM_PX,
) -> "anndata.AnnData":
    """
    將 Visium HD 2µm bins 分配給 Cellpose 細胞，回傳 cells × genes AnnData。

    Parameters
    ----------
    adata_path   : adata_002um.h5ad 路徑
    mask_path    : segmentation_masks.npy 路徑
    roi_x_px     : ROI 左上角 X（fullres px，來自 pipeline.yaml）
    roi_y_px     : ROI 左上角 Y（fullres px）
    pixel_size_um: fullres 像素尺寸（µm/px）
    """
    import anndata as ad
    import scanpy as sc
    import scipy.sparse as sp
    from scipy.ndimage import center_of_mass, labeled_comprehension

    logger.info(f"載入 bins AnnData：{adata_path}")
    adata = sc.read_h5ad(str(adata_path))
    logger.info(f"  bins n_obs={adata.n_obs}, n_vars={adata.n_vars}")

    logger.info(f"載入 Cellpose 遮罩：{mask_path}")
    seg_mask = np.load(str(mask_path))
    H, W = seg_mask.shape
    logger.info(f"  遮罩大小：{W}×{H}，最大 cell ID：{seg_mask.max()}")

    # ── 1. 取得 bin 中心座標（全域 fullres px）─────────────────────────────
    if "spatial" in adata.obsm:
        coords = adata.obsm["spatial"]  # shape (N, 2)：[col/x, row/y]
    elif "pxl_col_in_fullres" in adata.obs.columns and "pxl_row_in_fullres" in adata.obs.columns:
        logger.info("  obsm['spatial'] 不存在，從 obs 欄位建立座標")
        coords = np.stack([
            adata.obs["pxl_col_in_fullres"].values.astype(float),
            adata.obs["pxl_row_in_fullres"].values.astype(float),
        ], axis=1)
    else:
        raise ValueError("adata_002um.h5ad 缺少 obsm['spatial'] 與 obs 空間座標欄位")

    # ── 2. 換算至 ROI 局部像素座標 ───────────────────────────────────────
    roi_col = (coords[:, 0] - roi_x_px).round().astype(int)   # x → col
    roi_row = (coords[:, 1] - roi_y_px).round().astype(int)   # y → row

    # ── 3. 篩選有效範圍（在 ROI 內）──────────────────────────────────────
    valid = (
        (roi_col >= 0) & (roi_col < W) &
        (roi_row >= 0) & (roi_row < H)
    )
    n_valid = valid.sum()
    n_out = len(adata) - n_valid
    logger.info(f"  有效 bins：{n_valid}（ROI 內），{n_out} 個超出 ROI 範圍（忽略）")
    if n_out > 0 and n_out / len(adata) > 0.3:
        logger.warning(
            f"⚠️ 超出 ROI 範圍的 bins 超過 30%（{n_out}/{len(adata)}），"
            f"請確認 ROI 偏移 (x={roi_x_px}, y={roi_y_px}) 是否正確"
        )

    # ── 4. 查詢每個 bin 的細胞 ID ────────────────────────────────────────
    bin_cell_ids = np.zeros(len(adata), dtype=np.int32)
    valid_idx = np.where(valid)[0]
    bin_cell_ids[valid_idx] = seg_mask[roi_row[valid_idx], roi_col[valid_idx]]

    # 只保留被分配到細胞（cell_id > 0）的 bins
    in_cell_mask = bin_cell_ids > 0
    n_assigned = in_cell_mask.sum()
    logger.info(f"  分配至細胞的 bins：{n_assigned}（{n_assigned/len(adata)*100:.1f}%）")

    # ── 5. 取得所有唯一細胞 ID（按 mask 中出現順序排序）─────────────────
    unique_cells = np.unique(seg_mask[seg_mask > 0])
    n_cells = len(unique_cells)
    logger.info(f"  Cellpose 細胞數：{n_cells}")
    cell_id_to_row = {int(cid): i for i, cid in enumerate(unique_cells)}

    # ── 6. 建立分配矩陣 A（n_cells × n_bins），稀疏 ──────────────────────
    assigned_bin_idx = np.where(in_cell_mask)[0]
    assigned_cell_ids = bin_cell_ids[in_cell_mask]
    row_idx = np.array([cell_id_to_row[int(cid)] for cid in assigned_cell_ids], dtype=np.int32)

    A = sp.csr_matrix(
        (np.ones(len(row_idx), dtype=np.float32), (row_idx, assigned_bin_idx)),
        shape=(n_cells, adata.n_obs),
    )

    # ── 7. 矩陣乘法：(n_cells × n_bins) @ (n_bins × n_genes) ────────────
    logger.info("  計算 gene count 矩陣...")
    X_source = adata.X if sp.issparse(adata.X) else sp.csr_matrix(adata.X)
    X_cells = A @ X_source  # → sparse (n_cells × n_genes)
    if not sp.issparse(X_cells):
        X_cells = sp.csr_matrix(X_cells)

    # ── 8. 計算每細胞 bin 數量 ───────────────────────────────────────────
    n_bins_per_cell = np.array(A.sum(axis=1)).ravel().astype(int)

    # ── 9. 計算細胞重心（ROI 局部座標）──────────────────────────────────
    logger.info("  計算細胞重心...")
    centroids_yx = center_of_mass(
        seg_mask > 0,
        labels=seg_mask,
        index=unique_cells.tolist(),
    )
    centroids = np.array(centroids_yx, dtype=float)  # (n_cells, 2): (row, col)
    cx_px = centroids[:, 1]   # col = x（ROI 局部）
    cy_px = centroids[:, 0]   # row = y（ROI 局部）
    cx_um = cx_px * pixel_size_um
    cy_um = cy_px * pixel_size_um

    # ── 10. 計算細胞面積 ─────────────────────────────────────────────────
    area_px = labeled_comprehension(
        np.ones_like(seg_mask, dtype=np.int32),
        seg_mask,
        unique_cells.tolist(),
        func=np.sum,
        out_dtype=float,
        default=0.0,
    )
    area_um2 = area_px * (pixel_size_um ** 2)

    # ── 11. 組建輸出 AnnData ─────────────────────────────────────────────
    import pandas as pd

    obs = pd.DataFrame(
        {
            "cell_id":         unique_cells.astype(int),
            "n_bins":          n_bins_per_cell,
            "centroid_x_px":   cx_px,
            "centroid_y_px":   cy_px,
            "centroid_x_um":   cx_um,
            "centroid_y_um":   cy_um,
            "cell_area_px":    area_px,
            "cell_area_um2":   area_um2,
        },
        index=[f"cell_{int(cid)}" for cid in unique_cells],
    )

    result = ad.AnnData(X=X_cells, obs=obs, var=adata.var.copy())

    # obsm['spatial']：ROI 局部 µm 座標，[x_um, y_um]
    result.obsm["spatial"] = np.stack([cx_um, cy_um], axis=1)

    logger.info(
        f"  完成：{n_cells} 個細胞，"
        f"中位 genes/cell = {np.median((X_cells > 0).sum(axis=1)).astype(int)}，"
        f"中位 counts/cell = {int(np.median(np.array(X_cells.sum(axis=1)).ravel()))}"
    )
    return result


# ──────────────────────────────────────────────────────────────────────────────
# 多 ROI Pipeline（供 API 層呼叫）
# ──────────────────────────────────────────────────────────────────────────────

def run_counting_pipeline(config: dict, roi_name: Optional[str] = None):
    """
    對指定 ROI（或全部 ROI）執行 Cellpose RNA 計數，儲存 cellpose_cells.h5ad。

    Parameters
    ----------
    config   : pipeline.yaml 設定字典
    roi_name : 若指定，只跑單一 ROI；否則跑全部
    """
    from backend.src.utils.config import resolve_path

    paths = config.get("paths", {})
    rois  = config.get("rois", [])
    out_base = resolve_path(paths.get("output_dir", "results/analysis")) / "roi"

    if roi_name:
        rois = [r for r in rois if r.get("name") == roi_name]
        if not rois:
            raise ValueError(f"找不到 ROI：{roi_name}")

    if not rois:
        raise ValueError("pipeline.yaml 未設定任何 ROI")

    for roi in rois:
        rn = roi.get("name", "")
        if not rn:
            logger.warning("跳過無名稱的 ROI 設定")
            continue

        roi_dir = out_base / rn
        adata_path = roi_dir / "adata_002um.h5ad"
        mask_path  = roi_dir / "segmentation_masks.npy"
        out_path   = roi_dir / "cellpose_cells.h5ad"

        if not adata_path.exists():
            logger.warning(f"[{rn}] 找不到 adata_002um.h5ad，跳過（請先完成 Stage 0）")
            continue
        if not mask_path.exists():
            logger.warning(f"[{rn}] 找不到 segmentation_masks.npy，跳過（請先完成 Stage 1）")
            continue

        # 取得 ROI 裁切偏移
        roi_x_px = int(roi.get("x", 0))
        roi_y_px = int(roi.get("y", 0))
        pixel_size_um = float(roi.get("pixel_size_um", VISIUM_UM_PX))

        logger.info(f"[{rn}] 開始計數（ROI 偏移 x={roi_x_px}, y={roi_y_px}）...")
        adata_cells = count_rna_per_cell(
            adata_path=adata_path,
            mask_path=mask_path,
            roi_x_px=roi_x_px,
            roi_y_px=roi_y_px,
            pixel_size_um=pixel_size_um,
        )

        roi_dir.mkdir(parents=True, exist_ok=True)
        adata_cells.write_h5ad(str(out_path))
        logger.info(f"[{rn}] 儲存完成：{out_path}（{adata_cells.n_obs} cells）")

    logger.info("所有 ROI 計數完成")
