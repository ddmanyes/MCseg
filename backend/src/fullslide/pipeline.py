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

from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:  # pragma: no cover - 僅供型別檢查，避免匯入期拉進重量級套件
    import anndata as ad
    import pandas as pd


def bin_attribution(
    mask: np.ndarray,
    tp_path: str | Path,
    crop_y0: int,
    crop_x0: int,
    out_path: str | Path | None = None,
) -> "pd.DataFrame":
    """
    將 Visium HD 2µm bins 對應到 MCseg v2 細胞遮罩。

    Parameters
    ----------
    mask : np.ndarray
        (H, W) int 遮罩，像素值 = 細胞 ID，0 為背景。座標原點為裁切左上角。
    tp_path : str | Path
        `tissue_positions.parquet` 路徑。
    crop_y0, crop_x0 : int
        裁切左上角在**原始影像 fullres 座標系**中的位置。bin 的
        `pxl_row/col_in_fullres` 會扣除此原點後才查表。
    out_path : str | Path | None
        給定時額外寫出 parquet 快取。

    Returns
    -------
    pd.DataFrame
        欄位 `barcode`、`cell_id`；只含 `in_tissue == 1` 且落在細胞內
        （`cell_id > 0`）的 bins。
    """
    import pandas as pd

    tp = pd.read_parquet(
        str(tp_path),
        columns=["barcode", "in_tissue", "pxl_row_in_fullres", "pxl_col_in_fullres"],
    )
    tp = tp[tp["in_tissue"] == 1].copy()

    h, w = mask.shape
    row_local = (tp["pxl_row_in_fullres"].values - crop_y0).astype(np.int32).clip(0, h - 1)
    col_local = (tp["pxl_col_in_fullres"].values - crop_x0).astype(np.int32).clip(0, w - 1)
    tp["cell_id"] = mask[row_local, col_local]

    attr = tp[tp["cell_id"] > 0][["barcode", "cell_id"]].reset_index(drop=True)

    if out_path is not None:
        attr.to_parquet(str(out_path), index=False)
    return attr


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
