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
