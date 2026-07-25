"""
對位誤差的量化（bin 密度光柵化 + 次像素位移估計）
================================================

**定位（P0.5 之後）**：主要的座標對位由 `alignment.py` 的組合 homography 負責
（幾何正確、來源可追溯）。本模組是**驗證手段** —— 回答「JSON 套下去之後，
RNA 與影像還剩多少殘餘偏移？」20-50 px 的殘餘位移不會觸發
`bin_attribution` 的 30% 越界警告，卻足以讓 RNA 落到隔壁細胞。
"""
from __future__ import annotations

import logging
from pathlib import Path

import numpy as np

logger = logging.getLogger("pipeline.registration.align")


def render_bin_density(
    tp_path: str | Path,
    full_shape: tuple[int, int],
    downsample: int,
    transform: np.ndarray | None = None,
) -> np.ndarray:
    """
    把 Visium bin 質心光柵化成低解析密度圖。

    bin 只落在組織上（`in_tissue == 1`），因此密度圖的輪廓就是組織輪廓 ——
    正好可以和 H&E 縮圖做相位相關。

    Parameters
    ----------
    tp_path : str | Path
        `tissue_positions.parquet` 路徑。
    full_shape : tuple[int, int]
        目標影像的 `(height, width)`（fullres px）。
    downsample : int
        輸出解析度為 `(H // downsample, W // downsample)`。
    transform : np.ndarray | None
        3×3 homography，SR fullres → 影像 px（見 `alignment.compose_bin_to_image`）。
        給定時先套用再光柵化 —— 這樣估出的殘餘位移才是「JSON 套完之後還差多少」。

    Returns
    -------
    np.ndarray
        float32 密度圖；每格的值 = 落在該格的 bin 數。超出 `full_shape` 的 bin 丟棄
        （不夾邊，否則邊界會堆出假的高密度帶而污染相位相關）。
    """
    import pandas as pd

    tp = pd.read_parquet(
        str(tp_path),
        columns=["in_tissue", "pxl_row_in_fullres", "pxl_col_in_fullres"],
    )
    tp = tp[tp["in_tissue"] == 1]

    y = tp["pxl_row_in_fullres"].values.astype(float)
    x = tp["pxl_col_in_fullres"].values.astype(float)
    if transform is not None:
        from backend.src.fullslide.pipeline import _apply_homography

        x, y = _apply_homography(np.asarray(transform, dtype=float), x, y)

    h, w = full_shape
    out_h, out_w = max(1, h // downsample), max(1, w // downsample)
    ry = np.floor(y / downsample).astype(np.int64)
    rx = np.floor(x / downsample).astype(np.int64)

    keep = (ry >= 0) & (ry < out_h) & (rx >= 0) & (rx < out_w)
    dens = np.zeros((out_h, out_w), dtype=np.float32)
    np.add.at(dens, (ry[keep], rx[keep]), 1.0)
    return dens
