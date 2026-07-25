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
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

logger = logging.getLogger("pipeline.registration.align")


@dataclass
class AffineAlignment:
    """
    2×3 仿射變換，作用於 **(x, y)** 座標（不是 (row, col)）。

    軸序在此類程式中是頭號靜默錯誤來源 —— `phase_cross_correlation` 回傳
    `(row, col)`，而座標與仿射一律 `(x, y)`。本類別的介面**只接受 (x, y)**，
    轉換由 `shift_to_matrix` 一處負責。
    """

    matrix: list[list[float]] = field(
        default_factory=lambda: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]
    )
    source: str = "spaceranger_fullres"
    target: str = "raw_btf"
    estimated_error: float | None = None

    @classmethod
    def identity(cls) -> "AffineAlignment":
        return cls()

    @property
    def array(self) -> np.ndarray:
        return np.asarray(self.matrix, dtype=float).reshape(2, 3)

    def apply(self, coords_xy: np.ndarray) -> np.ndarray:
        """對 (N, 2) 的 `(x, y)` 座標套用變換。"""
        pts = np.asarray(coords_xy, dtype=float)
        m = self.array
        return pts @ m[:, :2].T + m[:, 2]

    def to_3x3(self) -> np.ndarray:
        """轉為 3×3 homography，供 `bin_attribution(transform=...)` 使用。"""
        return np.vstack([self.array, [0.0, 0.0, 1.0]])

    def is_identity(self, atol: float = 1e-12) -> bool:
        return bool(np.allclose(self.array, [[1, 0, 0], [0, 1, 0]], atol=atol))

    def to_dict(self) -> dict:
        return {
            "matrix": [[float(v) for v in row] for row in self.array],
            "source": self.source,
            "target": self.target,
            "estimated_error": self.estimated_error,
        }

    @classmethod
    def from_dict(cls, d: dict) -> "AffineAlignment":
        return cls(
            matrix=[[float(v) for v in row] for row in d.get("matrix", [[1, 0, 0], [0, 1, 0]])],
            source=str(d.get("source", "spaceranger_fullres")),
            target=str(d.get("target", "raw_btf")),
            estimated_error=(
                None if d.get("estimated_error") is None else float(d["estimated_error"])
            ),
        )


def shift_to_matrix(
    dy: float,
    dx: float,
    downsample: int = 1,
    estimated_error: float | None = None,
) -> AffineAlignment:
    """
    把 `estimate_shift` 的 `(dy, dx)` 轉成作用於 `(x, y)` 的平移矩陣。

    **兩個轉換都在這裡一次做完**，因為它們都不會報錯、只會靜默出錯：

    1. **軸序**：`phase_cross_correlation` 回傳 `(row, col)` = `(dy, dx)`，
       而仿射與座標一律 `(x, y)` = `(col, row)` → 必須交換。
    2. **單位**：估計值在 downsample 圖上，需乘回倍率才是 fullres px。
    """
    tx, ty = dx * downsample, dy * downsample
    return AffineAlignment(
        matrix=[[1.0, 0.0, float(tx)], [0.0, 1.0, float(ty)]],
        estimated_error=estimated_error,
    )


def estimate_shift(
    ref: np.ndarray, mov: np.ndarray, upsample_factor: int = 10
) -> tuple[float, float, float]:
    """
    以相位相關估計 `mov` 相對 `ref` 的次像素位移。

    Returns
    -------
    tuple[float, float, float]
        `(dy, dx, error)`，單位為**輸入圖的像素**。語意為
        `mov ≈ np.roll(ref, (dy, dx))` —— 即 mov 的內容相對 ref 往
        下 `dy`、右 `dx` 移動了多少。

    Notes
    -----
    兩張圖先各自 z-score 標準化。H&E 灰階與 bin 密度的量級差好幾個數量級，
    未標準化時直流分量會主導頻譜，讓相關峰失準。
    """
    from skimage.registration import phase_cross_correlation

    a, b = _zscore(ref), _zscore(mov)
    shift, error, _ = phase_cross_correlation(a, b, upsample_factor=upsample_factor)
    # skimage 回傳「要把 mov 移多少才能對上 ref」→ 取負號即 mov 相對 ref 的位移
    return float(-shift[0]), float(-shift[1]), float(error)


def _zscore(img: np.ndarray) -> np.ndarray:
    arr = np.asarray(img, dtype=np.float32)
    std = float(arr.std())
    if std == 0:
        return np.zeros_like(arr)
    return (arr - float(arr.mean())) / std


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
