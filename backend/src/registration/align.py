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
    x, y = load_bin_xy(tp_path, transform)
    h, w = full_shape
    return rasterize_xy(x, y, 0, 0, w, h, downsample)


def load_bin_xy(
    tp_path: str | Path, transform: np.ndarray | None = None
) -> tuple[np.ndarray, np.ndarray]:
    """讀出 `in_tissue==1` 的 bin 質心 `(x, y)`（影像 px；給定 transform 時已套用）。"""
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
    return x, y


def rasterize_xy(
    x: np.ndarray,
    y: np.ndarray,
    x0: int,
    y0: int,
    w: int,
    h: int,
    downsample: int,
) -> np.ndarray:
    """把 `(x, y)` 點集累加成 `((h//ds), (w//ds))` 密度圖，原點為 `(x0, y0)`。"""
    out_h, out_w = max(1, h // downsample), max(1, w // downsample)
    ry = np.floor((y - y0) / downsample).astype(np.int64)
    rx = np.floor((x - x0) / downsample).astype(np.int64)

    keep = (ry >= 0) & (ry < out_h) & (rx >= 0) & (rx < out_w)
    dens = np.zeros((out_h, out_w), dtype=np.float32)
    np.add.at(dens, (ry[keep], rx[keep]), 1.0)
    return dens


# ── 兩階段（coarse-to-fine）全片位移估計 ────────────────────────────────────


def tissue_gray(rgb: np.ndarray) -> np.ndarray:
    """H&E RGB → 「組織為高值」的單通道圖。

    必須反相：組織在 H&E 上是**暗**的，而 bin 密度是組織處**高**。相位相關
    對整體反相很敏感（相關峰會翻成負峰而抓不到），所以在進相關前先對齊極性。
    """
    arr = np.asarray(rgb)
    gray = arr.astype(np.float32).mean(axis=2) if arr.ndim == 3 else arr.astype(np.float32)
    return 255.0 - gray


def block_mean(img: np.ndarray, ds: int) -> np.ndarray:
    """整數倍降採樣（區塊平均）；邊緣不足一格者裁掉。"""
    if ds <= 1:
        return img.astype(np.float32)
    h, w = img.shape[:2]
    h2, w2 = (h // ds) * ds, (w // ds) * ds
    if h2 == 0 or w2 == 0:
        return np.zeros((0, 0), dtype=np.float32)
    return img[:h2, :w2].astype(np.float32).reshape(h // ds, ds, w // ds, ds).mean(axis=(1, 3))


def _pick_dense_windows(
    dens: np.ndarray, coarse_ds: int, n: int, window: int, full_shape: tuple[int, int]
) -> list[tuple[int, int]]:
    """從粗density圖挑 n 個高密度且彼此分離的窗格左上角（fullres px）。"""
    work = dens.copy()
    radius = max(1, window // coarse_ds)
    h, w = full_shape
    picks: list[tuple[int, int]] = []

    for _ in range(n):
        if work.max() <= 0:
            break
        ry, rx = np.unravel_index(int(np.argmax(work)), work.shape)
        cx, cy = (rx + 0.5) * coarse_ds, (ry + 0.5) * coarse_ds
        x0 = int(np.clip(cx - window / 2, 0, max(0, w - window)))
        y0 = int(np.clip(cy - window / 2, 0, max(0, h - window)))
        picks.append((x0, y0))
        # 抑制已選區域，強迫下一個窗格落在別處（三個窗互相獨立才有離群判斷力）
        work[
            max(0, ry - radius):ry + radius + 1,
            max(0, rx - radius):rx + radius + 1,
        ] = 0.0
    return picks


def estimate_shift_fullres(
    btf_path: str | Path,
    tp_path: str | Path,
    full_shape: tuple[int, int] | None = None,
    *,
    coarse_ds: int | None = None,
    fine_ds: int = 4,
    window: int = 1024,
    n_windows: int = 3,
    transform: np.ndarray | None = None,
) -> tuple[float, float, float]:
    """
    兩階段估計「bin 相對 H&E 影像」的位移，單位為 **fullres px**。

    ① **粗估**：在 1/`coarse_ds` 全片縮圖上跑一次相位相關，定出量級與窗格位置。
    ② **精修**：在 `n_windows` 個高 bin 密度區域各取 `window`×`window`，於
       1/`fine_ds` 重跑，取**中位數**（抗離群），`spread` = 三窗估計的最大差。

    **為何必須兩階段**：`coarse_ds=32` + `upsample_factor=10` 的理論極限是
    0.1 縮圖 px ＝ 3.2 fullres px；而傷害最大的 20-50 px 偏移在縮圖上只有
    0.6-1.6 px，SNR 太低。單階段會給出「看起來有數字但不可信」的結果 ——
    比沒有更危險。ds=32 時主導訊號是組織輪廓（適合粗對位），細節紋理已被
    抹平（不適合精修）。

    Returns
    -------
    tuple[float, float, float]
        `(dy, dx, spread)`。`spread` 為可信度指標：三窗估計愈一致愈可信；
        精修失敗時回退粗估，此時 `spread` 為 `inf`。
    """
    import tifffile

    from backend.src.roi.extractor import read_btf_crop
    from backend.src.roi.tile_server import THUMB_SCALE, _load_or_build_thumb

    coarse_ds = THUMB_SCALE if coarse_ds is None else coarse_ds
    btf_path = Path(btf_path)

    if full_shape is None:
        with tifffile.TiffFile(str(btf_path)) as tf:
            page = tf.pages[0]
            full_shape = (int(page.imagelength), int(page.imagewidth))

    x, y = load_bin_xy(tp_path, transform)
    if len(x) == 0:
        raise ValueError("tissue_positions 內沒有 in_tissue==1 的 bin，無法估計位移")

    # ── ① 粗估 ──────────────────────────────────────────────────────────
    thumb = tissue_gray(_load_or_build_thumb(btf_path, coarse_ds))
    dens = rasterize_xy(x, y, 0, 0, full_shape[1], full_shape[0], coarse_ds)
    hh = min(thumb.shape[0], dens.shape[0])
    ww = min(thumb.shape[1], dens.shape[1])
    dy_c, dx_c, _ = estimate_shift(thumb[:hh, :ww], dens[:hh, :ww])
    dy_c, dx_c = dy_c * coarse_ds, dx_c * coarse_ds
    logger.info(f"粗估位移（1/{coarse_ds}）：dy={dy_c:.1f}, dx={dx_c:.1f} fullres px")

    # ── ② 精修 ──────────────────────────────────────────────────────────
    ests: list[tuple[float, float]] = []
    for x0, y0 in _pick_dense_windows(dens[:hh, :ww], coarse_ds, n_windows, window, full_shape):
        try:
            crop, ax0, ay0 = read_btf_crop(btf_path, x0, y0, window, window)
        except (OSError, ValueError, NotImplementedError) as e:
            logger.warning(f"窗格 ({x0}, {y0}) 讀取失敗，略過：{e}")
            continue
        if crop.size == 0:
            continue

        ch, cw = crop.shape[:2]
        # 座標一律以 read_btf_crop 回傳的 actual origin 為準 —— tile 對齊可能讓
        # 實際原點與請求值不同，用請求值會整批偏移。
        ref = block_mean(tissue_gray(crop), fine_ds)
        mov = rasterize_xy(x, y, ax0, ay0, cw, ch, fine_ds)
        if ref.shape != mov.shape or mov.sum() == 0 or ref.std() == 0:
            continue

        dy_f, dx_f, _ = estimate_shift(ref, mov)
        ests.append((dy_f * fine_ds, dx_f * fine_ds))

    if not ests:
        logger.warning("精修階段無可用窗格，回退粗估（可信度未知）")
        return dy_c, dx_c, float("inf")

    arr = np.asarray(ests, dtype=float)
    dy, dx = float(np.median(arr[:, 0])), float(np.median(arr[:, 1]))
    spread = float(max(np.ptp(arr[:, 0]), np.ptp(arr[:, 1]))) if len(arr) > 1 else float("inf")
    logger.info(
        f"精修位移（{len(ests)} 窗，1/{fine_ds}）：dy={dy:.1f}, dx={dx:.1f}，spread={spread:.1f} px"
    )
    return dy, dx, spread
