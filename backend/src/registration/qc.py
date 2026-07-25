"""
對位 QC 疊圖
============

把 H&E 影像與 Visium bin 質心疊在同一張 patch 上輸出 PNG —— 對位有沒有問題，
肉眼一看便知，不必先相信任何估計數值。

判讀方式：bin 散點應**均勻落在組織上**。若散點整體偏離組織一小段距離，就是
殘餘位移（可用 `align.estimate_shift_fullres` 量化）；若散點壓縮或拉伸，
則是縮放／變換錯誤，該回頭檢查對位 JSON（`alignment.py`）。
"""
from __future__ import annotations

import logging
from pathlib import Path

import numpy as np

logger = logging.getLogger("pipeline.registration.qc")

QC_SUBDIR = Path("qc") / "alignment"


def render_overlay_patches(
    btf_path: str | Path,
    tp_path: str | Path,
    out_dir: str | Path,
    n: int = 3,
    size: int = 512,
    seed: int = 0,
    full_shape: tuple[int, int] | None = None,
    transform: np.ndarray | None = None,
    dpi: int = 300,
) -> list[Path]:
    """
    產生 n 張「H&E ＋ bin 質心散點」疊圖 PNG。

    Parameters
    ----------
    n, size : int
        取樣 patch 數與邊長（fullres px）。
    seed : int
        固定亂數種子 —— QC 圖必須可重現，否則兩次執行看到不同區域無從比較。
    transform : np.ndarray | None
        3×3 homography（`alignment.compose_bin_to_image`）。給定時散點畫的是
        **校正後**的位置，也就是實際用於計數的座標。

    Returns
    -------
    list[Path]
        產出的 PNG 路徑（`{out_dir}/patch_{i}.png`）。
    """
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import tifffile

    from backend.src.registration.align import load_bin_xy
    from backend.src.roi.extractor import read_btf_crop

    btf_path = Path(btf_path)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if full_shape is None:
        with tifffile.TiffFile(str(btf_path)) as tf:
            page = tf.pages[0]
            full_shape = (int(page.imagelength), int(page.imagewidth))

    x, y = load_bin_xy(tp_path, transform)
    if len(x) == 0:
        raise ValueError("tissue_positions 內沒有 in_tissue==1 的 bin，無法產生 QC 疊圖")

    rng = np.random.default_rng(seed)
    h, w = full_shape
    # 以 bin 質心為中心取樣 —— 隨機取全片座標多半落在空白區，看不出對位好壞
    idx = rng.choice(len(x), size=min(n, len(x)), replace=False)

    written: list[Path] = []
    for i, k in enumerate(idx):
        x0 = int(np.clip(x[k] - size / 2, 0, max(0, w - size)))
        y0 = int(np.clip(y[k] - size / 2, 0, max(0, h - size)))
        try:
            crop, ax0, ay0 = read_btf_crop(btf_path, x0, y0, size, size)
        except (OSError, ValueError, NotImplementedError) as e:
            logger.warning(f"QC patch ({x0}, {y0}) 讀取失敗，略過：{e}")
            continue
        if crop.size == 0:
            continue

        ch, cw = crop.shape[:2]
        # 散點座標一律以 read_btf_crop 回傳的 actual origin 為基準：tile 對齊
        # 可能讓實際原點不等於請求值，用請求值會讓整批散點偏移。
        sel = (x >= ax0) & (x < ax0 + cw) & (y >= ay0) & (y < ay0 + ch)

        fig, ax = plt.subplots(figsize=(6, 6), dpi=dpi)
        ax.imshow(crop)
        ax.scatter(
            x[sel] - ax0, y[sel] - ay0,
            s=1.0, c="#00e5ff", alpha=0.45, linewidths=0,
        )
        ax.set_title(
            f"patch {i}  fullres ({ax0}, {ay0})  {int(sel.sum()):,} bins", fontsize=8
        )
        ax.set_xlim(0, cw)
        ax.set_ylim(ch, 0)
        ax.axis("off")

        path = out_dir / f"patch_{i}.png"
        fig.savefig(str(path), dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        written.append(path)

    return written
