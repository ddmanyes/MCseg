"""
gen_fig3a_crc_comparison.py
============================
Fig. 3a — CRC Representative ROI Comparison (4-column)
ROI2 & ROI4 × H&E | V12 | P3 | SR

Output: manuscript/figures/04_crc_tas/fig3a_crc_comparison.png
"""

from __future__ import annotations

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import cv2
from pathlib import Path

# ── paths ────────────────────────────────────────────────────────────────────

HE_DIR   = Path("/Volumes/SSD/plan_a/crc_he_seg/results/rois")
MASK_DIR = Path("/Volumes/SSD/plan_a/crc_transcript_attribution/results/masks")
OUT_PATH = Path("/Volumes/SSD/plan_a/manuscript/figures/04_crc_tas/fig3a_crc_comparison.png")

SELECTED_ROIS = ["roi2", "roi4"]

# method colours (RGB 0-255)
COLOR_V2  = (70,  130, 220)   # blue
COLOR_V1  = (60,  170,  80)   # green
COLOR_SR  = (210,  90,   0)   # orange

# ── style ────────────────────────────────────────────────────────────────────

plt.rcParams.update({
    "font.family": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size":   8,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "savefig.facecolor": "white",
})

# ── helpers ──────────────────────────────────────────────────────────────────

def mask_to_rgba_overlay(mask: np.ndarray, color_rgb: tuple,
                          alpha: float = 0.35,
                          boundary_alpha: float = 0.88) -> np.ndarray:
    h, w = mask.shape
    overlay = np.zeros((h, w, 4), dtype=np.float32)
    fill = mask > 0
    overlay[fill, :3] = np.array(color_rgb) / 255.0
    overlay[fill, 3]  = alpha

    m32    = mask.astype(np.int32)
    binary = fill.astype(np.uint8)
    eroded = cv2.erode(binary, np.ones((3, 3), np.uint8), iterations=1)
    outer  = (binary - eroded).astype(bool)
    r_diff = np.zeros((h, w), bool); r_diff[:, :-1] = m32[:, 1:] != m32[:, :-1]
    d_diff = np.zeros((h, w), bool); d_diff[:-1, :]  = m32[1:, :] != m32[:-1, :]
    boundary = outer | r_diff | d_diff
    overlay[boundary, :3] = 1.0
    overlay[boundary, 3]  = boundary_alpha
    return overlay


def blend(base: np.ndarray, over: np.ndarray) -> np.ndarray:
    bg = base.astype(np.float32) / 255.0
    fg, a = over[..., :3], over[..., 3:4]
    return np.clip((fg * a + bg * (1.0 - a)) * 255, 0, 255).astype(np.uint8)


# ── main ─────────────────────────────────────────────────────────────────────

def main() -> None:
    n_rows = len(SELECTED_ROIS)
    n_cols = 4   # H&E | V12 | P3 | SR

    fig, axes = plt.subplots(n_rows, n_cols,
                             figsize=(n_cols * 3.5, n_rows * 3.6),
                             dpi=200)
    plt.subplots_adjust(hspace=0.04, wspace=0.03)

    col_titles = ["H&E", "MCseg v2", "MCseg v1", "SR"]
    col_colors = ["black", "#2196F3", "#4CAF50", "#FF5722"]

    for i, roi in enumerate(SELECTED_ROIS):
        # load H&E
        he_raw = np.load(HE_DIR / f"{roi}_he.npy")
        if he_raw.dtype != np.uint8:
            he_raw = (np.clip(he_raw, 0, 1) * 255).astype(np.uint8)

        # load masks
        v12 = np.load(MASK_DIR / f"v12_{roi}.npy")
        p3  = np.load(MASK_DIR / f"p3_{roi}.npy")
        sr  = np.load(MASK_DIR / f"sr_{roi}.npy")

        # composite images
        img_v2 = blend(he_raw, mask_to_rgba_overlay(v2, COLOR_V2))
        img_v1  = blend(he_raw, mask_to_rgba_overlay(v1,  COLOR_V1))
        img_sr  = blend(he_raw, mask_to_rgba_overlay(sr,  COLOR_SR))

        imgs = [he_raw, img_v2, img_v1, img_sr]

        for j, (img, ax) in enumerate(zip(imgs, axes[i])):
            ax.imshow(img, origin="upper", interpolation="antialiased")
            ax.set_xticks([]); ax.set_yticks([])

            # column title (first row only)
            if i == 0:
                ax.set_title(col_titles[j], fontsize=11, fontweight="bold",
                             color=col_colors[j], pad=5)

            # row label (first column only)
            if j == 0:
                ax.set_ylabel(roi.upper(), fontsize=10, fontweight="bold",
                              rotation=0, labelpad=32, va="center")

            # panel label "a" — top-left of first panel
            if i == 0 and j == 0:
                ax.text(-0.18, 1.06, "a", transform=ax.transAxes,
                        fontsize=11, fontweight="bold", va="top", ha="right")

    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PATH, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"✓ Saved: {OUT_PATH}")


if __name__ == "__main__":
    main()
