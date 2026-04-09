"""
gen_suppfig_s6_crc_strip.py
============================
Supp. Fig. S6 / S6b — CRC 7-ROI benchmarking strip
4 columns per row: H&E | V12 | P3 | SR

Outputs:
  supplementary/SuppFigS6.png   — with 200px zoom inset
  supplementary/SuppFigS6b.png  — no zoom (full-field overview)
"""

from __future__ import annotations

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import cv2
from pathlib import Path

# ── paths ─────────────────────────────────────────────────────────────────────

HE_DIR   = Path("/Volumes/SSD/plan_a/crc_he_seg/results/rois")
MASK_DIR = Path("/Volumes/SSD/plan_a/crc_transcript_attribution/results/masks")
OUT_DIR  = Path("/Volumes/SSD/plan_a/manuscript/supplementary")

ROIS = ["roi1", "roi2", "roi3", "roi4", "roi5", "roi6", "roi7"]

# method colours (RGB 0-255)
COLOR_V12 = (70,  130, 220)   # blue
COLOR_P3  = (60,  170,  80)   # green
COLOR_SR  = (210,  90,   0)   # orange

# ── style ──────────────────────────────────────────────────────────────────────

plt.rcParams.update({
    "font.family": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size":   8,
    "savefig.dpi": 300,
    "savefig.facecolor": "white",
})

COL_TITLES  = ["H&E", "V12", "P3", "SR"]
COL_COLORS  = ["black", "#2196F3", "#4CAF50", "#FF5722"]

# ── helpers ───────────────────────────────────────────────────────────────────

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


def load_data(roi: str):
    he = np.load(HE_DIR / f"{roi}_he.npy")
    if he.dtype != np.uint8:
        he = (np.clip(he, 0, 1) * 255).astype(np.uint8)
    v12 = np.load(MASK_DIR / f"v12_{roi}.npy")
    p3  = np.load(MASK_DIR / f"p3_{roi}.npy")
    sr  = np.load(MASK_DIR / f"sr_{roi}.npy")
    return he, v12, p3, sr


# ── no-zoom strip (S6b) ───────────────────────────────────────────────────────

def make_no_zoom():
    n_rows = len(ROIS)
    n_cols = 4
    fig, axes = plt.subplots(n_rows, n_cols,
                             figsize=(n_cols * 3.2, n_rows * 3.2),
                             dpi=200)
    plt.subplots_adjust(hspace=0.03, wspace=0.02,
                        left=0.06, right=0.99, top=0.97, bottom=0.01)

    for i, roi in enumerate(ROIS):
        he, v12, p3, sr = load_data(roi)
        img_v12 = blend(he, mask_to_rgba_overlay(v12, COLOR_V12))
        img_p3  = blend(he, mask_to_rgba_overlay(p3,  COLOR_P3))
        img_sr  = blend(he, mask_to_rgba_overlay(sr,  COLOR_SR))
        imgs = [he, img_v12, img_p3, img_sr]

        for j, (img, ax) in enumerate(zip(imgs, axes[i])):
            ax.imshow(img, origin="upper", interpolation="antialiased")
            ax.set_xticks([]); ax.set_yticks([])

            # column title (first row only)
            if i == 0:
                ax.set_title(COL_TITLES[j], fontsize=11, fontweight="bold",
                             color=COL_COLORS[j], pad=5)
            # row label (first column only)
            if j == 0:
                ax.set_ylabel(roi.upper(), fontsize=9, fontweight="bold",
                              rotation=0, labelpad=34, va="center")

    # legend
    handles = [
        mpatches.Patch(color=COL_COLORS[1], label="V12 (this work)"),
        mpatches.Patch(color=COL_COLORS[2], label="P3 (visiumHD_pipeline_3)"),
        mpatches.Patch(color=COL_COLORS[3], label="SR (Space Ranger reference)"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=3,
               fontsize=8, framealpha=0.9, bbox_to_anchor=(0.5, 0.0))

    out = OUT_DIR / "SuppFigS6b.png"
    fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"✓ Saved: {out}")


# ── zoom strip (S6) ───────────────────────────────────────────────────────────

ZOOM_SZ = 200   # px crop size for inset
ZOOM_POS = (350, 350)  # (y, x) top-left of crop in 1000×1000 image

def add_zoom_inset(ax, img: np.ndarray, y0: int, x0: int, sz: int):
    crop = img[y0:y0+sz, x0:x0+sz]
    h_img, w_img = img.shape[:2]
    # inset axes: bottom-right corner
    ax_in = ax.inset_axes([0.60, 0.0, 0.40, 0.40])
    ax_in.imshow(crop, origin="upper", interpolation="antialiased")
    ax_in.set_xticks([]); ax_in.set_yticks([])
    # rectangle on main image
    from matplotlib.patches import Rectangle
    rect = Rectangle((x0, y0), sz, sz,
                     linewidth=1.2, edgecolor="yellow", facecolor="none", zorder=5)
    ax.add_patch(rect)


def make_zoom():
    n_rows = len(ROIS)
    n_cols = 4
    fig, axes = plt.subplots(n_rows, n_cols,
                             figsize=(n_cols * 3.2, n_rows * 3.2),
                             dpi=200)
    plt.subplots_adjust(hspace=0.03, wspace=0.02,
                        left=0.06, right=0.99, top=0.97, bottom=0.04)

    for i, roi in enumerate(ROIS):
        he, v12, p3, sr = load_data(roi)
        img_v12 = blend(he, mask_to_rgba_overlay(v12, COLOR_V12))
        img_p3  = blend(he, mask_to_rgba_overlay(p3,  COLOR_P3))
        img_sr  = blend(he, mask_to_rgba_overlay(sr,  COLOR_SR))
        imgs = [he, img_v12, img_p3, img_sr]

        for j, (img, ax) in enumerate(zip(imgs, axes[i])):
            ax.imshow(img, origin="upper", interpolation="antialiased")
            ax.set_xticks([]); ax.set_yticks([])
            add_zoom_inset(ax, img, ZOOM_POS[0], ZOOM_POS[1], ZOOM_SZ)

            if i == 0:
                ax.set_title(COL_TITLES[j], fontsize=11, fontweight="bold",
                             color=COL_COLORS[j], pad=5)
            if j == 0:
                ax.set_ylabel(roi.upper(), fontsize=9, fontweight="bold",
                              rotation=0, labelpad=34, va="center")

    handles = [
        mpatches.Patch(color=COL_COLORS[1], label="V12 (this work)"),
        mpatches.Patch(color=COL_COLORS[2], label="P3 (visiumHD_pipeline_3)"),
        mpatches.Patch(color=COL_COLORS[3], label="SR (Space Ranger reference)"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=3,
               fontsize=8, framealpha=0.9, bbox_to_anchor=(0.5, 0.0))

    out = OUT_DIR / "SuppFigS6.png"
    fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"✓ Saved: {out}")


# ── main ──────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    make_no_zoom()
    make_zoom()
