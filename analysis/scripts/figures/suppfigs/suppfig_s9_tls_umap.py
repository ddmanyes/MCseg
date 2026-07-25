"""
gen_suppfig_s10_roi9_dust.py
============================
Supp. Fig. S10 — ROI9: Pigmented Alveolar Macrophages (Dust Cells)

Source: roi9_quenching_sidebyside.png (2×2 archive image, 2723×2423)
  Grid: horizontal split rows ~1200/1311, vertical split cols ~1300/1431
  Sub-panels extracted (titles cropped out):
    TL rows 101-1199, cols 20-1299   → (a) H&E
    TR rows 101-1199, cols 1431-2703 → (b) Xenium DAPI (quenching)
    BL rows 1311-2402, cols 20-1299  → (c) Xenium segmentation on H&E
    BR rows 1311-2402, cols 1431-2703→ (d) MCseg v2 segmentation on H&E

Corrected narrative:
  Both methods segment cells; key difference is TRANSCRIPT characterization —
  smFISH transcripts are quenched by pigment (visible in panel b),
  while Visium HD captures are unaffected (SPP1+ shown in Fig. 2).

Cell counts (verified):
  Xenium: 9,694  |  MCseg v2: 8,523

Output: manuscript/supplementary/SuppFigS10.png
"""

from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from PIL import Image

ARCHIVE = Path("/Volumes/SSD/plan_a/manuscript/figures/_archive/root")
OUT_DIR = Path("/Volumes/SSD/plan_a/manuscript/supplementary")
OUT_DIR.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({
    "font.family":    ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size":      7,
    "savefig.dpi":    300,
    "savefig.facecolor": "white",
})
MM_TO_IN = 1 / 25.4

def load_crop(path, top=0, bot=None, left=0, right=None):
    img = np.array(Image.open(path).convert("RGB"))
    h, w = img.shape[:2]
    return img[top:(bot or h), left:(right or w)]

def add_panel_label(ax, letter, x=-0.04, y=1.04):
    ax.text(x, y, letter, transform=ax.transAxes,
            fontsize=11, fontweight="bold", va="top", ha="right")


def main():
    src = ARCHIVE / "roi9_quenching_sidebyside.png"

    # Extract four sub-panels (crop embedded titles)
    panel_a = load_crop(src, top=101, bot=1199, left=20,   right=1299)  # H&E
    panel_b = load_crop(src, top=101, bot=1199, left=1431, right=2703)  # Xenium DAPI
    panel_c = load_crop(src, top=1311,           left=20,   right=1299)  # Xenium seg
    panel_d = load_crop(src, top=1311,           left=1431, right=2703)  # MCseg v2 seg

    fig_w = 183 * MM_TO_IN
    fig_h = 140 * MM_TO_IN

    fig, axes = plt.subplots(
        2, 2, figsize=(fig_w, fig_h),
        gridspec_kw=dict(hspace=0.10, wspace=0.03,
                         left=0.01, right=0.99, top=0.96, bottom=0.01)
    )

    # Panel a: H&E
    axes[0, 0].imshow(panel_a, aspect="auto")
    axes[0, 0].axis("off")
    add_panel_label(axes[0, 0], "a")
    axes[0, 0].set_title("H&E — Pigmented Alveolar Macrophages (ROI9)",
                         fontsize=7, fontweight="bold", pad=2)

    # Panel b: Xenium DAPI quenching
    axes[0, 1].imshow(panel_b, aspect="auto")
    axes[0, 1].axis("off")
    add_panel_label(axes[0, 1], "b")
    axes[0, 1].set_title("Xenium DAPI — Fluorescence Quenching\n"
                         "(pigment absorbs excitation light; signal void in dust cell zones)",
                         fontsize=7, fontweight="bold", pad=2)

    # Panel c: Xenium segmentation
    axes[1, 0].imshow(panel_c, aspect="auto")
    axes[1, 0].axis("off")
    add_panel_label(axes[1, 0], "c")
    axes[1, 0].set_title("Xenium Segmentation on H&E\n"
                         "(9,694 cells; DAPI-based — cells detected, transcripts quenched)",
                         fontsize=7, fontweight="bold", pad=2)

    # Panel d: MCseg v2 segmentation
    axes[1, 1].imshow(panel_d, aspect="auto")
    axes[1, 1].axis("off")
    add_panel_label(axes[1, 1], "d")
    axes[1, 1].set_title("MCseg v2 Segmentation on H&E\n"
                         "(8,523 cells; H&E morphology — enables Visium HD SPP1⁺ characterization)",
                         fontsize=7, fontweight="bold", pad=2)

    out = OUT_DIR / "SuppFigS10.png"
    fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"✓ Saved: {out}")


if __name__ == "__main__":
    main()
