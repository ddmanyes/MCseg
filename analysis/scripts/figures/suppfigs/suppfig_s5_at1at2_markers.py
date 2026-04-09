"""
gen_suppfig_s5_roi10_at1at2.py
================================
Supp. Fig. S5 — ROI10: AT2 Pneumocyte Detection (MCseg v1 vs v2)

Streamlined 2-panel layout (1 row × 2 cols):
  a  MCseg v1 cell-type map (ROI10, Normal Lung)
  b  AT2 detection-rate bar chart (MCseg v1 vs v2)

Note:
  - MCseg v2 cell-type map is shown in Fig. 2d (not duplicated here).
  - AT1 platform comparison (Xenium vs MCseg geometric vs RNA) is shown
    in Fig. 2e–f; AT1 is intentionally excluded from this panel to avoid
    inconsistency with the different computation method used there.

Verified data from AnnData (roi10_cellpose_dilate.h5ad / roi10_v12.h5ad):
  MCseg v1 (n=948):  AT2=38.9%
  MCseg v2 (n=1085): AT2=45.0%

Output: manuscript/supplementary/SuppFigS5.png
"""

from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from PIL import Image

# ── paths ──────────────────────────────────────────────────────────────────

ARCHIVE = Path("/Volumes/SSD/plan_a/manuscript/figures/_archive/root")
OUT_DIR = Path("/Volumes/SSD/plan_a/manuscript/supplementary")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ── verified detection rates ────────────────────────────────────────────────

DATA = {
    "MCseg v1": {"AT2 (any)": 38.9, "n_cells": 948},
    "MCseg v2": {"AT2 (any)": 45.0, "n_cells": 1085},
}

COLOR_V1  = "#6baed6"
COLOR_V2  = "#2171b5"

# ── style ───────────────────────────────────────────────────────────────────

plt.rcParams.update({
    "font.family":       ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size":         7,
    "savefig.dpi":       300,
    "savefig.facecolor": "white",
    "axes.linewidth":    0.6,
})

MM = 1 / 25.4   # mm → inches


# ── helpers ─────────────────────────────────────────────────────────────────

def load_crop(filename: str,
              top: int = 0, bot: int | None = None,
              left: int = 0, right: int | None = None) -> np.ndarray:
    img = np.array(Image.open(ARCHIVE / filename).convert("RGB"))
    h, w = img.shape[:2]
    return img[top : (bot or h), left : (right or w)]


def panel_label(ax, letter: str, x: float = -0.06, y: float = 1.06):
    ax.text(x, y, letter, transform=ax.transAxes,
            fontsize=11, fontweight="bold", va="top", ha="right",
            clip_on=False)


def draw_scale_bar(ax, img_w: int, um_per_px: float = 0.2738,
                   bar_um: float = 50, frac_from_right: float = 0.06,
                   frac_from_bot: float = 0.06, color: str = "white"):
    """Draw a scale bar on an image axis.  img_w in pixels."""
    bar_px = bar_um / um_per_px
    x1 = img_w * (1 - frac_from_right) - bar_px
    x2 = img_w * (1 - frac_from_right)
    y  = img_w * frac_from_bot          # approximate y using width for square-ish crops
    ax.plot([x1, x2], [y, y], lw=2, color=color, solid_capstyle="butt",
            transform=ax.transData)
    ax.text((x1 + x2) / 2, y - img_w * 0.025, f"{bar_um} µm",
            ha="center", va="bottom", fontsize=5.5, color=color,
            transform=ax.transData)


def plot_detection_bars(ax):
    """Bar chart: AT2 detection rate for MCseg v1 vs v2."""
    v1_at2 = DATA["MCseg v1"]["AT2 (any)"]
    v2_at2 = DATA["MCseg v2"]["AT2 (any)"]

    x     = np.array([0.0, 0.6])
    width = 0.40

    bars1 = ax.bar(x[0], v1_at2, width,
                   color=COLOR_V1, alpha=0.90,
                   edgecolor="white", linewidth=0.5,
                   label="MCseg v1")
    bars2 = ax.bar(x[1], v2_at2, width,
                   color=COLOR_V2, alpha=0.90,
                   edgecolor="white", linewidth=0.5,
                   label="MCseg v2")

    # value labels on bars
    for bar in [bars1[0], bars2[0]]:
        h = bar.get_height()
        ax.text(bar.get_x() + bar.get_width() / 2, h + 0.5,
                f"{h:.1f}%", ha="center", va="bottom",
                fontsize=6, fontweight="bold", color="#333333")

    # Δ annotation
    delta = v2_at2 - v1_at2
    ypos  = v2_at2 + 3.0
    ax.annotate("", xy=(x[1], v2_at2 + 1.0), xytext=(x[0], v1_at2 + 1.0),
                arrowprops=dict(arrowstyle="-", color="#555555", lw=0.8))
    ax.text((x[0] + x[1]) / 2, ypos, f"+{delta:.1f} pp",
            ha="center", va="bottom", fontsize=6.5,
            color=COLOR_V2, fontweight="bold")

    ax.set_xticks(x)
    ax.set_xticklabels(["MCseg v1", "MCseg v2"], fontsize=7)
    ax.set_ylabel("AT2 Pneumocytes Positive (%)", fontsize=7)
    ax.set_ylim(0, 56)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis="both", labelsize=6.5)
    ax.grid(axis="y", alpha=0.20, lw=0.5, zorder=0)

    # n labels below x-axis ticks
    for xi, key in zip(x, ["MCseg v1", "MCseg v2"]):
        ax.text(xi, -4.5, f"n={DATA[key]['n_cells']}", ha="center",
                va="top", fontsize=5.5, color="#555555")


# ── main ────────────────────────────────────────────────────────────────────

def main():

    # ── load & crop source images ────────────────────────────────────────

    # cp3  (2823×935): left = H&E, right = MCseg v1 cell map
    cp3 = load_crop("roi10_panel_A_winner_cp3.png", top=205, bot=897)
    panel_a = cp3[:, cp3.shape[1] // 2:]           # v1 cell map only  (1412×692)

    # ── figure: 183mm wide × 90mm tall ──────────────────────────────────
    # Single row: a (left, wider) | b (right, bar chart)
    # Panel b = bar chart; MCseg v2 cell map is in Fig. 2d (not duplicated here)
    fig_w = 183 * MM
    fig_h = 90 * MM

    fig = plt.figure(figsize=(fig_w, fig_h))

    # explicit width_ratios: image gets 1.6x the chart width
    gs = fig.add_gridspec(
        1, 2,
        width_ratios=[1.6, 1.0],
        wspace=0.35,
        left=0.04, right=0.97,
        top=0.87, bottom=0.10,
    )

    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])

    # ── panel a: MCseg v1 cell map ───────────────────────────────────────
    ax_a.imshow(panel_a, aspect="auto")
    ax_a.axis("off")
    panel_label(ax_a, "a")
    ax_a.set_title("MCseg v1 — Cell-type Map  (ROI10, Normal Lung)\n"
                   "948 cells · nuclei + radial dilation",
                   fontsize=7, fontweight="bold", pad=3)

    # ── panel b: detection-rate bar chart ────────────────────────────────
    panel_label(ax_b, "b")
    ax_b.set_title("AT2 Pneumocyte Detection Rate:\nMCseg v1 vs MCseg v2",
                   fontsize=7, fontweight="bold", pad=3)
    plot_detection_bars(ax_b)

    # ── save ─────────────────────────────────────────────────────────────
    out = OUT_DIR / "SuppFigS5.png"
    fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"✓ Saved: {out}")

    print("\nVerified detection rates:")
    for method, d in DATA.items():
        print(f"  {method} (n={d['n_cells']}): AT2={d['AT2 (any)']:.1f}%")


if __name__ == "__main__":
    main()
