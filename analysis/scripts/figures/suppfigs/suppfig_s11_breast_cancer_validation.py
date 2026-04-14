"""
gen_suppfig_s11_roi10_at1at2.py
================================
Supp. Fig. S11 — ROI10: AT1/AT2 Pneumocyte Discrimination

Assembles four panels comparing MCseg v1 vs MCseg v2 for normal lung:
  a  MCseg v1 cell-type map (cellpose_dilate)
  b  MCseg v2 cell-type map (v12)
  c  Key gene markers on MCseg v2 cells (SFTPC / AT2; AGER / AT1)
  d  Quantitative comparison bar chart: AT2 & AT1 detection rates (v1 vs v2)

Verified numbers from AnnData (roi10_cellpose_dilate.h5ad / roi10_v12.h5ad):
  MCseg v1: AT2 (any SFTPC/SFTPB/SFTPA1/SFTPA2 > 0) = 38.9%, AT1 (any) = 5.1%
            SFTPC alone = 29.3%, AGER alone = 1.2%; n_cells = 948
  MCseg v2: AT2 (any) = 45.0%, AT1 (any) = 5.7%
            SFTPC alone = 34.1%, AGER alone = 1.6%; n_cells = 1085

Source images (archive, verified from real data):
  roi10_panel_A_winner_cp3.png  — H&E + MCseg v1 cell map (2823×935)
  roi10_panel_A_winner_v12.png  — H&E + MCseg v2 cell map (2800×1000)
  roi10_panel_B_markers_v12.png — 6 gene marker maps on MCseg v2 (3000×1800)

Output: manuscript/supplementary/SuppFigS11.png
"""

from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from pathlib import Path
from PIL import Image

# ── paths ──────────────────────────────────────────────────────────────────

ARCHIVE = Path("/Volumes/SSD/plan_a/manuscript/figures/_archive/root")
OUT_DIR = Path("/Volumes/SSD/plan_a/manuscript/supplementary")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ── verified detection rates (from AnnData, any-marker-positive) ──────────

DATA = {
    "MCseg v1": {
        "AT2 (any)": 38.9,
        "AT1 (any)": 5.1,
        "SFTPC":     29.3,
        "AGER":       1.2,
        "n_cells":   948,
    },
    "MCseg v2": {
        "AT2 (any)": 45.0,
        "AT1 (any)": 5.7,
        "SFTPC":     34.1,
        "AGER":       1.6,
        "n_cells":   1085,
    },
}

COLOR_V1 = "#6baed6"   # mid-blue  (MCseg v1)
COLOR_V2 = "#2171b5"   # dark-blue (MCseg v2)
COLOR_AT2 = "#2166AC"  # AT2 blue (matches cell map legend)
COLOR_AT1 = "#B2182B"  # AT1 red

# ── style ──────────────────────────────────────────────────────────────────

plt.rcParams.update({
    "font.family":    ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size":      7,
    "savefig.dpi":    300,
    "savefig.facecolor": "white",
})

MM_TO_IN = 1 / 25.4


def load_crop(filename: str, top: int = 0, bot: int | None = None,
              left: int = 0, right: int | None = None) -> np.ndarray:
    img = np.array(Image.open(ARCHIVE / filename).convert("RGB"))
    h, w = img.shape[:2]
    bot   = bot   if bot   is not None else h
    right = right if right is not None else w
    return img[top:bot, left:right]


def add_panel_label(ax, letter, x=-0.04, y=1.04):
    ax.text(x, y, letter, transform=ax.transAxes,
            fontsize=11, fontweight="bold", va="top", ha="right")


def plot_detection_bars(ax):
    """Panel d: grouped bar chart MCseg v1 vs v2 for AT2 and AT1 detection."""
    # Two groups: AT2 | AT1
    # Each group: v1 bar, v2 bar
    groups   = ["AT2 Pneumocyte\n(any marker)", "AT1 Pneumocyte\n(any marker)"]
    v1_vals  = [DATA["MCseg v1"]["AT2 (any)"], DATA["MCseg v1"]["AT1 (any)"]]
    v2_vals  = [DATA["MCseg v2"]["AT2 (any)"], DATA["MCseg v2"]["AT1 (any)"]]

    x     = np.array([0, 1.2])
    width = 0.45
    bars1 = ax.bar(x - width / 2, v1_vals, width, color=COLOR_V1,
                   alpha=0.85, label="MCseg v1", edgecolor="white", linewidth=0.5)
    bars2 = ax.bar(x + width / 2, v2_vals, width, color=COLOR_V2,
                   alpha=0.85, label="MCseg v2", edgecolor="white", linewidth=0.5)

    # value labels
    for bar in list(bars1) + list(bars2):
        h = bar.get_height()
        ax.text(bar.get_x() + bar.get_width() / 2, h + 0.6,
                f"{h:.1f}%", ha="center", va="bottom",
                fontsize=6.5, fontweight="bold")

    # Δ annotations (placed just above the taller bar)
    for xi, dv1, dv2 in zip(x, v1_vals, v2_vals):
        delta = dv2 - dv1
        ypos  = max(dv1, dv2) + 2.0
        ax.text(xi, ypos, f"+{delta:.1f}pp",
                ha="center", va="bottom", fontsize=6.5,
                color="#C0392B", fontweight="bold")

    ax.set_xticks(x)
    ax.set_xticklabels(groups, fontsize=7)
    ax.set_ylabel("% Cells Positive\n(v1: n=948, v2: n=1085)", fontsize=6.5)
    ax.set_ylim(0, 58)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", alpha=0.25, lw=0.5)

    # legend
    handles = [
        mpatches.Patch(color=COLOR_V1, label="MCseg v1\n(Cellpose + radial dilation)"),
        mpatches.Patch(color=COLOR_V2, label="MCseg v2\n(Ensemble + Voronoi)"),
    ]
    ax.legend(handles=handles, fontsize=6, loc="upper right",
              framealpha=0.9, edgecolor="#ccc")

    # physics limit annotation for AT1
    ax.axhline(7, color="#B2182B", lw=0.7, ls=":", alpha=0.6)
    ax.text(1.15, 7.0, "~5–6%\n(physics limit)",
            ha="center", va="bottom", fontsize=5.5, color="#B2182B", style="italic")


def main():
    # ── load and crop panels ────────────────────────────────────────────────

    # cp3: (2823, 935), image content starts at row 205
    # Left half = H&E, right half = MCseg v1 cell map
    cp3_full  = load_crop("roi10_panel_A_winner_cp3.png", top=205, bot=897)
    w_cp3     = cp3_full.shape[1]
    panel_a   = cp3_full[:, w_cp3 // 2:]     # MCseg v1 cell map (right half)

    # v12: (2800, 1000), content [158:851]
    # Left half = H&E, right half = MCseg v2 cell map
    v12_full  = load_crop("roi10_panel_A_winner_v12.png", top=158, bot=851)
    w_v12     = v12_full.shape[1]
    # Use v12's H&E for context (shared ROI), plus v12 cell map
    panel_b_he  = v12_full[:, :w_v12 // 2]   # H&E (same for both methods)
    panel_b_map = v12_full[:, w_v12 // 2:]   # MCseg v2 cell map

    # Markers: (3000, 1800), content [206:1603]
    # 2×3 grid: [SFTPC | SFTPB | SPP1 ] / [PECAM1 | AGER | RTKN2]
    markers_full = load_crop("roi10_panel_B_markers_v12.png", top=206, bot=1603)
    h_m, w_m = markers_full.shape[:2]
    col_w = w_m // 3
    row_h = h_m // 2
    panel_c_sftpc = markers_full[:row_h, :col_w]              # SFTPC (AT2)
    panel_c_ager  = markers_full[row_h:, col_w : col_w * 2]   # AGER  (AT1)

    # ── figure layout ───────────────────────────────────────────────────────
    # 2 rows × 2 cols; use aspect='auto' to fill axes cells cleanly
    fig_w = 183 * MM_TO_IN
    fig_h = 140 * MM_TO_IN

    fig = plt.figure(figsize=(fig_w, fig_h))
    gs  = fig.add_gridspec(
        2, 2,
        height_ratios=[1.0, 1.0],
        hspace=0.12, wspace=0.04,
        left=0.03, right=0.99,
        top=0.96, bottom=0.05,
    )
    ax_a  = fig.add_subplot(gs[0, 0])
    ax_b  = fig.add_subplot(gs[0, 1])
    ax_c  = fig.add_subplot(gs[1, 0])
    ax_d  = fig.add_subplot(gs[1, 1])

    # Panel a: MCseg v1 cell map
    ax_a.imshow(panel_a, aspect="auto")
    ax_a.axis("off")
    add_panel_label(ax_a, "a")
    ax_a.set_title("MCseg v1 — Cell-type Map\n(948 cells, Cellpose + radial dilation)",
                   fontsize=7, fontweight="bold", pad=2)

    # Panel b: MCseg v2 cell map with H&E context
    combo_b = np.concatenate([panel_b_he, panel_b_map], axis=1)
    ax_b.imshow(combo_b, aspect="auto")
    ax_b.axis("off")
    add_panel_label(ax_b, "b")
    ax_b.set_title("MCseg v2 — H&E + Cell-type Map\n(1085 cells, Ensemble + Voronoi)",
                   fontsize=7, fontweight="bold", pad=2)
    for xi, label in zip([0.25, 0.75], ["H&E", "Cell-type map"]):
        ax_b.text(xi, -0.03, label, transform=ax_b.transAxes,
                  fontsize=6, ha="center", va="top", style="italic")

    # Panel c: SFTPC (AT2) + AGER (AT1) gene maps side-by-side
    min_h = min(panel_c_sftpc.shape[0], panel_c_ager.shape[0])
    combo_c = np.concatenate([panel_c_sftpc[:min_h], panel_c_ager[:min_h]], axis=1)
    ax_c.imshow(combo_c, aspect="auto")
    ax_c.axis("off")
    add_panel_label(ax_c, "c")
    ax_c.set_title("Gene Markers on MCseg v2 Cells", fontsize=7, fontweight="bold", pad=2)
    for xi, label, color in [
        (0.25, "SFTPC (AT2 Marker)", COLOR_AT2),
        (0.75, "AGER (AT1 — Physics Limit)", COLOR_AT1),
    ]:
        ax_c.text(xi, -0.03, label, transform=ax_c.transAxes,
                  fontsize=6, ha="center", va="top",
                  color=color, fontweight="bold", style="italic")

    # Panel d: bar chart
    add_panel_label(ax_d, "d")
    ax_d.set_title("Cell-type Detection Rate: MCseg v1 vs v2",
                   fontsize=7, fontweight="bold", pad=2)
    plot_detection_bars(ax_d)

    # ── save ────────────────────────────────────────────────────────────────
    out = OUT_DIR / "SuppFigS11.png"
    fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"✓ Saved: {out}")

    # Print data summary for verification
    print("\nVerified detection rates (from AnnData):")
    for method, d in DATA.items():
        print(f"  {method} (n={d['n_cells']}): "
              f"AT2={d['AT2 (any)']:.1f}%, AT1={d['AT1 (any)']:.1f}%, "
              f"SFTPC={d['SFTPC']:.1f}%, AGER={d['AGER']:.1f}%")


if __name__ == "__main__":
    main()
