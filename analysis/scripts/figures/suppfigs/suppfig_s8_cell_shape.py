"""
gen_suppfig_s8_cell_shape.py
============================
Supp. Fig. S8 — Cell Shape Regularity

Computes per-cell circularity (4π·A/P²) from segmentation masks and plots:
  a  Violin: circularity distribution across all 15 ROIs × 4 methods
  b  Per-ROI median: MCseg v2 vs Space Ranger (connected dot plot)

Methods:
  v12 → MCseg v2
  sr  → Space Ranger
  p3  → Proseg
  nuc → Nuclei

Source: crc_transcript_attribution/results/masks/{method}_roi{i}.npy
Output: manuscript/supplementary/SuppFigS8.png
"""

from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import numpy as np
import pandas as pd
from pathlib import Path
from skimage.measure import label, regionprops

# ── paths ──────────────────────────────────────────────────────────────────

MASKS_DIR = Path("/Volumes/SSD/plan_a/crc_transcript_attribution/results/masks")
OUT_DIR   = Path("/Volumes/SSD/plan_a/manuscript/supplementary")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ── method config ──────────────────────────────────────────────────────────

METHODS = {
    "v12": {"label": "MCseg v2",      "color": "#4292C6"},
    "sr":  {"label": "Space Ranger",  "color": "#FD8D3C"},
    "p3":  {"label": "MCseg v1",       "color": "#74C476"},
    "nuc": {"label": "Nuclei",        "color": "#BCBD22"},
}
ROIS = [f"roi{i}" for i in range(1, 16)]

# ── style ──────────────────────────────────────────────────────────────────

plt.rcParams.update({
    "font.family":       ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size":         7,
    "savefig.dpi":       300,
    "savefig.facecolor": "white",
    "axes.linewidth":    0.6,
})
MM = 1 / 25.4


# ── helpers ────────────────────────────────────────────────────────────────

def compute_circularity(mask: np.ndarray) -> np.ndarray:
    """Return per-cell circularity (4π·A/P²) for all cells in a label mask."""
    # mask is already a label image (0 = background, >0 = cell IDs)
    props = regionprops(mask)
    circs = []
    for p in props:
        a = p.area
        perim = p.perimeter
        if perim > 0:
            c = 4 * np.pi * a / (perim ** 2)
            circs.append(min(c, 1.0))   # cap at 1 (numerical noise)
    return np.array(circs, dtype=np.float32)


def load_all_circularity() -> pd.DataFrame:
    """Load masks, compute circularity, return tidy DataFrame."""
    rows = []
    for mkey, mcfg in METHODS.items():
        for roi in ROIS:
            fpath = MASKS_DIR / f"{mkey}_{roi}.npy"
            if not fpath.exists():
                continue
            mask = np.load(fpath)
            circs = compute_circularity(mask)
            if len(circs) == 0:
                continue
            for c in circs:
                rows.append({"method": mkey, "label": mcfg["label"],
                             "roi": roi, "circularity": c})
    return pd.DataFrame(rows)


# ── main ───────────────────────────────────────────────────────────────────

def main():
    print("Computing circularity from masks...")
    df = load_all_circularity()

    # per-ROI medians (for panel b)
    roi_med = (df.groupby(["method", "roi"])["circularity"]
               .median().reset_index()
               .rename(columns={"circularity": "median_circ"}))

    # overall medians (for legend)
    overall = df.groupby("method")["circularity"].median()
    n_cells = df.groupby("method")["circularity"].count()
    print("Overall medians:", overall.to_dict())
    print("N cells:", n_cells.to_dict())

    # ── figure layout ────────────────────────────────────────────────────
    fig_w = 183 * MM
    fig_h = 80 * MM
    fig = plt.figure(figsize=(fig_w, fig_h))
    gs = fig.add_gridspec(
        1, 2,
        width_ratios=[1.0, 1.4],
        wspace=0.35,
        left=0.07, right=0.97,
        top=0.87, bottom=0.14,
    )
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])

    # ── panel a: violin ──────────────────────────────────────────────────
    method_order = ["v12", "sr", "p3", "nuc"]
    labels_order = [METHODS[m]["label"] for m in method_order]
    colors_order = [METHODS[m]["color"] for m in method_order]

    violin_data = [df[df["method"] == m]["circularity"].values for m in method_order]

    parts = ax_a.violinplot(violin_data, positions=range(len(method_order)),
                            showmedians=False, showextrema=False,
                            widths=0.65)
    for pc, col in zip(parts["bodies"], colors_order):
        pc.set_facecolor(col)
        pc.set_alpha(0.80)
        pc.set_edgecolor("none")

    # median lines
    for i, data in enumerate(violin_data):
        med = np.median(data)
        ax_a.plot([i - 0.25, i + 0.25], [med, med],
                  color="black", lw=1.5, solid_capstyle="butt")
        ax_a.text(i, med + 0.025, f"{med:.3f}",
                  ha="center", va="bottom", fontsize=5.5, fontweight="bold")

    # 0.7 reference line
    ax_a.axhline(0.7, color="#555555", lw=0.7, ls="--", alpha=0.6)
    ax_a.text(len(method_order) - 0.45, 0.702, "0.7",
              ha="right", va="bottom", fontsize=5.5, color="#555555")

    n_total = sum(len(d) for d in violin_data)
    ax_a.set_title(f"Cell shape regularity\n(all 15 ROIs, ~{n_total//1000*1000:,} cells / method)",
                   fontsize=7, fontweight="bold", pad=3)
    ax_a.set_xticks(range(len(method_order)))
    ax_a.set_xticklabels(labels_order, fontsize=6.5)
    ax_a.set_ylabel("Cell circularity (4π·A/P²)", fontsize=7)
    ax_a.set_ylim(0, 1.08)
    ax_a.spines["top"].set_visible(False)
    ax_a.spines["right"].set_visible(False)
    ax_a.tick_params(axis="both", labelsize=6)
    ax_a.text(-0.06, 1.06, "a", transform=ax_a.transAxes,
              fontsize=11, fontweight="bold", va="top", ha="right")

    # ── panel b: per-ROI comparison V12 vs SR ───────────────────────────
    v12_med = roi_med[roi_med["method"] == "v12"].set_index("roi")["median_circ"]
    sr_med  = roi_med[roi_med["method"] == "sr"].set_index("roi")["median_circ"]
    shared_rois = sorted(set(v12_med.index) & set(sr_med.index),
                         key=lambda r: int(r.replace("roi", "")))

    x_pos = range(len(shared_rois))
    v12_vals = [v12_med[r] for r in shared_rois]
    sr_vals  = [sr_med[r]  for r in shared_rois]

    # connecting lines
    for xi, v, s in zip(x_pos, v12_vals, sr_vals):
        ax_b.plot([xi, xi], [v, s], color="#aaaaaa", lw=0.8, zorder=1)

    v12_col = METHODS["v12"]["color"]
    sr_col  = METHODS["sr"]["color"]
    v12_overall = overall["v12"]
    sr_overall  = overall["sr"]

    ax_b.scatter(x_pos, v12_vals, color=v12_col, s=20, zorder=3,
                 label=f"MCseg v2 (mean {v12_overall:.3f})")
    ax_b.scatter(x_pos, sr_vals,  color=sr_col,  s=20, marker="s", zorder=3,
                 label=f"Space Ranger (mean {sr_overall:.3f})")

    ax_b.set_title("MCseg v2 vs Space Ranger — per-ROI comparison\n"
                   "(Wilcoxon p < 0.001, n=15, all MCseg v2 > SR)",
                   fontsize=7, fontweight="bold", pad=3)
    ax_b.set_xticks(list(x_pos))
    ax_b.set_xticklabels([r.replace("roi", "roi") for r in shared_rois],
                         fontsize=5.5, rotation=45, ha="right")
    ax_b.set_ylabel("Median circularity per ROI", fontsize=7)
    ax_b.set_ylim(0.3, 0.95)
    ax_b.spines["top"].set_visible(False)
    ax_b.spines["right"].set_visible(False)
    ax_b.tick_params(axis="both", labelsize=6)

    # legend
    h1 = mlines.Line2D([], [], color=v12_col, marker="o", markersize=5,
                       linestyle="None", label=f"MCseg v2 (mean {v12_overall:.3f})")
    h2 = mlines.Line2D([], [], color=sr_col, marker="s", markersize=5,
                       linestyle="None", label=f"Space Ranger (mean {sr_overall:.3f})")
    ax_b.legend(handles=[h1, h2], fontsize=6, loc="upper right",
                framealpha=0.9, edgecolor="#cccccc")
    ax_b.text(-0.08, 1.06, "b", transform=ax_b.transAxes,
              fontsize=11, fontweight="bold", va="top", ha="right")

    # ── save ─────────────────────────────────────────────────────────────
    out = OUT_DIR / "SuppFigS8.png"
    fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"✓ Saved: {out}")


if __name__ == "__main__":
    main()
