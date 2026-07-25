"""
fig_suppS13_enact_summary.py
Supplementary Figure S13 — ENACT CRC Benchmark Summary

Panel A: FTC vs GT Coverage scatter (coverage efficiency)
Panel B: Per-cell-type F1 grouped bar + subset wF1 summary
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from pathlib import Path

RESULTS_DIR = Path("/Volumes/SSD/plan_a/submission_bioinformatics/results/enact_method_comparison")
OUT_DIR     = Path("/Volumes/SSD/plan_a/submission_bioinformatics/supplementary")

# ── Shared style ───────────────────────────────────────────────────────────────
method_colors = {
    "StarDist": "#aaaaaa",
    "MCseg":    "#2196F3",
    "SR":       "#FF9800",
    "NUC":      "#9C27B0",
    "ProSeg":   "#F44336",
}
method_markers = {
    "StarDist": "x",
    "MCseg":    "o",
    "SR":       "s",
    "NUC":      "^",
    "ProSeg":   "D",
}

# ── Load data ──────────────────────────────────────────────────────────────────
intrinsic = pd.read_csv(RESULTS_DIR / "intrinsic_metrics.csv")
f1_df     = pd.read_csv(RESULTS_DIR / "method_comparison_f1.csv")

intrinsic["efficiency"] = intrinsic["gt_coverage"] / intrinsic["ftc"]

method_map = {
    "StarDist": "StarDist+WBA",
    "MCseg":    "MCseg+WBA",
    "SR":       "SR+WBA",
    "NUC":      "NUC+WBA",
    "ProSeg":   "ProSeg+WBA",
}
intrinsic["method_f1_key"] = intrinsic["method"].map(method_map)
merged = intrinsic.merge(
    f1_df[["method", "subset_wf1"]],
    left_on="method_f1_key", right_on="method", suffixes=("", "_f1")
)

# Panel B — exclude StarDist (ref)
f1_cmp = f1_df[~f1_df["method"].str.contains(r"\(ref\)", regex=True)].copy()
name_map = {
    "StarDist+WBA": "StarDist",
    "MCseg+WBA":    "MCseg",
    "SR+WBA":       "SR",
    "NUC+WBA":      "NUC",
    "ProSeg+WBA":   "ProSeg",
}
f1_cmp["label"] = f1_cmp["method"].map(name_map)
methods_order = ["StarDist", "MCseg", "SR", "NUC", "ProSeg"]
f1_cmp = f1_cmp.set_index("label").loc[methods_order].reset_index()

# ── Figure layout ──────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(14, 5))
gs  = gridspec.GridSpec(1, 3, figure=fig, width_ratios=[1.1, 1.6, 0.5],
                        wspace=0.38)

ax_a  = fig.add_subplot(gs[0])
ax_b1 = fig.add_subplot(gs[1])
ax_b2 = fig.add_subplot(gs[2])

fig.text(0.01, 0.97, "A", fontsize=14, fontweight="bold", va="top")
fig.text(0.40, 0.97, "B", fontsize=14, fontweight="bold", va="top")

# ── Panel A: FTC vs GT Coverage scatter ───────────────────────────────────────
x_ref = np.linspace(0, 1.05, 100)
ax_a.plot(x_ref, x_ref, "k--", alpha=0.25, linewidth=1,
          label="efficiency = 1 (baseline)")

label_offsets = {
    "StarDist": ( 0.01, -0.05),
    "MCseg":    ( 0.01,  0.02),
    "SR":       (-0.01, -0.05),
    "NUC":      ( 0.01,  0.02),
    "ProSeg":   (-0.10, -0.05),
}
for _, row in merged.iterrows():
    m     = row["method"]
    alpha = 0.45 if m == "StarDist" else 1.0
    ax_a.scatter(row["ftc"], row["gt_coverage"],
                 color=method_colors[m], marker=method_markers[m],
                 s=110, zorder=3, alpha=alpha)
    label = f"{m}\n(self-ref)" if m == "StarDist" else \
            f"{m}\n(eff={row['efficiency']:.2f})"
    dx, dy = label_offsets.get(m, (0.01, 0.02))
    ax_a.annotate(label, (row["ftc"] + dx, row["gt_coverage"] + dy),
                  fontsize=7.5, color=method_colors[m])

ax_a.set_xlabel("FTC (Fraction of Tissue Covered)", fontsize=10)
ax_a.set_ylabel("GT Coverage (centroid match rate)", fontsize=10)
ax_a.set_title("Coverage Efficiency\n(above diagonal = efficient)", fontsize=9.5)
ax_a.set_xlim(0, 1.1)
ax_a.set_ylim(0, 1.1)
ax_a.legend(fontsize=7.5, loc="lower right")
ax_a.grid(True, alpha=0.25)

# ── Panel B1: grouped bar (3 cell types × 5 methods) ─────────────────────────
cell_types = ["Epithelial", "Stromal", "Immune"]
f1_cols    = ["epithelial_f1", "stromal_f1", "immune_f1"]

n_methods = len(methods_order)
group_w   = 0.72
bar_w     = group_w / n_methods
group_gap = 0.32
x_centers = np.arange(len(cell_types)) * (group_w + group_gap)

for i, row in enumerate(f1_cmp.itertuples()):
    m     = row.label
    x_pos = x_centers + (i - n_methods / 2 + 0.5) * bar_w
    vals  = [getattr(row, col) for col in f1_cols]
    bars  = ax_b1.bar(x_pos, vals, width=bar_w * 0.88,
                      color=method_colors[m], label=m,
                      edgecolor="white", linewidth=0.5)
    for bar, val in zip(bars, vals):
        if val > 0.22:
            ax_b1.text(bar.get_x() + bar.get_width() / 2,
                       bar.get_height() + 0.01,
                       f"{val:.2f}", ha="center", va="bottom",
                       fontsize=6, color=method_colors[m])

ax_b1.set_xticks(x_centers)
ax_b1.set_xticklabels(cell_types, fontsize=10)
ax_b1.set_ylabel("F1 Score", fontsize=10)
ax_b1.set_ylim(0, 1.0)
ax_b1.set_title("Per-cell-type F1\n(+WBA methods)", fontsize=9.5)
ax_b1.legend(title="Method", fontsize=7.5, title_fontsize=8,
             loc="upper right", framealpha=0.85)
ax_b1.grid(True, alpha=0.22, axis="y")
ax_b1.axhline(0.5, color="grey", linestyle=":", linewidth=0.8, alpha=0.5)
ax_b1.annotate(
    "Immune F1 < 0.37\n(all methods)",
    xy=(x_centers[2], 0.36), xytext=(x_centers[2] + 0.42, 0.52),
    fontsize=7.5, color="#555555",
    arrowprops=dict(arrowstyle="->", color="#888888", lw=0.8),
)

# ── Panel B2: subset wF1 summary bar ─────────────────────────────────────────
methods_rev    = methods_order[::-1]
wf1_rev        = f1_cmp["subset_wf1"].values[::-1]
bar_colors_rev = [method_colors[m] for m in methods_rev]

bars2 = ax_b2.barh(methods_rev, wf1_rev,
                   color=bar_colors_rev, height=0.5, edgecolor="white")
for bar, val, m in zip(bars2, wf1_rev, methods_rev):
    ax_b2.text(val + 0.005, bar.get_y() + bar.get_height() / 2,
               f"{val:.3f}", va="center", fontsize=8.5,
               color=method_colors[m])

ax_b2.set_xlim(0.55, 0.82)
ax_b2.set_xlabel("Subset Weighted F1", fontsize=10)
ax_b2.set_title("Overall\nsubset wF1", fontsize=9.5)
ax_b2.axvline(0.7, color="grey", linestyle=":", linewidth=0.8, alpha=0.5)
ax_b2.grid(True, alpha=0.22, axis="x")

# ── Save ───────────────────────────────────────────────────────────────────────
out_path = OUT_DIR / "SuppFigS13.png"
plt.savefig(out_path, dpi=200, bbox_inches="tight")
plt.close()
print(f"Saved: {out_path}")
