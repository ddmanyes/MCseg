"""
fig_celltype_f1_breakdown.py
Per-cell-type F1 breakdown: epithelial / stromal / immune for each method.

Two panels:
  Left  — grouped bar chart (3 cell types × 5 methods)
  Right — subset_wF1 summary bar (same methods, for reference)
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

RESULTS_DIR = Path("/Volumes/SSD/plan_a/submission_bioinformatics/results/enact_method_comparison")

# ── Load ───────────────────────────────────────────────────────────────────────
df = pd.read_csv(RESULTS_DIR / "method_comparison_f1.csv")

# Exclude StarDist (ref) row — different match rate makes it incomparable
df = df[~df["method"].str.contains(r"\(ref\)", regex=True)].copy()

name_map = {
    "StarDist+WBA": "StarDist",
    "MCseg+WBA":    "MCseg",
    "SR+WBA":       "SR",
    "NUC+WBA":      "NUC",
    "ProSeg+WBA":   "ProSeg",
}
df["label"] = df["method"].map(name_map).fillna(df["method"])

methods  = ["StarDist", "MCseg", "SR", "NUC", "ProSeg"]
df = df.set_index("label").loc[methods].reset_index()

cell_types = ["Epithelial", "Stromal", "Immune"]
f1_cols    = ["epithelial_f1", "stromal_f1", "immune_f1"]

method_colors = {
    "StarDist": "#aaaaaa",
    "MCseg":    "#2196F3",
    "SR":       "#FF9800",
    "NUC":      "#9C27B0",
    "ProSeg":   "#F44336",
}

# ── Figure ─────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(13, 5),
                         gridspec_kw={"width_ratios": [3, 1]})
fig.suptitle("ENACT CRC Benchmark — Per-Cell-Type Annotation F1 (+WBA methods)",
             fontsize=12, fontweight="bold", y=1.01)

# ── Panel A: grouped bar chart ─────────────────────────────────────────────────
ax = axes[0]

n_methods = len(methods)
group_w   = 0.75
bar_w     = group_w / n_methods
group_gap = 0.35
x_centers = np.arange(len(cell_types)) * (group_w + group_gap)

for i, row in enumerate(df.itertuples()):
    m     = row.label
    x_pos = x_centers + (i - n_methods / 2 + 0.5) * bar_w
    vals  = [getattr(row, col) for col in f1_cols]
    bars  = ax.bar(x_pos, vals, width=bar_w * 0.88,
                   color=method_colors[m], label=m,
                   edgecolor="white", linewidth=0.5)
    for bar, val in zip(bars, vals):
        if val > 0.22:
            ax.text(bar.get_x() + bar.get_width() / 2,
                    bar.get_height() + 0.008,
                    f"{val:.2f}", ha="center", va="bottom",
                    fontsize=6.5, color=method_colors[m])

ax.set_xticks(x_centers)
ax.set_xticklabels(cell_types, fontsize=12)
ax.set_ylabel("F1 Score", fontsize=11)
ax.set_ylim(0, 1.0)
ax.set_title("Per-cell-type F1 by segmentation method", fontsize=10)
ax.legend(title="Method", fontsize=8.5, title_fontsize=9,
          loc="upper right", framealpha=0.85)
ax.grid(True, alpha=0.25, axis="y")
ax.axhline(0.5, color="black", linestyle=":", linewidth=0.8, alpha=0.35)
ax.text(x_centers[-1] + group_w / 2 + 0.04, 0.505,
        "F1 = 0.5", fontsize=7.5, color="grey", va="bottom")

# annotation: immune is universally low
ax.annotate(
    "Immune F1 < 0.37\nacross all methods",
    xy=(x_centers[2], 0.37), xytext=(x_centers[2] + 0.48, 0.54),
    fontsize=8.5, color="#555555",
    arrowprops=dict(arrowstyle="->", color="#888888", lw=0.9),
)

# ── Panel B: subset_wF1 horizontal bar ────────────────────────────────────────
ax2 = axes[1]

methods_rev    = methods[::-1]
wf1_vals_rev   = df["subset_wf1"].values[::-1]
bar_colors_rev = [method_colors[m] for m in methods_rev]

bars2 = ax2.barh(methods_rev, wf1_vals_rev,
                 color=bar_colors_rev, height=0.5, edgecolor="white")
for bar, val, m in zip(bars2, wf1_vals_rev, methods_rev):
    ax2.text(val + 0.005, bar.get_y() + bar.get_height() / 2,
             f"{val:.3f}", va="center", fontsize=9, color=method_colors[m])

ax2.set_xlim(0.55, 0.82)
ax2.set_xlabel("Subset Weighted F1", fontsize=10)
ax2.set_title("Overall\nsubset wF1", fontsize=10)
ax2.axvline(0.7, color="black", linestyle=":", linewidth=0.8, alpha=0.35)
ax2.grid(True, alpha=0.25, axis="x")

plt.tight_layout()
out_path = RESULTS_DIR / "fig_celltype_f1_breakdown.png"
plt.savefig(out_path, dpi=200, bbox_inches="tight")
plt.close()
print(f"Saved: {out_path}")
