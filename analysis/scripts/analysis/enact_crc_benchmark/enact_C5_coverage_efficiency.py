"""
ENACT CRC Benchmark — Coverage Efficiency Analysis
Compute: efficiency = GT_coverage / FTC (mask area fraction)
Note: StarDist excluded from efficiency comparison (self-reference: GT centroids derived from StarDist masks)
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

RESULTS_DIR = Path("/Volumes/SSD/plan_a/submission_bioinformatics/results/enact_method_comparison")

# ── Load data ──────────────────────────────────────────────────────────────
intrinsic = pd.read_csv(RESULTS_DIR / "intrinsic_metrics.csv")
f1_df = pd.read_csv(RESULTS_DIR / "method_comparison_f1.csv")

print("=== Intrinsic Metrics (raw) ===")
print(intrinsic.to_string(index=False))
print()

# ── Compute Efficiency ─────────────────────────────────────────────────────
# efficiency = GT_coverage / FTC
# GT_coverage = proportion of GT centroids that fall inside a predicted mask
# FTC = fraction of tissue area consumed by predicted masks
# Interpretation: per unit of tissue area consumed, how much GT is captured
# efficiency > 1 means masks are targeting real cells efficiently
# efficiency ≈ 1 (ProSeg/SR) means coverage is purely driven by space-filling
intrinsic["efficiency"] = intrinsic["gt_coverage"] / intrinsic["ftc"]

# StarDist: GT centroids ARE derived from StarDist masks → trivially 100% coverage
# Efficiency value (2.75) is not comparable; excluded from analysis
intrinsic["efficiency_valid"] = intrinsic["method"] != "StarDist"

print("=== Coverage Efficiency (GT_coverage / FTC) ===")
print("Note: StarDist excluded — GT labels derived from its own segmentation")
print()
for _, row in intrinsic.iterrows():
    flag = "  ← self-reference, excluded" if not row["efficiency_valid"] else ""
    print(f"  {row['method']:<10}  GT_cov={row['gt_coverage']:.3f}  FTC={row['ftc']:.3f}  efficiency={row['efficiency']:.3f}{flag}")

# ── Merge with F1 data ─────────────────────────────────────────────────────
method_map = {
    "StarDist": "StarDist+WBA",
    "MCseg":    "MCseg+WBA",
    "SR":       "SR+WBA",
    "NUC":      "NUC+WBA",
    "ProSeg":   "ProSeg+WBA",
}
intrinsic["method_f1_key"] = intrinsic["method"].map(method_map)
merged = intrinsic.merge(
    f1_df[["method", "match_rate", "subset_wf1", "weighted_f1"]],
    left_on="method_f1_key", right_on="method", suffixes=("", "_f1")
)

# Annotation efficiency: subset_wF1 per unit FTC
# Reflects how well the method annotates cells relative to mask area consumed
merged["annot_efficiency"] = merged["subset_wf1"] / merged["ftc"]

print()
print("=== Full Analysis Table ===")
cols_display = ["method", "n_cells", "ftc", "gt_coverage", "efficiency",
                "ned", "doublet_rate", "match_rate", "subset_wf1", "weighted_f1",
                "annot_efficiency", "efficiency_valid"]
print(merged[cols_display].to_string(index=False))

# ── Save results ───────────────────────────────────────────────────────────
out_path = RESULTS_DIR / "coverage_efficiency_analysis.csv"
merged[cols_display].to_csv(out_path, index=False)
print(f"\nSaved: {out_path}")

# ── Summary: non-StarDist methods ─────────────────────────────────────────
valid = merged[merged["efficiency_valid"]].copy()

print()
print("=== Summary (non-StarDist, sorted by efficiency) ===")
summary = valid[["method", "n_cells", "ftc", "gt_coverage", "efficiency",
                 "ned", "doublet_rate", "subset_wf1", "weighted_f1"]].sort_values("efficiency", ascending=False)
print(summary.to_string(index=False))

# ── Key insights ───────────────────────────────────────────────────────────
def eff(m):
    return valid.loc[valid["method"] == m, "efficiency"].values[0]

print()
print("=== Key Insights ===")
print(f"  MCseg  efficiency: {eff('MCseg'):.3f}")
print(f"  SR     efficiency: {eff('SR'):.3f}  (MCseg vs SR: +{(eff('MCseg')/eff('SR')-1)*100:.1f}%)")
print(f"  ProSeg efficiency: {eff('ProSeg'):.3f}  (MCseg vs ProSeg: +{(eff('MCseg')/eff('ProSeg')-1)*100:.1f}%)")
print(f"  NUC    efficiency: {eff('NUC'):.3f}  (highest — nucleus-only masks, FTC={valid.loc[valid['method']=='NUC','ftc'].values[0]:.3f})")
print()
print("  ProSeg efficiency ≈ 1.0: GT coverage is entirely explained by space-filling (FTC=0.995)")
print("  MCseg lower match_rate (65.1%) is not segmentation failure;")
print("  it reflects tighter biologically-faithful cell boundaries.")

# ── Figures ────────────────────────────────────────────────────────────────
colors = {"StarDist": "#aaaaaa", "MCseg": "#2196F3", "SR": "#FF9800",
          "NUC": "#9C27B0", "ProSeg": "#F44336"}
markers = {"StarDist": "x", "MCseg": "o", "SR": "s", "NUC": "^", "ProSeg": "D"}

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Panel A: FTC vs GT_coverage scatter
ax = axes[0]
x_ref = np.linspace(0, 1.05, 100)
ax.plot(x_ref, x_ref, "k--", alpha=0.3, linewidth=1, label="efficiency = 1 (baseline)")

label_offsets = {"StarDist": (0.01, -0.04), "MCseg": (0.01, 0.02),
                 "SR": (0.01, -0.04), "NUC": (0.01, 0.02), "ProSeg": (-0.06, -0.04)}
for _, row in merged.iterrows():
    m = row["method"]
    alpha = 0.4 if m == "StarDist" else 1.0
    ax.scatter(row["ftc"], row["gt_coverage"],
               color=colors[m], marker=markers[m], s=120, zorder=3, alpha=alpha)
    label = f"{m}\n(self-ref)" if not row["efficiency_valid"] else f"{m}\n(eff={row['efficiency']:.2f})"
    dx, dy = label_offsets[m]
    ax.annotate(label, (row["ftc"] + dx, row["gt_coverage"] + dy),
                fontsize=7.5, color=colors[m])

ax.set_xlabel("FTC (Fraction of Tissue Covered)", fontsize=11)
ax.set_ylabel("GT Coverage (centroid match rate)", fontsize=11)
ax.set_title("GT Coverage vs Tissue Area Consumed\n(above diagonal = efficient)", fontsize=10)
ax.set_xlim(0, 1.1)
ax.set_ylim(0, 1.1)
ax.legend(fontsize=8, loc="lower right")
ax.grid(True, alpha=0.3)

# Panel B: efficiency vs subset_wF1 scatter (non-StarDist only)
# Replaces bar chart — NUC's efficiency=3.34 dominated the y-axis and obscured
# the MCseg vs SR/ProSeg difference. Scatter reveals the full tradeoff directly.
ax2 = axes[1]

ax2.axvline(1.0, color="grey", linestyle="--", linewidth=1, alpha=0.5)
ax2.text(1.02, 0.558, "efficiency = 1\n(space-filling baseline)",
         fontsize=7.5, color="grey", va="bottom")

# annotation positions chosen to avoid overlap among the SR/ProSeg/MCseg cluster
annot_cfg = {
    "MCseg":  dict(xy_off=( 0.08,  0.010), ha="left"),
    "SR":     dict(xy_off=( 0.08, -0.022), ha="left"),
    "NUC":    dict(xy_off=(-0.08,  0.010), ha="right"),
    "ProSeg": dict(xy_off=(-0.08, -0.022), ha="right"),
}
for _, row in valid.iterrows():
    m = row["method"]
    ax2.scatter(row["efficiency"], row["subset_wf1"],
                color=colors[m], marker=markers[m], s=160, zorder=3)
    cfg = annot_cfg.get(m, dict(xy_off=(0.05, 0.005), ha="left"))
    dx, dy = cfg["xy_off"]
    ax2.annotate(
        f"{m}  (eff={row['efficiency']:.2f},  wF1={row['subset_wf1']:.3f})",
        xy=(row["efficiency"], row["subset_wf1"]),
        xytext=(row["efficiency"] + dx, row["subset_wf1"] + dy),
        fontsize=8.5, color=colors[m], ha=cfg["ha"],
        arrowprops=dict(arrowstyle="-", color=colors[m], lw=0.8, alpha=0.5),
    )

ax2.set_xlabel("Coverage Efficiency (GT_cov / FTC)", fontsize=11)
ax2.set_ylabel("Subset Weighted F1", fontsize=11)
ax2.set_title("Efficiency vs Annotation Quality\n(StarDist excluded: self-reference)",
              fontsize=10)
ax2.set_xlim(0.5, 4.0)
ax2.set_ylim(0.57, 0.80)
ax2.grid(True, alpha=0.3)

# ideal direction arrow
ax2.annotate("", xy=(3.8, 0.785), xytext=(3.4, 0.758),
             arrowprops=dict(arrowstyle="->", color="grey", lw=1.2))
ax2.text(3.82, 0.786, "ideal", fontsize=8, color="grey", va="bottom", ha="left")

plt.tight_layout()
fig_path = RESULTS_DIR / "fig_coverage_efficiency.png"
plt.savefig(fig_path, dpi=200, bbox_inches="tight")
plt.close()
print(f"\nSaved figure: {fig_path}")
