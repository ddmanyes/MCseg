"""
fig_suppS11_with_pvals.py
Supplementary Figure S11 (updated) — ENACT CRC benchmark: 3-condition comparison
with bootstrap 95% CI and McNemar's test between MCseg+Lookup and MCseg+WBA.

Conditions:
  1. StarDist+WBA  — Lotfollahi et al. reference (aggregate only, no CI)
  2. MCseg+Lookup  — per-cell data available → bootstrap CI
  3. MCseg+WBA     — per-cell data available → bootstrap CI

Statistical tests:
  - Bootstrap 95% CI for micro F1 and weighted F1 (MCseg conditions)
  - McNemar's test on cells covered by BOTH conditions (correct/incorrect)
  - One-sided bootstrap permutation p-value for WBA > Lookup

Output: supplementary/SuppFigS11.png

Usage:
  cd /Volumes/SSD/plan_a
  uv run python submission_bioinformatics/scripts/analysis/enact_crc_benchmark/fig_suppS11_with_pvals.py
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
from sklearn.metrics import f1_score, confusion_matrix
import seaborn as sns
from scipy.stats import chi2 as chi2_dist
from matplotlib.gridspec import GridSpec

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE       = Path("/Volumes/SSD/plan_a/submission_bioinformatics")
LOOKUP_CSV = BASE / "results" / "enact_crc_f1"     / "gt_matched.csv"
WBA_CSV    = BASE / "results" / "enact_crc_f1_wba" / "gt_matched_wba.csv"
OUT_DIR    = BASE / "supplementary"

# StarDist reference values (Lotfollahi et al. 2025, full-denominator n=20,991)
STARDIST_MICRO    = 0.708
STARDIST_WEIGHTED = 0.758

CLASSES       = ["epithelial cells", "stromal cells", "immune cells"]
CLASS_LABELS  = ["Epithelial", "Stromal", "Immune"]
N_BOOTSTRAP   = 5000
RNG_SEED      = 42

COLORS = {
    "StarDist": "#aaaaaa",
    "Lookup":   "#64B5F6",
    "WBA":      "#1565C0",
}


# ── Statistical helpers ───────────────────────────────────────────────────────

def bootstrap_f1_ci(
    y_true: np.ndarray,
    y_pred: np.ndarray,
    average: str = "weighted",
    n: int = N_BOOTSTRAP,
    seed: int = RNG_SEED,
) -> tuple[float, float, float]:
    """Return (point, lower_95, upper_95) via percentile bootstrap."""
    rng = np.random.default_rng(seed)
    size = len(y_true)
    scores = []
    for _ in range(n):
        idx = rng.integers(0, size, size=size)
        scores.append(
            f1_score(y_true[idx], y_pred[idx], average=average,
                     labels=CLASSES, zero_division=0)
        )
    scores = np.array(scores)
    point = f1_score(y_true, y_pred, average=average,
                     labels=CLASSES, zero_division=0)
    return point, float(np.percentile(scores, 2.5)), float(np.percentile(scores, 97.5))


def mcnemar_test(
    y_true: np.ndarray,
    y_pred_a: np.ndarray,
    y_pred_b: np.ndarray,
) -> tuple[float, float]:
    """McNemar's test (with continuity correction): b wins over a.
    Returns (chi2, p_value)."""
    correct_a = y_true == y_pred_a
    correct_b = y_true == y_pred_b
    n_10 = int(np.sum(correct_a & ~correct_b))
    n_01 = int(np.sum(~correct_a & correct_b))
    if (n_10 + n_01) == 0:
        return 0.0, 1.0
    chi2_val = (abs(n_10 - n_01) - 1) ** 2 / (n_10 + n_01)
    p = chi2_dist.sf(chi2_val, df=1)
    return float(chi2_val), float(p)


def pvalue_stars(p: float) -> str:
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"


# ── Load data ─────────────────────────────────────────────────────────────────
# pred_label has values: 'epithelial cells', 'stromal cells', 'immune cells',
# 'unmatched' (centroid outside any MCseg mask), 'other' (low-confidence label).
# Covered subset = rows where pred_label is one of the three GT classes.
# This gives Lookup=11,709 and WBA=11,948, matching f1_summary.csv n_valid.

lookup_raw = pd.read_csv(LOOKUP_CSV)
wba_raw    = pd.read_csv(WBA_CSV)

lookup_valid = lookup_raw[lookup_raw["pred_label"].isin(CLASSES)].reset_index(drop=True)
wba_valid    = wba_raw[wba_raw["pred_label"].isin(CLASSES)].reset_index(drop=True)

y_true_l = lookup_valid["gt_label"].values
y_pred_l = lookup_valid["pred_label"].values
y_true_w = wba_valid["gt_label"].values
y_pred_w = wba_valid["pred_label"].values

print(f"Lookup covered subset: {len(lookup_valid):,}  (expected 11,709)")
print(f"WBA    covered subset: {len(wba_valid):,}  (expected 11,948)")

# ── Bootstrap CI ──────────────────────────────────────────────────────────────

print("Computing bootstrap CIs (n=5,000)…")
micro_l,    micro_lo_l,    micro_hi_l    = bootstrap_f1_ci(y_true_l, y_pred_l, "micro")
micro_w,    micro_lo_w,    micro_hi_w    = bootstrap_f1_ci(y_true_w, y_pred_w, "micro")
weighted_l, weighted_lo_l, weighted_hi_l = bootstrap_f1_ci(y_true_l, y_pred_l, "weighted")
weighted_w, weighted_lo_w, weighted_hi_w = bootstrap_f1_ci(y_true_w, y_pred_w, "weighted")

print(f"  Lookup  micro={micro_l:.4f} [{micro_lo_l:.4f}, {micro_hi_l:.4f}]")
print(f"  WBA     micro={micro_w:.4f} [{micro_lo_w:.4f}, {micro_hi_w:.4f}]")
print(f"  Lookup  wF1  ={weighted_l:.4f} [{weighted_lo_l:.4f}, {weighted_hi_l:.4f}]")
print(f"  WBA     wF1  ={weighted_w:.4f} [{weighted_lo_w:.4f}, {weighted_hi_w:.4f}]")

# ── McNemar's test on cells covered by BOTH conditions ───────────────────────
# Both files have the same 20,991 GT rows in identical order.
# "Covered by both" = pred_label in CLASSES in both files (same row index).

both_covered_mask = (
    lookup_raw["pred_label"].isin(CLASSES) & wba_raw["pred_label"].isin(CLASSES)
)
shared_lookup_pred = lookup_raw.loc[both_covered_mask, "pred_label"].values
shared_wba_pred    = wba_raw.loc[both_covered_mask,    "pred_label"].values
shared_gt          = lookup_raw.loc[both_covered_mask, "gt_label"].values

print(f"Cells covered by both conditions: {both_covered_mask.sum():,}")

chi2_val, mcnemar_p = mcnemar_test(shared_gt, shared_lookup_pred, shared_wba_pred)
stars = pvalue_stars(mcnemar_p)
print(f"McNemar χ²={chi2_val:.2f}  p={mcnemar_p:.4f}  {stars}")

# ── Confusion matrix (MCseg+WBA covered subset) ───────────────────────────────

cm = confusion_matrix(y_true_w, y_pred_w, labels=CLASSES, normalize="true")

# ── Figure ────────────────────────────────────────────────────────────────────

fig = plt.figure(figsize=(14, 5.8))
fig.patch.set_facecolor("white")

gs = GridSpec(1, 2, figure=fig, width_ratios=[1.6, 1.0], wspace=0.40)
ax_bar = fig.add_subplot(gs[0])
ax_cm  = fig.add_subplot(gs[1])

# ── Panel A: grouped bar chart ────────────────────────────────────────────────

cond_labels = [
    "StarDist+WBA\n(ref · n=20,991)",
    "MCseg+Lookup\n(subset · n=11,709)",
    "MCseg+WBA\n(subset · n=11,948)",
]
x = np.arange(3)
w = 0.30

data = [
    # (micro,    micro_lo,    micro_hi,    weighted,    weighted_lo,    weighted_hi,    color_micro, alpha)
    (STARDIST_MICRO, None, None, STARDIST_WEIGHTED, None, None, COLORS["StarDist"], 0.55),
    (micro_l, micro_lo_l, micro_hi_l, weighted_l, weighted_lo_l, weighted_hi_l, COLORS["Lookup"], 1.0),
    (micro_w, micro_lo_w, micro_hi_w, weighted_w, weighted_lo_w, weighted_hi_w, COLORS["WBA"],    1.0),
]

for i, (mv, mlo, mhi, wv, wlo, whi, col, alpha) in enumerate(data):
    yerr_m = None if mlo is None else [[mv - mlo], [mhi - mv]]
    yerr_w = None if wlo is None else [[wv - wlo], [whi - wv]]

    ax_bar.bar(x[i] - w / 2, mv, w, color=col, alpha=alpha,
               yerr=yerr_m, capsize=4,
               error_kw={"elinewidth": 1.5, "ecolor": "#444"},
               label="Micro F1" if i == 0 else "_")
    top_m = mhi + 0.015 if mhi else mv + 0.015
    ax_bar.text(x[i] - w / 2, top_m, f"{mv:.3f}",
                ha="center", va="bottom", fontsize=8.5, color=col,
                fontweight="bold" if i > 0 else "normal")

    ax_bar.bar(x[i] + w / 2, wv, w,
               color=col if i > 0 else "#cccccc", alpha=alpha * 0.70,
               yerr=yerr_w, capsize=4,
               error_kw={"elinewidth": 1.5, "ecolor": "#444"},
               label="Weighted F1" if i == 0 else "_")
    top_w = whi + 0.015 if whi else wv + 0.015
    ax_bar.text(x[i] + w / 2, top_w, f"{wv:.3f}",
                ha="center", va="bottom", fontsize=8.5,
                color=col if i > 0 else "#999",
                fontweight="bold" if i > 0 else "normal")

# Significance bracket between Lookup and WBA
y_brk = 0.905
ax_bar.plot([x[1] + w / 2, x[1] + w / 2, x[2] - w / 2, x[2] - w / 2],
            [y_brk - 0.008, y_brk, y_brk, y_brk - 0.008],
            color="black", lw=1.2)
ax_bar.text((x[1] + x[2]) / 2, y_brk + 0.004, stars,
            ha="center", va="bottom", fontsize=12, color="black")
ax_bar.text((x[1] + x[2]) / 2, y_brk - 0.030,
            f"McNemar p={mcnemar_p:.4f}",
            ha="center", va="top", fontsize=7.5, color="#444")

# StarDist denominator note
ax_bar.text(x[0], 0.36,
            "★ Full-denominator reference\n(Lotfollahi et al. 2025)",
            ha="center", va="bottom", fontsize=7, color="#999", style="italic")

ax_bar.set_xticks(x)
ax_bar.set_xticklabels(cond_labels, fontsize=9.5)
ax_bar.set_ylabel("F1 Score", fontsize=11)
ax_bar.set_ylim(0.33, 0.96)
ax_bar.set_title(
    "ENACT CRC Benchmark — 3-Condition F1\n"
    "(error bars: 95% bootstrap CI, n=5,000 resamples)",
    fontsize=10.5,
)
ax_bar.legend(loc="lower right", fontsize=9, framealpha=0.85)
ax_bar.grid(True, alpha=0.22, axis="y")
ax_bar.axhline(0.708, color="#aaa", linestyle="--", linewidth=0.8, alpha=0.5)
ax_bar.text(x[-1] + 0.47, 0.712, "StarDist ref",
            fontsize=6.5, color="#aaa", va="bottom")

ax_bar.text(
    0.02, 0.02,
    "⚠  Denominators differ: StarDist = all 20,991 GT cells;\n"
    "    MCseg = covered subset only. Not directly comparable.",
    transform=ax_bar.transAxes,
    fontsize=7, color="#666", va="bottom",
    bbox=dict(boxstyle="round,pad=0.3", facecolor="#fffff8",
              edgecolor="#ccc", alpha=0.9),
)

# ── Panel B: confusion matrix ─────────────────────────────────────────────────

sns.heatmap(
    cm, annot=True, fmt=".2f", cmap="Blues",
    xticklabels=CLASS_LABELS, yticklabels=CLASS_LABELS,
    ax=ax_cm, linewidths=0.5, linecolor="#ddd",
    vmin=0, vmax=1,
    cbar_kws={"label": "Proportion (row-normalised)", "shrink": 0.85},
)
ax_cm.set_xlabel("Predicted", fontsize=10)
ax_cm.set_ylabel("True (GT)", fontsize=10)
ax_cm.set_title("Confusion Matrix\nMCseg+WBA (covered subset n=11,948)", fontsize=10)
ax_cm.set_xticklabels(CLASS_LABELS, rotation=30, ha="right", fontsize=9)
ax_cm.set_yticklabels(CLASS_LABELS, rotation=0, fontsize=9)

# ── Figure-level caption line ─────────────────────────────────────────────────

fig.text(
    0.5, -0.03,
    f"Statistical note: McNemar χ²={chi2_val:.2f}, p={mcnemar_p:.4f} ({stars}), "
    f"shared cells n={int(both_covered_mask.sum()):,}. "
    f"Bootstrap 95% CI: Lookup micro [{micro_lo_l:.3f}–{micro_hi_l:.3f}], "
    f"WBA micro [{micro_lo_w:.3f}–{micro_hi_w:.3f}].",
    ha="center", fontsize=7.5, color="#555",
)

plt.tight_layout()

out_path = OUT_DIR / "SuppFigS11.png"
plt.savefig(out_path, dpi=200, bbox_inches="tight")
plt.close()
print(f"\nSaved → {out_path}")
