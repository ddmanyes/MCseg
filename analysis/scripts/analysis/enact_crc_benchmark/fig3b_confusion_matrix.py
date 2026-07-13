"""
fig3b_confusion_matrix.py
Fig 3(b): MCseg+lookup confusion matrix on covered cells,
with MCseg vs StarDist+WBA recall comparison on shared-covered subset.

Layout: 1×2
  Left:  Row-normalised confusion matrix (MCseg+lookup, n=11,709)
  Right: Per-class recall bar chart — MCseg vs StarDist on
         both-covered cells (n=10,275)

Output:
  manuscript/figures/fig3/fig3b.png

Usage:
  cd /Volumes/SSD/plan_a
  uv run python submission_bioinformatics/scripts/analysis/enact_crc_benchmark/fig3b_confusion_matrix.py
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd
from sklearn.metrics import confusion_matrix, f1_score

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE        = Path("/Volumes/SSD/plan_a")
LOOKUP_CSV  = BASE / "submission_bioinformatics/results/enact_crc_f1/gt_matched.csv"
GT_EVAL_CSV = (
    BASE
    / "tissue sample/ENACT_supporting_files/public_data/human_colorectal"
    / "paper_results/chunks/weighted_by_area/celltypist_results/eval"
    / "cell_annotation_eval.csv"
)
# Figure 3 panels b + c（混淆矩陣 + recall 長條，合圖中裁成兩格）
OUT_PATH = BASE / "submission_bioinformatics/figures/fig3/fig3bc.png"
OUT_PATH.parent.mkdir(parents=True, exist_ok=True)

CLASSES      = ["epithelial cells", "stromal cells", "immune cells"]
CLASS_LABELS = ["Epithelial", "Stromal", "Immune"]

# ── Load & align data ─────────────────────────────────────────────────────────
lookup_raw = pd.read_csv(LOOKUP_CSV)
gt_eval    = pd.read_csv(GT_EVAL_CSV, usecols=["cell_x", "cell_y", "gt_label", "pred_label_clean"])

assert len(lookup_raw) == len(gt_eval) == 20991, "Row count mismatch"

combined = lookup_raw.copy()
combined["sd_pred"] = gt_eval["pred_label_clean"].values
combined["sd_pred"] = combined["sd_pred"].fillna("unmatched")
combined.loc[combined["sd_pred"] == "no label", "sd_pred"] = "unmatched"

# MCseg covered + labelled (3 classes only)
mcseg_valid = combined[
    (combined["mcseg_cell_id"] > 0)
    & combined["pred_label"].isin(CLASSES)
    & combined["gt_label"].isin(CLASSES)
].copy()

# Both covered: MCseg ∩ StarDist (fair per-class recall comparison)
both_valid = mcseg_valid[mcseg_valid["sd_pred"].isin(CLASSES)].copy()

n_mcseg = len(mcseg_valid)
n_both  = len(both_valid)

y_true_m   = mcseg_valid["gt_label"].values
y_pred_m   = mcseg_valid["pred_label"].values
y_true_b   = both_valid["gt_label"].values
y_pred_m_b = both_valid["pred_label"].values
y_pred_s_b = both_valid["sd_pred"].values

micro_m   = f1_score(y_true_m, y_pred_m, labels=CLASSES, average="micro", zero_division=0)
micro_m_b = f1_score(y_true_b, y_pred_m_b, labels=CLASSES, average="micro", zero_division=0)
micro_s_b = f1_score(y_true_b, y_pred_s_b, labels=CLASSES, average="micro", zero_division=0)

print(f"MCseg covered:  {n_mcseg:,}  micro F1 = {micro_m:.3f}")
print(f"Both covered:   {n_both:,}  MCseg F1 = {micro_m_b:.3f}  StarDist F1 = {micro_s_b:.3f}")

# ── Confusion matrix (MCseg+lookup) ──────────────────────────────────────────
cm = confusion_matrix(y_true_m, y_pred_m, labels=CLASSES, normalize="true")

# ── Per-class recall on shared subset ────────────────────────────────────────
def per_class_recall(y_true: np.ndarray, y_pred: np.ndarray) -> list[float]:
    return [
        float((y_pred[y_true == c] == c).sum() / (y_true == c).sum())
        if (y_true == c).sum() > 0 else 0.0
        for c in CLASSES
    ]

recall_m = per_class_recall(y_true_b, y_pred_m_b)
recall_s = per_class_recall(y_true_b, y_pred_s_b)

# ── Figure ────────────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(11, 4.5))
fig.patch.set_facecolor("white")
gs = gridspec.GridSpec(1, 2, figure=fig, width_ratios=[1.0, 1.0], wspace=0.45)

ax_cm  = fig.add_subplot(gs[0])
ax_bar = fig.add_subplot(gs[1])

# ── Left panel: row-normalised confusion matrix ───────────────────────────────
im = ax_cm.imshow(cm, vmin=0, vmax=1, cmap="Blues")
cbar = plt.colorbar(im, ax=ax_cm, fraction=0.046, pad=0.04)
cbar.ax.tick_params(labelsize=12)

ax_cm.set_xticks(range(3))
ax_cm.set_yticks(range(3))
ax_cm.set_xticklabels(CLASS_LABELS, fontsize=14)
ax_cm.set_yticklabels(CLASS_LABELS, fontsize=14)
ax_cm.set_xlabel("Predicted", fontsize=15)
ax_cm.set_ylabel("True label", fontsize=15)

for i in range(3):
    for j in range(3):
        text_color = "white" if cm[i, j] > 0.6 else "black"
        ax_cm.text(
            j, i, f"{cm[i, j]:.2f}",
            ha="center", va="center",
            fontsize=17, color=text_color, fontweight="bold",
        )

# ── Right panel: per-class recall comparison ─────────────────────────────────
x     = np.arange(3)
width = 0.34
c_mcseg = "#1565C0"
c_sd    = "#aaaaaa"

bars_m = ax_bar.bar(
    x - width / 2, recall_m, width,
    label="MCseg",
    color=c_mcseg, alpha=0.88,
)
bars_s = ax_bar.bar(
    x + width / 2, recall_s, width,
    label="StarDist+WBA",
    color=c_sd, alpha=0.88,
)

for bar in list(bars_m) + list(bars_s):
    h = bar.get_height()
    ax_bar.text(
        bar.get_x() + bar.get_width() / 2, h + 0.012,
        f"{h:.2f}", ha="center", va="bottom", fontsize=13,
    )

ax_bar.set_xticks(x)
ax_bar.set_xticklabels(CLASS_LABELS, fontsize=14)
ax_bar.tick_params(axis="y", labelsize=12)
ax_bar.set_ylabel("Recall", fontsize=15)
ax_bar.set_ylim(0, 1.18)
ax_bar.legend(fontsize=12, loc="upper right")
ax_bar.spines[["top", "right"]].set_visible(False)

plt.tight_layout()
plt.savefig(OUT_PATH, dpi=300, bbox_inches="tight", facecolor="white")
print(f"Saved: {OUT_PATH}")
