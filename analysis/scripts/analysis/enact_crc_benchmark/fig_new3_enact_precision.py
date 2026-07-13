"""
fig_new3_enact_precision.py
===========================
Figure 3 panel a — MCseg 在 ENACT CRC 資料集對 GT centroid 的覆蓋

左：全局 overview（8× downsample）——MCseg mask 輪廓 + GT centroid
    綠點 = matched，紅點 = unmatched（黃框為右圖放大範圍）
右：局部放大（腺體–間質界面，900×900 px）——同色系

僅保留組織影像；WBA confusion matrix 由 fig3b_confusion_matrix.py 產出
（Figure 3 panel b/c），此處不再重複繪製。matched/unmatched 比例會印到
stdout，供 figure caption 使用。

輸出：
  submission_bioinformatics/figures/fig3/fig3a.png（300 DPI）

執行：
  cd /Volumes/SSD/plan_a
  uv run python submission_bioinformatics/scripts/analysis/enact_crc_benchmark/fig_new3_enact_precision.py
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import tifffile
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
from skimage.segmentation import find_boundaries
from scipy.ndimage import binary_dilation
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "figures"))
from _panel_style import add_scale_bar

# ── Paths ──────────────────────────────────────────────────────────────────────
BASE     = Path("/Volumes/SSD/plan_a/submission_bioinformatics")
RES_DIR  = BASE / "results" / "enact_crc_f1"
WBA_DIR  = BASE / "results" / "enact_crc_f1_wba"
# Figure 3 panel a（組織 overview + zoom）。檔名對應合圖字母。
OUT_PATH = BASE / "figures" / "fig3" / "fig3a.png"

CROP_X0, CROP_Y0 = 5154, 4635   # ENACT local coord → mask pixel offset

# ── Colours ────────────────────────────────────────────────────────────────────
COL_MATCHED   = "#44CC77"
COL_UNMATCHED = "#FF4444"
CLASSES       = ["epithelial cells", "stromal cells", "immune cells"]
CLASS_LABELS  = ["Epithelial", "Stromal", "Immune"]

# ── Load data ──────────────────────────────────────────────────────────────────
print("Loading H&E crop…")
he   = tifffile.imread(str(RES_DIR / "he_crop.tif"))          # (H, W, 3) uint8
H, W = he.shape[:2]

print("Loading MCseg mask…")
mask = np.load(str(RES_DIR / "mcseg_mask.npy")).astype(np.int32)

print("Loading GT matched (lookup)…")
gt = pd.read_csv(RES_DIR / "gt_matched.csv")
gt = gt.dropna(subset=["cell_x", "cell_y"]).copy()
gt["col_local"] = (gt["cell_x"] - CROP_X0).astype(int)
gt["row_local"] = (gt["cell_y"] - CROP_Y0).astype(int)
valid = (
    (gt["col_local"] >= 0) & (gt["col_local"] < W) &
    (gt["row_local"] >= 0) & (gt["row_local"] < H)
)
gt       = gt[valid].copy()
matched  = gt[gt["mcseg_cell_id"] > 0]
unmatched= gt[gt["mcseg_cell_id"] == 0]
print(f"  GT in bounds: {len(gt):,}  matched: {len(matched):,}  unmatched: {len(unmatched):,}")

# ── Overview: 8× downsample + boundary overlay ────────────────────────────────
print("Building overview (8× downsample)…")
STEP  = 8
he_ov = he[::STEP, ::STEP].copy()

print("  computing mask boundaries…")
boundary = find_boundaries(mask, mode="outer")
boundary_ds = binary_dilation(boundary, iterations=1)[::STEP, ::STEP]
he_ov[boundary_ds] = [255, 255, 200]   # pale-yellow boundary (visible on pink H&E)

# ── Zoom region: upper-centre of crop (gland–stroma interface) ─────────────────
ZOOM_R0   = 1800
ZOOM_C0   = 8100
ZOOM_SIZE = 900   # pixels full-res ≈ 450 µm at 0.5 µm/px

he_zoom    = he[ZOOM_R0:ZOOM_R0+ZOOM_SIZE, ZOOM_C0:ZOOM_C0+ZOOM_SIZE].copy()
mask_zoom  = mask[ZOOM_R0:ZOOM_R0+ZOOM_SIZE, ZOOM_C0:ZOOM_C0+ZOOM_SIZE]
bound_zoom = binary_dilation(find_boundaries(mask_zoom, mode="outer"), iterations=1)
he_zoom[bound_zoom] = [255, 255, 180]   # pale-yellow

def in_zoom(df):
    return df[
        (df["row_local"] >= ZOOM_R0) & (df["row_local"] < ZOOM_R0 + ZOOM_SIZE) &
        (df["col_local"] >= ZOOM_C0) & (df["col_local"] < ZOOM_C0 + ZOOM_SIZE)
    ].copy()

zm_m = in_zoom(matched)
zm_u = in_zoom(unmatched)
print(f"  zoom matched: {len(zm_m):,}  unmatched: {len(zm_u):,}")

# ── Figure layout ──────────────────────────────────────────────────────────────
# Figure 3 panel a 只放組織影像：左 = overview、右 = zoom。
# 混淆矩陣改由 fig3bc.py 統一產出（panel b），此處不再重複繪製。
h_ov, w_ov = he_ov.shape[:2]
fig = plt.figure(figsize=(14, 7), dpi=300)
gs  = GridSpec(
    1, 2, figure=fig,
    left=0.03, right=0.97, top=0.93, bottom=0.03, wspace=0.10,
    width_ratios=[w_ov / h_ov, 1.0],   # 兩張圖等高（zoom 為正方形）
)

ax_ov   = fig.add_subplot(gs[0, 0])   # overview
ax_zoom = fig.add_subplot(gs[0, 1])   # zoom

# ── A1: Overview ───────────────────────────────────────────────────────────────
ax_ov.imshow(he_ov, origin="upper", interpolation="nearest")

um_c_ds = unmatched["col_local"].values / STEP
um_r_ds = unmatched["row_local"].values / STEP
m_c_ds  = matched["col_local"].values   / STEP
m_r_ds  = matched["row_local"].values   / STEP

ax_ov.scatter(um_c_ds, um_r_ds, s=1.5, c=COL_UNMATCHED, alpha=0.65,
              linewidths=0, rasterized=True, zorder=2)
ax_ov.scatter(m_c_ds,  m_r_ds,  s=1.5, c=COL_MATCHED,   alpha=0.65,
              linewidths=0, rasterized=True, zorder=3)

# zoom rectangle
zr0 = ZOOM_R0 / STEP;  zc0 = ZOOM_C0 / STEP
zr1 = (ZOOM_R0 + ZOOM_SIZE) / STEP;  zc1 = (ZOOM_C0 + ZOOM_SIZE) / STEP
ax_ov.add_patch(mpatches.Rectangle(
    (zc0, zr0), zc1 - zc0, zr1 - zr0,
    linewidth=1.5, edgecolor="yellow", facecolor="none", zorder=4
))

# scale bar: 1 mm，白底板 + 深色尺線（原本白字白線壓在組織上看不清）
# downsampled 像素邊長 = STEP × 0.5 µm = 4 µm
add_scale_bar(ax_ov, img_w_px=he_ov.shape[1], fontsize=15,
              scale_um=1000, um_per_px=STEP * 0.5, label="1 mm")

ax_ov.axis("off")

# ── A2: Zoom ───────────────────────────────────────────────────────────────────
ax_zoom.imshow(he_zoom, origin="upper", interpolation="nearest")
ax_zoom.scatter(
    zm_m["col_local"] - ZOOM_C0, zm_m["row_local"] - ZOOM_R0,
    s=9, c=COL_MATCHED,   alpha=0.75, linewidths=0, rasterized=True, zorder=2
)
ax_zoom.scatter(
    zm_u["col_local"] - ZOOM_C0, zm_u["row_local"] - ZOOM_R0,
    s=12, c=COL_UNMATCHED, alpha=0.90, linewidths=0.4,
    edgecolors="white", rasterized=True, zorder=3
)

# scale bar: 100 µm @ 0.5 µm/px（full-res zoom）
add_scale_bar(ax_zoom, img_w_px=ZOOM_SIZE, fontsize=15,
              scale_um=100, um_per_px=0.5)
ax_zoom.set_title("Zoom (gland–stroma interface)", fontsize=16, pad=5)
ax_zoom.axis("off")

# ── 圖說用數字（不畫進圖裡，改由 figure caption 敘述）─────────────────────────
match_rate = len(matched) / len(gt) * 100
type_stats = []
for lbl, short in zip(CLASSES, CLASS_LABELS):
    n  = (gt["gt_label"] == lbl).sum()
    nu = ((gt["gt_label"] == lbl) & (gt["mcseg_cell_id"] == 0)).sum()
    type_stats.append(f"{short} {nu/n*100:.0f}%")

print("\n── Caption 用數字 ─────────────────────────────")
print(f"  GT matched (綠, covered by MCseg):   {match_rate:.1f}%  n = {len(matched):,}")
print(f"  GT unmatched (紅, centroid in gap): {100-match_rate:.1f}%  n = {len(unmatched):,}")
print(f"  Unmatched rate by cell type: {' · '.join(type_stats)}")

# ── Save ───────────────────────────────────────────────────────────────────────
OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(str(OUT_PATH), dpi=300, bbox_inches="tight")
print(f"Saved → {OUT_PATH}")
