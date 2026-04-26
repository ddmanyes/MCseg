"""
視覺化：GT centroids 疊在 MCseg mask 上
- 灰底 = MCseg 分割細胞（每個細胞隨機顏色）
- 綠點 = GT centroid 已匹配（mcseg_cell_id > 0）
- 紅點 = GT centroid 未匹配（落在空隙或範圍外）
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

RESULT_DIR = Path("/Volumes/SSD/plan_a/submission_bioinformatics/results/enact_crc_f1")
CROP_X0, CROP_Y0 = 5154, 4635   # ENACT local coords of crop origin

# ── 載入資料 ────────────────────────────────────────────────
print("載入 mask...")
mask = np.load(RESULT_DIR / "mcseg_mask.npy")   # (H, W) int32

print("載入 GT matched...")
gt = pd.read_csv(RESULT_DIR / "gt_matched.csv")
# cell_x/y 是 ENACT local coords；轉為 mask 內像素座標
gt = gt.dropna(subset=["cell_x", "cell_y"])
gt["col_local"] = (gt["cell_x"] - CROP_X0).astype(int)
gt["row_local"] = (gt["cell_y"] - CROP_Y0).astype(int)

matched   = gt[gt["mcseg_cell_id"] > 0]
unmatched = gt[gt["mcseg_cell_id"] == 0]

H, W = mask.shape
print(f"mask: {H}×{W},  GT total: {len(gt)},  matched: {len(matched)},  unmatched: {len(unmatched)}")

# ── 建立 RGB 視覺化 ────────────────────────────────────────
print("建立 mask RGB...")
rng = np.random.default_rng(42)
n_cells = int(mask.max()) + 1
colors = rng.integers(60, 220, size=(n_cells, 3), dtype=np.uint8)
colors[0] = [20, 20, 20]   # background 深灰

# downsample 4x for speed
step = 4
mask_ds = mask[::step, ::step]
rgb = colors[mask_ds]   # (H//4, W//4, 3)

# ── 全圖 overview ──────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(20, 10), dpi=150)

ax = axes[0]
ax.imshow(rgb, origin="upper", interpolation="nearest")
col_m = matched["col_local"].values / step
row_m = matched["row_local"].values / step
col_u = unmatched["col_local"].values / step
row_u = unmatched["row_local"].values / step

ax.scatter(col_u, row_u, s=1.5, c="#ff4444", alpha=0.4, linewidths=0, label=f"Unmatched ({len(unmatched):,})")
ax.scatter(col_m, row_m, s=1.5, c="#44ff88", alpha=0.6, linewidths=0, label=f"Matched ({len(matched):,})")
ax.set_title(f"Full crop — GT centroids on MCseg mask\n(match rate {len(matched)/len(gt)*100:.1f}%)", fontsize=11)
ax.legend(markerscale=5, fontsize=8, loc="upper right")
ax.axis("off")

# ── 局部放大（左上 1/4）──────────────────────────────────
ax2 = axes[1]
H4, W4 = H // 4, W // 4
mask_zoom = mask[:H4, :W4]
rgb_zoom  = colors[mask_zoom]
ax2.imshow(rgb_zoom, origin="upper", interpolation="nearest")

in_zoom_m = matched[(matched["row_local"] < H4) & (matched["col_local"] < W4)]
in_zoom_u = unmatched[(unmatched["row_local"] < H4) & (unmatched["col_local"] < W4)]
ax2.scatter(in_zoom_u["col_local"], in_zoom_u["row_local"], s=4, c="#ff4444", alpha=0.6, linewidths=0, label=f"Unmatched ({len(in_zoom_u):,})")
ax2.scatter(in_zoom_m["col_local"], in_zoom_m["row_local"], s=4, c="#44ff88", alpha=0.8, linewidths=0, label=f"Matched ({len(in_zoom_m):,})")
ax2.set_title(f"Zoom: top-left quarter ({W4}×{H4} px)", fontsize=11)
ax2.legend(markerscale=3, fontsize=8, loc="upper right")
ax2.axis("off")

plt.tight_layout()
out = RESULT_DIR / "fig_gt_overlay.png"
fig.savefig(out, dpi=150, bbox_inches="tight")
print(f"Saved: {out}")
plt.close()
