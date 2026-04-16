#!/usr/bin/env python3
"""
viz_comparison.py (Minimalist 1x3 Version)
三格對比圖：H&E | MCseg v2 預測 | Xenium GT
設計：移除 Zoom View 與選取框，僅展示純淨的全景對照。
"""

import numpy as np
import cv2
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

ROOT    = Path("/Volumes/SSD/plan_a/autoresearch_seg")
DATA    = ROOT / "data"
RESULTS = ROOT / "results"
PRED_CACHE = RESULTS / "pred_mask_best.npy"
OUT_PNG    = Path("/Volumes/SSD/plan_a/manuscript/figures/fig1/fig1b_comparison_3panel.png")

# 1. 數據載入
img = np.load(DATA / "he_patch.npy")
gt  = np.load(DATA / "gt_mask.npy")
pred = np.load(PRED_CACHE)

BLUE   = (70,  130, 220)
ORANGE = (210,  90,   0)

def mask_to_rgba_overlay(mask, color_rgb, alpha=0.35, draw_boundary=True, boundary_alpha=0.9):
    h, w = mask.shape
    overlay = np.zeros((h, w, 4), dtype=np.float32)
    fill    = mask > 0
    overlay[fill, :3] = np.array(color_rgb) / 255.0
    overlay[fill, 3] = alpha
    if draw_boundary:
        m32 = mask.astype(np.int32)
        binary = (mask > 0).astype(np.uint8)
        eroded = cv2.erode(binary, np.ones((3,3), np.uint8), iterations=1)
        outer = (binary - eroded).astype(bool)
        r_diff = np.zeros((h, w), dtype=bool); r_diff[:, :-1] = m32[:, 1:] != m32[:, :-1]
        d_diff = np.zeros((h, w), dtype=bool); d_diff[:-1, :]  = m32[1:, :] != m32[:-1, :]
        inter_bound = (r_diff | d_diff) & fill
        overlay[outer | inter_bound, :3] = 1.0
        overlay[outer | inter_bound, 3] = boundary_alpha
    return overlay

def blend(base_rgb, overlay_rgba):
    bg  = base_rgb.astype(np.float32) / 255.0
    fg  = overlay_rgba[..., :3]
    a   = overlay_rgba[..., 3:4]
    out = fg * a + bg * (1.0 - a)
    return np.clip(out * 255, 0, 255).astype(np.uint8)

img_pred = blend(img, mask_to_rgba_overlay(pred, BLUE,   alpha=0.35))
img_gt   = blend(img, mask_to_rgba_overlay(gt,   ORANGE, alpha=0.45))

# Scale bar parameters (VisiumHD H&E fullres: 0.2737 µm/px)
PIXEL_SIZE_UM = 0.2737   # µm per pixel
SCALE_UM      = 50       # µm
scale_px = SCALE_UM / PIXEL_SIZE_UM  # ≈ 365 px

def add_scale_bar(ax, img_h, img_w):
    """Add a 100 µm white scale bar at bottom-right corner."""
    margin = img_w * 0.04          # 4% margin from edge
    bar_y  = img_h * 0.93          # 93% down the image
    x_end  = img_w - margin
    x_start = x_end - scale_px
    ax.plot([x_start, x_end], [bar_y, bar_y],
            color="white", linewidth=4, solid_capstyle="butt")
    ax.text((x_start + x_end) / 2, bar_y - img_h * 0.025,
            f"{SCALE_UM} µm",
            color="white", ha="center", va="bottom",
            fontsize=13, fontweight="bold",
            bbox=dict(facecolor="none", edgecolor="none", pad=0))

# 2. 繪圖：極簡 1x3 佈局
PANEL_TITLES = ["H&E Image", "MCseg", "Xenium GT"]
fig, axes = plt.subplots(1, 3, figsize=(21, 7), dpi=150)
plt.subplots_adjust(wspace=0.03)

for ax, title, panel in zip(axes, PANEL_TITLES, [img, img_pred, img_gt]):
    ax.imshow(panel)
    ax.set_title(title, fontsize=18, fontweight="bold", pad=15)
    ax.axis("off")

# Scale bar on all three panels
h, w = img.shape[:2]
for ax in axes:
    add_scale_bar(ax, h, w)

# 專業圖例 (放置於底部)
from matplotlib.patches import Patch
legend_elements = [Patch(facecolor=np.array(BLUE)/255, alpha=0.7, label="MCseg"),
                   Patch(facecolor=np.array(ORANGE)/255, alpha=0.7, label="Xenium GT")]
fig.legend(handles=legend_elements, loc="lower center", ncol=2, fontsize=14, frameon=False, bbox_to_anchor=(0.5, 0.05))

plt.savefig(OUT_PNG, bbox_inches="tight", facecolor="white")
print(f"✅ Minimalist 3-panel figure generated: {OUT_PNG}")
