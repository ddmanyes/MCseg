import numpy as np
import matplotlib.pyplot as plt
import json
from pathlib import Path
import cv2

# 設定路徑
ROOT = Path("/Volumes/SSD/plan_a/crc_he_seg")
ROIS_DIR = ROOT / "results" / "rois"
MASKS_DIR = ROOT / "results"
OUT_PATH = Path("/Volumes/SSD/plan_a/manuscript/figures/02_methods/fig2b_crc_selected.png")

# 僅精選 ROI2 與 ROI4
SELECTED_ROIS = ["roi2", "roi4"]
AP_SCORES = {"roi2": 0.5326, "roi4": 0.4885}

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
        overlay[outer | r_diff | d_diff, :3] = 1.0
        overlay[outer | r_diff | d_diff, 3] = boundary_alpha
    return overlay

def blend(base, over):
    bg = base.astype(np.float32) / 255.0
    fg, a = over[..., :3], over[..., 3:4]
    return np.clip((fg * a + bg * (1.0 - a)) * 255, 0, 255).astype(np.uint8)

print(f"Generating Main Text Fig 2b (ROI2 & ROI4 only)...")
fig, axes = plt.subplots(2, 3, figsize=(16, 11), dpi=200) # 更高的 DPI 保證正文品質
plt.subplots_adjust(hspace=0.05, wspace=0.05)

for i, roi_name in enumerate(SELECTED_ROIS):
    he = np.load(ROIS_DIR / f"{roi_name}_he.npy")
    if he.dtype != np.uint8: he = (np.clip(he, 0, 1) * 255).astype(np.uint8)
    pred = np.load(MASKS_DIR / f"{roi_name}_pred_mask.npy")
    gt   = np.load(MASKS_DIR / f"{roi_name}_gt_mask_rendered.npy")
    
    img_pred = blend(he, mask_to_rgba_overlay(pred, BLUE, alpha=0.35))
    img_gt   = blend(he, mask_to_rgba_overlay(gt, ORANGE, alpha=0.45))
    
    axes[i, 0].imshow(he)
    axes[i, 0].set_ylabel(f"{roi_name.upper()}\n(AP={AP_SCORES[roi_name]:.4f})", fontsize=14, fontweight='bold')
    if i == 0: axes[i, 0].set_title("H&E Image", fontsize=16, fontweight='bold')
    axes[i, 0].set_xticks([]); axes[i, 0].set_yticks([])

    axes[i, 1].imshow(img_pred)
    if i == 0: axes[i, 1].set_title("V12 Prediction", fontsize=16, fontweight='bold')
    axes[i, 1].set_xticks([]); axes[i, 1].set_yticks([])

    axes[i, 2].imshow(img_gt)
    if i == 0: axes[i, 2].set_title("SR GT (Official Reference)", fontsize=16, fontweight='bold', color='red')
    axes[i, 2].set_xticks([]); axes[i, 2].set_yticks([])

plt.savefig(OUT_PATH, bbox_inches='tight')
plt.close()
print(f"✅ Main Text Figure 2b generated: {OUT_PATH}")
