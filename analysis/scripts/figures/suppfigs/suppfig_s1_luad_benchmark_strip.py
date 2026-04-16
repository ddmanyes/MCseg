import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import cv2
from pathlib import Path

# ── 設定與路徑 ────────────────────────────────────────────────────────
MASK_DIR = Path("/Volumes/SSD/plan_a/xenium_he_seg/results/masks")
PRED_DIR = Path("/Volumes/SSD/plan_a/xenium_visiumhd_comparison/results/benchmark_v12")
OUT_DIR     = Path("/Volumes/SSD/plan_a/xenium_visiumhd_comparison/results/benchmark_v12/figures")
OUT_DIR_SUB = Path("/Volumes/SSD/plan_a/submission_bioinformatics/supplementary")
OUT_DIR.mkdir(parents=True, exist_ok=True)
OUT_DIR_SUB.mkdir(parents=True, exist_ok=True)

PIXEL_SIZE_UM = 0.2737   # LUAD VisiumHD H&E
SCALE_UM      = 50       # µm

# ROI 配置 (PQ 2026-03-25)
ROIS = [
    ("roi1", "roi1_pred.npy", 0.6615, "#2ecc71", "Tumor boundary",         800, 1600),
    ("6",    "6_pred.npy",    0.6097, "#2ecc71", "Tumor core",              400, 1200),
    ("2",    "2_pred.npy",    0.5534, "#2ecc71", "Tumor stroma",            400,  800),
    ("3",    "3_pred.npy",    0.5229, "#2ecc71", "Mixed tumor-stroma",        0, 1600),
    ("4",    "4_pred.npy",    0.5101, "#2ecc71", "Normal-tumor interface",  400, 1200),
    ("5",    "5_pred.npy",    0.4683, "#f39c12", "Alveolar region",         400,  800),
]

BLUE   = (70,  130, 220) # V12
ORANGE = (210,  90,   0) # Xenium GT
ZOOM_SZ = 600

# ── 幫助函式 ──────────────────────────────────────────────────────────
def mask_to_rgba_overlay(mask, color_rgb, alpha=0.32, draw_boundary=True, boundary_alpha=0.9):
    """將 instance mask 轉為帶有「白色邊界強化」的高清晰度 RGBA 疊層 (與 CRC 風格一致)"""
    h, w = mask.shape
    overlay = np.zeros((h, w, 4), dtype=np.float32)
    fill = mask > 0
    overlay[fill, :3] = np.array(color_rgb) / 255.0
    overlay[fill, 3] = alpha
    
    if draw_boundary:
        m32 = mask.astype(np.int32)
        binary = (mask > 0).astype(np.uint8)
        eroded = cv2.erode(binary, np.ones((3,3), np.uint8), iterations=1)
        outer = (binary - eroded).astype(bool)
        # 細胞對細胞邊界
        r_diff = np.zeros((h, w), dtype=bool); r_diff[:, :-1] = m32[:, 1:] != m32[:, :-1]
        d_diff = np.zeros((h, w), dtype=bool); d_diff[:-1, :]  = m32[1:, :] != m32[:-1, :]
        inter_bound = (r_diff | d_diff) & fill
        bound = outer | inter_bound
        # 邊界強化為「白色」確保對比度
        overlay[bound, :3] = 1.0 # White
        overlay[bound, 3] = boundary_alpha
    return overlay

def blend_images(base, overlay):
    bg = base.astype(np.float32) / 255.0
    fg, a = overlay[..., :3], overlay[..., 3:4]
    return np.clip((fg * a + bg * (1.0 - a)) * 255, 0, 255).astype(np.uint8)

def add_scale_bar(ax, img_h, img_w):
    """White 50 µm scale bar at bottom-right corner."""
    scale_px = SCALE_UM / PIXEL_SIZE_UM   # ≈ 183 px
    margin_x = img_w * 0.05
    margin_y = img_h * 0.06
    bar_y  = img_h - margin_y
    bar_x1 = img_w - margin_x
    bar_x0 = bar_x1 - scale_px
    ax.plot([bar_x0, bar_x1], [bar_y, bar_y],
            color="white", linewidth=3, solid_capstyle="butt", zorder=10)
    ax.text((bar_x0 + bar_x1) / 2, bar_y - img_h * 0.025,
            f"{SCALE_UM} µm",
            color="white", ha="center", va="bottom",
            fontsize=9, fontweight="bold", zorder=10)

def render_row(axes, rname, pq, desc, color_code, zy, zx):
    print(f"  Rendering {rname}...")
    he   = np.load(MASK_DIR / f"{rname}_he_crop.npy")
    gt   = np.load(MASK_DIR / f"{rname}_xenium_gt_mask.npy")
    pred = np.load(PRED_DIR / f"{rname}_pred.npy")
    
    if he.dtype != np.uint8: he = (np.clip(he, 0, 1) * 255).astype(np.uint8)
    h, w = he.shape[:2]
    
    # 疊層合成
    ov_v12 = mask_to_rgba_overlay(pred, BLUE,   alpha=0.35)
    ov_gt  = mask_to_rgba_overlay(gt,   ORANGE, alpha=0.45)
    
    img_v12 = blend_images(he, ov_v12)
    img_gt  = blend_images(he, ov_gt)
    
    # Col 0: H&E Raw
    axes[0].imshow(he)
    axes[0].set_ylabel(f"ROI {rname.replace('roi','')}\n{desc}\nPQ={pq:.3f}", fontsize=11, fontweight='bold', color=color_code)
    if rname == "roi1": axes[0].set_title("H&E Image", fontsize=15, fontweight='bold')
    axes[0].set_xticks([]); axes[0].set_yticks([])
    add_scale_bar(axes[0], h, w)

    # Col 1: MCseg (Sharp Boundaries)
    axes[1].imshow(img_v12)
    if rname == "roi1": axes[1].set_title("MCseg", fontsize=15, fontweight='bold')
    axes[1].set_xticks([]); axes[1].set_yticks([])
    add_scale_bar(axes[1], h, w)

    # Col 2: Xenium GT (Sharp Boundaries)
    axes[2].imshow(img_gt)
    if rname == "roi1": axes[2].set_title("Xenium GT", fontsize=15, fontweight='bold', color='red')
    axes[2].set_xticks([]); axes[2].set_yticks([])
    add_scale_bar(axes[2], h, w)

    # Col 3: Zoom in
    y0, x0 = max(0, zy), max(0, zx)
    y1, x1 = min(h, y0+ZOOM_SZ), min(w, x0+ZOOM_SZ)
    zoom_h  = y1 - y0
    zoom_w  = min(ZOOM_SZ, w - x0)
    zoom_v12 = img_v12[y0:y1, x0:x0+ZOOM_SZ]
    axes[3].imshow(zoom_v12)
    if rname == "roi1": axes[3].set_title("Zoom in", fontsize=15, fontweight='bold')
    axes[3].set_xticks([]); axes[3].set_yticks([])
    add_scale_bar(axes[3], zoom_h, zoom_w)

# ── 起跑渲染 ──────────────────────────────────────────────────────────
fig, axes = plt.subplots(6, 4, figsize=(20, 36), dpi=300) # 300 DPI Ultra-HD
plt.subplots_adjust(hspace=0.08, wspace=0.10)

for i, (name, fname, pq, color, desc, zy, zx) in enumerate(ROIS):
    render_row(axes[i], name, pq, desc, color, zy, zx)

out_path = OUT_DIR / "v12_vs_gt_strip_ULTRAHD.png"
plt.savefig(out_path, bbox_inches='tight', facecolor='white')
out_supp = OUT_DIR_SUB / "SuppFigS1.png"
plt.savefig(out_supp, bbox_inches='tight', facecolor='white')
plt.close()
print(f"✅ ULTRA-HD 6-ROI Strip generated at: {out_path}")
print(f"✅ SuppFigS1 saved at: {out_supp}")
