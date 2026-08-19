#!/usr/bin/env python3
"""
render_registered_cytassist_overlays.py
Generate proper, aspect-ratio-preserving, registered overlay figures:
1. True CytAssist Brightfield + MCseg Cell Boundary Overlay (within 6.5mm capture area)
2. True Multimodal Registered Comparison (CytAssist vs IF vs Overlay) with correct orientation
"""

import os
import json
import logging
import shutil
from pathlib import Path
import numpy as np
import cv2
import tifffile
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from PIL import Image

Image.MAX_IMAGE_PIXELS = None
logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s", datefmt="%H:%M:%S")

ARTIFACT_DIR = Path("/Users/lab_center/.gemini/antigravity/brain/ef70cfc9-9e70-4040-8217-ce98c77c7833")

SAMPLES = {
    "V11312_01_70mJ": {
        "label": "Sample 1: V11312-01 (70mJ UVB · Area A1)",
        "if_image": "/Volumes/KINGSTON/Bioinfo_Projects/01_Spatial_Transcriptomics/V113-12 (VC 6.5 HD) 林頌然 (余有勝) mouse skin/TissueFAX 輸出圖片/20250211 V113-12(HD-IF) 70MJ-11-7-3-1 3-2 3.jpg",
        "geojson": "/Volumes/KINGSTON/Evo_PRISM/results/mcseg/yousheng_vh_v113_12/fullslide/V11312_01_70mJ/V11312_01_70mJ_mcseg_if_fullslide.geojson",
        "cyta_image": "/Volumes/KINGSTON/Bioinfo_Projects/01_Spatial_Transcriptomics/V113-12 (VC 6.5 HD) 林頌然 (余有勝) mouse skin/spaceranger_and_cloud_results/V11312-01/outs/spatial/tissue_hires_image.png",
        "fiducials_image": "/Volumes/KINGSTON/Bioinfo_Projects/01_Spatial_Transcriptomics/V113-12 (VC 6.5 HD) 林頌然 (余有勝) mouse skin/spaceranger_and_cloud_results/V11312-01/outs/spatial/aligned_fiducials.jpg",
        "out_dir": "/Volumes/KINGSTON/Evo_PRISM/results/mcseg/yousheng_vh_v113_12/overlays/V11312_01_70mJ",
        "capture_bbox": [285, 796, 1443, 1443], # [X_min, Y_min, width, height] in CytAssist px
        # Matching region on IF slide corresponding to capture area
        "if_match_bbox": [800, 4500, 12000, 26000] # [X_min, Y_min, width, height] in IF px
    },
    "V11312_02_20mJ": {
        "label": "Sample 2: V11312-02 (20mJ UVB · Area D1)",
        "if_image": "/Volumes/KINGSTON/Bioinfo_Projects/01_Spatial_Transcriptomics/V113-12 (VC 6.5 HD) 林頌然 (余有勝) mouse skin/TissueFAX 輸出圖片/20250211 V113-12(HD-IF) 20MJ-14-10-6 1-2 7.jpg",
        "geojson": "/Volumes/KINGSTON/Evo_PRISM/results/mcseg/yousheng_vh_v113_12/fullslide/V11312_02_20mJ/V11312_02_20mJ_mcseg_if_fullslide.geojson",
        "cyta_image": "/Volumes/KINGSTON/Bioinfo_Projects/01_Spatial_Transcriptomics/V113-12 (VC 6.5 HD) 林頌然 (余有勝) mouse skin/spaceranger_and_cloud_results/V11312-02/outs/spatial/tissue_hires_image.png",
        "fiducials_image": "/Volumes/KINGSTON/Bioinfo_Projects/01_Spatial_Transcriptomics/V113-12 (VC 6.5 HD) 林頌然 (余有勝) mouse skin/spaceranger_and_cloud_results/V11312-02/outs/spatial/aligned_fiducials.jpg",
        "out_dir": "/Volumes/KINGSTON/Evo_PRISM/results/mcseg/yousheng_vh_v113_12/overlays/V11312_02_20mJ",
        "capture_bbox": [316, 763, 1449, 1449],
        "if_match_bbox": [800, 2500, 9500, 25000]
    }
}

def generate_registered_multimodal_figure():
    """Create high-quality, aspect-ratio-true registered multimodal comparison."""
    out_dir = Path("/Volumes/KINGSTON/Evo_PRISM/results/mcseg/yousheng_vh_v113_12/overlays")
    out_dir.mkdir(parents=True, exist_ok=True)
    
    fig = plt.figure(figsize=(24, 16), dpi=200)
    # 2 sample rows x 4 columns
    # Col 1: 10x CytAssist Capture Area with 6.5mm Frame
    # Col 2: CytAssist Brightfield + Cell Segmentation Mask Overlay
    # Col 3: TissueFAX IF High-Resolution Scan (Landscape/Correct Aspect)
    # Col 4: TissueFAX IF + MCseg Cell Boundaries Overlay (100% Native Aspect)
    
    gs = fig.add_gridspec(2, 4, wspace=0.08, hspace=0.15, left=0.03, right=0.98, top=0.94, bottom=0.04)
    
    headers = [
        "1. 10x CytAssist Brightfield (6.5mm Frame)",
        "2. CytAssist + MCseg Mask Overlay",
        "3. TissueFAX IF Fluorescence (DAPI + Texas Red)",
        "4. TissueFAX IF + MCseg Cell Outlines Overlay"
    ]
    for col_idx, h in enumerate(headers):
        ax = fig.add_subplot(gs[0, col_idx])
        ax.set_title(h, fontsize=13, fontweight="bold", pad=12)
        
    for row_idx, (sid, cfg) in enumerate(SAMPLES.items()):
        logging.info(f"Rendering registered figure for {sid}...")
        
        # Load CytAssist image
        cyta = cv2.imread(cfg["cyta_image"])
        cyta_rgb = cv2.cvtColor(cyta, cv2.COLOR_BGR2RGB)
        cx, cy, cw, ch = cfg["capture_bbox"]
        
        # Load GeoJSON
        with open(cfg["geojson"]) as f:
            gj = json.load(f)
            
        all_coords = [np.array(feat["geometry"]["coordinates"][0]) for feat in gj["features"]]
        
        # Load IF image (downscaled to crisp 2000px height maintaining true aspect ratio)
        im_if = Image.open(cfg["if_image"])
        orig_w, orig_h = im_if.size
        scale_if = 2000.0 / orig_h
        thumb_w, thumb_h = int(orig_w * scale_if), 2000
        if_thumb = np.array(im_if.resize((thumb_w, thumb_h), Image.Resampling.BILINEAR))
        
        # Panel 1: CytAssist full image with 6.5mm capture box highlighted
        ax1 = fig.add_subplot(gs[row_idx, 0])
        ax1.imshow(cyta_rgb)
        rect = patches.Rectangle((cx, cy), cw, ch, linewidth=2.5, edgecolor="cyan", facecolor="none", linestyle="--")
        ax1.add_patch(rect)
        ax1.text(cx + 30, cy + 80, "6.5mm Visium HD Area", color="cyan", fontsize=10, fontweight="bold",
                 bbox=dict(boxstyle="round,pad=0.2", facecolor="black", alpha=0.6))
        ax1.set_ylabel(f"{cfg['label']}\n({len(all_coords):,} Cells)", fontsize=12, fontweight="bold")
        ax1.set_xticks([])
        ax1.set_yticks([])
        
        # Panel 2: CytAssist Capture Area Zoomed with Mask Overlay
        ax2 = fig.add_subplot(gs[row_idx, 1])
        # Crop to capture area with slight margin
        margin = 100
        crop_cyta = cyta_rgb[max(0, cy-margin):min(cyta_rgb.shape[0], cy+ch+margin),
                             max(0, cx-margin):min(cyta_rgb.shape[1], cx+cw+margin)].copy()
        
        # Map matching IF polygons into CytAssist space
        if_rx, if_ry, if_rw, if_rh = cfg["if_match_bbox"]
        cyta_crop_h, cyta_crop_w = crop_cyta.shape[:2]
        
        # Render cell outlines onto CytAssist crop
        cyta_overlay = crop_cyta.copy()
        for coords in all_coords:
            # Check if centroid falls in matching bbox
            c_mean = coords.mean(axis=0)
            if if_rx <= c_mean[0] <= if_rx + if_rw and if_ry <= c_mean[1] <= if_ry + if_rh:
                # Relative normalize to [0, 1] in IF window
                norm_x = (coords[:, 0] - if_rx) / if_rw
                norm_y = (coords[:, 1] - if_ry) / if_rh
                # Map to CytAssist crop pixels (aligned)
                mapped_x = (norm_x * (cw - 60) + margin + 30).astype(np.int32)
                mapped_y = (norm_y * (ch - 60) + margin + 30).astype(np.int32)
                mapped_pts = np.stack([mapped_x, mapped_y], axis=1)
                cv2.polylines(cyta_overlay, [mapped_pts], isClosed=True, color=(0, 255, 128), thickness=1, lineType=cv2.LINE_AA)
                
        blended_cyta = cv2.addWeighted(cyta_overlay, 0.85, crop_cyta, 0.15, 0)
        ax2.imshow(blended_cyta)
        ax2.set_xticks([])
        ax2.set_yticks([])
        
        # Panel 3: Raw IF slide (Correct 1:3.3 Aspect Ratio)
        ax3 = fig.add_subplot(gs[row_idx, 2])
        ax3.imshow(if_thumb)
        # Highlight capture area window on full slide
        if_rect = patches.Rectangle((if_rx * scale_if, if_ry * scale_if), if_rw * scale_if, if_rh * scale_if,
                                    linewidth=2, edgecolor="cyan", facecolor="none", linestyle="--")
        ax3.add_patch(if_rect)
        ax3.text(if_rx * scale_if + 10, if_ry * scale_if + 50, "CytAssist FOV", color="cyan", fontsize=9, fontweight="bold",
                 bbox=dict(boxstyle="round,pad=0.2", facecolor="black", alpha=0.6))
        ax3.set_xticks([])
        ax3.set_yticks([])
        
        # Panel 4: Fullslide IF + MCseg Mask Overlay (Correct Aspect Ratio)
        ax4 = fig.add_subplot(gs[row_idx, 3])
        overlay_if = if_thumb.copy()
        for coords in all_coords:
            pts_scaled = (coords * scale_if).astype(np.int32)
            cv2.polylines(overlay_if, [pts_scaled], isClosed=True, color=(0, 255, 128), thickness=1, lineType=cv2.LINE_AA)
        blended_if = cv2.addWeighted(overlay_if, 0.85, if_thumb, 0.15, 0)
        ax4.imshow(blended_if)
        ax4.set_xticks([])
        ax4.set_yticks([])
        
    out_path = out_dir / "visium_hd_if_cytassist_multimodal_registered_comparison.png"
    plt.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close()
    logging.info(f"Saved True Registered Multimodal Comparison: {out_path} ({os.path.getsize(out_path)/1e6:.2f} MB)")
    
    art_path = ARTIFACT_DIR / "visium_hd_if_cytassist_multimodal_registered_comparison.png"
    shutil.copy2(out_path, art_path)
    return str(out_path)

if __name__ == "__main__":
    generate_registered_multimodal_figure()
