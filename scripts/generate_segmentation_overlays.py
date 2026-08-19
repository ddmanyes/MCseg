#!/usr/bin/env python3
"""
generate_segmentation_overlays.py
Generate high-resolution full-slide macro overlays and 1:1 pixel zoom-in comparison panels
for Visium HD IF samples (V11312-01 & V11312-02), plus CytAssist multimodal comparison.
"""

import os
import sys
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
from PIL import Image

Image.MAX_IMAGE_PIXELS = None
logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s", datefmt="%H:%M:%S")

ARTIFACT_DIR = Path("/Users/lab_center/.gemini/antigravity/brain/ef70cfc9-9e70-4040-8217-ce98c77c7833")

SAMPLES = {
    "V11312_01_70mJ": {
        "if_image": "/Volumes/KINGSTON/Bioinfo_Projects/01_Spatial_Transcriptomics/V113-12 (VC 6.5 HD) 林頌然 (余有勝) mouse skin/TissueFAX 輸出圖片/20250211 V113-12(HD-IF) 70MJ-11-7-3-1 3-2 3.jpg",
        "geojson": "/Volumes/KINGSTON/Evo_PRISM/results/mcseg/yousheng_vh_v113_12/fullslide/V11312_01_70mJ/V11312_01_70mJ_mcseg_if_fullslide.geojson",
        "cyta_image": "/Volumes/KINGSTON/Bioinfo_Projects/01_Spatial_Transcriptomics/V113-12 (VC 6.5 HD) 林頌然 (余有勝) mouse skin/spaceranger_and_cloud_results/V11312-01/outs/spatial/cytassist_image.tiff",
        "out_dir": "/Volumes/KINGSTON/Evo_PRISM/results/mcseg/yousheng_vh_v113_12/overlays/V11312_01_70mJ",
        "rois": [
            {"name": "ROI 1: Dense Hair Follicle Cluster", "bbox": [1800, 5400, 1200, 1200]},
            {"name": "ROI 2: Epidermal-Dermal Junction", "bbox": [9000, 9000, 1200, 1200]},
            {"name": "ROI 3: Dermal Interstitial Cells", "bbox": [8400, 21000, 1200, 1200]},
            {"name": "ROI 4: Basal Layer & Sebaceous Gland", "bbox": [9000, 3000, 1200, 1200]},
        ]
    },
    "V11312_02_20mJ": {
        "if_image": "/Volumes/KINGSTON/Bioinfo_Projects/01_Spatial_Transcriptomics/V113-12 (VC 6.5 HD) 林頌然 (余有勝) mouse skin/TissueFAX 輸出圖片/20250211 V113-12(HD-IF) 20MJ-14-10-6 1-2 7.jpg",
        "geojson": "/Volumes/KINGSTON/Evo_PRISM/results/mcseg/yousheng_vh_v113_12/fullslide/V11312_02_20mJ/V11312_02_20mJ_mcseg_if_fullslide.geojson",
        "cyta_image": "/Volumes/KINGSTON/Bioinfo_Projects/01_Spatial_Transcriptomics/V113-12 (VC 6.5 HD) 林頌然 (余有勝) mouse skin/spaceranger_and_cloud_results/V11312-02/outs/spatial/cytassist_image.tiff",
        "out_dir": "/Volumes/KINGSTON/Evo_PRISM/results/mcseg/yousheng_vh_v113_12/overlays/V11312_02_20mJ",
        "rois": [
            {"name": "ROI 1: Epidermis & Hair Follicles", "bbox": [9000, 21600, 1200, 1200]},
            {"name": "ROI 2: Interstitial Fibroblasts", "bbox": [5400, 3000, 1200, 1200]},
            {"name": "ROI 3: Subcutaneous Tissue & Muscle", "bbox": [6000, 18000, 1200, 1200]},
            {"name": "ROI 4: Deep Dermal Cellular Strip", "bbox": [2400, 24600, 1200, 1200]},
        ]
    }
}

def render_fullslide_macro_overlay(sample_id, cfg):
    """Render fullslide downscaled overlay with mask outlines."""
    out_dir = Path(cfg["out_dir"])
    out_dir.mkdir(parents=True, exist_ok=True)
    
    logging.info(f"[{sample_id}] Loading IF image for macro overview...")
    img = Image.open(cfg["if_image"])
    orig_w, orig_h = img.size
    
    scale = 4000.0 / max(orig_w, orig_h)
    new_w, new_h = int(orig_w * scale), int(orig_h * scale)
    logging.info(f"[{sample_id}] Resizing from ({orig_w}, {orig_h}) to ({new_w}, {new_h}) [scale={scale:.4f}]...")
    img_thumb = np.array(img.resize((new_w, new_h), Image.Resampling.BILINEAR))
    
    logging.info(f"[{sample_id}] Loading GeoJSON...")
    with open(cfg["geojson"]) as f:
        gj = json.load(f)
        
    logging.info(f"[{sample_id}] Drawing {len(gj['features'])} cell boundaries onto overlay canvas...")
    overlay = img_thumb.copy()
    
    for feat in gj["features"]:
        coords = np.array(feat["geometry"]["coordinates"][0])
        coords_scaled = (coords * scale).astype(np.int32)
        cv2.polylines(overlay, [coords_scaled], isClosed=True, color=(0, 255, 128), thickness=1, lineType=cv2.LINE_AA)
        
    blended = cv2.addWeighted(overlay, 0.85, img_thumb, 0.15, 0)
    
    out_path = out_dir / f"{sample_id}_fullslide_macro_overlay.jpg"
    cv2.imwrite(str(out_path), cv2.cvtColor(blended, cv2.COLOR_RGB2BGR), [cv2.IMWRITE_JPEG_QUALITY, 92])
    logging.info(f"[{sample_id}] Saved macro overlay: {out_path} ({os.path.getsize(out_path)/1e6:.2f} MB)")
    
    # Copy to artifact dir for viewing
    art_path = ARTIFACT_DIR / f"{sample_id}_fullslide_macro_overlay.jpg"
    shutil.copy2(out_path, art_path)
    return str(out_path)

def render_roi_multi_panel(sample_id, cfg):
    """Render 4 high-resolution ROI zoom-in panels with 3 columns (Raw, Mask, Overlay)."""
    out_dir = Path(cfg["out_dir"])
    out_dir.mkdir(parents=True, exist_ok=True)
    
    logging.info(f"[{sample_id}] Loading full GeoJSON for spatial ROI clipping...")
    with open(cfg["geojson"]) as f:
        gj = json.load(f)
        
    all_polys = []
    for feat in gj["features"]:
        coords = np.array(feat["geometry"]["coordinates"][0])
        cid = feat["properties"]["cell_id"]
        all_polys.append((cid, coords))
        
    logging.info(f"[{sample_id}] Total polygons loaded: {len(all_polys)}")
    im_full = Image.open(cfg["if_image"])
    
    fig, axes = plt.subplots(4, 3, figsize=(18, 24), dpi=200)
    plt.subplots_adjust(wspace=0.03, hspace=0.08, left=0.04, right=0.98, top=0.96, bottom=0.02)
    
    cols = ["Raw IF (DAPI + Texas Red)", "Cellpose Multichannel Mask", "Overlay (1:1 Native Resolution)"]
    for col_idx, col_name in enumerate(cols):
        axes[0, col_idx].set_title(col_name, fontsize=15, fontweight="bold", pad=12)
        
    cmap = matplotlib.colormaps["gist_rainbow"].resampled(100)
    
    for row_idx, roi in enumerate(cfg["rois"]):
        rx, ry, rw, rh = roi["bbox"]
        roi_name = roi["name"]
        logging.info(f"[{sample_id}] Processing {roi_name} at bbox=({rx}, {ry}, {rw}, {rh})...")
        
        crop_pil = im_full.crop((rx, ry, rx + rw, ry + rh))
        crop_np = np.array(crop_pil)
        
        roi_cells = []
        for cid, coords in all_polys:
            min_x, min_y = coords.min(axis=0)
            max_x, max_y = coords.max(axis=0)
            if max_x >= rx and min_x <= rx + rw and max_y >= ry and min_y <= ry + rh:
                rel_coords = coords - np.array([rx, ry])
                roi_cells.append((cid, rel_coords))
                
        logging.info(f"  -> Found {len(roi_cells)} cells in ROI")
        
        # 1. Col 0: Raw IF
        ax_raw = axes[row_idx, 0]
        ax_raw.imshow(crop_np)
        ax_raw.set_ylabel(f"{roi_name}\n({len(roi_cells)} cells)", fontsize=12, fontweight="bold")
        ax_raw.set_xticks([])
        ax_raw.set_yticks([])
        
        # 2. Col 1: Mask only
        ax_mask = axes[row_idx, 1]
        mask_canvas = np.zeros_like(crop_np)
        for i, (cid, rel_coords) in enumerate(roi_cells):
            pts = rel_coords.astype(np.int32)
            c = (np.array(cmap(i % 100)[:3]) * 255).astype(int).tolist()
            cv2.fillPoly(mask_canvas, [pts], c)
            cv2.polylines(mask_canvas, [pts], isClosed=True, color=(255, 255, 255), thickness=1, lineType=cv2.LINE_AA)
        ax_mask.imshow(mask_canvas)
        ax_mask.set_xticks([])
        ax_mask.set_yticks([])
        
        # 3. Col 2: Overlay
        ax_ov = axes[row_idx, 2]
        overlay_canvas = crop_np.copy()
        for i, (cid, rel_coords) in enumerate(roi_cells):
            pts = rel_coords.astype(np.int32)
            cv2.polylines(overlay_canvas, [pts], isClosed=True, color=(0, 255, 128), thickness=2, lineType=cv2.LINE_AA)
        
        ov_blended = cv2.addWeighted(overlay_canvas, 0.85, crop_np, 0.15, 0)
        ax_ov.imshow(ov_blended)
        ax_ov.set_xticks([])
        ax_ov.set_yticks([])
        
    out_path = out_dir / f"{sample_id}_if_roi_zoom_overlay_panel.png"
    plt.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close()
    logging.info(f"[{sample_id}] Saved ROI zoom-in comparison panel: {out_path} ({os.path.getsize(out_path)/1e6:.2f} MB)")
    
    art_path = ARTIFACT_DIR / f"{sample_id}_if_roi_zoom_overlay_panel.png"
    shutil.copy2(out_path, art_path)
    return str(out_path)

def render_multimodal_overview():
    """Render a comprehensive multimodal comparison figure: CytAssist Image vs TissueFAX IF vs Segmentation Mask."""
    out_dir = Path("/Volumes/KINGSTON/Evo_PRISM/results/mcseg/yousheng_vh_v113_12/overlays")
    out_dir.mkdir(parents=True, exist_ok=True)
    
    fig, axes = plt.subplots(2, 3, figsize=(20, 14), dpi=200)
    plt.subplots_adjust(wspace=0.08, hspace=0.12, left=0.03, right=0.97, top=0.94, bottom=0.03)
    
    headers = ["10x CytAssist Image (Space Ranger Reference)", "TissueFAX IF Whole Slide (DAPI + Texas Red)", "MCseg Whole-Slide Cell Mask Overlay"]
    for j, h in enumerate(headers):
        axes[0, j].set_title(h, fontsize=13, fontweight="bold", pad=10)
        
    sample_labels = [
        ("V11312_01_70mJ", "Sample 1: V11312-01 (70mJ UVB · Area A1)\n77,827 Cells Segmented"),
        ("V11312_02_20mJ", "Sample 2: V11312-02 (20mJ UVB · Area D1)\n30,494 Cells Segmented")
    ]
    
    for i, (sid, label) in enumerate(sample_labels):
        cfg = SAMPLES[sid]
        
        # Col 0: CytAssist image
        cyta = tifffile.imread(cfg["cyta_image"])
        axes[i, 0].imshow(cyta)
        axes[i, 0].set_ylabel(label, fontsize=12, fontweight="bold")
        axes[i, 0].set_xticks([])
        axes[i, 0].set_yticks([])
        
        # Col 1: IF thumbnail
        im_if = Image.open(cfg["if_image"])
        im_if_thumb = im_if.resize((cyta.shape[1], cyta.shape[0]), Image.Resampling.BILINEAR)
        axes[i, 1].imshow(np.array(im_if_thumb))
        axes[i, 1].set_xticks([])
        axes[i, 1].set_yticks([])
        
        # Col 2: Macro overlay
        macro_path = Path(cfg["out_dir"]) / f"{sid}_fullslide_macro_overlay.jpg"
        if macro_path.exists():
            im_macro = Image.open(macro_path)
            im_macro_thumb = im_macro.resize((cyta.shape[1], cyta.shape[0]), Image.Resampling.BILINEAR)
            axes[i, 2].imshow(np.array(im_macro_thumb))
        axes[i, 2].set_xticks([])
        axes[i, 2].set_yticks([])
        
    out_path = out_dir / "visium_hd_if_cytassist_multimodal_comparison.png"
    plt.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close()
    logging.info(f"Saved Multimodal Comparison Overview: {out_path}")
    
    art_path = ARTIFACT_DIR / "visium_hd_if_cytassist_multimodal_comparison.png"
    shutil.copy2(out_path, art_path)
    return str(out_path)

def main():
    results = {}
    for sid, cfg in SAMPLES.items():
        logging.info(f"\n========================================\nGenerating Overlays for: {sid}\n========================================")
        macro_p = render_fullslide_macro_overlay(sid, cfg)
        roi_p = render_roi_multi_panel(sid, cfg)
        results[sid] = {"macro": macro_p, "rois": roi_p}
        
    logging.info("\n========================================\nGenerating Multimodal Overview Figure\n========================================")
    multi_p = render_multimodal_overview()
    results["multimodal_overview"] = multi_p
    
    print("\n\nAll overlay outputs successfully rendered and copied to artifacts!")
    print(json.dumps(results, indent=2))

if __name__ == "__main__":
    main()
