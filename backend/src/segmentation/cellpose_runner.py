# ==========================================================
# Provenance: cell_segmentation/scripts/run_segmentation_v3.py
# Original Project: cellpose-he-segmentation
# Extracted: 2026-02-20
# Note: 提取核心函數為可匯入模組，配置改由 pipeline.yaml 驅動
# ==========================================================
"""
Cellpose 全圖分割執行器

核心功能：
- Logic A 雙尺寸合併策略 (Small + Large)
- Macenko 色彩標準化 → Hematoxylin 提取
- Eosin Watershed 細胞質擴張
- 分塊拼接 (Tiling + Stitching)
"""

import os
import sys
import numpy as np
import cv2
import tifffile
import torch
from tqdm import tqdm
from cellpose import models, core
from pathlib import Path
from skimage.segmentation import watershed

from .macenko import MacenkoNormalizer, apply_clahe


def get_best_tissue_patch(page_obj, step=50, grid_n=8):
    """
    Finds the densest tissue crop using a low-res thumbnail.
    Used for auto-calibrating Macenko normalization.
    """
    thumb_level = len(page_obj.pages) - 1
    thumb = page_obj.pages[thumb_level].asarray()

    # Handle multi-channel images (RGB, RGBA)
    if thumb.ndim == 3:
        if thumb.shape[-1] == 4:
            thumb = thumb[..., :3]
        if thumb.shape[-1] == 3:
            gray = cv2.cvtColor(thumb, cv2.COLOR_RGB2GRAY)
        else:
            # Fallback for unexpected number of channels
            gray = thumb[..., 0]
    else:
        # Already grayscale
        gray = thumb

    tissue = (gray < 220).astype(np.float32)

    h_t, w_t = tissue.shape
    best_score = -1
    best_rc = (0, 0)

    ph = h_t // grid_n
    pw = w_t // grid_n

    for r in range(0, h_t - ph, step):
        for c in range(0, w_t - pw, step):
            score = tissue[r:r+ph, c:c+pw].sum()
            if score > best_score:
                best_score = score
                best_rc = (r, c)

    scale = page_obj.pages[0].shape[0] / h_t
    ry, rx = best_rc
    full_y = int(ry * scale)
    full_x = int(rx * scale)
    full_h = int(ph * scale)
    full_w = int(pw * scale)

    return full_y, full_x, full_h, full_w


def run_eosin_watershed(nuclei_masks: np.ndarray, eosin_img: np.ndarray, bg_threshold: int):
    """Use Eosin channel watershed to expand cytoplasm."""
    markers = nuclei_masks.copy()
    fg = (eosin_img > bg_threshold).astype(np.uint8)
    fg_pct = fg.mean() * 100
    if fg_pct < 1.0:
        print(f"⚠️  Eosin watershed skipped: fg coverage {fg_pct:.1f}% < 1% (eosin max={eosin_img.max()}, threshold={bg_threshold})")
        return nuclei_masks
    print(f"  Eosin fg coverage: {fg_pct:.1f}%")
    return watershed(-eosin_img, markers, mask=fg)


def _merge_masks_logic_a(masks_small, masks_large, fragment_threshold=50):
    """
    Strategy LOGIC_A:
    1. If large cell covers > 1 small cell -> Keep small cells (Under-segmentation)
    2. If large cell covers <= 1 small cell -> Keep large cell (Over-segmentation)
    """
    merged = masks_small.copy()

    large_ids = np.unique(masks_large)
    large_ids = large_ids[large_ids > 0]

    next_id = masks_small.max() + 1

    for lid in large_ids:
        region = (masks_large == lid)
        small_ids_in_region = np.unique(masks_small[region])
        small_ids_in_region = small_ids_in_region[small_ids_in_region > 0]

        significant_small = []
        for sid in small_ids_in_region:
            if np.sum(masks_small == sid) >= fragment_threshold:
                significant_small.append(sid)

        if len(significant_small) <= 1:
            merged[region] = next_id
            next_id += 1

    return merged


def run_segmentation(config: dict):
    """
    執行全圖分割的主函數。
    """
    # 支援傳入完整配置或僅傳入 segmentation 區塊
    if "segmentation" in config and "paths" in config:
        # 傳入的是完整配置
        paths = config["paths"]
        input_path = paths.get("he_image")
        output_dir = paths.get("masks_dir", "results/masks")
        config = config["segmentation"].copy()
    else:
        # 傳入的是 segmentation 區塊
        input_path = config.get("input_path")
        output_dir = config.get("output_dir", "results/masks")

    if input_path is None:
        raise ValueError("未指定輸入影像路徑 (he_image or input_path is None)")

    model_config = config.get("cellpose_model", {})
    strategy = config.get("strategy", {})
    pp_config = config.get("postprocessing", {})
    tile_config = config.get("tiling", {})
    prep_config = config.get("preprocessing", {})
    out_config = config.get("output", {})

    model_type = model_config.get("model_type", "nuclei")
    use_gpu = model_config.get("use_gpu", True) and core.use_gpu()
    batch_size = model_config.get("batch_size", 16)

    dia_small = strategy.get("dia_small", 30.0)
    dia_large = strategy.get("dia_large", 60.0)
    flow_thresh = strategy.get("flow_threshold", 0.4)
    cellprob_thresh = strategy.get("cellprob_threshold", -1.0)
    frag_thresh = strategy.get("fragment_threshold", 200)

    block_size = tile_config.get("block_size", 2048)
    overlap = tile_config.get("overlap", 256)

    normalize_stains = prep_config.get("normalize_stains", True)

    save_flows = out_config.get("save_flows", True)
    mask_filename = out_config.get("mask_filename", "segmentation_masks.npy")
    mask_tif_filename = out_config.get("mask_tif_filename", "segmentation_masks.tif")

    os.makedirs(output_dir, exist_ok=True)

    print(f"🔬 Loading image: {input_path}")
    tif = tifffile.TiffFile(input_path)
    page0 = tif.pages[0]
    img_full = page0.asarray()

    # Standardize image shape/channels
    is_grayscale = False
    if img_full.ndim == 2:
        is_grayscale = True
        print("  Info: Input image is grayscale")
    elif img_full.ndim == 3:
        if img_full.shape[-1] == 4:
            img_full = img_full[..., :3]
        if img_full.shape[-1] == 1:
            img_full = img_full[..., 0]
            is_grayscale = True
            print("  Info: Input image is grayscale (1-channel)")

    if is_grayscale:
        normalize_stains = False
        print("  Info: Skipping stain normalization for grayscale image")

    h_full, w_full = img_full.shape[:2]
    print(f"  Image size: {w_full} x {h_full}")

    # Macenko 校正
    normalizer = MacenkoNormalizer()
    if normalize_stains:
        calib_roi = prep_config.get("calibration_roi")
        if calib_roi:
            cy, cx, ch, cw = calib_roi
            calib_patch = img_full[cy:cy+ch, cx:cx+cw]
        else:
            fy, fx, fh, fw = get_best_tissue_patch(tif)
            calib_patch = img_full[fy:fy+fh, fx:fx+fw]

        success = normalizer.fit(calib_patch)
        if success:
            print("✅ Macenko calibration successful")
        else:
            print("⚠️ Macenko calibration failed, using grayscale fallback")

    # Cellpose Model
    model = models.CellposeModel(gpu=use_gpu, pretrained_model=model_type)
    print(f"🧠 Cellpose model: {model_type}, GPU: {use_gpu}")

    # Tiling
    final_masks = np.zeros((h_full, w_full), dtype=np.int32)
    global_id = 0

    ny = max(1, (h_full - overlap) // (block_size - overlap) + 1)
    nx = max(1, (w_full - overlap) // (block_size - overlap) + 1)
    total_tiles = ny * nx

    print(f"🧩 Processing {total_tiles} tiles ({nx}x{ny})")

    for iy in tqdm(range(ny), desc="Rows"):
        for ix in range(nx):
            y0 = iy * (block_size - overlap)
            x0 = ix * (block_size - overlap)
            y1 = min(y0 + block_size, h_full)
            x1 = min(x0 + block_size, w_full)

            tile = img_full[y0:y1, x0:x1]

            # Preprocess
            if normalize_stains and normalizer.stain_matrix is not None:
                gray = normalizer.extract_hematoxylin(tile)
            elif tile.ndim == 3 and tile.shape[-1] >= 3:
                gray = cv2.cvtColor(tile[..., :3], cv2.COLOR_RGB2GRAY)
            else:
                # Already grayscale or single channel
                gray = tile.squeeze()

            gray = apply_clahe(gray)
            input_img = np.stack([gray, gray, gray], axis=-1)

            # Dual inference
            masks_s, _, _ = model.eval(
                input_img, diameter=dia_small,
                flow_threshold=flow_thresh,
                cellprob_threshold=cellprob_thresh,
                batch_size=batch_size
            )

            masks_l, _, _ = model.eval(
                input_img, diameter=dia_large,
                flow_threshold=flow_thresh,
                cellprob_threshold=cellprob_thresh,
                batch_size=batch_size
            )

            # Merge
            merged = _merge_masks_logic_a(masks_s, masks_l, frag_thresh)

            # Stitch into global
            inner_y0 = overlap // 2 if iy > 0 else 0
            inner_x0 = overlap // 2 if ix > 0 else 0
            inner_y1 = merged.shape[0] - (overlap // 2 if iy < ny - 1 else 0)
            inner_x1 = merged.shape[1] - (overlap // 2 if ix < nx - 1 else 0)

            inner_mask = merged[inner_y0:inner_y1, inner_x0:inner_x1]

            ids = np.unique(inner_mask)
            ids = ids[ids > 0]

            for old_id in ids:
                global_id += 1
                final_masks[y0+inner_y0:y0+inner_y1, x0+inner_x0:x0+inner_x1][inner_mask == old_id] = global_id

    print(f"🔢 Total cells: {global_id}")

    # Eosin Expansion
    if pp_config.get("enable_eosin_watershed", False):
        print("🔬 Running Eosin Watershed expansion...")
        # Approximate Eosin from RGB
        if is_grayscale:
            # For grayscale H&E, tissue is DARK (low value), background is LIGHT (high value).
            # We need to invert it so that tissue becomes the 'peak' for watershed.
            eosin_approx = 255.0 - img_full.astype(np.float32)
        else:
            eosin_approx = img_full[:, :, 0].astype(np.float32) - img_full[:, :, 2].astype(np.float32)

        eosin_approx = np.clip(eosin_approx, 0, 255).astype(np.uint8)
        bg_thresh = pp_config.get("eosin_bg_threshold", 40)
        final_masks = run_eosin_watershed(final_masks, eosin_approx, bg_thresh)

    # Save
    npy_path = os.path.join(output_dir, mask_filename)
    tif_path = os.path.join(output_dir, mask_tif_filename)

    np.save(npy_path, final_masks)
    print(f"💾 Saved: {npy_path}")

    tifffile.imwrite(tif_path, final_masks.astype(np.uint16), compression='zlib')
    print(f"💾 Saved: {tif_path}")

    if save_flows:
        print("ℹ️ Flow saving skipped in pipeline mode (available in original script)")

    return final_masks


# ── Per-ROI segmentation (Stage 1 主要執行路徑) ──────────────────────────────


def run_segmentation_rois(config: dict, progress_callback=None):
    """對所有 ROI 的 he_crop.tif 執行分割，結果分別存至各 ROI 目錄。"""
    paths      = config.get("paths", {})
    output_dir = paths.get("output_dir", "results/analysis")
    rois       = config.get("rois", [])
    seg_cfg    = config.get("segmentation", {})
    roi_base   = Path(output_dir) / "roi"

    # 按 config rois 順序收集
    roi_paths: list[tuple[str, Path]] = []
    for roi in rois:
        roi_name = roi.get("name", "")
        he_crop  = roi_base / roi_name / "he_crop.tif"
        if he_crop.exists():
            roi_paths.append((roi_name, he_crop))

    # 掃描目錄補充未在 config 的 ROI
    known = {r[0] for r in roi_paths}
    if roi_base.exists():
        for d in sorted(roi_base.iterdir()):
            if d.is_dir() and d.name not in known:
                he_crop = d / "he_crop.tif"
                if he_crop.exists():
                    roi_paths.append((d.name, he_crop))

    if not roi_paths:
        raise ValueError("找不到 he_crop.tif，請先在 Stage 0 執行 ROI 裁切")

    n = len(roi_paths)
    print(f"🗂️ 找到 {n} 個 ROI 待分割")

    for i, (roi_name, he_crop_path) in enumerate(roi_paths):
        if progress_callback:
            progress_callback(i / n, f"ROI {i+1}/{n}: {roi_name}")
        print(f"\n{'='*50}")
        print(f"📂 處理 ROI: {roi_name} ({i+1}/{n})")
        _run_single_roi_segmentation(he_crop_path, roi_name, seg_cfg)

    print(f"\n✅ 所有 {n} 個 ROI 分割完成")


def _run_single_roi_segmentation(he_crop_path: Path, roi_name: str, seg_cfg: dict):
    """對單一 he_crop.tif 執行 Cellpose 分割，結果存至同目錄。"""
    model_config      = seg_cfg.get("cellpose_model", {})
    strategy          = seg_cfg.get("strategy", {})
    pp_config         = seg_cfg.get("postprocessing", {})
    tile_config       = seg_cfg.get("tiling", {})
    prep_config       = seg_cfg.get("preprocessing", {})
    out_config        = seg_cfg.get("output", {})

    model_type        = model_config.get("model_type", "cyto2")
    use_gpu           = model_config.get("use_gpu", True) and core.use_gpu()
    batch_size        = model_config.get("batch_size", 4)
    dia_small         = strategy.get("dia_small", 30.0)
    dia_large         = strategy.get("dia_large", 60.0)
    flow_thresh       = strategy.get("flow_threshold", 0.4)
    cellprob_thresh   = strategy.get("cellprob_threshold", -1.0)
    frag_thresh       = strategy.get("fragment_threshold", 200)
    block_size        = tile_config.get("block_size", 2048)
    overlap           = tile_config.get("overlap", 256)
    normalize_stains  = prep_config.get("normalize_stains", True)
    mask_filename     = out_config.get("mask_filename", "masks.npy")
    mask_tif_filename = out_config.get("mask_tif_filename", "masks.tif")

    output_dir = he_crop_path.parent

    print(f"🔬 Loading: {he_crop_path}")
    img_full = tifffile.imread(str(he_crop_path))
    if img_full.ndim == 3 and img_full.shape[-1] == 4:
        img_full = img_full[..., :3]

    is_grayscale = img_full.ndim == 2 or (img_full.ndim == 3 and img_full.shape[-1] == 1)
    if is_grayscale:
        normalize_stains = False
        if img_full.ndim == 3:
            img_full = img_full[..., 0]

    h_full, w_full = img_full.shape[:2]
    print(f"  Image size: {w_full} x {h_full}")

    normalizer = MacenkoNormalizer()
    if normalize_stains:
        success = normalizer.fit(img_full)
        print("✅ Macenko calibration successful" if success else "⚠️ Macenko fallback to grayscale")

    model = models.CellposeModel(gpu=use_gpu, pretrained_model=model_type)
    print(f"🧠 Model: {model_type}, GPU: {use_gpu}")

    final_masks = np.zeros((h_full, w_full), dtype=np.int32)
    global_id   = 0

    # 小於 block_size 的影像直接以 1 tile 處理
    if h_full <= block_size and w_full <= block_size:
        ny, nx = 1, 1
    else:
        ny = max(1, (h_full - overlap) // (block_size - overlap) + 1)
        nx = max(1, (w_full - overlap) // (block_size - overlap) + 1)

    print(f"🧩 Processing {ny * nx} tiles ({nx}x{ny})")

    for iy in tqdm(range(ny), desc=roi_name):
        for ix in range(nx):
            y0 = iy * (block_size - overlap) if ny > 1 else 0
            x0 = ix * (block_size - overlap) if nx > 1 else 0
            y1 = min(y0 + block_size, h_full)
            x1 = min(x0 + block_size, w_full)

            tile = img_full[y0:y1, x0:x1]

            if normalize_stains and normalizer.stain_matrix is not None:
                gray = normalizer.extract_hematoxylin(tile)
            elif tile.ndim == 3 and tile.shape[-1] >= 3:
                gray = cv2.cvtColor(tile[..., :3], cv2.COLOR_RGB2GRAY)
            else:
                gray = tile.squeeze()

            gray = apply_clahe(gray)
            input_img = np.stack([gray, gray, gray], axis=-1)

            masks_s, _, _ = model.eval(input_img, diameter=dia_small,
                flow_threshold=flow_thresh, cellprob_threshold=cellprob_thresh,
                batch_size=batch_size)
            masks_l, _, _ = model.eval(input_img, diameter=dia_large,
                flow_threshold=flow_thresh, cellprob_threshold=cellprob_thresh,
                batch_size=batch_size)

            merged = _merge_masks_logic_a(masks_s, masks_l, frag_thresh)

            inner_y0 = overlap // 2 if iy > 0 else 0
            inner_x0 = overlap // 2 if ix > 0 else 0
            inner_y1 = merged.shape[0] - (overlap // 2 if iy < ny - 1 else 0)
            inner_x1 = merged.shape[1] - (overlap // 2 if ix < nx - 1 else 0)

            inner_mask = merged[inner_y0:inner_y1, inner_x0:inner_x1]
            ids = np.unique(inner_mask)
            ids = ids[ids > 0]

            for old_id in ids:
                global_id += 1
                final_masks[
                    y0+inner_y0:y0+inner_y1,
                    x0+inner_x0:x0+inner_x1
                ][inner_mask == old_id] = global_id

    print(f"🔢 Total cells: {global_id}")

    if pp_config.get("enable_eosin_watershed", False) and not is_grayscale:
        print("🔬 Running Eosin Watershed expansion...")
        eosin_approx = img_full[:, :, 0].astype(np.float32) - img_full[:, :, 2].astype(np.float32)
        eosin_approx = np.clip(eosin_approx, 0, 255).astype(np.uint8)
        bg_thresh = pp_config.get("eosin_bg_threshold", 40)
        final_masks = run_eosin_watershed(final_masks, eosin_approx, bg_thresh)

    npy_path = output_dir / mask_filename
    tif_path = output_dir / mask_tif_filename

    np.save(str(npy_path), final_masks)
    print(f"💾 Saved: {npy_path}")

    tifffile.imwrite(str(tif_path), final_masks.astype(np.uint16), compression='zlib')
    print(f"💾 Saved: {tif_path}")

    return final_masks
