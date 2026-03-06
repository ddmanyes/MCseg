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

import logging
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
from scipy.ndimage import binary_dilation

from .macenko import MacenkoNormalizer, apply_clahe

logger = logging.getLogger("pipeline.segmentation")


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


def run_eosin_watershed(nuclei_masks: np.ndarray, eosin_img: np.ndarray, bg_threshold: int, brightness_mode: bool = False):
    """Use tissue detection to expand cytoplasm via watershed.
    brightness_mode=True: tissue mask = max(R,G,B) < (255 - bg_threshold)
    brightness_mode=False (legacy): tissue mask = R-B > bg_threshold
    """
    markers = nuclei_masks.copy()
    if brightness_mode:
        if eosin_img.ndim == 3:
            brightness = eosin_img[:, :, :3].astype(np.float32).max(axis=2)
        else:
            brightness = eosin_img.astype(np.float32)
        fg = (brightness < (255 - bg_threshold)).astype(np.uint8)  # 亮度低 = 組織
        fg_pct = fg.mean() * 100
        if fg_pct < 1.0:
            logger.warning(f"Tissue mask skipped: coverage {fg_pct:.1f}% < 1%")
            return nuclei_masks
        logger.info(f"Tissue coverage: {fg_pct:.1f}%")
        # 使用亮度反轉之後做 watershed marker
        return watershed(brightness, markers, mask=fg.astype(bool))
    else:
        fg = (eosin_img > bg_threshold).astype(np.uint8)
        fg_pct = fg.mean() * 100
        if fg_pct < 1.0:
            logger.warning(f"Eosin watershed skipped: fg coverage {fg_pct:.1f}% < 1% (eosin max={eosin_img.max()}, threshold={bg_threshold})")
            return nuclei_masks
        logger.info(f"Eosin fg coverage: {fg_pct:.1f}%")
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
            # 清除所有含入的小細胞像素（包含延伸到大細胞外的殘留像素）
            for sid in small_ids_in_region:
                merged[merged == sid] = 0
            merged[region] = next_id
            next_id += 1

    return merged


def merge_enclosed_cells(masks: np.ndarray,
                         dilation_px: int = 6,
                         coverage_threshold: float = 0.55) -> np.ndarray:
    """
    將被單一較大細胞大面積包圍的細胞合併到該細胞。

    判斷方法（覆蓋率 + 四象限）：
    1. 對每個細胞做 dilation_px px dilation ring
    2. 若 ring 中同一個較大細胞佔比 >= coverage_threshold → 可能被包圍
    3. 再確認該較大細胞出現在 ring 的所有 4 個象限（排除「只在一側的鄰居」）
    4. 通過後合併

    此方法容許 Cellpose 在細胞間留下的背景間隙（通常 < 40% ring 面積），
    同時透過四象限檢查避免將相鄰而非包圍的細胞錯誤合併。
    """
    result = masks.copy()
    struct = np.ones((3, 3), dtype=bool)  # 8-connectivity

    cell_ids = np.unique(result)
    cell_ids = cell_ids[cell_ids > 0]
    sizes = {cid: int(np.sum(result == cid)) for cid in cell_ids}
    cell_ids_sorted = sorted(cell_ids, key=lambda cid: sizes[cid])

    n_merged = 0
    for cid in cell_ids_sorted:
        cell_mask = (result == cid)
        if not cell_mask.any():
            continue  # 已被前一輪合併

        # ── dilation ring ──────────────────────────────────────────────────
        dilated = binary_dilation(cell_mask, structure=struct, iterations=dilation_px)
        ring = dilated & ~cell_mask
        ring_vals = result[ring]
        if ring_vals.size == 0:
            continue

        # ── 找出 ring 中佔比最高的細胞（可含背景）──────────────────────────
        unique_all, counts_all = np.unique(ring_vals, return_counts=True)
        id_to_count = dict(zip(unique_all.tolist(), counts_all.tolist()))

        # 只考慮比自己大的細胞
        cell_counts = {u: c for u, c in id_to_count.items()
                       if u != 0 and sizes.get(int(u), 0) > sizes.get(cid, 0)}
        if not cell_counts:
            continue

        parent = int(max(cell_counts, key=cell_counts.get))
        coverage = cell_counts[parent] / ring_vals.size  # 含背景的覆蓋率
        if coverage < coverage_threshold:
            continue

        # ── 四象限檢查：parent 必須出現在細胞質心的四個象限 ─────────────────
        ys, xs = np.where(cell_mask)
        cy, cx = float(ys.mean()), float(xs.mean())
        parent_ring = np.zeros_like(ring, dtype=bool)
        parent_ring[ring] = (ring_vals == parent)
        py, px = np.where(parent_ring)
        if py.size == 0:
            continue
        q1 = bool(np.any((py < cy) & (px < cx)))   # top-left
        q2 = bool(np.any((py < cy) & (px >= cx)))  # top-right
        q3 = bool(np.any((py >= cy) & (px < cx)))  # bottom-left
        q4 = bool(np.any((py >= cy) & (px >= cx))) # bottom-right
        if not (q1 and q2 and q3 and q4):
            continue  # parent 未出現在所有象限 → 只是鄰居，非包圍

        result[cell_mask] = parent
        n_merged += 1
        logger.debug(f"Enclosed cell {cid}({sizes[cid]}px) → {parent}({sizes[parent]}px) cov={coverage:.2f}")

    if n_merged:
        logger.info(f"merge_enclosed_cells: {n_merged} cells merged")
    return result


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

    logger.info(f"Loading image: {input_path}")
    tif = tifffile.TiffFile(input_path)
    page0 = tif.pages[0]
    img_full = page0.asarray()

    # Standardize image shape/channels
    is_grayscale = False
    if img_full.ndim == 2:
        is_grayscale = True
        logger.info("Input image is grayscale")
    elif img_full.ndim == 3:
        if img_full.shape[-1] == 4:
            img_full = img_full[..., :3]
        if img_full.shape[-1] == 1:
            img_full = img_full[..., 0]
            is_grayscale = True
            logger.info("Input image is grayscale (1-channel)")

    if is_grayscale:
        normalize_stains = False
        logger.info("Skipping stain normalization for grayscale image")

    h_full, w_full = img_full.shape[:2]
    logger.info(f"Image size: {w_full} x {h_full}")

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
            logger.info("Macenko calibration successful")
        else:
            logger.warning("Macenko calibration failed, using grayscale fallback")

    tif.close()

    # Cellpose Model
    model = models.CellposeModel(gpu=use_gpu, pretrained_model=model_type)
    logger.info(f"Cellpose model: {model_type}, GPU: {use_gpu}")

    # Tiling
    final_masks = np.zeros((h_full, w_full), dtype=np.int32)
    global_id = 0

    ny = max(1, (h_full - overlap) // (block_size - overlap) + 1)
    nx = max(1, (w_full - overlap) // (block_size - overlap) + 1)
    total_tiles = ny * nx

    logger.info(f"Processing {total_tiles} tiles ({nx}x{ny})")

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

    logger.info(f"Total cells: {global_id}")

    # 合併被完全封閉在其他細胞內的子細胞
    if pp_config.get("enable_merge_enclosed", True):
        logger.info("Running merge_enclosed_cells...")
        final_masks = merge_enclosed_cells(final_masks)

    # ⚠️ 注意：Eosin Watershed 不修改 final_masks（分割結果應保持純 LOGIC_A 輸出）
    # cyto_mask.npy 另外在下方獨立計算，供 Proseg 使用

    # Save
    npy_path = os.path.join(output_dir, mask_filename)
    tif_path = os.path.join(output_dir, mask_tif_filename)

    np.save(npy_path, final_masks)
    logger.info(f"Saved: {npy_path}")

    tifffile.imwrite(tif_path, final_masks.astype(np.uint16), compression='zlib')
    logger.info(f"Saved: {tif_path}")

    if save_flows:
        logger.info("Flow saving skipped in pipeline mode (available in original script)")

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
    logger.info(f"找到 {n} 個 ROI 待分割")

    for i, (roi_name, he_crop_path) in enumerate(roi_paths):
        if progress_callback:
            progress_callback(i / n, f"ROI {i+1}/{n}: {roi_name}")
        logger.info("=" * 50)
        logger.info(f"處理 ROI: {roi_name} ({i+1}/{n})")
        _run_single_roi_segmentation(he_crop_path, roi_name, seg_cfg)

    logger.info(f"所有 {n} 個 ROI 分割完成")


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
    clahe_clip_limit  = float(prep_config.get("clahe_clip_limit", 2.0))
    mask_filename     = out_config.get("mask_filename", "masks.npy")
    mask_tif_filename = out_config.get("mask_tif_filename", "masks.tif")

    output_dir = he_crop_path.parent

    logger.info(f"Loading: {he_crop_path}")
    img_full = tifffile.imread(str(he_crop_path))
    if img_full.ndim == 3 and img_full.shape[-1] == 4:
        img_full = img_full[..., :3]

    is_grayscale = img_full.ndim == 2 or (img_full.ndim == 3 and img_full.shape[-1] == 1)
    if is_grayscale:
        normalize_stains = False
        if img_full.ndim == 3:
            img_full = img_full[..., 0]

    h_full, w_full = img_full.shape[:2]
    logger.info(f"Image size: {w_full} x {h_full}")

    normalizer = MacenkoNormalizer()
    if normalize_stains:
        success = normalizer.fit(img_full)
        if success:
            logger.info("Macenko calibration successful")
        else:
            logger.warning("Macenko fallback to grayscale")

    model = models.CellposeModel(gpu=use_gpu, pretrained_model=model_type)
    logger.info(f"Model: {model_type}, GPU: {use_gpu}")

    final_masks = np.zeros((h_full, w_full), dtype=np.int32)
    global_id   = 0
    flow_canvas = None  # (H, W, 3) uint8

    # 小於 block_size 的影像直接以 1 tile 處理
    if h_full <= block_size and w_full <= block_size:
        ny, nx = 1, 1
    else:
        ny = max(1, (h_full - overlap) // (block_size - overlap) + 1)
        nx = max(1, (w_full - overlap) // (block_size - overlap) + 1)

    logger.info(f"Processing {ny * nx} tiles ({nx}x{ny})")

    for iy in tqdm(range(ny), desc=roi_name):
        for ix in range(nx):
            y0 = iy * (block_size - overlap) if ny > 1 else 0
            x0 = ix * (block_size - overlap) if nx > 1 else 0
            y1 = min(y0 + block_size, h_full)
            x1 = min(x0 + block_size, w_full)

            tile = img_full[y0:y1, x0:x1]

            if normalize_stains and normalizer.stain_matrix is not None:
                H, E = normalizer.extract_he_channels(tile)
            elif tile.ndim == 3 and tile.shape[-1] >= 3:
                H = cv2.cvtColor(tile[..., :3], cv2.COLOR_RGB2GRAY)
                E = H  # 無法分離時退化為灰階
            else:
                H = tile.squeeze()
                E = H

            H = apply_clahe(H, clip_limit=clahe_clip_limit)
            E = apply_clahe(E, clip_limit=clahe_clip_limit)

            # 根據模型類型選擇通道策略：
            #   nuclei          → [H, H, H]（僅核）
            #   cyto / cyto2 / cyto3 → [E, H, 0]（細胞質 + 核雙通道）
            if model_type in ("cyto", "cyto2", "cyto3"):
                zeros = np.zeros_like(H)
                input_img = np.stack([E, H, zeros], axis=-1)
            else:
                input_img = np.stack([H, H, H], axis=-1)

            masks_s, flows_s, _ = model.eval(input_img, diameter=dia_small,
                flow_threshold=flow_thresh, cellprob_threshold=cellprob_thresh,
                batch_size=batch_size)
            masks_l, _, _ = model.eval(input_img, diameter=dia_large,
                flow_threshold=flow_thresh, cellprob_threshold=cellprob_thresh,
                batch_size=batch_size)

            # 將 dP flows 拼接到全圖 canvas (flows_s[0] is (H, W, 3) RGB uint8 or float)
            if flows_s and len(flows_s) > 0:
                dp = flows_s[0]  
                if flow_canvas is None:
                    flow_canvas = np.zeros((h_full, w_full, 3), dtype=np.uint8)
                if dp.ndim == 3 and dp.shape[-1] == 3:
                    dp_u8 = np.clip(dp, 0, 255).astype(np.uint8)
                    flow_canvas[y0:y1, x0:x1, :] = dp_u8[:y1-y0, :x1-x0, :]

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

    logger.info(f"Total cells: {global_id}")

    # 合併被完全封閉在其他細胞內的子細胞
    if pp_config.get("enable_merge_enclosed", True):
        logger.info("Running merge_enclosed_cells...")
        final_masks = merge_enclosed_cells(final_masks)

    # ── Flow 視覺化（小尺寸 dP），用於分割品質檢查 ----------------------------------
    if flow_canvas is not None:
        try:
            from PIL import Image as _Image
            import io as _io
            flow_rgb = flow_canvas
            flow_u8 = _Image.fromarray(flow_rgb)
            flow_buf = _io.BytesIO()
            flow_u8.save(flow_buf, "JPEG", quality=85)
            flows_preview_path = output_dir / "flows_preview.jpg"
            flows_preview_path.write_bytes(flow_buf.getvalue())
            logger.info(f"Saved Flow Preview: {flows_preview_path}")
        except Exception as e:
            logger.warning(f"Flow visualization failed: {e}")

    if pp_config.get("enable_eosin_watershed", False) and not is_grayscale:
        logger.info("Generating Eosin Cytoplasm Mask (Brightness Method)...")
        bg_thresh = pp_config.get("eosin_bg_threshold", 40)
        # 亮度法：白色空直背景被排除，組織區域保留
        brightness = img_full[:, :, :3].astype(np.float32).max(axis=2)
        is_background = (brightness > (255 - bg_thresh))
        # 產生純細胞質遮罩 (Cyto Mask) 供 Proseg 約束使用 (0=背景, 1=組織)
        cyto_mask = (~is_background).astype(np.int32)
        cyto_npy_path = output_dir / "cyto_mask.npy"
        np.save(str(cyto_npy_path), cyto_mask)
        logger.info(f"Saved Cyto Mask: {cyto_npy_path}")

    npy_path = output_dir / mask_filename
    tif_path = output_dir / mask_tif_filename

    np.save(str(npy_path), final_masks)
    logger.info(f"Saved: {npy_path}")

    tifffile.imwrite(str(tif_path), final_masks.astype(np.uint16), compression='zlib')
    logger.info(f"Saved: {tif_path}")

    return final_masks
