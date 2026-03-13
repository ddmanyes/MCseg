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


def reconcile_stitched_labels(masks: np.ndarray, min_flat_length: int = 4) -> np.ndarray:
    """
    偵測並修補 Cellpose 拼縫偽影（直線邊界）。
    如果兩個相鄰細胞在拼接線上共享一段水平或垂直的直線邊界，則將其合併。
    """
    from collections import defaultdict
    result = masks.copy()
    
    # 統計所有在相鄰像素點但 ID 不同的對
    # (pair, depth) -> [indices along the boundary]
    h_adj = defaultdict(list) # ((l1, l2), row) -> [cols...]
    v_adj = defaultdict(list) # ((l1, l2), col) -> [rows...]
    
    # 1. 水平邊界掃描 (標記物件 L1 在 L2 上方且接觸)
    diff_h = (result[:-1, :] != result[1:, :]) & (result[:-1, :] > 0) & (result[1:, :] > 0)
    rows_h, cols_h = np.where(diff_h)
    for r, c in zip(rows_h, cols_h):
        l1, l2 = result[r, c], result[r+1, c]
        pair = tuple(sorted((l1, l2)))
        h_adj[(pair, r)].append(c)
            
    # 2. 垂直邊界掃描
    diff_v = (result[:, :-1] != result[:, 1:]) & (result[:, :-1] > 0) & (result[:, 1:] > 0)
    rows_v, cols_v = np.where(diff_v)
    for r, c in zip(rows_v, cols_v):
        l1, l2 = result[r, c], result[r, c+1]
        pair = tuple(sorted((l1, l2)))
        v_adj[(pair, c)].append(r)

    merge_map = {}
    def find_root(i):
        root = i
        while root in merge_map:
            root = merge_map[root]
        # Path compression
        curr = i
        while curr in merge_map:
            next_p = merge_map[curr]
            merge_map[curr] = root
            curr = next_p
        return root

    # 3. 判定與合併：如果某對 label 在同一條像素線上連續接觸超過閾值
    # 這代表它們是被 Tiling Grid 硬生生切斷的
    for (pair, r), cols in h_adj.items():
        if len(cols) >= min_flat_length:
            cols.sort()
            # 檢查是否為連續直線段
            if cols[-1] - cols[0] + 1 == len(cols):
                u, v = find_root(pair[0]), find_root(pair[1])
                if u != v: merge_map[u] = v

    for (pair, c), rows in v_adj.items():
        if len(rows) >= min_flat_length:
            rows.sort()
            if rows[-1] - rows[0] + 1 == len(rows):
                u, v = find_root(pair[0]), find_root(pair[1])
                if u != v: merge_map[u] = v

    if not merge_map:
        return result
        
    logger.info(f"Detected tiling artifacts: reconciling {len(merge_map)} cell fragments...")
    
    # 4. 套用合併映射 (優化：使用 lookup table 避免效能瓶頸)
    id_map = {cid: find_root(cid) for cid in np.unique(result) if cid > 0}
    max_id = int(result.max())
    lookup = np.arange(max_id + 1, dtype=result.dtype)
    for src, dst in id_map.items():
        lookup[src] = dst
        
    result = lookup[result]
    return result


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
            logger.warning("Macenko calibration failed (possible lack of tissue or poor contrast), falling back to grayscale.")

    tif.close()

    # Cellpose Model
    model = models.CellposeModel(gpu=use_gpu, pretrained_model=model_type)
    logger.info(f"Cellpose model: {model_type}, GPU: {use_gpu}")

    # ==========================
    # 執行分割 (內建分塊拼接 Tiling)
    # ==========================
    # 優點：Cellpose 內建 Tiling 會在 Flow/Probability 階段進行平均融合，徹底消除拼縫直線
    
    # 1. 全域預處理
    if normalize_stains and normalizer.stain_matrix is not None:
        logger.info("Applying Macenko Color Deconvolution (Hematoxylin)...")
        gray_full = normalizer.extract_hematoxylin(img_full)
    elif img_full.ndim == 3 and img_full.shape[-1] >= 3:
        gray_full = cv2.cvtColor(img_full[..., :3], cv2.COLOR_RGB2GRAY)
    else:
        gray_full = img_full.squeeze()
        
    gray_full = apply_clahe(gray_full)
    input_full = np.stack([gray_full, gray_full, gray_full], axis=-1)
    
    # 2. 雙直徑評估 (Dual Inference)
    # 智能決定是否需要分塊：若 ROI 尺寸小於 block_size，則強制不切圖以避免接縫
    # [優化] 預設 block_size 提高到 3072，處理中型 ROI 更穩健
    effective_block_size = max(block_size, 3072)
    do_tile = (h_full > effective_block_size) or (w_full > effective_block_size)
    tile_overlap_ratio = 0.25 if do_tile else 0.1 # 調高重疊率以利縫合
    
    eval_kwargs = {
        "channels":           [0, 0], # 強制 Cellpose 識別為 2D 影像
        "do_3D":              False,
        "flow_threshold":     flow_thresh,
        "cellprob_threshold": cellprob_thresh,
        "batch_size":         batch_size,
        "bsize":              256, # 恢復為 Cellpose 4.0 預設的 256 (Transformer 模型硬性要求)
        "tile_overlap":       tile_overlap_ratio,
        "resample":           True,
        "stitch_threshold":   0.0, # 設為 0 以避免誤啟動 3D 縫合邏輯導致 z_axis 報錯
    }

    if do_tile:
        logger.info(f"Running Cellpose Tiling (overlap={tile_overlap_ratio:.2f}, stitch_thresh=0.5)...")
    else:
        logger.info("ROI is small enough, running without tiling to ensure no artifacts.")

    logger.info(f"Step 1/2: Evaluating small diameter ({dia_small})...")
    masks_s, _, _ = model.eval(input_full, diameter=dia_small, **eval_kwargs)
    
    logger.info(f"Step 2/2: Evaluating large diameter ({dia_large})...")
    masks_l, _, _ = model.eval(input_full, diameter=dia_large, **eval_kwargs)
    
    # 資源釋放：input_full 為 3-通道疊加，佔用大量記憶體。
    # 既然預測已完成，主動將其清空以利後續合併運算。
    input_full = None
    import gc
    gc.collect()

    # 3. 合併雙尺寸結果
    logger.info("Merging small and large diameter masks...")
    final_masks = _merge_masks_logic_a(masks_s, masks_l, frag_thresh)
    
    # [新增] 拼縫修補：將被 Tiling 切斷的細胞片段重新融合
    logger.info("Running reconcile_stitched_labels to fix tiling artifacts...")
    final_masks = reconcile_stitched_labels(final_masks)
    
    global_id = int(final_masks.max())
    logger.info(f"Stitching finished. Total cells: {global_id}")

    # 合併被完全封閉在其他細胞內的子細胞
    if pp_config.get("enable_merge_enclosed", True):
        logger.info("Running merge_enclosed_cells...")
        final_masks = merge_enclosed_cells(final_masks)


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

# 可被 ROI 覆寫的欄位 → (seg_cfg section, key)
_ROI_OVERRIDE_FIELD_MAP: dict[str, tuple[str, str]] = {
    "model_type":         ("cellpose_model", "model_type"),
    "dia_small":          ("strategy",       "dia_small"),
    "dia_large":          ("strategy",       "dia_large"),
    "flow_threshold":     ("strategy",       "flow_threshold"),
    "cellprob_threshold": ("strategy",       "cellprob_threshold"),
    "fragment_threshold": ("strategy",       "fragment_threshold"),
    "block_size":         ("tiling",         "block_size"),
    "overlap":            ("tiling",         "overlap"),
}


def _merge_roi_params(global_seg_cfg: dict, roi_overrides: dict) -> dict:
    """將 ROI 特定覆寫合併進全域分割設定（深複製，不修改原始設定）。"""
    import copy
    cfg = copy.deepcopy(global_seg_cfg)
    for key, (section, field) in _ROI_OVERRIDE_FIELD_MAP.items():
        if key in roi_overrides and roi_overrides[key] is not None:
            cfg.setdefault(section, {})[field] = roi_overrides[key]
    return cfg


def run_segmentation_rois(config: dict, progress_callback=None,
                          roi_overrides: dict | None = None,
                          target_roi: str | None = None):
    """對所有（或指定）ROI 的 he_crop.tif 執行分割，結果分別存至各 ROI 目錄。

    Args:
        roi_overrides: {roi_name: {field: value}} 形式的 ROI 個別參數覆寫。
                       未指定的欄位沿用全域 config 設定。
        target_roi: 若指定，只重跑此 ROI（單 ROI 重做模式）。
    """
    paths      = config.get("paths", {})
    output_dir = paths.get("output_dir", "results/analysis")
    rois       = config.get("rois", [])
    seg_cfg    = config.get("segmentation", {})
    roi_base   = Path(output_dir) / "roi"
    overrides  = roi_overrides or {}

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

    # 單 ROI 重做模式：過濾只留指定的 ROI
    if target_roi:
        filtered = [(name, path) for name, path in roi_paths if name == target_roi]
        if not filtered:
            raise ValueError(f"找不到 ROI '{target_roi}' 的 he_crop.tif")
        roi_paths = filtered
        logger.info(f"單 ROI 重做模式：只處理 {target_roi}")

    n = len(roi_paths)
    logger.info(f"找到 {n} 個 ROI 待分割")

    for i, (roi_name, he_crop_path) in enumerate(roi_paths):
        if progress_callback:
            progress_callback(i / n, f"ROI {i+1}/{n}: {roi_name}")
        logger.info("=" * 50)
        logger.info(f"處理 ROI: {roi_name} ({i+1}/{n})")

        roi_specific = overrides.get(roi_name, {})
        effective_seg_cfg = _merge_roi_params(seg_cfg, roi_specific) if roi_specific else seg_cfg
        if roi_specific:
            logger.info(f"  套用個別參數覆寫：{roi_specific}")

        _run_single_roi_segmentation(he_crop_path, roi_name, effective_seg_cfg)

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

    # ==========================
    # 執行分割 (原生 Tiling + 雙尺寸合併)
    # ==========================
    # 1. 影像預處理 (全圖)
    if normalize_stains and normalizer.stain_matrix is not None:
        H, E = normalizer.extract_he_channels(img_full)
    elif img_full.ndim == 3 and img_full.shape[-1] >= 3:
        H = cv2.cvtColor(img_full[..., :3], cv2.COLOR_RGB2GRAY)
        E = H
    else:
        H = img_full.squeeze()
        E = H

    H = apply_clahe(H, clip_limit=clahe_clip_limit)
    E = apply_clahe(E, clip_limit=clahe_clip_limit)

    # 選擇通道策略
    if model_type in ("cyto", "cyto2", "cyto3"):
        zeros = np.zeros_like(H)
        input_img = np.stack([E, H, zeros], axis=-1)
    else:
        input_img = np.stack([H, H, H], axis=-1)

    # 2. 評估參數
    effective_block_size = max(block_size, 3072)
    do_tile = (h_full > effective_block_size) or (w_full > effective_block_size)
    tile_overlap_ratio = 0.25 if do_tile else 0.1
    
    eval_kwargs = {
        "channels":           [0, 0],
        "do_3D":              False,
        "flow_threshold":     flow_thresh,
        "cellprob_threshold": cellprob_thresh,
        "batch_size":         batch_size,
        "bsize":              256,
        "tile_overlap":       tile_overlap_ratio,
        "resample":           True,
        "stitch_threshold":   0.0,
    }

    logger.info(f"Evaluating small diameter ({dia_small})...")
    masks_s, flows_s, _ = model.eval(input_img, diameter=dia_small, **eval_kwargs)
    
    logger.info(f"Evaluating large diameter ({dia_large})...")
    masks_l, _, _ = model.eval(input_img, diameter=dia_large, **eval_kwargs)

    # 用於儲存 Flow 預覽
    flow_canvas = None
    if flows_s and len(flows_s) > 0:
        dp = flows_s[0]
        if dp.ndim == 3 and dp.shape[-1] == 3:
            flow_canvas = np.clip(dp, 0, 255).astype(np.uint8)

    # 資源中途清理
    input_img = None
    import gc; gc.collect()

    # 3. 合併與修正
    logger.info("Merging masks and reconciling tiling artifacts...")
    final_masks = _merge_masks_logic_a(masks_s, masks_l, frag_thresh)
    final_masks = reconcile_stitched_labels(final_masks)

    global_id = int(final_masks.max())
    logger.info(f"Segmentation finished. Total cells: {global_id}")

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

    if pp_config.get("enable_eosin_watershed", True) and not is_grayscale:
        # 使用多尺度判斷：高亮度 為背景
        bg_thresh = pp_config.get("eosin_bg_threshold", 50)
        brightness = img_full[:, :, :3].astype(np.float32).max(axis=2)
        is_background = (brightness > (255 - bg_thresh))
        
        # 產生組織遮罩 (Cyto Mask)
        cyto_mask_raw = (~is_background).astype(np.uint8)
        
        # [極度放寬] 暫停形態學操作，避免削減邊緣
        cyto_mask = (~is_background).astype(np.int32)
        cyto_npy_path = output_dir / "segmentation_masks_cyto.npy"
        np.save(str(cyto_npy_path), cyto_mask)
        logger.info(f"Saved Cleaned Cyto Mask for Proseg: {cyto_npy_path}")
        # [優化] 釋放組織遮罩記憶體
        cyto_mask = None
        is_background = None
        brightness = None

    npy_path = output_dir / mask_filename
    tif_path = output_dir / mask_tif_filename

    np.save(str(npy_path), final_masks)
    logger.info(f"Saved: {npy_path}")

    tifffile.imwrite(str(tif_path), final_masks.astype(np.uint16), compression='zlib')
    logger.info(f"Saved: {tif_path}")

    # [優化] 在回傳前徹底清理大型 Array
    import gc
    final_masks = None
    gc.collect()

    return None # 改為回傳 None，因為呼叫端通常直接讀磁碟，減少記憶體對象傳遞
