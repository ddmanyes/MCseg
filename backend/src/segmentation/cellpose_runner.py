# ==========================================================
# MSseg — MCseg v2 分割引擎
# 基於 autoresearch_seg/segment_best.py（V12 Voronoi 集成）
# 改編自 visiumHD_pipeline_3，移除 Proseg，整合 MCseg v2
# ==========================================================
"""
MCseg v2 多模型集成分割器

核心演算法（V12）：
  1. 預處理：CLAHE + 組織遮罩 + Ruifrok H&E 色彩分離
  2. 多模型推論：cyto3 × 3 直徑 + 可選 hematoxylin pass + 可選 cpsam × 3
  3. 非重疊合併（merge_masks_fast）
  4. 轉錄本密度補救（可選，需 vhd_csv）
  5. Voronoi 擴張（防止重疊）
  6. 清理 + 重新編號

保留 pipeline_3 的架構：
  - per-ROI 分割（run_segmentation_rois）
  - ROI 個別參數覆寫（roi_overrides）
  - 單 ROI 重做模式（target_roi）
"""

from __future__ import annotations

import gc
import logging
import time
from pathlib import Path

import cv2
import numpy as np
import tifffile
from scipy import ndimage

logger = logging.getLogger("pipeline.segmentation")


def _clear_gpu_cache() -> None:
    """釋放 CUDA / MPS 顯存碎片（CLAUDE.md §13）。"""
    try:
        import torch
        if torch.cuda.is_available():
            torch.cuda.empty_cache()
        elif hasattr(torch.backends, "mps") and torch.backends.mps.is_available():
            torch.mps.empty_cache()
    except ImportError:
        pass


# ─────────────────────────────────────────────────────────
# 預處理工具
# ─────────────────────────────────────────────────────────

def apply_clahe(img: np.ndarray, clip_limit: float = 3.0,
                tile_size: int = 8) -> np.ndarray:
    """對 RGB 或灰階影像套用 CLAHE。"""
    if img.ndim == 3 and img.shape[-1] >= 3:
        lab = cv2.cvtColor(img[..., :3], cv2.COLOR_RGB2LAB)
        l_ch, a_ch, b_ch = cv2.split(lab)
        clahe = cv2.createCLAHE(clipLimit=clip_limit,
                                tileGridSize=(tile_size, tile_size))
        return cv2.cvtColor(cv2.merge((clahe.apply(l_ch), a_ch, b_ch)),
                            cv2.COLOR_LAB2RGB)
    clahe = cv2.createCLAHE(clipLimit=clip_limit,
                            tileGridSize=(tile_size, tile_size))
    return clahe.apply(img.astype(np.uint8))


def create_tissue_mask(img: np.ndarray) -> np.ndarray:
    """從 H&E 影像建立組織遮罩（排除白色背景）。"""
    gray = (cv2.cvtColor(img[..., :3], cv2.COLOR_RGB2GRAY)
            if img.ndim == 3 else img)
    tissue = (gray < 220).astype(np.uint8)
    kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (11, 11))
    tissue = cv2.morphologyEx(tissue, cv2.MORPH_CLOSE, kernel)
    kernel2 = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (5, 5))
    tissue = cv2.dilate(tissue, kernel2, iterations=2)
    return tissue.astype(bool)


def color_deconvolution_he(img: np.ndarray) -> np.ndarray:
    """Ruifrok & Johnston H&E 色彩分離，提取 Hematoxylin 通道（uint8）。"""
    img_f = img[..., :3].astype(np.float64) + 1.0
    od = -np.log(img_f / 256.0)

    he_matrix = np.array([
        [0.6500286, 0.7041680, 0.2860126],   # Hematoxylin
        [0.0728940, 0.9904310, 0.1155140],   # Eosin
        [0.2688350, 0.5706770, 0.7768750],   # DAB (residual)
    ], dtype=np.float64)
    for i in range(3):
        norm = np.linalg.norm(he_matrix[i])
        if norm > 0:
            he_matrix[i] /= norm

    stains = (od.reshape(-1, 3) @ np.linalg.inv(he_matrix).T
              ).reshape(img.shape[:2] + (3,))
    hema = np.clip(stains[:, :, 0], 0, None)
    h_max = np.percentile(hema, 99.5)
    if h_max > 0:
        hema = np.clip(hema / h_max, 0, 1)

    hema_u8 = (hema * 255).astype(np.uint8)
    clahe = cv2.createCLAHE(clipLimit=2.0, tileGridSize=(8, 8))
    return clahe.apply(hema_u8)


# ─────────────────────────────────────────────────────────
# 後處理工具
# ─────────────────────────────────────────────────────────

def voronoi_expand(mask: np.ndarray, max_distance: int,
                   tissue_mask: np.ndarray | None = None) -> np.ndarray:
    """
    Voronoi 約束擴張：每個背景像素分配給最近細胞，距離上限 max_distance px。
    優於 expand_labels：不產生重疊，細胞自然填滿可用空間。
    """
    binary = mask > 0
    if not binary.any():
        return mask.copy()

    dist, nearest_idx = ndimage.distance_transform_edt(~binary,
                                                        return_indices=True)
    expanded = mask[nearest_idx[0], nearest_idx[1]]
    expanded[dist > max_distance] = 0
    if tissue_mask is not None:
        expanded[~tissue_mask] = 0
    return expanded.astype(np.int32)


def merge_masks_fast(base_mask: np.ndarray, new_mask: np.ndarray,
                     max_overlap_ratio: float = 0.15,
                     min_size: int = 15) -> tuple[np.ndarray, int]:
    """
    將 new_mask 中不重疊的細胞併入 base_mask（原地修改）。
    回傳 (base_mask, n_added)。

    使用 regionprops 取 bounding box，在小 patch 上做 overlap 計算，
    避免舊版 O(n_cells × H×W) 全圖布林掃描。
    """
    from skimage.measure import regionprops

    next_id = int(base_mask.max()) + 1
    added = 0
    base_occupied = base_mask > 0

    for prop in regionprops(new_mask):
        pixel_count = prop.area
        if pixel_count < min_size:
            continue
        r0, c0, r1, c1 = prop.bbox
        nid = prop.label
        crop_new  = new_mask[r0:r1, c0:c1] == nid
        crop_base = base_occupied[r0:r1, c0:c1]
        overlap = int((crop_base & crop_new).sum())
        if overlap / pixel_count < max_overlap_ratio:
            empty = crop_new & ~crop_base
            if int(empty.sum()) >= min_size:
                rows_e, cols_e = np.where(empty)
                base_mask[rows_e + r0, cols_e + c0] = next_id
                base_occupied[rows_e + r0, cols_e + c0] = True
                next_id += 1
                added += 1

    return base_mask, added


def clean_mask(mask: np.ndarray, min_size: int = 20,
               max_size: int = 6000) -> np.ndarray:
    """移除面積過小或過大的細胞（原地修改）。"""
    labels_arr, counts_arr = np.unique(mask, return_counts=True)
    remove = labels_arr[
        ((counts_arr < min_size) | (counts_arr > max_size)) & (labels_arr > 0)
    ]
    if len(remove) > 0:
        mask[np.isin(mask, remove)] = 0
    return mask


def relabel_sequential(mask: np.ndarray) -> np.ndarray:
    """重新編號為從 1 開始的連續 ID（LUT 向量化，O(max_id)）。"""
    unique_labels = np.unique(mask)
    unique_labels = unique_labels[unique_labels > 0]
    if len(unique_labels) == 0:
        return mask
    lut = np.zeros(int(mask.max()) + 1, dtype=mask.dtype)
    lut[unique_labels] = np.arange(1, len(unique_labels) + 1, dtype=mask.dtype)
    return lut[mask]


# ─────────────────────────────────────────────────────────
# 轉錄本密度補救（可選）
# ─────────────────────────────────────────────────────────

def find_transcript_seeds(
    vhd_csv: str,
    img_shape: tuple[int, int, int],
    existing_mask: np.ndarray,
    tissue_mask: np.ndarray,
) -> tuple[np.ndarray, int]:
    """
    從 Visium HD 轉錄本密度尋找 Cellpose 遺漏的細胞位置。
    vhd_csv 需含 'x', 'y' 欄位（影像像素座標）。
    回傳 (seed_mask, n_added)。
    """
    import pandas as pd

    try:
        df = pd.read_csv(vhd_csv)
    except Exception:
        return np.zeros(img_shape[:2], dtype=np.int32), 0

    h, w = img_shape[:2]
    density = np.zeros((h, w), dtype=np.float32)
    x_c = np.clip(df["x"].values.astype(int), 0, w - 1)
    y_c = np.clip(df["y"].values.astype(int), 0, h - 1)
    np.add.at(density, (y_c, x_c), 1)

    density_smooth = cv2.GaussianBlur(density, (0, 0), sigmaX=5.0)
    density_max = ndimage.maximum_filter(density_smooth, size=15)
    local_max = (density_smooth == density_max) & (density_smooth > 2.0)

    exist_dil = cv2.dilate(
        (existing_mask > 0).astype(np.uint8),
        cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (11, 11)),
    )
    local_max = local_max & (exist_dil == 0) & tissue_mask

    peak_labels, n_peaks = ndimage.label(local_max)
    seed_mask = np.zeros((h, w), dtype=np.int32)
    next_id = int(existing_mask.max()) + 1
    added = 0

    if n_peaks > 0:
        centroids = ndimage.center_of_mass(
            local_max, peak_labels, range(1, n_peaks + 1)
        )
        for cy, cx in centroids:
            cy_i, cx_i = int(round(cy)), int(round(cx))
            if 0 <= cy_i < h and 0 <= cx_i < w:
                region = existing_mask[
                    max(0, cy_i - 5) : min(h, cy_i + 6),
                    max(0, cx_i - 5) : min(w, cx_i + 6),
                ]
                if (region > 0).sum() / max(region.size, 1) < 0.1:
                    cv2.circle(seed_mask, (cx_i, cy_i), 5, int(next_id), -1)
                    next_id += 1
                    added += 1

    seed_mask[~tissue_mask] = 0
    return seed_mask, added


# ─────────────────────────────────────────────────────────
# MCseg v2 核心：多模型集成
# ─────────────────────────────────────────────────────────

def run_mcseg_v2(
    img: np.ndarray,
    cfg: dict,
    vhd_csv: str | None = None,
) -> np.ndarray:
    """
    MCseg v2 主分割函數（V12 Voronoi 集成）。

    Args:
        img:     RGB H&E 影像 (H, W, 3)
        cfg:     mcseg_v2 設定 dict（來自 config segmentation.mcseg_v2）
        vhd_csv: 轉錄本密度 CSV 路徑（可選；不存在時靜默跳過）

    Returns:
        int32 細胞分割遮罩 (H, W)
    """
    from cellpose import core, models

    t0 = time.time()

    # ── 參數 ────────────────────────────────────────────────
    use_gpu          = bool(cfg.get("use_gpu", True)) and core.use_gpu()
    batch_size       = int(cfg.get("batch_size", 4))
    dia_small        = float(cfg.get("dia_small", 13.0))
    dia_mid          = float(cfg.get("dia_mid", 17.0))
    dia_large        = float(cfg.get("dia_large", 22.0))
    use_hematoxylin  = bool(cfg.get("use_hematoxylin", True))
    use_cpsam        = bool(cfg.get("use_cpsam", False))
    voronoi_dist     = int(cfg.get("voronoi_distance", 9))
    flow_thresh      = float(cfg.get("flow_threshold", 0.4))
    cellprob_thresh  = float(cfg.get("cellprob_threshold", -2.0))
    min_size         = int(cfg.get("min_size", 20))
    max_size         = int(cfg.get("max_size", 6000))
    use_rescue       = bool(cfg.get("use_transcript_rescue", True))
    clahe_clip       = float(cfg.get("clahe_clip_limit", 3.0))
    # ── cpsam 獨立直徑與 cellprob（論文 7-pass 規格）──────────────────────
    # dia_cpsam_auto=0 → Cellpose 自動偵測直徑（論文 Pass 5/7 30 auto）
    # dia_cpsam_small=16 → 論文 Pass 6 固定 16px
    dia_cpsam_auto   = float(cfg.get("dia_cpsam_auto",   0.0))   # 0 = auto
    dia_cpsam_small  = float(cfg.get("dia_cpsam_small",  16.0))  # fixed 16px
    cellprob_cpsam_auto  = float(cfg.get("cellprob_cpsam_auto",  -1.0))
    cellprob_cpsam_small = float(cfg.get("cellprob_cpsam_small", -3.0))
    cellprob_cpsam_hema  = float(cfg.get("cellprob_cpsam_hema",  -1.0))

    logger.info("[MCseg v2] === V12 Voronoi 集成分割 ===")
    logger.info(
        f"  cyto3 dia: {dia_small}/{dia_mid}/{dia_large} | "
        f"cpsam={use_cpsam} (dia_auto={dia_cpsam_auto or 'auto'}, dia_small={dia_cpsam_small}) | "
        f"voronoi_d={voronoi_dist}"
    )

    # ── 1. 預處理 ────────────────────────────────────────────
    enhanced = apply_clahe(img, clip_limit=clahe_clip, tile_size=8)
    tissue_mask = create_tissue_mask(img)

    hema: np.ndarray | None = None
    if use_hematoxylin:
        hema = color_deconvolution_he(img)

    # ── 2. 多模型推論 ────────────────────────────────────────
    results: dict[str, np.ndarray] = {}

    eval_base = dict(
        channels=[0, 0],
        flow_threshold=flow_thresh,
        cellprob_threshold=cellprob_thresh,
        min_size=10,
        batch_size=batch_size,
    )

    logger.info(f"  [{time.time()-t0:.0f}s] 載入 cyto3...")
    cyto3 = models.CellposeModel(model_type="cyto3", gpu=use_gpu)

    logger.info(f"  [{time.time()-t0:.0f}s] cyto3 dia={dia_mid} (RGB)...")
    m, _, _ = cyto3.eval(enhanced, diameter=dia_mid,
                         augment=True, resample=True, **eval_base)
    results["cyto3_mid"] = m
    logger.info(f"    → {m.max()} cells")

    logger.info(f"  [{time.time()-t0:.0f}s] cyto3 dia={dia_small} (small)...")
    m, _, _ = cyto3.eval(
        enhanced, diameter=dia_small, augment=False, resample=True,
        **{**eval_base, "cellprob_threshold": cellprob_thresh - 1.0},
    )
    results["cyto3_small"] = m
    logger.info(f"    → {m.max()} cells")

    logger.info(f"  [{time.time()-t0:.0f}s] cyto3 dia={dia_large} (large)...")
    m, _, _ = cyto3.eval(
        enhanced, diameter=dia_large, augment=False, resample=True,
        **{**eval_base, "cellprob_threshold": cellprob_thresh + 1.0},
    )
    results["cyto3_large"] = m
    logger.info(f"    → {m.max()} cells")

    if use_hematoxylin and hema is not None:
        hema_rgb = np.stack([hema, hema, hema], axis=-1)
        logger.info(f"  [{time.time()-t0:.0f}s] cyto3 dia={dia_mid} (hematoxylin)...")
        m, _, _ = cyto3.eval(hema_rgb, diameter=dia_mid,
                             augment=True, resample=True, **eval_base)
        results["cyto3_hema"] = m
        logger.info(f"    → {m.max()} cells")
        del hema_rgb

    del cyto3
    gc.collect()
    _clear_gpu_cache()

    if use_cpsam:
        logger.info(f"  [{time.time()-t0:.0f}s] 載入 cpsam...")
        try:
            cpsam = models.CellposeModel(model_type="cpsam", gpu=use_gpu)
            cpsam_base = dict(
                channels=[0, 0],
                flow_threshold=flow_thresh,
                min_size=10,
                batch_size=batch_size,
                augment=False,
                resample=False,
            )

            # Pass 5（論文）：cpsam CLAHE-RGB，dia=auto（0 → Cellpose 自偵測 ~30px），cellprob=-1.0
            _dia_auto = dia_cpsam_auto if dia_cpsam_auto > 0 else None  # None = Cellpose auto
            logger.info(f"  [{time.time()-t0:.0f}s] cpsam Pass5 dia={'auto' if _dia_auto is None else _dia_auto} (RGB, cellprob={cellprob_cpsam_auto})...")
            m, _, _ = cpsam.eval(
                enhanced, diameter=_dia_auto,
                **{**cpsam_base, "cellprob_threshold": cellprob_cpsam_auto},
            )
            results["cpsam_auto"] = m
            logger.info(f"    → {m.max()} cells")

            # Pass 6（論文）：cpsam CLAHE-RGB，dia=16px，cellprob=-3.0
            logger.info(f"  [{time.time()-t0:.0f}s] cpsam Pass6 dia={dia_cpsam_small} (RGB, cellprob={cellprob_cpsam_small})...")
            m, _, _ = cpsam.eval(
                enhanced, diameter=float(dia_cpsam_small),
                **{**cpsam_base, "cellprob_threshold": cellprob_cpsam_small},
            )
            results["cpsam_small"] = m
            logger.info(f"    → {m.max()} cells")

            # Pass 7（論文）：cpsam Hematoxylin，dia=auto，cellprob=-1.0
            if use_hematoxylin and hema is not None:
                hema_rgb2 = np.stack([hema, hema, hema], axis=-1)
                logger.info(f"  [{time.time()-t0:.0f}s] cpsam Pass7 dia={'auto' if _dia_auto is None else _dia_auto} (hema, cellprob={cellprob_cpsam_hema})...")
                m, _, _ = cpsam.eval(
                    hema_rgb2, diameter=_dia_auto,
                    **{**cpsam_base, "cellprob_threshold": cellprob_cpsam_hema},
                )
                results["cpsam_hema"] = m
                logger.info(f"    → {m.max()} cells")
                del hema_rgb2

            del cpsam
            gc.collect()
            _clear_gpu_cache()
        except Exception as e:
            logger.warning(f"  cpsam 失敗（跳過）：{e}")

    del hema
    gc.collect()

    logger.info(f"  [{time.time()-t0:.0f}s] 所有模型完成，開始合併...")

    # ── 3. 集成合併（以 cyto3_mid 為基底）────────────────────
    base_mask = results["cyto3_mid"].copy().astype(np.int32)
    target_shape = base_mask.shape
    for key, mask in results.items():
        if key == "cyto3_mid":
            continue
        m = mask.astype(np.int32)
        if m.shape != target_shape:
            import cv2 as _cv2
            m = _cv2.resize(m, (target_shape[1], target_shape[0]),
                            interpolation=_cv2.INTER_NEAREST)
        base_mask, n_added = merge_masks_fast(base_mask, m)
        logger.info(f"    合併 {key}: +{n_added} cells")
    logger.info(f"  合併後：{base_mask.max()} cells")

    # ── 4. 轉錄本密度補救（可選）────────────────────────────
    if use_rescue and vhd_csv and Path(vhd_csv).exists():
        logger.info(f"  [{time.time()-t0:.0f}s] 轉錄本密度補救...")
        rescue_mask, n_rescued = find_transcript_seeds(
            vhd_csv, img.shape, base_mask, tissue_mask
        )
        if n_rescued > 0:
            next_rescue_id = int(base_mask.max()) + 1
            for rid in np.unique(rescue_mask)[1:]:
                rpix = rescue_mask == rid
                if (base_mask[rpix] > 0).sum() == 0:
                    base_mask[rpix] = next_rescue_id
                    next_rescue_id += 1
            logger.info(f"  補救 {n_rescued} cells")
    elif use_rescue and vhd_csv:
        logger.info(f"  vhd_csv 不存在，跳過轉錄本補救：{vhd_csv}")

    # ── 5. 清理 + Voronoi 擴張 ───────────────────────────────
    base_mask[~tissue_mask] = 0
    base_mask = clean_mask(base_mask, min_size=min_size, max_size=max_size)
    base_mask = relabel_sequential(base_mask)
    logger.info(f"  擴張前：{base_mask.max()} cells")

    final_mask = voronoi_expand(base_mask, max_distance=voronoi_dist,
                                tissue_mask=tissue_mask)
    final_mask = clean_mask(final_mask, min_size=min_size, max_size=max_size)

    n_final = int(len(np.unique(final_mask)) - 1)
    logger.info(f"  [{time.time()-t0:.0f}s] 最終：{n_final} cells")
    return final_mask.astype(np.int32)


# ─────────────────────────────────────────────────────────
# ROI 分割（維持 pipeline_3 架構）
# ─────────────────────────────────────────────────────────

# 可被 ROI 覆寫的欄位（均屬 mcseg_v2 section）
_ROI_OVERRIDE_FIELDS: frozenset[str] = frozenset({
    "use_gpu", "batch_size",
    "dia_small", "dia_mid", "dia_large",
    "use_hematoxylin", "use_cpsam",
    "dia_cpsam_auto", "dia_cpsam_small",
    "cellprob_cpsam_auto", "cellprob_cpsam_small", "cellprob_cpsam_hema",
    "voronoi_distance",
    "flow_threshold", "cellprob_threshold",
    "min_size", "max_size",
    "use_transcript_rescue",
    "clahe_clip_limit",
})


def _merge_roi_params(global_seg_cfg: dict, roi_overrides: dict) -> dict:
    """將 ROI 個別覆寫合併進全域分割設定（深複製，不改原始 dict）。"""
    import copy
    cfg = copy.deepcopy(global_seg_cfg)
    mcseg = cfg.setdefault("mcseg_v2", {})
    for key, val in roi_overrides.items():
        if key in _ROI_OVERRIDE_FIELDS and val is not None:
            mcseg[key] = val
    return cfg


def run_segmentation_rois(
    config: dict,
    progress_callback=None,
    roi_overrides: dict | None = None,
    target_roi: str | None = None,
) -> None:
    """
    對所有（或指定）ROI 的 he_crop.tif 執行 MCseg v2 分割。

    Args:
        config:            完整 pipeline config dict
        progress_callback: fn(progress: float, message: str)
        roi_overrides:     {roi_name: {field: value}}
        target_roi:        若指定，只重跑此 ROI
    """
    paths      = config.get("paths", {})
    output_dir = paths.get("output_dir", "results/analysis")
    rois       = config.get("rois", [])
    seg_cfg    = config.get("segmentation", {})
    roi_base   = Path(output_dir) / "roi"
    overrides  = roi_overrides or {}

    # 收集 ROI 路徑（先依 config 順序，再補掃描目錄）
    roi_paths: list[tuple[str, Path]] = []
    for roi in rois:
        roi_name = roi.get("name", "")
        he_crop  = roi_base / roi_name / "he_crop.tif"
        if he_crop.exists():
            roi_paths.append((roi_name, he_crop))

    known = {r[0] for r in roi_paths}
    if roi_base.exists():
        for d in sorted(roi_base.iterdir()):
            if d.is_dir() and d.name not in known:
                he_crop = d / "he_crop.tif"
                if he_crop.exists():
                    roi_paths.append((d.name, he_crop))

    if not roi_paths:
        raise ValueError("找不到 he_crop.tif，請先在 Stage 0 執行 ROI 裁切")

    if target_roi:
        filtered = [(n, p) for n, p in roi_paths if n == target_roi]
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
        effective_cfg = (
            _merge_roi_params(seg_cfg, roi_specific) if roi_specific else seg_cfg
        )
        if roi_specific:
            logger.info(f"  套用個別參數覆寫：{roi_specific}")

        _run_single_roi(he_crop_path, roi_name, effective_cfg)

    if progress_callback:
        progress_callback(1.0, f"全部 {n} 個 ROI 分割完成")
    logger.info(f"所有 {n} 個 ROI 分割完成")


def _run_single_roi(he_crop_path: Path, _roi_name: str,
                    seg_cfg: dict) -> None:
    """對單一 he_crop.tif 執行 MCseg v2，結果存至同目錄。"""
    mcseg_cfg = seg_cfg.get("mcseg_v2", {})
    out_cfg   = seg_cfg.get("output", {})

    mask_filename     = out_cfg.get("mask_filename",     "segmentation_masks.npy")
    mask_tif_filename = out_cfg.get("mask_tif_filename", "segmentation_masks.tif")
    output_dir        = he_crop_path.parent

    logger.info(f"讀取：{he_crop_path}")
    img = tifffile.imread(str(he_crop_path))
    if img.ndim == 3 and img.shape[-1] == 4:
        img = img[..., :3]
    if img.ndim == 2:
        img = np.stack([img, img, img], axis=-1)

    logger.info(f"影像尺寸：{img.shape[1]} × {img.shape[0]}")

    # 轉錄本密度 CSV（若存在則啟用補救）
    vhd_csv = str(output_dir / "vhd_pseudo_transcripts.csv")

    final_masks = run_mcseg_v2(img, mcseg_cfg, vhd_csv=vhd_csv)

    npy_path = output_dir / mask_filename
    tif_path = output_dir / mask_tif_filename

    np.save(str(npy_path), final_masks)
    logger.info(f"儲存：{npy_path}")

    tif_dtype = np.uint16 if final_masks.max() <= 65535 else np.uint32
    if tif_dtype == np.uint32:
        logger.warning(f"Cell count {final_masks.max()} exceeds uint16 range — saving TIF as uint32")
    tifffile.imwrite(str(tif_path), final_masks.astype(tif_dtype),
                     compression="zlib")
    logger.info(f"儲存：{tif_path}")

    del final_masks, img
    gc.collect()


# ─────────────────────────────────────────────────────────
# 全圖 Tiled 分割（MPS 安全版）
# ─────────────────────────────────────────────────────────

def run_tiled_mcseg_v2(
    img: np.ndarray,
    cfg: dict,
    tile_size: int = 1024,
    overlap: int = 128,
    progress_callback=None,
) -> np.ndarray:
    """
    全圖 tiled MCseg v2 分割。

    兩階段設計（避免 Voronoi 在 tile 邊界產生接縫）：
      Phase 1（per-tile）：Cellpose 多直徑推論 + merge_masks_fast
      Phase 2（全圖）：拼接 → Voronoi 擴張 → clean_mask

    MPS 安全設定：
      - tile_size=1024（預設），比 2048 佔用更少 GPU 記憶體
      - augment=False 全程停用（augment=True 在 MPS 上易 OOM）
      - batch_size 限制為 cfg 設定值，建議 ≤ 2
      - 捕捉 MPS RuntimeError 後自動 fallback 到 CPU

    Args:
        img:               (H, W, 3) uint8 RGB H&E 影像（全圖或大 ROI）
        cfg:               mcseg_v2 設定 dict
        tile_size:         每塊大小（px），MPS 安全建議 1024
        overlap:           相鄰塊重疊寬度（px）
        progress_callback: fn(progress: float, message: str)

    Returns:
        int32 細胞分割遮罩 (H, W)
    """
    from cellpose import core, models

    t0 = time.time()

    # ── 參數 ────────────────────────────────────────────────
    use_gpu         = bool(cfg.get("use_gpu", True)) and core.use_gpu()
    batch_size      = int(cfg.get("batch_size", 2))   # MPS 建議 ≤ 2
    dia_small       = float(cfg.get("dia_small", 13.0))
    dia_mid         = float(cfg.get("dia_mid", 17.0))
    dia_large       = float(cfg.get("dia_large", 22.0))
    use_hematoxylin = bool(cfg.get("use_hematoxylin", True))
    voronoi_dist    = int(cfg.get("voronoi_distance", 9))
    flow_thresh     = float(cfg.get("flow_threshold", 0.4))
    cellprob_thresh = float(cfg.get("cellprob_threshold", -2.0))
    min_size        = int(cfg.get("min_size", 20))
    max_size        = int(cfg.get("max_size", 6000))
    use_cpsam       = bool(cfg.get("use_cpsam", False))
    clahe_clip      = float(cfg.get("clahe_clip_limit", 3.0))
    # ── cpsam 獨立直徑與 cellprob（論文 7-pass 規格）──────────────────────
    dia_cpsam_auto       = float(cfg.get("dia_cpsam_auto",       0.0))   # 0 = auto
    dia_cpsam_small      = float(cfg.get("dia_cpsam_small",     16.0))
    cellprob_cpsam_auto  = float(cfg.get("cellprob_cpsam_auto",  -1.0))
    cellprob_cpsam_small = float(cfg.get("cellprob_cpsam_small", -3.0))
    cellprob_cpsam_hema  = float(cfg.get("cellprob_cpsam_hema",  -1.0))

    H, W = img.shape[:2]
    y_starts = list(range(0, H, tile_size))
    x_starts = list(range(0, W, tile_size))
    total_tiles = len(y_starts) * len(x_starts)

    logger.info(
        f"[Tiled MCseg v2] 全圖 {W}×{H}px  "
        f"tile={tile_size}px overlap={overlap}px  "
        f"tiles={total_tiles}  gpu={use_gpu}"
    )

    # ── 預處理（全圖一次性，節省重複計算）────────────────────
    enhanced     = apply_clahe(img, clip_limit=clahe_clip, tile_size=8)
    tissue_mask  = create_tissue_mask(img)
    hema_full: np.ndarray | None = None
    if use_hematoxylin:
        hema_full = color_deconvolution_he(img)

    # ── Phase 1：per-tile Cellpose ────────────────────────
    stitched = np.zeros((H, W), dtype=np.int32)
    current_max = 0

    logger.info(f"  載入 cyto3 模型 (gpu={use_gpu})")
    cyto3 = models.CellposeModel(model_type="cyto3", gpu=use_gpu)

    cpsam = None
    if use_cpsam:
        logger.info(f"  載入 cpsam 模型 (gpu={use_gpu})")
        try:
            cpsam = models.CellposeModel(model_type="cpsam", gpu=use_gpu)
            logger.info("  cpsam 載入成功")
        except Exception as e:
            logger.warning(f"  cpsam 載入失敗（跳過）：{e}")
            cpsam = None

    eval_base = dict(
        channels=[0, 0],
        flow_threshold=flow_thresh,
        cellprob_threshold=cellprob_thresh,
        min_size=10,
        batch_size=batch_size,
        augment=False,   # MPS 安全：停用 augment
        resample=True,
    )

    for ti, y in enumerate(y_starts):
        y0e = y - overlap if y > 0 else 0
        y1e = min(y + tile_size + overlap, H)

        for tj, x in enumerate(x_starts):
            tile_idx = ti * len(x_starts) + tj + 1
            x0e = x - overlap if x > 0 else 0
            x1e = min(x + tile_size + overlap, W)

            enh_tile  = enhanced[y0e:y1e, x0e:x1e]
            hema_tile = hema_full[y0e:y1e, x0e:x1e] if hema_full is not None else None

            msg = f"Tile {tile_idx}/{total_tiles} ({x},{y})"
            if progress_callback:
                progress_callback(tile_idx / total_tiles * 0.85, msg)
            logger.info(f"  [{time.time()-t0:.0f}s] {msg}")

            # per-tile 多直徑推論 + merge（不做 Voronoi）
            tile_results: dict[str, np.ndarray] = {}
            try:
                m, _, _ = cyto3.eval(enh_tile, diameter=dia_mid, **eval_base)
                tile_results["mid"] = m

                m, _, _ = cyto3.eval(
                    enh_tile, diameter=dia_small,
                    **{**eval_base, "cellprob_threshold": cellprob_thresh - 1.0},
                )
                tile_results["small"] = m

                m, _, _ = cyto3.eval(
                    enh_tile, diameter=dia_large,
                    **{**eval_base, "cellprob_threshold": cellprob_thresh + 1.0},
                )
                tile_results["large"] = m

                if hema_tile is not None:
                    hema_rgb = np.stack([hema_tile] * 3, axis=-1)
                    m, _, _ = cyto3.eval(hema_rgb, diameter=dia_mid, **eval_base)
                    tile_results["hema"] = m

                if cpsam is not None:
                    _dia_auto = dia_cpsam_auto if dia_cpsam_auto > 0 else None
                    cpsam_base = {**eval_base, "augment": False, "resample": False}
                    # Pass 5: cpsam RGB dia=auto, cellprob_cpsam_auto
                    m, _, _ = cpsam.eval(
                        enh_tile, diameter=_dia_auto,
                        **{**cpsam_base, "cellprob_threshold": cellprob_cpsam_auto},
                    )
                    tile_results["cpsam_auto"] = m
                    # Pass 6: cpsam RGB dia=16, cellprob_cpsam_small
                    m, _, _ = cpsam.eval(
                        enh_tile, diameter=float(dia_cpsam_small),
                        **{**cpsam_base, "cellprob_threshold": cellprob_cpsam_small},
                    )
                    tile_results["cpsam_small"] = m
                    # Pass 7: cpsam Hema dia=auto, cellprob_cpsam_hema
                    if hema_tile is not None:
                        hema_rgb2 = np.stack([hema_tile] * 3, axis=-1)
                        m, _, _ = cpsam.eval(
                            hema_rgb2, diameter=_dia_auto,
                            **{**cpsam_base, "cellprob_threshold": cellprob_cpsam_hema},
                        )
                        tile_results["cpsam_hema"] = m

            except RuntimeError as e:
                if "MPS" in str(e) or "out of memory" in str(e).lower():
                    logger.warning(f"  MPS OOM on tile {tile_idx}，fallback CPU")
                    gc.collect()
                    cpu_model = models.CellposeModel(model_type="cyto3", gpu=False)
                    m, _, _ = cpu_model.eval(enh_tile, diameter=dia_mid,
                                             **{**eval_base, "batch_size": 1})
                    tile_results["mid"] = m
                    del cpu_model
                else:
                    raise

            # merge tile results
            base = tile_results.get("mid", np.zeros_like(enh_tile[:, :, 0])).copy().astype(np.int32)
            target_h, target_w = base.shape
            for key, mask in tile_results.items():
                if key == "mid":
                    continue
                m = mask.astype(np.int32)
                if m.shape != (target_h, target_w):
                    m = cv2.resize(m, (target_w, target_h), interpolation=cv2.INTER_NEAREST)
                base, _ = merge_masks_fast(base, m)


            # 裁掉 overlap，只保留有效區域
            act_top   = y - y0e
            act_bot   = y1e - (y + tile_size)
            act_left  = x - x0e
            act_right = x1e - (x + tile_size)
            v0 = act_top
            v1 = base.shape[0] - act_bot if act_bot > 0 else base.shape[0]
            u0 = act_left
            u1 = base.shape[1] - act_right if act_right > 0 else base.shape[1]
            valid = base[v0:v1, u0:u1].copy()

            # ID offset 避免衝突
            valid[valid > 0] += current_max

            # 邊界 ID 對齊（上方 + 左方）
            mappings: dict[int, int] = {}
            if y > 0:
                prev_row = stitched[y - 1, x:x + valid.shape[1]]
                curr_row = valid[0, :len(prev_row)]
                mm = (prev_row > 0) & (curr_row > 0)
                for p, c in zip(prev_row[mm], curr_row[mm]):
                    mappings.setdefault(int(c), int(p))
            if x > 0:
                prev_col = stitched[y:y + valid.shape[0], x - 1]
                curr_col = valid[:len(prev_col), 0]
                mm = (prev_col > 0) & (curr_col > 0)
                for p, c in zip(prev_col[mm], curr_col[mm]):
                    mappings.setdefault(int(c), int(p))
            for c_lbl, p_lbl in mappings.items():
                valid[valid == c_lbl] = p_lbl

            prev_max = current_max
            current_max = max(current_max, int(valid.max()))
            tw = min(x + tile_size, W) - x
            th = min(y + tile_size, H) - y
            stitched[y:y + th, x:x + tw] = valid[:th, :tw]

            # 新增細胞 = ID 超過上一塊 max 者（邊界合併的細胞已被重映射到舊 ID）
            n_new_tile = int(np.unique(valid[valid > prev_max]).size)
            logger.info(f"    cells in tile (new): {n_new_tile}  total: {current_max}")

    del cyto3
    if cpsam is not None:
        del cpsam
    del enhanced
    if hema_full is not None:
        del hema_full
    gc.collect()
    _clear_gpu_cache()

    # ── Phase 2：全圖 Voronoi + 清理 ─────────────────────
    if progress_callback:
        progress_callback(0.90, "Voronoi 擴張（全圖）...")
    logger.info(f"  [{time.time()-t0:.0f}s] Phase 2：清理 + Voronoi 擴張")

    stitched[~tissue_mask] = 0
    stitched = clean_mask(stitched, min_size=min_size, max_size=max_size)
    stitched = relabel_sequential(stitched)
    logger.info(f"  擴張前：{stitched.max()} cells")

    final = voronoi_expand(stitched, max_distance=voronoi_dist, tissue_mask=tissue_mask)
    final = clean_mask(final, min_size=min_size, max_size=max_size)
    n_final = int(len(np.unique(final)) - 1)

    if progress_callback:
        progress_callback(1.0, f"完成：{n_final:,} 個細胞")
    logger.info(f"  [{time.time()-t0:.0f}s] 全圖分割完成：{n_final:,} cells")

    return final.astype(np.int32)


# ─────────────────────────────────────────────────────────
# Preview 用：單 patch 快速分割
# ─────────────────────────────────────────────────────────

def run_preview_patch(img_patch: np.ndarray, mcseg_cfg: dict) -> np.ndarray:
    """
    對小 patch 執行 MCseg v2 快速預覽。
    自動停用 cpsam 與 transcript rescue 以加快速度。
    """
    cfg = dict(mcseg_cfg)
    cfg["use_transcript_rescue"] = False
    cfg["use_cpsam"] = False
    return run_mcseg_v2(img_patch, cfg, vhd_csv=None)
