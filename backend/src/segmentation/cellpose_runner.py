# ==========================================================
# MSseg — MCseg v2 分割引擎
# 基於 autoresearch_seg/segment_best.py（V12 Voronoi 集成）
# 改編自 visiumHD_pipeline_3，移除 Proseg，整合 MCseg v2
# ==========================================================
"""
MCseg v2 多模型集成分割器

核心演算法（V12）：
  1. 預處理：CLAHE + 組織遮罩 + Ruifrok H&E 色彩分離
  2. 多模型推論：主模型 × 3 直徑 + 可選 hematoxylin pass + 可選 cpsam × 3
  3. 非重疊合併（merge_masks_fast）
  4. 轉錄本密度補救（可選，需 vhd_csv）
  5. Voronoi 擴張（防止重疊）
  6. 清理 + 重新編號

⚠️ 「cyto3」在本環境**實際上是 cpsam**
---------------------------------------
歷史文件把前 3-4 個 pass 稱為「cyto3 多直徑集成」，但 `cellpose 4.0.8` 已移除
`model_type` 參數（會印 `model_type argument is not used in v4.0.1+`）並一律載入
`pretrained_model` 預設值 `cpsam`。本機 `~/.cellpose/models/` 只有 `cpsam` 與
`cpsam_v2`，**cyto3 權重從未存在**。

`diameter` 仍然有效（用於把影像縮放到模型的 30px 細胞徑），所以「多直徑集成」
的機制沒有失效 —— 失準的只有模型名稱。既有結果全部是 cpsam 產生的。

`_load_primary_model` 會在載入後記錄實際權重路徑，避免這件事再次被文件蓋掉。

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


def _load_primary_model(use_gpu: bool):
    """
    載入多直徑集成的主模型，並記錄**實際**載入的權重。

    歷史程式碼寫 `model_type="cyto3"`，但 cellpose 4.0.1+ 已停用該參數（只會印
    一行 warning 然後忽略），實際載入的一律是 `pretrained_model` 預設的 `cpsam`。
    這裡不再傳那個無效參數，並把權重路徑記進 log —— 讓「跑的是哪個模型」
    永遠可從執行紀錄查證，而不是靠文件宣稱。
    """
    from cellpose import models

    model = models.CellposeModel(gpu=use_gpu)
    logger.info(f"  主模型權重：{getattr(model, 'pretrained_model', '未知')}")
    return model


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

    # Stain vectors are rows: OD = C @ M, so C = OD @ inv(M).
    stains = (od.reshape(-1, 3) @ np.linalg.inv(he_matrix)
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
# 分塊處理工具（全片規模：遮罩可能落在 memmap，不可整份進 RAM）
# ─────────────────────────────────────────────────────────

DEFAULT_BLOCK = 4096
# margin 必須 ≥ voronoi_distance；256 是其最大合理值（9）的 28 倍，留足餘裕
DEFAULT_MARGIN = 256


def iter_blocks(h: int, w: int, block: int = DEFAULT_BLOCK):
    """逐塊產生 `(y0, x0, y1, x1)`。"""
    for y0 in range(0, h, block):
        for x0 in range(0, w, block):
            yield y0, x0, min(y0 + block, h), min(x0 + block, w)


def global_label_sizes(labels, block: int = DEFAULT_BLOCK) -> np.ndarray:
    """
    逐塊累加各 label 的像素數，回傳 `counts[label]`。

    跨塊的細胞會被正確加總 —— 這是分塊處理仍能做**全域**面積過濾的前提。
    """
    h, w = labels.shape
    counts = np.zeros(1, dtype=np.int64)
    for y0, x0, y1, x1 in iter_blocks(h, w, block):
        c = np.bincount(np.asarray(labels[y0:y1, x0:x1]).ravel())
        if len(c) > len(counts):
            grown = np.zeros(len(c), dtype=np.int64)
            grown[: len(counts)] = counts
            counts = grown
        counts[: len(c)] += c
    return counts


def apply_lut_blocked(labels, lut: np.ndarray, block: int = DEFAULT_BLOCK) -> None:
    """就地套用 label → label 的查表（分塊，適用 memmap）。"""
    h, w = labels.shape
    for y0, x0, y1, x1 in iter_blocks(h, w, block):
        labels[y0:y1, x0:x1] = lut[np.asarray(labels[y0:y1, x0:x1])]


def clean_and_relabel_blocked(
    labels, min_size: int, max_size: int, block: int = DEFAULT_BLOCK
) -> int:
    """
    分塊版的「面積過濾 ＋ 連續重新編號」，回傳保留的細胞數。

    等價於 `clean_mask` ＋ `relabel_sequential`，但面積統計是**全域**的
    （見 `global_label_sizes`），因此跨塊細胞不會被誤判為過小。
    """
    sizes = global_label_sizes(labels, block)
    keep = (sizes >= min_size) & (sizes <= max_size)
    keep[0] = False                       # 背景
    n_keep = int(keep.sum())

    lut = np.zeros(len(sizes), dtype=np.int32)
    lut[keep] = np.arange(1, n_keep + 1, dtype=np.int32)
    apply_lut_blocked(labels, lut, block)
    return n_keep


def blocked_voronoi(
    labels,
    out,
    max_distance: int,
    tissue_fn=None,
    block: int = DEFAULT_BLOCK,
    margin: int = DEFAULT_MARGIN,
) -> None:
    """
    分塊 Voronoi 擴張：`labels` → `out`，避免一次建立全圖距離場。

    **為何是精確的**：`voronoi_expand` 只會把距離 ≤ `max_distance` 的背景像素
    指派出去，因此任一像素的結果只取決於半徑 `max_distance` 內的內容。只要
    `margin ≥ max_distance`，帶 margin 的塊內結果與全圖計算**完全一致**，
    跨塊細胞不會產生接縫。

    **必須寫到另一個陣列**：若原地覆寫，後續塊讀到的 margin 會是已擴張的值，
    等於把擴張結果再擴張一次。

    Parameters
    ----------
    tissue_fn : Callable[[int, int, int, int], np.ndarray] | None
        `(y0, x0, y1, x1) -> bool 陣列`，回傳該區域的組織遮罩。全片模式下
        組織遮罩不整份保存（21504×47104 bool ≈ 1 GB），改為按需重算。
    """
    if margin < max_distance:
        raise ValueError(f"margin ({margin}) 必須 ≥ voronoi 距離 ({max_distance})")

    h, w = labels.shape
    for y0, x0, y1, x1 in iter_blocks(h, w, block):
        ey0, ex0 = max(0, y0 - margin), max(0, x0 - margin)
        ey1, ex1 = min(h, y1 + margin), min(w, x1 + margin)

        sub = np.asarray(labels[ey0:ey1, ex0:ex1])
        tissue = tissue_fn(ey0, ex0, ey1, ex1) if tissue_fn is not None else None
        expanded = voronoi_expand(sub, max_distance=max_distance, tissue_mask=tissue)

        out[y0:y1, x0:x1] = expanded[y0 - ey0:y1 - ey0, x0 - ex0:x1 - ex0]


FULL_SEG_PROGRESS_FILENAME = "full_seg_progress.json"


def _save_seg_progress(path: Path, cfg_hash: str, done, current_max: int) -> None:
    """記錄已完成的 tile 與目前最大 label ID（供中斷續跑）。"""
    import json

    payload = {
        "config_hash": cfg_hash,
        "done_tiles": sorted([int(a), int(b)] for a, b in done),
        "current_max": int(current_max),
    }
    tmp = path.with_suffix(".tmp")
    tmp.write_text(json.dumps(payload), encoding="utf-8")
    tmp.replace(path)          # 原子替換：避免寫到一半被中斷而讀到半份進度


def _load_seg_progress(path: Path, cfg_hash: str) -> dict | None:
    """
    讀取續跑進度；**config hash 不符即視為無效**。

    參數改過就必須重跑 —— 否則會把兩組參數的 tile 拼在一起，產出一份看起來
    正常但實際錯誤的遮罩。
    """
    import json

    if not path.exists():
        return None
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (ValueError, OSError) as e:
        logger.warning(f"續跑進度檔損壞，改為重新開始：{e}")
        return None

    if data.get("config_hash") != cfg_hash:
        logger.info(
            f"分割參數已變更（{data.get('config_hash')} → {cfg_hash}），忽略舊進度重新開始"
        )
        return None
    return {
        "done_tiles": data.get("done_tiles", []),
        "current_max": int(data.get("current_max", 0)),
    }


def config_hash(cfg: dict, *extra) -> str:
    """
    分割設定的短雜湊，用於命名續跑用的暫存檔。

    參數不同的兩次執行必須拿到不同的 hash —— 否則中斷後重跑會把兩組參數的
    結果拼在一起，產出一份**看起來正常但實際錯誤**的遮罩。
    """
    import hashlib
    import json

    payload = json.dumps([cfg, *extra], sort_keys=True, default=str)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()[:12]


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
    from cellpose import core

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
        f"  主模型 dia: {dia_small}/{dia_mid}/{dia_large} | "
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

    logger.info(f"  [{time.time()-t0:.0f}s] 載入主模型...")
    primary = _load_primary_model(use_gpu)

    logger.info(f"  [{time.time()-t0:.0f}s] 主模型 dia={dia_mid} (RGB)...")
    m, _, _ = primary.eval(enhanced, diameter=dia_mid,
                         augment=True, resample=True, **eval_base)
    results["primary_mid"] = m
    logger.info(f"    → {m.max()} cells")

    logger.info(f"  [{time.time()-t0:.0f}s] 主模型 dia={dia_small} (small)...")
    m, _, _ = primary.eval(
        enhanced, diameter=dia_small, augment=False, resample=True,
        **{**eval_base, "cellprob_threshold": cellprob_thresh - 1.0},
    )
    results["primary_small"] = m
    logger.info(f"    → {m.max()} cells")

    logger.info(f"  [{time.time()-t0:.0f}s] 主模型 dia={dia_large} (large)...")
    m, _, _ = primary.eval(
        enhanced, diameter=dia_large, augment=False, resample=True,
        **{**eval_base, "cellprob_threshold": cellprob_thresh + 1.0},
    )
    results["primary_large"] = m
    logger.info(f"    → {m.max()} cells")

    if use_hematoxylin and hema is not None:
        hema_rgb = np.stack([hema, hema, hema], axis=-1)
        logger.info(f"  [{time.time()-t0:.0f}s] 主模型 dia={dia_mid} (hematoxylin)...")
        m, _, _ = primary.eval(hema_rgb, diameter=dia_mid,
                             augment=True, resample=True, **eval_base)
        results["primary_hema"] = m
        logger.info(f"    → {m.max()} cells")
        del hema_rgb

    del primary
    gc.collect()
    _clear_gpu_cache()

    if use_cpsam:
        logger.info(f"  [{time.time()-t0:.0f}s] 載入 cpsam...")
        try:
            cpsam = _load_primary_model(use_gpu)
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

    # ── 3. 集成合併（以 primary_mid 為基底）────────────────────
    base_mask = results["primary_mid"].copy().astype(np.int32)
    target_shape = base_mask.shape
    for key, mask in results.items():
        if key == "primary_mid":
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


def validate_roi_overrides(roi_overrides: dict) -> "tuple[dict, list[str]]":
    """檢查單一 ROI 的覆寫欄位。

    「哪些分割參數可以被 ROI 覆寫」是分割領域的知識，這裡是對外的唯一入口——
    呼叫端（API 驗證、參數合併）不需要也不應該直接碰 _ROI_OVERRIDE_FIELDS。

    Returns
    -------
    (clean, invalid)
        clean   : 可套用的覆寫（欄位名合法，值非 None）
        invalid : 未知欄位名清單，順序同輸入
    """
    clean = {
        k: v for k, v in roi_overrides.items()
        if k in _ROI_OVERRIDE_FIELDS and v is not None
    }
    invalid = [k for k in roi_overrides if k not in _ROI_OVERRIDE_FIELDS]
    return clean, invalid


def _merge_roi_params(global_seg_cfg: dict, roi_overrides: dict) -> dict:
    """將 ROI 個別覆寫合併進全域分割設定（深複製，不改原始 dict）。"""
    import copy
    cfg = copy.deepcopy(global_seg_cfg)
    mcseg = cfg.setdefault("mcseg_v2", {})
    clean, invalid = validate_roi_overrides(roi_overrides)
    if invalid:
        logger.warning(f"  忽略未知的 ROI 覆寫欄位：{invalid}")
    mcseg.update(clean)
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
    img: np.ndarray | None = None,
    cfg: dict | None = None,
    tile_size: int = 1024,
    overlap: int = 128,
    progress_callback=None,
    tile_reader=None,
    full_shape: tuple[int, int] | None = None,
    work_dir: str | Path | None = None,
    resume: bool = True,
) -> np.ndarray:
    """
    全圖 tiled MCseg v2 分割（支援串流讀取、memmap 落地、中斷續跑）。

    兩階段設計（避免 Voronoi 在 tile 邊界產生接縫）：
      Phase 1（per-tile）：讀圖 → CLAHE → Cellpose 多直徑推論 → merge_masks_fast
      Phase 2（分塊）：全域面積過濾 → 重新編號 → 分塊 Voronoi 擴張

    MPS 安全設定：
      - tile_size=1024（預設），比 2048 佔用更少 GPU 記憶體
      - augment=False 全程停用（augment=True 在 MPS 上易 OOM）
      - batch_size 限制為 cfg 設定值，建議 ≤ 2
      - 捕捉 MPS RuntimeError 後自動 fallback 到 CPU

    Args:
        img:               (H, W, 3) uint8 RGB 影像。與 `tile_reader` 二擇一。
        cfg:               mcseg_v2 設定 dict
        tile_size:         每塊大小（px），MPS 安全建議 1024
        overlap:           相鄰塊重疊寬度（px）
        progress_callback: fn(progress: float, message: str)
        tile_reader:       `fn(x, y, w, h) -> (h, w, 3) uint8`，串流讀取單塊影像。
                           給定時**整張影像不進 RAM**，全片（>6 GB）才跑得動。
        full_shape:        `(H, W)`；使用 `tile_reader` 時必填。
        work_dir:          給定時，標籤圖落在 `tmp_labels_{hash}.npy` memmap
                           （21504×47104 int32 ≈ 4 GB 走磁碟而非 RAM），
                           並啟用中斷續跑。
        resume:            `work_dir` 存在且 config hash 相符時跳過已完成的 tile。

    Returns:
        int32 細胞分割遮罩 (H, W)

    Notes:
        影像預處理（CLAHE、組織遮罩、Hematoxylin）改為 **per-tile（含 overlap）**，
        兩條路徑（ndarray / tile_reader）走同一份程式碼，結果逐位元一致。
        這同時讓全圖路徑與 ROI 路徑更可比 —— ROI 的 CLAHE 本來就是在一小塊
        影像上做的，而舊版全圖 CLAHE 的 8×8 網格落在整張片子上，等於近乎全域
        等化，與 ROI 結果不可直接比較。
    """
    from cellpose import core

    t0 = time.time()
    cfg = cfg or {}

    # ── 輸入來源：ndarray 或串流 reader ──────────────────────
    if (img is None) == (tile_reader is None):
        raise ValueError("`img` 與 `tile_reader` 必須且只能提供其中一個")

    if tile_reader is None:
        H, W = img.shape[:2]

        def tile_reader(x, y, w, h, _img=img):     # noqa: ARG001 - 統一介面
            return _img[y:y + h, x:x + w]
    else:
        if full_shape is None:
            raise ValueError("使用 `tile_reader` 時必須提供 `full_shape`")
        H, W = int(full_shape[0]), int(full_shape[1])

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

    y_starts = list(range(0, H, tile_size))
    x_starts = list(range(0, W, tile_size))
    total_tiles = len(y_starts) * len(x_starts)

    logger.info(
        f"[Tiled MCseg v2] 全圖 {W}×{H}px  "
        f"tile={tile_size}px overlap={overlap}px  "
        f"tiles={total_tiles}  gpu={use_gpu}  "
        f"{'串流' if img is None else 'in-memory'}"
    )

    # ── 標籤圖：memmap 落地（可續跑）或 RAM ──────────────────
    cfg_hash = config_hash(cfg, H, W, tile_size, overlap)
    labels_path: Path | None = None
    state = {"done_tiles": [], "current_max": 0}

    if work_dir is not None:
        work_dir = Path(work_dir)
        work_dir.mkdir(parents=True, exist_ok=True)
        labels_path = work_dir / f"tmp_labels_{cfg_hash}.npy"
        progress_path = work_dir / FULL_SEG_PROGRESS_FILENAME
        loaded = _load_seg_progress(progress_path, cfg_hash) if resume else None

        if loaded is not None and labels_path.exists():
            state = loaded
            stitched = np.lib.format.open_memmap(str(labels_path), mode="r+")
            logger.info(
                f"  續跑：已完成 {len(state['done_tiles'])}/{total_tiles} 個 tile"
                f"（config hash {cfg_hash}）"
            )
        else:
            stitched = np.lib.format.open_memmap(
                str(labels_path), mode="w+", dtype=np.int32, shape=(H, W)
            )
    else:
        progress_path = None
        stitched = np.zeros((H, W), dtype=np.int32)

    done = {tuple(t) for t in state["done_tiles"]}
    current_max = int(state["current_max"])

    logger.info(f"  載入主模型 (gpu={use_gpu})")
    primary = _load_primary_model(use_gpu)

    cpsam = None
    if use_cpsam:
        logger.info(f"  載入 cpsam 模型 (gpu={use_gpu})")
        try:
            cpsam = _load_primary_model(use_gpu)
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

    # ── Phase 1：per-tile 讀圖 → 前處理 → Cellpose ──────────
    for ti, y in enumerate(y_starts):
        y0e = y - overlap if y > 0 else 0
        y1e = min(y + tile_size + overlap, H)

        for tj, x in enumerate(x_starts):
            tile_idx = ti * len(x_starts) + tj + 1
            if (ti, tj) in done:
                logger.info(f"  [SKIP] Tile {tile_idx}/{total_tiles} ({x},{y}) 已完成")
                continue

            x0e = x - overlap if x > 0 else 0
            x1e = min(x + tile_size + overlap, W)

            msg = f"Tile {tile_idx}/{total_tiles} ({x},{y})"
            if progress_callback:
                progress_callback(tile_idx / total_tiles * 0.85, msg)
            logger.info(f"  [{time.time()-t0:.0f}s] {msg}")

            raw_tile = np.ascontiguousarray(
                tile_reader(x0e, y0e, x1e - x0e, y1e - y0e)
            )
            enh_tile = apply_clahe(raw_tile, clip_limit=clahe_clip, tile_size=8)
            tissue_tile = create_tissue_mask(raw_tile)
            hema_tile = color_deconvolution_he(raw_tile) if use_hematoxylin else None

            # per-tile 多直徑推論 + merge（不做 Voronoi）
            tile_results: dict[str, np.ndarray] = {}
            try:
                m, _, _ = primary.eval(enh_tile, diameter=dia_mid, **eval_base)
                tile_results["mid"] = m

                m, _, _ = primary.eval(
                    enh_tile, diameter=dia_small,
                    **{**eval_base, "cellprob_threshold": cellprob_thresh - 1.0},
                )
                tile_results["small"] = m

                m, _, _ = primary.eval(
                    enh_tile, diameter=dia_large,
                    **{**eval_base, "cellprob_threshold": cellprob_thresh + 1.0},
                )
                tile_results["large"] = m

                if hema_tile is not None:
                    hema_rgb = np.stack([hema_tile] * 3, axis=-1)
                    m, _, _ = primary.eval(hema_rgb, diameter=dia_mid, **eval_base)
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
                    cpu_model = _load_primary_model(False)
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

            # 組織遮罩就地套用：全片模式下不保留整份 tissue mask（bool 全圖 ≈ 1 GB）
            if tissue_tile.shape == base.shape:
                base[~tissue_tile] = 0

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
                prev_row = np.asarray(stitched[y - 1, x:x + valid.shape[1]])
                curr_row = valid[0, :len(prev_row)]
                mm = (prev_row > 0) & (curr_row > 0)
                for p, c in zip(prev_row[mm], curr_row[mm]):
                    mappings.setdefault(int(c), int(p))
            if x > 0:
                prev_col = np.asarray(stitched[y:y + valid.shape[0], x - 1])
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
            done.add((ti, tj))

        # 每完成一列就存檔：中斷後可自此列續跑
        if progress_path is not None:
            if hasattr(stitched, "flush"):
                stitched.flush()
            _save_seg_progress(progress_path, cfg_hash, done, current_max)

    del primary
    if cpsam is not None:
        del cpsam
    gc.collect()
    _clear_gpu_cache()

    # ── Phase 2：分塊清理 + Voronoi ──────────────────────────
    if progress_callback:
        progress_callback(0.90, "Voronoi 擴張（分塊）...")
    logger.info(f"  [{time.time()-t0:.0f}s] Phase 2：清理 + Voronoi 擴張（分塊）")

    n_before = clean_and_relabel_blocked(stitched, min_size, max_size)
    logger.info(f"  擴張前：{n_before} cells")

    def _tissue_fn(by0, bx0, by1, bx1):
        return create_tissue_mask(
            np.ascontiguousarray(tile_reader(bx0, by0, bx1 - bx0, by1 - by0))
        )

    if labels_path is not None:
        expanded_path = labels_path.with_name(f"tmp_expanded_{cfg_hash}.npy")
        final = np.lib.format.open_memmap(
            str(expanded_path), mode="w+", dtype=np.int32, shape=(H, W)
        )
    else:
        expanded_path = None
        final = np.zeros((H, W), dtype=np.int32)

    blocked_voronoi(stitched, final, max_distance=voronoi_dist, tissue_fn=_tissue_fn)
    n_final = clean_and_relabel_blocked(final, min_size, max_size)

    if progress_callback:
        progress_callback(1.0, f"完成：{n_final:,} 個細胞")
    logger.info(f"  [{time.time()-t0:.0f}s] 全圖分割完成：{n_final:,} cells")

    if expanded_path is not None:
        final.flush()
        result = np.array(final)     # 呼叫端需要一份可自由使用的陣列
        del final, stitched
        gc.collect()
        # 暫存檔只在最終遮罩確實產出後才刪 —— 提早刪會讓中斷的執行無從續跑
        for p in (labels_path, expanded_path):
            try:
                p.unlink()
            except OSError as e:
                logger.warning(f"  暫存檔清理失敗（不影響結果）：{p.name} — {e}")
        if progress_path is not None:
            try:
                progress_path.unlink()
            except OSError:
                pass
        return result.astype(np.int32)

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
