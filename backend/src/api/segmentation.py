"""Stage 1: 細胞分割 API"""
import asyncio
import base64
import logging
from pathlib import Path
from typing import Optional
from fastapi import APIRouter, BackgroundTasks
from pydantic import BaseModel

from backend.src.utils.config import load_config, resolve_path
from backend.src.utils.logging import set_current_stage

router = APIRouter()
logger = logging.getLogger("pipeline.api.segmentation")

_task_status = {"status": "idle", "progress": 0.0, "message": ""}
_task_lock = asyncio.Lock()


class SegmentationParams(BaseModel):
    """前端傳入的分割參數覆寫（所有欄位皆可選）"""
    mode: str = "roi"                         # "roi"（各 ROI he_crop.tif）或 "full"（完整 BTF 全圖）
    model_type: Optional[str] = None          # cyto2 / cyto3 / nuclei
    use_gpu: Optional[bool] = None
    batch_size: Optional[int] = None
    dia_small: Optional[float] = None
    dia_large: Optional[float] = None
    flow_threshold: Optional[float] = None
    cellprob_threshold: Optional[float] = None
    fragment_threshold: Optional[int] = None
    normalize_stains: Optional[bool] = None
    clahe_clip_limit: Optional[float] = None
    enable_eosin_watershed: Optional[bool] = None
    eosin_bg_threshold: Optional[int] = None
    block_size: Optional[int] = None
    overlap: Optional[int] = None


class PreviewRequest(BaseModel):
    """快速 Patch 預覽請求（對 he_crop.tif 取小塊跑 Cellpose）"""
    roi_name: Optional[str] = None   # 若為空則自動選第一個有 he_crop.tif 的 ROI
    x: int = 0
    y: int = 0
    patch_size: int = 512
    # param overrides（繼承自 SegmentationParams，不含 block/overlap/watershed）
    model_type: Optional[str] = None
    use_gpu: Optional[bool] = None
    dia_small: Optional[float] = None
    dia_large: Optional[float] = None
    flow_threshold: Optional[float] = None
    cellprob_threshold: Optional[float] = None
    fragment_threshold: Optional[int] = None
    normalize_stains: Optional[bool] = None
    clahe_clip_limit: Optional[float] = None  # 即時 CLAHE 預覽
    enable_eosin_watershed: Optional[bool] = None   # ← 新增
    eosin_bg_threshold: Optional[int] = None        # ← 新增


def _apply_overrides(config: dict, p: SegmentationParams) -> dict:
    """將前端傳入的參數覆寫到 config dict（只覆寫非 None 的欄位）"""
    seg = config.setdefault("segmentation", {})

    model = seg.setdefault("cellpose_model", {})
    if p.model_type is not None:
        model["model_type"] = p.model_type
    if p.use_gpu is not None:
        model["use_gpu"] = p.use_gpu
    if p.batch_size is not None:
        model["batch_size"] = p.batch_size

    strategy = seg.setdefault("strategy", {})
    if p.dia_small is not None:
        strategy["dia_small"] = p.dia_small
    if p.dia_large is not None:
        strategy["dia_large"] = p.dia_large
    if p.flow_threshold is not None:
        strategy["flow_threshold"] = p.flow_threshold
    if p.cellprob_threshold is not None:
        strategy["cellprob_threshold"] = p.cellprob_threshold
    if p.fragment_threshold is not None:
        strategy["fragment_threshold"] = p.fragment_threshold

    preproc = seg.setdefault("preprocessing", {})
    if p.normalize_stains is not None:
        preproc["normalize_stains"] = p.normalize_stains
    if p.clahe_clip_limit is not None:
        preproc["clahe_clip_limit"] = p.clahe_clip_limit

    postproc = seg.setdefault("postprocessing", {})
    if p.enable_eosin_watershed is not None:
        postproc["enable_eosin_watershed"] = p.enable_eosin_watershed
    if p.eosin_bg_threshold is not None:
        postproc["eosin_bg_threshold"] = p.eosin_bg_threshold

    tiling = seg.setdefault("tiling", {})
    if p.block_size is not None:
        tiling["block_size"] = p.block_size
    if p.overlap is not None:
        tiling["overlap"] = p.overlap

    return config


def _run_preview_sync(he_crop_path: Path, req: PreviewRequest, config: dict) -> dict:
    """同步執行 patch 預覽（在 executor 中呼叫）"""
    import io
    import numpy as np
    import cv2
    import tifffile
    from PIL import Image
    from skimage.segmentation import find_boundaries
    from cellpose import models, core

    seg_cfg  = config.get("segmentation", {})
    mc       = seg_cfg.get("cellpose_model", {})
    st       = seg_cfg.get("strategy", {})
    prep     = seg_cfg.get("preprocessing", {})

    model_type    = req.model_type     or mc.get("model_type",    "cyto2")
    use_gpu       = (req.use_gpu if req.use_gpu is not None else mc.get("use_gpu", True)) and core.use_gpu()
    dia_small     = req.dia_small      if req.dia_small      is not None else st.get("dia_small",      30.0)
    dia_large     = req.dia_large      if req.dia_large      is not None else st.get("dia_large",      60.0)
    flow_thresh   = req.flow_threshold if req.flow_threshold is not None else st.get("flow_threshold",   0.4)
    cellprob_thresh = req.cellprob_threshold if req.cellprob_threshold is not None else st.get("cellprob_threshold", -1.0)
    frag_thresh   = req.fragment_threshold   if req.fragment_threshold  is not None else st.get("fragment_threshold",  200)
    normalize     = req.normalize_stains     if req.normalize_stains     is not None else prep.get("normalize_stains", True)

    # ── 讀取 patch ────────────────────────────────────────────────────────────
    img_full = tifffile.imread(str(he_crop_path))
    if img_full.ndim == 3 and img_full.shape[-1] == 4:
        img_full = img_full[..., :3]
    h, w = img_full.shape[:2]
    x0 = max(0, req.x);  x1 = min(w, x0 + req.patch_size)
    y0 = max(0, req.y);  y1 = min(h, y0 + req.patch_size)
    patch = img_full[y0:y1, x0:x1]
    if patch.size == 0:
        return {"status": "error", "message": f"座標超出影像範圍（{w}×{h}）"}

    # ── Macenko → Hematoxylin（與 cellpose_runner 相同流程）────────────────────
    from backend.src.segmentation.macenko import MacenkoNormalizer, apply_clahe
    normalizer = MacenkoNormalizer()
    if normalize and patch.ndim == 3 and patch.shape[-1] == 3:
        if normalizer.fit(patch):
            gray = normalizer.extract_hematoxylin(patch)
        else:
            gray = cv2.cvtColor(patch[..., :3], cv2.COLOR_RGB2GRAY)
    else:
        gray = cv2.cvtColor(patch[..., :3], cv2.COLOR_RGB2GRAY) if patch.ndim == 3 else patch.squeeze()
    # 優先使用 req 傳入值，fallback 到 config
    clip_limit = float(req.clahe_clip_limit if req.clahe_clip_limit is not None else prep.get("clahe_clip_limit", 2.0))
    gray = apply_clahe(gray, clip_limit=clip_limit)
    input_img = np.stack([gray, gray, gray], axis=-1)

    # ── Cellpose 雙尺寸推論 + LOGIC_A 合併 ────────────────────────────────────
    model = models.CellposeModel(gpu=use_gpu, pretrained_model=model_type)
    masks_s, flows_s, _ = model.eval(input_img, diameter=dia_small,
                                     flow_threshold=flow_thresh, cellprob_threshold=cellprob_thresh)
    masks_l, _, _ = model.eval(input_img, diameter=dia_large,
                               flow_threshold=flow_thresh, cellprob_threshold=cellprob_thresh)
    from backend.src.segmentation.cellpose_runner import _merge_masks_logic_a
    merged = _merge_masks_logic_a(masks_s, masks_l, frag_thresh)

    # ── 疊圖（H&E + 綠色細胞邊界）─────────────────────────────────────────────
    overlay = patch[..., :3].copy() if patch.ndim == 3 else np.stack([patch] * 3, axis=-1)
    boundaries = find_boundaries(merged, mode='thick')
    overlay[boundaries] = [50, 255, 80]   # 亮綠邊界

    n_cells = int(len(np.unique(merged)) - 1)
    cv2.putText(overlay, f"n={n_cells}", (8, 22),
                cv2.FONT_HERSHEY_SIMPLEX, 0.6, (255, 255, 0), 2)

    buf = io.BytesIO()
    Image.fromarray(overlay.astype(np.uint8)).save(buf, "JPEG", quality=85)

    # ── Macenko 前處理圖（Cellpose 實際看到的輸入）────────────────────────────
    mac_buf = io.BytesIO()
    Image.fromarray(input_img.astype(np.uint8)).save(mac_buf, "JPEG", quality=85)

    # ── Flow Direction（HSV 光流方向圖）─────────────────────────────────────
    # Cellpose 4.x flows: [0]=XY viz (H,W,3), [1]=dP (2,H,W), [2]=cellprob (H,W)
    dP = flows_s[1]
    dy, dx = dP[0].astype(np.float32), dP[1].astype(np.float32)
    angle = np.arctan2(dy, dx)                                   # -π ~ π
    mag   = np.sqrt(dy ** 2 + dx ** 2)
    hue   = ((angle + np.pi) / (2 * np.pi) * 179).astype(np.uint8)
    sat   = np.clip(mag / (mag.max() + 1e-6) * 255, 0, 255).astype(np.uint8)
    val   = np.full_like(hue, 220)
    flow_rgb = cv2.cvtColor(np.stack([hue, sat, val], axis=-1), cv2.COLOR_HSV2RGB)
    # 把 merged 邊界疊在 flow 圖上（白色）方便對照細胞位置
    flow_rgb[boundaries] = [255, 255, 255]
    flow_buf = io.BytesIO()
    Image.fromarray(flow_rgb).save(flow_buf, "JPEG", quality=85)

    # ── Cyto Mask (Eosin Watershed) ───────────────────────────────────────────
    cyto_b64 = None
    enable_cyto = req.enable_eosin_watershed if hasattr(req, 'enable_eosin_watershed') and getattr(req, 'enable_eosin_watershed') is not None else prep.get("enable_eosin_watershed", False)
    if enable_cyto and patch.ndim == 3 and patch.shape[-1] == 3:
        bg_thresh = req.eosin_bg_threshold if hasattr(req, 'eosin_bg_threshold') and getattr(req, 'eosin_bg_threshold') is not None else prep.get("eosin_bg_threshold", 40)
        # 兒素主要標記組織在：用最大通道亮度判斷被刑物（純白空直 = 所有通道都露忧）
        # 背景（空直）：Brightness 高 → 廢止 Proseg 擴展
        # 組織（細胞質+核）：Brightness 中/低 → 允許 Proseg 在此擴展
        brightness = patch.astype(np.float32).max(axis=2)  # (H, W) max(R,G,B)
        is_background = (brightness > (255 - bg_thresh))   # 白色區域
        cyto_mask = ~is_background                          # 組織區域
        # 視覺化：背景變暗，組織保留
        cyto_overlay = patch.astype(np.float32).copy()
        cyto_overlay[is_background] = cyto_overlay[is_background] * 0.2 + np.array([20, 20, 20]) * 0.8
        from skimage.segmentation import find_boundaries
        border = find_boundaries(cyto_mask.astype(np.uint8), mode='outer')
        cyto_overlay[border] = [0, 255, 200]
        cyto_overlay = np.clip(cyto_overlay, 0, 255).astype(np.uint8)
        cyto_buf = io.BytesIO()
        Image.fromarray(cyto_overlay).save(cyto_buf, "JPEG", quality=85)
        cyto_b64 = base64.b64encode(cyto_buf.getvalue()).decode()

    return {
        "status": "ok",
        "data": {
            "image_b64":   base64.b64encode(buf.getvalue()).decode(),
            "macenko_b64": base64.b64encode(mac_buf.getvalue()).decode(),
            "flows_b64":   base64.b64encode(flow_buf.getvalue()).decode(),
            "cyto_b64":    cyto_b64,
            "n_cells":     n_cells,
            "roi_name":    he_crop_path.parent.name,
            "patch_info":  f"x={x0},y={y0} size={x1-x0}×{y1-y0}",
        },
    }


@router.get("/status")
async def get_status():
    return _task_status


async def _run_segmentation(config: dict):
    global _task_status
    set_current_stage("segmentation")
    _task_status = {"status": "running", "progress": 0.0, "message": "啟動 Cellpose..."}

    def _progress(progress: float, message: str):
        _task_status.update({"progress": progress, "message": message})

    mode = config.get("_mode", "roi")
    try:
        if mode == "full":
            from backend.src.segmentation.cellpose_runner import run_segmentation
            await asyncio.get_event_loop().run_in_executor(None, run_segmentation, config)
        else:
            from backend.src.segmentation.cellpose_runner import run_segmentation_rois
            await asyncio.get_event_loop().run_in_executor(None, run_segmentation_rois, config, _progress)
        _task_status = {"status": "done", "progress": 1.0, "message": "分割完成"}
    except Exception as e:
        logger.error(f"分割失敗：{e}")
        _task_status = {"status": "error", "progress": 0.0, "message": str(e)}


@router.post("/run")
async def run_segmentation(background_tasks: BackgroundTasks, params: SegmentationParams = SegmentationParams()):
    async with _task_lock:
        if _task_status["status"] == "running":
            return {"status": "error", "message": "任務執行中"}
        config = load_config()
        config = _apply_overrides(config, params)
        config["_mode"] = params.mode
        background_tasks.add_task(_run_segmentation, config)
    return {"status": "ok", "message": "細胞分割已啟動"}


@router.post("/preview_preproc")
async def preview_preproc(req: PreviewRequest):
    """
    ⚡ 極速前處理預覽（不跑 Cellpose）
    只執行 Macenko + CLAHE，回傳原始 H&E 與前處理後灰階圖並排比較。
    通常 < 1 秒完成，適合快速調整 clip_limit 參數。
    """
    import io
    import numpy as np
    import cv2
    import tifffile
    from PIL import Image

    config = load_config()
    output_dir = config.get("paths", {}).get("output_dir", "results/analysis")
    rois = config.get("rois", [])

    # 找 he_crop.tif（與 run_preview 相同邏輯）
    he_crop_path: Optional[Path] = None
    if req.roi_name:
        candidate = Path(output_dir) / "roi" / req.roi_name / "he_crop.tif"
        if candidate.exists():
            he_crop_path = candidate
    if he_crop_path is None:
        for roi in rois:
            candidate = Path(output_dir) / "roi" / roi["name"] / "he_crop.tif"
            if candidate.exists():
                he_crop_path = candidate
                break
    if he_crop_path is None:
        roi_base = Path(output_dir) / "roi"
        if roi_base.exists():
            for d in sorted(roi_base.iterdir()):
                candidate = d / "he_crop.tif"
                if candidate.exists():
                    he_crop_path = candidate
                    break
    if he_crop_path is None:
        return {"status": "error", "message": "找不到 he_crop.tif"}

    try:
        from backend.src.segmentation.macenko import MacenkoNormalizer, apply_clahe

        seg_cfg = config.get("segmentation", {})
        prep = seg_cfg.get("preprocessing", {})

        img_full = tifffile.imread(str(he_crop_path))
        if img_full.ndim == 3 and img_full.shape[-1] == 4:
            img_full = img_full[..., :3]
        h, w = img_full.shape[:2]

        # 取 patch
        x0 = max(0, req.x);  x1 = min(w, x0 + req.patch_size)
        y0 = max(0, req.y);  y1 = min(h, y0 + req.patch_size)
        patch = img_full[y0:y1, x0:x1]
        if patch.size == 0:
            return {"status": "error", "message": f"座標超出影像範圍（{w}×{h}）"}

        # Macenko + CLAHE
        normalize = req.normalize_stains if req.normalize_stains is not None else prep.get("normalize_stains", True)
        clip_limit = float(req.clahe_clip_limit if req.clahe_clip_limit is not None else prep.get("clahe_clip_limit", 2.0))

        normalizer = MacenkoNormalizer()
        if normalize and patch.ndim == 3 and patch.shape[-1] == 3:
            success = normalizer.fit(patch)
            if success:
                gray = normalizer.extract_hematoxylin(patch)
                method = f"Macenko + CLAHE({clip_limit})"
            else:
                gray = cv2.cvtColor(patch[..., :3], cv2.COLOR_RGB2GRAY)
                method = f"Grayscale fallback + CLAHE({clip_limit})"
        else:
            gray = cv2.cvtColor(patch[..., :3], cv2.COLOR_RGB2GRAY) if patch.ndim == 3 else patch.squeeze()
            method = f"Grayscale + CLAHE({clip_limit})"
            method_combined_text = f"Grayscale + CLAHE({clip_limit})"

        gray_clahe = apply_clahe(gray, clip_limit=clip_limit)

        # 建立並排比較圖（原始 H&E | 前處理灰階）
        ph, pw = patch.shape[:2]
        
        # ── 計算 Cyto Mask (如果開啟 Eosin Watershed) ─────────────────────────────────
        cyto_b64 = None
        enable_cyto = req.enable_eosin_watershed if hasattr(req, 'enable_eosin_watershed') and getattr(req, 'enable_eosin_watershed') is not None else prep.get("enable_eosin_watershed", False)
        if enable_cyto and patch.ndim == 3 and patch.shape[-1] == 3:
            bg_thresh = req.eosin_bg_threshold if hasattr(req, 'eosin_bg_threshold') and getattr(req, 'eosin_bg_threshold') is not None else prep.get("eosin_bg_threshold", 40)
            brightness = patch.astype(np.float32).max(axis=2)
            is_background = (brightness > (255 - bg_thresh))
            cyto_mask = ~is_background
            cyto_overlay = patch.astype(np.float32).copy()
            cyto_overlay[is_background] = cyto_overlay[is_background] * 0.2 + np.array([20, 20, 20]) * 0.8
            from skimage.segmentation import find_boundaries
            border = find_boundaries(cyto_mask.astype(np.uint8), mode='outer')
            cyto_overlay[border] = [0, 255, 200]
            cyto_overlay = np.clip(cyto_overlay, 0, 255).astype(np.uint8)
            cyto_buf = io.BytesIO()
            Image.fromarray(cyto_overlay).save(cyto_buf, "JPEG", quality=85)
            cyto_b64 = base64.b64encode(cyto_buf.getvalue()).decode()

        # 原圖 (Patch) 為對照
        raw_buf = io.BytesIO()
        Image.fromarray(patch.astype(np.uint8)).save(raw_buf, "JPEG", quality=85)

        # Macenko 前處理圖（Cellpose 實際看到的輸入）
        input_img_for_mac_buf = np.stack([gray_clahe, gray_clahe, gray_clahe], axis=-1)
        mac_buf = io.BytesIO()
        Image.fromarray(input_img_for_mac_buf.astype(np.uint8)).save(mac_buf, "JPEG", quality=85)

        # Determine method string for the return dictionary
        method_return_dict = "Macenko Hematoxylin + CLAHE" if (normalize and patch.ndim==3 and patch.shape[-1]==3 and normalizer.fit(patch)) else "RGB2GRAY + CLAHE"

        # The original combined image logic for display
        he_rgb = patch[..., :3] if patch.ndim == 3 else np.stack([patch] * 3, axis=-1)
        gray_rgb = np.stack([gray_clahe, gray_clahe, gray_clahe], axis=-1)

        # 分隔線
        sep = np.full((ph, 4, 3), 80, dtype=np.uint8)
        combined = np.concatenate([he_rgb, sep, gray_rgb], axis=1)

        # 標注文字
        font = cv2.FONT_HERSHEY_SIMPLEX
        cv2.putText(combined, "H&E 原始", (8, 22), font, 0.55, (255,255,80), 1)
        cv2.putText(combined, method, (pw + 12, 22), font, 0.45, (80,255,160), 1)

        buf = io.BytesIO()
        Image.fromarray(combined.astype(np.uint8)).save(buf, "JPEG", quality=88)

        return {
            "status": "ok",
            "data": {
                "image_b64": base64.b64encode(raw_buf.getvalue()).decode(),  # 原圖
                "macenko_b64": base64.b64encode(mac_buf.getvalue()).decode(), # 灰階圖
                "cyto_b64": cyto_b64, # 新增 Cyto Mask
                "method": method_return_dict,
                "patch_info": f"x={x0},y={y0} size={x1-x0}×{y1-y0}",
                "clip_limit": clip_limit,
            }
        }
    except Exception as e:
        logger.error(f"前處理預覽失敗：{e}", exc_info=True)
        return {"status": "error", "message": str(e)}


@router.post("/run_preview")
async def run_preview(req: PreviewRequest):
    """對指定 ROI 的 he_crop.tif 取 patch 執行快速 Cellpose，回傳疊圖 base64。"""
    config = load_config()
    output_dir = config.get("paths", {}).get("output_dir", "results/analysis")
    rois       = config.get("rois", [])

    # 找 he_crop.tif
    he_crop_path: Optional[Path] = None
    if req.roi_name:
        candidate = Path(output_dir) / "roi" / req.roi_name / "he_crop.tif"
        if candidate.exists():
            he_crop_path = candidate

    if he_crop_path is None:
        for roi in rois:
            candidate = Path(output_dir) / "roi" / roi["name"] / "he_crop.tif"
            if candidate.exists():
                he_crop_path = candidate
                break

    if he_crop_path is None:
        # 掃描所有 roi 子目錄
        roi_base = Path(output_dir) / "roi"
        if roi_base.exists():
            for d in sorted(roi_base.iterdir()):
                candidate = d / "he_crop.tif"
                if candidate.exists():
                    he_crop_path = candidate
                    break

    if he_crop_path is None:
        return {"status": "error", "message": "找不到 he_crop.tif，請先在 Stage 0 執行 ROI 裁切"}

    logger.info(f"Preview patch: {he_crop_path} x={req.x} y={req.y} size={req.patch_size}")
    try:
        result = await asyncio.get_event_loop().run_in_executor(
            None, _run_preview_sync, he_crop_path, req, config
        )
        return result
    except Exception as e:
        logger.error(f"Preview 失敗：{e}", exc_info=True)
        return {"status": "error", "message": str(e)}


@router.get("/preview")
async def get_preview(roi_name: Optional[str] = None):
    """
    回傳指定 ROI 的分割遮罩疊加在 H&E 上的預覽圖（base64 JPEG）。
    roi_name 省略時取第一個有遮罩的 ROI。
    同時回傳所有可用 ROI 列表（供前端下拉選擇）。
    """
    import numpy as np
    import tifffile
    import cv2
    from skimage.segmentation import find_boundaries

    config = load_config()
    paths  = config.get("paths", {})
    output_dir        = resolve_path(paths.get("output_dir", "results/analysis"))
    mask_tif_filename = config.get("segmentation", {}).get("output", {}).get("mask_tif_filename", "segmentation_masks.tif")

    roi_base = output_dir / "roi"
    # 掃描所有有遮罩的 ROI
    available_rois = []
    if roi_base.exists():
        for d in sorted(roi_base.iterdir()):
            if d.is_dir() and not d.name.startswith(".") and (d / mask_tif_filename).exists():
                available_rois.append(d.name)

    if not available_rois:
        return {"status": "error", "message": "尚未執行分割，找不到遮罩檔案"}

    # 選定 ROI
    target = roi_name if roi_name in available_rois else available_rois[0]
    roi_dir   = roi_base / target
    mask_path = roi_dir / mask_tif_filename
    he_path   = roi_dir / "he_crop.tif"

    try:
        mask = tifffile.imread(str(mask_path)).astype(np.int32)

        # 嘗試讀取 H&E 背景
        if he_path.exists():
            he = tifffile.imread(str(he_path))
            if he.ndim == 2:
                he = np.stack([he, he, he], axis=-1)
            elif he.shape[-1] == 4:
                he = he[..., :3]
            # 若尺寸不符（不應發生），fallback 到純色背景
            if he.shape[:2] != mask.shape[:2]:
                he = np.full((*mask.shape, 3), 240, dtype=np.uint8)
        else:
            he = np.full((*mask.shape, 3), 240, dtype=np.uint8)

        # 疊加：繪製細胞邊界（亮綠色）+ 細胞核填色（半透明）
        overlay = he.copy().astype(np.float32)
        boundaries = find_boundaries(mask, mode="thick")
        overlay[boundaries] = [50, 255, 80]   # 亮綠邊界

        overlay = overlay.astype(np.uint8)
        n_cells = int(len(np.unique(mask)) - 1)  # 扣除背景 0

        # 縮放至長邊 1200px 輸出
        h, w = overlay.shape[:2]
        max_dim = 1200
        scale = min(max_dim / h, max_dim / w, 1.0)
        if scale < 1.0:
            overlay = cv2.resize(overlay, (int(w * scale), int(h * scale)))

        cv2.putText(overlay, f"n={n_cells} cells", (8, 28),
                    cv2.FONT_HERSHEY_SIMPLEX, 0.7, (255, 255, 0), 2)

        _, buf = cv2.imencode(".jpg", cv2.cvtColor(overlay, cv2.COLOR_RGB2BGR),
                              [cv2.IMWRITE_JPEG_QUALITY, 82])
        img_b64 = base64.b64encode(buf.tobytes()).decode()

        # ── Cyto Mask（如果存在 cyto_mask.npy）────────────────────────────────
        cyto_b64 = None
        cyto_npy_path = roi_dir / "cyto_mask.npy"
        if cyto_npy_path.exists() and he_path.exists():
            try:
                import io as _io
                from PIL import Image as _Image
                cyto_mask_arr = np.load(str(cyto_npy_path))  # 0=背景, 1=組織
                he_cyto = tifffile.imread(str(he_path))
                if he_cyto.ndim == 2:
                    he_cyto = np.stack([he_cyto]*3, axis=-1)
                elif he_cyto.shape[-1] == 4:
                    he_cyto = he_cyto[..., :3]
                # 確保 cyto_mask 與 he_cyto 尺寸一致
                if he_cyto.shape[:2] == cyto_mask_arr.shape[:2]:
                    is_bg = cyto_mask_arr == 0
                    cyto_vis = he_cyto.astype(np.float32)
                    cyto_vis[is_bg] = cyto_vis[is_bg] * 0.2 + np.array([20, 20, 20]) * 0.8
                    from skimage.segmentation import find_boundaries as _fb
                    border = _fb((~is_bg).astype(np.uint8), mode='outer')
                    cyto_vis[border] = [0, 255, 200]
                    cyto_vis = np.clip(cyto_vis, 0, 255).astype(np.uint8)
                    # 縮放
                    if scale < 1.0:
                        cyto_vis = cv2.resize(cyto_vis, (int(he_cyto.shape[1]*scale), int(he_cyto.shape[0]*scale)))
                    cb = _io.BytesIO()
                    _Image.fromarray(cyto_vis).save(cb, "JPEG", quality=82)
                    cyto_b64 = base64.b64encode(cb.getvalue()).decode()
            except Exception:
                pass

        # ── Flows Preview（如果存在 flows_preview.jpg）────────────────────────────────
        flows_b64 = None
        flows_jpg_path = roi_dir / "flows_preview.jpg"
        if flows_jpg_path.exists():
            try:
                with open(flows_jpg_path, "rb") as f:
                    flows_b64 = base64.b64encode(f.read()).decode()
            except Exception:
                pass

        return {
            "status": "ok",
            "data": {
                "image_b64":      img_b64,
                "cyto_b64":       cyto_b64,
                "flows_b64":      flows_b64,
                "roi":            target,
                "n_cells":        n_cells,
                "available_rois": available_rois,
            },
        }
    except Exception as e:
        return {"status": "error", "message": str(e)}
