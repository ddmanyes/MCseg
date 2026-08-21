"""Stage 1: MCseg v2 細胞分割 API"""
from __future__ import annotations

import asyncio
import base64
import logging
from pathlib import Path
from typing import Optional

import cv2
import numpy as np
from fastapi import APIRouter, BackgroundTasks
from pydantic import BaseModel, Field

from backend.src.fullslide.pipeline import (
    _SENTINEL_FULL,
    _upper_is_open,
    run_full_slide_segmentation,
)
from backend.src.utils.config import load_config, load_state, resolve_path, save_state_key
from backend.src.utils.logging import set_current_stage

router = APIRouter()
logger = logging.getLogger("pipeline.api.segmentation")

_task_status: dict = {"status": "idle", "progress": 0.0, "message": ""}
_task_lock = asyncio.Lock()

_full_status: dict = {"status": "idle", "progress": 0.0, "message": ""}
_full_lock   = asyncio.Lock()

_PREVIEW_JPEG_QUALITY = 85


class SegmentationParams(BaseModel):
    """前端傳入的 MCseg v2 參數覆寫（所有欄位皆可選）。"""

    mode: str = "roi"

    use_gpu: Optional[bool] = None
    batch_size: Optional[int] = None

    dia_small: Optional[float] = None
    dia_mid: Optional[float] = None
    dia_large: Optional[float] = None

    use_hematoxylin: Optional[bool] = None
    use_cpsam: Optional[bool] = None

    dia_cpsam_auto: Optional[float] = None
    dia_cpsam_small: Optional[float] = None
    cellprob_cpsam_auto: Optional[float] = None
    cellprob_cpsam_small: Optional[float] = None
    cellprob_cpsam_hema: Optional[float] = None

    voronoi_distance: Optional[int] = None
    min_size: Optional[int] = None
    max_size: Optional[int] = None

    flow_threshold: Optional[float] = None
    cellprob_threshold: Optional[float] = None
    clahe_clip_limit: Optional[float] = None

    use_transcript_rescue: Optional[bool] = None

    roi_overrides: Optional[dict] = None
    target_roi: Optional[str] = None


class PreviewRequest(BaseModel):
    """快速 Patch 預覽請求"""

    roi_name: Optional[str] = None
    x: int = Field(0, ge=0)
    y: int = Field(0, ge=0)
    patch_size: int = Field(512, ge=64, le=2048)

    use_gpu: Optional[bool] = None
    batch_size: Optional[int] = None
    dia_small: Optional[float] = None
    dia_mid: Optional[float] = None
    dia_large: Optional[float] = None
    use_hematoxylin: Optional[bool] = None
    use_cpsam: Optional[bool] = None
    voronoi_distance: Optional[int] = None
    flow_threshold: Optional[float] = None
    cellprob_threshold: Optional[float] = None
    clahe_clip_limit: Optional[float] = None


class FullSegParams(BaseModel):
    """全圖分割請求：裁切窗格與 cpsam 開關（全部可選）。

    裁切座標為**原始影像 fullres px**。`None` 或 `-1` 代表該邊取影像邊界，
    因此不帶任何參數即為「整張影像」（維持既有行為）。
    """

    crop_x0: Optional[int] = None
    crop_x1: Optional[int] = None
    crop_y0: Optional[int] = None
    crop_y1: Optional[int] = None
    use_cpsam: Optional[bool] = None


def validate_crop(p: FullSegParams) -> Optional[str]:
    """驗證裁切座標；合法回傳 None，否則回傳錯誤訊息。

    只做「與影像尺寸無關」的檢查（負值、上下界順序）；實際邊界夾制在
    `resolve_crop_window` 完成，因為那時才知道影像尺寸。
    """
    for name, v0 in (("crop_x0", p.crop_x0), ("crop_y0", p.crop_y0)):
        if v0 is not None and v0 < 0:
            return f"{name} 不可為負（收到 {v0}）"
    for name, v1 in (("crop_x1", p.crop_x1), ("crop_y1", p.crop_y1)):
        if v1 is not None and v1 < 0 and v1 != _SENTINEL_FULL:
            return f"{name} 不可為負（收到 {v1}）"

    if not _upper_is_open(p.crop_x1):
        x0 = p.crop_x0 or 0
        if p.crop_x1 <= x0:
            return f"crop_x1 必須大於 crop_x0（收到 {x0} → {p.crop_x1}）"
    if not _upper_is_open(p.crop_y1):
        y0 = p.crop_y0 or 0
        if p.crop_y1 <= y0:
            return f"crop_y1 必須大於 crop_y0（收到 {y0} → {p.crop_y1}）"
    return None


def _maybe_set(d: dict, key: str, val) -> None:
    if val is not None:
        d[key] = val


def _apply_overrides(config: dict, p: SegmentationParams) -> dict:
    seg   = config.setdefault("segmentation", {})
    mcseg = seg.setdefault("mcseg_v2", {})
    _maybe_set(mcseg, "use_gpu",               p.use_gpu)
    _maybe_set(mcseg, "batch_size",            p.batch_size)
    _maybe_set(mcseg, "dia_small",             p.dia_small)
    _maybe_set(mcseg, "dia_mid",               p.dia_mid)
    _maybe_set(mcseg, "dia_large",             p.dia_large)
    _maybe_set(mcseg, "use_hematoxylin",       p.use_hematoxylin)
    _maybe_set(mcseg, "use_cpsam",             p.use_cpsam)
    _maybe_set(mcseg, "dia_cpsam_auto",        p.dia_cpsam_auto)
    _maybe_set(mcseg, "dia_cpsam_small",       p.dia_cpsam_small)
    _maybe_set(mcseg, "cellprob_cpsam_auto",   p.cellprob_cpsam_auto)
    _maybe_set(mcseg, "cellprob_cpsam_small",  p.cellprob_cpsam_small)
    _maybe_set(mcseg, "cellprob_cpsam_hema",   p.cellprob_cpsam_hema)
    _maybe_set(mcseg, "voronoi_distance",      p.voronoi_distance)
    _maybe_set(mcseg, "min_size",              p.min_size)
    _maybe_set(mcseg, "max_size",              p.max_size)
    _maybe_set(mcseg, "flow_threshold",        p.flow_threshold)
    _maybe_set(mcseg, "cellprob_threshold",    p.cellprob_threshold)
    _maybe_set(mcseg, "clahe_clip_limit",      p.clahe_clip_limit)
    _maybe_set(mcseg, "use_transcript_rescue", p.use_transcript_rescue)
    return config


def _build_preview_mcseg_cfg(config: dict, req: PreviewRequest) -> dict:
    base = dict(config.get("segmentation", {}).get("mcseg_v2", {}))
    _maybe_set(base, "use_gpu",            req.use_gpu)
    _maybe_set(base, "batch_size",         req.batch_size)
    _maybe_set(base, "dia_small",          req.dia_small)
    _maybe_set(base, "dia_mid",            req.dia_mid)
    _maybe_set(base, "dia_large",          req.dia_large)
    _maybe_set(base, "use_hematoxylin",    req.use_hematoxylin)
    _maybe_set(base, "use_cpsam",          req.use_cpsam)
    _maybe_set(base, "voronoi_distance",   req.voronoi_distance)
    _maybe_set(base, "flow_threshold",     req.flow_threshold)
    _maybe_set(base, "cellprob_threshold", req.cellprob_threshold)
    _maybe_set(base, "clahe_clip_limit",   req.clahe_clip_limit)
    return base


def _find_he_crop(config: dict, roi_name: str | None) -> Path | None:
    output_dir = config.get("paths", {}).get("output_dir", "results/analysis")
    roi_base   = resolve_path(output_dir) / "roi"

    if roi_name:
        c = roi_base / roi_name / "he_crop.tif"
        if c.exists():
            return c

    for roi in config.get("rois", []):
        c = roi_base / roi.get("name", "") / "he_crop.tif"
        if c.exists():
            return c

    if roi_base.exists():
        for d in sorted(roi_base.iterdir()):
            c = d / "he_crop.tif"
            if c.exists():
                return c
    return None


async def _run_segmentation(config: dict) -> None:
    global _task_status
    set_current_stage("segmentation")
    _task_status = {"status": "running", "progress": 0.0, "message": "啟動 MCseg v2..."}

    def _progress(progress: float, message: str) -> None:
        _task_status.update({"progress": progress, "message": message})

    roi_overrides: dict = config.pop("_roi_overrides", None) or load_state().get(
        "roi_seg_overrides", {}
    )
    target_roi: str | None = config.pop("_target_roi", None)

    try:
        from backend.src.segmentation.cellpose_runner import run_segmentation_rois
        import functools

        fn = functools.partial(
            run_segmentation_rois, config, _progress, roi_overrides, target_roi
        )
        await asyncio.get_running_loop().run_in_executor(None, fn)
        _task_status = {"status": "done", "progress": 1.0, "message": "分割完成"}
    except Exception as e:
        logger.error(f"分割失敗：{e}", exc_info=True)
        _task_status = {"status": "error", "progress": 0.0, "message": "分割失敗，請查閱 log"}


@router.get("/status")
async def get_status():
    return _task_status


@router.post("/run")
async def run_segmentation(
    background_tasks: BackgroundTasks,
    params: SegmentationParams = SegmentationParams(),
):
    async with _task_lock:
        if _task_status["status"] == "running":
            return {"status": "error", "message": "任務執行中"}
        config = load_config()
        config = _apply_overrides(config, params)
        overrides = params.roi_overrides or {}
        save_state_key("roi_seg_overrides", overrides)
        config["_roi_overrides"] = overrides
        config["_target_roi"] = params.target_roi
        _task_status["status"] = "running"
        _task_status["progress"] = 0.0
        _task_status["message"] = "初始化..."
        background_tasks.add_task(_run_segmentation, config)
    return {"status": "ok", "message": "MCseg v2 分割已啟動"}


@router.get("/full_seg_status")
async def get_full_seg_status():
    return _full_status


@router.post("/run_full")
async def run_full_segmentation(
    background_tasks: BackgroundTasks,
    params: FullSegParams = FullSegParams(),
):
    """對 BTF 全圖（或指定裁切窗格）執行 tiled MCseg v2 分割（MPS 安全模式）。"""
    err = validate_crop(params)
    if err:
        return {"status": "error", "message": err}

    async with _full_lock:
        if _full_status["status"] == "running":
            return {"status": "error", "message": "全圖分割任務執行中"}
        config = load_config()
        _full_status["status"]   = "running"
        _full_status["progress"] = 0.0
        _full_status["message"]  = "初始化..."
        background_tasks.add_task(_run_full_segmentation, config, params)
    return {"status": "ok", "message": "全圖分割已啟動"}


async def _run_full_segmentation(
    config: dict, params: FullSegParams | None = None
) -> None:
    """薄層：executor 呼叫 `run_full_slide_segmentation` + 狀態轉譯。

    編排邏輯（tile_reader 組裝、OOM 防護、MPS 鉗制、metadata 寫入）已下沉到
    `fullslide.pipeline.run_full_slide_segmentation`（架構深化 P8），CLI 全片
    路徑共用同一份實作與安全防護。
    """
    global _full_status
    set_current_stage("segmentation")
    _full_status = {"status": "running", "progress": 0.0, "message": "讀取全圖影像..."}

    def _progress(p: float, msg: str) -> None:
        _full_status.update({"progress": p, "message": msg})

    params = params or FullSegParams()

    try:
        import functools

        loop = asyncio.get_running_loop()
        result = await loop.run_in_executor(
            None,
            functools.partial(
                run_full_slide_segmentation,
                config,
                crop_x0=params.crop_x0, crop_y0=params.crop_y0,
                crop_x1=params.crop_x1, crop_y1=params.crop_y1,
                use_cpsam=params.use_cpsam,
                progress=_progress,
            ),
        )

        scope = "全圖" if result.is_full_image else f"窗格 ({result.crop_x0},{result.crop_y0})"
        _full_status = {
            "status": "done", "progress": 1.0,
            "message": f"{scope}分割完成：{result.n_cells:,} 個細胞  →  {result.mask_path.name}",
            "n_cells": result.n_cells,
            "output": str(result.mask_path),
        }
    except Exception as e:
        logger.error(f"全圖分割失敗：{e}", exc_info=True)
        if isinstance(e, MemoryError):
            # MemoryError 訊息含有尺寸資訊但不含路徑，可安全回傳
            safe_msg = str(e).split("，請縮小")[0] + "，請縮小裁切範圍或改用 ROI 模式。"
        else:
            safe_msg = "全圖分割失敗，請查閱 log"
        _full_status = {"status": "error", "progress": 0.0, "message": safe_msg}


@router.get("/roi_overrides")
async def get_roi_overrides():
    return {"status": "ok", "data": load_state().get("roi_seg_overrides", {})}


@router.put("/roi_overrides")
async def put_roi_overrides(body: dict):
    from backend.src.segmentation.cellpose_runner import validate_roi_overrides

    # 驗證 ROI 名稱：只允許已存在於 config 的 ROI，防止路徑穿越攻擊
    config = load_config()
    valid_names = {r.get("name", "") for r in config.get("rois", [])}
    invalid_names = [k for k in body if k not in valid_names]
    if invalid_names:
        return {"status": "error", "message": f"未知的 ROI 名稱：{invalid_names}"}

    # 驗證每個 ROI 的覆寫欄位名稱（哪些欄位可覆寫由分割領域決定）
    for roi_name, overrides in body.items():
        if not isinstance(overrides, dict):
            return {"status": "error", "message": f"ROI '{roi_name}' 的覆寫值必須是 dict"}
        _, invalid_fields = validate_roi_overrides(overrides)
        if invalid_fields:
            return {"status": "error", "message": f"ROI '{roi_name}' 包含未知欄位：{invalid_fields}"}

    save_state_key("roi_seg_overrides", body)
    return {"status": "ok"}


def _run_preview_sync(he_crop_path: Path, req: PreviewRequest, config: dict) -> dict:
    import io
    import tifffile
    from PIL import Image
    from skimage.segmentation import find_boundaries
    from backend.src.segmentation.cellpose_runner import run_preview_patch, apply_clahe

    mcseg_cfg = _build_preview_mcseg_cfg(config, req)

    img_full = tifffile.imread(str(he_crop_path))
    if img_full.ndim == 3 and img_full.shape[-1] == 4:
        img_full = img_full[..., :3]
    if img_full.ndim == 2:
        img_full = np.stack([img_full, img_full, img_full], axis=-1)

    h, w = img_full.shape[:2]
    x0 = max(0, req.x);  x1 = min(w, x0 + req.patch_size)
    y0 = max(0, req.y);  y1 = min(h, y0 + req.patch_size)
    patch = img_full[y0:y1, x0:x1]
    if patch.size == 0:
        return {"status": "error", "message": f"座標超出影像範圍（{w}×{h}）"}

    merged = run_preview_patch(patch, mcseg_cfg)

    clip = float(mcseg_cfg.get("clahe_clip_limit", 3.0))
    enhanced = apply_clahe(patch, clip_limit=clip)
    mac_buf = io.BytesIO()
    Image.fromarray(enhanced.astype(np.uint8)).save(mac_buf, "JPEG",
                                                     quality=_PREVIEW_JPEG_QUALITY)

    overlay = patch.copy()
    boundaries = find_boundaries(merged, mode="thick")
    overlay[boundaries] = [50, 255, 80]
    n_cells = int(len(np.unique(merged)) - 1)
    cv2.putText(overlay, f"n={n_cells}", (8, 22),
                cv2.FONT_HERSHEY_SIMPLEX, 0.6, (255, 255, 0), 2)

    buf = io.BytesIO()
    Image.fromarray(overlay.astype(np.uint8)).save(buf, "JPEG",
                                                    quality=_PREVIEW_JPEG_QUALITY)

    return {
        "status": "ok",
        "data": {
            "image_b64":   base64.b64encode(buf.getvalue()).decode(),
            "clahe_b64": base64.b64encode(mac_buf.getvalue()).decode(),
            "flows_b64":   None,
            "n_cells":     n_cells,
            "roi_name":    he_crop_path.parent.name,
            "patch_info":  f"x={x0},y={y0} size={x1-x0}×{y1-y0}",
        },
    }


@router.post("/run_preview")
async def run_preview(req: PreviewRequest):
    config = load_config()
    he_crop_path = _find_he_crop(config, req.roi_name)
    if he_crop_path is None:
        return {"status": "error",
                "message": "找不到 he_crop.tif，請先在 Stage 0 執行 ROI 裁切"}
    logger.info(f"Preview: {he_crop_path} x={req.x} y={req.y} size={req.patch_size}")
    try:
        result = await asyncio.get_running_loop().run_in_executor(
            None, _run_preview_sync, he_crop_path, req, config
        )
        return result
    except Exception as e:
        logger.error(f"Preview 失敗：{e}", exc_info=True)
        return {"status": "error", "message": "預覽失敗，請查閱 log"}


@router.post("/preview_preproc")
async def preview_preproc(req: PreviewRequest):
    """極速前處理預覽（不跑分割）：原圖 + CLAHE 增強 + Hematoxylin 通道。"""
    import io
    import tifffile
    from PIL import Image
    from backend.src.segmentation.cellpose_runner import (
        apply_clahe,
        color_deconvolution_he,
    )

    config = load_config()
    he_crop_path = _find_he_crop(config, req.roi_name)
    if he_crop_path is None:
        return {"status": "error", "message": "找不到 he_crop.tif"}

    try:
        mcseg_cfg = _build_preview_mcseg_cfg(config, req)
        clip = float(mcseg_cfg.get("clahe_clip_limit", 3.0))

        img_full = tifffile.imread(str(he_crop_path))
        if img_full.ndim == 3 and img_full.shape[-1] == 4:
            img_full = img_full[..., :3]
        h, w = img_full.shape[:2]
        x0 = max(0, req.x);  x1 = min(w, x0 + req.patch_size)
        y0 = max(0, req.y);  y1 = min(h, y0 + req.patch_size)
        patch = img_full[y0:y1, x0:x1]
        if patch.size == 0:
            return {"status": "error",
                    "message": f"座標超出影像範圍（{w}×{h}）"}
        if patch.ndim == 2:
            patch = np.stack([patch, patch, patch], axis=-1)

        raw_buf = io.BytesIO()
        Image.fromarray(patch.astype(np.uint8)).save(raw_buf, "JPEG",
                                                      quality=_PREVIEW_JPEG_QUALITY)

        enhanced = apply_clahe(patch, clip_limit=clip)
        mac_buf = io.BytesIO()
        Image.fromarray(enhanced.astype(np.uint8)).save(mac_buf, "JPEG",
                                                         quality=_PREVIEW_JPEG_QUALITY)

        hema = color_deconvolution_he(patch)
        hema_buf = io.BytesIO()
        Image.fromarray(
            np.stack([hema, hema, hema], axis=-1).astype(np.uint8)
        ).save(hema_buf, "JPEG", quality=_PREVIEW_JPEG_QUALITY)

        return {
            "status": "ok",
            "data": {
                "image_b64":   base64.b64encode(raw_buf.getvalue()).decode(),
                "clahe_b64": base64.b64encode(mac_buf.getvalue()).decode(),
                "hema_b64":    base64.b64encode(hema_buf.getvalue()).decode(),
                "method":      f"CLAHE(clip={clip}) + Ruifrok Hematoxylin",
                "patch_info":  f"x={x0},y={y0} size={x1-x0}×{y1-y0}",
                "clip_limit":  clip,
            },
        }
    except Exception as e:
        logger.error(f"前處理預覽失敗：{e}", exc_info=True)
        return {"status": "error", "message": "前處理預覽失敗，請查閱 log"}


@router.get("/preview")
async def get_preview(roi_name: Optional[str] = None):
    """回傳指定 ROI 分割遮罩疊加 H&E 的預覽圖（base64 JPEG）。"""
    import tifffile

    config     = load_config()
    paths      = config.get("paths", {})
    output_dir = resolve_path(paths.get("output_dir", "results/analysis"))
    mask_tif   = (
        config.get("segmentation", {}).get("output", {})
        .get("mask_tif_filename", "segmentation_masks.tif")
    )

    roi_base = output_dir / "roi"
    available: list[str] = []
    if roi_base.exists():
        for d in sorted(roi_base.iterdir()):
            if d.is_dir() and not d.name.startswith(".") and (d / mask_tif).exists():
                available.append(d.name)

    if not available:
        return {"status": "error", "message": "尚未執行分割，找不到遮罩檔案"}

    target    = roi_name if roi_name in available else available[0]
    roi_dir   = roi_base / target

    try:
        mask = tifffile.imread(str(roi_dir / mask_tif)).astype(np.int32)
        he_path = roi_dir / "he_crop.tif"
        # 灰底是「疊不上去」的降級結果，不是正常畫面 —— 必須說明原因，
        # 否則使用者只看到孤零零的輪廓，無從得知遮罩其實已經過期
        warning: Optional[str] = None
        if he_path.exists():
            he = tifffile.imread(str(he_path))
            if he.ndim == 2:
                he = np.stack([he, he, he], axis=-1)
            elif he.shape[-1] == 4:
                he = he[..., :3]
            if he.shape[:2] != mask.shape[:2]:
                warning = (
                    f"遮罩 {mask.shape[1]}×{mask.shape[0]} 與 H&E crop "
                    f"{he.shape[1]}×{he.shape[0]} 尺寸不符 —— 遮罩是舊 ROI 範圍留下的，"
                    f"請重新執行 Stage 1 分割。目前以灰底顯示，非真實組織。"
                )
                logger.warning(f"ROI '{target}' {warning}")
                he = np.full((*mask.shape, 3), 240, dtype=np.uint8)
        else:
            warning = "找不到 he_crop.tif（請先在 Stage 0 執行 ROI 裁切），目前以灰底顯示"
            logger.warning(f"ROI '{target}' {warning}")
            he = np.full((*mask.shape, 3), 240, dtype=np.uint8)

        from skimage.segmentation import find_boundaries
        overlay = he.copy()
        overlay[find_boundaries(mask, mode="thick")] = [50, 255, 80]
        n_cells = int(len(np.unique(mask)) - 1)

        h, w = overlay.shape[:2]
        scale = min(1200 / h, 1200 / w, 1.0)
        if scale < 1.0:
            overlay = cv2.resize(overlay, (int(w * scale), int(h * scale)))

        cv2.putText(overlay, f"n={n_cells} cells", (8, 28),
                    cv2.FONT_HERSHEY_SIMPLEX, 0.7, (255, 255, 0), 2)

        _, buf = cv2.imencode(
            ".jpg", cv2.cvtColor(overlay, cv2.COLOR_RGB2BGR),
            [cv2.IMWRITE_JPEG_QUALITY, 82],
        )

        return {
            "status": "ok",
            "data": {
                "image_b64":      base64.b64encode(buf.tobytes()).decode(),
                "flows_b64":      None,
                "roi":            target,
                "n_cells":        n_cells,
                "available_rois": available,
                "orig_w":         w,
                "orig_h":         h,
                "warning":        warning,
            },
        }
    except Exception as e:
        logger.error(f"取得分割預覽失敗：{e}", exc_info=True)
        return {"status": "error", "message": "取得分割預覽失敗，請查閱 log"}
