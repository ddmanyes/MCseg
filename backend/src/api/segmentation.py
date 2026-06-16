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

from backend.src.utils.config import load_config, load_state, resolve_path, save_state
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
        save_state({"roi_seg_overrides": overrides})
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
async def run_full_segmentation(background_tasks: BackgroundTasks):
    """對 BTF 全圖執行 tiled MCseg v2 分割（MPS 安全模式）。"""
    async with _full_lock:
        if _full_status["status"] == "running":
            return {"status": "error", "message": "全圖分割任務執行中"}
        config = load_config()
        _full_status["status"]   = "running"
        _full_status["progress"] = 0.0
        _full_status["message"]  = "初始化..."
        background_tasks.add_task(_run_full_segmentation, config)
    return {"status": "ok", "message": "全圖分割已啟動"}


async def _run_full_segmentation(config: dict) -> None:
    global _full_status
    set_current_stage("segmentation")
    _full_status = {"status": "running", "progress": 0.0, "message": "讀取全圖影像..."}

    def _progress(p: float, msg: str) -> None:
        _full_status.update({"progress": p, "message": msg})

    try:
        import gc
        import numpy as np
        import tifffile
        import zarr as _zarr
        from backend.src.segmentation.cellpose_runner import run_tiled_mcseg_v2

        paths      = config.get("paths", {})
        output_dir = resolve_path(paths["output_dir"])
        btf_path   = paths.get("he_image", "")
        seg_cfg    = config.get("segmentation", {}).get("mcseg_v2", {})

        if not btf_path or not Path(btf_path).exists():
            raise FileNotFoundError(f"找不到 BTF/TIFF：{btf_path}  請在 config paths.he_image 指定")

        _progress(0.02, "讀取 BTF 全圖（tile-based）...")
        with tifffile.TiffFile(str(btf_path)) as tif:
            store = tif.aszarr()
            z = _zarr.open(store, mode="r")
            arr = z if not isinstance(z, _zarr.Group) else z[0]
            shape = arr.shape
            h_img, w_img = shape[0], shape[1]
            n_ch = shape[2] if len(shape) > 2 else 1
            estimated_gb = h_img * w_img * n_ch / 1024 ** 3
            if estimated_gb > 6.0:
                raise MemoryError(
                    f"全圖 {w_img}×{h_img}px ≈ {estimated_gb:.1f} GB，"
                    f"超過安全載入上限（6 GB）。"
                    f"請改用 ROI 裁切模式（Stage 0 + Stage 1）。"
                )
            img = np.array(arr)
        if img.ndim == 3 and img.shape[-1] == 4:
            img = img[..., :3]

        _progress(0.05, f"全圖尺寸 {img.shape[1]}×{img.shape[0]}px，開始 tiled 分割...")

        # MPS 安全設定：tile_size=1024, batch_size≤2
        tile_size = int(config.get("full_seg", {}).get("tile_size", 1024))
        overlap   = int(config.get("full_seg", {}).get("overlap", 128))
        seg_cfg_safe = dict(seg_cfg)
        seg_cfg_safe["batch_size"] = min(int(seg_cfg_safe.get("batch_size", 2)), 2)
        seg_cfg_safe["use_cpsam"] = False   # cpsam 在全圖模式太耗記憶體

        loop = asyncio.get_running_loop()
        import functools
        final_mask = await loop.run_in_executor(
            None,
            functools.partial(run_tiled_mcseg_v2, img, seg_cfg_safe,
                              tile_size, overlap, _progress),
        )
        del img
        gc.collect()

        out_path = output_dir / "full_image_segmentation_masks.npy"
        output_dir.mkdir(parents=True, exist_ok=True)
        np.save(str(out_path), final_mask)
        n_cells = int(len(np.unique(final_mask)) - 1)

        _full_status = {
            "status": "done", "progress": 1.0,
            "message": f"全圖分割完成：{n_cells:,} 個細胞  →  {out_path.name}",
            "n_cells": n_cells,
            "output": str(out_path),
        }
    except Exception as e:
        logger.error(f"全圖分割失敗：{e}", exc_info=True)
        if isinstance(e, MemoryError):
            # MemoryError 訊息含有尺寸資訊但不含路徑，可安全回傳
            safe_msg = str(e).split("，請改")[0] + "，請改用 ROI 裁切模式。"
        else:
            safe_msg = "全圖分割失敗，請查閱 log"
        _full_status = {"status": "error", "progress": 0.0, "message": safe_msg}


@router.get("/roi_overrides")
async def get_roi_overrides():
    return {"status": "ok", "data": load_state().get("roi_seg_overrides", {})}


@router.put("/roi_overrides")
async def put_roi_overrides(body: dict):
    from backend.src.segmentation.cellpose_runner import _ROI_OVERRIDE_FIELDS

    # 驗證 ROI 名稱：只允許已存在於 config 的 ROI，防止路徑穿越攻擊
    config = load_config()
    valid_names = {r.get("name", "") for r in config.get("rois", [])}
    invalid_names = [k for k in body if k not in valid_names]
    if invalid_names:
        return {"status": "error", "message": f"未知的 ROI 名稱：{invalid_names}"}

    # 驗證每個 ROI 的覆寫欄位名稱
    for roi_name, overrides in body.items():
        if not isinstance(overrides, dict):
            return {"status": "error", "message": f"ROI '{roi_name}' 的覆寫值必須是 dict"}
        invalid_fields = [k for k in overrides if k not in _ROI_OVERRIDE_FIELDS]
        if invalid_fields:
            return {"status": "error", "message": f"ROI '{roi_name}' 包含未知欄位：{invalid_fields}"}

    save_state({"roi_seg_overrides": body})
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
        if he_path.exists():
            he = tifffile.imread(str(he_path))
            if he.ndim == 2:
                he = np.stack([he, he, he], axis=-1)
            elif he.shape[-1] == 4:
                he = he[..., :3]
            if he.shape[:2] != mask.shape[:2]:
                he = np.full((*mask.shape, 3), 240, dtype=np.uint8)
        else:
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
            },
        }
    except Exception as e:
        logger.error(f"取得分割預覽失敗：{e}", exc_info=True)
        return {"status": "error", "message": "取得分割預覽失敗，請查閱 log"}
