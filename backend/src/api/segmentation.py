"""Stage 1: 細胞分割 API"""
import asyncio
import base64
import logging
from fastapi import APIRouter, BackgroundTasks

from backend.src.utils.config import load_config
from backend.src.utils.logging import set_current_stage

router = APIRouter()
logger = logging.getLogger("pipeline.api.segmentation")

_task_status = {"status": "idle", "progress": 0.0, "message": ""}


@router.get("/status")
async def get_status():
    return _task_status


async def _run_segmentation(config: dict):
    global _task_status
    set_current_stage("segmentation")
    _task_status = {"status": "running", "progress": 0.0, "message": "啟動 Cellpose..."}
    try:
        from backend.src.segmentation.cellpose_runner import run_segmentation
        await asyncio.get_event_loop().run_in_executor(None, run_segmentation, config)
        _task_status = {"status": "done", "progress": 1.0, "message": "分割完成"}
    except Exception as e:
        logger.error(f"分割失敗：{e}")
        _task_status = {"status": "error", "progress": 0.0, "message": str(e)}


@router.post("/run")
async def run_segmentation(background_tasks: BackgroundTasks):
    if _task_status["status"] == "running":
        return {"status": "error", "message": "任務執行中"}
    config = load_config()
    background_tasks.add_task(_run_segmentation, config)
    return {"status": "ok", "message": "細胞分割已啟動"}


@router.get("/preview")
async def get_preview():
    """回傳分割遮罩預覽圖（base64 PNG）"""
    from pathlib import Path
    import numpy as np
    config = load_config()
    mask_path = Path(config["paths"]["masks_dir"]) / config["segmentation"]["output"]["mask_tif_filename"]
    if not mask_path.exists():
        return {"status": "error", "message": "尚未執行分割"}
    try:
        import tifffile
        import cv2
        mask = tifffile.imread(str(mask_path))
        # 將 label map 轉為 RGB 預覽
        from skimage.color import label2rgb
        preview = (label2rgb(mask, bg_label=0) * 255).astype(np.uint8)
        preview_small = cv2.resize(preview, (800, int(800 * preview.shape[0] / preview.shape[1])))
        _, buf = cv2.imencode(".jpg", cv2.cvtColor(preview_small, cv2.COLOR_RGB2BGR), [cv2.IMWRITE_JPEG_QUALITY, 80])
        img_b64 = base64.b64encode(buf.tobytes()).decode()
        return {"status": "ok", "data": {"image_b64": img_b64}}
    except Exception as e:
        return {"status": "error", "message": str(e)}
