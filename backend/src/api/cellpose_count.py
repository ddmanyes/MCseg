"""Stage 2：Cellpose RNA 計數 API"""
import asyncio
import logging
from fastapi import APIRouter, BackgroundTasks
from pydantic import BaseModel
from typing import Optional

from backend.src.utils.config import load_config
from backend.src.utils.logging import set_current_stage

router = APIRouter()
logger = logging.getLogger("pipeline.api.cellpose_count")

_status = {"status": "idle", "progress": 0.0, "message": ""}
_lock   = asyncio.Lock()


class CountParams(BaseModel):
    roi_name: Optional[str] = None   # None = 全部 ROI


@router.get("/status")
async def get_status():
    global _status
    # 記憶體 idle 時，查磁碟是否已有 cellpose_cells.h5ad
    if _status["status"] == "idle":
        try:
            config = load_config()
            from backend.src.utils.config import resolve_path
            out_base = resolve_path(config["paths"]["output_dir"]) / "roi"
            rois = config.get("rois", [])
            if rois and any((out_base / r.get("name", "") / "cellpose_cells.h5ad").exists() for r in rois):
                _status = {"status": "done", "progress": 1.0, "message": "Count complete (restored from disk)"}
        except Exception:
            pass
    return _status


async def _run_count(config: dict, roi_name: Optional[str]):
    global _status
    set_current_stage("count")
    _status = {"status": "running", "progress": 0.0, "message": "分配 RNA 至 Cellpose 細胞..."}
    try:
        from backend.src.cellpose_counter.counter import run_counting_pipeline
        await asyncio.get_running_loop().run_in_executor(
            None, run_counting_pipeline, config, roi_name
        )
        _status = {"status": "done", "progress": 1.0, "message": "RNA counting complete"}
    except Exception as e:
        logger.error(f"Cellpose 計數失敗：{e}")
        _status = {"status": "error", "progress": 0.0, "message": str(e)}


@router.post("/run")
async def run_count(background_tasks: BackgroundTasks, params: Optional[CountParams] = None):
    async with _lock:
        if _status["status"] == "running":
            return {"status": "error", "message": "任務執行中"}
        _status["status"] = "running"

    config = load_config()
    roi_name = params.roi_name if params else None
    background_tasks.add_task(_run_count, config, roi_name)
    return {"status": "ok", "message": "RNA 計數已啟動"}


@router.get("/available_rois")
async def get_available_rois():
    """列出所有已有 cellpose_cells.h5ad 的 ROI"""
    try:
        config = load_config()
        from backend.src.utils.config import resolve_path
        out_base = resolve_path(config["paths"]["output_dir"]) / "roi"
        rois = config.get("rois", [])
        result = []
        for roi in rois:
            name = roi.get("name", "")
            has_mask  = (out_base / name / "segmentation_masks.npy").exists()
            has_count = (out_base / name / "cellpose_cells.h5ad").exists()
            result.append({
                "name":      name,
                "has_mask":  has_mask,
                "has_count": has_count,
            })
        return {"status": "ok", "data": result}
    except Exception as e:
        return {"status": "error", "message": str(e)}
