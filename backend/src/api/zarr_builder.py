"""Stage 2: Zarr 建構 API"""
import asyncio
import logging
from fastapi import APIRouter, BackgroundTasks

from backend.src.utils.config import load_config
from backend.src.utils.logging import set_current_stage

router = APIRouter()
logger = logging.getLogger("pipeline.api.zarr")

_task_status = {"status": "idle", "progress": 0.0, "message": ""}


@router.get("/status")
async def get_status():
    return _task_status


async def _run_zarr(config: dict):
    global _task_status
    set_current_stage("zarr")
    _task_status = {"status": "running", "progress": 0.0, "message": "建構 Zarr..."}
    try:
        from backend.src.zarr_builder.builder import build_zarr
        await asyncio.get_event_loop().run_in_executor(None, build_zarr, config)
        _task_status = {"status": "done", "progress": 1.0, "message": "Zarr 建構完成"}
    except Exception as e:
        logger.error(f"Zarr 建構失敗：{e}")
        _task_status = {"status": "error", "progress": 0.0, "message": str(e)}


@router.post("/build")
async def build_zarr(background_tasks: BackgroundTasks):
    if _task_status["status"] == "running":
        return {"status": "error", "message": "任務執行中"}
    config = load_config()
    background_tasks.add_task(_run_zarr, config)
    return {"status": "ok", "message": "Zarr 建構已啟動"}
