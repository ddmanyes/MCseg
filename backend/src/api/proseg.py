"""Stage 3: Proseg 完整執行 API"""
import asyncio
import logging
from fastapi import APIRouter, BackgroundTasks

from backend.src.utils.config import load_config
from backend.src.utils.logging import set_current_stage

router = APIRouter()
logger = logging.getLogger("pipeline.api.proseg")

_task_status = {"status": "idle", "progress": 0.0, "message": ""}


@router.get("/status")
async def get_status():
    return _task_status


async def _run_proseg(config: dict):
    global _task_status
    set_current_stage("proseg")
    _task_status = {"status": "running", "progress": 0.0, "message": "啟動 Proseg..."}
    try:
        from backend.src.proseg.pipeline import ProsegPipeline
        pipeline = ProsegPipeline(config)
        await asyncio.get_event_loop().run_in_executor(None, pipeline.run_full)
        _task_status = {"status": "done", "progress": 1.0, "message": "Proseg 完成"}
    except Exception as e:
        logger.error(f"Proseg 失敗：{e}")
        _task_status = {"status": "error", "progress": 0.0, "message": str(e)}


@router.post("/run")
async def run_proseg(background_tasks: BackgroundTasks):
    if _task_status["status"] == "running":
        return {"status": "error", "message": "任務執行中"}
    config = load_config()
    background_tasks.add_task(_run_proseg, config)
    return {"status": "ok", "message": "Proseg 已啟動"}
