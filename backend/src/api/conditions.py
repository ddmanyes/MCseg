"""Stage 2.5: Proseg 參數條件測試 API"""
import asyncio
import logging
from typing import Any
from fastapi import APIRouter, BackgroundTasks
from pydantic import BaseModel

from backend.src.utils.config import load_config
from backend.src.utils.logging import set_current_stage

router = APIRouter()
logger = logging.getLogger("pipeline.api.conditions")

_task_status = {"status": "idle", "progress": 0.0, "message": "", "completed": 0, "total": 0}
_results: list[dict] = []


class ConditionGridRequest(BaseModel):
    max_dist: list[float] = [20, 40]
    compactness: list[float] = [0.03, 0.06]
    dilation: list[int] = [10, 20]
    roi_name: str = ""              # 空字串 = 使用 config 第一個 ROI
    quick_mode: bool = True


@router.get("/status")
async def get_status():
    return _task_status


@router.get("/results")
async def get_results():
    return {"status": "ok", "data": _results}


@router.get("/recommend")
async def recommend():
    if not _results:
        return {"status": "error", "message": "尚無結果，請先執行測試"}
    import math
    # 簡單啟發式：最大化 n_cells * median_genes
    best = max(_results, key=lambda r: r.get("n_cells", 0) * r.get("median_genes", 0))
    # 清除 NaN/inf 以符合 JSON 規範
    safe_best = {k: (None if isinstance(v, float) and (math.isnan(v) or math.isinf(v)) else v)
                 for k, v in best.items()}
    return {"status": "ok", "data": safe_best}


async def _run_conditions(config: dict, request: ConditionGridRequest):
    global _task_status, _results
    set_current_stage("conditions")
    from backend.src.proseg.condition_tester import ConditionTester
    tester = ConditionTester(config)
    grid = {
        "max_dist": request.max_dist,
        "compactness": request.compactness,
        "dilation": request.dilation,
    }
    total = len(request.max_dist) * len(request.compactness) * len(request.dilation)
    _task_status = {"status": "running", "progress": 0.0, "message": "開始條件測試...", "completed": 0, "total": total}
    _results = []

    def on_progress(completed: int, result: dict):
        _results.append(result)
        _task_status["completed"] = completed
        _task_status["progress"] = completed / total
        _task_status["message"] = f"完成 {completed}/{total} 條件"

    try:
        await asyncio.get_event_loop().run_in_executor(
            None, tester.run_grid, grid, request.roi_name, on_progress
        )
        _task_status["status"] = "done"
        _task_status["message"] = f"所有 {total} 個條件測試完成"
    except Exception as e:
        logger.error(f"條件測試失敗：{e}")
        _task_status = {"status": "error", "progress": 0.0, "message": str(e), "completed": 0, "total": 0}


@router.get("/thumbnail/{condition_idx}")
async def get_thumbnail(condition_idx: int):
    """回傳指定條件的 HE + 細胞輪廓疊圖縮圖（base64 JPEG）"""
    import base64
    config = load_config()
    from backend.src.utils.config import resolve_path
    cond_dir = resolve_path(config["paths"]["conditions_dir"]) / f"cond_{condition_idx:02d}"
    preview_path = cond_dir / "preview.jpg"
    if not preview_path.exists():
        return {"status": "error", "message": "縮圖尚未生成，請先執行條件測試"}
    img_b64 = base64.b64encode(preview_path.read_bytes()).decode()
    return {"status": "ok", "data": {"image_b64": img_b64}}


@router.get("/thumbnail_hd/{condition_idx}")
async def get_thumbnail_hd(condition_idx: int):
    """回傳高畫質 zoom 縮圖（200px 原圖裁切 → 800px，base64 JPEG）"""
    import base64
    config = load_config()
    from backend.src.utils.config import resolve_path
    cond_dir = resolve_path(config["paths"]["conditions_dir"]) / f"cond_{condition_idx:02d}"
    preview_path = cond_dir / "preview_hd.jpg"
    if not preview_path.exists():
        return {"status": "error", "message": "HD 縮圖尚未生成，請先執行條件測試"}
    img_b64 = base64.b64encode(preview_path.read_bytes()).decode()
    return {"status": "ok", "data": {"image_b64": img_b64}}


@router.post("/run")
async def run_conditions(request: ConditionGridRequest, background_tasks: BackgroundTasks):
    if _task_status["status"] == "running":
        return {"status": "error", "message": "任務執行中"}
    config = load_config()
    background_tasks.add_task(_run_conditions, config, request)
    return {"status": "ok", "message": "條件測試已啟動"}
