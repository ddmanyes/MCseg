"""Stage 2.5：Proseg RNA 重分配 API"""
import asyncio
import logging
from fastapi import APIRouter, BackgroundTasks
from pydantic import BaseModel
from typing import Optional

from backend.src.utils.config import load_config
from backend.src.utils.logging import set_current_stage

router = APIRouter()
logger = logging.getLogger("pipeline.api.proseg_rna")

_status = {"status": "idle", "progress": 0.0, "message": ""}
_lock   = asyncio.Lock()


class ProsegRNAParams(BaseModel):
    roi_name: Optional[str] = None   # None = 全部 ROI


@router.get("/status")
async def get_status():
    global _status
    # 記憶體 idle 時，查磁碟是否已有 proseg_cells.h5ad
    if _status["status"] == "idle":
        try:
            config = load_config()
            from backend.src.utils.config import resolve_path
            out_base = resolve_path(config["paths"]["output_dir"]) / "roi"
            rois = config.get("rois", [])
            if rois and any(
                (out_base / r.get("name", "") / "proseg_cells.h5ad").exists()
                for r in rois
            ):
                _status = {
                    "status": "done",
                    "progress": 1.0,
                    "message": "Proseg RNA 重分配已完成（從磁碟恢復）",
                }
        except Exception:
            pass
    return _status


async def _run_proseg_rna(config: dict, roi_name: Optional[str]):
    global _status
    set_current_stage("proseg_rna")
    _status = {
        "status": "running",
        "progress": 0.0,
        "message": "Proseg RNA 重分配執行中...",
    }
    try:
        from backend.src.proseg.runner import run_proseg_rna_pipeline
        await asyncio.get_running_loop().run_in_executor(
            None, run_proseg_rna_pipeline, config, roi_name
        )
        _status = {
            "status": "done",
            "progress": 1.0,
            "message": "Proseg RNA 重分配完成",
        }
    except Exception as e:
        logger.error(f"Proseg RNA 重分配失敗：{e}", exc_info=True)
        _status = {"status": "error", "progress": 0.0, "message": str(e)}


@router.post("/run")
async def run_proseg_rna(
    background_tasks: BackgroundTasks,
    params: Optional[ProsegRNAParams] = None,
):
    async with _lock:
        if _status["status"] == "running":
            return {"status": "error", "message": "任務執行中"}

    config   = load_config()
    roi_name = params.roi_name if params else None
    background_tasks.add_task(_run_proseg_rna, config, roi_name)
    return {"status": "ok", "message": "Proseg RNA 重分配已啟動"}


@router.get("/available_rois")
async def get_available_rois():
    """列出所有已有必要輸入（adata + mask）的 ROI 及 proseg_cells.h5ad 狀態"""
    try:
        config = load_config()
        from backend.src.utils.config import resolve_path
        out_base = resolve_path(config["paths"]["output_dir"]) / "roi"
        rois = config.get("rois", [])
        result = []
        for roi in rois:
            name = roi.get("name", "")
            try:
                has_adata  = (out_base / name / "adata_002um.h5ad").exists()
                has_mask   = (out_base / name / "segmentation_masks.npy").exists()
                has_proseg = (out_base / name / "proseg_cells.h5ad").exists()
                result.append({
                    "name":       name,
                    "has_adata":  has_adata,
                    "has_mask":   has_mask,
                    "has_proseg": has_proseg,
                })
            except (PermissionError, OSError) as e:
                logger.warning(f"跳過 {name}：{e}")
                continue
        return {"status": "ok", "data": result}
    except Exception as e:
        return {"status": "error", "message": str(e)}
