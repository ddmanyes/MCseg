"""Stage 4: 下游分析 API"""
import asyncio
import base64
import logging
from pathlib import Path
from fastapi import APIRouter, BackgroundTasks
from pydantic import BaseModel
from typing import Optional

from backend.src.utils.config import load_config, save_config
from backend.src.utils.logging import set_current_stage

router = APIRouter()
logger = logging.getLogger("pipeline.api.analysis")

class AnalysisParams(BaseModel):
    min_genes: Optional[int] = None
    min_counts: Optional[int] = None
    max_genes: Optional[int] = None
    min_cells: Optional[int] = None
    max_pct_mito: Optional[float] = None
    target_sum: Optional[int] = None
    n_top_genes: Optional[int] = None
    n_pcs: Optional[int] = None
    n_neighbors: Optional[int] = None
    resolution: Optional[float] = None
    min_dist: Optional[float] = None

_task_status = {"status": "idle", "progress": 0.0, "message": ""}


@router.get("/status")
async def get_status():
    return _task_status


async def _run_analysis(config: dict):
    global _task_status
    set_current_stage("analysis")
    _task_status = {"status": "running", "progress": 0.0, "message": "執行聚類分析..."}
    try:
        from backend.src.analysis.pipeline import run_analysis_pipeline
        await asyncio.get_event_loop().run_in_executor(None, run_analysis_pipeline, config)
        _task_status = {"status": "done", "progress": 1.0, "message": "分析完成"}
    except Exception as e:
        logger.error(f"分析失敗：{e}")
        _task_status = {"status": "error", "progress": 0.0, "message": str(e)}


@router.post("/run")
async def run_analysis(background_tasks: BackgroundTasks, params: Optional[AnalysisParams] = None):
    if _task_status["status"] == "running":
        return {"status": "error", "message": "任務執行中"}
    
    config = load_config()
    
    if params:
        # Patch configuration
        cfg = config.setdefault("analysis", {})
        
        pre = cfg.setdefault("preprocessing", {})
        cellular = pre.setdefault("cellular", {})
        if params.min_genes is not None: cellular["min_genes"] = params.min_genes
        if params.min_counts is not None: cellular["min_counts"] = params.min_counts
        if params.max_genes is not None: cellular["max_genes"] = params.max_genes
        if params.min_cells is not None: cellular["min_cells"] = params.min_cells
        if params.max_pct_mito is not None: cellular["max_pct_mito"] = params.max_pct_mito
        
        norm = pre.setdefault("normalization", {})
        if params.target_sum is not None: norm["target_sum"] = params.target_sum
        
        hvg = pre.setdefault("hvg", {})
        if params.n_top_genes is not None: hvg["n_top_genes"] = params.n_top_genes
        
        clus = cfg.setdefault("clustering", {})
        if params.n_pcs is not None: clus["n_pcs"] = params.n_pcs
        if params.n_neighbors is not None: clus["n_neighbors"] = params.n_neighbors
        if params.resolution is not None: clus["resolution"] = params.resolution
        if params.min_dist is not None: clus["min_dist"] = params.min_dist

        save_config(config)

    background_tasks.add_task(_run_analysis, config)
    return {"status": "ok", "message": "分析已啟動"}


@router.get("/umap")
async def get_umap():
    """回傳最新的 UMAP 圖（base64 PNG）"""
    config = load_config()
    fig_dir = Path(config["paths"]["figure_dir"])
    umap_path = fig_dir / "umap.png"
    if not umap_path.exists():
        return {"status": "error", "message": "UMAP 圖尚未產生"}
    img_b64 = base64.b64encode(umap_path.read_bytes()).decode()
    return {"status": "ok", "data": {"image_b64": img_b64}}
