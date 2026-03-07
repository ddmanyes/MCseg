"""Stage 4: 下游分析 API（支援三步驟分段執行 + 舊版整合執行）"""
import asyncio
import base64
import logging
from pathlib import Path
from fastapi import APIRouter, BackgroundTasks
from fastapi.responses import FileResponse
from pydantic import BaseModel
from typing import Optional

from backend.src.utils.config import load_config, save_state
from backend.src.utils.logging import set_current_stage

router = APIRouter()
logger = logging.getLogger("pipeline.api.analysis")

# ─────────────────────── Pydantic Models ───────────────────────────

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


class QCParams(BaseModel):
    min_genes: Optional[int] = None
    min_counts: Optional[int] = None
    max_genes: Optional[int] = None
    min_cells: Optional[int] = None
    max_pct_mito: Optional[float] = None
    n_top_genes: Optional[int] = None
    n_pcs: Optional[int] = None
    merge_rois: Optional[bool] = None   # None = 使用 config 預設
    roi_name: Optional[str] = None       # 單一 ROI 模式時指定（None = 取第一個）


class UMAPExploreParams(BaseModel):
    n_pcs: Optional[int] = None
    n_neighbors: Optional[int] = None
    min_dist: Optional[float] = None
    resolutions: list[float] = [0.3, 0.5, 0.8]


class HeatmapParams(BaseModel):
    resolution: float
    n_top_genes: int = 20       # dotplot：每 cluster 幾個 marker gene
    n_heatmap_genes: int = 50   # heatmap：顯示方差最高的幾個 HVGs


class AnnotateParams(BaseModel):
    resolution: float
    model_name: str = "Human_Colorectal_Cancer.pkl"


class ApplyLabelsParams(BaseModel):
    resolution: float
    labels: dict[str, str]     # {cluster_id: cell_type_name}


# ─────────────────────── 舊版整合 status ──────────────────────────

_task_status = {"status": "idle", "progress": 0.0, "message": ""}
_task_lock = asyncio.Lock()

# ────────────────── 三步驟各自 status / lock ───────────────────────

_qc_status    = {"status": "idle", "progress": 0.0, "message": ""}
_umap_status  = {"status": "idle", "progress": 0.0, "message": ""}
_heat_status  = {"status": "idle", "progress": 0.0, "message": ""}
_annot_status = {"status": "idle", "progress": 0.0, "message": ""}

_qc_lock    = asyncio.Lock()
_umap_lock  = asyncio.Lock()
_heat_lock  = asyncio.Lock()
_annot_lock = asyncio.Lock()

# 暫存最近一次的圖表資料（避免每次查詢都讀磁碟）
_qc_images:    dict[str, str] = {}
_umap_images:  dict[str, str] = {}
_heatmap_images: dict[str, str] = {}
_annot_suggestions: dict[str, str] = {}   # 最近一次 CellTypist 建議


# ──────────────────── 舊版整合執行 /run ───────────────────────────

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
    async with _task_lock:
        if _task_status["status"] == "running":
            return {"status": "error", "message": "任務執行中"}

    config = load_config()
    if params:
        _patch_config_from_analysis_params(config, params)
        save_state({"analysis": config.get("analysis", {})})

    background_tasks.add_task(_run_analysis, config)
    return {"status": "ok", "message": "分析已啟動"}


@router.get("/umap")
async def get_umap():
    """回傳最新的 UMAP 圖（base64 PNG）— 舊版相容"""
    config = load_config()
    from backend.src.utils.config import resolve_path
    fig_dir = resolve_path(config["paths"]["figure_dir"])
    umap_path = fig_dir / "umap.png"
    if not umap_path.exists():
        return {"status": "error", "message": "UMAP 圖尚未產生"}
    img_b64 = base64.b64encode(umap_path.read_bytes()).decode()
    return {"status": "ok", "data": {"image_b64": img_b64}}


# ─────────────────────── helpers ──────────────────────────────────

def _get_fig_dir() -> Path:
    from backend.src.utils.config import resolve_path
    config = load_config()
    return resolve_path(config["paths"]["figure_dir"])


def _get_output_dir() -> Path:
    from backend.src.utils.config import resolve_path
    config = load_config()
    return resolve_path(config["paths"]["output_dir"])


def _load_disk_images(keys_paths: list[tuple[str, Path]]) -> dict[str, str]:
    """從磁碟載入存在的圖檔，回傳 {key: base64}"""
    result: dict[str, str] = {}
    for key, path in keys_paths:
        if path.exists():
            result[key] = base64.b64encode(path.read_bytes()).decode()
    return result


# ─────────────────────── Step 1: QC ───────────────────────────────

@router.get("/qc_status")
async def get_qc_status():
    global _qc_status
    # 記憶體為 idle 時，查磁碟判斷是否已有產出
    if _qc_status["status"] == "idle":
        try:
            output_dir = _get_output_dir()
            if (output_dir / "qc_preprocessed.h5ad").exists():
                _qc_status = {"status": "done", "progress": 1.0, "message": "QC 完成（從磁碟恢復）"}
        except Exception:
            pass
    return _qc_status


@router.get("/qc_images")
async def get_qc_images():
    """回傳最新的 QC 圖（violin / scatter / elbow），記憶體空時自動查磁碟"""
    global _qc_images
    if not _qc_images:
        try:
            fig_dir = _get_fig_dir()
            _qc_images = _load_disk_images([
                ("violin",   fig_dir / "qc_violin.png"),
                ("scatter",  fig_dir / "qc_scatter.png"),
                ("elbow",    fig_dir / "pca_elbow.png"),
                ("pre_qc",   fig_dir / "overlay_pre_qc.png"),
                ("post_qc",  fig_dir / "overlay_post_qc.png"),
            ])
        except Exception:
            pass
    if not _qc_images:
        return {"status": "error", "message": "QC 圖尚未產生"}
    return {"status": "ok", "data": _qc_images}


@router.get("/overlay_hd/{name}")
async def download_overlay_hd(name: str):
    """下載 HD 疊圖（300 DPI）。name = pre_qc | post_qc"""
    allowed = {"pre_qc": "overlay_pre_qc_hd.png", "post_qc": "overlay_post_qc_hd.png"}
    if name not in allowed:
        from fastapi import HTTPException
        raise HTTPException(status_code=404, detail=f"未知圖名：{name}")
    fig_dir = _get_fig_dir()
    path = fig_dir / allowed[name]
    if not path.exists():
        from fastapi import HTTPException
        raise HTTPException(status_code=404, detail="HD 疊圖尚未產生，請先執行 QC")
    return FileResponse(str(path), media_type="image/png", filename=allowed[name])


async def _run_qc(config: dict):
    global _qc_status, _qc_images
    set_current_stage("analysis")
    _qc_status = {"status": "running", "progress": 0.0, "message": "執行 QC 前處理..."}
    try:
        from backend.src.analysis.pipeline import run_qc_step
        result = await asyncio.get_event_loop().run_in_executor(None, run_qc_step, config)
        _qc_images = result
        _qc_status = {"status": "done", "progress": 1.0, "message": f"QC 完成，已產生 {len(result)} 張圖表"}
    except Exception as e:
        logger.error(f"QC Step 失敗：{e}")
        _qc_status = {"status": "error", "progress": 0.0, "message": str(e)}


@router.post("/run_qc")
async def run_qc(background_tasks: BackgroundTasks, params: Optional[QCParams] = None):
    async with _qc_lock:
        if _qc_status["status"] == "running":
            return {"status": "error", "message": "QC 任務執行中"}

    config = load_config()
    if params:
        _patch_config_from_qc_params(config, params)
        save_state({"analysis": config.get("analysis", {})})
        # 分析來源選擇（in-memory 覆寫，不寫回 state）
        if params.merge_rois is not None:
            config.setdefault("analysis", {})["merge_rois"] = params.merge_rois
        if params.roi_name is not None and not params.merge_rois:
            # 將指定 ROI 移至列表第一位，供 run_qc_step 的 rois[0] 讀取
            rois = config.get("rois", [])
            target = next((r for r in rois if r.get("name") == params.roi_name), None)
            if target:
                config["rois"] = [target] + [r for r in rois if r.get("name") != params.roi_name]

    background_tasks.add_task(_run_qc, config)
    return {"status": "ok", "message": "QC 前處理已啟動"}


@router.get("/available_rois")
async def get_available_rois():
    """列出所有已有 proseg_cells.h5ad 的 ROI（供 Stage 4 來源選擇）"""
    try:
        config = load_config()
        from backend.src.utils.config import resolve_path
        out_base = resolve_path(config["paths"]["output_dir"]) / "roi"
        rois = config.get("rois", [])
        result = []
        for roi in rois:
            name = roi.get("name", "")
            h5ad = out_base / name / "proseg_cells.h5ad"
            result.append({"name": name, "available": h5ad.exists()})
        return {"status": "ok", "data": result}
    except Exception as e:
        return {"status": "error", "message": str(e)}


# ─────────────────────── Step 2: UMAP ─────────────────────────────

@router.get("/umap_status")
async def get_umap_status():
    global _umap_status
    if _umap_status["status"] == "idle":
        try:
            output_dir = _get_output_dir()
            if (output_dir / "umap_computed.h5ad").exists():
                _umap_status = {"status": "done", "progress": 1.0, "message": "UMAP 完成（從磁碟恢復）"}
        except Exception:
            pass
    return _umap_status


@router.get("/umap_images")
async def get_umap_images():
    """回傳 UMAP 圖（各 resolution 一張 + grid），記憶體空時自動查磁碟"""
    global _umap_images
    if not _umap_images:
        try:
            fig_dir = _get_fig_dir()
            # 掃描所有 umap_res*.png
            disk_imgs = _load_disk_images([("grid", fig_dir / "umap_grid.png")])
            for p in sorted(fig_dir.glob("umap_res*.png")):
                key = p.stem.replace("umap_res", "")  # e.g. "0.5"
                disk_imgs[key] = base64.b64encode(p.read_bytes()).decode()
            _umap_images = disk_imgs
        except Exception:
            pass
    if not _umap_images:
        return {"status": "error", "message": "UMAP 圖尚未產生"}
    return {"status": "ok", "data": _umap_images}


async def _run_umap(config: dict, p: UMAPExploreParams):
    global _umap_status, _umap_images
    set_current_stage("analysis")
    _umap_status = {"status": "running", "progress": 0.0, "message": "計算 UMAP..."}
    try:
        from backend.src.analysis.pipeline import run_umap_step
        clus_cfg = config.get("analysis", {}).get("clustering", {})
        result = await asyncio.get_event_loop().run_in_executor(
            None,
            run_umap_step,
            config,
            p.resolutions,
            p.n_pcs or clus_cfg.get("n_pcs", 30),
            p.n_neighbors or clus_cfg.get("n_neighbors", 15),
            p.min_dist if p.min_dist is not None else clus_cfg.get("min_dist", 0.3),
        )
        _umap_images = result
        _umap_status = {
            "status": "done",
            "progress": 1.0,
            "message": f"UMAP 完成，{len(p.resolutions)} 個解析度",
            "resolutions": p.resolutions,
        }
    except Exception as e:
        logger.error(f"UMAP Step 失敗：{e}")
        _umap_status = {"status": "error", "progress": 0.0, "message": str(e)}


@router.post("/run_umap")
async def run_umap_explore(background_tasks: BackgroundTasks, params: Optional[UMAPExploreParams] = None):
    async with _umap_lock:
        if _umap_status["status"] == "running":
            return {"status": "error", "message": "UMAP 任務執行中"}

    p = params or UMAPExploreParams()
    config = load_config()

    # 更新 UMAP 參數至設定
    clus = config.setdefault("analysis", {}).setdefault("clustering", {})
    if p.n_pcs is not None:       clus["n_pcs"] = p.n_pcs
    if p.n_neighbors is not None: clus["n_neighbors"] = p.n_neighbors
    if p.min_dist is not None:    clus["min_dist"] = p.min_dist
    save_state({"analysis": {"clustering": clus}})

    background_tasks.add_task(_run_umap, config, p)
    return {"status": "ok", "message": "UMAP 探索已啟動"}


# ─────────────────────── Step 3: Heatmap ──────────────────────────

@router.get("/heatmap_status")
async def get_heatmap_status():
    global _heat_status
    if _heat_status["status"] == "idle":
        try:
            fig_dir = _get_fig_dir()
            if (fig_dir / "heatmap.png").exists() or (fig_dir / "dotplot.png").exists():
                _heat_status = {"status": "done", "progress": 1.0, "message": "圖表完成（從磁碟恢復）"}
        except Exception:
            pass
    return _heat_status


@router.get("/heatmap")
async def get_heatmap():
    global _heatmap_images
    if not _heatmap_images:
        try:
            fig_dir = _get_fig_dir()
            for key, fname in [("heatmap", "heatmap.png"), ("dotplot", "dotplot.png")]:
                p = fig_dir / fname
                if p.exists():
                    _heatmap_images[key] = base64.b64encode(p.read_bytes()).decode()
        except Exception:
            pass
    if not _heatmap_images:
        return {"status": "error", "message": "圖表尚未產生"}
    return {"status": "ok", "data": _heatmap_images}


async def _run_heatmap(config: dict, p: HeatmapParams):
    global _heat_status, _heatmap_images
    set_current_stage("analysis")
    _heatmap_images = {}  # 清空舊快取，避免前端在新圖完成前取到舊圖
    _heat_status = {"status": "running", "progress": 0.0, "message": f"產生熱圖 + 點圖 (res={p.resolution})..."}
    try:
        from backend.src.analysis.pipeline import run_heatmap_step
        result = await asyncio.get_event_loop().run_in_executor(
            None, run_heatmap_step, config, p.resolution, p.n_top_genes, p.n_heatmap_genes
        )
        _heatmap_images = result  # dict: {"heatmap": b64, "dotplot": b64}
        _heat_status = {"status": "done", "progress": 1.0, "message": "熱圖 + 點圖完成"}
    except Exception as e:
        logger.error(f"Heatmap Step 失敗：{e}")
        _heat_status = {"status": "error", "progress": 0.0, "message": str(e)}


@router.post("/run_heatmap")
async def run_heatmap(background_tasks: BackgroundTasks, params: HeatmapParams):
    async with _heat_lock:
        if _heat_status["status"] == "running":
            return {"status": "error", "message": "熱圖任務執行中"}

    config = load_config()
    background_tasks.add_task(_run_heatmap, config, params)
    return {"status": "ok", "message": "熱圖產生已啟動"}


# ─────────────────────── Config Helpers ───────────────────────────

def _patch_config_from_analysis_params(config: dict, params: AnalysisParams):
    cfg = config.setdefault("analysis", {})
    pre = cfg.setdefault("preprocessing", {})
    cellular = pre.setdefault("cellular", {})
    if params.min_genes is not None:     cellular["min_genes"] = params.min_genes
    if params.min_counts is not None:    cellular["min_counts"] = params.min_counts
    if params.max_genes is not None:     cellular["max_genes"] = params.max_genes
    if params.min_cells is not None:     cellular["min_cells"] = params.min_cells
    if params.max_pct_mito is not None:  cellular["max_pct_mito"] = params.max_pct_mito
    norm = pre.setdefault("normalization", {})
    if params.target_sum is not None:    norm["target_sum"] = params.target_sum
    hvg = pre.setdefault("hvg", {})
    if params.n_top_genes is not None:   hvg["n_top_genes"] = params.n_top_genes
    clus = cfg.setdefault("clustering", {})
    if params.n_pcs is not None:         clus["n_pcs"] = params.n_pcs
    if params.n_neighbors is not None:   clus["n_neighbors"] = params.n_neighbors
    if params.resolution is not None:    clus["resolution"] = params.resolution
    if params.min_dist is not None:      clus["min_dist"] = params.min_dist


def _patch_config_from_qc_params(config: dict, params: QCParams):
    cfg = config.setdefault("analysis", {})
    pre = cfg.setdefault("preprocessing", {})
    cellular = pre.setdefault("cellular", {})
    if params.min_genes is not None:     cellular["min_genes"] = params.min_genes
    if params.min_counts is not None:    cellular["min_counts"] = params.min_counts
    if params.max_genes is not None:     cellular["max_genes"] = params.max_genes
    if params.min_cells is not None:     cellular["min_cells"] = params.min_cells
    if params.max_pct_mito is not None:  cellular["max_pct_mito"] = params.max_pct_mito
    hvg = pre.setdefault("hvg", {})
    if params.n_top_genes is not None:   hvg["n_top_genes"] = params.n_top_genes
    clus = cfg.setdefault("clustering", {})
    if params.n_pcs is not None:         clus["n_pcs"] = params.n_pcs


# ─────────────────── Step 3: CellTypist 標註 endpoints ────────────

@router.get("/cluster_info")
async def get_cluster_info(resolution: float):
    """取得指定 resolution 的 cluster ID 列表，以及已套用的標籤（若有）。"""
    config = load_config()
    try:
        from backend.src.analysis.pipeline import get_cluster_ids
        ids, existing_labels = get_cluster_ids(config, resolution)
        return {"status": "ok", "data": {"cluster_ids": ids, "existing_labels": existing_labels}}
    except Exception as e:
        return {"status": "error", "message": str(e)}


@router.get("/celltypist_models")
async def get_celltypist_models():
    """回傳支援的 CellTypist 模型清單。"""
    from backend.src.analysis.pipeline import CELLTYPIST_MODELS
    return {"status": "ok", "data": CELLTYPIST_MODELS}


async def _run_annotate(config: dict, p: AnnotateParams):
    global _annot_status, _annot_suggestions
    set_current_stage("analysis")
    _annot_suggestions = {}
    _annot_status = {
        "status": "running", "progress": 0.0,
        "message": f"CellTypist 標註中 (res={p.resolution}, model={p.model_name})...",
    }
    try:
        from backend.src.analysis.pipeline import run_celltypist_annotation
        result = await asyncio.get_event_loop().run_in_executor(
            None, run_celltypist_annotation, config, p.resolution, p.model_name
        )
        _annot_suggestions = result
        _annot_status = {
            "status": "done", "progress": 1.0,
            "message": f"CellTypist 標註完成（{len(result)} 個 cluster）",
        }
    except Exception as e:
        logger.error(f"CellTypist 標註失敗：{e}")
        _annot_status = {"status": "error", "progress": 0.0, "message": str(e)}


@router.post("/annotate")
async def run_annotation(background_tasks: BackgroundTasks, params: AnnotateParams):
    async with _annot_lock:
        if _annot_status["status"] == "running":
            return {"status": "error", "message": "標註任務執行中"}
    config = load_config()
    background_tasks.add_task(_run_annotate, config, params)
    return {"status": "ok", "message": "CellTypist 標註已啟動"}


@router.get("/annotate_status")
async def get_annotate_status():
    if _annot_status.get("status") == "done" and _annot_suggestions:
        return {**_annot_status, "suggestions": _annot_suggestions}
    return _annot_status


@router.get("/annotate_suggestions")
async def get_annotate_suggestions():
    return {"status": "ok", "data": _annot_suggestions}


@router.post("/apply_labels")
async def apply_labels_endpoint(params: ApplyLabelsParams):
    """同步套用標籤（操作快，不需要背景任務）。"""
    config = load_config()
    try:
        from backend.src.analysis.pipeline import apply_cluster_labels
        apply_cluster_labels(config, params.resolution, params.labels)
        return {"status": "ok", "message": f"已套用 {len(params.labels)} 個 cluster 標籤"}
    except Exception as e:
        logger.error(f"套用標籤失敗：{e}")
        return {"status": "error", "message": str(e)}
