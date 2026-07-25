"""Stage 4: Browser 格式匯出 API（Pipeline 3 版本，使用 Cellpose mask 轉多邊形）"""
import asyncio
import functools
import logging
from fastapi import APIRouter, BackgroundTasks
from pydantic import BaseModel

from backend.src.utils.config import load_config, resolve_path
from backend.src.utils.logging import set_current_stage

router = APIRouter()
logger = logging.getLogger("pipeline.api.export")

_xenium_status = {"status": "idle", "progress": 0.0, "message": ""}
_loupe_status  = {"status": "idle", "progress": 0.0, "message": ""}
_result_status: dict = {"status": "idle", "progress": 0.0, "message": ""}
_result_images: dict = {}
_xenium_lock   = asyncio.Lock()
_loupe_lock    = asyncio.Lock()
_result_lock   = asyncio.Lock()


class ExportRequest(BaseModel):
    input_h5ad: str = ""   # 空字串 = 使用 config 預設輸出路徑
    output_dir: str = ""


@router.get("/status/xenium")
async def xenium_status():
    return _xenium_status


@router.get("/status/loupe")
async def loupe_status():
    return _loupe_status


# ──────────────────────────────────────────────────────────────────────────────
# 匯出背景任務（編排在 backend/src/export/jobs.py，這裡只負責狀態與執行緒調度）
# ──────────────────────────────────────────────────────────────────────────────

def _report_xenium(fraction: float, message: str) -> None:
    """給 export job 的進度回報 callback（在 executor thread 內被呼叫）。"""
    global _xenium_status
    _xenium_status = {"status": "running", "progress": fraction, "message": message}


def _report_loupe(fraction: float, message: str) -> None:
    global _loupe_status
    _loupe_status = {"status": "running", "progress": fraction, "message": message}


async def _run_xenium(config: dict, req: ExportRequest):
    """背景任務：把匯出丟進 thread pool，狀態回報寫進模組級 dict。"""
    global _xenium_status
    set_current_stage("export")
    _xenium_status = {"status": "running", "progress": 0.0, "message": "匯出至 Xenium Explorer..."}
    try:
        from backend.src.export.jobs import run_xenium_export

        result = await asyncio.get_running_loop().run_in_executor(
            None,
            functools.partial(
                run_xenium_export,
                config,
                req.input_h5ad,
                req.output_dir,
                progress=_report_xenium,
            ),
        )
        message = (
            f"Xenium 匯出完成（{len(result.output_dirs)} 個 ROI bundle）"
            if result.is_merged else "Xenium 匯出完成"
        )
        _xenium_status = {"status": "done", "progress": 1.0, "message": message}
    except Exception as e:
        logger.error(f"Xenium 匯出失敗：{e}", exc_info=True)
        _xenium_status = {"status": "error", "progress": 0.0, "message": "Xenium 匯出失敗，請查閱 log"}


async def _run_loupe(config: dict, req: ExportRequest):
    """背景任務：把匯出丟進 thread pool，狀態回報寫進模組級 dict。"""
    global _loupe_status
    set_current_stage("export")
    _loupe_status = {"status": "running", "progress": 0.0, "message": "匯出至 Loupe Browser..."}
    try:
        from backend.src.export.jobs import run_loupe_export

        await asyncio.get_running_loop().run_in_executor(
            None,
            functools.partial(
                run_loupe_export,
                config,
                req.input_h5ad,
                req.output_dir,
                progress=_report_loupe,
            ),
        )
        _loupe_status = {"status": "done", "progress": 1.0, "message": "Loupe 匯出完成"}
    except Exception as e:
        logger.error(f"Loupe 匯出失敗：{e}", exc_info=True)
        _loupe_status = {"status": "error", "progress": 0.0, "message": "Loupe 匯出失敗，請查閱 log"}


@router.post("/xenium")
async def export_xenium(req: ExportRequest, background_tasks: BackgroundTasks):
    global _xenium_status
    async with _xenium_lock:
        if _xenium_status["status"] == "running":
            return {"status": "error", "message": "任務執行中"}
        config = load_config()
        _xenium_status = {"status": "running", "progress": 0.0, "message": "準備匯出..."}
        background_tasks.add_task(_run_xenium, config, req)
    return {"status": "ok", "message": "Xenium 匯出已啟動"}


@router.post("/loupe")
async def export_loupe(req: ExportRequest, background_tasks: BackgroundTasks):
    global _loupe_status
    async with _loupe_lock:
        if _loupe_status["status"] == "running":
            return {"status": "error", "message": "任務執行中"}
        config = load_config()
        _loupe_status = {"status": "running", "progress": 0.0, "message": "準備匯出..."}
        background_tasks.add_task(_run_loupe, config, req)
    return {"status": "ok", "message": "Loupe 匯出已啟動"}


# ──────────────────────────────────────────────────────────────────────────────
# Result Visualizations（標註後視覺化）
# ──────────────────────────────────────────────────────────────────────────────

@router.get("/result_status")
async def get_result_status():
    return _result_status


@router.get("/result_images")
async def get_result_images():
    global _result_status, _result_images
    if not _result_images:
        try:
            fig_dir = resolve_path(load_config()["paths"]["figure_dir"])
            loaded: dict[str, str] = {}
            import base64
            # 固定圖
            for key, fname in [
                ("result_umap",     "result_umap.png"),
                ("result_dotplot",  "result_dotplot.png"),
                ("result_heatmap",  "result_heatmap.png"),
            ]:
                p = fig_dir / fname
                if p.exists():
                    loaded[key] = base64.b64encode(p.read_bytes()).decode()
            # 空間圖（支援單/多 ROI）
            import re as _re
            for p in sorted(fig_dir.glob("result_spatial*.png")):
                m = _re.match(r"(result_spatial(?:_filled)?(?:_.+)?)\.png$", p.name)
                if m:
                    loaded[m.group(1)] = base64.b64encode(p.read_bytes()).decode()
            if loaded:
                _result_images = loaded
                _result_status = {"status": "done", "progress": 1.0, "message": "已從磁碟載入結果圖"}
        except Exception as e:
            logger.warning(f"從磁碟載入結果圖失敗：{e}")
    if not _result_images:
        return {"status": "error", "message": "結果圖尚未產生，請先執行「生成結果圖」"}
    return {"status": "ok", "data": _result_images}


@router.post("/generate_result")
async def generate_result(background_tasks: BackgroundTasks):
    global _result_status
    async with _result_lock:
        if _result_status.get("status") == "running":
            return {"status": "running", "message": "已在執行中"}
        config = load_config()
        _result_status = {"status": "running", "progress": 0.0, "message": "生成結果視覺化中..."}
        background_tasks.add_task(_run_generate_result, config)
    return {"status": "started"}


async def _run_generate_result(config: dict):
    global _result_status, _result_images
    set_current_stage("export")
    try:
        from backend.src.analysis.pipeline import run_result_visualizations
        result = await asyncio.get_running_loop().run_in_executor(
            None, run_result_visualizations, config
        )
        _result_images = result
        _result_status = {
            "status": "done",
            "progress": 1.0,
            "message": f"結果圖生成完成（{len(result)} 張）",
        }
    except Exception as e:
        import traceback
        logger.error(f"結果視覺化失敗：{e!r}\n{traceback.format_exc()}")
        _result_status = {"status": "error", "progress": 0.0, "message": "結果視覺化失敗，請查閱 log"}
