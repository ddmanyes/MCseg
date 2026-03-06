"""Stage 5: Browser 格式匯出 API"""
import asyncio
import logging
from pathlib import Path
from fastapi import APIRouter, BackgroundTasks
from pydantic import BaseModel

from backend.src.utils.config import load_config, resolve_path
from backend.src.utils.logging import set_current_stage

router = APIRouter()
logger = logging.getLogger("pipeline.api.export")

_xenium_status = {"status": "idle", "progress": 0.0, "message": ""}
_loupe_status  = {"status": "idle", "progress": 0.0, "message": ""}


class ExportRequest(BaseModel):
    input_h5ad: str = ""   # 空字串 = 使用 config 預設輸出路徑
    output_dir: str = ""


@router.get("/status/xenium")
async def xenium_status():
    return _xenium_status


@router.get("/status/loupe")
async def loupe_status():
    return _loupe_status


async def _run_xenium(config: dict, req: ExportRequest):
    global _xenium_status
    set_current_stage("export")
    _xenium_status = {"status": "running", "progress": 0.0, "message": "匯出至 Xenium Explorer..."}
    try:
        from backend.src.export.xenium_exporter import XeniumExporter, generate_combined_geojson
        import json

        paths = config.get("paths", {})
        rois = config.get("rois", [{}])
        roi_name = rois[0].get("name", "") if rois else ""

        zarr_dir = resolve_path(paths.get("zarr_dir", "results/zarr"))
        output_dir_base = resolve_path(paths.get("output_dir", "results/analysis"))
        export_dir = resolve_path(paths.get("export_dir", "results/export"))

        zarr_path = zarr_dir / roi_name / "proseg_integrated.zarr" if roi_name else zarr_dir / "proseg_integrated.zarr"

        # ROI 輸出目錄（實際存放資料的位置）
        roi_out_dir = output_dir_base / "roi" / roi_name if roi_name else output_dir_base
        tile_proseg_dir = roi_out_dir / "proseg_tiles"

        # 合併多邊形 GeoJSON（若尚未生成則從 tile 重建）
        combined_poly_path = roi_out_dir / "combined_proseg_results_qc.json"
        if not combined_poly_path.exists() and tile_proseg_dir.exists() and zarr_path.exists():
            logger.info("未找到合併 GeoJSON，從 tile 自動生成...")
            combined = generate_combined_geojson(tile_proseg_dir, zarr_path, config)
            combined_poly_path.parent.mkdir(parents=True, exist_ok=True)
            with open(combined_poly_path, "w", encoding="utf-8") as f:
                json.dump(combined, f)
            logger.info(f"合併 GeoJSON 已儲存：{combined_poly_path}（{len(combined['features'])} 個多邊形）")

        # 決定 h5ad 路徑：優先使用請求參數，其次搜尋標準輸出路徑
        if req.input_h5ad:
            h5ad_path = Path(req.input_h5ad)
        else:
            h5ad_path = None
            for candidate in ["clustered_final.h5ad", "umap_computed.h5ad", "qc_preprocessed.h5ad", "proseg_cells.h5ad"]:
                for search_dir in [output_dir_base, roi_out_dir]:
                    p = search_dir / candidate
                    if p.exists():
                        h5ad_path = p
                        break
                if h5ad_path:
                    break
            if h5ad_path is None:
                raise FileNotFoundError(
                    f"找不到分析結果 h5ad，請先執行 Stage 4 或指定路徑。"
                    f"搜尋位置：{output_dir_base}、{roi_out_dir}"
                )

        out_dir = Path(req.output_dir) if req.output_dir else export_dir / "xenium"

        exporter = XeniumExporter(
            zarr_path=zarr_path if zarr_path.exists() else None,
            poly_json_path=combined_poly_path if combined_poly_path.exists() else None,
            transcripts_csv_path=None,  # tile-level transcripts 暫不合並
        )
        await asyncio.get_event_loop().run_in_executor(
            None, exporter.export, h5ad_path, out_dir,
        )
        _xenium_status = {"status": "done", "progress": 1.0, "message": "Xenium 匯出完成"}
    except Exception as e:
        logger.error(f"Xenium 匯出失敗：{e}")
        _xenium_status = {"status": "error", "progress": 0.0, "message": str(e)}


async def _run_loupe(config: dict, req: ExportRequest):
    global _loupe_status
    set_current_stage("export")
    _loupe_status = {"status": "running", "progress": 0.0, "message": "匯出至 Loupe Browser..."}
    try:
        from backend.src.export.loupe_exporter import LoupeExporter

        paths = config.get("paths", {})
        rois = config.get("rois", [{}])
        roi_name = rois[0].get("name", "") if rois else ""

        proseg_dir = resolve_path(paths.get("proseg_dir", "results/proseg"))
        output_dir_base = resolve_path(paths.get("output_dir", "results/analysis"))
        export_dir = resolve_path(paths.get("export_dir", "results/export"))
        whitelist = config.get("export", {}).get("loupe", {}).get("whitelist_path", "")

        poly_json_path = proseg_dir / roi_name / "combined_proseg_results_qc.json" if roi_name else proseg_dir / "combined_proseg_results_qc.json"

        # 決定 h5ad 路徑
        if req.input_h5ad:
            h5ad_path = Path(req.input_h5ad)
        else:
            for candidate in ["clustered_final.h5ad", "umap_computed.h5ad", "qc_preprocessed.h5ad"]:
                p = output_dir_base / candidate
                if p.exists():
                    h5ad_path = p
                    break
            else:
                raise FileNotFoundError(f"找不到分析結果 h5ad，請先執行 Stage 4 或指定路徑。搜尋位置：{output_dir_base}")

        out_dir = Path(req.output_dir) if req.output_dir else export_dir / "loupe"

        exporter = LoupeExporter(
            poly_json_path=poly_json_path if poly_json_path.exists() else None,
            whitelist_path=resolve_path(whitelist) if whitelist else None,
        )
        await asyncio.get_event_loop().run_in_executor(
            None, exporter.export, h5ad_path, out_dir,
        )
        _loupe_status = {"status": "done", "progress": 1.0, "message": "Loupe 匯出完成"}
    except Exception as e:
        logger.error(f"Loupe 匯出失敗：{e}")
        _loupe_status = {"status": "error", "progress": 0.0, "message": str(e)}


@router.post("/xenium")
async def export_xenium(req: ExportRequest, background_tasks: BackgroundTasks):
    if _xenium_status["status"] == "running":
        return {"status": "error", "message": "任務執行中"}
    config = load_config()
    background_tasks.add_task(_run_xenium, config, req)
    return {"status": "ok", "message": "Xenium 匯出已啟動"}


@router.post("/loupe")
async def export_loupe(req: ExportRequest, background_tasks: BackgroundTasks):
    if _loupe_status["status"] == "running":
        return {"status": "error", "message": "任務執行中"}
    config = load_config()
    background_tasks.add_task(_run_loupe, config, req)
    return {"status": "ok", "message": "Loupe 匯出已啟動"}
