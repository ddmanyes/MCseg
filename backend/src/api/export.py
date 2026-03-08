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

        zarr_dir = resolve_path(paths.get("zarr_dir", "results/zarr"))
        output_dir_base = resolve_path(paths.get("output_dir", "results/analysis"))
        export_dir = resolve_path(paths.get("export_dir", "results/export"))

        # ── 偵測是否為合併模式 ──────────────────────────────────────
        # Stage 4 合併模式的輸出（umap_computed.h5ad）位於 output_dir_base 而非 roi 子目錄
        merged_h5ad = output_dir_base / "umap_computed.h5ad"
        is_merged_mode = merged_h5ad.exists() and len(rois) > 1

        if is_merged_mode:
            logger.info(f"偵測到合併模式（{len(rois)} 個 ROI），合併所有 GeoJSON...")
            combined_poly_path = output_dir_base / "combined_all_rois.json"

            if not combined_poly_path.exists():
                all_features: list = []
                for roi in rois:
                    rn = roi.get("name", "")
                    if not rn:
                        continue
                    roi_out_dir = output_dir_base / "roi" / rn
                    roi_poly_path = roi_out_dir / "combined_proseg_results_qc.json"
                    roi_zarr = zarr_dir / rn / "proseg_integrated.zarr"

                    # 若 ROI 的合併 GeoJSON 不存在，從 tile 重建
                    if not roi_poly_path.exists():
                        tile_dir = roi_out_dir / "proseg_tiles"
                        if tile_dir.exists() and roi_zarr.exists():
                            logger.info(f"  [{rn}] 從 tile 重建 GeoJSON...")
                            roi_combined = generate_combined_geojson(tile_dir, roi_zarr, config)
                            roi_poly_path.parent.mkdir(parents=True, exist_ok=True)
                            with open(roi_poly_path, "w", encoding="utf-8") as f:
                                json.dump(roi_combined, f)
                        else:
                            logger.warning(f"  [{rn}] 找不到 proseg_tiles 或 zarr，跳過")
                            continue

                    # 載入 GeoJSON，為 full_id 加上 ROI 前綴（對應 merge_all_rois 的 obs_names 格式）
                    with open(roi_poly_path, "r", encoding="utf-8") as f:
                        roi_geo = json.load(f)
                    for feat in roi_geo.get("features", []):
                        orig_id = feat["properties"].get("full_id", "")
                        feat["properties"]["full_id"] = f"{rn}__{orig_id}"
                        all_features.append(feat)
                    logger.info(f"  [{rn}] 加入 {len(roi_geo.get('features', []))} 個多邊形")

                combined = {"type": "FeatureCollection", "features": all_features}
                with open(combined_poly_path, "w", encoding="utf-8") as f:
                    json.dump(combined, f)
                logger.info(f"合併 GeoJSON 已儲存：{combined_poly_path}（{len(all_features)} 個多邊形）")

            h5ad_path = merged_h5ad if not req.input_h5ad else Path(req.input_h5ad)
            zarr_path = None  # 合併模式無單一 zarr 影像
            roi_pixel_size_um = rois[0].get("pixel_size_um", 0.2737) if rois else 0.2737

        else:
            # ── 單 ROI 模式（原有邏輯）──────────────────────────────
            roi_name = rois[0].get("name", "") if rois else ""
            zarr_path = zarr_dir / roi_name / "proseg_integrated.zarr" if roi_name else zarr_dir / "proseg_integrated.zarr"
            roi_out_dir = output_dir_base / "roi" / roi_name if roi_name else output_dir_base
            tile_proseg_dir = roi_out_dir / "proseg_tiles"

            combined_poly_path = roi_out_dir / "combined_proseg_results_qc.json"
            if not combined_poly_path.exists() and tile_proseg_dir.exists() and zarr_path.exists():
                logger.info("未找到合併 GeoJSON，從 tile 自動生成...")
                combined = generate_combined_geojson(tile_proseg_dir, zarr_path, config)
                combined_poly_path.parent.mkdir(parents=True, exist_ok=True)
                with open(combined_poly_path, "w", encoding="utf-8") as f:
                    json.dump(combined, f)
                logger.info(f"合併 GeoJSON 已儲存：{combined_poly_path}（{len(combined['features'])} 個多邊形）")

            if req.input_h5ad:
                h5ad_path = Path(req.input_h5ad)
            else:
                h5ad_path = None
                for candidate in ["umap_computed.h5ad", "qc_preprocessed.h5ad", "proseg_cells.h5ad"]:
                    for search_dir in [roi_out_dir, output_dir_base]:
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

            roi_pixel_size_um = rois[0].get("pixel_size_um", 0.2737) if rois else 0.2737

        out_dir = Path(req.output_dir) if req.output_dir else export_dir / "xenium"

        exporter = XeniumExporter(
            zarr_path=zarr_path if zarr_path and zarr_path.exists() else None,
            poly_json_path=combined_poly_path if combined_poly_path.exists() else None,
            transcripts_csv_path=None,
            pixel_size_um=roi_pixel_size_um,
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
