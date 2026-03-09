"""Stage 4: Browser 格式匯出 API（Pipeline 3 版本，使用 Cellpose mask 轉多邊形）"""
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
_xenium_lock   = asyncio.Lock()
_loupe_lock    = asyncio.Lock()


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
# 核心：Cellpose mask → GeoJSON 多邊形
# ──────────────────────────────────────────────────────────────────────────────

def _mask_to_geojson(
    mask_path: Path,
    pixel_size_um: float,
    min_area_px: int = 20,
) -> dict:
    """
    將 segmentation_masks.npy 轉換為 GeoJSON FeatureCollection。

    座標：ROI 局部 µm（原點 = ROI 左上角），與 cellpose_cells.h5ad obsm['spatial'] 一致。
    """
    import numpy as np
    from skimage import measure

    seg_mask = np.load(str(mask_path))
    unique_ids = np.unique(seg_mask)
    unique_ids = unique_ids[unique_ids > 0]   # 去掉背景 0

    features = []
    for cid in unique_ids:
        cell_mask = (seg_mask == cid).astype(np.uint8)
        if cell_mask.sum() < min_area_px:
            continue

        # 找輪廓（skimage：row/col，padded 以捕捉邊緣細胞）
        padded = np.pad(cell_mask, 1, mode="constant")
        contours = measure.find_contours(padded, 0.5)
        if not contours:
            continue

        # 取最大輪廓
        contour = max(contours, key=len)
        # 還原 padding offset，再轉成 (x, y) µm
        # skimage contour = (row, col)，去掉 pad 偏移 1
        xy_um = np.column_stack([
            (contour[:, 1] - 1) * pixel_size_um,   # col → x
            (contour[:, 0] - 1) * pixel_size_um,   # row → y
        ])

        # 確保多邊形閉合
        if not np.allclose(xy_um[0], xy_um[-1]):
            xy_um = np.vstack([xy_um, xy_um[0]])

        # 簡化：點數 > 200 時間隔取樣（避免 Xenium Explorer 過慢）
        if len(xy_um) > 200:
            step = max(1, len(xy_um) // 100)
            xy_um = np.vstack([xy_um[::step], xy_um[-1:]])

        features.append({
            "type": "Feature",
            "geometry": {
                "type": "Polygon",
                "coordinates": [xy_um.tolist()],
            },
            "properties": {
                "full_id":   f"cell_{int(cid)}",
                "cell_id":   int(cid),
            },
        })

    logger.info(f"  生成 {len(features)} 個 Cellpose 多邊形")
    return {"type": "FeatureCollection", "features": features}


# ──────────────────────────────────────────────────────────────────────────────
# Xenium 匯出
# ──────────────────────────────────────────────────────────────────────────────

async def _run_xenium(config: dict, req: ExportRequest):
    global _xenium_status
    set_current_stage("export")
    _xenium_status = {"status": "running", "progress": 0.0, "message": "匯出至 Xenium Explorer..."}
    try:
        from backend.src.export.xenium_exporter import XeniumExporter
        import json

        paths = config.get("paths", {})
        rois  = config.get("rois", [{}])

        output_dir_base = resolve_path(paths.get("output_dir", "results/analysis"))
        export_dir      = resolve_path(paths.get("export_dir", "results/export"))

        # 初始化，避免 UnboundLocalError
        combined_poly_path: "Path | None" = None

        # 單 ROI vs 合併模式
        merged_h5ad = output_dir_base / "umap_computed.h5ad"
        is_merged_mode = merged_h5ad.exists() and len(rois) > 1

        if is_merged_mode:
            logger.info(f"合併模式（{len(rois)} 個 ROI）")
            all_features: list = []
            for roi in rois:
                rn = roi.get("name", "")
                if not rn:
                    continue
                roi_out_dir   = output_dir_base / "roi" / rn
                mask_path     = roi_out_dir / "segmentation_masks.npy"
                pixel_size_um = float(roi.get("pixel_size_um", 0.2737))

                if not mask_path.exists():
                    logger.warning(f"  [{rn}] 找不到遮罩，跳過")
                    continue

                roi_geo = _mask_to_geojson(mask_path, pixel_size_um)

                # 加入全域座標偏移
                roi_x_um = roi.get("x", 0) * roi.get("pixel_size_um", 0.2737)
                roi_y_um = roi.get("y", 0) * roi.get("pixel_size_um", 0.2737)
                for feat in roi_geo.get("features", []):
                    orig_id = feat["properties"].get("full_id", "")
                    feat["properties"]["full_id"] = f"{rn}__{orig_id}"
                    _shift_geojson_coords(feat, roi_x_um, roi_y_um)
                    all_features.append(feat)
                logger.info(f"  [{rn}] {len(roi_geo['features'])} 個多邊形")

            combined_poly_path = output_dir_base / "combined_cellpose_polygons.json"
            with open(combined_poly_path, "w", encoding="utf-8") as f:
                json.dump({"type": "FeatureCollection", "features": all_features}, f)

            h5ad_path = merged_h5ad if not req.input_h5ad else Path(req.input_h5ad)
            roi_pixel_size_um = rois[0].get("pixel_size_um", 0.2737) if rois else 0.2737

            he_image_path = resolve_path(paths.get("he_image", "")) if paths.get("he_image") else None
            he_crop_bounds = None
            if he_image_path and rois:
                _x0 = min(r.get("x", 0) for r in rois)
                _y0 = min(r.get("y", 0) for r in rois)
                _x1 = max(r.get("x", 0) + r.get("width_px", 0) for r in rois)
                _y1 = max(r.get("y", 0) + r.get("height_px", 0) for r in rois)
                he_crop_bounds = (_x0, _y0, _x1, _y1)

        else:
            roi_name    = rois[0].get("name", "") if rois else ""
            roi_out_dir = output_dir_base / "roi" / roi_name if roi_name else output_dir_base
            mask_path   = roi_out_dir / "segmentation_masks.npy"
            pixel_size_um = float(rois[0].get("pixel_size_um", 0.2737)) if rois else 0.2737

            if not mask_path.exists():
                raise FileNotFoundError(f"找不到 {roi_name} 的 segmentation_masks.npy，請先完成 Stage 1")

            logger.info(f"單 ROI 模式（{roi_name}），從 Cellpose 遮罩生成多邊形...")
            roi_geo = _mask_to_geojson(mask_path, pixel_size_um)
            combined_poly_path = roi_out_dir / "cellpose_polygons.json"
            with open(combined_poly_path, "w", encoding="utf-8") as f:
                json.dump(roi_geo, f)

            h5ad_path = None
            if req.input_h5ad:
                h5ad_path = Path(req.input_h5ad)
            else:
                for candidate in ["umap_computed.h5ad", "qc_preprocessed.h5ad", "cellpose_cells.h5ad"]:
                    for search_dir in [roi_out_dir, output_dir_base]:
                        p = search_dir / candidate
                        if p.exists():
                            h5ad_path = p
                            break
                    if h5ad_path:
                        break
            if h5ad_path is None:
                raise FileNotFoundError(
                    f"找不到分析結果 h5ad，請先執行 Stage 3 分析或指定路徑。"
                    f"搜尋位置：{output_dir_base}、{roi_out_dir}"
                )

            roi_pixel_size_um = pixel_size_um
            he_image_path     = None
            he_crop_bounds    = None

        out_dir = Path(req.output_dir) if req.output_dir else export_dir / "xenium"

        exporter = XeniumExporter(
            zarr_path=None,
            poly_json_path=combined_poly_path if (combined_poly_path and combined_poly_path.exists()) else None,
            transcripts_csv_path=None,
            pixel_size_um=roi_pixel_size_um,
            he_image_path=he_image_path,
            he_crop_bounds=he_crop_bounds,
        )
        await asyncio.get_running_loop().run_in_executor(
            None, exporter.export, h5ad_path, out_dir,
        )
        _xenium_status = {"status": "done", "progress": 1.0, "message": "Xenium 匯出完成"}
    except Exception as e:
        logger.error(f"Xenium 匯出失敗：{e}")
        _xenium_status = {"status": "error", "progress": 0.0, "message": str(e)}


# ──────────────────────────────────────────────────────────────────────────────
# Loupe 匯出
# ──────────────────────────────────────────────────────────────────────────────

async def _run_loupe(config: dict, req: ExportRequest):
    global _loupe_status
    set_current_stage("export")
    _loupe_status = {"status": "running", "progress": 0.0, "message": "匯出至 Loupe Browser..."}
    try:
        from backend.src.export.loupe_exporter import LoupeExporter

        paths    = config.get("paths", {})
        rois     = config.get("rois", [{}])
        roi_name = rois[0].get("name", "") if rois else ""

        output_dir_base = resolve_path(paths.get("output_dir", "results/analysis"))
        export_dir      = resolve_path(paths.get("export_dir", "results/export"))
        whitelist       = config.get("export", {}).get("loupe", {}).get("whitelist_path", "")

        # 多邊形 JSON（Cellpose mask 轉出）
        roi_out_dir    = output_dir_base / "roi" / roi_name if roi_name else output_dir_base
        poly_json_path = roi_out_dir / "cellpose_polygons.json"
        if not poly_json_path.exists():
            # 嘗試合併版路徑
            alt = output_dir_base / "combined_cellpose_polygons.json"
            poly_json_path = alt if alt.exists() else None

        # h5ad
        if req.input_h5ad:
            h5ad_path = Path(req.input_h5ad)
        else:
            for candidate in ["clustered_final.h5ad", "umap_computed.h5ad", "qc_preprocessed.h5ad"]:
                p = output_dir_base / candidate
                if p.exists():
                    h5ad_path = p
                    break
            else:
                raise FileNotFoundError(f"找不到分析結果 h5ad，請先執行 Stage 3 分析。搜尋位置：{output_dir_base}")

        out_dir = Path(req.output_dir) if req.output_dir else export_dir / "loupe"

        exporter = LoupeExporter(
            poly_json_path=poly_json_path if poly_json_path and poly_json_path.exists() else None,
            whitelist_path=resolve_path(whitelist) if whitelist else None,
        )
        await asyncio.get_running_loop().run_in_executor(
            None, exporter.export, h5ad_path, out_dir,
        )
        _loupe_status = {"status": "done", "progress": 1.0, "message": "Loupe 匯出完成"}
    except Exception as e:
        logger.error(f"Loupe 匯出失敗：{e}")
        _loupe_status = {"status": "error", "progress": 0.0, "message": str(e)}


@router.post("/xenium")
async def export_xenium(req: ExportRequest, background_tasks: BackgroundTasks):
    async with _xenium_lock:
        if _xenium_status["status"] == "running":
            return {"status": "error", "message": "任務執行中"}
    config = load_config()
    background_tasks.add_task(_run_xenium, config, req)
    return {"status": "ok", "message": "Xenium 匯出已啟動"}


@router.post("/loupe")
async def export_loupe(req: ExportRequest, background_tasks: BackgroundTasks):
    async with _loupe_lock:
        if _loupe_status["status"] == "running":
            return {"status": "error", "message": "任務執行中"}
    config = load_config()
    background_tasks.add_task(_run_loupe, config, req)
    return {"status": "ok", "message": "Loupe 匯出已啟動"}


# ──────────────────────────────────────────────────────────────────────────────
# 工具函式
# ──────────────────────────────────────────────────────────────────────────────

def _shift_geojson_coords(feat: dict, dx: float, dy: float):
    """In-place 平移 GeoJSON feature 的座標。"""
    def _shift(coords):
        if not coords:
            return coords
        if isinstance(coords[0], (int, float)):
            return [coords[0] + dx, coords[1] + dy] + list(coords[2:])
        return [_shift(c) for c in coords]

    geom = feat.get("geometry", {})
    if geom and geom.get("coordinates") is not None:
        geom["coordinates"] = _shift(geom["coordinates"])
