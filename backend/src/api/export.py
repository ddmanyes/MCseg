"""Stage 4: Browser 格式匯出 API（Pipeline 3 版本，使用 Cellpose mask 轉多邊形）"""
import asyncio
import logging
from pathlib import Path
from fastapi import APIRouter, BackgroundTasks
from pydantic import BaseModel

from backend.src.utils.config import load_config, resolve_path
from backend.src.utils.constants import VISIUM_UM_PX
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
    mask_source: str = "auto"  # "auto", "cellpose", "proseg"


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
                "full_id":   str(int(cid) - 1),   # Cellpose cid 從 1 起，h5ad obs_names 從 '0' 起
                "cell_id":   int(cid),
            },
        })

    logger.info(f"  生成 {len(features)} 個 Cellpose 多邊形")
    return {"type": "FeatureCollection", "features": features}


# ──────────────────────────────────────────────────────────────────────────────
# Xenium 匯出
# ──────────────────────────────────────────────────────────────────────────────

def _read_proseg_geojson(path: Path) -> dict:
    """讀取 proseg_results.json（Proseg 輸出為 gzip 壓縮，副檔名仍是 .json）"""
    import gzip, json
    with open(str(path), "rb") as fh:
        magic = fh.read(2)
    if magic == b"\x1f\x8b":
        with gzip.open(str(path), "rt", encoding="utf-8") as fh:
            return json.load(fh)
    with open(str(path), "r", encoding="utf-8") as fh:
        return json.load(fh)


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
        data_root_dir   = resolve_path(paths.get("data_root", "."))

        # Find h5ad_path first
        h5ad_path = None
        if req.input_h5ad:
            p = Path(req.input_h5ad)
            if p.is_absolute():
                h5ad_path = p
            else:
                # 嘗試從 output_dir 與 data_root 兩處查找
                candidate_a = output_dir_base / p
                candidate_b = data_root_dir / p
                if candidate_a.exists():
                    h5ad_path = candidate_a
                elif candidate_b.exists():
                    h5ad_path = candidate_b
                else:
                    raise FileNotFoundError(f"找不到指定的 h5ad：{req.input_h5ad}\n搜尋位置：{candidate_a} 或 {candidate_b}")
        else:
            for candidate in ["umap_computed.h5ad", "qc_preprocessed.h5ad", "cellpose_cells.h5ad"]:
                p = output_dir_base / candidate
                if p.exists():
                    h5ad_path = p
                    break
        if h5ad_path is None:
            raise FileNotFoundError(f"找不到分析結果 h5ad，請先執行 Stage 3 分析。\n搜尋位置：{output_dir_base}")

        # Check if the h5ad is merged mode by inspecting obs_names
        import scanpy as sc
        adata_head = sc.read_h5ad(str(h5ad_path), backed="r")
        first_obs = adata_head.obs_names[0] if len(adata_head) > 0 else ""
        is_merged_mode = "__" in first_obs
        active_roi = adata_head.uns.get("active_roi", None) if "active_roi" in adata_head.uns else None
        del adata_head # Release backed file handle

        combined_poly_path: "Path | None" = None

        if is_merged_mode:
            logger.info(f"合併模式（{len(rois)} 個 ROI），根據 h5ad 的 obs_names 判定")
            all_features: list = []
            for roi in rois:
                rn = roi.get("name", "")
                if not rn:
                    continue
                roi_out_dir   = output_dir_base / "roi" / rn
                # Proseg 資料從 data_root 或 output_dir 查找
                proseg_json = data_root_dir / "roi" / rn / "_proseg_work" / "proseg_results.json"
                if not proseg_json.exists():
                    proseg_json = roi_out_dir / "_proseg_work" / "proseg_results.json"
                mask_path   = roi_out_dir / "segmentation_masks.npy"
                pixel_size_um = float(roi.get("pixel_size_um", VISIUM_UM_PX))

                if req.mask_source in ["auto", "proseg"] and proseg_json.exists():
                    logger.info(f"  [{rn}] 發現 Proseg 結果，使用 Proseg 擴展邊界")
                    roi_geo = _read_proseg_geojson(proseg_json)
                    for feat in roi_geo.get("features", []):
                        if "full_id" not in feat["properties"]:
                            cid = feat["properties"].get("cell") or feat["properties"].get("cell_id")
                            if cid is not None:
                                # proseg_cells.h5ad 的 obs_names 是純數字字串 '0','1','2'...
                                feat["properties"]["full_id"] = str(int(cid))
                elif req.mask_source in ["auto", "cellpose"] and mask_path.exists():
                    logger.info(f"  [{rn}] 使用 Cellpose 遮罩生成多邊形")
                    roi_geo = _mask_to_geojson(mask_path, pixel_size_um)
                else:
                    logger.warning(f"  [{rn}] 找不到選擇的 {req.mask_source} 遮罩來源，跳過")
                    continue

                # 加入全域座標偏移
                roi_x_um = roi.get("x", 0) * roi.get("pixel_size_um", VISIUM_UM_PX)
                roi_y_um = roi.get("y", 0) * roi.get("pixel_size_um", VISIUM_UM_PX)
                for feat in roi_geo.get("features", []):
                    orig_id = feat["properties"].get("full_id", "")
                    feat["properties"]["full_id"] = f"{rn}__{orig_id}"
                    _shift_geojson_coords(feat, roi_x_um, roi_y_um)
                    all_features.append(feat)
                logger.info(f"  [{rn}] {len(roi_geo['features'])} 個多邊形")

            combined_poly_path = output_dir_base / "combined_cellpose_polygons.json"
            with open(combined_poly_path, "w", encoding="utf-8") as f:
                json.dump({"type": "FeatureCollection", "features": all_features}, f)

            roi_pixel_size_um = rois[0].get("pixel_size_um", VISIUM_UM_PX) if rois else VISIUM_UM_PX

            he_image_path = resolve_path(paths.get("he_image", "")) if paths.get("he_image") else None
            he_crop_bounds = None
            if he_image_path and rois:
                _x0 = min(r.get("x", 0) for r in rois)
                _y0 = min(r.get("y", 0) for r in rois)
                _x1 = max(r.get("x", 0) + r.get("width_px", 0) for r in rois)
                _y1 = max(r.get("y", 0) + r.get("height_px", 0) for r in rois)
                he_crop_bounds = (_x0, _y0, _x1, _y1)

        else:
            roi_name    = active_roi or (rois[-1].get("name", "") if rois else "")
            roi_out_dir = output_dir_base / "roi" / roi_name if roi_name else output_dir_base
            mask_path   = roi_out_dir / "segmentation_masks.npy"
            # Proseg 資料從 data_root 或 output_dir 查找
            proseg_json = data_root_dir / "roi" / roi_name / "_proseg_work" / "proseg_results.json"
            if not proseg_json.exists():
                proseg_json = roi_out_dir / "_proseg_work" / "proseg_results.json"
            
            # Find the correct pixel_size_um for this ROI
            roi_cfg = next((r for r in rois if r.get("name") == roi_name), {})
            pixel_size_um = float(roi_cfg.get("pixel_size_um", VISIUM_UM_PX))

            if req.mask_source in ["auto", "proseg"] and proseg_json.exists():
                logger.info(f"單 ROI 模式（{roi_name}），從 Proseg 結果生成多邊形...")
                roi_geo = _read_proseg_geojson(proseg_json)
                for feat in roi_geo.get("features", []):
                    if "full_id" not in feat["properties"]:
                        cid = feat["properties"].get("cell") or feat["properties"].get("cell_id")
                        if cid is not None:
                            feat["properties"]["full_id"] = str(int(cid))
            elif req.mask_source in ["auto", "cellpose"] and mask_path.exists():
                logger.info(f"單 ROI 模式（{roi_name}），從 Cellpose 遮罩生成多邊形...")
                roi_geo = _mask_to_geojson(mask_path, pixel_size_um)
            else:
                raise FileNotFoundError(f"找不到 {roi_name} 的 {req.mask_source} 遮罩來源或檔案不存在")
            combined_poly_path = roi_out_dir / "cellpose_polygons.json"
            with open(combined_poly_path, "w", encoding="utf-8") as f:
                json.dump(roi_geo, f)

            roi_pixel_size_um = pixel_size_um
            he_image_path     = roi_out_dir / "he_crop.tif"
            if not he_image_path.exists():
                logger.warning(f"單 ROI 模式找不到 he_crop.tif，將不用底圖匯出")
                he_image_path = None
            he_crop_bounds    = None

        if req.output_dir:
            out_dir = Path(req.output_dir)
        else:
            out_dir = export_dir / "xenium" if is_merged_mode else roi_out_dir / "export_xenium"

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
        data_root_dir   = resolve_path(paths.get("data_root", "."))
        whitelist       = config.get("export", {}).get("loupe", {}).get("whitelist_path", "")

        # Find h5ad_path first
        h5ad_path = None
        if req.input_h5ad:
            p = Path(req.input_h5ad)
            if p.is_absolute():
                h5ad_path = p
            else:
                candidate_a = output_dir_base / p
                candidate_b = data_root_dir / p
                if candidate_a.exists():
                    h5ad_path = candidate_a
                elif candidate_b.exists():
                    h5ad_path = candidate_b
                else:
                    raise FileNotFoundError(f"找不到指定的 h5ad：{req.input_h5ad}")
        else:
            for candidate in ["clustered_final.h5ad", "umap_computed.h5ad", "qc_preprocessed.h5ad"]:
                p = output_dir_base / candidate
                if p.exists():
                    h5ad_path = p
                    break
            if h5ad_path is None:
                raise FileNotFoundError(f"找不到分析結果 h5ad，請先執行 Stage 3 分析。搜尋位置：{output_dir_base}")

        # Check if the h5ad is merged mode by inspecting obs_names
        import scanpy as sc
        adata_head = sc.read_h5ad(str(h5ad_path), backed="r")
        first_obs = adata_head.obs_names[0] if len(adata_head) > 0 else ""
        is_merged_mode = "__" in first_obs
        active_roi = adata_head.uns.get("active_roi", None) if "active_roi" in adata_head.uns else None
        del adata_head

        poly_json_path: "Path | None" = None
        import json

        if is_merged_mode:
            logger.info("Loupe 匯出：合併模式，產生 combined_cellpose_polygons.json")
            all_features: list = []
            for roi in rois:
                rn = roi.get("name", "")
                if not rn:
                    continue
                roi_out_dir   = output_dir_base / "roi" / rn
                # Proseg 資料從 data_root 或 output_dir 查找
                proseg_json = data_root_dir / "roi" / rn / "_proseg_work" / "proseg_results.json"
                if not proseg_json.exists():
                    proseg_json = roi_out_dir / "_proseg_work" / "proseg_results.json"
                mask_path     = roi_out_dir / "segmentation_masks.npy"
                pixel_size_um = float(roi.get("pixel_size_um", VISIUM_UM_PX))

                if req.mask_source in ["auto", "proseg"] and proseg_json.exists():
                    roi_geo = _read_proseg_geojson(proseg_json)
                    for feat in roi_geo.get("features", []):
                        if "full_id" not in feat["properties"]:
                            cid = feat["properties"].get("cell") or feat["properties"].get("cell_id")
                            if cid is not None:
                                feat["properties"]["full_id"] = str(int(cid))
                elif req.mask_source in ["auto", "cellpose"] and mask_path.exists():
                    roi_geo = _mask_to_geojson(mask_path, pixel_size_um)
                else:
                    continue

                roi_x_um = roi.get("x", 0) * roi.get("pixel_size_um", VISIUM_UM_PX)
                roi_y_um = roi.get("y", 0) * roi.get("pixel_size_um", VISIUM_UM_PX)
                for feat in roi_geo.get("features", []):
                    orig_id = feat["properties"].get("full_id", "")
                    feat["properties"]["full_id"] = f"{rn}__{orig_id}"
                    _shift_geojson_coords(feat, roi_x_um, roi_y_um)
                    all_features.append(feat)

            poly_json_path = output_dir_base / "combined_cellpose_polygons.json"
            with open(poly_json_path, "w", encoding="utf-8") as f:
                json.dump({"type": "FeatureCollection", "features": all_features}, f)
        else:
            roi_name    = active_roi or (rois[-1].get("name", "") if rois else "")
            roi_out_dir = output_dir_base / "roi" / roi_name if roi_name else output_dir_base
            mask_path   = roi_out_dir / "segmentation_masks.npy"
            # Proseg 資料從 data_root 或 output_dir 查找
            proseg_json = data_root_dir / "roi" / roi_name / "_proseg_work" / "proseg_results.json"
            if not proseg_json.exists():
                proseg_json = roi_out_dir / "_proseg_work" / "proseg_results.json"
            
            roi_cfg = next((r for r in rois if r.get("name") == roi_name), {})
            pixel_size_um = float(roi_cfg.get("pixel_size_um", VISIUM_UM_PX))

            if req.mask_source in ["auto", "proseg"] and proseg_json.exists():
                logger.info(f"Loupe 匯出：單 ROI 模式（{roi_name}），發現 Proseg 結果")
                roi_geo = _read_proseg_geojson(proseg_json)
                for feat in roi_geo.get("features", []):
                    if "full_id" not in feat["properties"]:
                        cid = feat["properties"].get("cell") or feat["properties"].get("cell_id")
                        if cid is not None:
                            feat["properties"]["full_id"] = str(int(cid))
                poly_json_path = roi_out_dir / "cellpose_polygons.json"
                with open(poly_json_path, "w", encoding="utf-8") as f:
                    json.dump(roi_geo, f)
            elif req.mask_source in ["auto", "cellpose"] and mask_path.exists():
                logger.info(f"Loupe 匯出：單 ROI 模式（{roi_name}），從 Cellpose 遮罩生成")
                roi_geo = _mask_to_geojson(mask_path, pixel_size_um)
                poly_json_path = roi_out_dir / "cellpose_polygons.json"
                with open(poly_json_path, "w", encoding="utf-8") as f:
                    json.dump(roi_geo, f)
            else:
                logger.warning(f"找不到 {roi_name} 的 {req.mask_source} 遮罩來源或 proseg，將不匯出多邊形層")

        if req.output_dir:
            out_dir = Path(req.output_dir)
        else:
            out_dir = export_dir / "loupe" if is_merged_mode else roi_out_dir / "export_loupe"

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
