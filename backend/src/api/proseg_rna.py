"""Stage 2.5：Proseg RNA 重分配 API"""
import asyncio
import logging
from fastapi import APIRouter, BackgroundTasks, HTTPException
from pydantic import BaseModel
from typing import Optional, Any
import base64
import json
import numpy as np
import tifffile
import cv2
from pathlib import Path
from skimage.segmentation import find_boundaries

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
@router.get("/comparison/{roi_name}")
async def get_comparison(roi_name: str):
    """
    產生三層疊圖供前端比較：
    1. HE 底圖 (base64)
    2. Cellpose 輪廓 (base64 PNG, 透明背景)
    3. Proseg 輪廓 (base64 PNG, 透明背景)
    """
    config = load_config()
    from backend.src.utils.config import resolve_path
    from backend.src.api.export import _mask_to_geojson, _read_proseg_geojson
    from backend.src.utils.constants import VISIUM_UM_PX

    out_base = resolve_path(config["paths"]["output_dir"]) / "roi"
    roi_dir = out_base / roi_name
    he_path = roi_dir / "he_crop.tif"
    mask_path = roi_dir / "segmentation_masks.npy"
    proseg_json_path = roi_dir / "_proseg_work" / "proseg_results.json"

    if not he_path.exists():
        raise HTTPException(status_code=404, detail="找不到 H&E 影像")

    try:
        # 1. 讀取 HE 並轉為 base64
        he = tifffile.imread(str(he_path))
        if he.ndim == 3 and he.shape[-1] == 4:
            he = he[..., :3]
        H, W = he.shape[:2]
        _, buffer = cv2.imencode(".jpg", cv2.cvtColor(he, cv2.COLOR_RGB2BGR), [cv2.IMWRITE_JPEG_QUALITY, 85])
        he_b64 = base64.b64encode(buffer).decode()

        # 取得 ROI 參數
        rois = config.get("rois", [])
        roi_cfg = next((r for r in rois if r.get("name") == roi_name), {})
        pixel_size_um = float(roi_cfg.get("pixel_size_um", VISIUM_UM_PX))

        # 2. 生成 Cellpose 輪廓圖層
        cellpose_b64 = ""
        if mask_path.exists():
            overlay = np.zeros((H, W, 4), dtype=np.uint8) # BGRA
            seg_mask = np.load(str(mask_path))
            
            # 使用 find_boundaries 確保 1:1 像素精確度，避免 findContours 的 0.5px 位移
            # mode='thick' 對應 1px 像素邊界
            boundaries = find_boundaries(seg_mask, mode='thick')
            
            # Cyan: B=255, G=255, R=0, A=180
            overlay[boundaries] = [255, 255, 0, 180] 
            
            _, buffer = cv2.imencode(".png", overlay)
            cellpose_b64 = base64.b64encode(buffer).decode()

        # 3. 生成 Proseg 輪廓圖層
        proseg_b64 = ""
        if proseg_json_path.exists():
            overlay = np.zeros((H, W, 4), dtype=np.uint8) # BGRA
            try:
                proseg_geo = _read_proseg_geojson(proseg_json_path)
                for feat in proseg_geo.get("features", []):
                    geometry = feat.get("geometry", {})
                    g_type = geometry.get("type", "Polygon")
                    
                    coordinates = []
                    if g_type == "Polygon":
                        coordinates = geometry.get("coordinates", [])
                    elif g_type == "MultiPolygon":
                        for poly in geometry.get("coordinates", []):
                            coordinates.extend(poly)
                    
                    for ring in coordinates:
                        # µm -> px
                        pts = (np.array(ring) / pixel_size_um).astype(np.int32)
                        if pts.size >= 6:
                            # Red: B=0, G=0, R=255, A=200
                            cv2.polylines(overlay, [pts], True, (0, 0, 255, 200), 1)
                
                _, buffer = cv2.imencode(".png", overlay)
                proseg_b64 = base64.b64encode(buffer).decode()
            except Exception as e:
                logger.warning(f"解析 Proseg GeoJSON 失敗供比較視圖：{e}")

        return {
            "status": "ok",
            "data": {
                "he": he_b64,
                "cellpose": cellpose_b64,
                "proseg": proseg_b64,
                "width": W,
                "height": H
            }
        }
    except Exception as e:
        logger.error(f"產生比較視圖失敗：{e}", exc_info=True)
        return {"status": "error", "message": str(e)}
