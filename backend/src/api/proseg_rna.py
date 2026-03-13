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
                has_cyto   = (out_base / name / "segmentation_masks_cyto.npy").exists()
                has_proseg = (out_base / name / "proseg_cells.h5ad").exists()
                result.append({
                    "name":          name,
                    "has_adata":     has_adata,
                    "has_mask":      has_mask,
                    "has_cyto_mask": has_cyto,
                    "has_proseg":    has_proseg,
                })
            except (PermissionError, OSError) as e:
                logger.warning(f"跳過 {name}：{e}")
                continue
        return {"status": "ok", "data": result}
    except Exception as e:
        return {"status": "error", "message": str(e)}
@router.get("/comparison/{roi_name}")
async def get_proseg_comparison(
    roi_name: str,
    show_he: bool = True,
    show_cellpose: bool = True,
    show_proseg: bool = True
):
    """
    獲取 ROI 的對照影像。
    為了徹底解決 CSS 疊圖位移問題，改由後端進行合併渲染，回傳單一完美的疊圖。
    """
    config = load_config()
    from backend.src.utils.config import resolve_path
    from backend.src.api.export import _read_proseg_geojson
    from backend.src.utils.constants import VISIUM_UM_PX

    cyto_active = config.get("proseg", {}).get("stage25", {}).get("cyto_protection", False)

    out_base = resolve_path(config["paths"]["output_dir"]) / "roi"
    roi_dir = out_base / roi_name
    he_path = roi_dir / "he_crop.tif"
    mask_path = roi_dir / "segmentation_masks.npy"
    proseg_json_path = roi_dir / "_proseg_work" / "proseg_results.json"

    if not he_path.exists():
        raise HTTPException(status_code=404, detail="找不到 H&E 影像")

    try:
        # 1. 讀取基礎影像
        he = tifffile.imread(str(he_path))
        if he.ndim == 3 and he.shape[-1] == 4:
            he = he[..., :3]
        H, W = he.shape[:2]

        # 建立畫布層
        if show_he:
            # 必須複製，避免原地修改緩存中的影像
            canvas = he.copy().astype(np.uint8)
        else:
            canvas = np.zeros((H, W, 3), dtype=np.uint8)

        # 2. 疊加 Cellpose (亮青色)
        if show_cellpose and mask_path.exists():
            seg_mask = np.load(str(mask_path))
            boundaries = find_boundaries(seg_mask, mode='thick')
            # 青色: Cyan [R=0, G=255, B=255]
            canvas[boundaries] = [0, 255, 255]

        # 3. 疊加 Proseg (鮮紅色)
        if show_proseg and proseg_json_path.exists():
            rois = config.get("rois", [])
            roi_cfg = next((r for r in rois if r.get("name") == roi_name), {})
            pixel_size_um = float(roi_cfg.get("pixel_size_um", VISIUM_UM_PX))
            
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
                            # 紅色: Red [R=255, G=0, B=0]
                            cv2.polylines(canvas, [pts], True, (255, 0, 0), 1)
            except Exception as e:
                logger.warning(f"解析 Proseg GeoJSON 失敗：{e}")

        # 4. 轉為 base64 JPEG
        _, buffer = cv2.imencode(".jpg", cv2.cvtColor(canvas, cv2.COLOR_RGB2BGR), [cv2.IMWRITE_JPEG_QUALITY, 85])
        img_b64 = base64.b64encode(buffer).decode()

        return {
            "status": "ok",
            "data": {
                "combined_b64": img_b64,
                "width": W,
                "height": H
            }
        }
    except Exception as e:
        logger.error(f"產生比較視圖失敗：{e}", exc_info=True)
        return {"status": "error", "message": str(e)}
