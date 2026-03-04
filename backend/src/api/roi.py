"""Stage 0: ROI 定義與裁切 API"""
import asyncio
import json
import logging
from pathlib import Path
from typing import Optional
from fastapi import APIRouter, BackgroundTasks
from fastapi.responses import Response as FastAPIResponse
from pydantic import BaseModel

from backend.src.utils.config import load_config, save_config
from backend.src.utils.logging import set_current_stage

router = APIRouter()
logger = logging.getLogger("pipeline.api.roi")

# 任務狀態追蹤
_task_status = {"status": "idle", "progress": 0.0, "message": ""}

# DZI tile server singleton（避免重複開啟 BTF）
_tile_server = None


def _get_tile_server():
    global _tile_server
    if _tile_server is None:
        from backend.src.roi.tile_server import DZITileServer
        config   = load_config()
        he_path  = config.get('paths', {}).get('he_image', '')
        b002     = config.get('paths', {}).get('binned_002', '')
        hires    = str(Path(b002) / 'spatial' / 'tissue_hires_image.png') if b002 else None
        scalef   = 0.1
        if b002:
            sf_file = Path(b002) / 'spatial' / 'scalefactors_json.json'
            if sf_file.exists():
                with open(sf_file) as f:
                    scalef = json.load(f).get('tissue_hires_scalef', 0.1)
        _tile_server = DZITileServer(he_path, hires, scalef)
    return _tile_server


class RoiConfig(BaseModel):
    name: str
    tissue: str
    # 格式 A：fullres pixel
    x: Optional[int] = None
    y: Optional[int] = None
    width_px: Optional[int] = None
    height_px: Optional[int] = None
    pixel_size_um: float = 0.2737
    # 格式 B：µm
    x_um: Optional[float] = None
    y_um: Optional[float] = None
    width_um: Optional[float] = None
    height_um: Optional[float] = None


@router.get("/status")
async def get_status():
    return _task_status


@router.get("/list")
async def list_rois():
    config = load_config()
    return {"status": "ok", "data": config.get("rois", [])}


@router.post("/add")
async def add_roi(roi: RoiConfig):
    config = load_config()
    rois = config.get("rois", [])
    rois = [r for r in rois if r.get("name") != roi.name]  # 去重
    rois.append(roi.model_dump(exclude_none=True))
    config["rois"] = rois
    save_config(config)
    return {"status": "ok", "message": f"ROI '{roi.name}' 已新增"}


@router.delete("/{roi_name}")
async def delete_roi(roi_name: str):
    config = load_config()
    rois = [r for r in config.get("rois", []) if r.get("name") != roi_name]
    config["rois"] = rois
    save_config(config)
    return {"status": "ok", "message": f"ROI '{roi_name}' 已刪除"}


@router.get("/overview")
async def get_overview():
    """取得 ROI 預覽基礎縮圖與 metadata"""
    set_current_stage("roi")
    config = load_config()
    try:
        from backend.src.roi.extractor import get_overview as ext_get_overview
        data = ext_get_overview(config)
        return {"status": "ok", "data": data}
    except Exception as e:
        logger.error(f"取得 overview 失敗：{e}")
        return {"status": "error", "message": str(e)}


async def _run_extract(config: dict):
    global _task_status
    set_current_stage("roi")
    _task_status = {"status": "running", "progress": 0.0, "message": "開始裁切..."}
    try:
        from backend.src.roi.extractor import RoiExtractor
        extractor = RoiExtractor(config)
        await asyncio.get_event_loop().run_in_executor(None, extractor.run_all)
        _task_status = {"status": "done", "progress": 1.0, "message": "裁切完成"}
    except Exception as e:
        logger.error(f"ROI 裁切失敗：{e}")
        _task_status = {"status": "error", "progress": 0.0, "message": str(e)}


@router.post("/extract")
async def run_extract(background_tasks: BackgroundTasks):
    if _task_status["status"] == "running":
        return {"status": "error", "message": "任務執行中"}
    config = load_config()
    background_tasks.add_task(_run_extract, config)
    return {"status": "ok", "message": "ROI 裁切已啟動"}


# ── OpenSeadragon DZI tile server ─────────────────────────────────────────────

@router.get("/dzi")
async def roi_dzi():
    """DZI XML descriptor for OpenSeadragon viewer."""
    try:
        ts = _get_tile_server()
        return FastAPIResponse(content=ts.get_dzi(), media_type='application/xml')
    except Exception as e:
        logger.error(f"DZI descriptor 失敗：{e}")
        return FastAPIResponse(content=f"<error>{e}</error>", status_code=500,
                               media_type='application/xml')


@router.get("/tiles/{level}/{tile_name}")
async def roi_tile(level: int, tile_name: str):
    """Serve a single DZI tile as JPEG (tile_name format: '{tx}_{ty}.jpg')."""
    try:
        ts  = _get_tile_server()
        xy  = tile_name.replace('.jpg', '').replace('.jpeg', '').split('_')
        tx, ty = int(xy[0]), int(xy[1])
        jpg = ts.get_tile(level, tx, ty)
        return FastAPIResponse(
            content=jpg,
            media_type='image/jpeg',
            headers={'Cache-Control': 'public, max-age=3600'},
        )
    except Exception as e:
        logger.error(f"Tile {level}/{tile_name} 失敗：{e}")
        return FastAPIResponse(content=b'', status_code=500, media_type='image/jpeg')


@router.get("/dzi_files/{level}/{tile_name}")
async def roi_dzi_files(level: int, tile_name: str):
    """OSD derives tile URL as '{dzi_url}_files/{level}/{tx}_{ty}.jpeg' — alias to roi_tile."""
    return await roi_tile(level, tile_name)
