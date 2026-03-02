"""資料設定 API：掃描資料目錄、套用路徑配置"""
import logging
from typing import Optional

from fastapi import APIRouter
from pydantic import BaseModel

from backend.src.utils.config import load_config, save_config
from backend.src.utils.discovery import scan_data_root

router = APIRouter()
logger = logging.getLogger("pipeline.api.data")


class ScanRequest(BaseModel):
    data_root: str


class ApplyRequest(BaseModel):
    he_image: Optional[str] = None
    binned_002: Optional[str] = None
    binned_008: Optional[str] = None
    xenium_outs: Optional[str] = None


@router.post("/scan")
async def scan_directory(req: ScanRequest):
    """掃描指定目錄，回傳發現的檔案清單"""
    try:
        result = scan_data_root(req.data_root)
        return {"status": "ok", "data": result.to_dict()}
    except Exception as e:
        logger.error(f"掃描失敗：{e}")
        return {"status": "error", "message": str(e)}


@router.post("/apply")
async def apply_paths(req: ApplyRequest):
    """將發現的路徑寫入 pipeline.yaml"""
    try:
        config = load_config()
        paths = config.setdefault("paths", {})
        updates = req.model_dump(exclude_none=True)

        for key, value in updates.items():
            if value:  # 只更新非空值
                paths[key] = value

        # 同時保存 data_root（方便下次掃描）
        save_config(config)
        logger.info(f"已套用 {len(updates)} 項路徑設定")
        return {"status": "ok", "message": f"已更新 {len(updates)} 項路徑", "data": paths}
    except Exception as e:
        logger.error(f"套用失敗：{e}")
        return {"status": "error", "message": str(e)}


@router.get("/status")
async def get_data_status():
    """取得目前 paths 配置狀態（哪些已填、哪些為空）"""
    config = load_config()
    paths = config.get("paths", {})
    required_keys = ["he_image", "binned_002", "binned_008", "xenium_outs"]
    status = {}
    for key in required_keys:
        val = paths.get(key, "")
        status[key] = {
            "path": val,
            "configured": bool(val),
        }
    return {"status": "ok", "data": status}
