"""資料設定 API：掃描資料目錄、套用路徑配置、目錄瀏覽"""
import logging
import os
from pathlib import Path
from typing import Optional

from fastapi import APIRouter, Query
from pydantic import BaseModel

from backend.src.utils.config import load_config, resolve_path, save_config
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


@router.get("/disk-status")
async def get_disk_status():
    """
    掃描磁碟，回傳各 Stage 實際完成狀態（重啟後前端可據此恢復 UI 狀態）。
    """
    config = load_config()
    paths = config.get("paths", {})
    output_dir = resolve_path(paths.get("output_dir", "results/analysis"))
    zarr_dir   = resolve_path(paths.get("zarr_dir",   "results/zarr"))

    roi_base = output_dir / "roi"

    # Stage 0: ROI — 有任何 he_crop.tif
    roi_dirs = [d for d in roi_base.iterdir() if d.is_dir() and not d.name.startswith(".")] \
               if roi_base.exists() else []
    roi_done = any((d / "he_crop.tif").exists() for d in roi_dirs)
    roi_names = [d.name for d in roi_dirs if (d / "he_crop.tif").exists()]

    # Stage 1: Segmentation — 有任何 segmentation_masks.npy
    seg_done = any((d / "segmentation_masks.npy").exists() for d in roi_dirs)

    # Stage 2: Zarr — zarr 目錄下有 *.zarr
    zarr_files = list(zarr_dir.glob("*.zarr")) if zarr_dir.exists() else []
    zarr_done = len(zarr_files) > 0

    # Stage 3: Proseg — 有任何 proseg_cells.h5ad 或 cells.geojson
    proseg_done = any(
        (d / "proseg_cells.h5ad").exists() or (d / "cells.geojson").exists()
        for d in roi_dirs
    )

    # Stage 4: Analysis — 有任何 clustering.h5ad
    analysis_done = any((d / "clustering.h5ad").exists() for d in roi_dirs)

    return {
        "status": "ok",
        "data": {
            "roi":         {"done": roi_done,      "roi_names": roi_names},
            "segmentation":{"done": seg_done},
            "zarr":        {"done": zarr_done,     "files": [str(f) for f in zarr_files]},
            "proseg":      {"done": proseg_done},
            "analysis":    {"done": analysis_done},
        },
    }


@router.get("/browse")
async def browse_directory(path: str = Query("~", description="要瀏覽的目錄路徑")):
    """
    瀏覽本機目錄結構。回傳指定路徑下的子目錄與大型檔案。
    前端可用此 API 實現點擊式的資料夾選擇器。
    """
    try:
        target = Path(os.path.expanduser(path)).resolve()
        if not target.exists():
            return {"status": "error", "message": f"路徑不存在：{target}"}
        if not target.is_dir():
            return {"status": "error", "message": f"不是目錄：{target}"}

        items = []
        try:
            entries = sorted(target.iterdir(), key=lambda e: (not e.is_dir(), e.name.lower()))
        except PermissionError:
            return {"status": "error", "message": f"無權限存取：{target}"}

        for entry in entries:
            try:
                # 跳過隱藏檔案和系統目錄
                if entry.name.startswith(".") or entry.name.startswith("._"):
                    continue
                if entry.name in ("node_modules", "__pycache__", ".git"):
                    continue

                if entry.is_dir():
                    try:
                        child_count = sum(1 for c in entry.iterdir() if not c.name.startswith("."))
                    except PermissionError:
                        child_count = 0
                    items.append({
                        "name": entry.name,
                        "path": str(entry),
                        "type": "dir",
                        "children": child_count,
                    })
                elif entry.is_file():
                    try:
                        size = entry.stat().st_size
                    except OSError:
                        size = 0
                    ext = entry.suffix.lower()
                    if ext in (".btf", ".tif", ".tiff", ".h5", ".h5ad", ".parquet", ".zarr", ".yaml", ".json"):
                        units = ["B", "KB", "MB", "GB"]
                        s = float(size)
                        u = 0
                        while s >= 1024 and u < 3:
                            s /= 1024
                            u += 1
                        items.append({
                            "name": entry.name,
                            "path": str(entry),
                            "type": "file",
                            "size": size,
                            "size_human": f"{s:.1f} {units[u]}",
                        })
            except (PermissionError, OSError) as entry_err:
                logger.warning(f"跳過無法存取的項目 {entry.name}：{entry_err}")
                continue


        return {
            "status": "ok",
            "data": {
                "current": str(target),
                "parent": str(target.parent) if target.parent != target else None,
                "items": items,
            },
        }
    except Exception as e:
        logger.error(f"瀏覽目錄失敗：{e}")
        return {"status": "error", "message": str(e)}

