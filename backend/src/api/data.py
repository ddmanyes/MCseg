"""資料設定 API：掃描資料目錄、套用路徑配置、目錄瀏覽"""
import logging
import os
import string
import sys
from pathlib import Path
from typing import Optional

from fastapi import APIRouter, Query
from pydantic import BaseModel

from backend.src.utils.config import load_config, resolve_path, save_state
from backend.src.utils.discovery import scan_data_root

router = APIRouter()
logger = logging.getLogger("pipeline.api.data")


class ScanRequest(BaseModel):
    data_root: str


class ApplyRequest(BaseModel):
    data_root: Optional[str] = None
    he_image: Optional[str] = None
    binned_002: Optional[str] = None
    binned_008: Optional[str] = None
    output_dir: Optional[str] = None
    pixel_size_um: Optional[float] = None


@router.post("/scan")
async def scan_directory(req: ScanRequest):
    """掃描指定目錄，回傳發現的檔案清單"""
    try:
        result = scan_data_root(req.data_root)
        return {"status": "ok", "data": result.to_dict()}
    except Exception as e:
        logger.error(f"掃描失敗：{e}", exc_info=True)
        return {"status": "error", "message": "掃描失敗，請查閱 log"}


@router.post("/apply")
async def apply_paths(req: ApplyRequest):
    """將發現的路徑寫入 pipeline.yaml"""
    try:
        config = load_config()
        paths = config.setdefault("paths", {})
        updates = req.model_dump(exclude_none=True)

        pixel_size_um = updates.pop("pixel_size_um", None)

        for key, value in updates.items():
            if value is not None and value != "":
                paths[key] = value

        # 只有當 request 明確帶入 data_root 時才自動設定輸出目錄
        data_root_in_request = updates.get("data_root")
        if data_root_in_request:
            result_base = Path(os.path.expanduser(data_root_in_request)) / "MCseg_result"
            paths["output_dir"] = str(result_base / "analysis")
            paths["export_dir"] = str(result_base / "export")
            paths["figure_dir"] = str(result_base / "figures")
            logger.info(f"輸出目錄設為：{result_base}")
        elif "output_dir" in updates and updates["output_dir"]:
            logger.info(f"手動指定 output_dir：{updates['output_dir']}")

        state_update: dict = {"paths": config["paths"]}
        if pixel_size_um is not None:
            config.setdefault("global", {})["pixel_size_um"] = pixel_size_um
            state_update["global"] = config["global"]
            logger.info(f"已更新 pixel_size_um = {pixel_size_um}")

        # 先存 config，再建立目錄
        save_state(state_update)

        if data_root_in_request:
            for subdir in ("analysis", "export", "figures"):
                (result_base / subdir).mkdir(parents=True, exist_ok=True)
        elif "output_dir" in updates and updates["output_dir"]:
            resolve_path(updates["output_dir"]).mkdir(parents=True, exist_ok=True)
        logger.info(f"已套用 {len(updates)} 項路徑設定")
        return {"status": "ok", "message": f"已更新 {len(updates)} 項路徑", "data": paths}
    except Exception as e:
        logger.error(f"套用失敗：{e}", exc_info=True)
        return {"status": "error", "message": "套用路徑失敗，請查閱 log"}


@router.get("/output-dir")
async def get_output_dir():
    """取得目前輸出目錄設定"""
    try:
        config = load_config()
        paths = config.get("paths", {})
        output_dir = paths.get("output_dir", "results/analysis")
        resolved = str(resolve_path(output_dir))
        return {"status": "ok", "data": {"output_dir": output_dir, "resolved": resolved}}
    except Exception as e:
        logger.error(f"取得輸出目錄失敗：{e}", exc_info=True)
        return {"status": "error", "message": "取得輸出目錄失敗，請查閱 log"}


@router.get("/status")
async def get_data_status():
    """取得目前 paths 配置狀態（哪些已填、哪些為空）"""
    try:
        config = load_config()
        paths = config.get("paths", {})
        required_keys = ["he_image", "binned_002", "binned_008"]
        status = {}
        for key in required_keys:
            val = paths.get(key, "")
            status[key] = {
                "path": val,
                "configured": bool(val),
            }
        return {"status": "ok", "data": status}
    except Exception as e:
        logger.error(f"取得資料狀態失敗：{e}", exc_info=True)
        return {"status": "error", "message": "取得資料狀態失敗，請查閱 log"}


@router.get("/disk-status")
async def get_disk_status():
    """
    掃描磁碟，回傳各 Stage 實際完成狀態（重啟後前端可據此恢復 UI 狀態）。
    """
    try:
      config = load_config()
    except Exception as e:
        logger.error(f"disk-status 載入設定失敗：{e}", exc_info=True)
        return {"status": "error", "message": "磁碟狀態查詢失敗，請查閱 log"}
    paths = config.get("paths", {})
    output_dir = resolve_path(paths.get("output_dir", "results/analysis"))

    roi_base = output_dir / "roi"

    # Stage 0: ROI — 有任何 he_crop.tif
    if roi_base.exists():
        try:
            roi_dirs = [d for d in roi_base.iterdir() if d.is_dir() and not d.name.startswith(".")]
        except PermissionError:
            roi_dirs = []
    else:
        roi_dirs = []
    roi_done = any((d / "he_crop.tif").exists() for d in roi_dirs)
    roi_names = [d.name for d in roi_dirs if (d / "he_crop.tif").exists()]

    # Stage 1: Segmentation — 有任何 segmentation_masks.npy
    seg_done = any((d / "segmentation_masks.npy").exists() for d in roi_dirs)

    # Stage 2: RNA 計數 — 有任何 cellpose_cells.h5ad
    count_done = any((d / "cellpose_cells.h5ad").exists() for d in roi_dirs)

    # Stage 2.5: Proseg RNA 重分配 — 有任何 proseg_cells.h5ad
    proseg_rna_done = any((d / "proseg_cells.h5ad").exists() for d in roi_dirs)

    # Stage 3: Analysis — 有任何 qc_preprocessed.h5ad 或 umap_computed.h5ad
    analysis_done = (
        (output_dir / "umap_computed.h5ad").exists() or
        (output_dir / "qc_preprocessed.h5ad").exists() or
        any((d / "clustering.h5ad").exists() for d in roi_dirs)
    )

    return {
        "status": "ok",
        "data": {
            "roi":         {"done": roi_done,  "roi_names": roi_names},
            "segmentation":{"done": seg_done},
            "count":       {"done": count_done},
            "proseg_rna":  {"done": proseg_rna_done},
            "analysis":    {"done": analysis_done},
        },
    }


def _build_allowed_roots() -> list[Path]:
    home = Path(os.path.expanduser("~"))
    if sys.platform == "win32":
        drives = [Path(f"{d}:/") for d in string.ascii_uppercase if Path(f"{d}:/").exists()]
        return [home, *drives]
    return [home, Path("/Volumes"), Path("/tmp")]

_BROWSE_ALLOWED_ROOTS = _build_allowed_roots()


@router.get("/browse")
async def browse_directory(path: str = Query("~", description="要瀏覽的目錄路徑")):
    """
    瀏覽本機目錄結構。回傳指定路徑下的子目錄與大型檔案。
    前端可用此 API 實現點擊式的資料夾選擇器。
    僅限 ~/、/Volumes/、/tmp/ 底下的路徑。
    """
    try:
        target = Path(os.path.expanduser(path)).resolve()
        if not any(str(target).startswith(str(r)) for r in _BROWSE_ALLOWED_ROOTS):
            return {"status": "error", "message": "路徑不在允許範圍內"}
        if not target.exists():
            return {"status": "error", "message": "路徑不存在"}
        if not target.is_dir():
            return {"status": "error", "message": "不是目錄"}

        items = []
        try:
            entries = sorted(target.iterdir(), key=lambda e: (not e.is_dir(), e.name.lower()))
        except PermissionError:
            return {"status": "error", "message": "無權限存取此目錄"}

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
        logger.error(f"瀏覽目錄失敗：{e}", exc_info=True)
        return {"status": "error", "message": "瀏覽目錄失敗，請查閱 log"}

