"""Stage 2：Cellpose RNA 計數 API"""
import asyncio
import logging
from fastapi import APIRouter, BackgroundTasks
from pydantic import BaseModel
from typing import Optional

from backend.src.utils.config import load_config
from backend.src.utils.logging import set_current_stage

router = APIRouter()
logger = logging.getLogger("pipeline.api.cellpose_count")

_status = {"status": "idle", "progress": 0.0, "message": ""}
_lock   = asyncio.Lock()

_full_status = {"status": "idle", "progress": 0.0, "message": ""}
_full_lock   = asyncio.Lock()

FULL_COUNT_SUBDIR = "fullslide"
FULL_COUNT_FILENAME = "cells.h5ad"


class CountParams(BaseModel):
    roi_name: Optional[str] = None   # None = 全部 ROI


@router.get("/status")
async def get_status():
    global _status
    # 記憶體 idle 時，查磁碟是否已有 cellpose_cells.h5ad
    if _status["status"] == "idle":
        try:
            config = load_config()
            from backend.src.utils.config import resolve_path
            out_base = resolve_path(config["paths"]["output_dir"]) / "roi"
            rois = config.get("rois", [])
            if rois and any((out_base / r.get("name", "") / "cellpose_cells.h5ad").exists() for r in rois):
                _status = {"status": "done", "progress": 1.0, "message": "Count complete (restored from disk)"}
        except Exception:
            pass
    return _status


async def _run_count(config: dict, roi_name: Optional[str]):
    global _status
    set_current_stage("count")
    _status = {"status": "running", "progress": 0.0, "message": "分配 RNA 至 Cellpose 細胞..."}
    try:
        from backend.src.cellpose_counter.counter import run_counting_pipeline
        await asyncio.get_running_loop().run_in_executor(
            None, run_counting_pipeline, config, roi_name
        )
        _status = {"status": "done", "progress": 1.0, "message": "RNA counting complete"}
    except Exception as e:
        logger.error(f"Cellpose 計數失敗：{e}")
        _status = {"status": "error", "progress": 0.0, "message": str(e)}


@router.post("/run")
async def run_count(background_tasks: BackgroundTasks, params: Optional[CountParams] = None):
    async with _lock:
        if _status["status"] == "running":
            return {"status": "error", "message": "任務執行中"}
        _status["status"] = "running"

    config = load_config()
    roi_name = params.roi_name if params else None
    background_tasks.add_task(_run_count, config, roi_name)
    return {"status": "ok", "message": "RNA 計數已啟動"}


@router.get("/full_status")
async def get_full_count_status():
    return _full_status


@router.post("/run_full")
async def run_full_count(background_tasks: BackgroundTasks):
    """
    對全圖分割遮罩執行 RNA 計數，輸出 `{output_dir}/fullslide/cells.h5ad`。

    走 `fullslide.pipeline` 共用流程（bin_attribution → aggregate_cells →
    add_centroids），與 CLI 全片流程同一份實作。
    """
    from backend.src.fullslide.pipeline import resolve_full_count_inputs

    config = load_config()
    inputs, err = resolve_full_count_inputs(config)
    if err:
        return {"status": "error", "message": err}

    async with _full_lock:
        if _full_status["status"] == "running":
            return {"status": "error", "message": "全圖計數任務執行中"}
        _full_status.update({"status": "running", "progress": 0.0, "message": "初始化..."})
        background_tasks.add_task(_run_full_count, inputs)
    return {"status": "ok", "message": "全圖 RNA 計數已啟動"}


def _full_count_sync(inputs: dict, progress) -> tuple[int, int, dict]:
    """同步執行全圖計數（在 executor 中跑），回傳 (n_cells, n_genes, coverage)。"""
    import numpy as np

    from backend.src.fullslide.pipeline import (
        add_centroids,
        aggregate_cells,
        bin_attribution,
    )

    progress(0.05, "載入全圖遮罩...")
    mask = np.load(str(inputs["mask_path"]))

    # 等距擴張填補 Voronoi 間隙；與 ROI 路徑（counter.py）用同一個設定值，
    # 否則全圖與 ROI 的計數結果不可比。實測此步對 bin 命中率影響極大
    # （dpcp01 全片：20.9% → 35.3%）。
    dilation_px = inputs["dilation_px"]
    if dilation_px > 0:
        from skimage.segmentation import expand_labels

        um = dilation_px * inputs["pixel_size_um"]
        progress(0.15, f"等距擴張 {dilation_px}px（≈{um:.2f}µm）...")
        mask = expand_labels(mask, distance=dilation_px)

    progress(0.25, f"對應 2µm bins 至細胞（{inputs['transform_source']}）...")
    origin_x, origin_y = inputs["origin_xy"]
    attribution = bin_attribution(
        mask,
        inputs["tp_path"],
        crop_y0=origin_y,
        crop_x0=origin_x,
        scale=inputs["scale"],
        transform=inputs.get("transform"),
    )
    if attribution.empty:
        raise ValueError(
            "沒有任何 bin 落在細胞內 —— 請檢查影像與 Visium 資料是否對位"
        )

    coverage = attribution.attrs.get("coverage", {})
    frac_out = coverage.get("frac_out_of_image", 0.0)
    if frac_out:
        logger.info(
            f"bin 涵蓋率：{coverage['n_in_bounds']:,}/{coverage['n_total']:,} 落在影像內"
            f"（{frac_out:.1%} 在影像範圍外），其中 {coverage['n_assigned']:,} 個落在細胞上"
        )

    progress(0.55, f"聚合 {len(attribution):,} 個 bins 為細胞矩陣...")
    cells = aggregate_cells(attribution, inputs["h5_path"])

    progress(0.85, "計算細胞重心...")
    add_centroids(
        cells, mask, inputs["pixel_size_um"], origin_xy=inputs["origin_xy"]
    )

    out_dir = inputs["output_dir"] / FULL_COUNT_SUBDIR
    out_dir.mkdir(parents=True, exist_ok=True)
    cells.write_h5ad(str(out_dir / FULL_COUNT_FILENAME))
    return cells.n_obs, cells.n_vars, coverage


async def _run_full_count(inputs: dict) -> None:
    global _full_status
    set_current_stage("count")

    def _progress(p: float, msg: str) -> None:
        _full_status.update({"progress": p, "message": msg})

    try:
        if inputs["meta_missing"]:
            logger.warning(
                "找不到 full_image_segmentation_meta.json，裁切原點視為 (0, 0)。"
                "若分割時有指定裁切窗格，座標將會偏移。"
            )
        n_cells, n_genes, coverage = await asyncio.get_running_loop().run_in_executor(
            None, _full_count_sync, inputs, _progress
        )
        frac_out = coverage.get("frac_out_of_image", 0.0)
        # 落在影像外是資料的真實限制（例如高解析圖沒掃到 SR 畫布右緣），
        # 必須明說而非靜默吞掉 —— 使用者才能判斷該不該補掃或換圖。
        oob = f"；{frac_out:.1%} bins 落在影像範圍外" if frac_out >= 0.001 else ""
        _full_status = {
            "status": "done", "progress": 1.0,
            "message": f"全圖計數完成：{n_cells:,} cells × {n_genes:,} genes"
                       f"（{inputs['transform_source']}{oob}）"
                       f"  →  {FULL_COUNT_SUBDIR}/{FULL_COUNT_FILENAME}",
            "n_cells": n_cells,
            "coverage": coverage,
        }
    except Exception as e:
        logger.error(f"全圖 RNA 計數失敗：{e}", exc_info=True)
        # ValueError 為流程可預期的診斷訊息（不含路徑），可安全回傳
        safe_msg = str(e) if isinstance(e, ValueError) else "全圖 RNA 計數失敗，請查閱 log"
        _full_status = {"status": "error", "progress": 0.0, "message": safe_msg}


@router.post("/coverage_qc")
async def run_coverage_qc(
    grid_px: int = 2048,
    # None = 依網格面積自動推算（固定 200 在 2048px 網格只佔滿格的 0.26%，
    # 等於不過濾，會讓四分之一的網格被標記）
    min_bins: int | None = None,
    low_ratio: float = 0.3,
):
    """
    全片分割覆蓋率 QC：逐網格比對 bin 密度 vs 細胞密度，標記分割失敗的區域。

    回答總命中率看不出來的問題 —— 一個區域整片漏掉只讓總數低幾個百分點
    （真實案例：Day0/Day3 區只分出 380 顆碎片，該區有 31,969 bins）。

    同步執行：只讀 memmap 遮罩與 parquet，全片約數十秒，不需背景任務。
    """
    from backend.src.fullslide.coverage_qc import run_coverage_qc_from_config

    try:
        sections, grid = run_coverage_qc_from_config(
            load_config(), grid_px=grid_px, min_bins=min_bins, low_ratio=low_ratio
        )
    except ValueError as e:
        # resolve_full_count_inputs 的診斷訊息（不含絕對路徑），可安全回傳
        return {"status": "error", "message": str(e)}
    except Exception as e:
        logger.error(f"覆蓋率 QC 失敗：{e}", exc_info=True)
        return {"status": "error", "message": "覆蓋率 QC 失敗，請查閱 log"}

    return {
        "status": "ok",
        "data": {
            # 切片層級是主要結論（網格層級抓不到整片缺口，見 coverage_qc 模組說明）
            "sections": sections.to_dict(orient="records"),
            "sections_summary": sections.attrs.get("summary", {}),
            "grid_summary": grid.attrs["summary"],
            # 網格只回傳被標記的（全片可有數百格），用途是在切片內定位
            "grid_flagged": grid[grid["flagged"]].to_dict(orient="records"),
        },
    }


@router.get("/available_rois")
async def get_available_rois():
    """列出所有已有 cellpose_cells.h5ad 的 ROI"""
    try:
        config = load_config()
        from backend.src.utils.config import resolve_path
        out_base = resolve_path(config["paths"]["output_dir"]) / "roi"
        rois = config.get("rois", [])
        result = []
        for roi in rois:
            name = roi.get("name", "")
            has_mask  = (out_base / name / "segmentation_masks.npy").exists()
            has_count = (out_base / name / "cellpose_cells.h5ad").exists()
            result.append({
                "name":      name,
                "has_mask":  has_mask,
                "has_count": has_count,
            })
        return {"status": "ok", "data": result}
    except Exception as e:
        return {"status": "error", "message": str(e)}
