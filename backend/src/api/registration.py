"""對位配準 API（估計殘餘位移、產生 QC 疊圖）

**定位**：主要對位由 Stage 2 全圖計數的組合 homography 負責（`registration.alignment`）。
本 API 是**驗證與診斷**入口 —— 回答「對位到底有沒有問題、還差多少」。
"""
from __future__ import annotations

import base64
import logging
from pathlib import Path
from typing import Optional

from fastapi import APIRouter
from pydantic import BaseModel

from backend.src.utils.config import load_config, resolve_path, save_state

router = APIRouter()
logger = logging.getLogger("pipeline.api.registration")

QC_SUBDIR = "qc/alignment"


class QCPatchParams(BaseModel):
    n: int = 3
    size: int = 512
    seed: int = 0


class ApplyParams(BaseModel):
    roi_name: Optional[str] = None
    downsample: int = 32
    # 預設 false：壞估計被靜默套用比不修正更糟，一律要求使用者看過 residual 再啟用
    enable: bool = False


def _resolve_inputs(config: dict) -> tuple[Path, Path, tuple[int, int]] | None:
    """取 (he_image, tissue_positions, full_shape)；缺件回 None。"""
    import tifffile

    paths = config.get("paths", {})
    he = Path(paths.get("he_image", ""))
    binned = paths.get("binned_002", "")
    tp = Path(binned) / "spatial" / "tissue_positions.parquet" if binned else Path()
    if not he.exists() or not tp.exists():
        return None

    with tifffile.TiffFile(str(he)) as tf:
        page = tf.pages[0]
        shape = (int(page.imagelength), int(page.imagewidth))
    return he, tp, shape


def _resolve_transform(config: dict, full_shape: tuple[int, int]):
    """取用於計數的 bin→影像 變換（與 Stage 2 同一份決策邏輯）。"""
    from backend.src.fullslide.pipeline import resolve_bin_to_image_transform

    transform, scale, source = resolve_bin_to_image_transform(config, full_shape)
    if transform is None:
        import numpy as np

        # 近似縮放也表達成 3×3，讓下游只需處理一種型別
        transform = np.diag([scale[0], scale[1], 1.0])
    return transform, source


@router.get("/estimate")
async def estimate_alignment(roi_name: Optional[str] = None):
    """
    估計 bin 相對 H&E 影像的**殘餘**位移（已套用對位變換之後）。

    回傳 `spread` 作為可信度指標：三個窗格的估計愈一致愈可信。
    `spread` 偏大代表不是單純平移（可能有旋轉／縮放，或組織有形變）。
    """
    config = load_config()
    resolved = _resolve_inputs(config)
    if resolved is None:
        return {"status": "error", "message": "找不到 H&E 影像或 tissue_positions.parquet"}
    he, tp, full_shape = resolved

    try:
        from backend.src.registration.align import estimate_shift_fullres

        transform, source = _resolve_transform(config, full_shape)
        dy, dx, spread = estimate_shift_fullres(
            he, tp, full_shape=full_shape, transform=transform
        )
    except (ValueError, OSError, NotImplementedError) as e:
        logger.warning(f"位移估計失敗：{e}")
        return {"status": "error", "message": f"位移估計失敗：{e}"}

    return {
        "status": "ok",
        "data": {
            "dy": round(dy, 2),
            "dx": round(dx, 2),
            "spread": None if spread == float("inf") else round(spread, 2),
            "transform_source": source,
            "roi_name": roi_name,
            "unit": "fullres_px",
        },
    }


@router.post("/qc_patches")
async def make_qc_patches(params: Optional[QCPatchParams] = None):
    """產生 QC 疊圖 PNG（同步執行；n=3、size=512 實測數秒內完成）。"""
    params = params or QCPatchParams()
    config = load_config()
    resolved = _resolve_inputs(config)
    if resolved is None:
        return {"status": "error", "message": "找不到 H&E 影像或 tissue_positions.parquet"}
    he, tp, full_shape = resolved

    out_dir = resolve_path(config["paths"]["output_dir"]) / QC_SUBDIR
    try:
        from backend.src.registration.qc import render_overlay_patches

        transform, source = _resolve_transform(config, full_shape)
        written = render_overlay_patches(
            he, tp, out_dir,
            n=params.n, size=params.size, seed=params.seed,
            full_shape=full_shape, transform=transform,
        )
    except (ValueError, OSError, NotImplementedError) as e:
        logger.warning(f"QC 疊圖產生失敗：{e}")
        return {"status": "error", "message": f"QC 疊圖產生失敗：{e}"}

    return {
        "status": "ok",
        "data": {"n": len(written), "transform_source": source},
        "message": f"已產生 {len(written)} 張對位 QC 疊圖",
    }


@router.get("/qc_images")
async def get_qc_images():
    """回傳既有 QC 疊圖的 base64（前端直接 <img src=...> 顯示）。"""
    out_dir = resolve_path(load_config()["paths"]["output_dir"]) / QC_SUBDIR
    if not out_dir.exists():
        return {"status": "ok", "data": []}

    images = []
    for p in sorted(out_dir.glob("patch_*.png")):
        try:
            b64 = base64.b64encode(p.read_bytes()).decode("ascii")
        except OSError as e:
            logger.warning(f"跳過 {p.name}：{e}")
            continue
        images.append({"name": p.stem, "data": f"data:image/png;base64,{b64}"})
    return {"status": "ok", "data": images}


@router.post("/apply")
async def apply_alignment(params: Optional[ApplyParams] = None):
    """
    估計仿射修正並寫入 `state.json` 的 `alignment`。

    **不自動啟用**：`enable` 預設 false，寫入的 `enabled` 亦為 false。前端顯示
    矩陣與 residual 之後，由使用者按「套用」才以 `enable=true` 再呼叫一次 ——
    估計錯誤而被靜默套用，比完全不修正更糟。
    """
    params = params or ApplyParams()
    config = load_config()
    resolved = _resolve_inputs(config)
    if resolved is None:
        return {"status": "error", "message": "找不到 H&E 影像或 tissue_positions.parquet"}
    he, tp, full_shape = resolved

    try:
        from backend.src.registration.align import (
            estimate_affine,
            render_bin_density,
            tissue_gray,
        )
        from backend.src.roi.tile_server import _load_or_build_thumb

        transform, source = _resolve_transform(config, full_shape)
        ds = params.downsample
        ref = tissue_gray(_load_or_build_thumb(he, ds))
        mov = render_bin_density(tp, full_shape, ds, transform)
        hh, ww = min(ref.shape[0], mov.shape[0]), min(ref.shape[1], mov.shape[1])

        affine = estimate_affine(ref[:hh, :ww], mov[:hh, :ww])
    except (ValueError, OSError, NotImplementedError) as e:
        logger.warning(f"仿射估計失敗：{e}")
        return {"status": "error", "message": f"仿射估計失敗：{e}"}

    # 估計是在 1/ds 縮圖上做的 → 平移項需乘回 fullres px（線性部分與尺度無關）
    m = affine.array.copy()
    m[:, 2] *= ds
    record = {
        "enabled": bool(params.enable),
        "matrix": [[float(v) for v in row] for row in m],
        "source": "spaceranger_fullres",
        "target": "raw_btf",
        "estimated_error": affine.estimated_error,
    }
    save_state({"alignment": record})

    return {
        "status": "ok",
        "data": {
            **record,
            "downsample": ds,
            "transform_source": source,
            "translation_px": [round(float(m[0, 2]), 2), round(float(m[1, 2]), 2)],
        },
        "message": (
            "已寫入 state.json（尚未啟用，請確認 residual 後再套用）"
            if not params.enable else "已寫入並啟用對位修正"
        ),
    }
