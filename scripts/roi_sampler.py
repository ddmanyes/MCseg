"""
組織感知 ROI 隨機抽樣器
讀取 tissue_positions.parquet，隨機抽取 3–5 個有效 ROI，寫入 results/qc_rois.json。
"""
from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

ROOT = Path(__file__).parent.parent
PX_SIZE           = 0.2737
ROI_UM            = 274
ROI_PX            = int(ROI_UM / PX_SIZE)
N_TARGET          = 4
COVERAGE_THRESH   = 0.60
COVERAGE_FALLBACK = 0.40
MAX_TRIES         = 50


def _build_tissue_index(pos: pd.DataFrame, stride: int) -> set[tuple[int, int]]:
    rows = (pos["pxl_row_in_fullres"].values.astype(int) // stride)
    cols = (pos["pxl_col_in_fullres"].values.astype(int) // stride)
    return set(zip(rows.tolist(), cols.tolist()))


def _coverage(x0: int, y0: int, tissue_set: set, stride: int) -> float:
    r0 = y0 // stride
    r1 = (y0 + ROI_PX) // stride + 1
    c0 = x0 // stride
    c1 = (x0 + ROI_PX) // stride + 1
    total = (r1 - r0) * (c1 - c0)
    if total == 0:
        return 0.0
    hit = sum(1 for r in range(r0, r1) for c in range(c0, c1) if (r, c) in tissue_set)
    return hit / total


def sample_rois(
    binned_dir: Path | str | None = None,
    n_target: int = N_TARGET,
    coverage_thresh: float = COVERAGE_THRESH,
    coverage_fallback: float = COVERAGE_FALLBACK,
    max_tries: int = MAX_TRIES,
    seed: int = 42,
    out_path: Path | str | None = None,
) -> dict:
    if binned_dir is None:
        cfg = yaml.safe_load((ROOT / "config" / "pipeline.yaml").read_text(encoding="utf-8"))
        binned_dir = Path(cfg["paths"]["binned_002"])
    binned_dir = Path(binned_dir)

    if out_path is None:
        out_path = ROOT / "results" / "qc_rois.json"
    out_path = Path(out_path)

    pos_path = binned_dir / "tissue_positions.parquet"
    pos      = pd.read_parquet(pos_path)
    tissue   = pos[pos["in_tissue"] == 1].copy()

    x_min = int(tissue["pxl_col_in_fullres"].min())
    x_max = int(tissue["pxl_col_in_fullres"].max())
    y_min = int(tissue["pxl_row_in_fullres"].min())
    y_max = int(tissue["pxl_row_in_fullres"].max())

    stride     = max(8, ROI_PX // 64)
    tissue_set = _build_tissue_index(tissue, stride)

    np.random.seed(seed)
    rois: list[dict] = []
    tries     = 0
    threshold = coverage_thresh

    while len(rois) < n_target and tries < max_tries:
        x0 = int(np.random.randint(x_min, max(x_min + 1, x_max - ROI_PX)))
        y0 = int(np.random.randint(y_min, max(y_min + 1, y_max - ROI_PX)))
        cov = _coverage(x0, y0, tissue_set, stride)
        if cov >= threshold:
            rois.append({
                "name":          f"qc_roi_{len(rois) + 1}",
                "x":             x0,
                "y":             y0,
                "width_px":      ROI_PX,
                "height_px":     ROI_PX,
                "pixel_size_um": PX_SIZE,
                "coverage":      round(float(cov), 3),
            })
        tries += 1
        if tries == 40 and len(rois) < 2:
            print(f"⚠️  覆蓋率閾值降至 {coverage_fallback}")
            threshold = coverage_fallback

    result = {
        "timestamp":      datetime.now().isoformat(),
        "threshold_used": threshold,
        "rois":           rois,
    }

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(result, indent=2), encoding="utf-8")

    print(f"✅ 抽樣完成：{len(rois)} 個 ROI，閾值={threshold}")
    for r in rois:
        print(f"  {r['name']}  x={r['x']} y={r['y']}  coverage={r['coverage']:.1%}")

    return result


if __name__ == "__main__":
    sample_rois()
