"""
全圖 AnnData 組裝器
讀取 handoff_report.json，對所有 MCseg 遮罩執行 RNA 計數，
合併後寫入 results/analysis/cellpose_cells.h5ad。
"""
from __future__ import annotations

import gc
import json
import logging
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import yaml


def _find_root(start: Path) -> Path:
    for p in [start, start.parent, start.parent.parent]:
        if (p / "pyproject.toml").exists():
            return p
    return start


ROOT   = _find_root(Path(__file__).resolve().parent)
logger = logging.getLogger(__name__)


def build_adata(
    handoff_json: Path | str | None = None,
    rois_json:    Path | str | None = None,
    out_path:     Path | str | None = None,
    dilation_px:  int = 6,
) -> ad.AnnData:
    from backend.src.cellpose_counter.counter import count_rna_per_cell

    if handoff_json is None:
        handoff_json = ROOT / "results" / "handoff_report.json"
    if rois_json is None:
        rois_json = ROOT / "results" / "qc_rois.json"
    if out_path is None:
        out_path = ROOT / "results" / "analysis" / "cellpose_cells.h5ad"

    report    = json.loads(Path(handoff_json).read_text(encoding="utf-8"))
    masks_dir = Path(report["masks_dir"])
    binned_dir = Path(report["binned_dir"])

    rois = json.loads(Path(rois_json).read_text(encoding="utf-8"))["rois"]

    adata_path = binned_dir / "adata_002um.h5ad"
    if not adata_path.exists():
        raise FileNotFoundError(
            f"找不到 bin 矩陣：{adata_path}\n"
            "請先執行 Stage 0 產出 adata_002um.h5ad，或確認 binned_dir 路徑正確。"
        )

    adatas = []
    for roi in rois:
        name = roi["name"]
        mask_path = masks_dir / f"{name}_mcseg.npy"
        if not mask_path.exists():
            logger.warning(f"跳過 {name}：遮罩不存在 ({mask_path})")
            continue

        try:
            a = count_rna_per_cell(
                adata_path   = adata_path,
                mask_path    = mask_path,
                roi_x_px     = roi["x"],
                roi_y_px     = roi["y"],
                pixel_size_um= roi.get("pixel_size_um", 0.2737),
                dilation_px  = dilation_px,
            )
        except Exception as e:
            logger.warning(f"跳過 {name}：計數失敗 — {e}")
            continue

        # 負座標驗證
        if "spatial" in a.obsm:
            coords = a.obsm["spatial"]
            neg_mask = coords[:, 0] < 0
            if neg_mask.any():
                logger.warning(f"⚠️  {name} 偏移後出現 {neg_mask.sum()} 個負座標")

        a.obs["roi_name"] = name
        adatas.append(a)
        gc.collect()

    if not adatas:
        raise RuntimeError("所有 ROI 均失敗，無法建立 AnnData")

    combined = ad.concat(adatas, merge="same") if len(adatas) > 1 else adatas[0]

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    combined.write(out_path)
    print(f"[OK] AnnData saved: {out_path}  ({combined.n_obs} cells x {combined.n_vars} genes)")
    return combined


if __name__ == "__main__":
    build_adata()
