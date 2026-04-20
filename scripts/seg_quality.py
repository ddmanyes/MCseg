"""
NUC + MCseg 雙模式分割品質評估包裝器
讀取 results/qc_rois.json，對每個 ROI 執行 NUC（基準）與 MCseg（完整）分割，
儲存 results/qc/{roi_name}_nuc.npy 與 results/qc/{roi_name}_mcseg.npy。
"""
from __future__ import annotations

import json
import logging
from pathlib import Path

import numpy as np
import tifffile
import yaml

ROOT   = Path(__file__).parent.parent
logger = logging.getLogger(__name__)

NUC_CFG = {
    "use_gpu":             False,
    "batch_size":          4,
    "dia_small":           15.0,
    "dia_mid":             15.0,
    "dia_large":           15.0,
    "use_hematoxylin":     False,
    "use_cpsam":           False,
    "voronoi_distance":    0,
    "clahe_clip_limit":    0.0,
    "min_size":            20,
    "max_size":            6000,
    "flow_threshold":      0.4,
    "cellprob_threshold":  -2.0,
    "use_transcript_rescue": False,
}


def _read_roi_crop(he_path: Path, x0: int, y0: int, w: int, h: int) -> np.ndarray:
    """BTF tile-based 裁切，禁止全圖載入。"""
    with tifffile.TiffFile(str(he_path)) as tif:
        page = tif.pages[0]
        full = page.asarray()
        crop = full[y0:y0 + h, x0:x0 + w]
    if crop.ndim == 2:
        crop = np.stack([crop] * 3, axis=-1)
    return crop.astype(np.uint8)


def run_seg_quality(
    rois_json: Path | str | None = None,
    he_path:   Path | str | None = None,
    mcseg_cfg: dict | None = None,
    out_dir:   Path | str | None = None,
) -> None:
    """對每個 QC ROI 執行 NUC 與 MCseg 分割，儲存 .npy 遮罩。"""
    from backend.src.segmentation.cellpose_runner import run_mcseg_v2

    if rois_json is None:
        rois_json = ROOT / "results" / "qc_rois.json"
    if out_dir is None:
        out_dir = ROOT / "results" / "qc"
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    rois_data = json.loads(Path(rois_json).read_text(encoding="utf-8"))
    rois = rois_data["rois"]

    if he_path is None:
        cfg_raw = yaml.safe_load((ROOT / "config" / "pipeline.yaml").read_text(encoding="utf-8"))
        he_path = Path(cfg_raw["paths"]["he_image"])

    if mcseg_cfg is None:
        cfg_raw = yaml.safe_load((ROOT / "config" / "pipeline.yaml").read_text(encoding="utf-8"))
        mcseg_cfg = cfg_raw["segmentation"]["mcseg_v2"]

    for roi in rois:
        name = roi["name"]
        x0, y0 = roi["x"], roi["y"]
        w,  h  = roi["width_px"], roi["height_px"]

        try:
            img = _read_roi_crop(Path(he_path), x0, y0, w, h)
        except Exception as e:
            logger.warning(f"跳過 {name}：BTF 裁切失敗 — {e}")
            continue

        # NUC 基準線
        try:
            nuc_mask = run_mcseg_v2(img, NUC_CFG)
            np.save(out_dir / f"{name}_nuc.npy", nuc_mask)
            print(f"  {name} NUC  cells={nuc_mask.max()}")
        except Exception as e:
            logger.warning(f"跳過 {name} NUC：{e}")
            continue

        # MCseg 完整流程
        try:
            mcseg_mask = run_mcseg_v2(img, mcseg_cfg)
            np.save(out_dir / f"{name}_mcseg.npy", mcseg_mask)
            print(f"  {name} MCseg cells={mcseg_mask.max()}")
        except Exception as e:
            logger.warning(f"跳過 {name} MCseg：{e}")

    print(f"✅ seg_quality 完成，輸出至 {out_dir}")


if __name__ == "__main__":
    run_seg_quality()
