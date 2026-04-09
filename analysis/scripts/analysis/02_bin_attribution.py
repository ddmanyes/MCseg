"""
01_bin_attribution.py
=====================
將 Visium HD 2µm bins 歸屬到細胞遮罩：
  mask[row_local, col_local] → cell_id per bin

輸出：results/attribution/{method}_{roi}.parquet
欄位：barcode, pxl_row_in_fullres, pxl_col_in_fullres, cell_id
（cell_id=0 表示未歸屬到任何細胞）
"""

from __future__ import annotations

import sys
import json
import numpy as np
import pandas as pd
from pathlib import Path
from tqdm import tqdm

ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT))

import yaml
cfg = yaml.safe_load((ROOT / "config.yaml").read_text())
PATHS = cfg["paths"]
DATA  = cfg["data"]

MASKS_DIR      = ROOT / PATHS["masks_dir"]
ATTRIBUTION_DIR = ROOT / PATHS["attribution_dir"]
ATTRIBUTION_DIR.mkdir(parents=True, exist_ok=True)

with open(PATHS["roi_info"]) as f:
    ROI_INFO = json.load(f)

METHODS = DATA["methods"]


def load_tissue_positions() -> pd.DataFrame:
    """讀取 Visium HD 2µm tissue_positions.parquet（只保留 in_tissue bins）。"""
    tp_path = Path(PATHS["visiumhd_positions"])
    print(f"讀取 tissue_positions: {tp_path.name}")
    tp = pd.read_parquet(tp_path, columns=[
        "barcode", "in_tissue",
        "pxl_row_in_fullres", "pxl_col_in_fullres"
    ])
    tp = tp[tp["in_tissue"] == 1].reset_index(drop=True)
    print(f"  in_tissue bins: {len(tp):,}")
    return tp


def attribute_roi(tp: pd.DataFrame, mask: np.ndarray,
                  x0: int, y0: int, x1: int, y1: int) -> pd.DataFrame:
    """
    向量化歸屬（無 Python 迴圈）：
      - 過濾 ROI bbox 內的 bins
      - 計算 ROI 局部座標
      - mask lookup → cell_id
    """
    # 過濾 ROI bbox
    in_roi = (
        (tp["pxl_col_in_fullres"] >= x0) & (tp["pxl_col_in_fullres"] < x1) &
        (tp["pxl_row_in_fullres"] >= y0) & (tp["pxl_row_in_fullres"] < y1)
    )
    tp_roi = tp[in_roi].copy()

    if tp_roi.empty:
        return tp_roi

    # ROI 局部座標
    row_local = (tp_roi["pxl_row_in_fullres"].values - y0).astype(np.int32)
    col_local = (tp_roi["pxl_col_in_fullres"].values - x0).astype(np.int32)

    # 邊界 clip（防止座標略超出）
    H, W = mask.shape
    row_local = row_local.clip(0, H - 1)
    col_local = col_local.clip(0, W - 1)

    tp_roi["cell_id"] = mask[row_local, col_local]
    return tp_roi[["barcode", "pxl_row_in_fullres", "pxl_col_in_fullres", "cell_id"]]


def run_attribution():
    tp = load_tissue_positions()

    for method in METHODS:
        print(f"\n[{method.upper()}] bin attribution...")

        for roi_name, roi in tqdm(ROI_INFO.items(), desc=f"{method}"):
            dst = ATTRIBUTION_DIR / f"{method}_{roi_name}.parquet"
            if dst.exists():
                print(f"  {roi_name} 已存在，跳過")
                continue

            mask_path = MASKS_DIR / f"{method}_{roi_name}.npy"
            if not mask_path.exists():
                print(f"  ⚠️  遮罩不存在：{mask_path.name}，跳過")
                continue

            mask = np.load(mask_path)
            x0, y0, x1, y1 = roi["x0"], roi["y0"], roi["x1"], roi["y1"]

            result = attribute_roi(tp, mask, x0, y0, x1, y1)
            result.to_parquet(dst, index=False)

            n_attributed = (result["cell_id"] > 0).sum()
            capture_rate = n_attributed / len(result) if len(result) > 0 else 0
            print(f"  {roi_name}: {len(result):,} bins，歸屬率 {capture_rate:.1%}")

    print("\n✅ 01_bin_attribution.py 完成")
    print(f"   輸出：{ATTRIBUTION_DIR}")


if __name__ == "__main__":
    run_attribution()
