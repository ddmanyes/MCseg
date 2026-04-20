"""
FTC / NED / Artificial Co-expression Rate 計算器
讀取 results/qc/ 下的 .npy 遮罩，輸出 results/qc_metrics.csv。
"""
from __future__ import annotations

import gc
import json
import logging
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import yaml
from skimage.morphology import dilation, footprint_rectangle
from skimage.segmentation import expand_labels


def _find_root(start: Path) -> Path:
    for p in [start, start.parent, start.parent.parent]:
        if (p / "pyproject.toml").exists():
            return p
    return start


ROOT   = _find_root(Path(__file__).resolve().parent)
logger = logging.getLogger(__name__)

IMPOSSIBLE_PAIRS_MAP = {
    "crc":  [["EPCAM", "CD3E"], ["MUC2", "NKG7"], ["ACTA2", "CD3E"], ["PECAM1", "EPCAM"]],
    "luad": [["EPCAM", "CD3E"], ["SFTPC", "NKG7"], ["ACTA2", "CD3E"], ["PECAM1", "EPCAM"]],
}


def _hellinger(p: np.ndarray, q: np.ndarray) -> float:
    sp_sum = p.sum()
    sq_sum = q.sum()
    if sp_sum == 0 or sq_sum == 0:
        return 0.0
    p = p / sp_sum
    q = q / sq_sum
    return float(0.5 * np.sum((np.sqrt(p) - np.sqrt(q)) ** 2))


def _build_count_matrix(
    mask: np.ndarray,
    adata_roi: "sc.AnnData",
    roi_x: int,
    roi_y: int,
    dilation_px: int = 6,
) -> sp.csr_matrix:
    if dilation_px > 0:
        mask = expand_labels(mask, distance=dilation_px)

    n_cells = int(mask.max())
    if n_cells == 0:
        return sp.csr_matrix((0, adata_roi.n_vars))

    col = (adata_roi.obs["pxl_col_in_fullres"].values.astype(int) - roi_x)
    row = (adata_roi.obs["pxl_row_in_fullres"].values.astype(int) - roi_y)
    h, w = mask.shape
    valid = (col >= 0) & (col < w) & (row >= 0) & (row < h)

    col_v, row_v = col[valid], row[valid]
    cell_ids = mask[row_v, col_v]
    assigned = cell_ids > 0

    if assigned.sum() == 0:
        return sp.csr_matrix((n_cells, adata_roi.n_vars))

    bin_idx = np.where(valid)[0][assigned]
    cell_idx = cell_ids[assigned] - 1

    X = adata_roi.X
    M = sp.csr_matrix(
        (np.ones(assigned.sum()), (cell_idx, bin_idx)),
        shape=(n_cells, adata_roi.n_obs),
    )
    return M @ (X if sp.issparse(X) else sp.csr_matrix(X))


def compute_metrics(
    rois_json:    Path | str | None = None,
    binned_dir:   Path | str | None = None,
    qc_dir:       Path | str | None = None,
    out_csv:      Path | str | None = None,
    tissue_profile: str | None = None,
    ned_sample:   int = 200,
) -> pd.DataFrame:
    if rois_json is None:
        rois_json = ROOT / "results" / "qc_rois.json"
    if out_csv is None:
        out_csv = ROOT / "results" / "qc_metrics.csv"
    if qc_dir is None:
        qc_dir = ROOT / "results" / "qc"

    cfg_raw = yaml.safe_load((ROOT / "config" / "pipeline.yaml").read_text(encoding="utf-8"))
    if binned_dir is None:
        binned_dir = Path(cfg_raw["paths"]["binned_002"])
    if tissue_profile is None:
        tissue_profile = cfg_raw["global"].get("tissue_profile", "crc")

    binned_dir = Path(binned_dir)
    qc_dir     = Path(qc_dir)

    impossible_pairs = IMPOSSIBLE_PAIRS_MAP.get(tissue_profile, IMPOSSIBLE_PAIRS_MAP["crc"])
    rois = json.loads(Path(rois_json).read_text(encoding="utf-8"))["rois"]

    # 空間座標（先讀，輕量）
    pos_path = binned_dir / "tissue_positions.parquet"
    pos_df   = pd.read_parquet(pos_path) if pos_path.exists() else None

    # 讀取 bin 矩陣：h5ad 用 backed='r' 避免全圖載入 RAM
    mtx_path = binned_dir / "filtered_feature_bc_matrix"
    if mtx_path.exists():
        adata_full = sc.read_10x_mtx(str(mtx_path), var_names="gene_symbols", cache=False)
        if pos_df is not None:
            adata_full.obs = adata_full.obs.join(
                pos_df.set_index("barcode")[["pxl_row_in_fullres", "pxl_col_in_fullres", "in_tissue"]]
            )
    else:
        h5_path    = binned_dir / "adata_002um.h5ad"
        adata_full = sc.read_h5ad(str(h5_path), backed='r')
        if pos_df is not None:
            pos_idx = pos_df.set_index("barcode")[["pxl_row_in_fullres", "pxl_col_in_fullres", "in_tissue"]]
            for col in pos_idx.columns:
                if col not in adata_full.obs.columns:
                    adata_full.obs[col] = pos_idx[col]

    records = []
    for roi in rois:
        name = roi["name"]
        x0, y0 = roi["x"], roi["y"]
        w,  h  = roi["width_px"], roi["height_px"]

        # 裁切 ROI bins
        col_vals = adata_full.obs.get("pxl_col_in_fullres")
        row_vals = adata_full.obs.get("pxl_row_in_fullres")
        if col_vals is None or row_vals is None:
            logger.warning(f"跳過 {name}：缺少空間座標欄位")
            continue

        mask_roi = (
            col_vals.between(x0, x0 + w) &
            row_vals.between(y0, y0 + h)
        )
        # 只將當前 ROI 的 bins 載入記憶體（backed 模式下節省 RAM）
        a_view = adata_full[mask_roi]
        a = a_view.to_memory() if adata_full.isbacked else a_view.copy()
        if a.n_obs == 0:
            logger.warning(f"跳過 {name}：ROI 內無 bins")
            continue

        for method in ("nuc", "mcseg"):
            mask_path = qc_dir / f"{name}_{method}.npy"
            if not mask_path.exists():
                logger.warning(f"跳過 {name}_{method}：遮罩不存在")
                continue

            mask = np.load(mask_path)
            n_cells = int(mask.max())
            if n_cells == 0:
                records.append({"roi": name, "method": method,
                                 "n_cells": 0, "ftc": 0.0, "ned": 0.0, "coexp_rate": 0.0})
                continue

            count_mat = _build_count_matrix(mask.copy(), a, x0, y0, dilation_px=6)

            # FTC
            exp_mask = expand_labels(mask, distance=6) if method == "mcseg" else mask
            col_v = (a.obs["pxl_col_in_fullres"].values.astype(int) - x0)
            row_v = (a.obs["pxl_row_in_fullres"].values.astype(int) - y0)
            mh, mw = exp_mask.shape
            in_bounds = (col_v >= 0) & (col_v < mw) & (row_v >= 0) & (row_v < mh)
            in_tissue = a.obs.get("in_tissue", pd.Series(1, index=a.obs.index)).values[in_bounds].astype(bool)
            cell_ids_v = exp_mask[row_v[in_bounds], col_v[in_bounds]]
            ftc = float((cell_ids_v > 0)[in_tissue].mean()) if in_tissue.sum() > 0 else 0.0

            # NED
            dilated = dilation(mask, footprint_rectangle((3, 3)))
            ned_vals = []
            sample_ids = np.random.choice(
                range(1, n_cells + 1), min(ned_sample, n_cells), replace=False
            )
            for cid in sample_ids:
                neighbors = np.unique(dilated[mask == cid])
                neighbors = neighbors[(neighbors > 0) & (neighbors != cid)]
                if len(neighbors) == 0:
                    continue
                p = np.asarray(count_mat[cid - 1].toarray()).flatten().astype(float)
                for nid in neighbors[:3]:
                    if nid - 1 >= count_mat.shape[0]:
                        continue
                    q = np.asarray(count_mat[nid - 1].toarray()).flatten().astype(float)
                    ned_vals.append(_hellinger(p, q))
            ned = float(np.mean(ned_vals)) if ned_vals else 0.0

            # Artificial co-expression rate
            coexp_rates = []
            for ga, gb in impossible_pairs:
                if ga not in a.var_names or gb not in a.var_names:
                    continue
                ia = a.var_names.get_loc(ga)
                ib = a.var_names.get_loc(gb)
                if ia >= count_mat.shape[1] or ib >= count_mat.shape[1]:
                    continue
                ca = np.asarray(count_mat[:, ia].toarray()).flatten() > 0
                cb = np.asarray(count_mat[:, ib].toarray()).flatten() > 0
                coexp_rates.append(float((ca & cb).mean()))
            coexp = float(np.mean(coexp_rates)) if coexp_rates else 0.0

            records.append({
                "roi": name, "method": method,
                "n_cells": n_cells,
                "ftc": round(ftc, 3),
                "ned": round(ned, 3),
                "coexp_rate": round(coexp, 4),
            })
            print(f"  {name} {method}: cells={n_cells} FTC={ftc:.3f} NED={ned:.3f} coexp={coexp:.4f}")
            del count_mat
            gc.collect()

        # ROI 完成後釋放 ROI 子矩陣
        del a
        gc.collect()

    if adata_full.isbacked:
        adata_full.file.close()

    df = pd.DataFrame(records)
    out_csv = Path(out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_csv, index=False)
    print(f"[OK] metrics saved: {out_csv}")
    return df


if __name__ == "__main__":
    compute_metrics()
