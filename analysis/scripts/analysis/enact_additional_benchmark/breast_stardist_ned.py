"""
breast_stardist_ned.py
======================
計算乳癌 StarDist+WBA 的 NED（Neighbour Expression Divergence）。

流程（與 CRC t01/02 + 04_metrics_de 一致）：
  1. 從 cells_df.csv 讀取 StarDist polygon geometry
  2. 光柵化到 crop 區域的局部 instance mask（分塊處理，避免記憶體爆炸）
  3. 對 crop 內所有 in-tissue bins 做 WBA 分派（bin → cell_id，weight=1）
  4. 從 filtered_feature_bc_matrix.h5 建 cell × gene AnnData
  5. 計算 NED（Hellinger distance，1000 HVGs，3000 neighbour pairs）

輸出：
  submission_bioinformatics/results/enact_breast_f1/stardist_ned_breast.txt
"""

from __future__ import annotations
import time
import numpy as np
import pandas as pd
import scipy.sparse as sp
from pathlib import Path
from shapely import from_wkt
from shapely.affinity import translate
from skimage.draw import polygon as sk_polygon

# ── paths ──────────────────────────────────────────────────────────────────────
PLAN_A       = Path("/Volumes/SSD/plan_a")
BSC_DIR      = PLAN_A / "tissue sample/bsc/binned_outputs_311/binned_outputs/square_002um"
H5_PATH      = BSC_DIR / "filtered_feature_bc_matrix.h5"
TP_PATH      = BSC_DIR / "spatial/tissue_positions.parquet"
CELLS_DF_ZIP = PLAN_A / "tissue sample/ENACT_additional_samples/human_breast/output_files/human_breast_cancer.zip"
CELLS_DF_MEMBER = "home/oneai/synthetic_data/human_breast_cancer/cells_df.csv"
OUT_PATH     = PLAN_A / "submission_bioinformatics/results/enact_breast_f1/stardist_ned_breast.txt"

# crop bounds (WSI fullres coords)
CROP_X0, CROP_X1 = 3588, 27772
CROP_Y0, CROP_Y1 = 1641, 23646
CROP_W = CROP_X1 - CROP_X0   # 24184
CROP_H = CROP_Y1 - CROP_Y0   # 22005

# ── helpers ────────────────────────────────────────────────────────────────────

def load_h5_matrix(h5_path: Path):
    import h5py
    with h5py.File(h5_path, "r") as f:
        grp = f["matrix"]
        barcodes   = grp["barcodes"][:].astype(str)
        gene_names = grp["features"]["name"][:].astype(str)
        data    = grp["data"][:]
        indices = grp["indices"][:]
        indptr  = grp["indptr"][:]
        shape   = tuple(grp["shape"][:])
        mat = sp.csc_matrix((data, indices, indptr), shape=shape).T.tocsr()
    return barcodes, gene_names, mat


def rasterize_tile(cells_tile: pd.DataFrame, tile_x0: int, tile_y0: int,
                   tile_w: int, tile_h: int) -> np.ndarray:
    mask = np.zeros((tile_h, tile_w), dtype=np.int32)
    for _, row in cells_tile.iterrows():
        geom   = translate(row["geom"], xoff=-tile_x0, yoff=-tile_y0)
        coords = np.array(geom.exterior.coords)
        cols   = np.clip(coords[:, 0], 0, tile_w - 1)
        rows   = np.clip(coords[:, 1], 0, tile_h - 1)
        rr, cc = sk_polygon(rows, cols, shape=(tile_h, tile_w))
        if len(rr):
            mask[rr, cc] = int(row["local_id"])
    return mask


def compute_ned(mask: np.ndarray, adata) -> float:
    from scipy.ndimage import grey_dilation
    if adata.n_obs < 5 or int(mask.max()) < 2:
        return np.nan
    X = adata.X
    if sp.issparse(X):
        X = np.asarray(X.todense(), dtype=np.float32)
    else:
        X = np.array(X, dtype=np.float32)
    n_hvgs = min(1000, X.shape[1])
    if X.shape[1] > n_hvgs:
        hvg_idx = np.argpartition(X.var(axis=0), -n_hvgs)[-n_hvgs:]
        X = X[:, hvg_idx]
    row_sums = np.maximum(X.sum(axis=1, keepdims=True), 1e-10)
    X_prob   = (X / row_sums).astype(np.float32)
    obs_ids  = adata.obs["cell_id"].values
    cid_to_row = {int(c): i for i, c in enumerate(obs_ids)}
    struct  = np.ones((3, 3), dtype=np.int32)
    dilated = grey_dilation(mask.astype(np.int32), footprint=struct)
    bnd     = (mask > 0) & (dilated != mask)
    ci = mask[bnd].astype(np.int32)
    cj = dilated[bnd].astype(np.int32)
    valid = (cj > 0) & (cj != ci)
    ci, cj = ci[valid], cj[valid]
    if len(ci) == 0:
        return np.nan
    pairs = np.unique(np.sort(np.stack([ci, cj], axis=1), axis=1), axis=0)
    known = set(cid_to_row.keys())
    ok = np.array([(int(a) in known and int(b) in known) for a, b in pairs])
    pairs = pairs[ok]
    if len(pairs) < 5:
        return np.nan
    if len(pairs) > 3000:
        rng   = np.random.default_rng(42)
        pairs = pairs[rng.choice(len(pairs), 3000, replace=False)]
    i_idx = np.array([cid_to_row[int(a)] for a in pairs[:, 0]])
    j_idx = np.array([cid_to_row[int(b)] for b in pairs[:, 1]])
    Xi, Xj = X_prob[i_idx], X_prob[j_idx]
    hell = np.sqrt(
        np.sum((np.sqrt(np.maximum(Xi, 0)) - np.sqrt(np.maximum(Xj, 0)))**2, axis=1) / 2
    )
    return float(np.clip(np.mean(hell), 0, 1))


# ── main ───────────────────────────────────────────────────────────────────────

def main():
    t0 = time.time()
    print("=== Breast StarDist NED ===\n")

    # [1] load cells_df geometry
    print("[1] Loading cells_df geometry...")
    import zipfile, io
    z     = zipfile.ZipFile(CELLS_DF_ZIP)
    cells = pd.read_csv(io.BytesIO(z.read(CELLS_DF_MEMBER)),
                        usecols=["id", "geometry", "cell_x", "cell_y"])
    cells = cells[
        (cells["cell_x"] >= 0) & (cells["cell_x"] < CROP_W) &
        (cells["cell_y"] >= 0) & (cells["cell_y"] < CROP_H) &
        cells["geometry"].notna()
    ].copy().reset_index(drop=True)
    cells["local_id"] = np.arange(1, len(cells) + 1, dtype=np.int32)
    print(f"  {len(cells):,} cells in crop region")

    print("  Parsing Shapely geometries...")
    t1 = time.time()
    cells["geom"] = from_wkt(cells["geometry"].values)
    print(f"  Done in {time.time()-t1:.1f}s")

    # [2] rasterize tile-by-tile
    print("\n[2] Rasterising StarDist polygons to instance mask...")
    TILE      = 2000
    mask_full = np.zeros((CROP_H, CROP_W), dtype=np.int32)
    xt = list(range(0, CROP_W, TILE)) + [CROP_W]
    yt = list(range(0, CROP_H, TILE)) + [CROP_H]
    n_tiles = (len(xt)-1) * (len(yt)-1)
    done = 0
    for yi in range(len(yt)-1):
        for xi in range(len(xt)-1):
            tx0, tx1 = xt[xi], xt[xi+1]
            ty0, ty1 = yt[yi], yt[yi+1]
            tw, th   = tx1 - tx0, ty1 - ty0
            buf = 60
            tile_cells = cells[
                (cells["cell_x"] >= tx0 - buf) & (cells["cell_x"] < tx1 + buf) &
                (cells["cell_y"] >= ty0 - buf) & (cells["cell_y"] < ty1 + buf)
            ]
            if len(tile_cells):
                tile_mask = rasterize_tile(tile_cells, tx0, ty0, tw, th)
                nnz = tile_mask > 0
                mask_full[ty0:ty1, tx0:tx1][nnz] = tile_mask[nnz]
            done += 1
            if done % 20 == 0 or done == n_tiles:
                print(f"  tile {done}/{n_tiles}  elapsed {time.time()-t0:.0f}s", flush=True)

    coverage = float((mask_full > 0).mean())
    print(f"  Mask coverage: {coverage:.3f}")

    # [3] load H5 + tissue positions
    print("\n[3] Loading H5 matrix and tissue positions...")
    barcodes, gene_names, mat = load_h5_matrix(H5_PATH)
    bc_to_idx = {b: i for i, b in enumerate(barcodes)}

    tp      = pd.read_parquet(TP_PATH)
    tp_crop = tp[
        (tp["in_tissue"] == 1) &
        (tp["pxl_col_in_fullres"] >= CROP_X0) & (tp["pxl_col_in_fullres"] <= CROP_X1) &
        (tp["pxl_row_in_fullres"] >= CROP_Y0) & (tp["pxl_row_in_fullres"] <= CROP_Y1)
    ].copy()
    tp_crop["local_col"] = (tp_crop["pxl_col_in_fullres"] - CROP_X0).astype(int)
    tp_crop["local_row"] = (tp_crop["pxl_row_in_fullres"] - CROP_Y0).astype(int)
    tp_crop = tp_crop[
        (tp_crop["local_col"] >= 0) & (tp_crop["local_col"] < CROP_W) &
        (tp_crop["local_row"] >= 0) & (tp_crop["local_row"] < CROP_H)
    ]
    print(f"  {len(tp_crop):,} in-tissue crop bins")

    # [4] attribution: bin → StarDist cell via mask lookup
    print("\n[4] Attribution via rasterised mask...")
    locs = tp_crop[["local_row", "local_col", "barcode"]].copy()
    locs["cell_id"] = mask_full[
        locs["local_row"].values, locs["local_col"].values
    ]
    wba = locs[locs["cell_id"] > 0].copy()
    wba["weight"] = 1.0
    print(f"  Attributed: {len(wba):,} / {len(tp_crop):,} ({len(wba)/len(tp_crop):.1%})")
    print(f"  Unique cells: {wba['cell_id'].nunique():,}")

    # [5] build cell × gene AnnData
    print("\n[5] Building cell × gene AnnData...")
    import anndata as ad
    wba_valid  = wba[wba["barcode"].isin(bc_to_idx)].copy()
    cell_ids   = np.sort(wba_valid["cell_id"].unique().astype(int))
    cid_to_row = {c: i for i, c in enumerate(cell_ids)}
    n_cells    = len(cell_ids)
    print(f"  {n_cells:,} cells × {mat.shape[1]:,} genes")

    wba_valid["bin_row"]  = wba_valid["barcode"].map(bc_to_idx).astype(int)
    wba_valid["cell_row"] = wba_valid["cell_id"].astype(int).map(cid_to_row)
    wba_valid = wba_valid.dropna(subset=["cell_row"])
    wba_valid["cell_row"] = wba_valid["cell_row"].astype(int)

    unique_bin_rows  = wba_valid["bin_row"].unique()
    print(f"  Fetching {len(unique_bin_rows):,} unique bin rows from H5...")
    bin_mat          = sp.csr_matrix(mat[unique_bin_rows, :], dtype=np.float32)
    bin_row_to_local = {r: i for i, r in enumerate(unique_bin_rows)}
    wba_valid["local_bin"] = wba_valid["bin_row"].map(bin_row_to_local).astype(int)

    local_bins = wba_valid["local_bin"].values.astype(int)
    cell_rows  = wba_valid["cell_row"].values.astype(int)
    weights    = wba_valid["weight"].values.astype(np.float32)

    W_diag        = sp.diags(weights, format="csr", dtype=np.float32)
    weighted_bins = W_diag.dot(bin_mat[local_bins, :])
    agg_mat       = sp.csr_matrix(
        (np.ones(len(cell_rows), dtype=np.float32),
         (cell_rows, np.arange(len(cell_rows)))),
        shape=(n_cells, len(cell_rows))
    ).dot(weighted_bins)

    obs   = pd.DataFrame({"cell_id": cell_ids}, index=[str(c) for c in cell_ids])
    adata = ad.AnnData(X=agg_mat, obs=obs, var=pd.DataFrame(index=gene_names))
    print(f"  AnnData: {adata.n_obs:,} cells × {adata.n_vars:,} genes")

    # [6] Median UMI/cell
    print("\n[6] Computing median UMI/cell...")
    import scipy.sparse as sp_
    X = adata.X
    cell_umis = np.asarray(X.sum(axis=1) if sp_.issparse(X) else X.sum(axis=1)).ravel()
    umi_median = float(np.median(cell_umis))
    print(f"  Median UMI/cell StarDist+WBA: {umi_median:.1f}")

    # [7] NED
    print("\n[7] Computing NED (Hellinger, 3000 pairs, 1000 HVGs)...")
    ned = compute_ned(mask_full, adata)
    print(f"  NED StarDist+WBA (breast): {ned:.4f}")

    result = (
        f"NED StarDist+WBA (breast): {ned:.4f}\n"
        f"median_umi_per_cell: {umi_median:.1f}\n"
        f"n_cells: {adata.n_obs:,}\n"
        f"mask_coverage: {coverage:.4f}\n"
        f"elapsed: {time.time()-t0:.1f}s\n"
    )
    OUT_PATH.write_text(result)
    print(f"\n✓ Saved: {OUT_PATH}")
    print(f"  Total elapsed: {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
