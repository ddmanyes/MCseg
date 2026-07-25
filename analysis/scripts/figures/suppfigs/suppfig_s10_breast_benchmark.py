"""
suppfig_s10_breast_benchmark.py
================================
Supp. Fig. S10 — MCseg cross-cancer generalisation: Breast Cancer (ENACT)

Layout:
  a  overview H&E + dense inset + sparse inset (MCseg overlay + StarDist dots)
  b  FTC: MCseg vs StarDist+WBA
  c  Median UMI/cell: MCseg vs StarDist+WBA
  d  NED: MCseg vs StarDist+WBA (StarDist NED from polygon rasterisation)

Outputs:
  manuscript/supplementary/SuppFigS10.png
  submission_bioinformatics/supplementary/SuppFigS10.png
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scipy.sparse as sp
from pathlib import Path
import tifffile

# ── paths ──────────────────────────────────────────────────────────────────────

PLAN_A      = Path("/Volumes/SSD/plan_a")
RESULTS_DIR = PLAN_A / "submission_bioinformatics" / "results" / "enact_breast_f1"
BSC_DIR     = PLAN_A / "tissue sample" / "bsc" / "binned_outputs_311" / "binned_outputs" / "square_002um"
H5_PATH     = BSC_DIR / "filtered_feature_bc_matrix.h5"
TP_PATH     = BSC_DIR / "spatial" / "tissue_positions.parquet"
HE_PATH     = RESULTS_DIR / "he_crop_breast.tif"
MASK_PATH   = RESULTS_DIR / "mcseg_mask_breast.npy"
WBA_PATH    = RESULTS_DIR / "wba_attribution_breast.csv"
STARDIST_CSV = PLAN_A / "tissue sample" / "ENACT_additional_samples" / "human_breast" / "merged_results.csv"
SD_NED_PATH  = RESULTS_DIR / "stardist_ned_breast.txt"

OUT_DIR     = PLAN_A / "manuscript" / "supplementary"
OUT_DIR_SUB = PLAN_A / "submission_bioinformatics" / "supplementary"
OUT_DIR.mkdir(parents=True, exist_ok=True)
OUT_DIR_SUB.mkdir(parents=True, exist_ok=True)

# crop region coords (same as enact_breast_f1_wba.py)
CROP_X0, CROP_X1 = 3588, 27772
CROP_Y0, CROP_Y1 = 1641, 23646

# inset crops (800×800 each)
INSET_SZ = 800
# dense inset — highest cell density (density=0.671)
DENSE_Y0, DENSE_X0 = 17076, 7233
# sparse inset — moderate density, far from dense region
SPARSE_Y0, SPARSE_X0 = 8458, 10740

# overview crop — covers both insets (12000×9500)
OV_R0, OV_R1 = 6500, 18500
OV_C0, OV_C1 = 4500, 14000
# inset positions within overview
DENSE_IN_OV_R  = DENSE_Y0  - OV_R0   # 10576
DENSE_IN_OV_C  = DENSE_X0  - OV_C0   # 2733
SPARSE_IN_OV_R = SPARSE_Y0 - OV_R0   # 1958
SPARSE_IN_OV_C = SPARSE_X0 - OV_C0   # 6240

# ── style ──────────────────────────────────────────────────────────────────────

plt.rcParams.update({
    "font.family":       ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size":         7,
    "savefig.dpi":       300,
    "savefig.facecolor": "white",
    "axes.linewidth":    0.6,
})
MM = 1 / 25.4

COL_MCSEG    = "#4292C6"
COL_STARDIST = "#DCB400"

COLOR_MCSEG_FILL = (70, 130, 220)   # blue — matches fig3a MCseg v2 colour
COLOR_SD_DOT     = (60, 200, 80)    # green — matches fig3a style dots

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


def build_mcseg_adata(wba_crop: pd.DataFrame, gene_names, mat,
                      bc_to_idx: dict):
    """Build cell × gene AnnData for MCseg from WBA crop subset (vectorised)."""
    import anndata as ad

    wba_valid = wba_crop[wba_crop["cell_id"] > 0].copy()
    wba_valid = wba_valid[wba_valid["barcode"].isin(bc_to_idx)]

    cell_ids   = np.sort(wba_valid["cell_id"].unique().astype(int))
    cid_to_row = {c: i for i, c in enumerate(cell_ids)}
    n_cells    = len(cell_ids)
    n_genes    = mat.shape[1]
    print(f"  Building {n_cells:,} × {n_genes:,} matrix (vectorised)...")

    wba_valid = wba_valid.copy()
    wba_valid["bin_row"]  = wba_valid["barcode"].map(bc_to_idx)
    wba_valid["cell_row"] = wba_valid["cell_id"].astype(int).map(cid_to_row)
    wba_valid = wba_valid.dropna(subset=["bin_row", "cell_row"])
    wba_valid["bin_row"]  = wba_valid["bin_row"].astype(int)
    wba_valid["cell_row"] = wba_valid["cell_row"].astype(int)

    # fetch all unique bin rows at once
    unique_bin_rows  = wba_valid["bin_row"].unique()
    print(f"  Fetching {len(unique_bin_rows):,} unique bin rows from H5...")
    bin_mat          = sp.csr_matrix(mat[unique_bin_rows, :], dtype=np.float32)
    bin_row_to_local = {r: i for i, r in enumerate(unique_bin_rows)}
    wba_valid        = wba_valid.copy()
    wba_valid["local_bin"] = wba_valid["bin_row"].map(bin_row_to_local)

    local_bins = wba_valid["local_bin"].values.astype(int)
    cell_rows  = wba_valid["cell_row"].values.astype(int)
    weights    = wba_valid["weight"].values.astype(np.float32)

    # weighted_bins[i] = weight[i] * bin_mat[local_bins[i], :]
    W_diag       = sp.diags(weights, format="csr", dtype=np.float32)
    weighted_bins = W_diag.dot(bin_mat[local_bins, :])   # n_wba × n_genes

    # aggregate rows by cell_row
    agg_weights = sp.csr_matrix(
        (np.ones(len(cell_rows), dtype=np.float32),
         (cell_rows, np.arange(len(cell_rows)))),
        shape=(n_cells, len(cell_rows))
    )
    agg = agg_weights.dot(weighted_bins)   # n_cells × n_genes

    obs = pd.DataFrame({"cell_id": cell_ids}, index=[str(c) for c in cell_ids])
    obs["n_umis"]  = np.asarray(agg.sum(axis=1)).ravel()
    obs["n_genes"] = np.asarray((agg > 0).sum(axis=1)).ravel()

    return ad.AnnData(X=agg, obs=obs, var=pd.DataFrame(index=gene_names))


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
        gene_var = X.var(axis=0)
        hvg_idx  = np.argpartition(gene_var, -n_hvgs)[-n_hvgs:]
        X = X[:, hvg_idx]

    row_sums = np.maximum(X.sum(axis=1, keepdims=True), 1e-10)
    X_prob   = (X / row_sums).astype(np.float32)

    obs_cell_ids = adata.obs["cell_id"].values
    cid_to_row   = {int(c): i for i, c in enumerate(obs_cell_ids)}

    struct  = np.ones((3, 3), dtype=np.int32)
    dilated = grey_dilation(mask.astype(np.int32), footprint=struct)
    bnd     = (mask > 0) & (dilated != mask)
    ci      = mask[bnd].astype(np.int32)
    cj      = dilated[bnd].astype(np.int32)
    valid   = (cj > 0) & (cj != ci)
    ci, cj  = ci[valid], cj[valid]

    if len(ci) == 0:
        return np.nan

    pairs_arr = np.unique(np.sort(np.stack([ci, cj], axis=1), axis=1), axis=0)
    known     = set(cid_to_row.keys())
    ma = np.array([int(a) in known for a in pairs_arr[:, 0]])
    mb = np.array([int(b) in known for b in pairs_arr[:, 1]])
    pairs_arr = pairs_arr[ma & mb]

    if len(pairs_arr) < 5:
        return np.nan
    if len(pairs_arr) > 3000:
        rng = np.random.default_rng(42)
        idx = rng.choice(len(pairs_arr), 3000, replace=False)
        pairs_arr = pairs_arr[idx]

    i_idx = np.array([cid_to_row[int(a)] for a in pairs_arr[:, 0]])
    j_idx = np.array([cid_to_row[int(b)] for b in pairs_arr[:, 1]])
    Xi = X_prob[i_idx]
    Xj = X_prob[j_idx]
    sqrt_i = np.sqrt(np.maximum(Xi, 0))
    sqrt_j = np.sqrt(np.maximum(Xj, 0))
    hell = np.sqrt(np.sum((sqrt_i - sqrt_j) ** 2, axis=1) / 2)
    return float(np.clip(np.mean(hell), 0, 1))


def mask_to_rgba_overlay(mask: np.ndarray, color_rgb: tuple,
                         alpha: float = 0.25,
                         boundary_alpha: float = 0.88) -> np.ndarray:
    import cv2
    h, w = mask.shape
    overlay = np.zeros((h, w, 4), dtype=np.float32)
    fill = mask > 0
    overlay[fill, :3] = np.array(color_rgb) / 255.0
    overlay[fill, 3]  = alpha
    m32    = mask.astype(np.int32)
    binary = fill.astype(np.uint8)
    eroded = cv2.erode(binary, np.ones((3, 3), np.uint8), iterations=1)
    outer  = (binary - eroded).astype(bool)
    r_diff = np.zeros((h, w), bool); r_diff[:, :-1] = m32[:, 1:] != m32[:, :-1]
    d_diff = np.zeros((h, w), bool); d_diff[:-1, :]  = m32[1:, :] != m32[:-1, :]
    boundary = outer | r_diff | d_diff
    overlay[boundary, :3] = 1.0
    overlay[boundary, 3]  = boundary_alpha
    return overlay


def blend(base: np.ndarray, over: np.ndarray) -> np.ndarray:
    bg = base.astype(np.float32) / 255.0
    fg, a = over[..., :3], over[..., 3:4]
    return np.clip((fg * a + bg * (1.0 - a)) * 255, 0, 255).astype(np.uint8)


# ── main ───────────────────────────────────────────────────────────────────────

def main():
    print("=== SuppFig S10: Breast Cancer Benchmark ===")

    # ── load shared resources ─────────────────────────────────────────────────
    print("\n[2] Loading H5 and tissue positions...")
    barcodes, gene_names, mat = load_h5_matrix(H5_PATH)
    bc_to_idx = {b: i for i, b in enumerate(barcodes)}

    tp = pd.read_parquet(TP_PATH)
    tp_crop = tp[
        (tp["in_tissue"] == 1) &
        (tp["pxl_col_in_fullres"] >= CROP_X0) & (tp["pxl_col_in_fullres"] <= CROP_X1) &
        (tp["pxl_row_in_fullres"] >= CROP_Y0) & (tp["pxl_row_in_fullres"] <= CROP_Y1)
    ]
    n_crop_bins = len(tp_crop)
    crop_bcs    = set(tp_crop["barcode"].values)
    print(f"  Crop in-tissue bins: {n_crop_bins:,}")

    wba     = pd.read_csv(WBA_PATH)
    wba_crop = wba[wba["barcode"].isin(crop_bcs)]

    # ── FTC ───────────────────────────────────────────────────────────────────
    print("\n[3] FTC...")
    ftc_mcseg = float((wba_crop["cell_id"] > 0).sum() / n_crop_bins)
    sd        = pd.read_csv(STARDIST_CSV)
    # merged_results cell_x/cell_y are crop-relative coords; filter to crop bounds
    # StarDist WBA: num_shared_bins are counted multiple times across adjacent cells
    # Use num_unique_bins only to avoid double-counting (equivalent to MCseg's FTC)
    ftc_sd    = float(sd["num_unique_bins"].sum() / n_crop_bins)
    print(f"  FTC — MCseg: {ftc_mcseg:.3f}  StarDist: {ftc_sd:.3f}")

    # ── median UMI/cell ───────────────────────────────────────────────────────
    print("\n[4] Median UMI/cell...")
    # StarDist: load pre-computed value from H5-based WBA (same pipeline as MCseg)
    sd_ned_lines = {l.split(": ")[0]: l.split(": ")[1]
                    for l in SD_NED_PATH.read_text().splitlines() if ": " in l}
    sd_umi_median = float(sd_ned_lines["median_umi_per_cell"])

    wba_valid = wba_crop[wba_crop["cell_id"] > 0].copy()
    wba_valid = wba_valid[wba_valid["barcode"].isin(bc_to_idx)]
    unique_bcs   = wba_valid["barcode"].unique()
    bc_indices   = np.array([bc_to_idx[b] for b in unique_bcs])
    bin_umi_vals = np.asarray(mat[bc_indices, :].sum(axis=1)).ravel()
    bc_umi_map   = dict(zip(unique_bcs, bin_umi_vals))
    wba_valid    = wba_valid.copy()
    wba_valid["bin_umi"] = wba_valid["barcode"].map(bc_umi_map).fillna(0)
    cell_umi     = (wba_valid.groupby("cell_id")
                   .apply(lambda g: (g["bin_umi"] * g["weight"]).sum(), include_groups=False)
                   .reset_index(name="total_umi"))
    mcseg_umi_median = float(cell_umi["total_umi"].median())
    print(f"  Median UMI — MCseg: {mcseg_umi_median:.1f}  StarDist: {sd_umi_median:.1f}")

    # ── NED ───────────────────────────────────────────────────────────────────
    print("\n[5] Building AnnData for NED (MCseg)...")
    mask_full = np.load(MASK_PATH)
    adata     = build_mcseg_adata(wba_crop, gene_names, mat, bc_to_idx)
    print(f"  AnnData: {adata.n_obs:,} cells × {adata.n_vars:,} genes")
    print("  Computing NED (Hellinger, 3000 neighbour pairs, 1000 HVGs)...")
    ned_mcseg = compute_ned(mask_full, adata)
    print(f"  NED MCseg: {ned_mcseg:.3f}")

    # load pre-computed StarDist NED (polygon rasterisation, same dict as UMI)
    ned_sd = float(sd_ned_lines["NED StarDist+WBA (breast)"])
    print(f"  NED StarDist+WBA (cached): {ned_sd:.4f}")

    # ── overlay images ────────────────────────────────────────────────────────
    print("\n[6] Overlay images...")
    he = tifffile.imread(str(HE_PATH))

    # overview: pure H&E (no overlay — too large to render fast), with boxes
    he_ov   = he[OV_R0:OV_R1, OV_C0:OV_C1]

    # dense inset: 800×800, MCseg overlay + StarDist dots
    he_dense   = he[DENSE_Y0:DENSE_Y0+INSET_SZ, DENSE_X0:DENSE_X0+INSET_SZ]
    mask_dense = mask_full[DENSE_Y0:DENSE_Y0+INSET_SZ, DENSE_X0:DENSE_X0+INSET_SZ]
    img_dense  = blend(he_dense, mask_to_rgba_overlay(mask_dense, COLOR_MCSEG_FILL))

    # sparse inset: 800×800, MCseg overlay + StarDist dots
    he_sparse   = he[SPARSE_Y0:SPARSE_Y0+INSET_SZ, SPARSE_X0:SPARSE_X0+INSET_SZ]
    mask_sparse = mask_full[SPARSE_Y0:SPARSE_Y0+INSET_SZ, SPARSE_X0:SPARSE_X0+INSET_SZ]
    img_sparse  = blend(he_sparse, mask_to_rgba_overlay(mask_sparse, COLOR_MCSEG_FILL))

    # StarDist cells in each inset (crop-relative: cell_x=col, cell_y=row)
    sd_all = pd.read_csv(STARDIST_CSV)

    def sd_in_window(y0, x0, sz):
        w = sd_all[
            (sd_all["cell_y"] >= y0) & (sd_all["cell_y"] < y0 + sz) &
            (sd_all["cell_x"] >= x0) & (sd_all["cell_x"] < x0 + sz)
        ].copy()
        w["local_x"] = w["cell_x"] - x0
        w["local_y"] = w["cell_y"] - y0
        return w

    sd_dense  = sd_in_window(DENSE_Y0,  DENSE_X0,  INSET_SZ)
    sd_sparse = sd_in_window(SPARSE_Y0, SPARSE_X0, INSET_SZ)
    print(f"  StarDist cells — dense: {len(sd_dense)}, sparse: {len(sd_sparse)}")

    # ── figure ────────────────────────────────────────────────────────────────
    print("\n[7] Building figure...")
    import matplotlib.lines as mlines
    import matplotlib.patches as mpatches

    OV_H = OV_R1 - OV_R0   # 12000
    OV_W = OV_C1 - OV_C0   # 9500
    ov_aspect = OV_W / OV_H  # ~0.79

    fig_w = 183 * MM
    fig_h = 130 * MM
    fig = plt.figure(figsize=(fig_w, fig_h))

    # 2-row layout: top=images, bottom=bars
    gs = fig.add_gridspec(
        2, 1,
        height_ratios=[1.6, 1.0],
        hspace=0.38,
        left=0.06, right=0.97,
        top=0.93, bottom=0.08,
    )
    # top row: overview + dense inset + sparse inset
    gs_top = gs[0].subgridspec(1, 3, wspace=0.06,
                                width_ratios=[ov_aspect, 1.0, 1.0])
    ax_ov     = fig.add_subplot(gs_top[0])
    ax_dense  = fig.add_subplot(gs_top[1])
    ax_sparse = fig.add_subplot(gs_top[2])

    # bottom row: 3 bar panels centred
    gs_bot = gs[1].subgridspec(1, 3, wspace=0.55)
    ax_b = fig.add_subplot(gs_bot[0])
    ax_c = fig.add_subplot(gs_bot[1])
    ax_d = fig.add_subplot(gs_bot[2])

    # ── overview ──────────────────────────────────────────────────────────────
    ax_ov.imshow(he_ov, origin="upper", interpolation="antialiased")
    BOX_COLORS = {"dense": "#FFD700", "sparse": "#00CFFF"}
    for label, ry, rx in [("dense",  DENSE_IN_OV_R,  DENSE_IN_OV_C),
                           ("sparse", SPARSE_IN_OV_R, SPARSE_IN_OV_C)]:
        rect = mpatches.Rectangle(
            (rx, ry), INSET_SZ, INSET_SZ,
            linewidth=1.5, edgecolor=BOX_COLORS[label], facecolor="none"
        )
        ax_ov.add_patch(rect)
        ax_ov.text(rx + INSET_SZ / 2, ry - OV_H * 0.01,
                   label, color=BOX_COLORS[label],
                   ha="center", va="bottom", fontsize=5.5, fontweight="bold")
    ax_ov.set_xticks([]); ax_ov.set_yticks([])
    ax_ov.set_title("Breast Cancer (no re-tuning)", fontsize=6.5, fontweight="bold", pad=3)
    ax_ov.text(-0.06, 1.02, "a", transform=ax_ov.transAxes,
               fontsize=11, fontweight="bold", va="bottom", ha="right")

    def draw_inset(ax, img, sd_win, box_color, label):
        ax.imshow(img, origin="upper", interpolation="antialiased")
        ax.scatter(sd_win["local_x"], sd_win["local_y"],
                   s=3, color=np.array(COLOR_SD_DOT)/255,
                   linewidths=0, alpha=0.85)
        ax.set_xticks([]); ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_edgecolor(box_color); spine.set_linewidth(1.5)
        ax.set_title(label, fontsize=6, color=box_color, fontweight="bold", pad=3)
        # scale bar: 100 px ≈ 50 µm
        bar_px = 100
        ax.plot([INSET_SZ - bar_px - 15, INSET_SZ - 15],
                [INSET_SZ - 22, INSET_SZ - 22],
                color="white", linewidth=1.5, solid_capstyle="butt")
        ax.text(INSET_SZ - 15 - bar_px / 2, INSET_SZ - 30,
                "50 µm", color="white", ha="center", va="bottom", fontsize=5)

    draw_inset(ax_dense,  img_dense,  sd_dense,  BOX_COLORS["dense"],  "dense")
    draw_inset(ax_sparse, img_sparse, sd_sparse, BOX_COLORS["sparse"], "sparse")

    # shared legend on sparse inset
    h_mc = mpatches.Patch(color=np.array(COLOR_MCSEG_FILL)/255, alpha=0.7, label="MCseg")
    h_sd = mlines.Line2D([], [], color=np.array(COLOR_SD_DOT)/255, marker="o",
                         linestyle="None", markersize=3, label="StarDist")
    ax_sparse.legend(handles=[h_mc, h_sd], fontsize=5, loc="upper left",
                     framealpha=0.80, edgecolor="#cccccc", handlelength=1.0)

    def bar2(ax, labels, vals, colors, ylabel, title, panel_letter,
             ylim=None, fmt=".3f"):
        bars = ax.bar(labels, vals, color=colors, width=0.55, edgecolor="none")
        for bar, val in zip(bars, vals):
            ax.text(bar.get_x() + bar.get_width()/2,
                    val + (ylim[1] * 0.02 if ylim else max(vals) * 0.03),
                    f"{val:{fmt}}", ha="center", va="bottom",
                    fontsize=6, fontweight="bold")
        if ylim:
            ax.set_ylim(*ylim)
        ax.set_ylabel(ylabel, fontsize=6.5)
        ax.set_title(title, fontsize=6, pad=3)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.tick_params(axis="both", labelsize=6)
        ax.text(-0.22, 1.08, panel_letter, transform=ax.transAxes,
                fontsize=11, fontweight="bold", va="top", ha="right")

    bar2(ax_b,
         ["StarDist\n+WBA", "MCseg\n+WBA"],
         [ftc_sd, ftc_mcseg],
         [COL_STARDIST, COL_MCSEG],
         "FTC ↑", "Transcript\ncapture rate",
         "b", ylim=(0, 1.05))

    bar2(ax_c,
         ["StarDist\n+WBA", "MCseg"],
         [sd_umi_median, mcseg_umi_median],
         [COL_STARDIST, COL_MCSEG],
         "Median UMI/cell ↑", "Median UMI\nper cell",
         "c", ylim=(0, 1700), fmt=".0f")

    bar2(ax_d,
         ["StarDist\n+WBA", "MCseg\n+WBA"],
         [ned_sd, ned_mcseg],
         [COL_STARDIST, COL_MCSEG],
         "NED ↑", "Neighbour\nExpression Div.",
         "d", ylim=(0, 0.75))

    # ── save ──────────────────────────────────────────────────────────────────
    for out_dir in [OUT_DIR, OUT_DIR_SUB]:
        out = out_dir / "SuppFigS10.png"
        fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
        print(f"✓ Saved: {out}")
    plt.close(fig)

    print("\n=== Summary ===")
    print(f"MCseg    — FTC: {ftc_mcseg:.3f}  UMI/cell: {mcseg_umi_median:.1f}  NED: {ned_mcseg:.3f}")
    print(f"StarDist — FTC: {ftc_sd:.3f}  UMI/cell: {sd_umi_median:.1f}  NED: {ned_sd:.4f}")


if __name__ == "__main__":
    main()
