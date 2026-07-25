"""
suppfig_s4_roi10_at2_comparison.py
====================================
SuppFig S4 — ROI10 (Normal Alveolar Region)
Side-by-side Leiden-based AT2 cell-type map:
  Left  (a): 2Cseg  (roi10_cellpose_dilate)
  Right (b): MCseg  (roi10_v12)

Same Leiden pipeline as fig2d_leiden_celltype_map.py.

Output:
  submission_bioinformatics/supplementary/SuppFigS4.png

Run:
  cd /Volumes/SSD/plan_a
  uv run python submission_bioinformatics/scripts/figures/suppfigs/suppfig_s4_roi10_at2_comparison.py
"""

from __future__ import annotations
import warnings
warnings.filterwarnings("ignore")
import os
os.environ["OMP_NUM_THREADS"] = "1"

import numpy as np
import scipy.sparse as sp
import anndata as ad
import tifffile
import zarr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path
from skimage.segmentation import find_boundaries
from sklearn.decomposition import TruncatedSVD
from sklearn.neighbors import NearestNeighbors
import leidenalg
import igraph

# ── Paths ─────────────────────────────────────────────────────────────────────
PROJECT_ROOT = Path("/Volumes/SSD/plan_a/xenium_he_seg")
BTF_PATH = (Path("/Volumes/SSD/plan_a/tissue sample/LUAD/visium") /
            "Visium_HD_Human_Lung_Cancer_post_Xenium_Prime_5K_Experiment2_tissue_image.btf")
MASK_DIR = PROJECT_ROOT / "results" / "masks"
H5AD_DIR = PROJECT_ROOT / "results" / "visiumhd" / "visiumhd_cells"
OUT_PATH = Path("/Volumes/SSD/plan_a/submission_bioinformatics/supplementary/SuppFigS4.png")

ROI10 = dict(x=7562, y=19440, w=3194, h=1587)

PALETTE = {
    "AT2 Pneumocyte": ((0.13, 0.40, 0.67), 0.72),
    "Unresolved":     ((0.73, 0.73, 0.73), 0.35),
}
PALETTE_HEX = {
    "AT2 Pneumocyte": "#2166AC",
    "Unresolved":     "#BBBBBB",
}

plt.rcParams.update({
    "font.family": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 8, "axes.titlesize": 9, "pdf.fonttype": 42,
    "axes.linewidth": 0.8,
})
MM = 1 / 25.4


# ── Helpers ───────────────────────────────────────────────────────────────────

def load_he_crop(roi: dict) -> np.ndarray:
    x, y, w, h = roi["x"], roi["y"], roi["w"], roi["h"]
    with tifffile.TiffFile(str(BTF_PATH)) as tif:
        store = tif.aszarr()
        z = zarr.open(store, mode="r")
        arr = z if not isinstance(z, zarr.Group) else z[0]
        crop = np.array(arr[y:y + h, x:x + w])
    if crop.ndim == 3 and crop.shape[2] == 4:
        crop = crop[:, :, :3]
    return crop


def normalize_log1p(X_sparse) -> np.ndarray:
    X = sp.csr_matrix(X_sparse)
    totals = np.array(X.sum(axis=1)).flatten()
    totals[totals == 0] = 1
    X = X.multiply(1e4 / totals[:, None]).tocsr()
    X.data = np.log1p(X.data)
    return X.toarray()


def leiden_cluster(X_log: np.ndarray, n_top_hvg: int = 2000,
                   n_pcs: int = 50, k: int = 15,
                   resolution: float = 0.4, seed: int = 42) -> np.ndarray:
    var_all = (X_log ** 2).mean(0) - X_log.mean(0) ** 2
    top_idx = np.argsort(var_all)[-n_top_hvg:]
    X_hvg = X_log[:, top_idx]
    actual_pcs = min(n_pcs, X_hvg.shape[1] - 1, X_hvg.shape[0] - 1)
    svd = TruncatedSVD(n_components=actual_pcs, random_state=seed)
    X_pca = svd.fit_transform(X_hvg)
    nn = NearestNeighbors(n_neighbors=k + 1, algorithm="auto", n_jobs=1)
    nn.fit(X_pca)
    _, idxs = nn.kneighbors(X_pca)
    n = X_pca.shape[0]
    edges = [(int(i), int(j)) for i in range(n) for j in idxs[i, 1:]]
    g = igraph.Graph(n=n, edges=edges, directed=False)
    g.simplify()
    part = leidenalg.find_partition(
        g, leidenalg.RBConfigurationVertexPartition,
        resolution_parameter=resolution, seed=seed)
    return np.array(part.membership)


def identify_at2_cluster(leiden_labels: np.ndarray, X_log: np.ndarray,
                          var_names: np.ndarray) -> int:
    clusters = np.unique(leiden_labels)
    at2_genes = ["SFTPC", "SFTPB", "SFTPA1", "SFTPA2"]
    gene_idx = {g: i for i, g in enumerate(var_names)}
    present = [g for g in at2_genes if g in gene_idx]
    if not present:
        raise RuntimeError("No AT2 marker genes found.")
    best_cl, best_score = -1, -1.0
    for cl in clusters:
        mask = leiden_labels == cl
        score = max(float(X_log[mask, gene_idx[g]].mean()) for g in present)
        if score > best_score:
            best_score = score
            best_cl = cl
    return int(best_cl)


def build_composite(adata, mask_npy: np.ndarray, he: np.ndarray,
                    leiden_labels: np.ndarray, at2_cl: int) -> tuple[np.ndarray, int, int, float]:
    cell_annotation = np.where(leiden_labels == at2_cl, "AT2 Pneumocyte", "Unresolved")
    cell_ids = adata.obs["cell_id"].values.astype(int)
    id_max = int(mask_npy.max())
    lut_rgb   = np.zeros((id_max + 1, 3), dtype=np.float32)
    lut_alpha = np.zeros(id_max + 1,      dtype=np.float32)
    for i, cid in enumerate(cell_ids):
        if 0 < cid <= id_max:
            rgb, alpha     = PALETTE[cell_annotation[i]]
            lut_rgb[cid]   = rgb
            lut_alpha[cid] = alpha
    he_f     = he.astype(np.float32) / 255.0
    cell_rgb = lut_rgb[mask_npy]
    cell_a   = lut_alpha[mask_npy, None]
    fg       = (mask_npy > 0)[:, :, None]
    blended  = np.where(fg, (1 - cell_a) * he_f + cell_a * cell_rgb, he_f)
    bounds   = find_boundaries(mask_npy, mode="thin")
    blended[bounds] = [0.15, 0.15, 0.15]
    comp = (blended * 255).clip(0, 255).astype(np.uint8)
    at2_n   = int((cell_annotation == "AT2 Pneumocyte").sum())
    total_n = adata.n_obs
    at2_pct = at2_n / total_n
    return comp, at2_n, total_n, at2_pct


def run_method(h5ad_name: str, mask_name: str, he: np.ndarray, label: str):
    print(f"\n[{label}] Loading {h5ad_name}...")
    adata = ad.read_h5ad(str(H5AD_DIR / h5ad_name))
    print(f"  Cells: {adata.n_obs}")
    X_log = normalize_log1p(adata.X)
    leiden_labels = leiden_cluster(X_log, seed=42)
    at2_cl = identify_at2_cluster(leiden_labels, X_log, np.array(adata.var_names))
    mask_npy = np.load(str(MASK_DIR / mask_name))
    comp, at2_n, total_n, at2_pct = build_composite(adata, mask_npy, he, leiden_labels, at2_cl)
    print(f"  AT2: n={at2_n} ({at2_pct:.1%}), Total: {total_n}")
    return comp, at2_n, total_n, at2_pct


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    print("[0] Loading H&E crop...")
    he = load_he_crop(ROI10)
    h_px, w_px = he.shape[:2]

    comp_2cseg, at2_n_2c, total_2c, pct_2c = run_method(
        "roi10_cellpose_dilate.h5ad", "vhd_roi10_cellpose_dilate.npy", he, "2Cseg")
    comp_mcseg, at2_n_mc, total_mc, pct_mc = run_method(
        "roi10_v12.h5ad", "vhd_roi10_v12.npy", he, "MCseg")

    print("\n[Plotting]...")
    fig_w = 183 * MM
    fig_h = fig_w * h_px / w_px * 0.55 + 14 * MM

    fig = plt.figure(figsize=(fig_w, fig_h))
    gs = fig.add_gridspec(1, 2, left=0.005, right=0.995,
                          top=0.88, bottom=0.05,
                          wspace=0.08, width_ratios=[3, 1])
    ax_map = fig.add_subplot(gs[0])
    ax_bar = fig.add_subplot(gs[1])

    # ── Panel a: 2Cseg spatial map ────────────────────────────────────────────
    ax_map.imshow(comp_2cseg, aspect="auto", interpolation="bilinear")
    ax_map.axis("off")
    ax_map.set_title("2Cseg — ROI 10 (Normal Alveolar Region)",
                     fontsize=9, fontweight="bold", pad=4)

    unresolved_n = total_2c - at2_n_2c
    leg_handles = [
        mpatches.Patch(color=PALETTE_HEX["AT2 Pneumocyte"],
                       label=f"AT2 Pneumocyte  (n={at2_n_2c}, {pct_2c:.0%})"),
        mpatches.Patch(color=PALETTE_HEX["Unresolved"],
                       label=f"Unresolved  (n={unresolved_n}, {1 - pct_2c:.0%})"),
    ]
    ax_map.legend(handles=leg_handles, fontsize=6.5, loc="lower right",
                  framealpha=0.88, edgecolor="#ccc", handlelength=1.2)
    ax_map.text(-0.005, 1.08, "a", transform=ax_map.transAxes,
                fontsize=12, fontweight="bold", va="top", ha="right")

    px_per_um = 1 / 0.2737
    scale_px  = 50 * px_per_um
    margin_x  = w_px * 0.04
    bar_y     = h_px - h_px * 0.05
    ax_map.plot([margin_x, margin_x + scale_px], [bar_y, bar_y],
                color="white", linewidth=3, solid_capstyle="butt", zorder=10)
    ax_map.text((margin_x + margin_x + scale_px) / 2, bar_y - h_px * 0.025,
                "50 µm", color="white", ha="center", va="bottom",
                fontsize=7, fontweight="bold", zorder=10)

    # ── Panel b: bar chart ────────────────────────────────────────────────────
    methods  = ["2Cseg", "MCseg"]
    pcts_pct = [pct_2c * 100, pct_mc * 100]
    totals_n = [total_2c, total_mc]
    bar_colors = ["#7BAFD4", PALETTE_HEX["AT2 Pneumocyte"]]
    bars = ax_bar.bar(methods, pcts_pct, color=bar_colors, width=0.55,
                      edgecolor="none")

    delta = pct_mc * 100 - pct_2c * 100
    ax_bar.annotate(
        f"+{delta:.1f} pp",
        xy=(0.5, (pct_2c + pct_mc) / 2 * 100),
        xytext=(0.5, pct_mc * 100 + 2.5),
        xycoords=("data", "data"), textcoords=("data", "data"),
        ha="center", fontsize=7.5, color="#333333",
    )

    ax_bar.set_ylabel("AT2 Pneumocytes (%)", fontsize=8)
    ax_bar.set_ylim(0, pct_mc * 100 * 1.4)
    ax_bar.set_title("AT2 Pneumocyte\nDetection Rate",
                     fontsize=8, fontweight="bold", pad=4)
    ax_bar.spines[["top", "right"]].set_visible(False)
    ax_bar.tick_params(axis="both", labelsize=8)
    for bar, n in zip(bars, totals_n):
        ax_bar.text(bar.get_x() + bar.get_width() / 2,
                    bar.get_height() + 0.3,
                    f"n={n}", ha="center", va="bottom", fontsize=6.5)
    ax_bar.text(-0.2, 1.08, "b", transform=ax_bar.transAxes,
                fontsize=12, fontweight="bold", va="top", ha="right")

    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PATH, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"\n✅ Saved: {OUT_PATH}")


if __name__ == "__main__":
    main()
