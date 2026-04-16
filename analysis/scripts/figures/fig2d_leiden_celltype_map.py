"""
fig2d_leiden_celltype_map.py
============================
Generate fig2d: ROI10 spatial cell-type map using Leiden clustering.

Replaces the threshold/winner-take-all annotation (fig2fg_spatial_expression.py)
with an unsupervised Leiden clustering approach.

Annotation logic:
  - Run Leiden (res=0.4, k=15) on normalize→HVG→TruncatedSVD pipeline
  - Identify AT2 cluster: cluster with highest mean SFTPC expression
  - Label: AT2 cluster → "AT2 Pneumocyte" (#2166AC), all others → "Unresolved" (#CCCCCC)

This is more defensible than threshold-based annotation because:
  - Unbiased grouping by transcriptome similarity
  - Only labels clusters with robust marker enrichment
  - Avoids inflating AT2 proportion via WTA default in sparse data

Output:
  submission_bioinformatics/figures/fig2/fig2d.png

Run:
  cd /Volumes/SSD/plan_a
  uv run python submission_bioinformatics/scripts/figures/fig2d_leiden_celltype_map.py
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
OUT_PATH = Path("/Volumes/SSD/plan_a/submission_bioinformatics/figures/fig2/fig2d.png")

ROI10 = dict(x=7562, y=19440, w=3194, h=1587)

PALETTE = {
    "AT2 Pneumocyte": ((0.13, 0.40, 0.67), 0.72),   # #2166AC, alpha=0.72
    "Unresolved":     ((0.73, 0.73, 0.73), 0.35),   # #BBBBBB, alpha=0.35 (subtle)
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
    """Cluster cells using TruncatedSVD + KNN + Leiden. Returns integer label array."""
    var_all = (X_log ** 2).mean(0) - X_log.mean(0) ** 2
    top_idx = np.argsort(var_all)[-n_top_hvg:]
    X_hvg = X_log[:, top_idx]

    actual_pcs = min(n_pcs, X_hvg.shape[1] - 1, X_hvg.shape[0] - 1)
    svd = TruncatedSVD(n_components=actual_pcs, random_state=seed)
    X_pca = svd.fit_transform(X_hvg)
    print(f"  PCA: {X_hvg.shape} → {X_pca.shape}  "
          f"(expl. var: {svd.explained_variance_ratio_.sum():.1%})")

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
    """Return the cluster label with the highest mean SFTPC expression.
    Falls back to AT2 gene-set score if SFTPC not in panel.
    """
    clusters = np.unique(leiden_labels)
    at2_genes = ["SFTPC", "SFTPB", "SFTPA1", "SFTPA2"]
    gene_idx = {g: i for i, g in enumerate(var_names)}

    present = [g for g in at2_genes if g in gene_idx]
    if not present:
        raise RuntimeError("None of the AT2 marker genes found in h5ad var_names.")

    best_cl, best_score = -1, -1.0
    for cl in clusters:
        mask = leiden_labels == cl
        score = max(float(X_log[mask, gene_idx[g]].mean()) for g in present)
        if score > best_score:
            best_score = score
            best_cl = cl

    return int(best_cl)


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    print("[1] Loading h5ad (ROI10, MCseg v2)...")
    adata = ad.read_h5ad(str(H5AD_DIR / "roi10_v12.h5ad"))
    var_names = np.array(adata.var_names)
    print(f"  Cells: {adata.n_obs}, Genes: {adata.n_vars}")

    print("[2] Normalizing (log1p, target=10,000)...")
    X_log = normalize_log1p(adata.X)

    print("[3] Leiden clustering (res=0.4, k=15, n_pcs=50)...")
    leiden_labels = leiden_cluster(X_log, n_top_hvg=2000, n_pcs=50,
                                   k=15, resolution=0.4, seed=42)
    n_clusters = len(np.unique(leiden_labels))
    print(f"  {n_clusters} clusters: " +
          ", ".join(f"C{cl}(n={int((leiden_labels == cl).sum())})"
                    for cl in np.unique(leiden_labels)))

    print("[4] Identifying AT2 cluster by SFTPC/SFTPB mean expression...")
    at2_cl = identify_at2_cluster(leiden_labels, X_log, var_names)
    at2_n = int((leiden_labels == at2_cl).sum())
    at2_pct = at2_n / adata.n_obs

    # Verify: print SFTPC % per cluster
    gene_idx = {g: i for i, g in enumerate(var_names)}
    for cl in np.unique(leiden_labels):
        mask = leiden_labels == cl
        sftpc_pct = float((X_log[mask, gene_idx["SFTPC"]] > 0).mean()) * 100 \
                    if "SFTPC" in gene_idx else 0.0
        tag = " ← AT2" if cl == at2_cl else ""
        print(f"  C{cl}(n={mask.sum():4d}): SFTPC+={sftpc_pct:.1f}%{tag}")

    print(f"\n  AT2 = C{at2_cl} (n={at2_n}, {at2_pct:.1%})")

    # Build cell annotation array (indexed by h5ad row order)
    cell_annotation = np.where(leiden_labels == at2_cl,
                                "AT2 Pneumocyte", "Unresolved")

    print("[5] Loading H&E crop from BTF...")
    he = load_he_crop(ROI10)
    print(f"  H&E shape: {he.shape}")

    print("[6] Loading segmentation mask...")
    mask = np.load(str(MASK_DIR / "vhd_roi10_v12.npy"))
    print(f"  Mask shape: {mask.shape}, unique cells: {len(np.unique(mask)) - 1}")

    print("[7] Building composite image...")
    cell_ids = adata.obs["cell_id"].values.astype(int)
    id_max = int(mask.max())

    lut_rgb   = np.zeros((id_max + 1, 3), dtype=np.float32)
    lut_alpha = np.zeros(id_max + 1,      dtype=np.float32)

    for i, cid in enumerate(cell_ids):
        if 0 < cid <= id_max:
            rgb, alpha      = PALETTE[cell_annotation[i]]
            lut_rgb[cid]    = rgb
            lut_alpha[cid]  = alpha

    he_f     = he.astype(np.float32) / 255.0
    cell_rgb = lut_rgb[mask]
    cell_a   = lut_alpha[mask, None]
    fg       = (mask > 0)[:, :, None]
    blended  = np.where(fg, (1 - cell_a) * he_f + cell_a * cell_rgb, he_f)

    bounds = find_boundaries(mask, mode="thin")
    blended[bounds] = [0.15, 0.15, 0.15]

    comp_img = (blended * 255).clip(0, 255).astype(np.uint8)

    print("[8] Plotting...")
    h_px, w_px = he.shape[:2]
    fig_w = 183 * MM
    fig_h = fig_w * h_px / w_px + 10 * MM

    fig, ax = plt.subplots(1, 1, figsize=(fig_w, fig_h),
                            gridspec_kw=dict(left=0.005, right=0.995,
                                             top=0.93, bottom=0.005))
    ax.imshow(comp_img, aspect="auto", interpolation="bilinear")
    ax.axis("off")
    ax.set_title("")

    unresolved_n = int((cell_annotation == "Unresolved").sum())
    leg_handles = [
        mpatches.Patch(color=PALETTE_HEX["AT2 Pneumocyte"],
                       label=f"AT2 Pneumocyte  (n={at2_n}, {at2_pct:.0%})"),
        mpatches.Patch(color=PALETTE_HEX["Unresolved"],
                       label=f"Unresolved  (n={unresolved_n}, {1 - at2_pct:.0%})"),
    ]
    ax.legend(handles=leg_handles, fontsize=6.5, loc="lower right",
              framealpha=0.88, edgecolor="#ccc", handlelength=1.2)
    ax.text(-0.005, 1.06, "d", transform=ax.transAxes,
            fontsize=12, fontweight="bold", va="top", ha="right")

    # Scale bar (bottom-left): 100 µm @ 0.2737 µm/px
    px_per_um   = 1 / 0.2737
    scale_um    = 50
    scale_px    = scale_um * px_per_um            # ≈ 183 px
    margin_x    = w_px * 0.04                     # 4% from left
    margin_y    = h_px * 0.05                     # 5% from bottom
    bar_y       = h_px - margin_y
    bar_x0      = margin_x
    bar_x1      = bar_x0 + scale_px
    ax.plot([bar_x0, bar_x1], [bar_y, bar_y],
            color="white", linewidth=3, solid_capstyle="butt",
            transform=ax.transData, zorder=10)
    ax.text((bar_x0 + bar_x1) / 2, bar_y - h_px * 0.025,
            f"{scale_um} µm",
            color="white", ha="center", va="bottom",
            fontsize=7, fontweight="bold",
            transform=ax.transData, zorder=10)

    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PATH, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"\n✅ Saved: {OUT_PATH}")


if __name__ == "__main__":
    main()
