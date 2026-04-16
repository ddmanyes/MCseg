"""
33_fig2_spatial_maps.py
========================
Generate:
  fig2de.png — d: ROI10 H&E + AT2/AT1/Other scatter map
               e: AT1/AT2 dotplot (7 genes × 3 groups)

  fig2fg.png — f: ROI9 H&E + Leiden cluster scatter map (5 clusters)
               g: Leiden cluster marker dotplot

Clustering uses sklearn + leidenalg to avoid macOS OpenBLAS deadlock.
Output: manuscript/figures/03_luad_benchmark/
"""

from __future__ import annotations
import warnings
warnings.filterwarnings("ignore")

import os
os.environ["OMP_NUM_THREADS"] = "1"   # prevent BLAS mutex deadlock on macOS

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from pathlib import Path
import tifffile
import zarr
import scipy.sparse as sp
import anndata as ad
from sklearn.decomposition import TruncatedSVD
from sklearn.neighbors import NearestNeighbors
from skimage.measure import regionprops
from skimage.segmentation import find_boundaries
from scipy import ndimage as ndi
import leidenalg
import igraph

# ── Paths ────────────────────────────────────────────────────────────────────
PROJECT_ROOT = Path("/Volumes/SSD/plan_a/xenium_he_seg")
BTF_PATH = (Path("/Volumes/SSD/plan_a/tissue sample/LUAD/visium") /
            "Visium_HD_Human_Lung_Cancer_post_Xenium_Prime_5K_Experiment2_tissue_image.btf")
MASK_DIR = PROJECT_ROOT / "results" / "masks"
H5AD_DIR = PROJECT_ROOT / "results" / "visiumhd" / "visiumhd_cells"
OUT_DIR  = Path("/Volumes/SSD/plan_a/manuscript/figures/03_luad_benchmark")
OUT_DIR.mkdir(parents=True, exist_ok=True)

ROI9  = dict(x=5869,  y=11834, w=3854, h=3315)
ROI10 = dict(x=7562,  y=19440, w=3194, h=1587)

plt.rcParams.update({
    "font.family": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 8, "axes.titlesize": 8, "axes.labelsize": 7.5,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "axes.linewidth": 0.8,
    "pdf.fonttype": 42,
})
MM = 1 / 25.4


# ── Helpers ──────────────────────────────────────────────────────────────────

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


def get_centroids(mask_path: str) -> dict[int, tuple[float, float]]:
    """Returns {label: (col_x, row_y)} in mask-local pixels."""
    mask = np.load(mask_path)
    return {p.label: (float(p.centroid[1]), float(p.centroid[0]))
            for p in regionprops(mask)}


def normalize_log1p(adata_raw: ad.AnnData) -> np.ndarray:
    """Return dense log1p-normalized matrix (cells × genes)."""
    X = sp.csr_matrix(adata_raw.X)
    totals = np.array(X.sum(axis=1)).flatten()
    X = X.multiply(1e4 / totals[:, None]).tocsr()
    X.data = np.log1p(X.data)
    return X.toarray()


def leiden_cluster(X_log: np.ndarray, n_top_hvg: int = 2000,
                    n_pcs: int = 30, k: int = 15,
                    resolution: float = 0.3, seed: int = 42) -> np.ndarray:
    """Cluster cells using TruncatedSVD + KNN + Leiden. Returns integer label array."""
    # Top HVG by variance
    mean_sq = X_log.mean(axis=0) ** 2
    var_all = ((X_log ** 2).mean(axis=0)) - mean_sq
    top_idx = np.argsort(var_all)[-n_top_hvg:]
    X_hvg   = X_log[:, top_idx]

    # PCA
    svd = TruncatedSVD(n_components=n_pcs, random_state=seed)
    X_pca = svd.fit_transform(X_hvg)

    # KNN
    nn = NearestNeighbors(n_neighbors=k + 1, algorithm="auto", n_jobs=1)
    nn.fit(X_pca)
    _, idxs = nn.kneighbors(X_pca)

    # igraph + Leiden
    n = X_pca.shape[0]
    edges = [(int(i), int(j)) for i in range(n) for j in idxs[i, 1:]]
    g = igraph.Graph(n=n, edges=edges, directed=False)
    g.simplify()
    part = leidenalg.find_partition(
        g, leidenalg.RBConfigurationVertexPartition,
        resolution_parameter=resolution, seed=seed)
    return np.array(part.membership)


def find_top_markers(X_log: np.ndarray, labels: np.ndarray,
                      var_names: np.ndarray, top_n: int = 3,
                      min_pct: float = 0.05,
                      force_genes: list[str] | None = None
                      ) -> tuple[list[str], dict]:
    """Find top marker genes per cluster by mean log1p fold change.
    Returns (selected_genes, {cluster: [genes]}).
    force_genes: additional genes always included at the end.
    """
    clusters = sorted(np.unique(labels), key=int)
    selected: list[str] = []
    cl_gene_map: dict = {}

    for cl in clusters:
        mask = labels == cl
        mean_in  = X_log[mask].mean(axis=0)
        pct_in   = (X_log[mask] > 0).mean(axis=0)
        mean_out = X_log[~mask].mean(axis=0)
        score    = (mean_in - mean_out) * (pct_in >= min_pct)
        top_idx  = np.argsort(-score)
        cl_genes = [var_names[i] for i in top_idx
                    if var_names[i] not in selected][:top_n]
        selected.extend(cl_genes)
        cl_gene_map[cl] = cl_genes

    if force_genes:
        extra = [g for g in force_genes if g not in selected]
        selected.extend(extra)
        cl_gene_map["_extra"] = extra

    return selected, cl_gene_map


def build_dotplot_arrays(X_log: np.ndarray, labels: np.ndarray,
                          genes: list[str], var_names: np.ndarray
                          ) -> tuple[np.ndarray, np.ndarray]:
    """Compute pct_positive and mean_log1p per (cluster × gene)."""
    clusters = sorted(np.unique(labels), key=int)
    gene_idx = {g: i for i, g in enumerate(var_names)}
    n_gr, n_g = len(clusters), len(genes)
    pct_arr  = np.zeros((n_gr, n_g))
    mean_arr = np.zeros((n_gr, n_g))

    for ci, cl in enumerate(clusters):
        mask = labels == cl
        for gi, gene in enumerate(genes):
            if gene not in gene_idx:
                continue
            expr = X_log[mask, gene_idx[gene]]
            pct_arr[ci, gi]  = (expr > 0).mean() * 100
            mean_arr[ci, gi] = expr.mean()
    return pct_arr, mean_arr


def render_dotplot(ax, fig, pct_arr: np.ndarray, mean_arr: np.ndarray,
                    x_labels: list[str], y_labels: list[str],
                    title: str = "", vmax_mean: float | None = None,
                    max_dot_s: float = 400, legends_outside: bool = False):
    """Render dotplot onto ax. pct_arr/mean_arr shape: (n_groups × n_genes).
    y_labels listed top→bottom.
    """
    n_gr, n_g = pct_arr.shape
    if vmax_mean is None:
        vmax_mean = max(mean_arr.max(), 0.05)
    norm = mcolors.Normalize(vmin=0, vmax=vmax_mean)
    cmap = plt.cm.YlOrRd

    xs, ys, ss, cs = [], [], [], []
    for gi in range(n_gr):
        y_pos = n_gr - 1 - gi
        for gei in range(n_g):
            xs.append(gei)
            ys.append(y_pos)
            ss.append(max((pct_arr[gi, gei] / 100) * max_dot_s, 3))
            cs.append(mean_arr[gi, gei])

    sc_obj = ax.scatter(xs, ys, s=ss, c=cs, cmap=cmap, norm=norm,
                        edgecolors="grey", linewidths=0.3, zorder=3)

    ax.set_xlim(-0.6, n_g - 0.4)
    ax.set_ylim(-0.6, n_gr - 0.4)
    ax.set_xticks(range(n_g))
    ax.set_xticklabels(x_labels, fontsize=6.5, rotation=50, ha="right")
    ax.set_yticks(range(n_gr))
    ax.set_yticklabels(y_labels[::-1], fontsize=7)
    ax.grid(alpha=0.18, lw=0.4)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if title:
        ax.set_title(title, fontsize=7, fontweight="bold", pad=4)

    # Size legend circles — same scale as scatter (s=area, markersize=diameter)
    # markersize = 2 * sqrt(s / π)  so legend visually matches scatter dots
    leg_h = [
        Line2D([0], [0], marker="o", color="w",
               markerfacecolor="#aaa", markeredgecolor="grey",
               markeredgewidth=0.4,
               markersize=np.sqrt(max((p / 100) * max_dot_s, 2) / np.pi) * 2,
               label=f"{p}%")
        for p in [25, 50, 75]
    ]

    if legends_outside:
        # Colorbar: true outside placement via figure axes
        ax_pos = ax.get_position()
        cbar_ax = fig.add_axes([
            ax_pos.x1 + 0.012,
            ax_pos.y0 + ax_pos.height * 0.42,
            0.013,
            ax_pos.height * 0.52,
        ])
        cb = fig.colorbar(sc_obj, cax=cbar_ax)
        cb.set_label("Mean log(1+UMI)", fontsize=5.5, labelpad=2)
        cb.ax.tick_params(labelsize=5.5, pad=1)
        cb.set_ticks(np.linspace(0, vmax_mean, 5))
        # Size legend below colorbar
        ax.legend(handles=leg_h, title="% positive",
                  bbox_to_anchor=(1.02, 0.36), loc="upper left",
                  fontsize=5.5, title_fontsize=5.5,
                  handletextpad=0.3, labelspacing=0.3,
                  framealpha=0.85, edgecolor="#ccc")
    else:
        cbar_ax = ax.inset_axes([0.88, 0.38, 0.026, 0.55])
        cb = fig.colorbar(sc_obj, cax=cbar_ax)
        cb.set_label("Mean log(1+UMI)", fontsize=5.5, labelpad=2)
        cb.ax.tick_params(labelsize=5.5, pad=1)
        cb.set_ticks(np.linspace(0, vmax_mean, 5))
        ax.legend(handles=leg_h, title="% positive",
                  bbox_to_anchor=(1.0, 0.33), loc="upper right",
                  fontsize=5.5, title_fontsize=5.5,
                  handletextpad=0.3, labelspacing=0.3,
                  framealpha=0.85, edgecolor="#ccc")


def get_expr_dense(X_log: np.ndarray, gene: str,
                    var_names: np.ndarray) -> np.ndarray:
    if gene not in var_names:
        return np.zeros(X_log.shape[0])
    idx = np.where(var_names == gene)[0][0]
    return X_log[:, idx]


# ROI9 cluster → cell type annotation (Leiden res=0.3, n=5 clusters)
# Cross-validated against Xenium celltype_summary.csv (ROI9)
# C0 CD74+FTH1+APOE+VIM          → Macrophage         (FTH1=0.80, APOE=0.24)
# C1 SCGB1A1+++BPIFB1+KRT8+EPCAM → Club Epithelial    (SCGB1A1=5.87)
# C2 IGKC+IGHG1+++JCHAIN+         → Plasma Cell        (IGKC=7.05, IGHG1=6.09)
# C3 IGKC+MS4A1+IGHA1+IGHM+       → B Cell             (MS4A1=0.34, PAX5+)
# C4 SPP1+FTH1(2.48)+CSF1R+LGMN   → SPP1+ Macrophage   (n=248 ≈ VHD 4.12% SPP1+Mac)
ROI9_CLUSTER_NAMES = {
    0: "Macrophage",
    1: "Club Epithelial",
    2: "Plasma Cell",
    3: "B Cell",
    4: "SPP1\u207a Macrophage",
}

# Curated marker genes per cluster (literature-based, ordered for dotplot display)
# C0: canonical macrophage markers detectable in Visium HD
#   CD74 (MHC-II invariant chain, antigen presentation), VIM (mesenchymal, macrophage-expressed),
#   CXCR4 (chemokine receptor, monocyte/macrophage homing)
# C1: airway secretory / Club cell markers (highly specific)
# C2: plasma cell immunoglobulin chains (IGHG1 distinguishes from B cells)
# C3: B cell markers (MS4A1=CD20, IGHM=IgM naive marker, CD79A=B cell receptor)
# C4: SPP1+ macrophage (SPP1=osteopontin, FTH1=ferritin heavy chain/iron-laden 55% in C4,
#     CD68=lysosomal marker, TREM2=immunomodulatory receptor 12x enriched in C4 vs background)
ROI9_CURATED_GENES = {
    0: ["CD74",    "VIM",    "CXCR4"  ],   # Macrophage
    1: ["SCGB1A1", "BPIFB1", "WFDC2"  ],   # Club Epithelial
    2: ["IGHG1",   "IGKC",   "JCHAIN" ],   # Plasma Cell
    3: ["MS4A1",   "IGHM",   "CD79A"  ],   # B Cell
    4: ["SPP1",    "FTH1",   "CD68",  "TREM2"],  # SPP1+ Macrophage
}

# ═════════════════════════════════════════════════════════════════════════════
# Shared ROI10 data loader
# ═════════════════════════════════════════════════════════════════════════════

LUNG_GENE_SETS = {
    "AT2 Pneumocyte":      ["SFTPC", "SFTPB", "SFTPA1", "SFTPA2"],
    "Alveolar Macrophage": ["MARCO", "FABP4", "MCEMP1", "SPP1"],
    "Endothelial":         ["PECAM1", "VWF", "CLDN5"],
    "AT1 Pneumocyte":      ["AGER", "RTKN2", "CAV1", "HOPX"],
}
# Colours matching 11_normal_lung_celltype.py  +  Unassigned
CELL_TYPE_PALETTE = {
    "AT2 Pneumocyte":      ((0.13, 0.40, 0.67), 0.70),  # Blue    #2166AC
    "Alveolar Macrophage": ((0.83, 0.52, 0.10), 0.72),  # Orange  #D4841A
    "Endothelial":         ((0.10, 0.45, 0.25), 0.70),  # Green   #1A7340
    "AT1 Pneumocyte":      ((0.70, 0.09, 0.17), 0.72),  # Red     #B2182B
    "Unassigned":          ((0.55, 0.55, 0.55), 0.30),  # Gray
}
# hex for legend patches
CELL_TYPE_HEX = {
    "AT2 Pneumocyte":      "#2166AC",
    "Alveolar Macrophage": "#D4841A",
    "Endothelial":         "#1A7340",
    "AT1 Pneumocyte":      "#B2182B",
    "Unassigned":          "#8C8C8C",
}


def _load_roi10_data():
    """Hybrid cell typing on ROI10.

    AT1 Pneumocyte: threshold-based (AGER > 0 OR RTKN2 > 0), applied first.
    All other types: winner-take-all among remaining cells.

    Returns (he, comp_img, type_counts)
      type_counts: dict {type_name: n_cells}
    """
    adata = ad.read_h5ad(str(H5AD_DIR / "roi10_v12.h5ad"))
    gene_names = np.array(adata.var_names)
    X_log = normalize_log1p(adata)

    n_cells   = adata.n_obs
    cell_type = np.full(n_cells, "Unassigned", dtype=object)

    # ── Step 1: AT1 threshold (AGER > 0 OR RTKN2 > 0) ───────────────────
    at1_mask = np.zeros(n_cells, dtype=bool)
    for g in ["AGER", "RTKN2"]:
        if g in gene_names:
            at1_mask |= (get_expr_dense(X_log, g, gene_names) > 0)
    cell_type[at1_mask] = "AT1 Pneumocyte"

    # ── Step 2: winner-take-all for remaining cells ───────────────────────
    remaining = ~at1_mask
    wta_types = ["AT2 Pneumocyte", "Alveolar Macrophage", "Endothelial"]
    wta_genes = {t: LUNG_GENE_SETS[t] for t in wta_types}

    scores = np.zeros((n_cells, len(wta_types)), dtype=np.float32)
    for ci, (tname, markers) in enumerate(wta_genes.items()):
        for g in markers:
            if g in gene_names:
                expr = get_expr_dense(X_log, g, gene_names)
                scores[:, ci] = np.maximum(scores[:, ci], expr)

    winner_idx = scores.argmax(axis=1)
    max_score  = scores.max(axis=1)
    wta_result = np.where(max_score > 0,
                          np.array(wta_types)[winner_idx],
                          "Unassigned")
    cell_type[remaining] = wta_result[remaining]

    cell_ids = adata.obs["cell_id"].values.astype(int)
    all_types = list(LUNG_GENE_SETS.keys()) + ["Unassigned"]
    type_counts = {t: int((cell_type == t).sum()) for t in all_types}

    print(f"  Cell type breakdown: { {k:v for k,v in type_counts.items() if v>0} }")

    # ── Build composite image ─────────────────────────────────────────────
    print("  Loading mask & building composite ROI10...")
    mask = np.load(str(MASK_DIR / "vhd_roi10_v12.npy"))
    he   = load_he_crop(ROI10)

    id_max    = int(mask.max())
    lut_rgb   = np.zeros((id_max + 1, 3), dtype=np.float32)
    lut_alpha = np.zeros(id_max + 1,      dtype=np.float32)

    for i, cid in enumerate(cell_ids):
        if 0 < cid <= id_max:
            rgb, alpha       = CELL_TYPE_PALETTE[cell_type[i]]
            lut_rgb[cid]     = rgb
            lut_alpha[cid]   = alpha

    he_f     = he.astype(np.float32) / 255.0
    cell_rgb = lut_rgb[mask]
    cell_a   = lut_alpha[mask, None]
    fg       = (mask > 0)[:, :, None]
    blended  = np.where(fg, (1 - cell_a) * he_f + cell_a * cell_rgb, he_f)

    bounds = find_boundaries(mask, mode="thin")
    blended[bounds] = [0.15, 0.15, 0.15]

    comp_img = (blended * 255).clip(0, 255).astype(np.uint8)
    return he, comp_img, type_counts


# ═════════════════════════════════════════════════════════════════════════════
# Fig. 2d — ROI10 cell-type composite map  (standalone, full width)
# ═════════════════════════════════════════════════════════════════════════════

def make_fig2d():
    print("=== fig2d: ROI10 cell-type map (5 types, winner-take-all) ===")
    he, comp_img, type_counts = _load_roi10_data()

    n_total = sum(type_counts.values())
    h_px, w_px = he.shape[:2]
    fig_w = 183 * MM
    fig_h = fig_w * h_px / w_px + 10 * MM

    fig, ax = plt.subplots(1, 1, figsize=(fig_w, fig_h),
                            gridspec_kw=dict(left=0.005, right=0.995,
                                             top=0.93, bottom=0.005))
    ax.imshow(comp_img, aspect="auto", interpolation="bilinear")
    ax.axis("off")
    ax.set_title("ROI10 — Cell-type Map  (Normal Alveolar, MCseg v2, n=1,085)",
                 fontsize=8, fontweight="bold", pad=3)

    leg_handles = []
    for tname, hex_col in CELL_TYPE_HEX.items():
        n = type_counts.get(tname, 0)
        if n == 0:
            continue
        pct = n / n_total
        leg_handles.append(
            mpatches.Patch(color=hex_col,
                           label=f"{tname}  (n={n}, {pct:.0%})"))
    ax.legend(handles=leg_handles, fontsize=6, loc="lower right",
              framealpha=0.88, edgecolor="#ccc")
    ax.text(-0.005, 1.04, "d", transform=ax.transAxes,
            fontsize=12, fontweight="bold", va="top", ha="right")

    out = OUT_DIR / "fig2d.png"
    fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  ✅ Saved: {out}")


# ═════════════════════════════════════════════════════════════════════════════
# Fig. 2e — AT2/AT1 dotplot  (standalone, narrower)
# ═════════════════════════════════════════════════════════════════════════════

def make_fig2e():
    """Fig. 2e: dotplot for 5 cell types × marker genes (ROI10, Visium HD).
    Cell type labels match fig2d (winner-take-all).
    """
    print("=== fig2e: 5-type dotplot ===")

    # ── Reproduce winner-take-all classification ──────────────────────────
    adata = ad.read_h5ad(str(H5AD_DIR / "roi10_v12.h5ad"))
    gene_names = np.array(adata.var_names)
    X_log = normalize_log1p(adata)

    type_names = list(LUNG_GENE_SETS.keys())
    scores = np.zeros((adata.n_obs, len(type_names)), dtype=np.float32)
    for ci, (_, markers) in enumerate(LUNG_GENE_SETS.items()):
        for g in markers:
            if g in gene_names:
                scores[:, ci] = np.maximum(scores[:, ci],
                                            get_expr_dense(X_log, g, gene_names))
    winner_idx = scores.argmax(axis=1)
    max_score  = scores.max(axis=1)
    cell_type  = np.where(max_score > 0,
                          np.array(type_names)[winner_idx], "Unassigned")

    # ── Marker genes to display (grouped by type) ─────────────────────────
    GENE_GROUPS = [
        ("AT2 Pneumocyte",      ["SFTPC", "SFTPB", "SFTPA1"]),
        ("Alveolar Macrophage", ["MARCO", "FABP4", "SPP1"]),
        ("Endothelial",         ["PECAM1", "VWF"]),
        ("AT1 Pneumocyte",      ["AGER", "CAV1", "HOPX"]),
    ]
    DOT_GENES = [g for _, gs in GENE_GROUPS for g in gs if g in gene_names]

    # Row order: AT2 → AM → Endo → AT1 → Unassigned (top → bottom)
    ROW_ORDER = type_names + ["Unassigned"]
    DOT_LABELS = ["AT2 Pneumocyte", "Alv. Macrophage",
                  "Endothelial", "AT1 Pneumocyte", "Unassigned"]

    # ── Compute pct & mean per group × gene ──────────────────────────────
    n_gr, n_g = len(ROW_ORDER), len(DOT_GENES)
    pct_arr  = np.zeros((n_gr, n_g))
    mean_arr = np.zeros((n_gr, n_g))
    gene_idx = {g: i for i, g in enumerate(gene_names)}

    for ri, rtype in enumerate(ROW_ORDER):
        mask_r = cell_type == rtype
        for gi, gene in enumerate(DOT_GENES):
            if gene not in gene_idx:
                continue
            expr = X_log[mask_r, gene_idx[gene]]
            pct_arr[ri, gi]  = (expr > 0).mean() * 100
            mean_arr[ri, gi] = expr.mean()

    # ── Tick label colours (match cell-type group colour) ─────────────────
    gene_to_type = {g: tname for tname, gs in GENE_GROUPS for g in gs}
    tick_cols = [CELL_TYPE_HEX.get(gene_to_type.get(g, ""), "#333333")
                 for g in DOT_GENES]

    # ── Vertical separators between gene groups ───────────────────────────
    sep_positions = []
    pos = 0
    for _, gs in GENE_GROUPS[:-1]:
        pos += sum(1 for g in gs if g in gene_names)
        sep_positions.append(pos - 0.5)

    # ── Figure ───────────────────────────────────────────────────────────
    fig_w = 130 * MM
    fig_h =  88 * MM
    fig, ax = plt.subplots(1, 1, figsize=(fig_w, fig_h),
                            gridspec_kw=dict(left=0.18, right=0.78,
                                             top=0.88, bottom=0.22))

    render_dotplot(ax, fig, pct_arr, mean_arr,
                   x_labels=DOT_GENES, y_labels=DOT_LABELS,
                   title="Cell-type Marker Genes — ROI10 (Visium HD)",
                   legends_outside=True)

    for tick, col in zip(ax.get_xticklabels(), tick_cols):
        tick.set_color(col)
    for sp in sep_positions:
        ax.axvline(sp, color="grey", lw=0.6, ls="--", alpha=0.40)

    ax.text(-0.22, 1.06, "e", transform=ax.transAxes,
            fontsize=12, fontweight="bold", va="top", ha="right")

    out = OUT_DIR / "fig2e.png"
    fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  ✅ Saved: {out}")


# ═════════════════════════════════════════════════════════════════════════════
# Shared ROI9 data loader
# ═════════════════════════════════════════════════════════════════════════════

def _load_roi9_data():
    """Leiden clustering + composite image for ROI9.

    Returns (he, comp_img, labels, clusters, X_log, gene_names, cmap_cl, selected, cl_gene_map)
    """
    adata = ad.read_h5ad(str(H5AD_DIR / "roi9_v12.h5ad"))
    gene_names = np.array(adata.var_names)
    print("  Normalizing...")
    X_log = normalize_log1p(adata)

    print("  Leiden clustering (res=0.3)...")
    labels   = leiden_cluster(X_log, resolution=0.3, seed=42)
    clusters = sorted(np.unique(labels), key=int)
    n_cl     = len(clusters)
    print(f"  {n_cl} clusters")

    selected, cl_gene_map = find_top_markers(
        X_log, labels, gene_names, top_n=3, min_pct=0.05,
        force_genes=["SPP1", "CD68"],
    )

    # Colormap
    try:
        cmap_cl = matplotlib.colormaps["tab10"].resampled(max(n_cl, 2))
    except AttributeError:
        cmap_cl = matplotlib.cm.get_cmap("tab10", max(n_cl, 2))

    # Build composite image (mask-based, same style as fig2d)
    cell_ids = adata.obs["cell_id"].values.astype(int)
    print("  Loading mask & building composite ROI9...")
    mask = np.load(str(MASK_DIR / "vhd_roi9_v12.npy"))
    he   = load_he_crop(ROI9)

    id_max    = int(mask.max())
    lut_rgb   = np.zeros((id_max + 1, 3), dtype=np.float32)
    lut_alpha = np.zeros(id_max + 1,      dtype=np.float32)

    for ci, cl in enumerate(clusters):
        color    = np.array(cmap_cl(ci / max(n_cl - 1, 1))[:3], dtype=np.float32)
        cl_cids  = cell_ids[labels == cl]
        valid_m  = (cl_cids > 0) & (cl_cids <= id_max)
        lut_rgb[cl_cids[valid_m]]   = color
        lut_alpha[cl_cids[valid_m]] = 0.68

    he_f     = he.astype(np.float32) / 255.0
    cell_rgb = lut_rgb[mask]
    cell_a   = lut_alpha[mask, None]
    fg       = (mask > 0)[:, :, None]
    blended  = np.where(fg, (1 - cell_a) * he_f + cell_a * cell_rgb, he_f)

    bounds = find_boundaries(mask, mode="thin")
    blended[bounds] = [0.15, 0.15, 0.15]

    comp_img = (blended * 255).clip(0, 255).astype(np.uint8)
    return he, comp_img, labels, clusters, X_log, gene_names, cmap_cl, selected, cl_gene_map


# ═════════════════════════════════════════════════════════════════════════════
# Fig. 2f — ROI9 Leiden composite map  (standalone, full width)
# ═════════════════════════════════════════════════════════════════════════════

def make_fig2f():
    print("=== fig2f: ROI9 Leiden composite map ===")
    he, comp_img, labels, clusters, X_log, gene_names, cmap_cl, selected, _ = _load_roi9_data()

    n_cl  = len(clusters)
    h_px, w_px = he.shape[:2]
    fig_w = 183 * MM
    fig_h = fig_w * h_px / w_px + 10 * MM

    fig, ax = plt.subplots(1, 1, figsize=(fig_w, fig_h),
                            gridspec_kw=dict(left=0.005, right=0.995,
                                             top=0.95, bottom=0.005))
    ax.imshow(comp_img, aspect="auto", interpolation="bilinear")
    ax.axis("off")
    ax.set_title("")

    leg_handles = [
        mpatches.Patch(color=cmap_cl(ci / max(n_cl - 1, 1)),
                       label=ROI9_CLUSTER_NAMES.get(int(cl), f"Cluster {cl}"))
        for ci, cl in enumerate(clusters)
    ]
    ax.legend(handles=leg_handles, fontsize=7, loc="lower right",
              framealpha=0.88, edgecolor="#ccc")
    ax.text(-0.005, 1.03, "f", transform=ax.transAxes,
            fontsize=12, fontweight="bold", va="top", ha="right")

    # Scale bar (bottom-left): 100 µm @ 0.2737 µm/px
    scale_px = 50 / 0.2737           # ≈ 183 px
    margin_x = w_px * 0.04
    margin_y = h_px * 0.04
    bar_y    = h_px - margin_y
    bar_x0   = margin_x
    bar_x1   = bar_x0 + scale_px
    ax.plot([bar_x0, bar_x1], [bar_y, bar_y],
            color="white", linewidth=4, solid_capstyle="butt", zorder=10)
    ax.text((bar_x0 + bar_x1) / 2, bar_y - h_px * 0.02,
            "50 µm",
            color="white", ha="center", va="bottom",
            fontsize=8, fontweight="bold", zorder=10)

    out = Path("/Volumes/SSD/plan_a/submission_bioinformatics/figures/fig2/fig2h.png")
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  ✅ Saved: {out}")


# ═════════════════════════════════════════════════════════════════════════════
# Fig. 2g — ROI9 cluster marker dotplot  (standalone, narrower)
# ═════════════════════════════════════════════════════════════════════════════

def make_fig2g():
    print("=== fig2g: ROI9 cluster dotplot ===")
    _, _, labels, clusters, X_log, gene_names, cmap_cl, selected, cl_gene_map = _load_roi9_data()

    n_cl      = len(clusters)
    cl_labels = [ROI9_CLUSTER_NAMES.get(int(c), f"Cluster {c}") for c in clusters]

    # ── Build curated gene list ordered by cluster ───────────────────────────
    curated_selected = []
    curated_cl_gene_map = {}   # cluster_int → [genes]
    for cl in clusters:
        genes = [g for g in ROI9_CURATED_GENES.get(int(cl), []) if g in gene_names]
        curated_cl_gene_map[int(cl)] = genes
        curated_selected.extend(genes)

    n_genes   = len(curated_selected)
    pct_arr, mean_arr = build_dotplot_arrays(X_log, labels, curated_selected, gene_names)
    vmax_g    = max(mean_arr.max(), 0.1)

    fig_w = 183 * MM
    fig_h = 112 * MM
    fig, ax = plt.subplots(1, 1, figsize=(fig_w, fig_h),
                            gridspec_kw=dict(left=0.17, right=0.76,
                                             top=0.96, bottom=0.34))

    render_dotplot(ax, fig, pct_arr, mean_arr,
                   x_labels=curated_selected, y_labels=cl_labels,
                   title="",          # no title — described in figure caption
                   vmax_mean=vmax_g, max_dot_s=420, legends_outside=True)

    # ── Override tick aesthetics ──────────────────────────────────────────────
    ax.set_xticklabels(curated_selected, fontsize=8, rotation=45, ha="right",
                       rotation_mode="anchor")
    ax.set_yticklabels(cl_labels[::-1], fontsize=8.5, fontweight="normal")
    ax.tick_params(axis="both", length=0)

    # ── Alternating row shading ───────────────────────────────────────────────
    for i in range(n_cl):
        if i % 2 == 0:
            ax.axhspan(i - 0.46, i + 0.46, color="#f5f5f5", zorder=0, lw=0)

    # ── Group separators between clusters ────────────────────────────────────
    ci_map = {int(cl): ci for ci, cl in enumerate(clusters)}
    gene_cl_map = {g: int(cl)
                   for cl in clusters
                   for g in curated_cl_gene_map.get(int(cl), [])}
    boundary = 0
    for cl in clusters[:-1]:
        boundary += len(curated_cl_gene_map.get(int(cl), []))
        ax.axvline(boundary - 0.5, color="#bbbbbb", lw=0.8, ls="-", alpha=0.6, zorder=1)

    # ── Cluster-coloured brackets above x-axis ────────────────────────────────
    boundary = 0
    for ci, cl in enumerate(clusters):
        n_g_cl = len(curated_cl_gene_map.get(int(cl), []))
        if n_g_cl == 0:
            continue
        x0, x1 = boundary - 0.35, boundary + n_g_cl - 0.65
        color = cmap_cl(ci / max(n_cl - 1, 1))
        ax.annotate("", xy=(x1, n_cl - 0.15), xytext=(x0, n_cl - 0.15),
                    xycoords="data", textcoords="data",
                    arrowprops=dict(arrowstyle="-", color=color, lw=2.5),
                    annotation_clip=False)
        boundary += n_g_cl

    # ── Tick label colours per cluster ────────────────────────────────────────
    for tick in ax.get_xticklabels():
        gene = tick.get_text()
        if gene in gene_cl_map:
            tick.set_color(cmap_cl(ci_map[gene_cl_map[gene]] / max(n_cl - 1, 1)))

    # ── Lighter grid ──────────────────────────────────────────────────────────
    ax.grid(alpha=0.10, lw=0.35, zorder=0)
    ax.set_axisbelow(True)

    # ── Panel label ──────────────────────────────────────────────────────────
    ax.text(-0.20, 1.03, "g", transform=ax.transAxes,
            fontsize=12, fontweight="bold", va="top", ha="right")

    out = OUT_DIR / "fig2g.png"
    fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  ✅ Saved: {out}")


# ── Main ──────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    print("Generating Fig. 2d, 2e, 2f, 2g ...")
    make_fig2d()
    make_fig2e()
    make_fig2f()
    make_fig2g()
    print("\nDone.")
    print("  fig2d.png  — ROI10 cell-type composite map (standalone)")
    print("  fig2e.png  — AT2/AT1 dotplot (standalone)")
    print("  fig2fg.png — ROI9 Leiden clusters + dotplot")
