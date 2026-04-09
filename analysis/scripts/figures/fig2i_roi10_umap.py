"""
35_fig2_roi10_umap.py
=====================
Cell clustering & 2D embedding for MCseg v2 ROI10 (n=1,085 cells, Visium HD).

Pipeline (numba-free, macOS safe):
  normalize_log1p → top HVG → TruncatedSVD (PCA)
  → sklearn KNN graph → igraph
  → leidenalg clustering
  → sklearn t-SNE (perplexity=30) for 2D embedding

Output:
  manuscript/figures/fig2/fig2_roi10_umap.png

Run:
  cd /Volumes/SSD/plan_a
  uv run python manuscript/scripts/35_fig2_roi10_umap.py
"""

from __future__ import annotations
import warnings
warnings.filterwarnings("ignore")
import os
os.environ["OMP_NUM_THREADS"] = "1"

import numpy as np
import scipy.sparse as sp
import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path
from sklearn.decomposition import TruncatedSVD
from sklearn.neighbors import NearestNeighbors
from sklearn.manifold import TSNE
from collections import Counter
import leidenalg, igraph

H5AD_PATH = Path("/Volumes/SSD/plan_a/xenium_he_seg/results/visiumhd/visiumhd_cells/roi10_v12.h5ad")
OUT_PATH  = Path("/Volumes/SSD/plan_a/manuscript/figures/fig2/fig2_roi10_umap.png")

GENE_SETS = {
    "AT2 Pneumocyte":      ["SFTPC", "SFTPB", "SFTPA1", "SFTPA2"],
    "Alveolar Macrophage": ["MARCO", "FABP4", "MCEMP1", "SPP1"],
    "Endothelial":         ["PECAM1", "VWF", "CLDN5"],
    "AT1 Pneumocyte":      ["AGER", "RTKN2"],
}
PALETTE = {
    "AT2 Pneumocyte":      "#2166AC",
    "Alveolar Macrophage": "#D4841A",
    "Endothelial":         "#1A7340",
    "AT1 Pneumocyte":      "#B2182B",
    "Unassigned":          "#CCCCCC",
}

plt.rcParams.update({
    "font.family": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 8, "axes.titlesize": 9, "pdf.fonttype": 42,
    "axes.linewidth": 0.6,
})
MM = 1 / 25.4


# ── Helpers ───────────────────────────────────────────────────────────────────

def normalize_log1p(X_sparse) -> np.ndarray:
    X = sp.csr_matrix(X_sparse)
    totals = np.array(X.sum(axis=1)).flatten()
    totals[totals == 0] = 1
    X = X.multiply(1e4 / totals[:, None]).tocsr()
    X.data = np.log1p(X.data)
    return X.toarray()


def hybrid_cell_type(X_log: np.ndarray, var_names: np.ndarray) -> np.ndarray:
    n = X_log.shape[0]
    gene_idx = {g: i for i, g in enumerate(var_names)}
    cell_type = np.full(n, "Unassigned", dtype=object)

    at1 = np.zeros(n, dtype=bool)
    for g in GENE_SETS["AT1 Pneumocyte"]:
        if g in gene_idx:
            at1 |= (X_log[:, gene_idx[g]] > 0)
    cell_type[at1] = "AT1 Pneumocyte"

    wta_types = ["AT2 Pneumocyte", "Alveolar Macrophage", "Endothelial"]
    scores = np.zeros((n, len(wta_types)), dtype=np.float32)
    for ci, t in enumerate(wta_types):
        for g in GENE_SETS[t]:
            if g in gene_idx:
                scores[:, ci] = np.maximum(scores[:, ci], X_log[:, gene_idx[g]])
    remaining = ~at1
    winner = scores.argmax(axis=1)
    max_s   = scores.max(axis=1)
    cell_type[remaining] = np.where(
        max_s[remaining] > 0,
        np.array(wta_types)[winner[remaining]],
        "Unassigned")
    return cell_type


def build_knn_graph(X_pca: np.ndarray, k: int = 15) -> igraph.Graph:
    nn = NearestNeighbors(n_neighbors=k + 1, n_jobs=1)
    nn.fit(X_pca)
    _, idxs = nn.kneighbors(X_pca)
    n = X_pca.shape[0]
    edges = [(int(i), int(j)) for i in range(n) for j in idxs[i, 1:]]
    g = igraph.Graph(n=n, edges=edges, directed=False)
    g.simplify()
    return g


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    print("[1] Loading h5ad...")
    adata    = ad.read_h5ad(str(H5AD_PATH))
    var_names = np.array(adata.var_names)
    X_log    = normalize_log1p(adata.X)
    print(f"  Shape: {adata.shape}")

    print("[2] Hybrid cell typing...")
    cell_type   = hybrid_cell_type(X_log, var_names)
    all_types   = list(GENE_SETS.keys()) + ["Unassigned"]
    type_counts = {t: int((cell_type == t).sum()) for t in all_types}
    print(f"  {type_counts}")

    print("[3] HVG + TruncatedSVD (PCA)...")
    n_top_hvg = min(2000, adata.n_vars)
    var_all   = (X_log ** 2).mean(0) - X_log.mean(0) ** 2
    top_idx   = np.argsort(var_all)[-n_top_hvg:]
    X_hvg     = X_log[:, top_idx]

    n_pcs = min(50, X_hvg.shape[1] - 1, adata.n_obs - 1)
    svd   = TruncatedSVD(n_components=n_pcs, random_state=42)
    X_pca = svd.fit_transform(X_hvg)
    print(f"  PCA: {X_hvg.shape} → {X_pca.shape}  "
          f"(expl. var: {svd.explained_variance_ratio_.sum():.1%})")

    print("[4] Building KNN graph (k=15)...")
    g = build_knn_graph(X_pca, k=15)
    print(f"  Graph: {g.vcount()} nodes, {g.ecount()} edges")

    print("[5] Leiden clustering (res=0.4)...")
    part = leidenalg.find_partition(
        g, leidenalg.RBConfigurationVertexPartition,
        resolution_parameter=0.4, seed=42)
    leiden_labels = np.array(part.membership)
    n_clusters    = len(np.unique(leiden_labels))
    print(f"  {n_clusters} clusters")

    print("[6] t-SNE 2D embedding (sklearn, numba-free)...")
    coords = TSNE(n_components=2, perplexity=30, random_state=42,
                  n_jobs=1).fit_transform(X_pca)
    print(f"  t-SNE done. coord range x:{coords[:,0].min():.1f}~{coords[:,0].max():.1f}")

    print("[7] Plotting...")
    marker_genes = ["SFTPC", "AGER", "RTKN2", "MARCO", "PECAM1", "LAMP3"]
    gene_idx_map = {g: i for i, g in enumerate(var_names)}

    fig_w = 183 * MM
    fig_h = 150 * MM
    fig, axes = plt.subplots(3, 3, figsize=(fig_w, fig_h),
                              gridspec_kw=dict(
                                  left=0.06, right=0.97,
                                  top=0.93, bottom=0.05,
                                  hspace=0.50, wspace=0.38))

    def style_ax(ax, title, xlabel="Dim 1", ylabel="Dim 2"):
        ax.set_title(title, fontsize=7.5, fontweight="bold")
        ax.set_xlabel(xlabel, fontsize=6)
        ax.set_ylabel(ylabel, fontsize=6)
        ax.tick_params(labelsize=5.5)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    # ── [0,0] Cell type ────────────────────────────────────────────────────
    ax = axes[0, 0]
    for ct in all_types[::-1]:
        mask  = cell_type == ct
        color = PALETTE[ct]
        s     = 6 if ct != "Unassigned" else 3
        alpha = 0.9 if ct != "Unassigned" else 0.25
        ax.scatter(coords[mask, 0], coords[mask, 1],
                   c=color, s=s, alpha=alpha, linewidths=0, rasterized=True)
    handles = [
        mpatches.Patch(color=PALETTE[t],
                       label=f"{t.replace(' Pneumocyte','').replace('Alveolar ','Alv.')}  n={type_counts[t]}")
        for t in all_types if type_counts.get(t, 0) > 0
    ]
    ax.legend(handles=handles, fontsize=4.5, loc="best",
              framealpha=0.85, edgecolor="#ccc",
              handlelength=1, borderpad=0.4, labelspacing=0.25)
    style_ax(ax, f"Cell type  (n={adata.n_obs:,})")

    # ── [0,1] Leiden clusters ──────────────────────────────────────────────
    ax = axes[0, 1]
    cmap_l = plt.cm.get_cmap("tab10", n_clusters)
    for ci in range(n_clusters):
        mask = leiden_labels == ci
        dominant = Counter(cell_type[mask]).most_common(1)[0][0]
        short = (dominant.replace(" Pneumocyte", "")
                         .replace("Alveolar ", "Alv.")
                         .replace("Unassigned", "—"))
        ax.scatter(coords[mask, 0], coords[mask, 1],
                   c=[cmap_l(ci)], s=5, alpha=0.85, linewidths=0,
                   label=f"C{ci} (n={mask.sum()}, {short})", rasterized=True)
    ax.legend(fontsize=4.5, loc="best", framealpha=0.85, edgecolor="#ccc",
              handlelength=1, borderpad=0.4, labelspacing=0.2)
    style_ax(ax, f"Leiden clusters (res=0.4, k={n_clusters})")

    # ── [0,2] AT1 highlight ────────────────────────────────────────────────
    ax = axes[0, 2]
    is_at1 = cell_type == "AT1 Pneumocyte"
    ax.scatter(coords[~is_at1, 0], coords[~is_at1, 1],
               c="#DDDDDD", s=3, alpha=0.3, linewidths=0,
               rasterized=True, label="Other")
    ax.scatter(coords[is_at1, 0], coords[is_at1, 1],
               c="#B2182B", s=20, alpha=0.95, linewidths=0.4,
               edgecolors="white", rasterized=True,
               label=f"AT1  n={is_at1.sum()}")
    ax.legend(fontsize=6, loc="best", framealpha=0.85, edgecolor="#ccc")
    style_ax(ax, "AT1 Pneumocyte highlight")

    # ── Panels [1-2]: marker genes ─────────────────────────────────────────
    for k_i, gene in enumerate(marker_genes):
        row = 1 + k_i // 3
        col = k_i % 3
        ax  = axes[row, col]

        if gene in gene_idx_map:
            expr    = X_log[:, gene_idx_map[gene]]
            pct_pos = (expr > 0).mean() * 100
            pos_vals = expr[expr > 0]
            vmax    = float(np.percentile(pos_vals, 95)) if len(pos_vals) > 0 else 1.0
        else:
            expr    = np.zeros(adata.n_obs)
            pct_pos = 0.0
            vmax    = 1.0

        sc_obj = ax.scatter(coords[:, 0], coords[:, 1],
                             c=expr, cmap="YlOrRd", vmin=0, vmax=vmax,
                             s=4, alpha=0.85, linewidths=0, rasterized=True)
        cb = plt.colorbar(sc_obj, ax=ax, fraction=0.04, pad=0.02)
        cb.ax.tick_params(labelsize=4.5)
        label = f"{gene}  ({pct_pos:.1f}% pos.)" if gene in gene_idx_map \
                else f"{gene} ✗ (not in panel)"
        style_ax(ax, label)

    fig.suptitle(
        "ROI10 — MCseg v2 Leiden clustering  "
        "(n=1,085 cells, Visium HD 2 µm/bin)\n"
        "Layout: t-SNE (perplexity=30) on 50-PC embedding",
        fontsize=8.5, fontweight="bold", y=0.975)

    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PATH, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"\n✅ Saved: {OUT_PATH}")

    # ── Summary ────────────────────────────────────────────────────────────
    print(f"\n── Cluster summary ──────────────────────────────")
    for ci in range(n_clusters):
        mask     = leiden_labels == ci
        dominant = Counter(cell_type[mask]).most_common(1)[0]
        print(f"  C{ci}: n={mask.sum():4d}, "
              f"dominant={dominant[0]} ({dominant[1]/mask.sum():.0%})")


if __name__ == "__main__":
    main()
