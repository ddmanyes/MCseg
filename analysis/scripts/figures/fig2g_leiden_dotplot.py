"""
fig2g_leiden_dotplot.py
=======================
Generate fig2g: ROI10 Leiden cluster marker gene dotplot.

Replaces the winner-take-all 5-type dotplot with one derived from the same
Leiden clustering used in fig2d (res=0.4, k=15, seed=42).

Rows: 3 Leiden clusters
  C1 → "AT2 Pneumocyte"  (identified by 100% SFTPC co-expression)
  C0 → "Unresolved (C0)"
  C2 → "Unresolved (C2)"

Columns: marker genes spanning AT2, AT1, Macrophage, Endothelial lineages.
  Confirms C1 AT2 identity; shows absence of other lineage markers in C0/C2.

Output:
  submission_bioinformatics/figures/fig2/fig2g.png

Run:
  cd /Volumes/SSD/plan_a
  uv run python submission_bioinformatics/scripts/figures/fig2g_leiden_dotplot.py
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
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D
from pathlib import Path
from sklearn.decomposition import TruncatedSVD
from sklearn.neighbors import NearestNeighbors
import leidenalg
import igraph

H5AD_PATH = Path("/Volumes/SSD/plan_a/xenium_he_seg/results/visiumhd/visiumhd_cells/roi10_v12.h5ad")
OUT_PATH  = Path("/Volumes/SSD/plan_a/submission_bioinformatics/figures/fig2/fig2g.png")

# Marker genes to display — same panel as old dotplot + key AT1 genes
# Ordered: AT2 markers → AT1 markers → Macrophage → Endothelial
GENE_GROUPS = [
    ("AT2 Pneumocyte",      ["SFTPC", "SFTPB", "SFTPA1"]),
    ("AT1 Pneumocyte",      ["AGER", "RTKN2", "CAV1"]),
    ("Alveolar Macrophage", ["MARCO", "FABP4", "SPP1"]),
    ("Endothelial",         ["PECAM1", "VWF"]),
]
GROUP_COLOURS = {
    "AT2 Pneumocyte":      "#2166AC",
    "AT1 Pneumocyte":      "#B2182B",
    "Alveolar Macrophage": "#D4841A",
    "Endothelial":         "#1A7340",
}

plt.rcParams.update({
    "font.family": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 8, "axes.titlesize": 8.5, "pdf.fonttype": 42,
    "axes.linewidth": 0.7,
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


def render_dotplot(ax, fig, pct_arr: np.ndarray, mean_arr: np.ndarray,
                   x_labels: list[str], y_labels: list[str],
                   vmax_mean: float | None = None, max_dot_s: float = 380):
    n_gr, n_g = pct_arr.shape
    if vmax_mean is None:
        vmax_mean = max(float(mean_arr.max()), 0.05)
    norm = mcolors.Normalize(vmin=0, vmax=vmax_mean)
    cmap = plt.cm.YlOrRd

    xs, ys, ss, cs = [], [], [], []
    for gi in range(n_gr):
        y_pos = n_gr - 1 - gi
        for gei in range(n_g):
            xs.append(gei)
            ys.append(y_pos)
            ss.append(max((pct_arr[gi, gei] / 100) * max_dot_s, 2))
            cs.append(mean_arr[gi, gei])

    sc_obj = ax.scatter(xs, ys, s=ss, c=cs, cmap=cmap, norm=norm,
                        edgecolors="grey", linewidths=0.3, zorder=3)

    ax.set_xlim(-0.6, n_g - 0.4)
    ax.set_ylim(-0.6, n_gr - 0.4)
    ax.set_xticks(range(n_g))
    ax.set_xticklabels(x_labels, fontsize=7.5, rotation=45, ha="right",
                       rotation_mode="anchor")
    ax.set_yticks(range(n_gr))
    ax.set_yticklabels(y_labels[::-1], fontsize=8.5)
    ax.grid(alpha=0.12, lw=0.35, zorder=0)
    ax.set_axisbelow(True)
    ax.tick_params(axis="both", length=0)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Alternating row shading
    for i in range(n_gr):
        if i % 2 == 0:
            ax.axhspan(i - 0.46, i + 0.46, color="#f5f5f5", zorder=0, lw=0)

    # Colorbar (outside right)
    ax_pos = ax.get_position()
    cbar_ax = fig.add_axes([ax_pos.x1 + 0.013, ax_pos.y0 + ax_pos.height * 0.42,
                             0.014, ax_pos.height * 0.52])
    cb = fig.colorbar(sc_obj, cax=cbar_ax)
    cb.set_label("Mean log(1+UMI)", fontsize=6, labelpad=2)
    cb.ax.tick_params(labelsize=5.5, pad=1)
    cb.set_ticks(np.linspace(0, vmax_mean, 5))

    # Size legend
    leg_h = [
        Line2D([0], [0], marker="o", color="w",
               markerfacecolor="#aaa", markeredgecolor="grey",
               markeredgewidth=0.4,
               markersize=np.sqrt(max((p / 100) * max_dot_s, 2) / np.pi) * 2,
               label=f"{p}%")
        for p in [25, 50, 75]
    ]
    ax.legend(handles=leg_h, title="% positive",
              bbox_to_anchor=(1.02, 0.38), loc="upper left",
              fontsize=6, title_fontsize=6,
              handletextpad=0.3, labelspacing=0.3,
              framealpha=0.85, edgecolor="#ccc")

    return sc_obj


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    print("[1] Loading h5ad...")
    adata = ad.read_h5ad(str(H5AD_PATH))
    var_names = np.array(adata.var_names)
    print(f"  {adata.n_obs} cells, {adata.n_vars} genes")

    print("[2] Normalizing...")
    X_log = normalize_log1p(adata.X)

    print("[3] Leiden clustering (res=0.4, k=15)...")
    leiden_labels = leiden_cluster(X_log, resolution=0.4, k=15, seed=42)
    n_clusters = len(np.unique(leiden_labels))
    print(f"  {n_clusters} clusters")

    # Identify AT2 cluster by SFTPC mean expression
    gene_idx = {g: i for i, g in enumerate(var_names)}
    sftpc_mean = [float(X_log[leiden_labels == cl, gene_idx["SFTPC"]].mean())
                  if "SFTPC" in gene_idx else 0.0
                  for cl in range(n_clusters)]
    at2_cl = int(np.argmax(sftpc_mean))

    # Cluster sizes
    for cl in range(n_clusters):
        n = int((leiden_labels == cl).sum())
        sftpc_pct = float((X_log[leiden_labels == cl, gene_idx.get("SFTPC", 0)] > 0).mean()) * 100 \
                    if "SFTPC" in gene_idx else 0.0
        tag = " ← AT2" if cl == at2_cl else ""
        print(f"  C{cl}(n={n:4d}): SFTPC+={sftpc_pct:.1f}%{tag}")

    # Build cluster labels (AT2 cluster named, others as "Unresolved (Cx)")
    cluster_labels = []
    for cl in range(n_clusters):
        if cl == at2_cl:
            n = int((leiden_labels == cl).sum())
            cluster_labels.append(f"AT2 Pneumocyte  (C{cl}, n={n})")
        else:
            n = int((leiden_labels == cl).sum())
            cluster_labels.append(f"Unresolved  (C{cl}, n={n})")

    print("[4] Building dotplot arrays...")
    all_genes = [g for _, gs in GENE_GROUPS for g in gs if g in gene_idx]
    n_gr = n_clusters
    n_g  = len(all_genes)
    pct_arr  = np.zeros((n_gr, n_g))
    mean_arr = np.zeros((n_gr, n_g))

    for cl in range(n_clusters):
        mask = leiden_labels == cl
        for gi, gene in enumerate(all_genes):
            expr = X_log[mask, gene_idx[gene]]
            pct_arr[cl, gi]  = float((expr > 0).mean() * 100)
            mean_arr[cl, gi] = float(expr.mean())

    vmax = max(float(mean_arr.max()), 0.1)

    print("[5] Plotting...")
    fig_w = 145 * MM
    fig_h =  90 * MM
    fig, ax = plt.subplots(1, 1, figsize=(fig_w, fig_h),
                            gridspec_kw=dict(left=0.32, right=0.78,
                                             top=0.88, bottom=0.30))

    render_dotplot(ax, fig, pct_arr, mean_arr,
                   x_labels=all_genes, y_labels=cluster_labels,
                   vmax_mean=vmax, max_dot_s=380)

    # Colour x-tick labels by lineage group
    gene_to_group = {g: grp for grp, gs in GENE_GROUPS for g in gs}
    for tick in ax.get_xticklabels():
        gene = tick.get_text()
        col = GROUP_COLOURS.get(gene_to_group.get(gene, ""), "#333333")
        tick.set_color(col)

    # Group separator verticals
    boundary = 0
    for _, gs in GENE_GROUPS[:-1]:
        boundary += sum(1 for g in gs if g in gene_idx)
        ax.axvline(boundary - 0.5, color="#bbbbbb", lw=0.8, ls="-", alpha=0.55, zorder=1)

    # Coloured brackets above x-axis labelling gene groups
    boundary = 0
    for grp_name, gs in GENE_GROUPS:
        n_g_grp = sum(1 for g in gs if g in gene_idx)
        if n_g_grp == 0:
            continue
        x0 = boundary - 0.35
        x1 = boundary + n_g_grp - 0.65
        col = GROUP_COLOURS.get(grp_name, "#666666")
        ax.annotate("", xy=(x1, n_clusters - 0.12), xytext=(x0, n_clusters - 0.12),
                    xycoords="data", textcoords="data",
                    arrowprops=dict(arrowstyle="-", color=col, lw=2.2),
                    annotation_clip=False)
        ax.text((x0 + x1) / 2, n_clusters - 0.0,
                grp_name.replace(" Pneumocyte", "").replace("Alveolar ", "Alv."),
                ha="center", va="bottom", fontsize=5.5, color=col,
                fontweight="bold", transform=ax.transData)
        boundary += n_g_grp

    ax.set_title("")
    ax.text(-0.44, 1.07, "g", transform=ax.transAxes,
            fontsize=12, fontweight="bold", va="top", ha="right")

    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PATH, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"\n✅ Saved: {OUT_PATH}")


if __name__ == "__main__":
    main()
