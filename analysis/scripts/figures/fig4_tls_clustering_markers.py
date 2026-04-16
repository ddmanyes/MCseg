"""
30_fig_candidate_e.py
======================
生成 Candidate E (ROI 15, top island) 的等效分析圖：
  - Fig_E_clustering: H&E + V12 cell types + SR clusters (仿 fig5d)
  - Fig_E_markers: V12 vs SR marker purity (仿 fig5e)

使用 KMeans 取代 leidenalg (macOS 掛起問題)。
Candidate E: x0=54800, y0=600, x1=55800, y1=1600 (roi15)
"""

from __future__ import annotations

import os
os.environ["OMP_NUM_THREADS"] = "1"   # prevent leiden macOS deadlock

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import Normalize
import scanpy as sc
from scipy.sparse import issparse
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA as skPCA
from skimage.segmentation import find_boundaries
from pathlib import Path

# ── paths ──────────────────────────────────────────────────────────────────
ROOT     = Path("/Volumes/SSD/plan_a/crc_transcript_attribution")
HE_DIR   = Path("/Volumes/SSD/plan_a/crc_he_seg/results/rois")
MASK_DIR = ROOT / "results/masks"
AD_DIR   = ROOT / "results/anndata"
OUT_DIR      = Path("/Volumes/SSD/plan_a/manuscript/figures/05_validation")
OUT_DIR_SUB  = Path("/Volumes/SSD/plan_a/submission_bioinformatics/figures/fig4")
OUT_DIR.mkdir(parents=True, exist_ok=True)
OUT_DIR_SUB.mkdir(parents=True, exist_ok=True)

PIXEL_SIZE_UM = 0.2738   # CRC VisiumHD H&E: 0.2738 µm/px
SCALE_UM      = 50

def add_scale_bar(ax, h: int, w: int) -> None:
    """White 50 µm scale bar at bottom-right corner."""
    scale_px = SCALE_UM / PIXEL_SIZE_UM   # ≈ 183 px
    margin_x = w * 0.05
    margin_y = h * 0.06
    bar_y  = h - margin_y
    bar_x1 = w - margin_x
    bar_x0 = bar_x1 - scale_px
    ax.plot([bar_x0, bar_x1], [bar_y, bar_y],
            color="white", linewidth=3, solid_capstyle="butt", zorder=10)
    ax.text((bar_x0 + bar_x1) / 2, bar_y - h * 0.025,
            f"{SCALE_UM} µm",
            color="white", ha="center", va="bottom",
            fontsize=9, fontweight="bold", zorder=10)

ROI = "roi15"

# ── marker genes (figs 5e) ─────────────────────────────────────────────────
MARKERS       = ["JCHAIN", "CEACAM5", "LYZ", "VIM", "MT-CO3", "CD79A"]
MARKER_CMAPS  = {"JCHAIN": "Blues", "CEACAM5": "Reds", "LYZ": "YlOrRd",
                 "VIM": "Greens", "MT-CO3": "Purples", "CD79A": "RdPu"}


# ── helpers ────────────────────────────────────────────────────────────────

def load_data():
    print("Loading data...")
    he_img   = plt.imread(HE_DIR / f"{ROI}_he.png")
    v12_mask = np.load(MASK_DIR / f"v12_{ROI}.npy")
    sr_mask  = np.load(MASK_DIR / f"sr_{ROI}.npy")
    v12_ad   = sc.read_h5ad(AD_DIR / f"v12_{ROI}.h5ad")
    sr_ad    = sc.read_h5ad(AD_DIR  / f"sr_{ROI}.h5ad")

    if he_img.max() <= 1.0:
        he_img = (he_img * 255).astype(np.uint8)

    print(f"  V12: {v12_ad.shape}, SR: {sr_ad.shape}")
    return he_img, v12_mask, sr_mask, v12_ad, sr_ad


def preprocess(adata):
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    return adata


def kmeans_cluster(adata, n_clusters=5, n_pca=30, seed=42):
    """KMeans clustering on PCA embedding (no leidenalg needed)."""
    X = adata.X.toarray() if issparse(adata.X) else adata.X
    n_pca = min(n_pca, X.shape[1] - 1, X.shape[0] - 1)
    pca = skPCA(n_components=n_pca, random_state=seed)
    emb = pca.fit_transform(X)
    km  = KMeans(n_clusters=n_clusters, random_state=seed, n_init=10)
    labels = km.fit_predict(emb)
    return labels.astype(str)


_MARKERS_MAP = {
    "Plasma / B Cell":      ["JCHAIN", "CD79A", "MS4A1", "IGKC"],
    "Well-diff. Tumor":     ["CEACAM5", "EPCAM", "KRT20"],
    "Myeloid / Macrophage": ["LYZ", "CD68", "CSF1R"],
    "Stroma / Fibroblast":  ["VIM", "ACTA2", "DCN"],
}


def assign_v12_celltypes(adata, cluster_labels):
    """Marker-score–based cell type assignment (normalised by gene mean + per-cell max)."""
    adata.obs["cluster"] = cluster_labels
    var_names = list(adata.var_names)
    X_all = adata.X.toarray() if issparse(adata.X) else np.array(adata.X)
    gene_mean = X_all.mean(axis=0) + 1e-6

    ct_map = {}
    for cl in np.unique(cluster_labels):
        mask_cl = adata.obs["cluster"] == cl
        X_sub   = X_all[mask_cl.values]
        scores  = {}
        for ct, genes in _MARKERS_MAP.items():
            avail_idx = [var_names.index(g) for g in genes if g in var_names]
            if avail_idx:
                normed = X_sub[:, avail_idx] / gene_mean[avail_idx]
                scores[ct] = float(normed.max(axis=1).mean())
            else:
                scores[ct] = 0.0
        ct_map[cl] = max(scores, key=scores.get)
        print(f"    cluster {cl}: {ct_map[cl]}  "
              f"scores={dict(sorted(scores.items(), key=lambda x: -x[1]))}")
    return ct_map


V12_COLORS = {
    "Plasma / B Cell":      "#3498DB",
    "Well-diff. Tumor":     "#E74C3C",
    "Myeloid / Macrophage": "#F39C12",
    "Stroma / Fibroblast":  "#2ECC71",
    "Unknown":              "#BDC3C7",
}
SR_COLORS  = ["#E91E63", "#2196F3", "#4CAF50", "#FFC107", "#9C27B0"]


def create_cluster_panel(ax, he, mask, adata, cell_type_col, color_dict, title):
    ax.imshow(he, alpha=0.35)
    overlay = np.zeros((*mask.shape, 4), dtype=np.float32)
    for cid, ct in adata.obs[cell_type_col].items():
        mid = int(str(cid).replace("cell_", ""))
        cell_mask = mask == mid
        hex_c = color_dict.get(ct, "#808080")
        r, g, b = [int(hex_c[i:i+2], 16) / 255.0 for i in (1, 3, 5)]
        overlay[cell_mask] = [r, g, b, 0.75]
    ax.imshow(overlay)
    bounds = find_boundaries(mask, mode="outer")
    edge = np.zeros((*mask.shape, 4))
    edge[bounds] = [0, 0, 0, 0.9]
    ax.imshow(edge)
    ax.set_title(title, fontweight="bold", fontsize=12)
    ax.axis("off")


def render_marker(ax, he, mask, adata, gene):
    ax.imshow(he, alpha=0.35)
    overlay = np.zeros((*mask.shape, 4), dtype=np.float32)
    if gene in adata.var_names:
        expr = adata[:, gene].X.toarray().flatten() if issparse(adata[:, gene].X) \
               else adata[:, gene].X.flatten()
        vmax = np.percentile(expr[expr > 0], 95) if (expr > 0).any() else 1.0
        norm = Normalize(vmin=0, vmax=vmax)
        cmap = plt.get_cmap(MARKER_CMAPS.get(gene, "Oranges"))
        for cid, val in zip(adata.obs.index, expr):
            if val <= 0:
                continue
            mid = int(str(cid).replace("cell_", ""))
            cell_m = mask == mid
            c = cmap(norm(val))
            overlay[cell_m] = [c[0], c[1], c[2], 0.8]
    ax.imshow(overlay)
    bounds = find_boundaries(mask, mode="outer")
    edge = np.zeros((*mask.shape, 4))
    edge[bounds] = [0, 0, 0, 0.8]
    ax.imshow(edge)
    ax.set_title(gene, fontsize=10, fontweight="bold")
    ax.axis("off")
    if gene in adata.var_names:
        return norm, plt.get_cmap(MARKER_CMAPS.get(gene, "Oranges"))
    return Normalize(0, 1), plt.get_cmap("Oranges")


# ══════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════

def main():
    he_img, v12_mask, sr_mask, v12_ad, sr_ad = load_data()

    # Resize H&E to mask size
    from skimage.transform import resize as sk_resize
    he_v12 = sk_resize(he_img, v12_mask.shape + (3,),
                        anti_aliasing=True, preserve_range=True).astype(np.uint8)
    he_sr  = sk_resize(he_img, sr_mask.shape + (3,),
                        anti_aliasing=True, preserve_range=True).astype(np.uint8)

    print("Preprocessing...")
    preprocess(v12_ad)
    preprocess(sr_ad)

    print("KMeans clustering (V12, k=4)...")
    labels_v12 = kmeans_cluster(v12_ad, n_clusters=4)
    ct_map = assign_v12_celltypes(v12_ad, labels_v12)
    v12_ad.obs["cell_type"] = pd.Series(labels_v12, index=v12_ad.obs.index).map(ct_map).fillna("Unknown")
    print("  V12 cell types:", v12_ad.obs["cell_type"].value_counts().to_dict())

    print("KMeans clustering (SR, k=3)...")
    labels_sr = kmeans_cluster(sr_ad, n_clusters=3)
    ct_map_sr = assign_v12_celltypes(sr_ad, labels_sr)
    sr_ad.obs["cell_type"] = pd.Series(labels_sr, index=sr_ad.obs.index).map(ct_map_sr).fillna("Unknown")
    print("  SR cell types:", sr_ad.obs["cell_type"].value_counts().to_dict())

    # ── Zoom region (auto-detected: bottom-right SR-only + acellular H&E area) ─
    # In SR mask coordinates: rows 749–1000, cols 622–905
    ZOOM = dict(r0=749, c0=622, r1=1000, c1=905)

    # ── Fig E_clustering (仿 fig5d) ─────────────────────────────────────────
    print("\nGenerating Fig E_clustering (3-panel)...")
    from matplotlib.patches import Rectangle

    fig, axes = plt.subplots(1, 3, figsize=(20, 7), gridspec_kw={"wspace": 0.05})

    axes[0].imshow(he_v12, alpha=0.95)
    axes[0].set_title("H&E", fontweight="bold", fontsize=12)
    axes[0].axis("off")
    add_scale_bar(axes[0], he_v12.shape[0], he_v12.shape[1])

    # V12 panel — colour by cell type (same palette as SR for direct comparison)
    axes[1].imshow(he_v12, alpha=0.70)
    overlay_v12 = np.zeros((*v12_mask.shape, 4), dtype=np.float32)
    for cid, ct in v12_ad.obs["cell_type"].items():
        mid   = int(str(cid).replace("cell_", ""))
        cm    = v12_mask == mid
        hex_c = V12_COLORS.get(ct, "#808080")
        r, g, b = [int(hex_c[i:i+2], 16) / 255.0 for i in (1, 3, 5)]
        overlay_v12[cm] = [r, g, b, 0.45]
    axes[1].imshow(overlay_v12)
    bounds_v12 = find_boundaries(v12_mask, mode="outer")
    edge_v12 = np.zeros((*v12_mask.shape, 4))
    edge_v12[bounds_v12] = [0, 0, 0, 0.70]
    axes[1].imshow(edge_v12)
    axes[1].set_title("MCseg", fontweight="bold", fontsize=12)
    axes[1].axis("off")
    add_scale_bar(axes[1], he_v12.shape[0], he_v12.shape[1])
    v12_ct_counts = v12_ad.obs["cell_type"].value_counts()
    v12_handles = [
        mpatches.Patch(color=V12_COLORS.get(ct, "#808080"), label=f"{ct} (n={n})")
        for ct, n in v12_ct_counts.items()
    ]
    axes[1].legend(handles=v12_handles, loc="lower center",
                   bbox_to_anchor=(0.5, -0.12), fontsize=7.5, ncol=2)

    # ── Red rectangle on V12 panel (rescale ZOOM coords to v12 mask space) ──
    scale_r = v12_mask.shape[0] / sr_mask.shape[0]
    scale_c = v12_mask.shape[1] / sr_mask.shape[1]
    zr0 = int(ZOOM["r0"] * scale_r); zr1 = int(ZOOM["r1"] * scale_r)
    zc0 = int(ZOOM["c0"] * scale_c); zc1 = int(ZOOM["c1"] * scale_c)

    rect_v12 = Rectangle((zc0, zr0), zc1 - zc0, zr1 - zr0,
                          linewidth=1.5, edgecolor="red", facecolor="none")
    axes[1].add_patch(rect_v12)

    # SR panel — colour by cell type (same palette as V12 for direct comparison)
    axes[2].imshow(he_sr, alpha=0.70)
    overlay_sr = np.zeros((*sr_mask.shape, 4), dtype=np.float32)
    for cid, ct in sr_ad.obs["cell_type"].items():
        mid = int(str(cid).replace("cell_", ""))
        cm  = sr_mask == mid
        hex_c = V12_COLORS.get(ct, "#808080")
        r, g, b = [int(hex_c[i:i+2], 16) / 255.0 for i in (1, 3, 5)]
        overlay_sr[cm] = [r, g, b, 0.45]
    axes[2].imshow(overlay_sr)
    bounds_sr = find_boundaries(sr_mask, mode="outer")
    edge_sr = np.zeros((*sr_mask.shape, 4))
    edge_sr[bounds_sr] = [0, 0, 0, 0.70]
    axes[2].imshow(edge_sr)
    axes[2].set_title("SR", fontweight="bold", fontsize=12)
    axes[2].axis("off")
    add_scale_bar(axes[2], he_sr.shape[0], he_sr.shape[1])
    sr_ct_counts = sr_ad.obs["cell_type"].value_counts()
    sr_handles = [
        mpatches.Patch(color=V12_COLORS.get(ct, "#808080"),
                       label=f"{ct} (n={n})")
        for ct, n in sr_ct_counts.items()
    ]
    axes[2].legend(handles=sr_handles, loc="lower center",
                   bbox_to_anchor=(0.5, -0.12), fontsize=7.5, ncol=2)

    # Red rectangle on SR main panel (same ZOOM coords, already in SR space)
    rect_sr = Rectangle((ZOOM["c0"], ZOOM["r0"]),
                         ZOOM["c1"] - ZOOM["c0"], ZOOM["r1"] - ZOOM["r0"],
                         linewidth=1.5, edgecolor="red", facecolor="none")
    axes[2].add_patch(rect_sr)

    out_d = OUT_DIR / "fig4a_tls_clustering.png"
    plt.savefig(out_d, dpi=350, bbox_inches="tight")
    out_sub = OUT_DIR_SUB / "fig4a.png"
    plt.savefig(out_sub, dpi=350, bbox_inches="tight")
    plt.close(fig)
    print(f"✓ Saved: {out_d}")
    print(f"✓ Saved: {out_sub}")

    # ── Fig E_markers (仿 fig5e) ─────────────────────────────────────────────
    from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec

    n_m = len(MARKERS)
    print(f"\nGenerating Fig E_markers (2×{n_m})...")
    fig2 = plt.figure(figsize=(4 * n_m, 13))
    # Outer grid: 2 row-blocks (MCseg v2 / SR), each with a label strip + panels
    outer = GridSpec(2, 1, figure=fig2, hspace=0.30,
                     top=0.97, bottom=0.08, left=0.02, right=0.98)

    inner0 = GridSpecFromSubplotSpec(2, n_m, subplot_spec=outer[0],
                                     wspace=0.05, hspace=0.04,
                                     height_ratios=[0.07, 1])
    ax_label0 = fig2.add_subplot(inner0[0, :])
    ax_label0.axis("off")
    ax_label0.text(0.5, 0.4, "MCseg v2", ha="center", va="center",
                   fontsize=13, fontweight="bold", transform=ax_label0.transAxes)
    axes_v12 = [fig2.add_subplot(inner0[1, c]) for c in range(n_m)]

    inner1 = GridSpecFromSubplotSpec(2, n_m, subplot_spec=outer[1],
                                     wspace=0.05, hspace=0.04,
                                     height_ratios=[0.07, 1])
    ax_label1 = fig2.add_subplot(inner1[0, :])
    ax_label1.axis("off")
    ax_label1.text(0.5, 0.4, "SR", ha="center", va="center",
                   fontsize=13, fontweight="bold", transform=ax_label1.transAxes)
    axes_sr = [fig2.add_subplot(inner1[1, c]) for c in range(n_m)]

    v12_norms = []
    for i, g in enumerate(MARKERS):
        norm, cmap = render_marker(axes_v12[i], he_v12, v12_mask, v12_ad, g)
        v12_norms.append((norm, cmap))

    for i, g in enumerate(MARKERS):
        render_marker(axes_sr[i], he_sr, sr_mask, sr_ad, g)
        norm, cmap = v12_norms[i]
        cax = fig2.add_axes([
            axes_sr[i].get_position().x0 + 0.005,
            axes_sr[i].get_position().y0 - 0.025,
            axes_sr[i].get_position().width - 0.01,
            0.012
        ])
        cb = fig2.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap),
                           cax=cax, orientation="horizontal")
        cb.set_label("Expr. Level", fontsize=8)
        cb.ax.tick_params(labelsize=7)

    out_e = OUT_DIR / "fig4b_tls_markers.png"
    plt.savefig(out_e, dpi=350, bbox_inches="tight")
    plt.close(fig2)
    print(f"✓ Saved: {out_e}")

    print("\n✅ Done.")


if __name__ == "__main__":
    main()
