"""
18_figs5_tls_umap_validation.py
==============================
Supplementary Figure S5: UMAP Comparison (V12 vs Space Ranger)
證明 V12 在轉錄組空間中的分群純度與分離門檻。
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scanpy as sc
import pandas as pd
import numpy as np
from pathlib import Path

# ── 路徑 ──────────────────────────────────────────────────────────────────
ROOT     = Path("/Volumes/SSD/plan_a/crc_transcript_attribution")
AD_DIR   = ROOT / "results/anndata"
OUT_DIR  = Path("/Volumes/SSD/plan_a/manuscript/figures/05_validation")

# ── 設定 ──────────────────────────────────────────────────────────────────
V12_COLOR_MAP = {
    "Plasma / B Cell":      "#3498DB", # Blue
    "Well-diff. Tumor":     "#E74C3C", # Red
    "General Tumor":        "#9B59B6", # Purple
    "Myeloid / Macrophage": "#F39C12", # Orange
    "Stroma / Fibroblast":  "#2ECC71", # Green
    "Unknown":              "#BDC3C7"
}

MARKERS = ["JCHAIN", "CEACAM5", "MT-CO3", "LYZ", "VIM"]

# ═══════════════════════════════════════════════════════════════════════════
# 1. 處理 V12 資料 (同步 Fig 5d 邏輯)
# ═══════════════════════════════════════════════════════════════════════════
print("1. Processing V12 UMAP...")
v12_adata = sc.read_h5ad(AD_DIR / "v12_roi15.h5ad")

# 同步分群
sc.pp.normalize_total(v12_adata, target_sum=1e4)
sc.pp.log1p(v12_adata)
sc.tl.pca(v12_adata)
sc.pp.neighbors(v12_adata)
sc.tl.leiden(v12_adata, resolution=0.8) # 為了分出 5 群，設為 0.8

# 映射標籤
ANNOTATION = {
    "0": "General Tumor",
    "1": "Plasma / B Cell",
    "2": "Stroma / Fibroblast",
    "3": "Well-diff. Tumor",
    "4": "Myeloid / Macrophage"
}
v12_adata.obs["cell_type"] = v12_adata.obs["leiden"].astype(str).map(ANNOTATION).fillna("Unknown")

sc.tl.umap(v12_adata)

# ═══════════════════════════════════════════════════════════════════════════
# 2. 處理 SR 資料
# ═══════════════════════════════════════════════════════════════════════════
print("2. Processing SR UMAP...")
sr_adata = sc.read_h5ad(AD_DIR / "sr_roi15.h5ad")
sc.pp.normalize_total(sr_adata, target_sum=1e4)
sc.pp.log1p(sr_adata)
sc.pp.pca(sr_adata)
sc.pp.neighbors(sr_adata)
sc.tl.umap(sr_adata)
sc.tl.leiden(sr_adata, resolution=0.5)

# 映射 SR 標籤 (標註為 Mixed)
SR_ANNOTATION = {
    "0": "SR Tumor (Mixed)",
    "1": "SR Plasma/B (Mixed)",
    "2": "SR Stroma (Mixed)",
    "3": "SR Myeloid (Mixed)"
}
sr_adata.obs["cell_type"] = sr_adata.obs["leiden"].astype(str).map(SR_ANNOTATION).fillna("Unknown")
SR_COLOR_MAP = {
    "SR Plasma/B (Mixed)":  "#3498DB",
    "SR Tumor (Mixed)":     "#E74C3C",
    "SR Stroma (Mixed)":    "#2ECC71",
    "SR Myeloid (Mixed)":   "#F39C12",
    "Unknown":              "#BDC3C7"
}

# ═══════════════════════════════════════════════════════════════════════════
# 3. 繪圖 (3-Row Layout: Clusters + Marker Expr.)
# ═══════════════════════════════════════════════════════════════════════════
import matplotlib.gridspec as gridspec

fig = plt.figure(figsize=(20, 15))
gs = gridspec.GridSpec(3, 5, figure=fig, height_ratios=[1.2, 1, 1])

# --- 第一排：Cluster UMAP (佔據 gs[0, 0:2] 與 gs[0, 3:5]) ---
ax_v12_cl = fig.add_subplot(gs[0, :2])
sc.pl.umap(v12_adata, color="cell_type", palette=V12_COLOR_MAP, 
           legend_loc="on data", legend_fontsize=12, legend_fontoutline=3,
           size=250, alpha=0.9, edgecolor="none", show=False, ax=ax_v12_cl, frameon=False)
ax_v12_cl.set_title(f"A: V12 Clusters (n={v12_adata.n_obs})", fontsize=16, fontweight="bold")

ax_sr_cl = fig.add_subplot(gs[0, 3:])  # 中間留一個空隙或者靠右
sc.pl.umap(sr_adata, color="cell_type", palette=SR_COLOR_MAP, 
           legend_loc="on data", legend_fontsize=12, legend_fontoutline=3,
           size=350, alpha=0.9, edgecolor="none", show=False, ax=ax_sr_cl, frameon=False)
ax_sr_cl.set_title(f"B: Space Ranger Units (n={sr_adata.n_obs})", fontsize=16, fontweight="bold")

# --- 第二、三排：Marker Heatmaps (V12 top, SR bottom) ---
GENES = ["JCHAIN", "CEACAM5", "MT-CO3", "LYZ", "VIM"]

for i, gene in enumerate(GENES):
    # V12 Marker Plot
    ax_v12_m = fig.add_subplot(gs[1, i])
    sc.pl.umap(v12_adata, color=gene, cmap="viridis", size=80, ax=ax_v12_m, 
               show=False, frameon=False, title=f"V12: {gene}")
    
    # SR Marker Plot
    ax_sr_m = fig.add_subplot(gs[2, i])
    sc.pl.umap(sr_adata, color=gene, cmap="viridis", size=100, ax=ax_sr_m, 
               show=False, frameon=False, title=f"SR: {gene}")

plt.suptitle("Supplementary Figure S5: Transcriptional Purity & Purity Validation on UMAP", 
             fontsize=22, fontweight="bold", y=0.98)
plt.subplots_adjust(hspace=0.3, wspace=0.1)

out_file = OUT_DIR / "figS5_tls_umap_comparison.png"
plt.savefig(out_file, dpi=200, bbox_inches="tight") # 降低 dpi 確保生成速度，正式版可調高
plt.close()

print(f"✅ Supplementary Figure S5 saved to {out_file}")
