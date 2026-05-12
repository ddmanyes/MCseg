import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from sklearn.neighbors import BallTree
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# ── Paths ──────────────────────────────────────────────────────────────────
BASE    = Path('/Volumes/SSD/plan_a/submission_bioinformatics/analysis/treg_crc')
OUT     = BASE / 'results'
FIG_DIR = BASE / 'figures'
H5_PATH = Path('/Volumes/SSD/plan_a/tissue sample/ENACT_supporting_files/public_data/human_colorectal/input_files/filtered_feature_bc_matrix.h5')

# ── Load cell data ─────────────────────────────────────────────────────────
print("Loading cell centroids and labels...")
cells = pd.read_csv(OUT / 'highres_seg' / 'cell_centroids.csv')
labels = pd.read_csv(OUT / 'highres_seg' / 'celltypist_labels.csv')
cells = cells.merge(labels[['cell_id', 'celltypist_label', 'conf_score']], on='cell_id')

tregs = cells[cells['subtype'] == 'Treg'].reset_index(drop=True)
cd8_all = cells[cells['subtype'] == 'CD8'].reset_index(drop=True)

# Find nearest Treg for each CD8
print(f"Calculating distances from {len(cd8_all)} CD8+ T cells to {len(tregs)} Tregs...")
tree_treg = BallTree(tregs[['row_um', 'col_um']].values, metric='euclidean')
dist, _ = tree_treg.query(cd8_all[['row_um', 'col_um']].values, k=1)
cd8_all['dist_to_treg'] = dist.flatten()

# ── Extract RNA for CD8 cells ──────────────────────────────────────────────
print("Loading RNA data for CD8 cells...")
# To speed up, we rebuild the aggregation only for CD8 cell_ids
adata_full = sc.read_10x_h5(str(H5_PATH))
adata_full.var_names_make_unique()

# We need the bin attribution to know which barcodes belong to which CD8 cell
# Instead of reloading everything, we'll approximate by fetching the key genes
# But for accuracy, let's use the attribution logic if possible.
# Actually, the quickest way is to use the existing 'cell_centroids.csv' logic 
# but we need the counts. Let's do a mini-aggregation.

# Need attribution (re-run logic from 00_resegment)
mask = np.load(OUT / 'highres_seg' / 'mcseg_mask.npy')
tp = pd.read_parquet('/Volumes/SSD/plan_a/tissue sample/ENACT_supporting_files/public_data/human_colorectal/input_files/tissue_positions.parquet')
tp = tp[tp['in_tissue'] == 1]
# ... (standard attribution logic) ...
CROP_X0, CROP_Y0 = 5154, 4635
COL_OFFSET = 40598
row_l = (tp['pxl_row_in_fullres'].values - CROP_Y0).astype(np.int32)
col_l = (tp['pxl_col_in_fullres'].values - COL_OFFSET - CROP_X0).astype(np.int32)
H, W = mask.shape
valid_bins = (row_l >= 0) & (row_l < H) & (col_l >= 0) & (col_l < W)
tp['cell_id'] = 0
tp.loc[valid_bins, 'cell_id'] = mask[row_l[valid_bins], col_l[valid_bins]]

cd8_ids = set(cd8_all['cell_id'].values)
tp_cd8 = tp[tp['cell_id'].isin(cd8_ids)]
adata_cd8_bins = adata_full[adata_full.obs_names.isin(tp_cd8['barcode'])].copy()
adata_cd8_bins.obs['cell_id'] = tp_cd8.set_index('barcode')['cell_id'].reindex(adata_cd8_bins.obs_names).values

# Aggregate counts per CD8 cell
unique_cd8 = adata_cd8_bins.obs['cell_id'].unique()
cell_id_to_idx = {int(c): i for i, c in enumerate(unique_cd8)}
rows = np.array([cell_id_to_idx[int(c)] for c in adata_cd8_bins.obs['cell_id']])
cols = np.arange(len(adata_cd8_bins))
A = sp.csr_matrix((np.ones(len(cols)), (rows, cols)), shape=(len(unique_cd8), len(adata_cd8_bins)))
X_agg = A @ adata_cd8_bins.X

adata_cd8 = sc.AnnData(X=X_agg, var=adata_cd8_bins.var)
adata_cd8.obs_names = [str(c) for c in unique_cd8]
sc.pp.normalize_total(adata_cd8, target_sum=1e4)
sc.pp.log1p(adata_cd8)

# ── Merge expression back to cd8_all ───────────────────────────────────────
TARGET_GENES = ['GZMB', 'IFNG', 'PDCD1', 'LAG3', 'PRF1', 'CD8A']
for gene in TARGET_GENES:
    if gene in adata_cd8.var_names:
        expr = adata_cd8[:, gene].X.toarray().flatten()
        expr_df = pd.DataFrame({'cell_id': unique_cd8.astype(int), f'expr_{gene}': expr})
        cd8_all = cd8_all.merge(expr_df, on='cell_id', how='left')

# ── Analysis ───────────────────────────────────────────────────────────────
cd8_all['group'] = 'Far (>100um)'
cd8_all.loc[cd8_all['dist_to_treg'] < 100, 'group'] = 'Proximal (30-100um)'
cd8_all.loc[cd8_all['dist_to_treg'] < 30,  'group'] = 'Contact (<30um)'
group_order = ['Contact (<30um)', 'Proximal (30-100um)', 'Far (>100um)']

print("\n--- Mean Expression by Distance to Treg ---")
print(cd8_all.groupby('group')[[f'expr_{g}' for g in TARGET_GENES if f'expr_{g}' in cd8_all.columns]].mean())

# ── Plotting ───────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 3, figsize=(15, 5))
plot_genes = ['GZMB', 'PDCD1', 'IFNG']
palette = {'Contact (<30um)': '#e74c3c', 'Proximal (30-100um)': '#f39c12', 'Far (>100um)': '#3498db'}

for ax, gene in zip(axes, plot_genes):
    col = f'expr_{gene}'
    if col not in cd8_all.columns: continue
    sns.barplot(data=cd8_all, x='group', y=col, order=group_order, palette=palette, ax=ax, capsize=.1)
    ax.set_title(f'CD8 Expression: {gene}', fontweight='bold')
    ax.set_ylabel('log1p(CPM)')
    ax.set_xlabel('')
    ax.tick_params(axis='x', rotation=15)

plt.tight_layout()
plt.savefig(FIG_DIR / 'SuppFig_S13_treg_cd8_exhaustion.png', dpi=150, bbox_inches='tight')
print(f"\nSaved -> {FIG_DIR / 'SuppFig_S13_treg_cd8_exhaustion.png'}")
