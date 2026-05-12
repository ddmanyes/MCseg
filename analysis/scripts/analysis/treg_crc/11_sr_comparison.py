"""
Step 11: SR vs MCseg Treg spatial analysis comparison.

Runs the same Treg nearest-neighbour analysis on SR segmentation to test
whether MCseg's single-cell precision enables the analysis more sensitively.

Pipeline:
  1. Bin attribution using sr_mask.npy (same coordinate system as mcseg_mask)
  2. Aggregate counts per SR cell
  3. CellTypist annotation (Human_Colorectal_Cancer.pkl)
  4. Treg/CD8/SPP1+ nearest-neighbour analysis (same params as 02/05)
  5. Compare key metrics against MCseg results

Outputs:
  results/sr_celltypist_labels.csv
  results/sr_cell_centroids.csv
  results/sr_nn_stats.csv
  figures/panelSR_comparison.png
"""
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.neighbors import BallTree
from scipy.ndimage import find_objects
from pathlib import Path

BASE    = Path('/Volumes/SSD/plan_a/submission_bioinformatics/analysis/treg_crc')
OUT     = BASE / 'results'
FIG_DIR = BASE / 'figures'

SR_MASK = Path('/Volumes/SSD/plan_a/submission_bioinformatics/results/enact_method_comparison/sr_mask.npy')
H5_PATH = Path('/Volumes/SSD/plan_a/tissue sample/ENACT_supporting_files/public_data/human_colorectal/input_files/filtered_feature_bc_matrix.h5')
TP_PATH = Path('/Volumes/SSD/plan_a/tissue sample/ENACT_supporting_files/public_data/human_colorectal/input_files/tissue_positions.parquet')

CROP_X0    = 5154
CROP_Y0    = 4635
COL_OFFSET = 40598
UM_PER_PX  = 0.5   # sr_mask is in the same 0.5 µm/px space as mcseg_mask

N_PERM     = 1000
CONTACT_UM = 30
SEED       = 42

SUBTYPE_MAP = {
    'Regulatory T cells':        'Treg',
    'CD8+ T cells':              'CD8',
    'CD4+ T cells':              'CD4',
    'T helper 17 cells':         'Th17',
    'T follicular helper cells': 'Tfh',
    'NK cells':                  'NK',
    'Pro-inflammatory':          'ProInf_macro',
    'SPP1+':                     'SPP1_macro',
    'IgG+ Plasma':               'Plasma',
    'CMS1': 'Tumor', 'CMS2': 'Tumor', 'CMS3': 'Tumor',
    'Myofibroblasts': 'Stromal', 'Stromal 1': 'Stromal', 'Stromal 3': 'Stromal',
}

LABEL_MAP = {
    'CMS1': 'epithelial cells', 'CMS2': 'epithelial cells',
    'CMS3': 'epithelial cells', 'CMS4': 'epithelial cells',
    'Goblet cells': 'epithelial cells',
    'Mature Enterocytes type 1': 'epithelial cells',
    'Mature Enterocytes type 2': 'epithelial cells',
    'Stem-like/TA': 'epithelial cells', 'Intermediate': 'epithelial cells',
    'Proliferating': 'epithelial cells',
    'Myofibroblasts': 'stromal cells', 'Pericytes': 'stromal cells',
    'Smooth muscle cells': 'stromal cells', 'Stromal 1': 'stromal cells',
    'Stromal 2': 'stromal cells', 'Stromal 3': 'stromal cells',
    'Lymphatic ECs': 'stromal cells', 'Proliferative ECs': 'stromal cells',
    'Stalk-like ECs': 'stromal cells', 'Enteric glial cells': 'stromal cells',
    'CD19+CD20+ B': 'immune cells', 'CD4+ T cells': 'immune cells',
    'CD8+ T cells': 'immune cells', 'Regulatory T cells': 'immune cells',
    'NK cells': 'immune cells', 'IgA+ Plasma': 'immune cells',
    'IgG+ Plasma': 'immune cells', 'Mast cells': 'immune cells',
    'Pro-inflammatory': 'immune cells', 'SPP1+': 'immune cells',
}

# ── Step 1: SR mask centroids ──────────────────────────────────────────────
print('Loading SR mask...')
mask = np.load(str(SR_MASK))
print(f'  shape: {mask.shape}, cells: {int(mask.max()):,}')

ct_path        = OUT / 'sr_celltypist_labels.csv'
centroids_path = OUT / 'sr_cell_centroids.csv'

if not centroids_path.exists():
    print('Computing SR cell centroids...')
    slices = find_objects(mask)
    rows, cols, ids = [], [], []
    for i, sl in enumerate(slices):
        if sl is None:
            continue
        cid = i + 1
        sub = mask[sl]
        r, c = np.where(sub == cid)
        if len(r) == 0:
            continue
        rows.append(r.mean() + sl[0].start)
        cols.append(c.mean() + sl[1].start)
        ids.append(cid)
    centroids = pd.DataFrame({'cell_id': ids, 'row_px': rows, 'col_px': cols})
    centroids['row_um'] = centroids['row_px'] * UM_PER_PX
    centroids['col_um'] = centroids['col_px'] * UM_PER_PX
    print(f'  {len(centroids):,} centroids computed')
else:
    centroids = pd.read_csv(centroids_path)
    print(f'  Loaded {len(centroids):,} existing centroids')

# ── Step 2: Bin attribution + CellTypist ──────────────────────────────────
if not ct_path.exists():
    print('Bin attribution...')
    tp = pd.read_parquet(str(TP_PATH), columns=[
        'barcode', 'in_tissue', 'pxl_row_in_fullres', 'pxl_col_in_fullres'
    ])
    tp = tp[tp['in_tissue'] == 1]

    row_local = tp['pxl_row_in_fullres'].values - CROP_Y0
    col_local = tp['pxl_col_in_fullres'].values - COL_OFFSET - CROP_X0
    scale = 0.2738 / UM_PER_PX
    row_idx = (row_local * scale).astype(np.int32)
    col_idx = (col_local * scale).astype(np.int32)

    H, W = mask.shape
    valid = (row_idx >= 0) & (row_idx < H) & (col_idx >= 0) & (col_idx < W)
    tp = tp.copy()
    tp['cell_id'] = 0
    tp.loc[valid, 'cell_id'] = mask[row_idx[valid], col_idx[valid]]
    attribution = tp[tp['cell_id'] > 0][['barcode', 'cell_id']].reset_index(drop=True)
    print(f'  Attributed bins: {len(attribution):,}')

    print('Aggregating counts...')
    adata_full = sc.read_10x_h5(str(H5_PATH))
    adata_full.var_names_make_unique()
    adata_crop = adata_full[adata_full.obs_names.isin(attribution['barcode'].values)].copy()
    del adata_full

    barcode_to_cell = attribution.set_index('barcode')['cell_id']
    adata_crop.obs['cell_id'] = barcode_to_cell.reindex(adata_crop.obs_names).values
    cell_ids_v = adata_crop.obs['cell_id'].values.astype(np.int32)
    valid_mask = cell_ids_v > 0
    adata_valid = adata_crop[valid_mask]
    cell_ids_v  = cell_ids_v[valid_mask]
    unique_cells = np.unique(cell_ids_v)
    n_cells = len(unique_cells)
    print(f'  Unique SR cells with RNA: {n_cells:,}')

    cell_id_to_idx = {int(c): i for i, c in enumerate(unique_cells)}
    r_idx = np.array([cell_id_to_idx[int(c)] for c in cell_ids_v])
    c_idx = np.arange(len(cell_ids_v))
    A = sp.csr_matrix(
        (np.ones(len(cell_ids_v), dtype=np.float32), (r_idx, c_idx)),
        shape=(n_cells, adata_valid.n_obs),
    )
    X_agg = A @ adata_valid.X
    adata_cells = sc.AnnData(
        X=X_agg.tocsr() if sp.issparse(X_agg) else sp.csr_matrix(X_agg),
        var=adata_valid.var.copy(),
    )
    adata_cells.obs_names = [str(c) for c in unique_cells]
    del adata_crop, adata_valid, A

    sc.pp.normalize_total(adata_cells, target_sum=1e4)
    sc.pp.log1p(adata_cells)

    print('Running CellTypist...')
    import celltypist
    predictions = celltypist.annotate(
        adata_cells, model='Human_Colorectal_Cancer.pkl', majority_voting=False,
    )
    ct_labels   = predictions.predicted_labels['predicted_labels'].values
    conf_scores = predictions.probability_matrix.max(axis=1).values

    ct_df = pd.DataFrame({
        'cell_id':          unique_cells,
        'celltypist_label': ct_labels,
        'conf_score':       conf_scores,
        'broad_label':      [LABEL_MAP.get(lbl, 'other') for lbl in ct_labels],
    })
    ct_df.to_csv(ct_path, index=False)
    print(f'  Saved: {ct_path.name}')
    print(ct_df['broad_label'].value_counts().to_string())
else:
    print('  Loaded existing CellTypist labels')
    ct_df = pd.read_csv(ct_path)

# ── Merge + subtype ────────────────────────────────────────────────────────
df = centroids.merge(ct_df, on='cell_id', how='left')
df['subtype'] = df['celltypist_label'].map(SUBTYPE_MAP).fillna(df['broad_label'])
df.to_csv(centroids_path, index=False)

print('\nSR subtype counts:')
for s in ['Treg', 'CD8', 'SPP1_macro', 'Tumor']:
    print(f'  {s}: {(df.subtype == s).sum()}')

# ── Step 3: Nearest-neighbour analysis ────────────────────────────────────
print('\nNearest-neighbour analysis...')
tregs      = df[df['subtype'] == 'Treg'][['row_um', 'col_um']].values
cd8        = df[df['subtype'] == 'CD8'][['row_um', 'col_um']].values
all_coords = df[['row_um', 'col_um']].values

if len(tregs) == 0 or len(cd8) == 0:
    print('ERROR: no Tregs or CD8 found')
    raise SystemExit(1)

def nn_dist(query, target):
    tree = BallTree(target, metric='euclidean')
    d, _ = tree.query(query, k=1)
    return d.flatten()

obs_cd8 = nn_dist(tregs, cd8)
print(f'SR Treg->CD8 median: {np.median(obs_cd8):.1f} µm')
print(f'SR contact Tregs (<{CONTACT_UM}µm): {(obs_cd8 < CONTACT_UM).mean()*100:.1f}%')

rng = np.random.default_rng(SEED)
row_min, row_max = all_coords[:, 0].min(), all_coords[:, 0].max()
col_min, col_max = all_coords[:, 1].min(), all_coords[:, 1].max()
perm_cd8 = [
    np.median(nn_dist(np.column_stack([
        rng.uniform(row_min, row_max, len(tregs)),
        rng.uniform(col_min, col_max, len(tregs)),
    ]), cd8))
    for _ in range(N_PERM)
]
p_cd8       = (np.array(perm_cd8) <= np.median(obs_cd8)).mean()
contact_pct = (obs_cd8 < CONTACT_UM).mean() * 100

sr_stats = pd.DataFrame({
    'pair':          ['Treg-CD8'],
    'obs_median_um': [np.median(obs_cd8)],
    'perm_median_um':[np.median(perm_cd8)],
    'p_value':       [p_cd8],
    'contact_pct':   [contact_pct],
    'n_treg':        [len(tregs)],
    'n_cd8':         [len(cd8)],
})
sr_stats.to_csv(OUT / 'sr_nn_stats.csv', index=False)
print(sr_stats.to_string(index=False))

# ── Step 4: Comparison figure ──────────────────────────────────────────────
mc_stats  = pd.read_csv(OUT / 'nn_stats.csv')
mc_row    = mc_stats[mc_stats['pair'] == 'Treg-CD8'].iloc[0]
mc_ct     = pd.read_csv(OUT / 'highres_seg' / 'cell_centroids.csv')
mc_treg_n = (mc_ct.subtype == 'Treg').sum()
mc_cd8_n  = (mc_ct.subtype == 'CD8').sum()

fig, axes = plt.subplots(1, 2, figsize=(9, 4.5))

# Panel left: median distance
methods  = ['MCseg', 'SR']
obs_vals = [mc_row['obs_median_um'], np.median(obs_cd8)]
csr_vals = [mc_row['perm_median_um'], np.median(perm_cd8)]
x = np.arange(2)
w = 0.35
ax = axes[0]
b1 = ax.bar(x - w/2, obs_vals, w, label='Observed',
            color=['#c0392b', '#e67e22'], edgecolor='#2c3e50', linewidth=1.1)
b2 = ax.bar(x + w/2, csr_vals, w, label='CSR baseline',
            color='#bdc3c7', edgecolor='#2c3e50', linewidth=1.1)
for bar, val in zip(list(b1) + list(b2), obs_vals + csr_vals):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.8,
            f'{val:.0f}µm', ha='center', va='bottom', fontsize=8.5)
ax.set_xticks(x)
ax.set_xticklabels(methods, fontsize=11)
ax.set_ylabel('Treg → CD8 median distance (µm)')
ax.set_title('Nearest-neighbour distance\nTreg → CD8', fontweight='bold')
ax.legend(fontsize=8.5)
ax.spines[['top', 'right']].set_visible(False)
for xi, pv in zip(x, [mc_row['p_value'], p_cd8]):
    label = 'p < 0.001' if pv < 0.001 else f'p = {pv:.3f}'
    ax.text(xi, max(csr_vals) * 1.12, label, ha='center', fontsize=8, color='#2c3e50')

# Panel right: contact %
ax2 = axes[1]
contact_vals = [mc_row['contact_pct'], contact_pct]
bars = ax2.bar(methods, contact_vals,
               color=['#c0392b', '#e67e22'], edgecolor='#2c3e50',
               linewidth=1.1, width=0.45)
for bar, val in zip(bars, contact_vals):
    ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
             f'{val:.1f}%', ha='center', va='bottom', fontsize=10, fontweight='bold')
ax2.set_ylabel(f'Tregs within {CONTACT_UM}µm of CD8 (%)')
ax2.set_title(f'Contact-range Tregs\n(<{CONTACT_UM}µm from CD8)', fontweight='bold')
ax2.set_ylim(0, max(contact_vals) * 1.3)
ax2.spines[['top', 'right']].set_visible(False)
for xi, (nt, nc) in enumerate([(mc_treg_n, mc_cd8_n), (len(tregs), len(cd8))]):
    ax2.text(xi, -max(contact_vals) * 0.12,
             f'n={nt} Tregs\nn={nc} CD8', ha='center', fontsize=7.5, color='#7f8c8d')

fig.suptitle('MCseg vs Space Ranger — Treg spatial analysis',
             fontsize=12, fontweight='bold', y=1.01)
plt.tight_layout()
plt.savefig(FIG_DIR / 'panelSR_comparison.png', dpi=150, bbox_inches='tight')
plt.close()
print(f'\nSaved -> {FIG_DIR}/panelSR_comparison.png')
print('Done.')
