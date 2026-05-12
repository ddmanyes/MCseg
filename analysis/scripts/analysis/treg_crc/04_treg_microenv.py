"""
Step 4: Treg microenvironment composition analysis (Route 3).

For each Treg, count neighbour cell types within a 100um radius.
Split Tregs into contact (< 30um from CD8) vs non-contact groups.
Compare neighbourhood compositions between groups.

Outputs:
  results/treg_microenv.csv        — per-Treg neighbourhood counts
  results/microenv_stats.csv       — Mann-Whitney U + BH-FDR per cell type
  figures/panelC1_microenv_heatmap.png
  figures/panelC2_microenv_boxplot.png
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.neighbors import BallTree
from scipy import stats
from pathlib import Path

RADIUS_UM   = 100
CONTACT_UM  = 30
MIN_CELLS   = 5   # minimum cells with > 0 to include a cell type in stats

BASE    = Path('/Volumes/SSD/plan_a/submission_bioinformatics/analysis/treg_crc')
OUT     = BASE / 'results'
FIG_DIR = BASE / 'figures'

# ── Load data ──────────────────────────────────────────────────────────────
df = pd.read_csv(OUT / 'highres_seg' / 'cell_centroids.csv')

tregs = df[df['subtype'] == 'Treg'].reset_index(drop=True)
cd8   = df[df['subtype'] == 'CD8'][['row_um', 'col_um']].values
all_cells = df[['row_um', 'col_um']].values
all_subtypes = df['subtype'].values

print(f"Tregs: {len(tregs)}, CD8: {len(cd8)}, All cells: {len(df)}")

# ── Identify contact Tregs (< 30um from nearest CD8) ──────────────────────
tree_cd8 = BallTree(cd8, metric='euclidean')
d_cd8, _ = tree_cd8.query(tregs[['row_um', 'col_um']].values, k=1)
d_cd8 = d_cd8.flatten()
tregs = tregs.copy()
tregs['d_cd8_um'] = d_cd8
tregs['contact'] = d_cd8 < CONTACT_UM

n_contact = tregs['contact'].sum()
n_noncontact = (~tregs['contact']).sum()
print(f"Contact Tregs (<{CONTACT_UM}um from CD8): {n_contact} ({n_contact/len(tregs)*100:.1f}%)")
print(f"Non-contact Tregs:                        {n_noncontact}")

# ── Define analysis cell types (exclude very small groups) ────────────────
all_sub_counts = df['subtype'].value_counts()
analysis_subtypes = [s for s in all_sub_counts.index
                     if s != 'Treg' and all_sub_counts[s] >= 50]
print(f"\nAnalysis subtypes ({len(analysis_subtypes)}): {analysis_subtypes}")

# ── Build BallTree for all cells, query 100um neighbourhood ───────────────
tree_all = BallTree(all_cells, metric='euclidean')
treg_coords = tregs[['row_um', 'col_um']].values
neighbours = tree_all.query_radius(treg_coords, r=RADIUS_UM)

# ── Count neighbour cell types per Treg ───────────────────────────────────
records = []
for i, (treg_row, nbrs) in enumerate(zip(tregs.itertuples(), neighbours)):
    nbr_subtypes = all_subtypes[nbrs]
    nbr_subtypes = nbr_subtypes[nbr_subtypes != 'Treg']  # exclude Tregs from neighbourhood

    row = {
        'cell_id':   treg_row.cell_id,
        'row_um':    treg_row.row_um,
        'col_um':    treg_row.col_um,
        'd_cd8_um':  treg_row.d_cd8_um,
        'contact':   treg_row.contact,
        'n_total':   len(nbr_subtypes),
    }
    for sub in analysis_subtypes:
        count = (nbr_subtypes == sub).sum()
        row[f'n_{sub}'] = count
        row[f'frac_{sub}'] = count / max(len(nbr_subtypes), 1)
    records.append(row)

microenv = pd.DataFrame(records)
microenv.to_csv(OUT / 'treg_microenv.csv', index=False)
print(f"\nSaved microenv table: {len(microenv)} rows")
print(f"Median neighbourhood size: {microenv['n_total'].median():.0f} cells per Treg")

# ── Statistical comparison: contact vs non-contact ────────────────────────
contact_df    = microenv[microenv['contact']]
noncontact_df = microenv[~microenv['contact']]

stat_rows = []
for sub in analysis_subtypes:
    col = f'frac_{sub}'
    a = contact_df[col].values
    b = noncontact_df[col].values
    if (a > 0).sum() < MIN_CELLS and (b > 0).sum() < MIN_CELLS:
        continue
    u, p = stats.mannwhitneyu(a, b, alternative='two-sided')
    stat_rows.append({
        'subtype':           sub,
        'contact_median':    np.median(a),
        'noncontact_median': np.median(b),
        'contact_mean':      np.mean(a),
        'noncontact_mean':   np.mean(b),
        'U':                 u,
        'p_raw':             p,
        'n_contact':         len(a),
        'n_noncontact':      len(b),
    })

stat_df = pd.DataFrame(stat_rows)

# BH-FDR correction
from scipy.stats import rankdata
n_tests = len(stat_df)
ranks = rankdata(stat_df['p_raw'])
stat_df['p_fdr'] = np.minimum(stat_df['p_raw'] * n_tests / ranks, 1.0)
stat_df = stat_df.sort_values('p_raw').reset_index(drop=True)

stat_df.to_csv(OUT / 'microenv_stats.csv', index=False)
print("\n=== Microenvironment statistics (contact vs non-contact Tregs) ===")
print(stat_df[['subtype','contact_median','noncontact_median','p_raw','p_fdr']].to_string(index=False))

sig = stat_df[stat_df['p_fdr'] < 0.05]
print(f"\nSignificant after FDR < 0.05: {len(sig)} cell types")
if len(sig):
    print(sig[['subtype','contact_median','noncontact_median','p_fdr']].to_string(index=False))

# ── Figure C1: Heatmap of mean neighbourhood fractions ────────────────────
heat_data = pd.DataFrame({
    'Contact':     [contact_df[f'frac_{s}'].mean() for s in analysis_subtypes],
    'Non-contact': [noncontact_df[f'frac_{s}'].mean() for s in analysis_subtypes],
}, index=analysis_subtypes)

heat_data['diff'] = heat_data['Contact'] - heat_data['Non-contact']
heat_data = heat_data.sort_values('diff', ascending=False).drop(columns='diff')

label_map = {
    'Tumor': 'Tumor', 'CD8': 'CD8⁺ T', 'CD4': 'CD4⁺ T',
    'SPP1_macro': 'SPP1⁺ Macro', 'ProInf_macro': 'Pro-inflam Macro',
    'Plasma': 'Plasma cell', 'NK': 'NK', 'Th17': 'Th17',
    'Tfh': 'Tfh', 'Stromal': 'Stromal', 'stromal cells': 'Stromal (alt)',
    'other': 'Other immune', 'immune cells': 'Immune (misc)',
    'epithelial cells': 'Epithelial',
}
heat_data.index = [label_map.get(s, s) for s in heat_data.index]

fig, ax = plt.subplots(figsize=(4, 6))
sns.heatmap(
    heat_data * 100,
    annot=True, fmt='.1f', cmap='RdBu_r',
    center=0, linewidths=0.5, linecolor='#cccccc',
    cbar_kws={'label': 'Mean % of 100um neighbourhood'},
    ax=ax
)
ax.set_title('Treg 100um neighbourhood\ncomposition (%)', fontsize=11, fontweight='bold')
ax.set_xlabel('')
ax.set_ylabel('')
plt.tight_layout()
plt.savefig(FIG_DIR / 'panelC1_microenv_heatmap.png', dpi=150, bbox_inches='tight')
plt.close()
print(f'\nSaved -> {FIG_DIR}/panelC1_microenv_heatmap.png')

# ── Figure C2: Boxplot for top-4 most-different cell types ────────────────
stat_df2 = stat_df.copy()
stat_df2['abs_diff'] = abs(stat_df2['contact_median'] - stat_df2['noncontact_median'])
top_types = stat_df2.nlargest(4, 'abs_diff')['subtype'].tolist()

frac_cols = [f'frac_{s}' for s in top_types]
plot_df = microenv[['contact'] + frac_cols].copy()
plot_df['group'] = plot_df['contact'].map({True: 'Contact\nTreg', False: 'Non-contact\nTreg'})

plot_long = plot_df.melt(
    id_vars='group', value_vars=frac_cols,
    var_name='cell_type', value_name='fraction'
)
plot_long['cell_type'] = plot_long['cell_type'].str.replace('frac_', '', regex=False)
plot_long['cell_type'] = plot_long['cell_type'].map(lambda x: label_map.get(x, x))
plot_long['fraction'] *= 100

fig, axes = plt.subplots(1, len(top_types), figsize=(3 * len(top_types), 4), sharey=False)
if len(top_types) == 1:
    axes = [axes]

palette = {'Contact\nTreg': '#c0392b', 'Non-contact\nTreg': '#2980b9'}

for ax, sub in zip(axes, top_types):
    label = label_map.get(sub, sub)
    sub_data = plot_long[plot_long['cell_type'] == label]
    sns.boxplot(data=sub_data, x='group', y='fraction',
                hue='group', palette=palette, width=0.5,
                linewidth=1.2, fliersize=3, legend=False, ax=ax)
    ax.set_title(label, fontsize=10, fontweight='bold')
    ax.set_xlabel('')
    ax.set_ylabel('% of neighbourhood' if sub == top_types[0] else '')
    ax.spines[['top', 'right']].set_visible(False)

    row = stat_df[stat_df['subtype'] == sub]
    if len(row):
        p = row['p_fdr'].values[0]
        p_label = ('FDR < 0.001' if p < 0.001
                   else f'FDR = {p:.3f}' if p < 0.05
                   else f'FDR = {p:.2f}')
        ax.text(0.5, 0.97, p_label, transform=ax.transAxes,
                ha='center', va='top', fontsize=8.5,
                color='#c0392b' if p < 0.05 else '#7f8c8d')

fig.suptitle('Cell-type composition around contact vs non-contact Tregs (100um)',
             fontsize=11, y=1.02)
plt.tight_layout()
plt.savefig(FIG_DIR / 'panelC2_microenv_boxplot.png', dpi=150, bbox_inches='tight')
plt.close()
print(f'Saved -> {FIG_DIR}/panelC2_microenv_boxplot.png')
