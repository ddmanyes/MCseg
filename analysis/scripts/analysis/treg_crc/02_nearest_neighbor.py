"""
Step 2: Treg nearest-neighbour distance analysis (Direction B).
- Treg vs CD8, SPP1+ macro, Tumor
- CSR permutation baseline (n=1000)
- Output: results/nn_distances.csv, figures/panelB*.png
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.neighbors import BallTree
from pathlib import Path

BASE    = Path('/Volumes/SSD/plan_a/submission_bioinformatics/analysis/treg_crc')
OUT     = BASE / 'results'
FIG_DIR = BASE / 'figures'
N_PERM     = 1000
CONTACT_UM = 30
SEED = 42

df = pd.read_csv(OUT / 'highres_seg' / 'cell_centroids.csv')

tregs = df[df['subtype'] == 'Treg'][['row_um', 'col_um']].values
cd8   = df[df['subtype'] == 'CD8'][['row_um', 'col_um']].values
spp1  = df[df['subtype'] == 'SPP1_macro'][['row_um', 'col_um']].values
tumor = df[df['subtype'] == 'Tumor'][['row_um', 'col_um']].values
all_coords = df[['row_um', 'col_um']].values

print(f'Tregs: {len(tregs)}, CD8: {len(cd8)}, SPP1+: {len(spp1)}, Tumor: {len(tumor)}')

def nn_dist(query, target):
    tree = BallTree(target, metric='euclidean')
    d, _ = tree.query(query, k=1)
    return d.flatten()

obs_cd8   = nn_dist(tregs, cd8)
obs_spp1  = nn_dist(tregs, spp1)
obs_tumor = nn_dist(tregs, tumor)

print(f'\nObserved Treg->CD8   median: {np.median(obs_cd8):.1f} um')
print(f'Observed Treg->SPP1+ median: {np.median(obs_spp1):.1f} um')
print(f'Observed Treg->Tumor median: {np.median(obs_tumor):.1f} um')
print(f'\nContact Tregs (<{CONTACT_UM}um from CD8): '
      f'{(obs_cd8 < CONTACT_UM).sum()} / {len(tregs)} '
      f'({(obs_cd8 < CONTACT_UM).mean()*100:.1f}%)')

# CSR permutation
rng = np.random.default_rng(SEED)
row_min, row_max = all_coords[:, 0].min(), all_coords[:, 0].max()
col_min, col_max = all_coords[:, 1].min(), all_coords[:, 1].max()

perm_cd8, perm_spp1, perm_tumor = [], [], []
for _ in range(N_PERM):
    r = np.column_stack([rng.uniform(row_min, row_max, len(tregs)),
                         rng.uniform(col_min, col_max, len(tregs))])
    perm_cd8.append(np.median(nn_dist(r, cd8)))
    perm_spp1.append(np.median(nn_dist(r, spp1)))
    perm_tumor.append(np.median(nn_dist(r, tumor)))

p_cd8   = (np.array(perm_cd8)   <= np.median(obs_cd8)).mean()
p_spp1  = (np.array(perm_spp1)  <= np.median(obs_spp1)).mean()
p_tumor = (np.array(perm_tumor) <= np.median(obs_tumor)).mean()

print(f'\nPermutation p-values (one-sided, closer than random):')
print(f'  Treg->CD8:   p = {p_cd8:.4f}')
print(f'  Treg->SPP1+: p = {p_spp1:.4f}')
print(f'  Treg->Tumor: p = {p_tumor:.4f}')

# Save distances (observed + one random replicate for violin)
rng2 = np.random.default_rng(SEED + 1)
rand_pts = np.column_stack([rng2.uniform(row_min, row_max, len(tregs)),
                             rng2.uniform(col_min, col_max, len(tregs))])
records = []
for d, pair in [(obs_cd8, 'Treg-CD8'), (obs_spp1, 'Treg-SPP1+'), (obs_tumor, 'Treg-Tumor')]:
    records += [{'pair': pair, 'condition': 'Observed', 'distance_um': v} for v in d]
for d, pair in [(nn_dist(rand_pts, cd8), 'Treg-CD8'),
                (nn_dist(rand_pts, spp1), 'Treg-SPP1+'),
                (nn_dist(rand_pts, tumor), 'Treg-Tumor')]:
    records += [{'pair': pair, 'condition': 'Random', 'distance_um': v} for v in d]

dist_df = pd.DataFrame(records)
dist_df.to_csv(OUT / 'nn_distances.csv', index=False)

stats = pd.DataFrame({
    'pair':          ['Treg-CD8', 'Treg-SPP1+', 'Treg-Tumor'],
    'obs_median_um': [np.median(obs_cd8), np.median(obs_spp1), np.median(obs_tumor)],
    'perm_median_um':[np.median(perm_cd8), np.median(perm_spp1), np.median(perm_tumor)],
    'p_value':       [p_cd8, p_spp1, p_tumor],
    'contact_pct':   [(obs_cd8 < CONTACT_UM).mean() * 100, None, None],
})
stats.to_csv(OUT / 'nn_stats.csv', index=False)
print('\n', stats.to_string(index=False))

# ── Figure B1: violin ──────────────────────────────────────────────
fig, axes = plt.subplots(1, 3, figsize=(11, 4))
palette = {'Observed': '#c0392b', 'Random': '#bdc3c7'}
pairs   = ['Treg-CD8', 'Treg-SPP1+', 'Treg-Tumor']
pvals   = [p_cd8, p_spp1, p_tumor]

for ax, pair, pv in zip(axes, pairs, pvals):
    sub = dist_df[dist_df['pair'] == pair]
    sns.violinplot(data=sub, x='condition', y='distance_um',
                   hue='condition', palette=palette, inner='box',
                   linewidth=1.2, legend=False, cut=0, ax=ax)
    ax.set_title(pair.replace('-', ' → '), fontsize=11, fontweight='bold')
    ax.set_xlabel('')
    ax.set_ylabel('Distance (um)' if pair == 'Treg-CD8' else '')
    ax.set_ylim(bottom=0)
    label = f'p = {pv:.3f}' if pv >= 0.001 else 'p < 0.001'
    ax.text(0.5, 0.97, label, transform=ax.transAxes,
            ha='center', va='top', fontsize=9, color='#2c3e50')
    ax.spines[['top', 'right']].set_visible(False)

fig.suptitle('Treg nearest-neighbour distances: observed vs. random (CSR)',
             fontsize=12, y=1.02)
plt.tight_layout()
plt.savefig(FIG_DIR / 'panelB1_nn_violin.png', dpi=150, bbox_inches='tight')
plt.close()
print(f'Saved -> {FIG_DIR}/panelB1_nn_violin.png')

# ── Figure B2: contact Treg bar ────────────────────────────────────
rng3 = np.random.default_rng(SEED + 2)
contact_perm = np.array([
    (nn_dist(np.column_stack([rng3.uniform(row_min, row_max, len(tregs)),
                               rng3.uniform(col_min, col_max, len(tregs))]),
             cd8) < CONTACT_UM).mean() * 100
    for _ in range(500)
])
contact_obs = (obs_cd8 < CONTACT_UM).mean() * 100
p_contact = (contact_perm >= contact_obs).mean()

fig, ax = plt.subplots(figsize=(4, 4))
ax.bar(['Observed', 'Random\n(CSR)'],
       [contact_obs, np.median(contact_perm)],
       color=['#c0392b', '#bdc3c7'], edgecolor='#2c3e50', linewidth=1.2, width=0.5)
ax.errorbar(1, np.median(contact_perm),
            yerr=[[np.median(contact_perm) - np.percentile(contact_perm, 2.5)],
                  [np.percentile(contact_perm, 97.5) - np.median(contact_perm)]],
            fmt='none', color='#2c3e50', capsize=5, linewidth=1.5)
ax.set_ylabel(f'Tregs within {CONTACT_UM} um of CD8 (%)')
ax.set_title('Contact-range Tregs', fontweight='bold')
ax.spines[['top', 'right']].set_visible(False)
p_contact_label = f'p = {p_contact:.3f}' if p_contact >= 0.001 else 'p < 0.001'
ax.text(0.5, 0.95, p_contact_label, transform=ax.transAxes,
        ha='center', va='top', fontsize=10)
plt.tight_layout()
plt.savefig(FIG_DIR / 'panelB2_contact_bar.png', dpi=150, bbox_inches='tight')
plt.close()
print(f'Saved -> {FIG_DIR}/panelB2_contact_bar.png')
