"""
Step 5: SPP1+ macrophage nearest-neighbour analysis.

Completes the immunosuppressive triad:
  Treg <-> CD8    (already in 02_nearest_neighbor.py)
  Treg <-> SPP1+  (already in 02_nearest_neighbor.py)
  SPP1+ <-> CD8   ← NEW (missing edge)

Also generates a triangle summary figure showing all three cell-pair
distances side-by-side as the core "MCseg enables cell-contact analysis"
figure.

Outputs:
  results/spp1_nn_stats.csv
  figures/panelD_triad_triangle.png
  figures/panelD_violin.png
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from sklearn.neighbors import BallTree
from pathlib import Path

BASE    = Path('/Volumes/SSD/plan_a/submission_bioinformatics/analysis/treg_crc')
OUT     = BASE / 'results'
FIG_DIR = BASE / 'figures'
N_PERM     = 1000
CONTACT_UM = 30
SEED = 42

# ── Load data ─────────────────────────────────────────────────────────
df = pd.read_csv(OUT / 'highres_seg' / 'cell_centroids.csv')

tregs = df[df['subtype'] == 'Treg'][['row_um', 'col_um']].values
cd8   = df[df['subtype'] == 'CD8'][['row_um', 'col_um']].values
spp1  = df[df['subtype'] == 'SPP1_macro'][['row_um', 'col_um']].values
all_coords = df[['row_um', 'col_um']].values

print(f'Tregs: {len(tregs)}, CD8: {len(cd8)}, SPP1+: {len(spp1)}')

row_min, row_max = all_coords[:, 0].min(), all_coords[:, 0].max()
col_min, col_max = all_coords[:, 1].min(), all_coords[:, 1].max()

# ── Helper ─────────────────────────────────────────────────────────────
def nn_dist(query, target):
    tree = BallTree(target, metric='euclidean')
    d, _ = tree.query(query, k=1)
    return d.flatten()

def csr_permutation(query_n, target, n_perm, rng):
    medians = []
    for _ in range(n_perm):
        rand = np.column_stack([
            rng.uniform(row_min, row_max, query_n),
            rng.uniform(col_min, col_max, query_n),
        ])
        medians.append(np.median(nn_dist(rand, target)))
    return np.array(medians)

# ── SPP1+ as query ────────────────────────────────────────────────────
rng = np.random.default_rng(SEED)

obs_spp1_cd8  = nn_dist(spp1, cd8)
obs_spp1_treg = nn_dist(spp1, tregs)

print(f'\nSPP1+ -> CD8   median: {np.median(obs_spp1_cd8):.1f} um')
print(f'SPP1+ -> Treg  median: {np.median(obs_spp1_treg):.1f} um')

perm_spp1_cd8  = csr_permutation(len(spp1), cd8,   N_PERM, rng)
perm_spp1_treg = csr_permutation(len(spp1), tregs, N_PERM, rng)

p_spp1_cd8  = (perm_spp1_cd8  <= np.median(obs_spp1_cd8)).mean()
p_spp1_treg = (perm_spp1_treg <= np.median(obs_spp1_treg)).mean()

print(f'\nPermutation p-values (SPP1+ closer than random):')
print(f'  SPP1+ -> CD8:  p = {p_spp1_cd8:.4f}')
print(f'  SPP1+ -> Treg: p = {p_spp1_treg:.4f}')

# Contact: SPP1+ within CONTACT_UM of CD8
contact_spp1_cd8_obs = (obs_spp1_cd8 < CONTACT_UM).mean() * 100
rng2 = np.random.default_rng(SEED + 10)
contact_perm = np.array([
    (nn_dist(np.column_stack([rng2.uniform(row_min, row_max, len(spp1)),
                               rng2.uniform(col_min, col_max, len(spp1))]),
             cd8) < CONTACT_UM).mean() * 100
    for _ in range(500)
])
p_contact = (contact_perm >= contact_spp1_cd8_obs).mean()
print(f'\nSPP1+ within {CONTACT_UM}um of CD8: {contact_spp1_cd8_obs:.1f}% (CSR median: {np.median(contact_perm):.1f}%, p={p_contact:.4f})')

# ── Load existing Treg stats for triangle summary ─────────────────────
nn_stats = pd.read_csv(OUT / 'nn_stats.csv')
treg_cd8_row  = nn_stats[nn_stats['pair'] == 'Treg-CD8'].iloc[0]
treg_spp1_row = nn_stats[nn_stats['pair'] == 'Treg-SPP1+'].iloc[0]

# ── Save SPP1+ stats ───────────────────────────────────────────────────
spp1_stats = pd.DataFrame({
    'pair':           ['SPP1-CD8', 'SPP1-Treg'],
    'query':          ['SPP1+', 'SPP1+'],
    'target':         ['CD8', 'Treg'],
    'obs_median_um':  [np.median(obs_spp1_cd8), np.median(obs_spp1_treg)],
    'perm_median_um': [np.median(perm_spp1_cd8), np.median(perm_spp1_treg)],
    'p_value':        [p_spp1_cd8, p_spp1_treg],
})
spp1_stats.to_csv(OUT / 'spp1_nn_stats.csv', index=False)
print('\n', spp1_stats.to_string(index=False))

# ── Figure D1: violin (SPP1+ distances) ───────────────────────────────
rng3 = np.random.default_rng(SEED + 20)
rand_spp1 = np.column_stack([
    rng3.uniform(row_min, row_max, len(spp1)),
    rng3.uniform(col_min, col_max, len(spp1)),
])

records = []
for obs, rand, label in [
    (obs_spp1_cd8,  nn_dist(rand_spp1, cd8),   'SPP1+ → CD8'),
    (obs_spp1_treg, nn_dist(rand_spp1, tregs),  'SPP1+ → Treg'),
]:
    records += [{'pair': label, 'condition': 'Observed', 'distance_um': v} for v in obs]
    records += [{'pair': label, 'condition': 'Random',   'distance_um': v} for v in rand]
violin_df = pd.DataFrame(records)

palette = {'Observed': '#e67e22', 'Random': '#bdc3c7'}
fig, axes = plt.subplots(1, 2, figsize=(7.5, 4))
pairs_info = [
    ('SPP1+ → CD8',  p_spp1_cd8),
    ('SPP1+ → Treg', p_spp1_treg),
]
for ax, (pair, pv) in zip(axes, pairs_info):
    sub = violin_df[violin_df['pair'] == pair]
    sns.violinplot(data=sub, x='condition', y='distance_um',
                   hue='condition', palette=palette, inner='box',
                   linewidth=1.2, legend=False, cut=0, ax=ax)
    ax.set_title(pair, fontsize=11, fontweight='bold')
    ax.set_xlabel('')
    ax.set_ylabel('Distance (µm)' if pair == 'SPP1+ → CD8' else '')
    ax.set_ylim(bottom=0)
    label = f'p = {pv:.3f}' if pv >= 0.001 else 'p < 0.001'
    ax.text(0.5, 0.97, label, transform=ax.transAxes,
            ha='center', va='top', fontsize=9, color='#2c3e50')
    ax.spines[['top', 'right']].set_visible(False)

fig.suptitle('SPP1⁺ macrophage nearest-neighbour distances: observed vs. random (CSR)',
             fontsize=11, y=1.02)
plt.tight_layout()
plt.savefig(FIG_DIR / 'panelD_violin.png', dpi=150, bbox_inches='tight')
plt.close()
print(f'Saved -> {FIG_DIR}/panelD_violin.png')

# ── Figure D2: Triad triangle summary ─────────────────────────────────
# Collect all three edges
edges = [
    {
        'label':    'Treg ↔ CD8',
        'obs':      treg_cd8_row['obs_median_um'],
        'csr':      treg_cd8_row['perm_median_um'],
        'p':        treg_cd8_row['p_value'],
        'color':    '#c0392b',
    },
    {
        'label':    'SPP1⁺ ↔ CD8',
        'obs':      np.median(obs_spp1_cd8),
        'csr':      np.median(perm_spp1_cd8),
        'p':        p_spp1_cd8,
        'color':    '#e67e22',
    },
    {
        'label':    'SPP1⁺ ↔ Treg',
        'obs':      treg_spp1_row['obs_median_um'],
        'csr':      treg_spp1_row['perm_median_um'],
        'p':        treg_spp1_row['p_value'],
        'color':    '#8e44ad',
    },
]

fig, axes = plt.subplots(1, 3, figsize=(11, 4.2))

node_colors = {
    'Treg':  '#c0392b',
    'CD8':   '#2980b9',
    'SPP1⁺': '#e67e22',
}

for ax, edge in zip(axes, edges):
    cats = ['Observed', 'CSR baseline']
    vals = [edge['obs'], edge['csr']]
    colors = [edge['color'], '#bdc3c7']
    y_max = max(vals) * 1.40
    bars = ax.bar(cats, vals, color=colors, edgecolor='#2c3e50',
                  linewidth=1.2, width=0.5)
    ax.set_title(edge['label'], fontsize=11, fontweight='bold')
    ax.set_ylabel('Median distance (µm)')
    ax.set_ylim(0, y_max)
    ax.spines[['top', 'right']].set_visible(False)
    # p-value bracket (all in data coordinates)
    p_label = 'p < 0.001' if edge['p'] < 0.001 else f"p = {edge['p']:.3f}"
    bracket_y = max(vals) * 1.10
    ax.plot([0, 0, 1, 1],
            [max(vals)*1.03, bracket_y, bracket_y, max(vals)*1.03],
            color='#2c3e50', lw=1.2)
    ax.text(0.5, bracket_y + max(vals) * 0.04, p_label,
            ha='center', va='bottom', fontsize=9.5, color='#2c3e50')
    # value labels on bars
    for bar, val in zip(bars, vals):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + max(vals) * 0.02,
                f'{val:.0f} µm', ha='center', va='bottom', fontsize=8.5)

fig.suptitle('Immunosuppressive triad: SPP1⁺ TAM – Treg – CD8⁺ T cell\n'
             'Median nearest-neighbour distance vs. CSR random baseline',
             fontsize=11, y=1.02)
plt.tight_layout()
plt.savefig(FIG_DIR / 'panelD_triad_triangle.png', dpi=150, bbox_inches='tight')
plt.close()
print(f'Saved -> {FIG_DIR}/panelD_triad_triangle.png')

print('\nDone.')
