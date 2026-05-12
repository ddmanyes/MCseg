"""
Step 3: Treg/CD8 spatial density ratio map + Moran's I (Direction A).
Output: figures/panelA_density_morans.png, results/morans_results.csv
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

OUT     = Path('/Volumes/SSD/plan_a/submission_bioinformatics/analysis/treg_crc/results')
FIG_DIR = Path('/Volumes/SSD/plan_a/submission_bioinformatics/analysis/treg_crc/figures')

GRID_UM = 200
N_PERM  = 999
SEED    = 42

df = pd.read_csv(OUT / 'cell_centroids.csv')

row_min = df['row_um'].min(); row_max = df['row_um'].max()
col_min = df['col_um'].min(); col_max = df['col_um'].max()

row_edges = np.arange(row_min, row_max + GRID_UM, GRID_UM)
col_edges = np.arange(col_min, col_max + GRID_UM, GRID_UM)

tregs  = df[df['subtype'] == 'Treg']
cd8    = df[df['subtype'] == 'CD8']
t_cells = df[df['subtype'].isin(['Treg', 'CD8', 'CD4', 'Th17', 'Tfh', 'NK'])]

def bin_counts(sub):
    H, _, _ = np.histogram2d(sub['row_um'], sub['col_um'],
                              bins=[row_edges, col_edges])
    return H

H_treg  = bin_counts(tregs)
H_cd8   = bin_counts(cd8)
H_tcell = bin_counts(t_cells)

with np.errstate(invalid='ignore', divide='ignore'):
    ratio = np.where(H_tcell > 0, H_treg / (H_treg + H_cd8 + 1e-6), np.nan)

print(f'Grid: {GRID_UM}µm  shape: {H_treg.shape}')
print(f'Non-empty grids: {(H_tcell > 0).sum()}')
print(f'Ratio mean: {np.nanmean(ratio):.3f}  range: {np.nanmin(ratio):.3f}–{np.nanmax(ratio):.3f}')

# Moran's I (queen contiguity, manual permutation)
valid_idx = np.argwhere(~np.isnan(ratio))
vals = ratio[~np.isnan(ratio)]

def morans_i(v, idx):
    n = len(v)
    diffs = v - v.mean()
    W = np.zeros((n, n))
    for a, (r1, c1) in enumerate(idx):
        for b, (r2, c2) in enumerate(idx):
            if a != b and abs(r1 - r2) <= 1 and abs(c1 - c2) <= 1:
                W[a, b] = 1.0
    W_sum = W.sum()
    if W_sum == 0:
        return 0.0
    num = n * (diffs[:, None] * W * diffs[None, :]).sum()
    den = W_sum * (diffs ** 2).sum()
    return num / den if den != 0 else 0.0

print('Computing Moran\'s I...')
obs_I = morans_i(vals, valid_idx)
rng = np.random.default_rng(SEED)
perm_I = [morans_i(rng.permutation(vals), valid_idx) for _ in range(N_PERM)]
p_moran = (np.array(perm_I) >= obs_I).mean()

print(f"Global Moran's I = {obs_I:.4f},  permutation p = {p_moran:.4f}")

pd.DataFrame({'metric': ['Treg_CD8_ratio_morans_I'],
              'moran_i': [obs_I],
              'p_value': [p_moran]}).to_csv(OUT / 'morans_results.csv', index=False)

# ── Figure ─────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))

ax = axes[0]
im0 = ax.imshow(H_tcell.T, origin='lower', cmap='Blues',
                extent=[row_min, row_max, col_min, col_max], aspect='auto')
ax.set_title('T cell density\n(all T subtypes)', fontweight='bold')
ax.set_xlabel('Row (µm)'); ax.set_ylabel('Col (µm)')
plt.colorbar(im0, ax=ax, shrink=0.7, label='cells/grid')

ax = axes[1]
cmap = plt.cm.RdBu_r.copy()
cmap.set_bad('lightgrey')
vmax = np.nanpercentile(ratio, 95)
im1 = ax.imshow(np.ma.masked_invalid(ratio).T, origin='lower',
                cmap=cmap, vmin=0, vmax=vmax,
                extent=[row_min, row_max, col_min, col_max], aspect='auto')
ax.set_title(f'Treg:(Treg+CD8) ratio\n(grid={GRID_UM} µm)', fontweight='bold')
ax.set_xlabel('Row (µm)')
plt.colorbar(im1, ax=ax, shrink=0.7, label='Treg fraction')
ax.text(0.02, 0.97, f"Moran's I = {obs_I:.3f}\np = {p_moran:.3f}",
        transform=ax.transAxes, va='top', fontsize=8,
        bbox=dict(boxstyle='round', fc='white', alpha=0.8))

ax = axes[2]
ax.hist(perm_I, bins=40, color='#bdc3c7', edgecolor='white', linewidth=0.5)
ax.axvline(obs_I, color='#c0392b', linewidth=2, label=f'Observed I = {obs_I:.3f}')
ax.set_xlabel("Moran's I"); ax.set_ylabel('Permutation count')
ax.set_title("Moran's I permutation test\n(n=999)", fontweight='bold')
ax.legend(fontsize=9)
ax.spines[['top', 'right']].set_visible(False)

plt.tight_layout()
plt.savefig(FIG_DIR / 'panelA_density_morans.png', dpi=150, bbox_inches='tight')
plt.close()
print(f'Saved -> {FIG_DIR}/panelA_density_morans.png')
