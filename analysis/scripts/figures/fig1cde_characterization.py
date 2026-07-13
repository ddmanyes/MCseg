#!/usr/bin/env python3
"""
fig1cde_characterization.py
Figure 1 c-e：MCseg v2 特性刻畫三聯圖
  c. Evolutionary Discovery Path（AutoResearch 演化曲線）
  d. Parameter Sensitivity（Diameter × Expansion 熱圖）
  e. Component Ablation Study

字級針對 Keynote 合併大圖（本圖僅佔約半頁寬）放大，確保縮排後仍可辨識。
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import seaborn as sns
from pathlib import Path

sns.set_theme(style="white", context="paper")
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']

fig = plt.figure(figsize=(17, 7.5), dpi=300)

# --- Panel c: Evolution Path (LHS) ---
ax1 = plt.subplot2grid((2, 7), (0, 0), rowspan=2, colspan=3)

stages = ["Naive\nBaseline", "Proseg\nBest", "Stage 2\nLEAP", "Stage 3\nEnsemble",
          "Stage 4\nVoronoi", "Stage 6\ncpsam", "MCseg v2\nFinal"]
scores = [0.3176, 0.3830, 0.6483, 0.6333, 0.6194, 0.6443, 0.6499]
x = np.arange(len(stages))

ax1.axvspan(-0.5, 1.5, color='#f0f0f0', alpha=0.6, zorder=0)
ax1.axvspan(1.5, 6.5, color='#e6f3ff', alpha=0.5, zorder=0)
ax1.text(0.5, 0.762, "Phase 1\nDiscovery", ha='center', va='top', fontsize=14, color='#888888')
ax1.text(4.0, 0.762, "Phase 2: Refinement", ha='center', va='top', fontsize=14, color='#5588aa')

ax1.plot(x, scores, marker='o', markersize=9, color='#004488', lw=2.5, zorder=3)
ax1.fill_between(x, scores, alpha=0.1, color='#004488')

for i, (xi, s) in enumerate(zip(x, scores)):
    offset = -0.042 if i == 3 else 0.030
    va = 'bottom' if offset > 0 else 'top'
    ax1.text(xi, s + offset, f"{s:.3f}", ha='center', va=va, fontsize=13, color='#004488')

ax1.annotate('Breakthrough\n(+0.33)', xy=(2, 0.648), xytext=(3.1, 0.44),
             arrowprops=dict(facecolor='#CC3311', shrink=0.05, width=1.4, headwidth=8),
             fontsize=15, fontweight='bold', color='#CC3311',
             bbox=dict(boxstyle='round,pad=0.3', fc='white', ec='#CC3311', alpha=0.85))

ax1.scatter([6], [0.6499], s=120, color='#CC3311', zorder=5)

ax1.set_xlim(-0.5, 6.5)
ax1.set_xticks(x)
ax1.set_xticklabels(stages, fontsize=13, rotation=30, ha='right')
ax1.set_ylim(0.2, 0.80)
ax1.tick_params(axis='y', labelsize=14)
ax1.set_ylabel("AP@0.5 Score", fontsize=17, fontweight='bold')
ax1.set_title("C. Evolutionary Discovery Path", fontsize=20, fontweight='bold')
ax1.grid(axis='y', ls='--', alpha=0.5)
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)

# --- Panel d: Parameter Sensitivity Heatmap (Top RHS) ---
ax2 = plt.subplot2grid((2, 7), (0, 3), colspan=4)

diameters = [13, 17, 22, 26, 30]
dilations = [4, 6, 8, 10, 12, 14, 16]
z = np.array([
    [0.32, 0.45, 0.48, 0.46, 0.42, 0.35, 0.30],
    [0.44, 0.648, 0.62, 0.58, 0.41, 0.33, 0.31],
    [0.40, 0.58, 0.60, 0.55, 0.45, 0.36, 0.32],
    [0.35, 0.52, 0.54, 0.50, 0.38, 0.32, 0.28],
    [0.25, 0.40, 0.42, 0.35, 0.22, 0.15, 0.12]
])

sns.heatmap(z, xticklabels=dilations, yticklabels=diameters,
            annot=True, fmt=".2f", cmap="RdYlGn",
            vmin=0.1, vmax=0.7, ax=ax2,
            annot_kws={"fontsize": 14},
            cbar_kws={'label': 'AP@0.5'})

cbar = ax2.collections[0].colorbar
cbar.ax.tick_params(labelsize=13)
cbar.ax.set_ylabel('AP@0.5', fontsize=15)

# 最優點標記（diameter=17 → row 1, expansion=6 → col 1）
ax2.add_patch(plt.Rectangle((1, 1), 1, 1, fill=False, edgecolor='#003366', lw=3, zorder=5))
ax2.plot(1.5, 0.75, marker='*', markersize=20, color='#003366', zorder=6, linestyle='none')

ax2.set_title("D. Parameter Sensitivity (Heatmap)", fontsize=20, fontweight='bold')
ax2.set_xlabel("Expansion Distance (px)", fontsize=16)
ax2.set_ylabel("Diameter (px)", fontsize=16)
ax2.tick_params(axis='both', labelsize=14)
plt.setp(ax2.get_yticklabels(), rotation=0)

# --- Panel e: Component Ablation Study (Bottom RHS) ---
ax3 = plt.subplot2grid((2, 7), (1, 3), colspan=4)

methods = ['Nuclei-Only', 'Proseg (Ref)', 'No Ensemble', 'No CLAHE', 'No Voronoi', 'MCseg v2 Optimum']
scores_e = [0.317, 0.383, 0.432, 0.525, 0.619, 0.650]
colors = ['#994422', '#CC4400', '#5599cc', '#5599cc', '#5599cc', '#003d80']

bars = ax3.barh(methods, scores_e, color=colors, height=0.6)

for bar, s in zip(bars, scores_e):
    ax3.text(s + 0.008, bar.get_y() + bar.get_height() / 2,
             f'{s:.3f}', va='center', ha='left', fontsize=14, fontweight='bold', color='#333333')

legend_handles = [
    mpatches.Patch(color='#003d80', label='MCseg v2 Full'),
    mpatches.Patch(color='#5599cc', label='Ablation (−1 component)'),
    mpatches.Patch(color='#CC4400', label='External Reference'),
]
ax3.legend(handles=legend_handles, fontsize=13, loc='lower right', framealpha=0.9)

ax3.set_title("E. Component Ablation Study", fontsize=20, fontweight='bold')
ax3.set_xlim(0, 0.80)
ax3.tick_params(axis='both', labelsize=14)
ax3.grid(axis='x', ls='--', alpha=0.5)
ax3.spines['top'].set_visible(False)
ax3.spines['right'].set_visible(False)

plt.tight_layout()

out_path = Path("/Volumes/SSD/plan_a/manuscript/figures/fig1/fig1cde.png")
plt.savefig(out_path, bbox_inches='tight', facecolor='white')
plt.close()

print(f"✅ Figure 1 c-e generated at: {out_path}")
