import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from pathlib import Path

# 設置科學論文繪圖風格
sns.set_theme(style="white", context="paper")
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']
plt.rcParams['axes.labelweight'] = 'bold'
plt.rcParams['axes.titleweight'] = 'bold'

# V12 演化數據 (來自 optimization_history.md)
stages = [
    "Phase 1\nInitial Search",
    "Agent Best\n(Proseg)",
    "Stage 2\nCLAHE+Dilation\n(BREAKTHROUGH)",
    "Stage 3\nEnsemble",
    "Stage 4\nVoronoi",
    "Stage 6\ncpsam",
    "MCseg v2\nFinal"
]
scores = [0.3176, 0.3830, 0.6483, 0.6333, 0.6194, 0.6443, 0.6499]

# 繪圖
fig, ax = plt.subplots(figsize=(10, 5), dpi=300)

# 背景分區陰影
ax.axvspan(-0.5, 1.5, color='#f0f0f0', alpha=0.5, label='Phase 1: Discovery')
ax.axvspan(1.5, 6.5, color='#e6f3ff', alpha=0.5, label='Phase 2: Evolution')

# 繪製曲線與數據點
line, = ax.plot(stages, scores, marker='o', markersize=8, color='#004488', linewidth=2.5, linestyle='-', zorder=3)
ax.fill_between(range(len(stages)), scores, alpha=0.1, color='#004488')

# 重點標註：Stage 2 的躍遷
ax.annotate('Major Leap: +0.33 AP\n(CLAHE + Optimal Dilation)', 
            xy=(2, 0.6483), xytext=(2.2, 0.52),
            arrowprops=dict(facecolor='black', shrink=0.05, width=1.5, headwidth=8),
            fontsize=10, fontweight='bold', color='#CC3311',
            bbox=dict(boxstyle='round,pad=0.3', fc='white', ec='#CC3311', alpha=0.8))

# 最終 V12 標註
ax.annotate(f'MCseg v2 Optimum: {scores[-1]:.4f}', 
            xy=(6, scores[-1]), xytext=(4.5, 0.68),
            arrowprops=dict(facecolor='black', shrink=0.05, width=1.5, headwidth=8),
            fontsize=10, fontweight='bold', color='#004488')

# 優化軸標籤
ax.set_ylabel('Segmentation Performance (AP@0.5)', fontsize=12)
ax.set_ylim(0.2, 0.75)
ax.grid(axis='y', linestyle='--', alpha=0.7)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# 加入標題
plt.title('Autonomous Evolutionary Discovery Path of MCseg v2 Pipeline', fontsize=14, pad=20)

# 保存圖片
out_dir = Path("/Volumes/SSD/plan_a/manuscript/figures/02_methods")
out_dir.mkdir(parents=True, exist_ok=True)
out_path = out_dir / "fig_v12_evolution_path.png"
plt.tight_layout()
plt.savefig(out_path)
plt.close()

print(f"✅ MCseg v2 Evolution Plot generated at: {out_path}")
