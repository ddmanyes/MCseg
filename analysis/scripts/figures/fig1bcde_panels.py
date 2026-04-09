import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from pathlib import Path

# 設置科學論文風格
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica']

# 定義圖片路徑
base_dir = Path("/Volumes/SSD/plan_a/manuscript/figures")
path_a = base_dir / "01_intro" / "fig1_pipeline_schematic.png"
path_b = base_dir / "02_methods" / "fig_v12_comparison_3panel.png"
path_c = base_dir / "02_methods" / "fig_v12_sanity_strip.png"

# 檢查檔案是否存在
for p in [path_a, path_b, path_c]:
    if not p.exists():
        print(f"❌ Error: File not found {p}")

# 繪製總合圖 (1x2 佈局，頂部為 A，下方為 B, C)
fig = plt.figure(figsize=(16, 12), dpi=300)
fig.set_facecolor('white')

# --- Panel A: Pipeline (Top) ---
ax1 = plt.subplot2grid((3, 2), (0, 0), colspan=2)
img_a = mpimg.imread(path_a)
ax1.imshow(img_a)
ax1.axis('off')
ax1.text(0, 0.95, "A", transform=ax1.transAxes, fontsize=24, fontweight='bold', va='top')

# --- Panel B: LUAD Comparison (Bottom Left) ---
ax2 = plt.subplot2grid((3, 2), (1, 0), rowspan=2)
img_b = mpimg.imread(path_b)
ax2.imshow(img_b)
ax2.axis('off')
ax2.text(0, 0.95, "B", transform=ax2.transAxes, fontsize=24, fontweight='bold', va='top')

# --- Panel C: CRC Sanity Strip (Bottom Right) ---
ax3 = plt.subplot2grid((3, 2), (1, 1), rowspan=2)
img_c = mpimg.imread(path_c)
ax3.imshow(img_c)
ax3.axis('off')
ax3.text(0, 0.95, "C", transform=ax3.transAxes, fontsize=24, fontweight='bold', va='top')

plt.tight_layout()

# 保存
out_dir = Path("/Volumes/SSD/plan_a/manuscript/figures/01_intro")
out_path = out_dir / "fig1_composite_flagship.png"
plt.savefig(out_path, bbox_inches='tight', facecolor='white')
plt.close()

print(f"✅ Flagship Composite Figure 1 generated at: {out_path}")
