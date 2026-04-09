import matplotlib.pyplot as plt
import matplotlib.patches as patches

# 白色背景，解決黑屏問題
fig, ax = plt.subplots(figsize=(16, 9), facecolor='white')
ax.set_xlim(0, 24)
ax.set_ylim(0, 12)
ax.axis('off')

# 系統顏色設定
blue_colors = ['#E1F0FA', '#C1E0F5', '#A1D0F0', '#81C0EB']
orange_colors = ['#FFF3E0', '#FFE0B2', '#FFCC80']
gray_bg = '#F8F9FA'
stroke = '#2B2D42'
text_color = '#1A1A1A'

# 主標題
ax.text(12, 11.2, "MCseg v2: Spatial Transcriptomics Segmentation Pipeline", ha='center', va='center', fontsize=22, fontweight='black', color='#1D3557', family='sans-serif')

def draw_section(x, y, w, h, title, subtitle, color, ec):
    # 畫陰影
    shadow = patches.FancyBboxPatch((x+0.1, y-0.1), w, h, boxstyle="round,pad=0.3,rounding_size=0.3", facecolor='#d0d0d0', edgecolor='none', zorder=1)
    ax.add_patch(shadow)
    
    # 畫主方塊
    box = patches.FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.3,rounding_size=0.3", facecolor=color, edgecolor=ec, linewidth=2.5, zorder=2)
    ax.add_patch(box)
    
    # 段落標題
    ax.text(x + w/2, y + h - 0.7, title, ha='center', va='center', fontsize=16, fontweight='bold', color=text_color, zorder=3)
    if subtitle:
        ax.text(x + w/2, y + h - 1.4, subtitle, ha='center', va='center', fontsize=12, color='#555555', zorder=3)

# ================= 1. INPUT 區域 =================
draw_section(1, 2, 5, 8, "Input Data", "", '#FAFAFA', '#4A4E69')
y_inp = 7.6
for i, txt in enumerate(["High-Res H&E Image", "2µm Spatial Bins", "ROI Designation"]):
    bx = patches.FancyBboxPatch((1.6, y_inp), 3.8, 1.2, boxstyle="round,pad=0.1", facecolor='#E3E8ED', edgecolor=stroke, linewidth=1.5, zorder=3)
    ax.add_patch(bx)
    ax.text(3.5, y_inp + 0.6, txt, ha='center', va='center', fontsize=12, fontweight='bold', zorder=4)
    y_inp -= 1.8

# 左到中：大指示箭頭
ax.annotate('', xy=(7.3, 6), xytext=(6.2, 6), arrowprops=dict(arrowstyle="-|>,head_width=0.6,head_length=0.8", color=stroke, lw=4), zorder=0)

# ================= 2. 中間引擎區 (MCseg v2) =================
draw_section(7.5, 0.5, 9, 9.5, "MCseg v2 AI Ensemble", "", '#FFFDF7', '#D90429')

y_pos = 8.8
ax.text(12, y_pos, "Stage 1: Cyto3 Scaffold (Multi-Diameter)", ha='center', va='center', fontsize=12, fontweight='bold', color='#1D3557', zorder=4)
y_pos -= 0.6
for txt, c in zip([
    "Pass 1: Cyto3 Dia=17 | CLAHE-RGB",
    "Pass 2: Cyto3 Dia=17 | Hematoxylin",
    "Pass 3: Cyto3 Dia=22 | CLAHE-RGB",
    "Pass 4: Cyto3 Dia=13 | CLAHE-RGB"
], blue_colors):
    bx = patches.FancyBboxPatch((8.5, y_pos-0.35), 7, 0.7, boxstyle="round,pad=0.1", facecolor=c, edgecolor=stroke, linewidth=1.5, zorder=3)
    ax.add_patch(bx)
    ax.text(12, y_pos, txt, ha='center', va='center', fontsize=11, fontweight='bold', zorder=4)
    ax.annotate('', xy=(12, y_pos-0.45), xytext=(12, y_pos-0.35), arrowprops=dict(arrowstyle="-|>,head_width=0.2,head_length=0.2", color=stroke, lw=2), zorder=0)
    y_pos -= 0.9

y_pos -= 0.1
ax.text(12, y_pos, "Stage 2: CPSAM Attention Rescue", ha='center', va='center', fontsize=12, fontweight='bold', color='#D68C45', zorder=4)

y_pos -= 0.6
for txt, c in zip([
    "Pass 5: CPSAM Auto | CLAHE-RGB",
    "Pass 6: CPSAM Dia=16 | CLAHE-RGB",
    "Pass 7: CPSAM Auto | Hematoxylin"
], orange_colors):
    bx = patches.FancyBboxPatch((8.5, y_pos-0.35), 7, 0.7, boxstyle="round,pad=0.1", facecolor=c, edgecolor=stroke, linewidth=1.5, zorder=3)
    ax.add_patch(bx)
    ax.text(12, y_pos, txt, ha='center', va='center', fontsize=11, fontweight='bold', zorder=4)
    if txt != "Pass 7: CPSAM Auto | Hematoxylin":
        ax.annotate('', xy=(12, y_pos-0.45), xytext=(12, y_pos-0.35), arrowprops=dict(arrowstyle="-|>,head_width=0.2,head_length=0.2", color=stroke, lw=2), zorder=0)
    y_pos -= 0.9

# 針對 Fusion
ax.annotate('', xy=(12, y_pos+0.3), xytext=(12, y_pos+0.5), arrowprops=dict(arrowstyle="-|>,head_width=0.3,head_length=0.3", color='#D90429', lw=3), zorder=0)

y_pos -= 0.2
bx = patches.FancyBboxPatch((8.5, y_pos-0.4), 7, 1.0, boxstyle="round,pad=0.1", facecolor='#2A9D8F', edgecolor=stroke, linewidth=1.5, zorder=3)
ax.add_patch(bx)
ax.text(12, y_pos+0.2, "Spatial Non-Overlapping Merge", ha='center', va='center', fontsize=11, fontweight='bold', color='white', zorder=4)
ax.text(12, y_pos-0.2, "& Adaptive Voronoi Expansion", ha='center', va='center', fontsize=11, fontweight='bold', color='white', zorder=4)

# 中到右：大指示箭頭
ax.annotate('', xy=(17.8, 6), xytext=(16.7, 6), arrowprops=dict(arrowstyle="-|>,head_width=0.6,head_length=0.8", color=stroke, lw=4), zorder=0)

# ================= 3. OUTPUT 區域 =================
draw_section(18, 2, 5, 8, "Output & Validation", "", '#FAFAFA', '#4A4E69')

bx1 = patches.FancyBboxPatch((18.5, 7.8), 4, 1.2, boxstyle="round,pad=0.1", facecolor='#D8E2DC', edgecolor=stroke, linewidth=1.5, zorder=3)
ax.add_patch(bx1)
ax.text(20.5, 8.6, "High-Fidelity Masks", ha='center', va='center', fontsize=12, fontweight='bold', zorder=4)
ax.text(20.5, 8.1, "& Single-Cell RNA", ha='center', va='center', fontsize=12, fontweight='bold', zorder=4)

ax.annotate('', xy=(20.5, 7.3), xytext=(20.5, 7.8), arrowprops=dict(arrowstyle="-|>,head_width=0.4,head_length=0.4", color=stroke, lw=3), zorder=0)

ax.text(20.5, 6.8, "Evaluation Metrics", ha='center', va='center', fontsize=13, fontweight='bold', color='#E63946', zorder=4)
for i, txt in enumerate(["1. FTC (Tissue Capture)", "2. UMI Density", "3. NED (Divergence)"]):
    bx = patches.FancyBboxPatch((18.5, 5.2 - i*1.3), 4, 0.9, boxstyle="round,pad=0.1", facecolor='#FFF3E0', edgecolor=stroke, linewidth=1.5, zorder=3)
    ax.add_patch(bx)
    ax.text(20.5, 5.65 - i*1.3, txt, ha='center', va='center', fontsize=10, fontweight='bold', zorder=4)

# 儲存圖片，強制寫入白色背景以防止在 Markdown Dark 模式下文字消失！
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
plt.savefig('/Volumes/SSD/plan_a/manuscript/figures/fig1/fig1a_full_redesign_v2.svg', format='svg', facecolor=fig.get_facecolor(), transparent=False)
plt.savefig('/Volumes/SSD/plan_a/manuscript/figures/fig1/fig1a_full_redesign_v2.png', dpi=300, facecolor=fig.get_facecolor(), transparent=False)
