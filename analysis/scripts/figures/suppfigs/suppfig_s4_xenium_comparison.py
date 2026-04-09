import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from pathlib import Path

# 定義路徑
INPUT_DIR = Path("/Volumes/SSD/plan_a/crc_transcript_attribution/results/tls_exploration")
OUTPUT_PATH = Path("/Volumes/SSD/plan_a/manuscript/figures/supplementary/supp_fig_s4_tls_discovery_workflow.png")
OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)

# 定義圖片列表
IMAGES = [
    ("A: TLS Composite Score (Full Slide)", INPUT_DIR / "01_tls_composite_score_fullslice.png"),
    ("B: Marker Gene Panel (Full Slide)", INPUT_DIR / "02_tls_marker_panel_fullslice.png"),
    ("C: Local Moran's I Hotspots (p<0.01)", INPUT_DIR / "03_tls_local_morans_i.png"),
    ("D: Candidate 1 (Selected ROI 15)", INPUT_DIR / "tls_candidate_1_zoom.png"),
    ("E: Candidate 2 Zoom-in", INPUT_DIR / "tls_candidate_2_zoom.png"),
    ("F: Candidate 3 Zoom-in", INPUT_DIR / "tls_candidate_3_zoom.png"),
]

# 繪圖
fig, axes = plt.subplots(2, 3, figsize=(18, 11))
axes = axes.flatten()

for i, (title, img_path) in enumerate(IMAGES):
    if img_path.exists():
        img = mpimg.imread(img_path)
        axes[i].imshow(img)
        axes[i].set_title(title, fontsize=14, fontweight="bold")
    else:
        axes[i].text(0.5, 0.5, f"Missing:\n{img_path.name}", ha='center', va='center')
        print(f"⚠️ Warning: {img_path} not found")
    
    axes[i].axis("off")

plt.suptitle("Supplementary Figure S4: Objective ROI Discovery and Selection Workflow", 
             fontsize=20, fontweight="bold", y=0.98)
plt.tight_layout(rect=[0, 0, 1, 0.96])

plt.savefig(OUTPUT_PATH, dpi=200, bbox_inches="tight")
plt.close()

print(f"✅ Supplementary Figure S4 saved to {OUTPUT_PATH}")
