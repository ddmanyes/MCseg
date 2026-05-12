import warnings
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.ndimage import minimum_filter
from pathlib import Path
from PIL import Image

Image.MAX_IMAGE_PIXELS = None

ROOT_DIR = Path('/Volumes/SSD/plan_a/submission_bioinformatics/analysis/treg_crc/results/highres_seg')
FIG_DIR  = Path('/Volumes/SSD/plan_a/submission_bioinformatics/analysis/treg_crc/figures')

HE_PATH  = ROOT_DIR / 'he_crop.tif'
MSK_PATH = ROOT_DIR / 'mcseg_mask.npy'
CELL_PATH = ROOT_DIR / 'cell_centroids.csv'

ALPHA_MASK = 0.72   # fig2h style: H&E shows through as texture

SUBTYPE_COLORS = {
    'Tumor':            '#2ecc71',
    'epithelial cells': '#27ae60',
    'CD8':              '#2980b9',
    'CD4':              '#85c1e9',
    'Treg':             '#e74c3c',
    'Th17':             '#f39c12',
    'Tfh':              '#e67e22',
    'NK':               '#1abc9c',
    'Plasma':           '#9b59b6',
    'ProInf_macro':     '#f1c40f',
    'SPP1_macro':       '#e67e22',
    'Stromal':          '#95a5a6',
    'stromal cells':    '#aab7b8',
    'immune cells':     '#76d7c4',
    'other':            '#bdc3c7',
}
LEGEND_RENAME = {
    'SPP1_macro':       'SPP1⁺ Macrophage',
    'ProInf_macro':     'Pro-inflam Macro',
    'stromal cells':    'Stromal',
    'immune cells':     'Immune (misc)',
    'epithelial cells': 'Epithelial',
}
LEGEND_ORDER = ['Tumor', 'CD8', 'Treg', 'SPP1_macro', 'CD4', 'Th17',
                'NK', 'Plasma', 'ProInf_macro', 'Stromal', 'other']

ROIS = [
    {'name': 'ROI 1', 'tag': 'roi1_highres',
     'r0': 4991, 'r1': 7182, 'c0': 7118, 'c1': 9308},
    {'name': 'ROI 2', 'tag': 'roi2_highres',
     'r0': 7456, 'r1': 9646, 'c0': 7855, 'c1': 10045},
]

# ── Load shared data ───────────────────────────────────────────────────────
print("Loading H&E image...")
he_arr = np.array(Image.open(HE_PATH))
print(f"  shape: {he_arr.shape}")

print("Loading MCseg mask...")
mask_full = np.load(MSK_PATH)
print(f"  shape: {mask_full.shape}")

print("Loading cell centroids...")
cells       = pd.read_csv(CELL_PATH)
cell_lookup = cells.set_index('cell_id')['subtype'].to_dict()

def hex_to_rgb(h):
    return (int(h[1:3], 16)/255, int(h[3:5], 16)/255, int(h[5:7], 16)/255)

# ── Render each ROI ────────────────────────────────────────────────────────
for roi in ROIS:
    r0, c0 = roi['r0'], roi['c0']
    r1 = min(roi['r1'], mask_full.shape[0] - 1)
    c1 = min(roi['c1'], mask_full.shape[1] - 1)
    print(f"\n{roi['name']}  [{r0}:{r1}, {c0}:{c1}]")

    he_roi   = he_arr[r0:r1, c0:c1]
    mask_roi = mask_full[r0:r1, c0:c1]
    H, W     = he_roi.shape[:2]
    # Assuming 0.2737 um/px instead of 0.5 um/px
    pixel_size_um = 0.2737
    print(f"  {H}×{W} px  (~{H*pixel_size_um/1000:.2f}×{W*pixel_size_um/1000:.2f} mm)")

    # Build RGBA cell mask
    cell_ids = np.unique(mask_roi[mask_roi > 0])
    print(f"  Cells: {len(cell_ids)}")

    rgba = np.zeros((H, W, 4), dtype=np.float32)
    for cid in cell_ids:
        sub = cell_lookup.get(cid, 'other')
        rv, gv, bv = hex_to_rgb(SUBTYPE_COLORS.get(sub, '#bdc3c7'))
        px = mask_roi == cid
        rgba[px, 0] = rv
        rgba[px, 1] = gv
        rgba[px, 2] = bv
        rgba[px, 3] = ALPHA_MASK

    # 1-px dark border on cell boundaries
    eroded   = minimum_filter(mask_roi, size=2)
    boundary = (mask_roi != eroded) & (mask_roi > 0)
    rgba[boundary] = [0.05, 0.05, 0.05, 0.90]

    # ── Figure: single panel — H&E + cell mask ────────────────────────────
    fig, ax = plt.subplots(figsize=(7, 7))
    # 200 um scale bar
    scale_bar_px = int(200 / pixel_size_um)

    ax.imshow(he_roi, origin='upper')
    ax.imshow(rgba, origin='upper', interpolation='nearest')
    ax.set_title(f'{roi["name"]} — MCseg cell types (High Res)', fontsize=12, fontweight='bold')
    ax.set_xlabel('µm', fontsize=9)
    ax.set_ylabel('µm', fontsize=9)

    ticks_x = np.linspace(0, W, 5)
    ticks_y = np.linspace(0, H, 5)
    ax.set_xticks(ticks_x); ax.set_xticklabels([f'{int(t*pixel_size_um)}' for t in ticks_x], fontsize=7)
    ax.set_yticks(ticks_y); ax.set_yticklabels([f'{int(t*pixel_size_um)}' for t in ticks_y], fontsize=7)

    x0 = W * 0.05; y = H * 0.93
    ax.plot([x0, x0 + scale_bar_px], [y, y], color='white', linewidth=2.5, solid_capstyle='butt')
    ax.text(x0 + scale_bar_px/2, H*0.965, '200 µm',
            ha='center', va='top', fontsize=9, color='white', fontweight='bold')
    ax.spines[['top', 'right']].set_visible(False)

    present = {cell_lookup.get(c, 'other') for c in cell_ids}
    patches = [
        mpatches.Patch(facecolor=SUBTYPE_COLORS[s], edgecolor='#333',
                       linewidth=0.4, label=LEGEND_RENAME.get(s, s))
        for s in LEGEND_ORDER if s in present and s in SUBTYPE_COLORS
    ]
    ax.legend(handles=patches, fontsize=7.5, loc='lower right',
              framealpha=0.90, handlelength=1.2, ncol=2, columnspacing=0.5)

    plt.tight_layout()
    out_path = FIG_DIR / f'panelDE_{roi["tag"]}_hemask.png'
    plt.savefig(out_path, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"  Saved -> {out_path.name}")

print("\nDone.")
