"""
Assemble Supplementary Figure S11.

Layout (4 panels):
  Row 0: [  a: Treg NN violin (3 pairs)  ] [ b: Contact bar ]
  Row 1: [       c: Microenv boxplot              ]
  Row 2: [       d: Immunosuppressive triad       ]
"""
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from pathlib import Path

FIG_DIR = Path('/Volumes/SSD/plan_a/submission_bioinformatics/analysis/treg_crc/figures')

img_a = mpimg.imread(FIG_DIR / 'panelB1_nn_violin.png')
img_b = mpimg.imread(FIG_DIR / 'panelB2_contact_bar.png')
img_c = mpimg.imread(FIG_DIR / 'panelC2_microenv_boxplot.png')
img_d = mpimg.imread(FIG_DIR / 'panelD_triad_triangle.png')

fig = plt.figure(figsize=(13, 14), constrained_layout=True)
gs = fig.add_gridspec(3, 3, height_ratios=[1, 1, 1])

ax_a = fig.add_subplot(gs[0, :2])
ax_a.imshow(img_a)
ax_a.axis('off')
ax_a.text(-0.04, 1.06, 'a', transform=ax_a.transAxes, fontsize=18, fontweight='bold', va='top')

ax_b = fig.add_subplot(gs[0, 2])
ax_b.imshow(img_b)
ax_b.axis('off')
ax_b.text(-0.08, 1.06, 'b', transform=ax_b.transAxes, fontsize=18, fontweight='bold', va='top')

ax_c = fig.add_subplot(gs[1, :])
ax_c.imshow(img_c)
ax_c.axis('off')
ax_c.text(-0.02, 1.06, 'c', transform=ax_c.transAxes, fontsize=18, fontweight='bold', va='top')

ax_d = fig.add_subplot(gs[2, :])
ax_d.imshow(img_d)
ax_d.axis('off')
ax_d.text(-0.02, 1.06, 'd', transform=ax_d.transAxes, fontsize=18, fontweight='bold', va='top')

plt.savefig(FIG_DIR / 'SuppFig_S11_treg_spatial.png', dpi=300, bbox_inches='tight')
print(f"Saved -> {FIG_DIR / 'SuppFig_S11_treg_spatial.png'}")
