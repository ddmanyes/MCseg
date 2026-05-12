import numpy as np
import pandas as pd
from scipy.ndimage import find_objects
from pathlib import Path

ROOT = Path('/Volumes/SSD/plan_a/submission_bioinformatics/analysis/treg_crc/results/highres_seg')
OUT  = Path('/Volumes/SSD/plan_a/submission_bioinformatics/analysis/treg_crc/results/highres_seg')

mask = np.load(ROOT / 'mcseg_mask.npy')
ct   = pd.read_csv(ROOT / 'celltypist_labels.csv')

print(f'Mask shape: {mask.shape}  unique cells: {len(np.unique(mask[mask>0]))}')

slices = find_objects(mask)
rows, cols, ids = [], [], []
for i, sl in enumerate(slices):
    if sl is None:
        continue
    cid = i + 1
    sub = mask[sl]
    r, c = np.where(sub == cid)
    if len(r) == 0:
        continue
    rows.append(r.mean() + sl[0].start)
    cols.append(c.mean() + sl[1].start)
    ids.append(cid)

centroids = pd.DataFrame({'cell_id': ids, 'row_px': rows, 'col_px': cols})
# Assuming the full resolution is 0.2737 um/px (Visium HD standard fullres pixel size)
centroids['row_um'] = centroids['row_px'] * 0.2737
centroids['col_um'] = centroids['col_px'] * 0.2737

df = centroids.merge(ct, on='cell_id', how='left')

SUBTYPE_MAP = {
    'Regulatory T cells':        'Treg',
    'CD8+ T cells':              'CD8',
    'CD4+ T cells':              'CD4',
    'T helper 17 cells':         'Th17',
    'T follicular helper cells': 'Tfh',
    'NK cells':                  'NK',
    'Pro-inflammatory':          'ProInf_macro',
    'SPP1+':                     'SPP1_macro',
    'IgG+ Plasma':               'Plasma',
    'CMS1':                      'Tumor',
    'CMS2':                      'Tumor',
    'CMS3':                      'Tumor',
    'Myofibroblasts':            'Stromal',
    'Stromal 1':                 'Stromal',
    'Stromal 3':                 'Stromal',
}
df['subtype'] = df['celltypist_label'].map(SUBTYPE_MAP).fillna(df['broad_label'])

out_path = OUT / 'cell_centroids.csv'
df.to_csv(out_path, index=False)
print(f'Saved {len(df)} rows -> {out_path}')

print('\nSubtype counts:')
print(df['subtype'].value_counts().head(15).to_string())
