# Supplementary Note 1: MCseg v2 Algorithm Specification and AutoResearch Sandbox

## A. AutoResearch Sandbox Configuration

MCseg v2 was developed through automated agent-driven optimisation using the AutoResearch framework. The agent operated within a strictly sandboxed environment:

- **Evaluation metric**: AP@0.5 (average precision at IoU ≥ 0.5) against Xenium-derived ground-truth masks
- **Cycle budget**: ~80 overnight cycles, each completing within <5 minutes on a CPU-only macOS workstation (no GPU)
- **Modifiable scope**: a single `segment.py` script, preventing unintended side effects on other pipeline components
- **Available primitives**: Cellpose models (cyto2, cyto3, nuclei, cpsam), image preprocessing (CLAHE, HED stain deconvolution, Gaussian smoothing), post-processing (watershed, Voronoi tessellation, morphological operations)
- **Search space**: model selection and combination, diameter values, CLAHE parameters (clip limit, tile grid), ensemble overlap thresholds, Voronoi seed distance, area filter bounds

The agent autonomously discovered the 7-pass ensemble architecture and adaptive Voronoi expansion strategy. AP@0.5 improved from 0.32 (single-model baseline) to 0.65 (final MCseg v2), a 2× relative gain.

## B. MCseg v2 Pipeline Specification

The pipeline processes a single H&E-stained ROI image and produces a labelled cell mask.

### Step 1 — Image preprocessing

- CLAHE contrast enhancement: clip limit = 3.0, tile grid = 8 × 8
- Hematoxylin channel extraction via HED colour deconvolution (Macenko method)
- Outputs: CLAHE-RGB composite and hematoxylin single-channel image

### Step 2 — 7-pass multi-model cell detection

| Pass | Model | Input | Diameter (px) |
|------|-------|-------|:-------------:|
| 1 | cyto3 | CLAHE-RGB | 13 |
| 2 | cyto3 | CLAHE-RGB | 17 |
| 3 | cyto3 | CLAHE-RGB | 22 |
| 4 | cyto3 | Hematoxylin | 13 |
| 5 | cpsam | CLAHE-RGB | auto |
| 6 | cpsam | CLAHE-RGB | 16 |
| 7 | cpsam | Hematoxylin | auto |

All passes use `flow_threshold = 0.4`, `cellprob_threshold = 0.0`.

### Step 3 — Ensemble merging

Masks from all 7 passes are merged using greedy non-maximum suppression. Any two candidate masks with pixel overlap > 15% are resolved by retaining the higher-confidence mask. This preserves complementary detections from different model configurations while eliminating redundant proposals.

### Step 4 — Voronoi boundary expansion

Cell centroids from the merged ensemble serve as Voronoi seeds. Unclaimed pixels (not assigned to any detected cell) are partitioned by nearest-centroid assignment using a fixed seed distance *d* = 8 px in deployment mode. This step recovers cytoplasmic RNA signal that would otherwise fall within inter-mask gaps.

### Step 5 — Quality filtering

Final masks are filtered by area: cells with area < 20 px² (sub-nuclear fragments) or > 6000 px² (stitching artefacts or tissue debris) are discarded.
