# AutoResearch Agent Task Specification: H&E Cell Segmentation

This file is the task description provided to the AI agent at the start of each AutoResearch session. Adapt it to your own tissue type, image specifications, and available tools.

---

## Objective

Your sole objective is to maximise the **AP@0.5** (average precision at IoU ≥ 0.5) score between your predicted cell instance masks and a reference ground truth, by freely modifying `segment.py`.

## Core Principle: No Boundaries

You are **not** restricted to tuning hyperparameters. You may:

1. **Invent new algorithms**: custom watershed, polygon stripping, K-means clustering.
2. **Build multi-model ensembles**: run Cellpose `cyto2`, `cyto3`, `nuclei`, `cpsam` simultaneously and merge predictions by union, intersection, or non-maximum suppression.
3. **Apply image science**: HED stain deconvolution, Sobel edge enhancement, morphological dilation/erosion, adaptive thresholding.
4. **Abandon Cellpose entirely**: if classical OpenCV + transcript density thresholding can achieve a higher score, implement it.

## Rules

1. **Sandbox scope**: you may **only modify** `segment.py` in the same directory. The function signature must be preserved:
   ```python
   def build_and_predict(img, vhd_csv, gt_mask=None) -> np.ndarray:
       """Returns 2D int32 array: 0=background, >0=cell instance ID"""
   ```
2. **Time budget**: each run must complete within **300 seconds**. Models that time out count as failures.
3. **Available packages**: `cellpose`, `scikit-image`, `opencv-python`, `scipy`, `pandas`, `numpy` (and optionally `proseg` binary via `subprocess`).

## Evaluation

Run scoring with:
```bash
python segment.py
```

Observe the `Validation Score` output (AP@0.5 against ground truth).

**Target: AP@0.5 > 0.70**

## Scientific Context (adapt to your tissue)

- **Image**: H&E RGB numpy array, ~1460×1460 px, 0.2737 µm/px
- **Tissue**: LUAD lung adenocarcinoma (tumour boundary region)
- **Ground truth**: Xenium single-molecule cell boundaries (cytoplasm-inclusive, ~5 µm radius beyond nucleus)
- **Key challenge**: cells range from small lymphocytes (~8 px diameter) to large pleomorphic tumour cells (~30 px); any single diameter will miss one extreme

## Known Insights (update as experiments accumulate)

- CLAHE contrast enhancement (clip=3.0, tile=8) consistently improves detection
- Voronoi-constrained expansion outperforms `expand_labels` for inter-cell space allocation
- `cyto3` at multiple diameters (13, 17, 22 px) captures the size spectrum better than a single pass
- `cpsam` (Cellpose-SAM) adds complementary detections when used alongside `cyto3`
- Per-cell morphological loops cause timeouts — use vectorised operations only
- Watershed (`cv2.watershed`, `skimage.watershed`) consistently exceeds the time budget

## Response Format (strict)

```
APPROACH: [one-line description of your strategy]
RATIONALE: [2-3 sentences of scientific reasoning]
[complete Python code starting with import statements]
```
