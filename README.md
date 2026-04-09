# MCseg

**MCseg** is a no-code, browser-based platform for high-fidelity cell segmentation
of [10x Genomics Visium HD](https://www.10xgenomics.com/products/visium-hd-spatial-gene-expression)
spatial transcriptomics data.  
Its segmentation engine was discovered through AI-autonomous experimentation and runs
entirely in the browser — no installation, no server, no code.

---

## Features

| Feature | Description |
|---------|-------------|
| 🔬 **H&E Stain Deconvolution** | Isolates the hematoxylin channel via Ruifrok & Johnston (2001) stain vectors |
| ⚙️ **Marker-Controlled Expansion** | Nucleus seeds grown to cell bodies using BFS nearest-nucleus expansion |
| 📍 **Spatial Barcode Assignment** | Maps Visium HD barcodes to the nearest cell centroid |
| 🎨 **Interactive Visualization** | Zoom, pan, and four overlay modes (Original / Overlay / Borders / Mask) |
| ⬇️ **Export** | Download segmentation masks (PNG) and cell tables (CSV) |
| 🌙 **Dark / Light Theme** | Respects `prefers-color-scheme`; toggle in the header |

---

## Quick Start

1. **Open `index.html`** in a modern web browser (Chrome ≥ 90, Firefox ≥ 88, Edge ≥ 90, Safari ≥ 15).
2. **Upload an H&E image** — drag-and-drop or click the upload zone.  
   Accepted formats: PNG, JPEG, TIFF.
3. **(Optional) Upload Spatial Barcodes** — Visium HD `tissue_positions.csv` or any CSV with `x`/`y` columns.
4. **Adjust Parameters** in the sidebar:
   - *Min / Max Nucleus Area* — filter nuclei by area (px²)
   - *Cell Expansion* — how many pixels to grow cell bodies beyond nuclei
   - *Nucleus Sensitivity* — 0 = conservative (fewer detections), 100 = aggressive
5. **Click ▶ Run MCseg** and watch the progress bar.
6. **Explore results** with the toolbar view buttons.
7. **Export** the mask PNG or cell CSV with the export buttons.

---

## Segmentation Pipeline

```
H&E Image (RGB)
      │
      ▼
1. H&E Optical-Density Decomposition
   Project RGB → OD space; project onto hematoxylin stain vector
   (Ruifrok & Johnston, 2001) → normalised hematoxylin channel
      │
      ▼
2. Separable Gaussian Smoothing  (σ = 1.5)
      │
      ▼
3. Otsu Global Threshold  +  Sensitivity Adjustment
      │
      ▼
4. Morphological Opening  (erosion r=2 → dilation r=2)
   Removes small noise / debris
      │
      ▼
5. Union-Find Connected-Component Labelling  (4-connectivity)
      │
      ▼
6. Area Filtering  (min / max nucleus area)
      │
      ▼
7. BFS Nearest-Nucleus Expansion  (cell body approximation)
      │
      ▼
8. Per-Cell Statistics  (area, centroid, mean RGB)
      │
      ▼
9. Optional Spatial-Barcode Assignment
   Each barcode → nearest cell centroid within 60 px
      │
      ▼
Segmentation Result  (label map, cell table, PNG mask, CSV export)
```

---

## Supported Input Formats

### Tissue Images
| Format | Extension |
|--------|-----------|
| PNG | `.png` |
| JPEG | `.jpg`, `.jpeg` |
| TIFF | `.tif`, `.tiff` *(requires browser TIFF support)* |

Images larger than 3 000 × 3 000 px are automatically downsampled before processing.

### Spatial Barcodes (Visium HD)
CSV or TSV with the following columns (column names are case-insensitive):

| Column | Accepted names |
|--------|---------------|
| Barcode ID | `barcode`, `barcodes`, `barcode_id` |
| X coordinate | `x`, `col`, `pxl_col_in_fullres` |
| Y coordinate | `y`, `row`, `pxl_row_in_fullres` |

---

## Export Formats

| Export | Description |
|--------|-------------|
| **Mask PNG** | RGBA image; each cell drawn in a unique colour, background transparent |
| **Cells CSV** | `cell_id, cx, cy, area, mean_r, mean_g, mean_b, barcodes` |

---

## Project Structure

```
MCseg/
├── index.html            # Single-page application
├── css/
│   └── style.css         # Responsive dark/light-theme styles
├── js/
│   ├── utils.js          # Shared helpers (CSV parser, download, colour)
│   ├── segmentation.js   # MCsegEngine — core segmentation pipeline
│   ├── visualization.js  # Visualizer — canvas renderer, zoom/pan
│   └── app.js            # Application controller
├── LICENSE
└── README.md
```

---

## Browser Requirements

- Chrome ≥ 90 · Firefox ≥ 88 · Edge ≥ 90 · Safari ≥ 15
- JavaScript and the Canvas API enabled
- No installation or server required — runs entirely client-side

---

## License

[MIT](LICENSE) © 2026 ddmanyes
