# MCseg — High-Fidelity Visium HD Cell Segmentation

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/)

**MCseg** is a no-code, end-to-end analysis platform for 10x Genomics **Visium HD** (2 µm resolution) spatial transcriptomics data. Its core segmentation engine, **MCseg v2**, uses a multi-pass ensemble (cyto3 at three diameters + optional hematoxylin and cpsam passes, up to 7 passes) with adaptive Voronoi boundary expansion to achieve high-fidelity cell segmentation. GPU is optional; CPU fallback is supported.

> MCseg v2 achieves **PQ = 0.554 ± 0.064** on LUAD tissue (vs 0.432 ± 0.037 for single-model Cellpose baseline, **+28% relative improvement**), validated against Xenium Prime ground-truth masks.

<p align="center">
  <img src="docs/fig1a_pipeline.png" width="820" alt="MCseg v2 pipeline overview">
</p>

---

## Contents

[Quick Start](#quick-start) · [Pipeline Overview](#pipeline-overview) · [Interface Tour](#interface-tour) · [Example Results](#example-results) · [Usage Guide](#usage-guide) · [Algorithm](#mcseg-v2-algorithm) · [Configuration](#configuration) · [Troubleshooting](#troubleshooting) · [Citation](#citation) · [License](#license)

---

## Quick Start

### System Requirements

| Component | Minimum | Recommended | Notes |
|-----------|---------|-------------|-------|
| **OS** | macOS 12 (Monterey) | macOS 13+ | Linux (Ubuntu 20.04+) also supported |
| **CPU** | 4-core, any modern x86-64 or ARM | Apple Silicon (M1/M2/M3) | Apple Silicon provides MPS GPU acceleration |
| **RAM** | 16 GB | 32 GB | Cellpose loads full ROI crops into memory; large BTFs need more |
| **Storage** | 15 GB free | 30 GB+ free | ~8 GB for Python env (torch, cellpose); remainder for data & results |
| **Python** | 3.10 | 3.11 | Managed by `uv`; do not use system Python |
| **Node.js** | v18 | v20 LTS | For frontend (Vite + React) |
| **GPU** | — (CPU fallback) | Apple MPS or NVIDIA (CUDA 12.4) | GPU cuts per-ROI segmentation from ~10 min to ~2 min |

> **No GPU?** CPU mode works but is slow for large ROIs. A single 1500 × 1200 px ROI takes ~8–12 min on an Apple M2 CPU vs ~2 min with MPS.

### Prerequisites

**macOS (Homebrew recommended):**
```bash
# Install Homebrew if not present
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

# Install Node.js (latest LTS)
brew install node
```

**Linux (Ubuntu/Debian):**
```bash
curl -fsSL https://deb.nodesource.com/setup_20.x | sudo -E bash -
sudo apt-get install -y nodejs
```

### Installation

```bash
# 1. Install uv (Python package manager)
curl -LsSf https://astral.sh/uv/install.sh | sh
source ~/.zshrc   # or restart your terminal — required for uv to be in PATH

# 2. Clone and install
git clone https://github.com/ddmanyes/MCseg.git
cd MCseg
uv sync           # skip this step if your drive is ExFAT — see note below

# 3. Install frontend dependencies
cd frontend && npm install && cd ..

# 4. Launch (also handles Python env setup)
bash start.sh
```

Open **<http://localhost:3000>** in your browser.

> **Note for bash users:** replace `source ~/.zshrc` with `source ~/.bashrc`.  
> **ExFAT / external drive users:** skip `uv sync` in step 2 and run `bash start.sh` directly — it creates `.venv` as a symlink to `~/.venvs/msseg` (APFS) before installing dependencies, avoiding resource-fork corruption.

---

## Pipeline Overview

| Stage | Function | Key Output |
|-------|----------|------------|
| Data Setup | Auto-scan and validate raw data | `state.json` |
| Stage 0: ROI Extract | Crop ROI from Gigapixel BTF | `he_crop.tif`, `adata_002um.h5ad` |
| Stage 1: MCseg v2 | Multi-pass ensemble segmentation (4–7 passes) + Voronoi expansion | `segmentation_masks.npy` |
| Stage 2: RNA Count | Assign Visium HD bins to cells | `cellpose_cells.h5ad` |
| Stage 3: Analysis | QC → Normalise → PCA → UMAP → Leiden | `umap_computed.h5ad` |
| Stage 3.5: Explorer | Interactive spatial gene expression viewer | PNG export |
| Stage 4: Export | Xenium Explorer / Loupe Browser format | GeoJSON, CSV |

---

## Interface Tour

<table>
  <tr>
    <td align="center" width="50%">
      <img src="docs/sample/Operation%20interface/stage1_steup.png" width="400" alt="Data Setup"><br>
      <sub><b>① Data Setup</b> — scan BTF + binned matrices, set output dir</sub>
    </td>
    <td align="center" width="50%">
      <img src="docs/sample/Operation%20interface/stage2_ROI.png" width="400" alt="ROI Definition"><br>
      <sub><b>② ROI Definition</b> — draw regions on H&E overview</sub>
    </td>
  </tr>
  <tr>
    <td align="center" width="50%">
      <img src="docs/sample/Operation%20interface/stage3_seg.png" width="400" alt="MCseg v2 Segmentation"><br>
      <sub><b>③ MCseg v2 Segmentation</b> — multi-pass ensemble + preview</sub>
    </td>
    <td align="center" width="50%">
      <img src="docs/sample/Operation%20interface/stage4_count.png" width="400" alt="RNA Counting"><br>
      <sub><b>④ RNA Counting</b> — assign Visium HD bins to cells</sub>
    </td>
  </tr>
  <tr>
    <td align="center" width="50%">
      <img src="docs/sample/Operation%20interface/stage5_umap.png" width="400" alt="UMAP Analysis"><br>
      <sub><b>⑤ UMAP / Leiden</b> — multi-resolution cluster explorer</sub>
    </td>
    <td align="center" width="50%">
      <img src="docs/sample/Operation%20interface/stage5_annotation.png" width="400" alt="Cell-type Annotation"><br>
      <sub><b>⑥ Cell-type Annotation</b> — Celltypist auto-labelling</sub>
    </td>
  </tr>
  <tr>
    <td align="center" width="50%">
      <img src="docs/sample/Operation%20interface/stage6_explore.png" width="400" alt="Spatial Explorer"><br>
      <sub><b>⑦ Spatial Explorer</b> — interactive gene expression viewer</sub>
    </td>
    <td align="center" width="50%">
      <img src="docs/sample/Operation%20interface/stage7_output.png" width="400" alt="Export"><br>
      <sub><b>⑧ Export</b> — Xenium Explorer / Loupe Browser formats</sub>
    </td>
  </tr>
</table>

---

## Example Results

### Cell-type mapping on Visium HD (LUAD, Tumor Boundary ROI)

<p align="center">
  <img src="docs/fig_celltype_map.png" width="700" alt="Cell-type map — ROI1 tumor boundary, n=5,021 cells, Xenium GT validated">
</p>

> 7 cell types resolved by MCseg v2 + Celltypist, validated against Xenium Prime ground-truth (n = 5,021 cells, ROI 1).

### Spatial AT2 Pneumocyte detection overlaid on H&E

<p align="center">
  <img src="docs/fig_spatial_at2.png" width="500" alt="AT2 Pneumocyte (blue outlines, n=326, 30%) on H&E">
</p>

> AT2 Pneumocytes (SFTPC+, blue outlines, n = 326, 30%) detected directly on the H&E image — no GPU required.

### QC filtering (Stage 3)

<p align="center">
  <img src="docs/sample/result/qc_violin.png" width="780" alt="QC violin plots: UMI, genes per cell, % mitochondrial">
</p>

> Violin plots showing per-cell QC metrics after MCseg v2 segmentation — dashed lines indicate configurable thresholds.

### UMAP, marker genes and spatial cell-type map

<table>
  <tr>
    <td align="center" width="33%">
      <img src="docs/sample/result/result_umap.png" width="260" alt="UMAP annotated"><br>
      <sub>UMAP coloured by Celltypist annotation</sub>
    </td>
    <td align="center" width="33%">
      <img src="docs/sample/result/result_dotplot.png" width="260" alt="Marker gene dotplot"><br>
      <sub>Marker gene dotplot per cluster</sub>
    </td>
    <td align="center" width="33%">
      <img src="docs/sample/result/result_spatial_filled_1.png" width="260" alt="Spatial cell-type map"><br>
      <sub>Spatial cell-type map overlaid on H&E</sub>
    </td>
  </tr>
</table>

---

## Usage Guide

After launching (`bash start.sh`), open **http://localhost:3000** and follow the steps below.

### Step 1 — Data Setup

1. Click **Browse** to select your Visium HD sample folder (the root containing `spatial/` and `binned_outputs/`).
2. Click **Scan** — MCseg auto-detects the H&E image (`.btf` / `.tif`), 2 µm and 8 µm binned matrices.
3. Verify that all three files are found (green checkmarks), then click **Apply** to register them.
4. Set the **Output Directory** where results (`roi/`, `analysis/`) will be written, then click **Save**.

> **Data layout expected:**
> ```
> <sample>/
> ├── spatial/
> │   └── tissue_hires_image.btf          ← gigapixel H&E
> └── binned_outputs/
>     ├── square_002um/filtered_feature_bc_matrix/
>     └── square_008um/filtered_feature_bc_matrix/
> ```

### Step 2 — Stage 0: ROI Extraction

1. In the **Add ROI** form, fill in:
   - **Name** — a unique identifier (e.g. `roi1`)
   - **Tissue** — `crc` or `luad` (sets the matching parameter profile for this ROI)
   - **x / y / width / height** — region in full-resolution pixels (1 px = 0.2737 µm)
2. Click **Add** to register the ROI; repeat for all regions of interest.
3. Click **Run ROI Extraction** — MCseg tile-reads the BTF and crops `he_crop.tif` + `adata_002um.h5ad` per ROI.

### Step 3 — Stage 1: MCseg v2 Segmentation

1. Review the default parameters (pre-filled from the tissue profile):

   | Parameter | Default | Notes |
   |-----------|---------|-------|
   | `dia_small / mid / large` | 13 / 17 / 22 px | cyto3 cell diameter sweep |
   | `voronoi_distance` | 9 px | Voronoi expansion cap |
   | `use_hematoxylin` | true | adds H-channel passes |
   | `use_cpsam` | false | enable for complex/dense tissue |
   | `use_transcript_rescue` | true | fills in cells missed by morphology |
   | `use_gpu` | true | MPS / CUDA; falls back to CPU |

2. (Optional) Expand **ROI Overrides** to tune parameters per individual ROI.
3. Click **Preview** on one ROI to verify cell outlines before committing to a full run.
4. Click **Run All ROIs** — outputs `segmentation_masks.npy` per ROI.

### Step 4 — Stage 2: RNA Counting

1. Check the ROI list — each row shows whether a segmentation mask and count result exist.
2. Click **Run All** (or per-ROI **Run**) — each 2 µm bin is assigned to the nearest cell mask with a 6 px dilation.
3. Output: `cellpose_cells.h5ad` (cells × genes sparse matrix).

### Step 5 — Stage 3: Analysis

The analysis stage runs four sequential sub-steps:

| Sub-step | Button | Output |
|----------|--------|--------|
| 1. QC | **Run QC** | QC histograms; filtered cells |
| 2. UMAP | **Run UMAP** | PCA → UMAP → Leiden clusters |
| 3. Heatmap | **Run Heatmap** | Top marker gene heatmap |
| 4. Annotate | **Run Annotate** (Celltypist) | Automated cell-type labels |

Run each sub-step in order; results are visualised inline. Click **Apply Labels** after annotation to write cluster names back to the h5ad.

### Step 6 — Spatial Explorer (`✦`)

Interactive spatial gene expression viewer — available after Stage 3 completes.

1. Select an ROI from the dropdown.
2. Search for a gene or choose a preset panel (Immune/Tumor, Hair Follicle, etc.).
3. Switch between **Contour** (cell outlines) and **Set** (dot overlay) modes.
4. Export the current view as PNG.

### Step 7 — Stage 4: Export

The export page provides both result visualisation and format conversion:

**Visualisation tabs** (review before exporting):

| Tab | Content |
|-----|---------|
| Spatial | Colour-coded cluster map overlaid on H&E |
| UMAP | Dimensionality reduction plot |
| Dotplot | Marker gene expression per cluster |
| Heatmap | Top gene heatmap |

**Export formats:**

| Target | Output | Use for |
|--------|--------|---------|
| Xenium Explorer | GeoJSON cell boundaries + transcript CSV | Spatial visualisation |
| Loupe Browser | Barcode CSV with cluster labels | 10x Genomics Loupe |

Files are saved to `<output_dir>/roi/<roi_name>/export/`.

---

## MCseg v2 Algorithm

```text
1. CLAHE preprocessing (clip=3.0, tile=8×8) + Hematoxylin extraction
2. Multi-pass multi-model detection (4–7 passes depending on options):
   · cyto3 @ 13/17/22 px on CLAHE-RGB (3 passes, always)
   · cyto3 @ 17 px on Hematoxylin channel (1 pass, use_hematoxylin=true by default)
   · cpsam @ auto / 16 px / hematoxylin (up to 3 passes, use_cpsam=false by default)
3. Ensemble merging (IoU overlap threshold < 15%)
4. Voronoi boundary expansion (default d=9 px; d=8 px used in paper benchmark)
5. Quality filtering (20–6000 px²)
```

See [Supplementary Note 1](analysis/supplementary/Supplementary_Note_1.md) for full algorithm specification.

---

## Configuration

All parameters are managed in `config/pipeline.yaml`. Switch tissue type with one line:

```yaml
global:
  tissue_profile: crc   # or: luad
```

---

## Testing

```bash
uv run pytest backend/tests/ -v
```

---

## Troubleshooting

| Issue | Cause | Solution |
| ----- | ----- | -------- |
| `uv: command not found` after install | Shell profile not reloaded | Run `source ~/.zshrc` (zsh) or `source ~/.bashrc` (bash), or restart terminal |
| Backend fails to start (`address in use`) | Previous process still running | `start.sh` auto-kills ports 8001/3000; or run `lsof -ti:8001,3000 \| xargs kill -9` manually |
| `uv sync` fails on ExFAT drive | Resource-fork file corruption | `start.sh` handles this automatically; if running manually: `rm -rf .venv && mkdir -p ~/.venvs/msseg && ln -s ~/.venvs/msseg .venv && uv sync` |
| Out-of-memory during segmentation | ROI too large for available RAM | Reduce ROI size, or decrease `batch_size` (default 4 → try 2 or 1) |
| Slow segmentation | CPU mode | Enable GPU: set `use_gpu: true` in Stage 1 UI or `pipeline.yaml` |
| Too few cells detected | `cellprob_threshold` too high | Lower to `-2.0` or `-3.0` in Stage 1 UI |
| Fragmented small cells | `min_size` too low | Increase `min_size` (e.g., 50 px²) in Stage 1 UI |
| Low bin assignment rate | Voronoi gaps not filled | Set `rna_counting.dilation_px: 6` in `pipeline.yaml` (default is 6) |
| macOS `._*` file errors | ExFAT external drive | Pipeline auto-filters; manually: `find . -name "._*" -delete` |

---

## Citation

If you use MCseg in your research, please cite:

> Chan, C.-R. (詹麒儒), et al. MCseg: End-to-End Visium HD Analysis with AI-Optimised Ensemble Cell Segmentation. *Bioinformatics* (under review), 2026.

---

## Reproducibility

Analysis scripts and data for the paper are provided in the [`analysis/`](analysis/) directory:

```text
analysis/
├── scripts/
│   ├── analysis/     # Core analysis pipeline (01–08)
│   └── figures/      # Figure generation scripts (fig1–fig4, suppfigs)
├── data/             # Per-ROI metrics CSV files
└── supplementary/    # Supplementary Note 1, Table S1, Table S2
```

> **Manuscript**: The full manuscript will be linked here upon publication. Preprint / DOI to be added.

### AI-Autonomous Discovery (AutoResearch)

MCseg v2 was developed by running an AI agent loop over ~80 overnight cycles. The agent iteratively proposed, implemented, and scored segmentation architectures against Xenium ground truth—converging on the multi-model ensemble without human intervention.

Templates for adapting this paradigm to your own segmentation problem are provided in [`docs/autoResearch/`](docs/autoResearch/):

| File | Description |
|------|-------------|
| [`README.md`](docs/autoResearch/README.md) | Overview and adaptation guide |
| [`program.md`](docs/autoResearch/program.md) | Agent task specification template |
| [`segment_template.py`](docs/autoResearch/segment_template.py) | Sandbox starter script (MCseg v2 helpers included) |
| [`run_agent.py`](docs/autoResearch/run_agent.py) | Agent runner using the Anthropic API |

### Data Availability

| Dataset | Source |
|---------|--------|
| LUAD (6 ROIs) | 10x Genomics public demo data + Xenium Prime co-registration |
| CRC (15 ROIs) | 10x Genomics + GEO [GSE280318](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE280318) |

---

## License

MIT License — © 2026 詹麒儒 (Chan Chi Ru). See [LICENSE](LICENSE).
