# MCseg — High-Fidelity Visium HD Cell Segmentation

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/)

**MCseg** is a no-code, end-to-end analysis platform for 10x Genomics **Visium HD** (2 µm resolution) spatial transcriptomics data. Its core segmentation engine, **MCseg v2**, uses a 7-pass multi-model ensemble (cyto3 at three diameters + optional cpsam passes) with adaptive Voronoi boundary expansion to achieve high-fidelity cell segmentation without GPU requirements.

> MCseg v2 achieves PQ = 0.554 ± 0.064 on LUAD tissue (vs 0.432 ± 0.037 for single-model Cellpose baseline, +28% relative improvement), validated against Xenium Prime ground-truth masks.

---

## Citation

If you use MCseg in your research, please cite:

> Chan, C.-R. (詹麒儒), et al. MCseg: High-Fidelity Visium HD Cell Segmentation via AI-Autonomous Pipeline Discovery. *Bioinformatics* (under review), 2026.

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

## Installation

### Requirements

| Tool | Version | Notes |
|------|---------|-------|
| macOS | 12+ | Apple Silicon or Intel |
| Python | 3.10+ | Managed by `uv` |
| Node.js | v18+ | For frontend (Vite) |
| GPU | Optional | Apple MPS / NVIDIA CUDA 12.4; auto-falls back to CPU |

### Quick Start

```bash
# 1. Install uv (Python package manager)
curl -LsSf https://astral.sh/uv/install.sh | sh

# 2. Clone and install
git clone https://github.com/ddmanyes/MCseg.git
cd MCseg
uv sync

# 3. Install frontend dependencies
cd frontend && npm install && cd ..

# 4. Launch
bash start.sh
```

Open **<http://localhost:3000>** in your browser.

---

## Pipeline Overview

| Stage | Function | Key Output |
|-------|----------|------------|
| Data Setup | Auto-scan and validate raw data | `state.json` |
| Stage 0: ROI Extract | Crop ROI from Gigapixel BTF | `he_crop.tif`, `adata_002um.h5ad` |
| Stage 1: MCseg v2 | 7-pass ensemble segmentation + Voronoi expansion | `segmentation_masks.npy` |
| Stage 2: RNA Count | Assign Visium HD bins to cells | `cellpose_cells.h5ad` |
| Stage 3: Analysis | QC → Normalise → PCA → UMAP → Leiden | `umap_computed.h5ad` |
| Stage 3.5: Explorer | Interactive spatial gene expression viewer | PNG export |
| Stage 4: Export | Xenium Explorer / Loupe Browser format | GeoJSON, CSV |

---

## MCseg v2 Algorithm

```text
1. CLAHE preprocessing (clip=3.0, tile=8×8) + Hematoxylin extraction
2. 7-pass multi-model detection:
   · cyto3 @ 13/17/22 px on CLAHE-RGB and Hematoxylin channel
   · cpsam @ auto and 16 px (optional)
3. Ensemble merging (overlap threshold < 15%)
4. Voronoi boundary expansion (d=8 px in deployment)
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
| Slow segmentation | CPU mode | Set `use_gpu: true` in `pipeline.yaml` |
| Too few cells detected | `cellprob_threshold` too high | Lower to `-2.0` or `-3.0` |
| Fragmented small cells | `min_size` too low | Increase `min_size` (e.g., 50 px²) |
| Low bin assignment rate | `dilation_px=0` | Set `rna_counting.dilation_px: 6` |
| macOS `._*` file pollution | ExFAT external drive | Pipeline filters these automatically |

---

## License

MIT License — © 2026 詹麒儒 (Chan Chi Ru). See [LICENSE](LICENSE).
