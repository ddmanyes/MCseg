# MSseg AI Skills

Two portable AI skills for Visium HD spatial transcriptomics analysis.
Copy the `skills/` directory to any analysis folder and use with any AI assistant.

## Skills

| Skill | Purpose | Input | Output |
|-------|---------|-------|--------|
| `msseg-segment.md` | Cell segmentation + QC evaluation | H&E image + Visium HD bins | `results/masks/` + `handoff_report.json` |
| `msseg-analyze.md` | Dimensionality reduction + export | `handoff_report.json` | UMAP, Leiden clusters, Xenium/h5ad export |

## Quick Start

### In MSseg project (uv environment)

```bash
# Check environment and GPU
uv run python scripts/setup_env.py

# Phase 1: Segmentation
# Paste msseg-segment.md as system prompt → run guided flow

# Phase 2: Analysis
# Paste msseg-analyze.md as system prompt → run guided flow
```

### Portable (copy to any analysis folder)

```bash
# Copy skills to your analysis folder
cp -r MSseg/skills/ /your/analysis/folder/

cd /your/analysis/folder
uv run python skills/scripts/setup_env.py
```

Replace `scripts/` with `skills/scripts/` in all commands within the skill guides.

## Directory Structure

```
skills/
├── README.md               # This file
├── msseg-segment.md        # Phase 1 skill: segmentation + QC
├── msseg-analyze.md        # Phase 2 skill: analysis + export
├── IMPLEMENTATION_PLAN.md  # Development history
└── scripts/                # Portable support scripts
    ├── setup_env.py        # Environment check + GPU detection
    ├── roi_sampler.py      # Tissue-aware ROI random sampling
    ├── seg_quality.py      # NUC baseline vs MCseg quality assessment
    ├── qc_metrics.py       # FTC / NED / Co-expression metrics
    ├── write_handoff.py    # Handoff report generator
    ├── build_full_adata.py # Multi-ROI AnnData assembly
    ├── run_analysis.py     # Scanpy QC → UMAP → Leiden CLI
    └── export_mcseg.py     # Xenium zarr+GeoJSON + h5ad export
```

## Requirements

```bash
uv add cellpose scanpy anndata scikit-image zarr tifffile pandas numpy scipy pyarrow pyyaml
```

GPU (recommended):
- Windows (NVIDIA): `uv add torch torchvision --index-url https://download.pytorch.org/whl/cu128`
- macOS: CPU/MPS auto-detected

## Compatible AI Systems

These skills are designed to work with any AI assistant that supports system prompts:
- Claude (claude.ai / Claude Code)
- OpenCode
- Hermes
- Any assistant with code execution capability

## Handoff Protocol

Phase 1 outputs `results/handoff_report.json` with:
- Quality metrics (NED / FTC / Co-expression)
- Recommended QC parameters for Phase 2
- Warnings if quality is below threshold

Phase 2 reads this report at STEP 0 and adjusts parameters accordingly.
