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
├── supplementary/    # Supplementary Note 1, Table S1, Table S2
└── manuscript/       # Manuscript source (manuscript.md)
```

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

## License

MIT License — © 2026 詹麒儒 (Chan Chi Ru). See [LICENSE](LICENSE).

---

## 中文說明

MCseg 是專為 10x Genomics **Visium HD**（2µm 解析度）空間轉錄體學資料設計的全流程分析平台。核心分割引擎採用 **MCseg v2**——以 cyto3 多直徑集成（13/17/22px）搭配可選 Hematoxylin 通道 pass 與 Voronoi 擴張，取代傳統單模型雙尺寸策略，大幅提升複雜腫瘤微環境的細胞邊界精度（LUAD PQ=0.554 vs cellpose_dilate 0.432，+28%）。

架構採用 **FastAPI 後端**（port 8001）搭配 **React + Vite 前端**（port 3000），支援瀏覽器內視覺化操作與 WebSocket 即時日誌追蹤。

### 安裝與啟動

請參閱上方 Quick Start 步驟（中英相同）。

### 設定系統

所有分析參數集中於 `config/pipeline.yaml`。支援**組織 Profile 系統**：

| Profile | 適用組織 | TME panels |
|---------|---------|------------|
| `crc.yaml` | 大腸直腸癌 | 8 panels |
| `luad.yaml` | 肺腺癌 | 4 panels |

### MCseg v2 關鍵參數

```yaml
segmentation:
  mcseg_v2:
    dia_small: 13.0
    dia_mid: 17.0
    dia_large: 22.0
    use_hematoxylin: true
    use_cpsam: false
    voronoi_distance: 9
    min_size: 20
    max_size: 6000
    flow_threshold: 0.4
    cellprob_threshold: -2.0
    clahe_clip_limit: 3.0
```

### 常見問題

| 問題 | 原因 | 解法 |
|------|------|------|
| MCseg v2 速度慢 | CPU 模式 | 設定 `use_gpu: true` |
| 細胞數偏少 | cellprob_threshold 過高 | 降低至 -2.0 或 -3.0 |
| 細胞過小/碎片化 | min_size 過低 | 提高 `min_size`（如 50 px²） |
| bins 指派率低 | dilation_px=0 | 設定 `rna_counting.dilation_px: 6` |
| macOS `._*` 污染 | ExFAT 外接硬碟 | Pipeline 已內建自動過濾 |
