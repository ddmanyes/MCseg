# MCseg：AI 最佳化集成式細胞分割的 Visium HD 空間轉錄組端到端分析平台

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/)

✅ 無需撰碼的網頁介面 · ✅ 從 Gigapixel BTF 自訂 ROI · ✅ 端到端分析（QC → UMAP → 細胞型別標注）· ✅ 多 ROI 合併 · ✅ 互動式空間基因表現探索器 · ✅ 匯出至 Xenium Explorer · ✅ GPU 可選

**MCseg** 是針對 10x Genomics **Visium HD**（2 µm 解析度）空間轉錄組資料的無需撰碼端到端分析平台。從原始 Gigapixel BTF 影像出發，MCseg 涵蓋完整工作流程：自訂 ROI 裁切、高精度細胞分割、RNA 計數、下游分析（QC → UMAP → 細胞型別標注），以及一鍵匯出至 Xenium Explorer 或 Loupe Browser——全程透過網頁介面操作，無需任何程式設計背景。

其核心分割引擎 **MCseg** 採用 **AutoResearch** 範式開發——以 AI 自主架構搜索在約 80 個評估循環中收斂——產出七輪多模型 Cellpose 集成搭配 Voronoi 約束邊界擴張。以 LUAD 組織的 Xenium Prime ground truth 為基準，MCseg 達到 **PQ = 0.554 ± 0.064**——較最佳雙直徑基準線 **2Cseg**（PQ 0.432 ± 0.037）提升 **+28%**。在 CRC 中，MCseg 的轉錄本捕捉能力與 Space Ranger 相當（UMI 密度 11.6 vs 11.7 UMI/µm²），同時保持更高的轉錄邊界純度（NED 0.727 vs 0.712，p = 0.026）。GPU 為可選項，支援完整 CPU 回退。

<p align="center">
  <img src="docs/fig1a_pipeline.png" width="820" alt="MCseg pipeline overview">
</p>

---

## 目錄

[快速開始](#快速開始) · [流程概覽](#流程概覽) · [CLI（無介面模式）](#cli無介面全切片流程) · [介面導覽](#介面導覽) · [範例結果](#範例結果) · [輸出結構](#輸出結構) · [使用指南](#使用指南) · [演算法](#mcseg-演算法) · [設定](#設定) · [疑難排解](#疑難排解) · [引用](#引用) · [授權](#授權)

---

## 快速開始

### 系統需求

| 元件          | 最低需求                              | 建議需求                              | 備註                                                                                                         |
| ------------- | ------------------------------------- | ------------------------------------- | ------------------------------------------------------------------------------------------------------------ |
| **作業系統**  | macOS 12、Windows 10、Ubuntu 20.04    | macOS 13+ / Windows 11 / Ubuntu 22.04 | 三個平台均完整支援                                                                                           |
| **CPU**       | 4 核，任何現代 x86-64 或 ARM          | Apple Silicon（M1/M2/M3）或 AMD/Intel | Apple Silicon → MPS；NVIDIA → CUDA GPU 加速                                                                 |
| **記憶體**    | 8 GB                                  | 16 GB+                                | Cellpose 會將完整 ROI 裁切載入記憶體；超大 ROI（>2000×2000 px）或多 ROI 執行建議 32 GB                      |
| **儲存空間**  | 15 GB 可用                            | 30 GB+ 可用                           | ~8 GB 用於 Python 環境（torch、cellpose）；其餘用於資料與結果                                               |
| **Python**    | 3.10                                  | 3.11                                  | 由 `uv` 管理；請勿使用系統 Python                                                                           |
| **Node.js**   | v18                                   | v20 LTS                               | 前端用（Vite + React）；CLI 模式不需要 Node.js                                                              |
| **GPU**       | —（CPU 回退）                        | NVIDIA CUDA 12.x 或 Apple MPS         | GPU 可大幅縮短分割時間：CPU 約 30 分鐘（4-pass）/ 55 分鐘（7-pass）→ CUDA/MPS 約 5–10 / 15–25 分鐘        |

### 前置需求

**macOS（建議使用 Homebrew）：**

```bash
# 若尚未安裝 Homebrew
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

# 安裝 Node.js（最新 LTS）
brew install node
```

**Linux（Ubuntu/Debian）：**

```bash
curl -fsSL https://deb.nodesource.com/setup_20.x | sudo -E bash -
sudo apt-get install -y nodejs
```

**Windows（以系統管理員身分執行 PowerShell）：**

```powershell
# 安裝 Node.js（winget，Windows 10/11 內建）
winget install OpenJS.NodeJS.LTS

# 安裝 uv
powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex"
```

### 安裝

**macOS / Linux：**

```bash
# 1. 安裝 uv（Python 套件管理器）
curl -LsSf https://astral.sh/uv/install.sh | sh
source ~/.zshrc   # 或重新開啟終端機，確保 uv 在 PATH 中

# 2. 複製並安裝
git clone https://github.com/ddmanyes/MCseg.git
cd MCseg
uv sync           # 若磁碟為 ExFAT，請跳過此步驟——詳見下方說明

# 3. 安裝前端相依套件
cd frontend && npm install && cd ..

# 4. 啟動（同時處理 Python 環境設定）
bash start.sh
```

**Windows（PowerShell）：**

```powershell
# 1. 安裝 uv 後請重新啟動 PowerShell，確保其在 PATH 中

# 2. 複製並安裝
git clone https://github.com/ddmanyes/MCseg.git
cd MCseg

# 若磁碟為 NTFS（C:\、D:\ 等）：
uv sync

# 若磁碟為 ExFAT（外接 SSD，例如 K:\）：
$env:UV_LINK_MODE = "copy"; uv sync

# 3. 安裝前端相依套件
cd frontend; npm install; cd ..

# 4. 啟動後端 + 前端（需兩個終端機視窗）
# 終端機 1：
uv run uvicorn backend.main:app --port 8001
# 終端機 2：
cd frontend; npm run dev
```

在瀏覽器開啟 **[http://localhost:3000](http://localhost:3000)**。

> [!NOTE]
> **bash 使用者：** 將 `source ~/.zshrc` 改為 `source ~/.bashrc`。

> [!IMPORTANT]
> **ExFAT / 外接磁碟使用者（僅限 macOS）：** 跳過步驟 2 的 `uv sync`，直接執行 `bash start.sh`——腳本會在安裝前將 `.venv` 建立為指向 `~/.venvs/msseg`（APFS）的 symlink，避免 resource-fork 損壞。
> **ExFAT / 外接磁碟使用者（Windows）：** 使用 `$env:UV_LINK_MODE = "copy"; uv sync` 而非單純的 `uv sync`，以防止在非 NTFS 磁碟區發生硬連結失敗。若 `uv sync` 後 `.venv` 顯示為 1 KB 檔案，請執行 `cmd /c "attrib -H .venv && del .venv"` 後以 `UV_LINK_MODE=copy` 重新執行。

---

## 流程概覽

| 階段                 | 功能                                                         | 主要輸出                              |
| -------------------- | ------------------------------------------------------------ | ------------------------------------- |
| 資料設定             | 自動掃描並驗證原始資料                                       | `state.json`                         |
| Stage 0：ROI 裁切    | 從 Gigapixel BTF 裁切 ROI                                    | `he_crop.tif`、`adata_002um.h5ad`   |
| Stage 1：MCseg       | 多輪集成分割（4–7 輪）+ Voronoi 擴張                        | `segmentation_masks.npy`             |
| Stage 2：RNA 計數    | 將 Visium HD bins 指派至細胞                                 | `cellpose_cells.h5ad`                |
| Stage 3：分析        | QC → 正規化 → PCA → UMAP → Leiden                          | `umap_computed.h5ad`                 |
| Stage 3.5：探索器    | 互動式空間基因表現檢視器                                     | PNG 匯出                              |
| Stage 4：匯出        | Xenium Explorer / Loupe Browser 格式                         | `experiment.xenium`、zarr archives   |

---

## CLI（無介面）全切片流程

針對批次處理、HPC 叢集或腳本化 pipeline，MSseg 提供**命令列介面（CLI）**，無需開啟網頁介面即可一行指令跑完整全切片流程——**分割 → RNA 計數 → 細胞型態標注**：

1. **裁切** 從原始 BTF 取出 H&E（或載入既有 `he_crop.tif`）
2. **分割** 以 tiled MCseg v2 對整片分割（4-pass，或加 `--cpsam` 為 7-pass）→ `mcseg_mask.npy`
3. **Bin attribution／RNA 計數**（提供 `--tp` + `--h5` 時執行）→ `bin_attribution.parquet`
4. **聚合** cells×genes 矩陣（含細胞重心）→ `cells.h5ad`
5. **CellTypist** 細胞型態標注（未加 `--skip-celltypist` 時）→ `celltypist_labels.csv`（並回寫至 `cells.h5ad`）
6. **Xenium Explorer** 匯出（加 `--export-xenium` 時）→ `xenium_explorer/`

步驟 3–5 在 `--tp`/`--h5` 齊全時自動執行；只給 `--btf`/`--out` 則僅做分割。

### 基本語法

執行 `uv sync` 後，`msseg-segment` 指令即可直接使用：

```powershell
# Windows（PowerShell）— 簡短形式
uv run msseg-segment `
    --btf  "K:\path\to\image.btf" `
    --tp   "K:\path\to\tissue_positions.parquet" `
    --h5   "K:\path\to\filtered_feature_bc_matrix.h5" `
    --out  "K:\path\to\output_dir\" `
    --tissue crc `
    --cpsam
```

```bash
# macOS / Linux — 簡短形式
uv run msseg-segment \
    --btf  "/Volumes/SSD/image.btf" \
    --tp   "/Volumes/SSD/tissue_positions.parquet" \
    --h5   "/Volumes/SSD/filtered_feature_bc_matrix.h5" \
    --out  "/Volumes/SSD/output/" \
    --tissue crc \
    --cpsam
```

> 亦可使用模組形式：`uv run python -m backend.src.cli.segment ...`

### 常見用途

| 任務 | 指令旗標 |
|------|---------|
| **CRC 7-pass**（含 cpsam） | `--tissue crc --cpsam` |
| **LUAD 4-pass**（快速） | `--tissue luad` |
| **跳過 BTF 裁切**（重用現有 he_crop.tif） | `--he-crop path/to/he_crop.tif` |
| **從 BTF 裁切子區域** | `--btf image.btf --crop-y0 4635 --crop-y1 18599 --btf-col0 45752 --btf-col1 55840` |
| **跳過 CellTypist** | `--skip-celltypist` |
| **匯出 Xenium Explorer** | `--export-xenium`（需 `--tp` + `--h5`） |
| **僅使用 CPU** | `--no-gpu` |
| **自訂直徑** | `--dia-small 11 --dia-mid 15 --dia-large 20` |

### 所有選項

```
uv run python -m backend.src.cli.segment --help

  --btf PATH            原始 BigTIFF (.btf) 路徑
  --he-crop PATH        已裁切的 he_crop.tif（略過 BTF 裁切步驟）

  --crop-y0 PX          裁切起始 row（BTF 全圖座標，預設 0）
  --crop-y1 PX          裁切結束 row（-1 = 全圖）
  --btf-col0 PX         裁切起始 col（BTF 全圖座標，預設 0）
  --btf-col1 PX         裁切結束 col（-1 = 全圖）

  --tp PATH             tissue_positions.parquet 路徑
  --h5 PATH             filtered_feature_bc_matrix.h5 路徑
  --out DIR             輸出目錄（必填）

  --tissue {crc,luad,default}   組織 preset（預設 crc）
  --cpsam               啟用 cpsam（7-pass，需更長時間）
  --no-gpu              強制 CPU
  --batch-size N        Cellpose batch size（預設 2）
  --tile-size PX        Tile 大小（預設 1024）
  --overlap PX          Tile 重疊（預設 128）
  --dia-small/mid/large PX      覆寫 cyto3 直徑
  --voronoi-d PX        覆寫 Voronoi 距離
  --cellprob THRESH     覆寫 cellprob_threshold

  --celltypist-model MODEL      CellTypist 模型（預設 Human_Colorectal_Cancer.pkl）
  --skip-celltypist     跳過 CellTypist

  --export-xenium       匯出 Xenium Explorer bundle（需 --tp + --h5）
```

### 輸出檔案

```
<out>/
├── he_crop.tif               ← 裁切後 H&E 影像
├── mcseg_mask.npy            ← MCseg v2 細胞遮罩（int32, H×W）
├── bin_attribution.parquet   ← barcode → cell_id 對應表
├── cells.h5ad                ← cells × genes 矩陣（原始 counts、重心、celltypist 標籤）
├── celltypist_labels.csv     ← cell_id → celltypist_label
└── xenium_explorer/          ← Xenium Explorer bundle（僅 --export-xenium 時）
```

> [!TIP]
> CLI 支援**斷點續跑**：每個輸出檔若已存在則自動跳過該步驟，可隨時中斷後重新執行。

> [!NOTE]
> CLI 和 Web UI 使用**完全相同的 `cellpose_runner.py` 引擎**，參數語義一致。Web UI 做的任何 ROI 參數覆寫都可以直接翻譯成 `--dia-mid` / `--voronoi-d` 等 CLI 旗標。

---

## 介面導覽

<table>
  <tr>
    <td align="center" width="50%">
      <img src="docs/sample/Operation%20interface/stage1_steup.png" width="400" alt="資料設定"><br>
      <sub><b>① 資料設定</b> — 掃描 BTF + 分箱矩陣，設定輸出目錄</sub>
    </td>
    <td align="center" width="50%">
      <img src="docs/sample/Operation%20interface/stage2_ROI.png" width="400" alt="ROI 定義"><br>
      <sub><b>② ROI 定義</b> — 在 H&E 概覽圖上標記區域</sub>
    </td>
  </tr>
  <tr>
    <td align="center" width="50%">
      <img src="docs/sample/Operation%20interface/stage3_seg.png" width="400" alt="MCseg 分割"><br>
      <sub><b>③ MCseg 分割</b> — 多輪集成分割 + 預覽</sub>
    </td>
    <td align="center" width="50%">
      <img src="docs/sample/Operation%20interface/stage4_count.png" width="400" alt="RNA 計數"><br>
      <sub><b>④ RNA 計數</b> — 將 Visium HD bins 指派至細胞</sub>
    </td>
  </tr>
  <tr>
    <td align="center" width="50%">
      <img src="docs/sample/Operation%20interface/stage5_umap.png" width="400" alt="UMAP 分析"><br>
      <sub><b>⑤ UMAP / Leiden</b> — 多解析度叢集探索器</sub>
    </td>
    <td align="center" width="50%">
      <img src="docs/sample/Operation%20interface/stage5_annotation.png" width="400" alt="細胞型別標注"><br>
      <sub><b>⑥ 細胞型別標注</b> — Celltypist 自動標記</sub>
    </td>
  </tr>
  <tr>
    <td align="center" width="50%">
      <img src="docs/sample/Operation%20interface/stage6_explore.png" width="400" alt="空間探索器"><br>
      <sub><b>⑦ 空間探索器</b> — 互動式基因表現檢視器</sub>
    </td>
    <td align="center" width="50%">
      <img src="docs/sample/Operation%20interface/stage7_output.png" width="400" alt="匯出"><br>
      <sub><b>⑧ 匯出</b> — 空間分析結果；Xenium Explorer / Loupe Browser 匯出於同一頁面</sub>
    </td>
  </tr>
</table>

---

## 範例結果

### Visium HD 細胞型別空間分佈圖（LUAD，腫瘤邊界 ROI）

<p align="center">
  <img src="docs/fig2h.png" width="700" alt="細胞型別空間圖 — LUAD 腫瘤邊界，MCseg + Celltypist">
</p>

> MCseg + Celltypist 在 LUAD 腫瘤邊界 ROI 解析的細胞型別——T/B 淋巴球、Club 上皮細胞、漿細胞、B 細胞、SPP1⁺ 巨噬細胞疊加於 H&E 影像。

### AT2 肺泡細胞空間偵測疊加於 H&E

<p align="center">
  <img src="docs/fig_spatial_at2.png" width="500" alt="AT2 肺泡細胞（藍色輪廓，n=326，30%）疊加於 H&E">
</p>

> AT2 肺泡細胞（SFTPC+，藍色輪廓，n = 326，30%）直接偵測於 H&E 影像——無需 GPU。

### CRC 中的轉錄本歸因

> 在 CRC（15 個 ROI）中，MCseg 的每細胞 RNA 捕捉量與 Space Ranger 相當（UMI 密度 11.6 vs 11.7 UMI/µm²），同時達到更高的轉錄邊界純度（NED 0.727 vs 0.712，p = 0.026）。在三級淋巴組織中，MCseg 解析出**四**個功能性免疫細胞群，而 Space Ranger 僅三個，且細胞數多 44%（636 vs 440）。

### QC 過濾（Stage 3）

<p align="center">
  <img src="docs/sample/result/qc_violin.png" width="780" alt="QC 小提琴圖：UMI、每細胞基因數、粒線體比例">
</p>

> MCseg 分割後每細胞 QC 指標的小提琴圖——虛線表示可設定的過濾閾值。

### UMAP、標記基因與空間細胞型別圖

<table>
  <tr>
    <td align="center" width="25%">
      <img src="docs/sample/result/result_umap.png" width="200" alt="UMAP 標注圖"><br>
      <sub>以 Celltypist 標注著色的 UMAP</sub>
    </td>
    <td align="center" width="25%">
      <img src="docs/sample/result/result_dotplot.png" width="200" alt="標記基因點圖"><br>
      <sub>各叢集標記基因點圖</sub>
    </td>
    <td align="center" width="25%">
      <img src="docs/sample/result/result_heatmap.png" width="200" alt="頂級標記基因熱圖"><br>
      <sub>頂級標記基因熱圖</sub>
    </td>
    <td align="center" width="25%">
      <img src="docs/sample/result/result_spatial_filled_1.png" width="200" alt="空間細胞型別圖"><br>
      <sub>細胞型別空間圖疊加於 H&E</sub>
    </td>
  </tr>
</table>

### 匯出至 Xenium Explorer

MCseg 輸出可直接載入的 Xenium Explorer 套件（`experiment.xenium` + zarr archives）。下方截圖顯示 MCseg 匯出後直接在 Xenium Explorer 4.1.1 中載入的 CRC 資料。

<table>
  <tr>
    <td align="center" width="50%">
      <img src="docs/sample/result/xenium_capture/xenium_2.png" width="400" alt="Xenium Explorer 中的 H&E 影像與 MCseg 細胞邊界"><br>
      <sub>H&E 影像與 MCseg 細胞邊界</sub>
    </td>
    <td align="center" width="50%">
      <img src="docs/sample/result/xenium_capture/xenium_5.png" width="400" alt="Xenium Explorer 中的細胞型別標注群組"><br>
      <sub>細胞型別群組（Celltypist）— 互動式細胞資訊彈窗</sub>
    </td>
  </tr>
  <tr>
    <td align="center" width="50%">
      <img src="docs/sample/result/xenium_capture/xenium_3.png" width="400" alt="Xenium Explorer 中的轉錄本點圖視覺化"><br>
      <sub>轉錄本點圖疊加（SCGB1A1）</sub>
    </td>
    <td align="center" width="50%">
      <img src="docs/sample/result/xenium_capture/xenium_4.png" width="400" alt="Xenium Explorer 中的基因特異性轉錄本密度"><br>
      <sub>基因特異性轉錄本密度</sub>
    </td>
  </tr>
</table>

> **匯出套件結構**（`<output_dir>/export/`）：
>
> ```
> experiment.xenium
> morphology.ome.tif
> cells.zarr.zip
> transcripts.zarr.zip
> cell_feature_matrix.zarr.zip
> analysis.zarr.zip
> analysis_summary.html
> ```

---

## 輸出結構

完整執行後，輸出目錄將包含：

```text
<output_dir>/
├── analysis/
│   ├── roi/
│   │   └── {roi_name}/
│   │       ├── he_crop.tif                  ← H&E 裁切圖（Stage 0）
│   │       ├── adata_002um.h5ad             ← 2 µm 分箱矩陣（Stage 0）
│   │       ├── segmentation_masks.npy       ← MCseg 細胞遮罩（Stage 1）
│   │       ├── segmentation_masks.tif       ← 視覺化疊加圖（Stage 1）
│   │       ├── cellpose_cells.h5ad          ← 細胞 × 基因矩陣（Stage 2）
│   │       ├── cellpose_polygons.json       ← 細胞邊界多邊形（Stage 2）
│   │       └── transcripts_roi.csv          ← 每細胞轉錄本表（Stage 2）
│   ├── merged_all_rois.h5ad                 ← 多 ROI 合併 AnnData（Stage 3，合併模式）
│   ├── qc_preprocessed.h5ad                ← QC 後 AnnData（Stage 3）
│   ├── umap_computed.h5ad                   ← UMAP + Leiden 叢集（Stage 3）
│   ├── combined_cellpose_polygons.json      ← 合併多邊形（Stage 3，合併模式）
│   └── combined_transcripts.csv            ← 合併轉錄本（Stage 3，合併模式）
└── export/
    └── xenium/
        └── {roi_name}/
            ├── experiment.xenium            ← 在 Xenium Explorer 中載入此檔案
            ├── morphology.ome.tif
            ├── cells.zarr.zip
            ├── transcripts.zarr.zip
            ├── cell_feature_matrix.zarr.zip
            ├── analysis.zarr.zip
            └── analysis_summary.html
```

---

## 使用指南

啟動後（`bash start.sh`），開啟 **[http://localhost:3000](http://localhost:3000)** 並按照以下步驟操作。

> *以下時間為估計值，以 **Apple M2 CPU、16 GB RAM**、ROI 約 1500 × 1200 px 為測量基準。GPU（Apple MPS 或 NVIDIA CUDA）可將 Stage 1 縮短至 ~2–3 分鐘/ROI。*

### 步驟一 — 資料設定

1. 點選 **Browse** 選取 Visium HD 樣本資料夾（包含 `spatial/` 和 `binned_outputs/` 的根目錄）。
2. 點選 **Scan**——MCseg 自動偵測 H&E 影像（`.btf` / `.tif`）、2 µm 和 8 µm 分箱矩陣。
3. 確認三個檔案均已找到（綠色勾選），點選 **Apply** 完成登記。
4. 設定結果（`roi/`、`analysis/`）的**輸出目錄**，點選 **Save**。

> **預期資料目錄結構：**
>
> ```
> <sample>/
> ├── spatial/
> │   └── tissue_hires_image.btf          ← Gigapixel H&E
> └── binned_outputs/
>     ├── square_002um/filtered_feature_bc_matrix/
>     └── square_008um/filtered_feature_bc_matrix/
> ```

### 步驟二 — Stage 0：ROI 裁切（~1 分鐘/ROI）

1. 在 **Add ROI** 表單填寫：
   - **名稱** — 唯一識別碼（例如 `roi1`）
   - **組織類型** — `crc` 或 `luad`（為此 ROI 套用對應的參數設定檔）
   - **x / y / width / height** — 全解析度像素座標（1 px = 0.2737 µm）
2. 點選 **Add** 登記 ROI；重複以上步驟新增所有感興趣區域。
3. 點選 **Run ROI Extraction**——MCseg 以 tile-by-tile 方式讀取 BTF，裁切每個 ROI 的 `he_crop.tif` + `adata_002um.h5ad`。

### 步驟三 — Stage 1：MCseg 分割（CPU ~30 分鐘/ROI · GPU ~2–3 分鐘 · 預設 4-pass）

1. 確認預設參數（從組織設定檔預填）：
   | 參數                        | 預設值          | 說明                                                              |
   | --------------------------- | --------------- | ----------------------------------------------------------------- |
   | `dia_small / mid / large` | 13 / 17 / 22 px | cyto3 細胞直徑掃描範圍                                            |
   | `voronoi_distance`        | 9 px            | Voronoi 擴張上限                                                  |
   | `use_hematoxylin`         | true            | 加入蘇木精通道輪次                                                |
   | `use_cpsam`               | false           | 啟用以處理複雜/緻密組織（+3 輪，CPU 約 50–60 分鐘）              |
   | `use_transcript_rescue`   | true            | 補救形態學遺漏的細胞                                              |
   | `use_gpu`                 | true            | MPS / CUDA；自動回退至 CPU                                        |

   啟用 `use_cpsam` 時，cpsam 7-pass 規格可獨立微調（論文 Pass 5/6/7）：

   | 參數                    | 預設     | 說明                                  |
   | ----------------------- | -------- | ------------------------------------- |
   | `dia_cpsam_auto`      | 0（auto）| Pass 5/7 直徑；0 = Cellpose 自動（~30 px）|
   | `dia_cpsam_small`     | 16 px    | Pass 6 固定直徑                       |
   | `cellprob_cpsam_auto` | -1.0     | Pass 5（CLAHE-RGB, auto dia）         |
   | `cellprob_cpsam_small`| -3.0     | Pass 6（CLAHE-RGB, dia=16）           |
   | `cellprob_cpsam_hema` | -1.0     | Pass 7（蘇木精, auto dia）            |

2. （選擇性）展開 **ROI Overrides** 為個別 ROI 微調參數（含上述全部參數與 cpsam 7-pass 規格）。
3. 對其中一個 ROI 點選 **Preview**，在正式執行前確認細胞輪廓。
4. 點選 **Run All ROIs**——每個 ROI 輸出 `segmentation_masks.npy`。

> **整片分割（不分 ROI）：** **Run Full Segmentation** 會以 tiled MCseg v2 對整片切片分割（MPS 安全：tile=1024、batch≤2、停用 cpsam），輸出 `full_image_segmentation_masks.npy`。內建 6 GB 記憶體上限保護過大切片——超過時請改用 ROI 模式（或 [CLI](#cli無介面全切片流程)，以 tile 方式讀取 BTF、無此上限）。此模式只產生遮罩，計數與分析在 UI 上仍以 ROI 為單位。

### 步驟四 — Stage 2：RNA 計數（~2–3 分鐘/ROI）

1. 查看 ROI 清單——每行顯示分割遮罩與計數結果是否存在。
2. 點選 **Run All**（或個別 ROI 的 **Run**）——每個 2 µm bin 以 6 px 膨脹指派至最近的細胞遮罩。
3. 輸出：`cellpose_cells.h5ad`（細胞 × 基因稀疏矩陣）。

### 步驟五 — Stage 3：分析（~3–5 分鐘）

分析階段包含四個依序執行的子步驟：

| 子步驟      | 按鈕                                | 輸出                           |
| ----------- | ----------------------------------- | ------------------------------ |
| 1. QC       | **Run QC**                    | QC 直方圖；過濾後細胞          |
| 2. UMAP     | **Run UMAP**                  | PCA → UMAP → Leiden 叢集      |
| 3. 熱圖     | **Run Heatmap**               | 頂級標記基因熱圖               |
| 4. 標注     | **Run Annotate**（Celltypist）| 自動細胞型別標記               |

依序執行每個子步驟；結果即時顯示於介面。標注完成後點選 **Apply Labels**，將叢集名稱寫回 h5ad 檔案。

### 步驟六 — 空間探索器（`✦`）

互動式空間基因表現檢視器——Stage 3 完成後可使用。

1. 從下拉選單選取 ROI。
2. 搜尋基因，或選擇預設面板（免疫/腫瘤、毛囊等）。
3. 切換 **Contour**（細胞輪廓）和 **Set**（點圖疊加）模式。
4. 將目前檢視匯出為 PNG。

### 步驟七 — Stage 4：匯出（~2–5 分鐘/ROI）

匯出頁面同時提供結果視覺化與格式轉換：

**視覺化頁籤**（匯出前確認）：

| 頁籤    | 內容                           |
| ------- | ------------------------------ |
| 空間圖  | 叢集顏色圖疊加於 H&E           |
| UMAP    | 降維圖                         |
| 點圖    | 每叢集標記基因表現             |
| 熱圖    | 頂級基因熱圖                   |

**匯出格式：**

| 目標            | 輸出                                                         | 用途                                |
| --------------- | ------------------------------------------------------------ | ----------------------------------- |
| Xenium Explorer | Xenium 原生套件（`experiment.xenium` + zarr archives）     | 直接在 Xenium Explorer 4+ 中載入    |
| Loupe Browser   | `.cloupe` 檔案 + 含叢集標記的 barcode CSV                  | 10x Genomics Loupe Browser          |

檔案儲存至 `<output_dir>/export/xenium/{roi_name}/`。

---

## MCseg 演算法

```text
1. CLAHE 前處理（clip=3.0, tile=8×8）+ 蘇木精提取
2. 多輪多模型偵測（4–7 輪，依選項而定）：
   · cyto3 @ 13/17/22 px，以 CLAHE-RGB 為輸入（3 輪，固定）
   · cyto3 @ 17 px，以蘇木精通道為輸入（1 輪，use_hematoxylin=true，預設啟用）
   · cpsam @ auto / 16 px / 蘇木精（最多 3 輪，use_cpsam=false，預設停用）
3. 集成合併（IoU 重疊閾值 < 15%）
4. Voronoi 邊界擴張（預設 d=9 px；論文基準測試使用 d=8 px）
5. 品質過濾（20–6000 px²）
```

完整演算法規格請見 [Supplementary Note 1](analysis/supplementary/Supplementary_Note_1.md)。

---

## 設定

所有參數均在 `config/pipeline.yaml` 中管理。切換組織類型只需修改一行：

```yaml
global:
  tissue_profile: crc   # 或：luad
```

---

## 測試

```bash
uv sync --extra dev            # 安裝 pytest-asyncio + httpx（API 測試所需）
uv run pytest backend/tests/ -v
```

> **ExFAT／外接硬碟（macOS）：** `uv run` 會重建環境並可能覆蓋 `.venv` symlink。
> 請先清除 resource-fork 雜訊，再直接用 venv 跑 pytest：
>
> ```bash
> find . -name '._*' -delete && find ~/.venvs/msseg -name '._*' -delete
> .venv/bin/python -m pytest backend/tests/ -v
> ```

---

## 疑難排解

| 問題                                        | 原因                            | 解決方法                                                                                                                                        |
| ------------------------------------------- | ------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------- |
| 安裝後 `uv: command not found`            | Shell 設定檔未重新載入          | 執行 `source ~/.zshrc`（zsh）或 `source ~/.bashrc`（bash），或重新啟動終端機                                                                 |
| 後端啟動失敗（`address in use`）          | 先前的程序仍在執行              | `start.sh` 自動終止 8001/3000 埠；或手動執行 `lsof -ti:8001,3000 \| xargs kill -9`                                                          |
| ExFAT 磁碟上 `uv sync` 失敗               | Resource-fork 檔案損壞          | `start.sh` 自動處理；手動執行：`rm -rf .venv && mkdir -p ~/.venvs/msseg && ln -s ~/.venvs/msseg .venv && uv sync`                           |
| 分割時記憶體不足                            | ROI 超出可用 RAM                | 縮小 ROI，或調低 `batch_size`（預設 4 → 試試 2 或 1）                                                                                       |
| 分割速度過慢                                | CPU 模式                        | 啟用 GPU：在 Stage 1 UI 或 `pipeline.yaml` 中設定 `use_gpu: true`                                                                            |
| 偵測到的細胞過少                            | `cellprob_threshold` 過高     | 在 Stage 1 UI 中降低至 `-2.0` 或 `-3.0`                                                                                                     |
| 細胞分裂成過多小片段                        | `min_size` 過低               | 在 Stage 1 UI 中提高 `min_size`（例如 50 px²）                                                                                               |
| Bin 指派率偏低                              | Voronoi 間隙未填補              | 在 `pipeline.yaml` 中設定 `rna_counting.dilation_px: 6`（預設即為 6）                                                                        |
| CLI：Windows 上的 `.venv` 檔案錯誤        | 來自 macOS 的 ExFAT symlink     | 執行 `cmd /c "attrib -H K:\...\MSseg\.venv && del K:\...\MSseg\.venv"` 後執行 `$env:UV_LINK_MODE="copy"; uv sync`                          |
| CLI：`zarr < 3 not supported`             | tifffile 版本衝突               | 在 MSseg venv 內執行 `uv pip install "tifffile==2023.12.9"`                                                                                 |
| macOS `._*` 檔案錯誤                      | ExFAT 外接磁碟                  | Pipeline 自動過濾；手動清除：`find . -name "._*" -delete`                                                                                     |

---

## 引用

若您在研究中使用 MCseg，請引用：

> Chan, C.-R.\*, Chang, N.-W.\*, Wang, C.-Y., Tan, H.-Y.†, Lin, S.-J.† MCseg: End-to-end Visium HD spatial transcriptomics analysis with AI-optimised ensemble-based cell segmentation. *Bioinformatics*（審稿中），2026。

---

## 可重現性

論文相關分析腳本與資料存放於 [`analysis/`](analysis/) 目錄：

```text
analysis/
├── scripts/
│   ├── analysis/     # 核心分析 pipeline（01–08）
│   └── figures/      # 圖表生成腳本（fig1–fig4, suppfigs）
├── data/             # 每 ROI 指標 CSV 檔案
└── supplementary/    # Supplementary Note 1、Table S1、Table S2
```

> **論文全文**：出版後將在此提供連結。預印本 / DOI 待補。

### AI 自主探索（AutoResearch）

MCseg 採用 **AutoResearch** 範式開發（[Karpathy, 2026](https://github.com/karpathy/autoresearch)）——一個 AI 自主架構搜索框架，由 AI 代理在約 80 個循環中反覆提案、實作並評分完整分割 pipeline，以 Xenium ground truth 為評估基準，在無人工介入的情況下收斂至多模型集成方案。候選架構使用 Anthropic Claude API（`claude-sonnet-4-5`）評估。就我們所知，MCseg 是首個透過 AI 自主架構搜索開發的細胞分割方法。

用於將此範式調適至您自身分割問題的模板存放於 [`docs/autoResearch/`](docs/autoResearch/)：

| 檔案                                                          | 說明                                       |
| ------------------------------------------------------------- | ------------------------------------------ |
| [`README.md`](docs/autoResearch/README.md)                     | 概覽與調適指南                             |
| [`program.md`](docs/autoResearch/program.md)                   | 代理任務規格模板                           |
| [`segment_template.py`](docs/autoResearch/segment_template.py) | 沙箱起始腳本（含 MCseg 輔助函數）          |
| [`run_agent.py`](docs/autoResearch/run_agent.py)               | 使用 Anthropic API 的代理執行程式          |

### 資料可用性

| 資料集        | 來源                                                                                     |
| ------------- | ---------------------------------------------------------------------------------------- |
| LUAD（6 ROI） | 10x Genomics 公開展示資料 + Xenium Prime 共配準                                          |
| CRC（15 ROI） | 10x Genomics + GEO [GSE280318](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE280318) |

---

## 授權

MIT 授權條款 — © 2026 詹麒儒（Chan Chi Ru）。詳見 [LICENSE](LICENSE)。
