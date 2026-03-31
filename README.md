# MSseg

MSseg 是專為 10x Genomics **Visium HD**（2µm 解析度）空間轉錄體學資料設計的全流程分析平台。核心分割引擎採用 **MCseg v2**——以 cyto3 多直徑集成（13/17/22px）搭配可選 Hematoxylin 通道 pass 與 Voronoi 擴張，取代傳統單模型雙尺寸策略，大幅提升複雜腫瘤微環境的細胞邊界精度（LUAD PQ=0.554 vs cellpose_dilate 0.432，+28%）。

架構採用 **FastAPI 後端**（port 8001）搭配 **React + Vite 前端**（port 3000），支援瀏覽器內視覺化操作與 WebSocket 即時日誌追蹤。

---

## 安裝與啟動

### 系統需求

| 工具 | 版本 | 備註 |
| --- | --- | --- |
| macOS | 12 以上 | Apple Silicon（MPS）或 Intel 均可 |
| Python | 3.10 以上 | 由 uv 自動管理，無需手動安裝 |
| Node.js | v18 以上 | 供前端 Vite 使用 |
| GPU | 選配 | Apple MPS / NVIDIA CUDA 12.4，無 GPU 自動回退 CPU |

---

### 步驟 1：安裝 uv（Python 套件管理器）

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```

安裝後重新開啟終端機，確認安裝成功：

```bash
uv --version
```

> 若已安裝可跳過此步驟。

---

### 步驟 2：安裝 Node.js（供前端使用）

建議使用 [nvm](https://github.com/nvm-sh/nvm) 管理版本：

```bash
# 安裝 nvm
curl -o- https://raw.githubusercontent.com/nvm-sh/nvm/v0.39.7/install.sh | bash

# 重新開啟終端機後安裝 Node.js 22
nvm install 22
nvm use 22
```

或直接從 [nodejs.org](https://nodejs.org) 下載安裝 v18 以上版本。

確認安裝成功：

```bash
node --version   # 應顯示 v18.x 以上
npm --version
```

---

### 步驟 3：下載專案

```bash
git clone <repository-url> MSseg
cd MSseg
```

---

### 步驟 4：安裝 Python 依賴

```bash
uv sync
```

`uv sync` 會自動讀取 `pyproject.toml`，建立 `.venv` 並安裝所有套件（包含 Cellpose、PyTorch、Scanpy 等），**無需手動 pip install**。

> **外接硬碟（ExFAT）注意**：ExFAT 格式的外接硬碟不支援 symlink，建議直接執行下方的一鍵啟動腳本，它會自動將 `.venv` 移至本機 `~/.venvs/msseg` 以避免問題。

---

### 步驟 5：安裝前端依賴

```bash
cd frontend
npm install
cd ..
```

---

### 啟動方式

#### 方式 A：一鍵啟動（推薦）

```bash
bash start.sh
```

腳本會自動：
- 處理 `.venv` ExFAT 問題（symlink 至本機）
- 清除佔用 port 8001/3000 的舊行程
- 依序啟動後端（等待健康確認）→ 前端
- `Ctrl+C` 同時停止所有服務

#### 方式 B：分開啟動（開發模式）

```bash
# 終端機 1：後端
uv run uvicorn backend.main:app --reload --port 8001

# 終端機 2：前端
cd frontend && npm run dev
```

---

啟動後前往 **<http://localhost:3000>**

| 服務 | 網址 |
| --- | --- |
| 前端 UI | http://localhost:3000 |
| 後端 API | http://localhost:8001 |
| API 文件 | http://localhost:8001/docs |

---

## 設定系統

所有分析參數集中於 `config/pipeline.yaml`，**禁止在程式碼中硬編碼任何參數**。支援**組織 Profile 系統**，換組織只需修改一行。

### Tissue Profile（多組織支援）

```yaml
# pipeline.yaml
global:
  tissue_profile: crc    # 切換至 luad 只需改這一行
```

| Profile | MCseg v2 dia 範圍 | TME panels | 適用組織 |
| --- | --- | --- | --- |
| `crc.yaml` | 13 / 17 / 22 px | 8 panels | 大腸直腸癌 |
| `luad.yaml` | 13 / 17 / 22 px | 4 panels | 肺腺癌 |

三層合併順序（後者覆寫前者）：

```text
config/profiles/{tissue_profile}.yaml  ← 組織基底
        ↓
config/pipeline.yaml                   ← 專案覆寫
        ↓
results/state.json                     ← 執行期覆寫（UI 互動）
```

### MCseg v2 關鍵參數

```yaml
segmentation:
  mcseg_v2:
    # GPU / 批次
    use_gpu: true
    batch_size: 4
    # cyto3 多直徑集成
    dia_small: 13.0          # small cells（淋巴細胞等）
    dia_mid: 17.0            # 主要 pass（集成基底）
    dia_large: 22.0          # 大型細胞（上皮、巨噬細胞）
    # 可選 pass
    use_hematoxylin: true    # cyto3 on Ruifrok Hematoxylin 通道
    use_cpsam: false         # cpsam 集成（可選，顯著增加時間）
    # Voronoi 擴張
    voronoi_distance: 9      # px，約 2.5µm
    # 細胞過濾
    min_size: 20             # px²
    max_size: 6000           # px²
    # Cellpose 品質控制
    flow_threshold: 0.4
    cellprob_threshold: -2.0
    clahe_clip_limit: 3.0
    # 轉錄本密度補救（選配）
    use_transcript_rescue: true
```

---

## 流程說明

MSseg 共 6 個步驟。左側選單顯示進度，前一步驟完成前後續步驟鎖定。

---

### 資料設定 (Data Setup)

**目的**：自動掃描並驗證原始資料完整性。

後端 `discovery.py` 掃描指定的組織樣本根目錄，自動識別：
- H&E 影像：OME-TIFF BTF（BigTIFF 格式，多解析度金字塔）
- Visium HD bins：`square_002um/filtered_feature_bc_matrix.h5` + `tissue_positions.parquet`
- 空間資訊：`square_002um/spatial/scalefactors_json.json`

結果寫入 `results/state.json`，後續所有步驟從此讀取路徑。

---

### Stage 0: ROI 裁切 (ROI Extract)

**目的**：從 Gigapixel BTF 中精準裁切感興趣區域，避免全圖載入耗盡記憶體。

頁面內嵌 **OpenSeadragon** 多解析度瀏覽器，直接從 BTF 動態載入 DZI tile，支援流暢縮放與平移。

**後端原理**（`extractor.py`）：
1. tile-by-tile 讀取 BTF，記憶體峰值 < 2 GB
2. 篩選 ROI 框內 `in_tissue == 1` 的 bins
3. 稀疏矩陣 slice 輸出 `adata_002um.h5ad`

**輸出**：
- `results/roi/{ROI_NAME}/he_crop.tif`
- `results/roi/{ROI_NAME}/adata_002um.h5ad`

---

### Stage 1: 細胞分割 (MCseg v2)

**目的**：以 MCseg v2 多模型 Voronoi 集成，在 H&E 影像上標定每個細胞的像素邊界。

#### MCseg v2 演算法

```text
1. CLAHE 前處理（LAB 色彩空間局部對比增強）
2. 多 pass Cellpose cyto3 推論：
   ・pass 1：dia_small=13px（淋巴細胞、小細胞核）
   ・pass 2：dia_mid=17px（主要推論，集成基底）
   ・pass 3：dia_large=22px（腸上皮、巨噬細胞）
   ・pass 4（選配）：cyto3 on Ruifrok Hematoxylin 通道
   ・pass 5-7（選配）：cpsam×3 補充偵測
3. merge_masks_fast 集成（IoU 閾值去重）
4. Voronoi 擴張（max_distance=voronoi_distance px）
   ─ 每個背景像素指派至最近細胞
   ─ 無重疊保證（vs expand_labels 可能重疊）
5. 轉錄本密度補救（選配）
   ─ 從 vhd_pseudo_transcripts.csv 高密度區域植入補漏細胞
6. 尺寸過濾（min_size / max_size px²）
```

**效能**（LUAD，n=6 ROI）：

| 指標 | cellpose_dilate | MCseg v2 | 提升 |
| --- | :---: | :---: | :---: |
| PQ@0.5（mean） | 0.432 | **0.554** | **+28%** |
| SQ | — | 0.777 | — |
| RQ | — | 0.711 | — |

**每 ROI 獨立參數覆寫**：Stage 1 UI 提供 ROI 級別的 `dia_mid`、`voronoi_distance` 等欄位覆寫，無需修改全域設定。

**輸出**：
- `results/roi/{ROI_NAME}/segmentation_masks.npy`（H×W int32，0=背景）
- `results/roi/{ROI_NAME}/segmentation_masks.tif`

---

### Stage 2: RNA 計數 (RNA Count)

**目的**：將 Visium HD 2µm bins 的 RNA 計數依細胞遮罩分配至細胞層級。

**稀疏矩陣法**（`counter.py`）：

```text
1. Dilation（expand_labels，dilation_px=6 = 1.64 µm）
   ─ 填補 Voronoi 擴張後的細胞間隙 bins
2. Bin 座標映射至遮罩像素
3. 稀疏矩陣聚合：A @ adata_002um.X（純矩陣乘法，無迴圈）
```

速度：LUAD 837,530 bins < 30 秒。

**輸出**：
- `results/roi/{ROI_NAME}/cellpose_cells.h5ad`（cells × 18K genes）

---

### Stage 3: 下游分析 (Analysis)

**Scanpy 流程**：

```text
原始計數矩陣
→ QC 過濾（min_genes / max_genes / max_pct_mito / min_complexity）
→ Normalize（target_sum=10,000）+ log1p
→ HVG 篩選（top 2,000，seurat_v3）
→ PCA（n_components=20）
→ kNN Graph（n_neighbors=15）
→ UMAP（min_dist=0.3）
→ Leiden 聚類（resolution=0.5）
→ Marker Gene 計算
→ TME Panel 分析（由 tissue profile YAML 定義）
```

**輸出**：
- `results/roi/{ROI_NAME}/qc_preprocessed.h5ad`
- `results/roi/{ROI_NAME}/umap_computed.h5ad`
- `figures/{ROI_NAME}/`：UMAP、Violin、Dotplot（300 DPI PNG）

---

### Stage 4: Browser 匯出 (Export)

**Xenium Explorer**：
- `skimage.measure.find_contours` 從遮罩提取細胞邊界多邊形（GeoJSON）
- 含 Leiden cluster 標籤、UMAP 座標、gene expression

**Loupe Browser**：
- 輸出 barcode whitelist 與 cluster assignment

---

## 技術亮點

1. **MCseg v2 多模型 Voronoi 集成**：cyto3 三直徑 + 可選 Hematoxylin pass，以 merge_masks_fast IoU 去重後套用 Voronoi 擴張（無重疊保證），LUAD PQ 從 0.432 提升至 0.554（+28%）。

2. **Voronoi 擴張取代 expand_labels**：背景像素以 Voronoi tessellation 指派至最近細胞，距離上限 `voronoi_distance` px，根本消除邊界重疊問題。

3. **轉錄本密度補救（選配）**：利用 `vhd_pseudo_transcripts.csv` 高密度區域植入漏偵測細胞，改善 RNA 稀疏組織的細胞召回率。

4. **Tissue Profile 系統**：三層 YAML 合併，換組織只需改一行，TME panels、分割參數全自動切換。

5. **全向量化稀疏運算**：RNA 計數（A @ X）以 scipy 稀疏矩陣實現，萬級細胞 × 萬級基因 < 30 秒。

6. **完全非同步**：所有耗時後端任務以 FastAPI `BackgroundTasks` 執行，WebSocket 即時串流 log 至前端 xterm.js Terminal。

7. **OpenSeadragon DZI Tile Server**：Stage 0 直接對 BTF 分層串流，瀏覽器內流暢操作 Gigapixel 影像。

8. **macOS 外接硬碟防護**：自動過濾 `._*` metadata 檔案，確保在外接硬碟上穩健運行。

---

## 輸出檔案結構

```text
results/
  roi/{ROI_NAME}/
    he_crop.tif                   # Stage 0：H&E 裁切影像
    adata_002um.h5ad              # Stage 0：ROI 2µm bins
    segmentation_masks.npy        # Stage 1：MCseg v2 遮罩（H×W int32）
    segmentation_masks.tif        # Stage 1：同上，TIF 格式
    cellpose_cells.h5ad           # Stage 2：RNA 計數（cells × 18K genes）
    qc_preprocessed.h5ad          # Stage 3：QC 後矩陣
    umap_computed.h5ad            # Stage 3：含 UMAP + Leiden 標籤

figures/
  {ROI_NAME}/
    umap_leiden.png
    violin_qc.png
    dotplot_markers.png

results/export/{ROI_NAME}/
  cellpose_polygons.json          # Xenium Explorer 格式（GeoJSON）
  visiumhd_transcripts.csv        # Visium HD 轉錄點（x, y, gene）
  loupe_clusters.csv              # Loupe Browser 格式
```

---

## 測試與除錯

```bash
# 後端單元測試
uv run pytest backend/tests/ -v

# API 文件（啟動後前往）
# http://localhost:8001/docs
```

### 常見問題

| 問題 | 原因 | 解法 |
| --- | --- | --- |
| MCseg v2 速度慢 | CPU 模式 | 設定 `use_gpu: true`，確認 CUDA 可用 |
| 細胞數偏少 | cellprob_threshold 過高 | 降低至 -2.0 或 -3.0 |
| 細胞過小/碎片化 | min_size 過低 | 提高 `min_size`（如 50 px²） |
| bins 指派率低 | dilation_px=0 | 設定 `rna_counting.dilation_px: 6` |
| macOS `._*` 污染 | ExFAT 外接硬碟 | Pipeline 已內建自動過濾 |
