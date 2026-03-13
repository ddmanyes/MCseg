# VisiumHD Pipeline 3

VisiumHD Pipeline 3 是一個為處理 10x Genomics **Visium HD** 空間轉錄體學資料而設計的架構優化版。相較於前代，本版本將核心路徑精簡為 **Cellpose 原生分配模式**，大幅縮短分析流程；同時**保留了 Proseg 作為選配 (Stage 2.5)**，提供分子層級的高精度重分配能力。

架構同樣採用 **FastAPI 後端** 搭配 **React + Vite 前端**，支援瀏覽器中的視覺化操作與即時日誌追蹤。

---

## 🚀 快速啟動

### 系統需求

- **uv**：極速的 Python 套件管理器 (`curl -LsSf https://astral.sh/uv/install.sh | sh`)
- **Node.js**：供前端使用 (建議 v18 以上)

### 一鍵執行

```bash
cd visiumHD_pipeline_3
bash start.sh
```

啟動後，請開啟瀏覽器並前往：**<http://localhost:3000>**

> **注意**：後端預設使用 **port 8001**（避免與 Pipeline 2 的 port 8000 衝突）。

### 分開啟動 (開發模式)

**終端機 1 (後端)：**

```bash
uv run uvicorn backend.main:app --reload --port 8001
```

**終端機 2 (前端)：**

```bash
cd frontend
npm install
npm run dev
```

---

## 🖥 網頁介面與操作說明

Pipeline 3 共 6 個步驟，透過左側選單（或頂部進度條）在頁面間導航。前一步驟完成前，後續步驟會鎖定。

### 📂 資料設定 (Data Setup)

- **用途**：自動掃描並載入所需的原始數據。
- **操作說明**：點擊「掃描」按鈕並選擇組織樣本根目錄，後端自動尋找 H&E BTF、`square_002um`、`square_008um` 資料夾，完全免去手動輸入路徑的麻煩。

### ✂️ Stage 0: ROI 裁切 (ROI Extract)

- **用途**：從 gigapixel 等級的大型組織影像中裁切感興趣區域（ROI），避免全圖載入耗盡記憶體。
- **操作說明**：頁面內嵌 **OpenSeadragon** 多解析度瀏覽器，直接從原始 BTF 動態載入 DZI tile，支援即時縮放與平移。拖曳框選即可自動填入 fullres pixel 座標，系統精準裁切 H&E 影像與 Visium HD `h5ad` 矩陣。
- **輸出**：`he_crop.tif`、`adata_002um.h5ad`（ROI 範圍內的 2µm bins）

- **用途**：利用高解析度 H&E 影像標定細胞核與細胞質範圍。
- **操作說明**：調用 **Cellpose**（ViT Transformer 架構）對影像進行分割。本版本新增了 **「無縫拼接技術 (Label Reconciliation)」** 與 **「原生 Tiling 融合」**，徹底解決了大圖分塊運算時產生的邊界接縫與綠色線條問題。
- **輸出**：`segmentation_masks.npy`（H×W 整數陣列，像素值 = cell ID）

### 🧬 Stage 2: RNA 計數 (RNA Count)

- **用途**：將 Visium HD 2µm bins 的 RNA 計數依 Cellpose 分割遮罩分配至細胞層級。
- **操作說明**：後端讀取 `adata_002um.h5ad` 的 bin 空間座標，對應至 `segmentation_masks.npy` 的像素位置，以稀疏矩陣法高效彙總。此為預設路徑，具有極高的運算速度。
- **輸出**：`cellpose_cells.h5ad`（或 `proseg_cells.h5ad`，若後續執行 Stage 2.5）

### 🧪 Stage 2.5: Proseg 重分配 (Proseg - Optional)

- **用途**：進階功能。使用概率模型將單個 RNA 分子重新指派給最可能的細胞。
- **操作說明**：若 Cellpose 分割結果在某些緻密區域不夠理想，可執行此步驟。系統會自動建立 Zarr 數據立方體並調用 Proseg 進行 MCMC 抽樣分配。
- **輸出**：`proseg_cells.h5ad`、`proseg_results.json`

### 📊 Stage 3: 下游分析 (Analysis)

- **用途**：執行單細胞層級的品質控制（QC）、降維與聚類（Clustering）。
- **操作說明**：基於 Scanpy 引擎，濾除低基因表現細胞與高粒線體比例細胞，進行 Normalize → HVG → PCA → UMAP → Leiden 聚類，並繪製分析圖（存於 `figures/`）。
- **輸入**：`cellpose_cells.h5ad`（Stage 2 輸出）

### 📤 Stage 4: Browser 匯出 (Export)

- **用途**：將分析與聚類結果轉換為可視化軟體相容格式。
- **操作說明**：使用 `skimage.measure.find_contours` 從 Cellpose 遮罩提取細胞輪廓多邊形，轉換為 GeoJSON 格式。支援一鍵匯出至 **10x Genomics Xenium Explorer** 與 **Loupe Browser**。

---

## 🛠 技術亮點

1. **混合多路徑架構**：提供「極速 Cellpose 直接計數」與「高精度 Proseg 分子分配」雙重路徑，適配不同研究需求。
2. **無縫拼接優化**：首創 `reconcile_stitched_labels` 演算法，在 GPU 分塊運算後自動偵測並修補斷裂的細胞核，徹底消除切割線。
3. **完全非同步 (Fully Async)**：後端耗時任務採用 FastAPI `BackgroundTasks` 與 WebSocket 即時回報 log，保證大型 ROI 運算時 UI 流暢。
4. **xterm.js Terminal**：前端使用高效能 canvas 渲染終端機日誌，支援即時彩色視圖與異常標記。
5. **OpenSeadragon DZI Tile Server**：Stage 0 直接對 BTF 進行分層串流，支援瀏覽器內流暢操作 Gigapixel 大型圖檔。
6. **macOS 外接硬碟友好 (`._` 防護)**：自動過濾 ExFAT / APFS 混用產生的 metadata 檔案，確保 Pipeline 穩健運行。

---

## 📁 輸出檔案結構

```
results/
  roi/{ROI_NAME}/
    he_crop.tif                 # Stage 0：H&E 裁切圖
    adata_002um.h5ad            # Stage 0：ROI 範圍內的 2µm bins
    segmentation_masks.npy      # Stage 1：Cellpose 遮罩（H×W int，0=背景）
    cellpose_cells.h5ad         # Stage 2：細胞級計數矩陣（cells × genes）
  analysis/{ROI_NAME}/
    adata_processed.h5ad        # Stage 3：Scanpy 分析結果
    figures/                    # UMAP、Violin plot 等
  export/{ROI_NAME}/
    cellpose_polygons.json      # Stage 4：GeoJSON 細胞輪廓（Xenium Explorer）
```

## 🔬 與 Pipeline 2 的差異

| 功能 | Pipeline 2 | Pipeline 3 |
|------|-----------|-----------|
| 核心路徑 | Proseg Only | **Cellpose 優先 (Proseg 選配)** |
| 處理速度 | 較慢 (需 Zarr/MCMC) | **極速 (直接計數) / 靈活 (可選 MCMC)** |
| 拼接問題 | 存在邊界接縫 | **已修復 (Seamless Stitching)** |
| RNA 計數方式 | Proseg 機率分配 | **直接 mask 查詢 或 Proseg 重分配** |
| 後端 Port | 8000 | **8001** |
| 環境管理 | Pip/Conda | **uv (推薦)** |

## 🛠 測試與除錯

```bash
uv run pytest backend/tests/ -v
```
