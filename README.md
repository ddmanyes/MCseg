# VisiumHD Pipeline 2

VisiumHD Pipeline 2 是一個為處理 10x Genomics **Visium HD** 與 **Xenium** 空間轉錄體學資料而設計的現代化全端應用程式。此專案整合了從原圖處理、細胞分割、H&E 光柵化、到 Proseg 分子分配及下游分析的一站式流程。

本專案將原本基於 PyQt6 桌面版的流程，升級為 **FastAPI 後端** 搭配 **React + Vite 前端**，不僅效能大增，更可在瀏覽器中輕鬆進行視覺化操作與即時日誌追蹤。

---

## 🚀 快速啟動

### 系統需求

- **uv**：極速的 Python 套件管理器 (`curl -LsSf https://astral.sh/uv/install.sh | sh`)
- **Node.js**：供前端使用 (建議 v18 以上)

### 一鍵執行

您可以使用我們提供的啟動腳本，它會自動啟動後端伺服器 (port 8000) 並啟動前端開發伺服器 (port 3000)：

```bash
cd visiumHD_pipeline_2
bash start.sh
```

啟動後，請開啟瀏覽器並前往：**<http://localhost:3000>**

### 分開啟動 (開發模式)

若你需要分別除錯，可以在兩個終端機中分開執行：

**終端機 1 (後端)：**

```bash
uv run uvicorn backend.main:app --reload --port 8000
```

**終端機 2 (前端)：**

```bash
cd frontend
npm install
npm run dev
```

---

## 🖥 網頁介面與操作說明

Pipeline 包含多個串聯的步驟，您可以透過左側選單在不同頁面間導航。

### 📂 資料設定 (Data Setup)

- **用途**：自動掃描並載入所需的原始數據。
- **操作說明**：點擊「掃描」按鈕並使用**本機資料夾瀏覽器** (Folder Browser) 選擇您的組織樣本根目錄 (例如 CRC 或 LUAD 資料夾)。後端會自動為您尋找 H&E BTF、`square_002um`、`square_008um` 以及 Xenium `outs` 資料夾，完全免去手動輸入路徑的麻煩。

### ✂️ Stage 0: ROI 裁切 (ROI Extract)

- **用途**：大型組織影像 (達 10GB 以上) 全局處理太佔記憶體。此模組允許您僅載入特定的感興趣區域 (ROI - Region of Interest)。
- **操作說明**：定義好 ROI 後，系統會精準裁切 H&E 影像、Visium HD `h5ad` 矩陣，並利用 PyArrow 從 Xenium `transcripts.parquet` 中快速切出所需基因座標，大幅節省後續 GPU 資源。

### 🦠 Stage 1: 細胞分割 (Segmentation)

- **用途**：利用高解析度 H&E 影像標定細胞核與細胞質範圍。
- **操作說明**：點選執行後，後端會自動調用 **Cellpose** (可選 `nuclei` 或是針對大腸癌的 `cyto2` 模型) 對影像進行 tile-based 分割，並套用 Eosin 染色分析與分水嶺演算法 (Watershed) 優化細胞邊界。

### 🧱 Stage 2: Zarr 建構 (Zarr Builder)

- **用途**：將零散的分析資料打包為具擴充性的高維空間資料格式 (Zarr)。
- **操作說明**：系統會自動抓取 Visium HD `.h5` 矩陣與前一步驟產生的 `segmentation_masks.npy`，打包成適用於 SpatialData API 及 Dask 讀取的 `proseg_integrated.zarr`。

### ⚙️ Stage 2.5: 條件測試 (Condition Test)

- **用途**：在執行完整 Proseg 之前，用小型 ROI 驗證哪些參數組合最適合您的組織種類。
- **操作說明**：設定 Proseg 參數網格 (例如：`max_dist` 或 `dilation` 不同的值)，後端將會並行處理這些測試條件，幫助您省去瞎猜參數的時間浪費。

### 🚀 Stage 3: Proseg 執行 (Run Proseg)

- **用途**：將轉錄本準確分配給分割出的單一細胞。
- **操作說明**：決定好黃金參數 (Golden Params) 後啟動此階段。此步驟對計算資源要求最高，後端會使用 tile-based 分塊處理技術 (含 overlap padding) 以避免記憶體溢位，並在前端提供即時 Terminal 串流日誌。

### 📊 Stage 4: 下游分析 (Analysis)

- **用途**：執行單細胞層級的品質控制 (QC)、降維與聚類 (Clustering)。
- **操作說明**：基於 Scanpy 引擎，自動濾除低基因表現細胞或高粒線體細胞。進行 Normalize、選取高變異基因 (HVG)、PCA 降維、UMAP 投影，最終使用 Leiden 演算法找出細胞聚類並繪製出分析圖 (存於 `figures/` 資料夾)。

### 📤 Stage 5: Browser 匯出 (Export)

- **用途**：將分析與聚類結果轉換為可視化軟體的相容格式，方便與醫學界分享。
- **操作說明**：支援一鍵匯出至 **10x Genomics Xenium Explorer** 與 **Loupe Browser**，讓研究人員可以直接在互動軟體中探索發現的新細胞類型與基因分佈特徵。

---

## 🛠 技術亮點

1. **完全非同步 (Fully Async)**：後端耗時任務採用 FastAPI `BackgroundTasks`，保證 UI 不卡頓，並透過 **WebSocket** 即時將命令列輸出串流推送到前端。
2. **PyArrow Predicate Pushdown**：在讀取 Xenium `transcripts.parquet` (通常 >10GB) 時，不讀入全表，而是將條件推演至底層只讀取被選取的 ROI，大幅減少 RAM 使用量。
3. **極度節省 VRAM 的 GPU 分塊**：不管是 Cellpose 影像分割或者是 Proseg 行為，全都是透過 tile/overlapping 技術讓一般顯卡也能跑得動。
4. **macOS 外接硬碟友好 (`._` 防護)**：遇到 ExFAT 產生的 `._` 資源分叉檔時，套件管理 (`uv` symlink) 以及 `discovery.py` 自動過濾，解決無法建置環境的痛點。

---

## 測試與除錯

本專案搭載了完整的 `pytest` 測試框架 (超過 50 項測試)：

```bash
uv run pytest backend/tests/ -v
```

測試涵蓋資料結構掃描、API 回應狀態與 Async HTTP 客戶端，甚至包括 Xenium `parquet` 與 `BTF` 高階影像標籤的 metadata 驗證。
