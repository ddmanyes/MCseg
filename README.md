# VisiumHD Pipeline 3

VisiumHD Pipeline 3 是一個為處理 10x Genomics **Visium HD** 空間轉錄體學資料而設計的精簡全端應用程式。相較於 Pipeline 2，本版本**移除 Zarr 建構與 Proseg**，改用 Cellpose 分割遮罩直接將 2µm bins 分配至細胞，大幅縮短分析流程並降低運算資源需求。

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

### 🦠 Stage 1: 細胞分割 (Segmentation)

- **用途**：利用高解析度 H&E 影像標定細胞核與細胞質範圍。
- **操作說明**：調用 **Cellpose**（可選 `nuclei` 或 `cyto2` 模型）對影像進行 tile-based 分割，並套用 Eosin 染色分析與分水嶺演算法優化邊界。
- **輸出**：`segmentation_masks.npy`（H×W 整數陣列，像素值 = cell ID）

### 🧬 Stage 2: RNA 計數 (RNA Count)

- **用途**：將 Visium HD 2µm bins 的 RNA 計數依 Cellpose 分割遮罩分配至細胞層級。
- **操作說明**：後端讀取 `adata_002um.h5ad` 的 bin 空間座標，對應至 `segmentation_masks.npy` 的像素位置，以稀疏矩陣乘法高效彙總每個細胞的基因計數。**不需要 Zarr 或 Proseg**。
- **輸出**：`cellpose_cells.h5ad`（cells × genes 稀疏矩陣，含 ROI local µm 座標）

### 📊 Stage 3: 下游分析 (Analysis)

- **用途**：執行單細胞層級的品質控制（QC）、降維與聚類（Clustering）。
- **操作說明**：基於 Scanpy 引擎，濾除低基因表現細胞與高粒線體比例細胞，進行 Normalize → HVG → PCA → UMAP → Leiden 聚類，並繪製分析圖（存於 `figures/`）。
- **輸入**：`cellpose_cells.h5ad`（Stage 2 輸出）

### 📤 Stage 4: Browser 匯出 (Export)

- **用途**：將分析與聚類結果轉換為可視化軟體相容格式。
- **操作說明**：使用 `skimage.measure.find_contours` 從 Cellpose 遮罩提取細胞輪廓多邊形，轉換為 GeoJSON 格式。支援一鍵匯出至 **10x Genomics Xenium Explorer** 與 **Loupe Browser**。

---

## 🛠 技術亮點

1. **零 Zarr/Proseg 架構**：直接稀疏矩陣計數（`scipy.sparse.csr_matrix` + `A @ adata.X`），比 Pipeline 2 少兩個耗時 Stage，適合快速迭代分析。
2. **完全非同步 (Fully Async)**：後端耗時任務採用 FastAPI `BackgroundTasks`，保證 UI 不卡頓，並透過 **WebSocket** 即時將命令列輸出串流推送到前端。
3. **xterm.js Terminal**：前端 Terminal 元件使用 xterm.js canvas 渲染，支援 ANSI 顏色（ERROR 紅、WARNING 黃、DEBUG 灰、INFO 綠）與增量寫入，不因大量 log 重新渲染整個 DOM。
4. **TanStack Query Polling**：所有 stage 狀態查詢改用 TanStack Query，`refetchInterval` 只在 `status === 'running'` 時啟動，unmount 自動清理，解決 setInterval 記憶體洩漏問題。
5. **OpenSeadragon DZI Tile Server**：Stage 0 直接對 BTF 進行 tile-by-tile 讀取，透過 Deep Zoom Image 協定將 gigapixel 組織影像分層串流至前端。
6. **Cellpose 原生平滑輪廓**：Cellpose 訓練於 flow field 梯度，輸出輪廓天然圓滑，使用 `skimage.measure.find_contours` 提取後無需額外多邊形平滑處理。
7. **macOS 外接硬碟友好 (`._` 防護)**：`discovery.py` 自動過濾 `._*` 與 `.DS_Store`，uv 套件管理使用 symlink 指向 SSD，確保 ExFAT 環境下穩定運作。

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
| Zarr 建構 | ✅ Stage 2 | ❌ 移除 |
| Proseg 條件測試 | ✅ Stage 2.5 | ❌ 移除 |
| Proseg 分子分配 | ✅ Stage 3 | ❌ 移除 |
| RNA 計數方式 | Proseg 機率分配 | **直接 mask 查詢 + 稀疏矩陣** |
| 細胞輪廓來源 | Proseg GeoJSON | **skimage.find_contours** |
| 後端 Port | 8000 | **8001** |
| 分析輸入 | `proseg_cells.h5ad` | `cellpose_cells.h5ad` |

## 🛠 測試與除錯

```bash
uv run pytest backend/tests/ -v
```
