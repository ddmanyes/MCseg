# visiumHD_pipeline_2 — 完整實作計畫

> 版本：v1.0
> 建立日期：2026-03-03
> 基於：`visiumhd_pipeline` + `xenium_visiumhd_comparison` + `Proseg-Zarr-Integration`

---

## 一、專案目標

在原 `visiumhd_pipeline` 基礎上，整合三項重大擴充：

| # | 功能 | 來源參考 |
|---|------|--------|
| 1 | **ROI 裁切前處理**（Stage 0） | `xenium_visiumhd_comparison` |
| 2 | **Proseg 參數條件測試**（Stage 2.5） | 全新設計 |
| 3 | **Browser 輸出匯出**（Stage 5） | `Proseg-Zarr-Integration/scripts/export_*` |
| 4 | **React 前端 UI**（替換 PyQt6） | 全新設計（FastAPI + React + Vite） |

---

## 二、系統架構

```
visiumHD_pipeline_2/
├── CLAUDE.md                    # 開發維護規範
├── plan.md                      # 本文件
├── README.md                    # 使用說明（後續建立）
├── pyproject.toml               # Python 依賴（uv 管理）
├── package.json                 # 前端 monorepo 根配置
│
├── config/
│   ├── pipeline.yaml            # 統一參數設定（含 ROI 定義）
│   └── roi_presets.yaml         # ROI 預設模板
│
├── backend/                     # Python FastAPI 後端
│   ├── main.py                  # FastAPI 應用程式入口
│   ├── src/
│   │   ├── api/                 # REST API 端點
│   │   │   ├── roi.py
│   │   │   ├── segmentation.py
│   │   │   ├── zarr_builder.py
│   │   │   ├── proseg.py
│   │   │   ├── analysis.py
│   │   │   └── export.py
│   │   ├── roi/                 # [新] Stage 0: ROI 裁切
│   │   │   ├── extractor.py     # H&E BTF 裁切、AnnData 裁切
│   │   │   ├── xenium_crop.py   # Xenium 轉錄本/邊界裁切
│   │   │   └── visualizer.py    # ROI 預覽產圖
│   │   ├── segmentation/        # Stage 1（沿用）
│   │   │   ├── cellpose_runner.py
│   │   │   └── macenko.py
│   │   ├── zarr_builder/        # Stage 2（沿用）
│   │   │   └── builder.py
│   │   ├── proseg/              # Stage 2.5 + Stage 3
│   │   │   ├── pipeline.py      # 完整 Proseg 執行（沿用）
│   │   │   ├── condition_tester.py  # [新] 參數網格搜索
│   │   │   ├── zarr_handler.py
│   │   │   └── metrics.py       # 評估指標計算
│   │   ├── analysis/            # Stage 4（沿用）
│   │   │   ├── preprocessing.py
│   │   │   ├── clustering.py
│   │   │   └── annotation.py
│   │   ├── export/              # [新] Stage 5: Browser 匯出
│   │   │   ├── xenium_exporter.py  # → Xenium Explorer 格式
│   │   │   └── loupe_exporter.py   # → Loupe Browser 格式
│   │   └── utils/
│   │       ├── config.py
│   │       ├── constants.py     # 所有物理常數集中管理
│   │       └── logging.py
│   └── scripts/                 # CLI 腳本（可獨立執行）
│       ├── 00_roi/
│       │   ├── define_roi.py
│       │   └── extract_roi.py
│       ├── 01_segmentation/
│       ├── 02_build_zarr/
│       ├── 02b_condition_test/  # [新] Proseg 條件測試
│       │   └── run_condition_grid.py
│       ├── 03_proseg/
│       ├── 04_analysis/
│       └── 05_export/
│           ├── export_xenium.py
│           └── export_loupe.py
│
└── frontend/                    # React 前端
    ├── package.json
    ├── vite.config.ts
    ├── tsconfig.json
    ├── index.html
    └── src/
        ├── main.tsx
        ├── App.tsx
        ├── api/                 # API 呼叫層
        │   └── client.ts
        ├── components/
        │   ├── layout/
        │   │   ├── Sidebar.tsx
        │   │   └── Header.tsx
        │   ├── roi/
        │   │   ├── HeImageViewer.tsx    # H&E 影像 + ROI 繪製
        │   │   └── RoiEditor.tsx
        │   ├── pipeline/
        │   │   ├── StageCard.tsx
        │   │   └── ProgressBar.tsx
        │   ├── proseg/
        │   │   └── ConditionComparison.tsx  # 參數比較視覺化
        │   └── shared/
        │       ├── Terminal.tsx         # 即時 log 輸出
        │       └── ConfigEditor.tsx
        └── pages/
            ├── Stage0_ROI.tsx
            ├── Stage1_Segmentation.tsx
            ├── Stage2_Zarr.tsx
            ├── Stage2b_ConditionTest.tsx
            ├── Stage3_Proseg.tsx
            ├── Stage4_Analysis.tsx
            └── Stage5_Export.tsx
```

---

## 三、各 Stage 詳細規格

### Stage 0：ROI 定義與裁切（全新）

**目標**：在正式分析前，讓使用者定義感興趣區域並裁切所有相關資料。

#### 0-A：ROI 定義方式

```yaml
# config/pipeline.yaml 中的 ROI 定義格式
rois:
  - name: "tumor_boundary"
    tissue: "CRC"
    # 方式 1：fullres pixel 座標（Visium HD 座標系）
    x: 43490
    y: 12515
    width_px: 6569
    height_px: 4791
    pixel_size_um: 0.2737
    # 方式 2：µm 座標（Xenium 原生座標系，當 Xenium 為獨立切片時）
    # x_um: 4356
    # y_um: 856
    # width_um: 1800
    # height_um: 1800
```

#### 0-B：裁切目標

| 資料類型 | 裁切方法 | 輸出 |
|---------|---------|------|
| H&E BTF/TIFF | tile-based 讀取（避免全圖載入） | `.npy` 或 `.tif` crop |
| Visium HD 8µm AnnData | `pxl_col/row_in_fullres` 篩選 | subset `.h5ad` |
| Visium HD 2µm AnnData | `pxl_col/row_in_fullres` 篩選 | subset `.h5ad` |
| Xenium transcripts | `x_location/y_location` 篩選（parquet pushdown） | subset parquet |
| Xenium nucleus_boundaries | `vertex_x/y` 篩選 + cv2 光柵化 | label mask |
| Xenium morphology TIFF | µm → px 轉換後裁切 | crop array |

#### 0-C：座標系統管理

```
Visium HD fullres px  ×0.2737  →  µm
Visium HD hires px    ×scalef  →  fullres px （scalef 讀自 scalefactors_json.json）
Xenium morphology px  ×0.2125  →  µm （LUAD & CRC 相同）
```

---

### Stage 1：細胞分割（沿用 visiumhd_pipeline）

- Cellpose + Macenko 色彩標準化
- Logic A 雙尺寸合併（nuclei dia + cyto dia）
- 分塊執行（block_size=2048, overlap=256）
- 輸出：細胞核遮罩 `.npy`

---

### Stage 2：Zarr 建構（沿用 visiumhd_pipeline）

- SpatialData OME-Zarr 格式
- 整合 2µm + 8µm + 細胞核遮罩
- 輸出：`proseg_integrated.zarr`

---

### Stage 2.5：Proseg 參數條件測試（全新）

**目標**：在正式 Proseg 跑全圖前，用小型 ROI 快速評估不同參數組合，選出最佳條件。

#### 參數網格設計

```python
CONDITION_GRID = {
    "max_dist":     [20, 30, 40, 50],   # µm（RNA 擴散最大距離）
    "compactness":  [0.03, 0.06, 0.1],  # 細胞緊密度
    "dilation":     [10, 20, 30],        # px（核遮罩擴張）
}
# 預設：3×3×3 = 27 種組合，但建議先做 2×2×2 = 8 種快速測試
```

#### 評估指標

| 指標 | 計算方法 | 意義 |
|------|---------|------|
| `n_cells` | `len(adata)` | 細胞偵測數量 |
| `median_genes` | `np.median(adata.obs.n_genes)` | 每細胞基因豐富度 |
| `median_counts` | `np.median(adata.obs.total_counts)` | 每細胞 UMI 數 |
| `fraction_assigned` | `assigned_tx / total_tx` | RNA 指派率 |
| `cell_area_cv` | `std(area) / mean(area)` | 細胞大小均一性 |
| `silhouette_score` | sklearn | 聚類分離度（可選） |

#### 輸出格式

```
results/02b_conditions/
├── condition_grid.csv           # 所有條件的評估指標表
├── condition_heatmap.png        # 指標熱力圖
├── recommended_params.yaml      # 建議最佳參數
└── per_condition/
    ├── cond_01_d20_c006_m30/
    │   ├── adata.h5ad
    │   └── metrics.json
    └── ...
```

#### UI 視覺化

- 互動式散點圖（`n_cells` vs `median_genes`）
- 參數熱力圖（2D grids）
- 自動高亮「帕累托最優」條件

---

### Stage 3：Proseg 完整執行（沿用 + 優化）

- 使用 Stage 2.5 選出的最佳參數
- 分塊執行、合併、背景過濾
- 輸出：`processed_proseg_cyto.h5ad` + GeoJSON 多邊形

#### 核心優化：Nuclear Shield (核保護器) 演算法

**背景問題**：Proseg 本質為純機率模型 (Gaussian Mixture Model)，給定 RNA 座標與初始分配種子後，它有可能因為局部 RNA 濃度分佈，導致預測的細胞多邊形（Polygon）**越界吞噬**相鄰的細胞核，或延伸至無組織的背景區域。

**解決方案**：結合 **硬約束 (Hard Constraint)**、**邊界調停 (Watershed)** 與 **柔性引導 (Soft Prior)** 的三重防護體系：

1. **防線一（Watershed 邊界調停與細胞質約束）**：
   - *Watershed 分配*：在處理重疊或緊密相連的擴張核遮罩時，透過距離轉換（Distance Transform）建立細胞間的絕對等距中線（楚河漢界），確保初始餵給 Proseg 的每個 RNA 點都有公平且不重疊的歸屬。
   - *細胞質過濾 (Cytoplasm Constraint)*：利用 Stage 1 獨有的 **Eosin BG Threshold** 機制。透過 Macenko 演算法將 H&E 影像中的 Eosin (粉紅色細胞質) 訊號分離轉為灰階後，設定閥值將影像「二值化」，產生一張非黑即白的「細胞質實體遮罩」。**在送入 Proseg 之前，將超出該實體遮罩（即背景空腔與無組織區）的擴張種子與游離轉錄點強制剃除（設定為 Background / Unassigned）**，從物理底層防堵模型演算法向空腔區域無限延伸多邊形的可能。
2. **防線二（Python 核心硬約束）**：在資料前處理階段分配初始 `cell_id` 時，無論分水嶺 (Watershed/Dilation) 如何擴張，只要 RNA 點絕對落在「原始未擴張的細胞核遮罩」內，代碼將強制覆寫並鎖死其原始 ID，不容任何擴張或分水嶺計算的誤差污染。
3. **防線三（Proseg 模型封印）**：執行 Proseg 時強制加入 CLI 參數 `--nuclear-reassignment-prob 0.01`（系統預設為 0.2），將演算法在內部迭代中「擅自重新分配核內 RNA 歸屬」的機率降至最低的 1%，徹底封印 Proseg 越界搶奪相鄰核的能力。

**最終效果**：完美確保每一個細胞核範圍皆為「神聖不可侵犯」，而在細胞核以外的細胞質區間，又能透過 Eosin 遮罩防堵空腔蔓延，並讓邊界根據真實的 RNA 濃度與 Watershed 的公平劃分自由探索（透過調整 `max_dist` 與 `compactness`），達到既不重疊又高度吻合細胞形態的精準分割。

---

### Stage 4：下游分析（沿用 visiumhd_pipeline）

- QC → normalize → HVG → PCA → UMAP → Leiden
- CellTypist 細胞型標注
- 輸出：`clustered_final.h5ad` + 圖表

---

### Stage 5：Browser 格式匯出（全新，來自 Proseg-Zarr-Integration）

#### 5-A：Xenium Explorer 格式

**所需輸出檔案**：

| 檔案 | 格式 | 說明 |
|-----|------|------|
| `experiment.xenium` | JSON | 實驗元數據（含 pixel_size = 0.2645833） |
| `morphology.ome.tif` | OME-TIFF 金字塔 | H&E 多解析度影像 |
| `transcripts.zarr.zip` | Zarr（壓縮） | 轉錄點位 + 基因名稱 |
| `cell_boundaries/` | Zarr | 細胞多邊形向量 |

**關鍵處理步驟**：

1. 多邊形座標：µm → px（÷ 0.2645833）
2. Shapely 簡化 + 平滑（去除 Watershed 毛邊）
3. 修補 `experiment.xenium` pixel_size Bug（`spatialdata_xenium_explorer` 硬編碼 0.2125）
4. 轉錄點 ID 重映射（Proseg global ID → AnnData local index）

**Bug 修復**（必須保留）：

```python
# spatialdata_xenium_explorer.write() 後強制修補
exp_data["pixel_size"] = PROSEG_SCALE_UM_PX  # 0.2645833
```

#### 5-B：Loupe Browser 格式

**所需輸出**：

| 檔案 | 說明 |
|-----|------|
| `.cloupe` | HDF5，含表達矩陣 + 元數據 + 聚類 |
| `.geojson` | 細胞邊界（空間疊圖用） |

**關鍵限制**：

- 細胞條碼必須符合 10X Genomics 白名單（16 bp）
- 分類欄位唯一值 ≤ 32,768
- 使用 `loupepy.create_loupe_from_anndata()`
- GeoJSON 使用 `ijson` 串流讀取（避免大型 JSON 記憶體爆炸）

---

## 四、React 前端架構

### 技術棧

| 層級 | 技術 | 選擇原因 |
|------|------|--------|
| 後端 | Python FastAPI | 與現有 scanpy 生態無縫整合 |
| API 通訊 | WebSocket + REST | WebSocket 用於即時 log 串流，REST 用於操作控制 |
| 前端框架 | React 18 + TypeScript | 現代化、強型別 |
| 建構工具 | Vite | 快速 HMR，適合開發 |
| UI 元件庫 | shadcn/ui + Tailwind CSS | 無頭元件，高度可客製 |
| 狀態管理 | Zustand | 輕量，適合流水線狀態 |
| 影像檢視 | OpenSeadragon | 支援 DZI 金字塔大圖（H&E 檢視） |
| 圖表 | Recharts / D3 | 條件比較視覺化 |
| 路由 | React Router v6 | SPA 導向 |

### UI 頁面規劃

```
┌─────────────────────────────────────────────────────────────────┐
│  VisiumHD Pipeline 2                                [⚙] [?]     │
├─────────┬───────────────────────────────────────────────────────┤
│         │                                                        │
│ Sidebar │  主要內容區域                                          │
│         │                                                        │
│ [0] ROI │  ┌─────────────────────────────────────────────────┐  │
│ [1] SEG │  │  Stage 0: ROI 定義與裁切                         │  │
│ [2] ZRR │  │                                                  │  │
│ [2.5]   │  │  ┌─────────────────┐  ┌─────────────────────┐  │  │
│  COND   │  │  │  H&E 影像預覽   │  │  ROI 設定面板        │  │  │
│ [3] PRO │  │  │  (OpenSeaDragon│  │  ─────────────────── │  │  │
│ [4] ANA │  │  │   + ROI 覆蓋)  │  │  名稱: [tumor_bd]    │  │  │
│ [5] EXP │  │  │               │  │  座標類型: [fullres] │  │  │
│         │  │  │  [拖拉繪製     │  │  X: [43490]          │  │  │
│         │  │  │   ROI 方框]   │  │  Y: [12515]          │  │  │
│         │  │  │               │  │  W: [6569]           │  │  │
│         │  │  │               │  │  H: [4791]           │  │  │
│         │  │  └─────────────────┘  │  [+ 新增 ROI]        │  │  │
│         │  │                       └─────────────────────┘  │  │
│         │  │  [執行裁切] ████████░░░░ 68%                    │  │
│         │  └─────────────────────────────────────────────────┘  │
│         │                                                        │
│         │  ┌─────── Stage Log ─────────────────────────────┐   │
│         │  │ [11:23:01] BTF tile-based reading...           │   │
│         │  │ [11:23:04] Subset AnnData: 36724 bins          │   │
│         │  └────────────────────────────────────────────────┘   │
└─────────┴───────────────────────────────────────────────────────┘
```

### 頁面 2.5：Proseg 條件測試（核心新功能）

```
Stage 2.5: Proseg 條件測試
─────────────────────────────────────────────────────────────────

參數網格設定:
┌──────────────────────────────────────────────────────────────┐
│  max_dist (µm):    [20] [30] [40] [50]  ☑全選               │
│  compactness:      [0.03] [0.06] [0.1]  ☑全選               │
│  dilation (px):    [10] [20] [30]        ☑全選               │
│                                                               │
│  測試 ROI 大小:  500×500 µm（快速）  [自訂...]              │
│  估算執行時間:   ~15 分鐘（8 種條件）                        │
│                                                               │
│  [開始測試]  [預覽條件清單]                                  │
└──────────────────────────────────────────────────────────────┘

結果比較:
┌─────────────────────────────┐  ┌─────────────────────────────┐
│  散點圖                      │  │  熱力圖                      │
│  X: n_cells  Y: median_genes │  │  max_dist vs compactness    │
│  ●最佳  ○其他               │  │  色階: n_cells               │
└─────────────────────────────┘  └─────────────────────────────┘

推薦條件: max_dist=40, compactness=0.06, dilation=20
[套用至 Stage 3] [詳細比較報告]
```

---

## 五、FastAPI 後端 API 設計

### REST API 端點

```
GET  /api/health
GET  /api/config                    # 讀取 pipeline.yaml
PUT  /api/config                    # 更新 pipeline.yaml

# Stage 0: ROI
POST /api/roi/preview               # 產生 H&E thumbnail + ROI 標注
POST /api/roi/extract               # 執行裁切
GET  /api/roi/status

# Stage 1: 分割
POST /api/segmentation/run
GET  /api/segmentation/status
GET  /api/segmentation/preview      # 回傳分割遮罩預覽圖

# Stage 2: Zarr
POST /api/zarr/build
GET  /api/zarr/status

# Stage 2.5: 條件測試
POST /api/conditions/run            # body: { conditions: [...] }
GET  /api/conditions/status
GET  /api/conditions/results        # 回傳所有條件的指標 JSON
GET  /api/conditions/recommend      # 回傳建議最佳條件

# Stage 3: Proseg
POST /api/proseg/run
GET  /api/proseg/status
GET  /api/proseg/preview

# Stage 4: 分析
POST /api/analysis/run
GET  /api/analysis/status
GET  /api/analysis/umap             # 回傳 UMAP 圖 base64

# Stage 5: 匯出
POST /api/export/xenium
POST /api/export/loupe
GET  /api/export/status/{job_id}
```

### WebSocket 即時 Log

```
WS   /ws/log/{stage}               # 即時串流各 stage 的 stdout/stderr
```

---

## 六、關鍵技術決策

### 決策 1：前端與後端分離

- 後端：Python FastAPI（port 8000）
- 前端：Vite dev server（port 3000）
- 生產部署：Vite build → FastAPI 靜態服務
- 好處：前端可獨立開發，後端計算不影響 UI 回應

### 決策 2：長任務管理

- 使用 Python `asyncio` + `subprocess`（非同步執行 Proseg 二進位）
- 任務狀態：`IDLE → RUNNING → DONE / ERROR`
- 前端透過 WebSocket 接收即時 log
- 支援任務取消（SIGTERM）

### 決策 3：ROI 影像預覽

- H&E BTF/TIFF 太大無法直接傳輸
- 後端產生 hires thumbnail（~2000px）+ ROI 矩形疊加
- 傳輸 JPEG（壓縮），前端 OpenSeadragon 顯示
- 選擇 ROI 後，後端執行實際的 tile-based 裁切

### 決策 4：Proseg 條件測試並行化

- 每個條件獨立執行（各用不同臨時目錄）
- `asyncio.gather()` 或 `concurrent.futures.ProcessPoolExecutor`
- 最多同時跑 4 個條件（避免 VRAM 耗盡）
- 進度透過 WebSocket 回報每個條件完成情況

### 決策 5：常數集中管理

```python
# backend/src/utils/constants.py
XENIUM_UM_PX      = 0.2125   # Xenium morphology 像素尺寸
PROSEG_UM_PX      = 0.2645833  # Proseg 輸出座標系
VISIUM_UM_PX      = 0.2737   # Visium HD fullres 像素尺寸
XENIUM_NM_PX      = 212.5    # Xenium nm/px

GOLDEN_PARAMS = {
    "dilation":   20,
    "max_dist":   40.0,
    "compactness": 0.06,
    "samples":    500,
    "recorded":   150,
    "watershed":  True,
}
```

---

## 七、實作優先順序

### Phase 1（核心基礎）

- [ ] 建立專案骨架（pyproject.toml、package.json、目錄結構）
- [ ] FastAPI 後端骨架 + WebSocket log 系統
- [ ] 從 `visiumhd_pipeline` 移植 Stage 1-4 核心邏輯
- [ ] React 基礎框架（路由、Sidebar、Stage 頁面骨架）

### Phase 2（新功能 1：ROI 裁切）

- [ ] 從 `xenium_visiumhd_comparison` 移植 `subset_to_roi()`、`_read_btf_crop()`、`rasterize_nucleus_mask()`、`load_transcripts_roi()`
- [ ] `backend/src/roi/extractor.py` — 統一 ROI 裁切介面
- [ ] API endpoint：`POST /api/roi/extract`
- [ ] React 頁面 `Stage0_ROI.tsx`：H&E thumbnail 顯示 + ROI 矩形繪製

### Phase 3（新功能 2：條件測試）

- [ ] `backend/src/proseg/condition_tester.py` — 參數網格執行引擎
- [ ] 評估指標計算（n_cells、median_genes、fraction_assigned、cell_area_cv）
- [ ] API endpoint：`POST /api/conditions/run`、`GET /api/conditions/results`
- [ ] React 頁面 `Stage2b_ConditionTest.tsx`：散點圖 + 熱力圖

### Phase 4（新功能 3：Browser 匯出）

- [ ] 從 `Proseg-Zarr-Integration` 移植 `export_to_xenium_full.py` → `backend/src/export/xenium_exporter.py`
- [ ] 移植 `export_to_loupe_merged.py` → `backend/src/export/loupe_exporter.py`
- [ ] 保留並整合所有已知 Bug 修復（pixel_size patch、Loupe 白名單條碼等）
- [ ] API endpoint：`POST /api/export/xenium`、`POST /api/export/loupe`
- [ ] React 頁面 `Stage5_Export.tsx`

### Phase 5（UI 優化）

- [ ] OpenSeadragon 整合（H&E 大圖縮放）
- [ ] 條件測試的互動式圖表（Recharts）
- [ ] 整體 UI 精修（深色主題、動畫、錯誤提示）
- [ ] 一鍵完整流程執行（Run All Stages）

---

## 八、依賴清單

### Python（後端）

```toml
[tool.uv.dependencies]
# Web 框架
fastapi = ">=0.111.0"
uvicorn = {extras = ["standard"], version = ">=0.30.0"}
websockets = ">=12.0"
python-multipart = ">=0.0.9"

# 空間轉錄組核心（沿用）
scanpy = ">=1.10.0"
anndata = ">=0.10.0"
squidpy = ">=1.3.0"

# 影像與分割（沿用）
cellpose = "*"
tifffile = ">=2025.5.10"
opencv-python-headless = ">=4.13.0"
torch = {version = "*", source = "pytorch-cu124"}

# Zarr 與空間格式（沿用）
zarr = ">=2.18.3,<3.0.0"
spatialdata = ">=0.2.6"
geopandas = "*"
shapely = "*"

# Browser 匯出（新增）
spatialdata-xenium-explorer = "*"
loupepy = "*"
ijson = "*"           # 串流讀取大型 JSON

# 資料處理
pandas = ">=2.0.0"
numpy = ">=1.24.0"
scipy = "*"
pyarrow = ">=14.0.0"

# 工具
pyyaml = "*"
tqdm = "*"
```

### Node.js（前端）

```json
{
  "dependencies": {
    "react": "^18.3.0",
    "react-dom": "^18.3.0",
    "react-router-dom": "^6.23.0",
    "zustand": "^4.5.0",
    "axios": "^1.7.0",
    "openseadragon": "^4.1.0",
    "@openseadragon/types": "*",
    "recharts": "^2.12.0",
    "tailwindcss": "^3.4.0",
    "@radix-ui/react-*": "*",
    "lucide-react": "^0.379.0"
  },
  "devDependencies": {
    "vite": "^5.2.0",
    "@vitejs/plugin-react": "^4.3.0",
    "typescript": "^5.4.0"
  }
}
```

---

## 九、已知風險與對策

| 風險 | 嚴重度 | 對策 |
|------|--------|------|
| `spatialdata_xenium_explorer` pixel_size Bug | 高 | 強制覆寫 `experiment.xenium`（已有成熟修復） |
| Loupe 需要 10X 白名單條碼 | 中 | 從 `Proseg-Zarr-Integration` 移植條碼生成邏輯 |
| Proseg 條件測試 VRAM 不足 | 中 | 限制最大並行數 4，並支援 CPU-only fallback |
| BTF 大檔案讀取速度 | 低 | tile-based 讀取（已驗證可行） |
| React 與 FastAPI WebSocket 相容性 | 低 | 使用標準 WebSocket，已有成熟方案 |
| macOS `._*` 檔案污染 Proseg | 中 | 沿用現有 `rglob("._*")` 清除邏輯 |

---

## 十、參考來源

| 功能 | 來源檔案 |
|------|---------|
| ROI fullres px 裁切 | `xenium_visiumhd_comparison/scripts/02_baseline/cluster_8um.py` |
| BTF tile-based 讀取 | `xenium_visiumhd_comparison/scripts/03_pipeline/segment_cellpose.py` |
| Xenium 轉錄本裁切 | `xenium_visiumhd_comparison/scripts/03_pipeline/segment_proseg_zarr.py` |
| nucleus_boundaries 光柵化 | `xenium_visiumhd_comparison/scripts/03_pipeline/segment_proseg_zarr.py` |
| Xenium Explorer 匯出 | `Proseg-Zarr-Integration/scripts/export_to_xenium_full.py` |
| Loupe Browser 匯出 | `Proseg-Zarr-Integration/scripts/export_to_loupe_merged.py` |
| Cellpose + Logic A | `visiumhd_pipeline/src/segmentation/cellpose_runner.py` |
| Proseg Pipeline | `visiumhd_pipeline/src/proseg/pipeline.py` |
| Scanpy 分析流程 | `visiumhd_pipeline/src/analysis/` |
