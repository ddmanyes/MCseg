# visiumHD_pipeline_2 — 進度摘要

> 最後更新：2026-03-03（Session UI 升級）

## 專案定位

基於 `visiumhd_pipeline`（PyQt6）重構的 **Web 版空間轉錄組分析平台**。
整合 `xenium_visiumhd_comparison`（ROI 裁切）與 `Proseg-Zarr-Integration`（Browser 匯出），並新增 Stage 2.5 Proseg 條件測試。

**後端**：FastAPI + uvicorn（port 8000）
**前端**：React 18 + Vite + Tailwind CSS（port 3000）

---

## 已完成功能

### 後端模組（34 個 Python 檔）

| Stage | 模組路徑 | 功能 | 狀態 |
|-------|---------|------|------|
| Setup | `backend/src/utils/discovery.py` | 資料目錄自動掃描 | ✅ |
| Setup | `backend/src/api/data.py` | /api/data 掃描與設定 | ✅ |
| 0 | `backend/src/roi/extractor.py` | BTF tile 讀取 + AnnData ROI 裁切 | ✅ |
| 0 | `backend/src/roi/tile_server.py` | **DZI tile server（OpenSeadragon 用）** | ✅ |
| 0 | `backend/src/api/roi.py` | ROI CRUD + 預覽 + 萃取 + **/dzi + /tiles 端點** | ✅ |
| 1 | `backend/src/segmentation/cellpose_runner.py` | Cellpose + Logic A + Eosin Watershed | ✅ |
| 1 | `backend/src/segmentation/macenko.py` | Macenko 色彩標準化 | ✅ |
| 1 | `backend/src/api/segmentation.py` | 分割 API（含前端參數覆寫）| ✅ |
| 2 | `backend/src/zarr_builder/builder.py` | Zarr 建構（含 macOS 汙染防護）| ✅ |
| 2 | `backend/src/api/zarr_builder.py` | Zarr API | ✅ |
| 2.5 | `backend/src/proseg/condition_tester.py` | 網格搜尋 + HE 疊圖縮圖生成 | ✅ |
| 2.5 | `backend/src/api/conditions.py` | 條件測試 API + /thumbnail/{idx} | ✅ |
| 3 | `backend/src/proseg/pipeline.py` | Proseg 完整 pipeline | ✅ |
| 3 | `backend/src/api/proseg.py` | Proseg API | ✅ |
| 4 | `backend/src/analysis/preprocessing.py` | QC + 正規化 + HVG | ✅ |
| 4 | `backend/src/analysis/clustering.py` | PCA + UMAP + Leiden | ✅ |
| 4 | `backend/src/analysis/pipeline.py` | 分析串流 + 圖表輸出 | ✅ |
| 5 | `backend/src/export/xenium_exporter.py` | Xenium Explorer 格式（含 pixel_size bug 修補）| ✅ |
| 5 | `backend/src/export/loupe_exporter.py` | Loupe Browser 格式（ijson 串流）| ✅ |
| 5 | `backend/src/api/export.py` | 匯出 API | ✅ |

### 前端頁面（8 個 React 頁面）

| 路由 | 頁面 | 功能 |
|------|------|------|
| `/data` | DataSetup | 資料路徑掃描與設定 |
| `/roi` | Stage0_ROI | ROI 定義（CRUD）+ **OpenSeadragon gigapixel 互動式選取** |
| `/segmentation` | Stage1_Segmentation | **完整參數面板**（模型/雙尺寸/前後處理/分塊 + 3 個快速預設）|
| `/zarr` | Stage2_Zarr | Zarr 建構進度 |
| `/conditions` | Stage2b_ConditionTest | **Top 3 疊圖縮圖 + 完整排序表格 + 散點圖**|
| `/proseg` | Stage3_Proseg | Proseg 執行（自動載入推薦參數）|
| `/analysis` | Stage4_Analysis | UMAP 視覺化 |
| `/export` | Stage5_Export | Xenium Explorer + Loupe Browser 雙格式匯出 |

### 本 Session UI 升級（Priority 1–4 全完成）

| 功能 | 改動 |
| --- | --- |
| TanStack Query | 所有 stage 移除 `setInterval`，改用 `useStageStatus` hook |
| PipelineStepper | 頂部水平 8 步進度條 + Sidebar 鎖定邏輯 |
| xterm.js Terminal | ANSI 顏色（ERROR/WARNING/DEBUG/INFO）+ 增量寫入 |
| OpenSeadragon ROI | BTF DZI tile server + OSD viewer + Pan/Draw 模式切換 |

### 測試套件（5 個模組）

- `test_01_infra.py` — config、constants、discovery
- `test_02_data.py` — BTF、binned 目錄、xenium outs
- `test_03_api.py` — FastAPI 端點健康檢查
- `test_04_roi.py` — BTF tiled TIFF、extractor
- `test_05_xenium.py` — transcripts/cells parquet schema

### 已產出的真實資料

```
results/analysis/roi/
├── CRC_tumor_boundary/
│   ├── adata_002um.h5ad   ← 2µm Visium HD ROI AnnData
│   ├── adata_008um.h5ad   ← 8µm 聚合 ROI AnnData
│   └── he_crop.tif        ← H&E 裁切影像
└── CRC_normal_colon/
    ├── adata_002um.h5ad
    ├── adata_008um.h5ad
    └── he_crop.tif
```

---

## 關鍵技術備忘

| 項目 | 說明 |
|------|------|
| 物理常數 | `XENIUM=0.2125`, `VISIUM=0.2737`, `PROSEG=0.2645833` µm/px |
| Proseg 黃金參數 | `dilation=20, max_dist=40, compactness=0.06` |
| Xenium pixel_size bug | `_patch_experiment_xenium()` 必須在每次 write 後執行 |
| BTF 讀取 | tile-by-tile，禁止全圖載入 |
| macOS 汙染防護 | `rglob("._*")` 在所有 zarr 操作前執行 |
| Dask monkey-patch | `dask.config.set({"dataframe.query-planning": True})` |
| Loupe 條碼限制 | 16bp ATCG + "-1"，類別數 < 32,768 |

---

## 啟動方式

```bash
cd /Volumes/SSD/plan_a/visiumHD_pipeline_2
bash start.sh
# 後端：http://localhost:8000  (API Docs: /docs)
# 前端：http://localhost:3000
```

---

## 下一步（選擇性）

- [ ] 執行 `uv run pytest backend/tests/` 確認測試全通過
- [ ] 在 Stage 0 測試真實 BTF tile 讀取流程
- [ ] 補充 Stage 2.5 縮圖：有真實 Proseg 輸出後驗證 GeoJSON 讀取路徑
- [ ] 加入 Stage 4 marker gene 視覺化
