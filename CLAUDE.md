# visiumHD_pipeline_3 — 開發維護規範

## 1. 專案定位

**VisiumHD Pipeline 3** 是 Pipeline 2 的精簡版本，**移除 Zarr 建構與 Proseg**，改用 Cellpose 分割遮罩直接分配 RNA bins，大幅縮短分析流程。

**Pipeline 架構**：資料設定 → ROI 裁切 → Cellpose 分割 → RNA 計數 → Scanpy 分析 → Browser 匯出

| Stage | 功能 | 關鍵點 |
|-------|------|--------|
| Setup | 資料自動掃描 | /api/data/scan |
| 0 | ROI 裁切（BTF tile-by-tile）| extractor.py |
| 1 | Cellpose 分割 | segmentation_masks.npy（H×W int） |
| 2 | RNA 計數（新增）| cellpose_counter/counter.py，稀疏矩陣 |
| 3 | Scanpy 分析 | QC + UMAP + Leiden → cellpose_cells.h5ad |
| 4 | Browser 匯出 | skimage.find_contours → GeoJSON |

**相較 Pipeline 2 移除**：Zarr 建構（Stage 2）、Proseg 條件測試（Stage 2.5）、Proseg 完整執行（Stage 3）

---

## 2. 環境管理

### Python 後端

```bash
# 使用 uv 管理（強制，禁止使用 pip）
uv sync                           # 安裝/同步依賴
uv add <package>                  # 新增套件
uv run python <script>            # 執行腳本
uv run uvicorn backend.main:app --reload --port 8000  # 啟動後端
```

> **⚠️ ExFAT 注意**：`.venv` 是指向 SSD 上 venv 的 symlink（`/Volumes/SSD/plan_a/visiumHD_pipeline_2_venv`）。
> 若 symlink 失效（重開機後 `/tmp` 被清空是舊問題，現已修正指向 SSD），請重建：
>
> ```bash
> rm .venv && ln -s /Volumes/SSD/plan_a/visiumHD_pipeline_2_venv .venv
> # 清理 ExFAT 的 ._* 垃圾（避免 uv sync 失敗）
> find /Volumes/SSD/plan_a/visiumHD_pipeline_2_venv -name '._*' -delete
> UV_LINK_MODE=copy uv sync
> ```

### Node.js 前端

```bash
cd frontend
npm install                       # 安裝依賴（或 pnpm install）
npm run dev                       # 開發模式（port 3000）
npm run build                     # 建構生產版本
```

### 一鍵啟動（開發模式）

```bash
# Terminal 1: 後端（使用 port 8001 避免與 pipeline_2 衝突）
uv run uvicorn backend.main:app --reload --port 8001

# Terminal 2: 前端
cd frontend && npm run dev
```

---

## 3. 目錄結構規範

```
visiumHD_pipeline_3/
├── CLAUDE.md           # 本文件
├── config/
│   ├── pipeline.yaml   # 所有參數來源（禁止硬編碼）
│   └── roi_presets.yaml
├── backend/
│   ├── main.py         # FastAPI 入口（port 8001）
│   └── src/
│       ├── api/
│       │   ├── data.py            # /api/data
│       │   ├── roi.py             # /api/roi
│       │   ├── segmentation.py    # /api/segmentation
│       │   ├── cellpose_count.py  # /api/count（新增）
│       │   ├── analysis.py        # /api/analysis
│       │   └── export.py          # /api/export
│       ├── cellpose_counter/      # RNA 計數模組（新增）
│       │   └── counter.py
│       ├── roi/                   # Stage 0: ROI 裁切
│       ├── segmentation/          # Stage 1: Cellpose
│       ├── analysis/              # Stage 3: Scanpy
│       ├── export/                # Stage 4: Xenium/Loupe
│       └── utils/
│           └── constants.py       # 物理常數集中於此
└── frontend/src/
    ├── api/client.ts              # API 客戶端（無 zarr/proseg）
    ├── stores/pipelineStore.ts    # Zustand（含 count stage）
    └── pages/
        ├── DataSetup.tsx
        ├── Stage0_ROI.tsx
        ├── Stage1_Segmentation.tsx
        ├── Stage2_Count.tsx       # 新增
        ├── Stage3_Analysis.tsx
        └── Stage4_Export.tsx
```

---

## 4. 物理常數（集中管理）

所有座標換算常數必須從 `backend/src/utils/constants.py` 引用，**嚴禁在腳本內硬編碼**：

```python
# backend/src/utils/constants.py
XENIUM_UM_PX       = 0.2125     # Xenium morphology 像素尺寸（µm/px）
XENIUM_NM_PX       = 212.5      # 同，nm/px
PROSEG_UM_PX       = 0.2645833  # Proseg 輸出座標系（µm/px）
VISIUM_UM_PX       = 0.2737     # Visium HD fullres（µm/px，LUAD & CRC 相同）

# Proseg 黃金參數（可被 config 覆寫）
GOLDEN_PARAMS = {
    "dilation":     20,
    "max_dist":     40.0,
    "compactness":  0.06,
    "samples":      500,
    "recorded":     150,
    "watershed":    True,
    "connectivity": True,
}
```

---

## 5. ROI 裁切規範（Stage 0）

### 座標格式

```yaml
# config/pipeline.yaml 中的 ROI 定義
rois:
  - name: "tumor_boundary"
    tissue: "CRC"
    # 選擇一種格式：
    # 格式 A：Visium HD fullres pixel（CRC/LUAD 標準）
    x: 43490
    y: 12515
    width_px: 6569
    height_px: 4791
    pixel_size_um: 0.2737
    # 格式 B：Xenium 原生 µm（CRC Xenium 獨立切片）
    # x_um: 4356
    # y_um: 856
    # width_um: 1800
    # height_um: 1800
```

### 座標轉換公式

| 轉換 | 公式 |
|------|------|
| Visium fullres px → µm | `× 0.2737` |
| Visium fullres px → hires px | `× scalef`（從 `scalefactors_json.json` 讀取）|
| Xenium µm → morphology px | `÷ 0.2125` |
| Proseg µm → px | `÷ 0.2645833` |

### BTF tile-based 讀取

大型 BTF/TIFF 必須使用 tile-based 讀取（來自 `xenium_visiumhd_comparison`），**嚴禁全圖載入**：

```python
# 正確做法
def read_btf_crop(btf_path, x0, y0, w, h):
    with tifffile.TiffFile(btf_path) as tf:
        page = tf.pages[0]
        TW, TH = page.tilewidth, page.tilelength
        # ... tile-by-tile 讀取
```

---

## 6. Stage 2 RNA 計數核心演算法

Pipeline 3 直接用 Cellpose 遮罩分配 bins，**無需 Zarr 或 Proseg**。

```python
# backend/src/cellpose_counter/counter.py 核心流程

import numpy as np
import scanpy as sc
from scipy.sparse import csr_matrix
from scipy.ndimage import center_of_mass

# 1. 載入資料
mask = np.load("segmentation_masks.npy")        # H×W, pixel value = cell ID (0=bg)
adata = sc.read_h5ad("adata_002um.h5ad")

# 2. bin 座標轉換（Visium fullres px → ROI local pixel）
coords = adata.obsm['spatial']                   # (n_bins, 2) = [x_fullres, y_fullres]
col = (coords[:, 0] - roi_x_px).astype(int)
row = (coords[:, 1] - roi_y_px).astype(int)

# 3. 邊界過濾
valid_mask = (row >= 0) & (row < mask.shape[0]) & (col >= 0) & (col < mask.shape[1])

# 4. 查詢 cell ID（0 = 背景，>0 = 細胞）
cell_ids = mask[row[valid_mask], col[valid_mask]]
valid_cell = cell_ids > 0

# 5. 稀疏矩陣聚合（cell × bin），然後乘以計數矩陣
n_cells = mask.max()
A = csr_matrix((np.ones(valid_cell.sum()),
                (cell_ids[valid_cell] - 1,
                 np.where(valid_mask)[0][valid_cell])),
               shape=(n_cells, len(adata)))
count_matrix = A @ adata.X                       # (n_cells, n_genes)

# 6. 計算細胞質心
labels = np.unique(mask[mask > 0])
centroids = center_of_mass(np.ones_like(mask), mask, labels)  # [(row, col), ...]
```

**輸出 `cellpose_cells.h5ad`**：
- `obs` 欄位：`cell_id, n_bins, centroid_x_um, centroid_y_um, cell_area_um2`
- `obsm['spatial']`：ROI local µm 座標 [x_um, y_um]（原點 = ROI 左上角）

---

## 7. Browser 匯出規範（Stage 4）

### Xenium Explorer 匯出

**Pipeline 3 使用 skimage 從遮罩生成 GeoJSON**（無需 Proseg GeoJSON）：

```python
from skimage.measure import find_contours
import numpy as np

def _mask_to_geojson(mask_path, pixel_size_um, min_area_px=20):
    mask = np.load(mask_path)
    # 邊緣 padding 避免邊界細胞遺漏
    padded = np.pad(mask, 1, mode='constant')
    features = []
    for cell_id in np.unique(padded)[1:]:           # 跳過 0（背景）
        cell_mask = (padded == cell_id).astype(np.uint8)
        contours = find_contours(cell_mask, 0.5)
        if not contours:
            continue
        contour = max(contours, key=len) - 1        # 去除 padding 偏移
        if len(contour) < 4:
            continue
        # (row, col) → (x_um, y_um)
        coords = [(c[1] * pixel_size_um, c[0] * pixel_size_um) for c in contour]
        if len(coords) > 200:                        # 縮減頂點避免渲染過慢
            step = len(coords) // 200
            coords = coords[::step]
        coords.append(coords[0])                     # GeoJSON 閉合
        features.append({"type": "Feature", "geometry": {"type": "Polygon", "coordinates": [coords]}, ...})
    return {"type": "FeatureCollection", "features": features}
```

**注意**：Cellpose 輪廓天然平滑（flow field），**無需額外 shapely 平滑處理**。

### Loupe Browser 匯出

**必須使用 10X Genomics 白名單條碼**：

```python
# 10x_whitelist.txt 存於 backend/src/export/
with open("10x_whitelist.txt") as f:
    barcodes = [line.strip() for line in f]
adata.obs_names = [f"{barcodes[i]}-1" for i in range(len(adata))]
```

**分類欄位限制**：唯一值 > 32,768 須轉為整數：

```python
for col in adata.obs.columns:
    if adata.obs[col].nunique() > 32_000:
        adata.obs[col] = adata.obs[col].astype("category").cat.codes
```

**大型 GeoJSON 串流讀取**（禁止一次性 `json.load` 大型 JSON）：

```python
import ijson
with open(geojson_path, "rb") as f:
    for feat in ijson.items(f, "features.item"):
        # 逐個處理
```

---

## 8. FastAPI 後端規範

### 長任務模式

所有耗時操作（Proseg、Cellpose、匯出）必須使用非同步背景任務：

```python
from fastapi import BackgroundTasks
import asyncio

@app.post("/api/proseg/run")
async def run_proseg(background_tasks: BackgroundTasks):
    background_tasks.add_task(proseg_pipeline.run_async)
    return {"status": "started"}
```

### WebSocket Log 串流

```python
@app.websocket("/ws/log/{stage}")
async def websocket_log(websocket: WebSocket, stage: str):
    await manager.connect(websocket, stage)
    # 串流 subprocess stdout/stderr
```

### API 回應格式統一

```python
# 成功
{"status": "ok", "data": {...}}
# 錯誤
{"status": "error", "message": "...", "detail": "..."}
# 進度
{"status": "running", "progress": 0.68, "message": "..."}
```

---

## 9. React 前端規範

### 狀態管理（Zustand）

每個 Stage 的狀態集中在 `src/stores/pipelineStore.ts`：

```typescript
interface StageState {
  status: "idle" | "running" | "done" | "error";
  progress: number;
  logs: string[];
}
```

### API 呼叫

統一使用 `src/api/client.ts` 中的封裝函數，禁止在元件內直接 fetch：

```typescript
// src/api/client.ts
export const runRoiExtract = (config: RoiConfig) =>
  axios.post("/api/roi/extract", config);
```

### WebSocket Log

使用自訂 Hook `useStageLog(stage: string)` 訂閱即時 log：

```typescript
const logs = useStageLog("proseg");  // 自動連接 /ws/log/proseg
```

---

## 10. 資料安全規範

- **禁止執行 `cat` 大型生信檔案**（h5ad、zarr、BTF）
- **禁止使用 `pip`**，所有 Python 套件透過 `uv add`
- **所有參數必須從 `config/pipeline.yaml` 讀取**，禁止硬編碼
- **`results/` 目錄**：存放分析輸出，不納入 git
- **macOS 清理**：提交前執行 `find . -name "._*" -delete && find . -name ".DS_Store" -delete`

### Code Review 衍生規則（2026-03-04 更新）

#### DRY 原則：公用函數必須定義於模組層級

跨函數重複使用的邏輯（如型別轉換、解碼、格式化）**必須提升至模組層級定義，禁止在函數內部重複定義**。

```python
# ✅ 正確：模組層級定義，所有函數共用
def _decode_bytes(v) -> str:
    """將 10x H5 的 bytes 安全解碼為 str。"""
    if isinstance(v, bytes):
        try:
            return v.decode("utf-8")
        except UnicodeDecodeError:
            return v.decode("latin-1")
    return v

# ❌ 錯誤：在每個函數內分別定義 _decode_bytes / _dec / _decode...
```

#### subprocess 安全使用

所有 `subprocess.run()` 呼叫**禁止單獨使用 `check=False`**，必須搭配 returncode 檢查或 stderr 捕捉：

```python
# ✅ 正確
result = subprocess.run(["rm", "-rf", path], check=False, capture_output=True, text=True)
if result.returncode != 0:
    logger.warning(f"rm -rf 失敗：{result.stderr.strip()}")

# ❌ 錯誤：靜默失敗
subprocess.run(["rm", "-rf", path], check=False)
```

#### ROI 座標偏移後必須驗證

座標偏移後若出現負值代表 `roi_crop` 設定錯誤，需 warning：

```python
coords[:, 0] -= roi_crop["x0"]
neg_mask = coords[:, 0] < 0
if neg_mask.any():
    logger.warning(f"⚠️ ROI 偏移後出現 {neg_mask.sum()} 個負座標")
```

#### API 容錯：entry 層級 try/except

遍歷目錄或列表時，**單一項目的錯誤不應導致整個請求失敗**，應在 entry 層級 catch 並 continue：

```python
for entry in entries:
    try:
        ...
    except (PermissionError, OSError) as e:
        logger.warning(f"跳過 {entry.name}：{e}")
        continue
```

---

## 11. 測試規範

```bash
# 後端單元測試
uv run pytest backend/tests/ -v

# 前端測試
cd frontend && npm run test

# 端到端測試（小型 ROI）
uv run python backend/scripts/00_roi/extract_roi.py --dry-run
```

---

## 12. 效能注意事項

### 記憶體保護

- Visium HD 2µm 全圖：>100 萬 bins，需 backed mode 或先裁切
- BTF 影像：10-80 GB，**必須** tile-based 讀取
- Proseg GeoJSON：可達數 GB，**必須** ijson 串流

### GPU 資源管理

- Cellpose batch_size ≤ 8（RTX 4090 VRAM 限制）
- Proseg 條件測試並行 ≤ 4
- 每個 Stage 結束後顯式清理 GPU 快取：`torch.cuda.empty_cache()`

### Dask 分塊策略

- Zarr 讀寫：chunk size = (1, 1024, 1024) for images
- 大型 DataFrame：使用 dask.dataframe 延遲計算

---

## 13. 參考來源速查

| 需要移植的功能 | 來源路徑 |
|--------------|---------|
| ROI AnnData 裁切 | `xenium_visiumhd_comparison/scripts/02_baseline/cluster_8um.py` → `subset_to_roi()` |
| BTF tile-based 讀取 | `xenium_visiumhd_comparison/scripts/03_pipeline/segment_cellpose.py` → `_read_btf_crop()` |
| nucleus_boundaries 光柵化 | `xenium_visiumhd_comparison/scripts/03_pipeline/segment_proseg_zarr.py` → `rasterize_nucleus_mask()` |
| Xenium 轉錄本裁切 | `xenium_visiumhd_comparison/scripts/03_pipeline/segment_proseg_zarr.py` → `load_transcripts_roi()` |
| Xenium Explorer 匯出 | `Proseg-Zarr-Integration/scripts/export_to_xenium_full.py` |
| Loupe Browser 匯出 | `Proseg-Zarr-Integration/scripts/export_to_loupe_merged.py` |
| Proseg Pipeline | `visiumhd_pipeline/src/proseg/pipeline.py` |
| Cellpose + Logic A | `visiumhd_pipeline/src/segmentation/cellpose_runner.py` |
| Zarr Builder | `visiumhd_pipeline/scripts/02_build_zarr/create_zarr.py` |
| Scanpy 分析 | `visiumhd_pipeline/src/analysis/` |
