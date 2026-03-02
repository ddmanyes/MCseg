# visiumHD_pipeline_2 — 開發維護規範

## 1. 專案定位

**VisiumHD Pipeline 2** 是 `visiumhd_pipeline` 的下一代版本，整合了四項重大升級：

1. **Stage 0：ROI 裁切前處理** — 在分析前裁切感興趣區域
2. **Stage 2.5：Proseg 參數條件測試** — 自動評估最佳參數組合
3. **Stage 5：Browser 格式匯出** — 輸出 Xenium Explorer / Loupe Browser 相容格式
4. **React 前端** — 替換 PyQt6，使用現代 Web UI（FastAPI + React + Vite）

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

### Node.js 前端

```bash
cd frontend
npm install                       # 安裝依賴（或 pnpm install）
npm run dev                       # 開發模式（port 3000）
npm run build                     # 建構生產版本
```

### 一鍵啟動（開發模式）

```bash
# Terminal 1: 後端
uv run uvicorn backend.main:app --reload --port 8000

# Terminal 2: 前端
cd frontend && npm run dev
```

---

## 3. 目錄結構規範

```
visiumHD_pipeline_2/
├── CLAUDE.md           # 本文件
├── plan.md             # 實作計畫
├── config/
│   ├── pipeline.yaml   # 所有參數來源（禁止硬編碼）
│   └── roi_presets.yaml
├── backend/
│   ├── main.py         # FastAPI 入口
│   └── src/
│       ├── api/        # 路由端點
│       ├── roi/        # Stage 0: ROI 裁切
│       ├── segmentation/   # Stage 1
│       ├── zarr_builder/   # Stage 2
│       ├── proseg/     # Stage 2.5 + Stage 3
│       ├── analysis/   # Stage 4
│       ├── export/     # Stage 5（Xenium / Loupe）
│       └── utils/
│           └── constants.py  # 所有物理常數集中於此
└── frontend/
    └── src/
        ├── api/        # API 呼叫層
        ├── components/ # 可重用元件
        └── pages/      # 各 Stage 頁面
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

## 6. Proseg 相關規範

### 條件測試（Stage 2.5）

- 測試用小型 ROI（建議 500µm × 500µm 或 1000µm × 1000µm）
- 最多同時並行 4 個條件（避免 VRAM 耗盡）
- 使用獨立臨時目錄（`results/02b_conditions/cond_{idx}/`）

### macOS `._*` 污染防護

每次讀取 Zarr 或 Parquet 前必須清理：

```python
import pathlib
for f in pathlib.Path(zarr_dir).rglob("._*"):
    f.unlink(missing_ok=True)
```

### Dask Monkey Patch

在任何 spatialdata 操作前必須確保 query-planning 啟用：

```python
import dask
dask.config.set({"dataframe.query-planning": True})
```

---

## 7. Browser 匯出規範（Stage 5）

### Xenium Explorer 匯出

**必須修補 pixel_size Bug**（`spatialdata_xenium_explorer` 硬編碼 0.2125）：

```python
# 每次 write() 後都要執行
import json
exp_file = Path(out_dir) / "experiment.xenium"
with open(exp_file, "r+") as f:
    data = json.load(f)
    data["pixel_size"] = PROSEG_UM_PX  # 0.2645833
    f.seek(0)
    json.dump(data, f, indent=2)
    f.truncate()
```

**多邊形平滑化**（必須，否則 Xenium Explorer 渲染慢）：

```python
from shapely.affinity import scale

poly_px = scale(poly_um, xfact=1/PROSEG_UM_PX, yfact=1/PROSEG_UM_PX, origin=(0,0))
poly_smooth = poly_px.simplify(0.4, preserve_topology=True)
poly_smooth = poly_smooth.buffer(1.0, join_style=1).buffer(-1.0, join_style=1)
```

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
