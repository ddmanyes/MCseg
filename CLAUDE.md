# MSseg — 開發維護規範

## 1. 專案定位

**MSseg** 是以 Pipeline 3 框架為基礎、整合 **MCseg v2** 分割引擎的 Visium HD 全流程分析平台。

MCseg v2 核心設計：**多模型 Voronoi 集成**——以 cyto3 三直徑（13/17/22px）+ 可選 Hematoxylin pass 進行多輪推論，透過 `merge_masks_fast` IoU 去重後套用 Voronoi 擴張，徹底取代舊版 LOGIC_A 雙尺寸策略。

**Pipeline 架構**：資料設定 → ROI 裁切 → MCseg v2 分割 → RNA 計數 → Scanpy 分析 → Browser 匯出

| Stage | 功能 | 關鍵點 |
|-------|------|--------|
| Setup | 資料自動掃描 | /api/data/scan |
| 0 | ROI 裁切（BTF tile-by-tile）| extractor.py |
| 1 | MCseg v2 分割 | segmentation_masks.npy（H×W int） |
| 2 | RNA 計數（Voronoi + 等距擴張）| `counter.py`，expand_labels(6px) → `cellpose_cells.h5ad` |
| 3 | Scanpy 分析 | QC + UMAP + Leiden |
| 4 | Browser 匯出 | skimage.find_contours → GeoJSON |

---

## 2. 環境管理

### Python 後端

```bash
# 使用 uv 管理（強制，禁止使用 pip）
uv sync                           # 安裝/同步依賴
uv add <package>                  # 新增套件
uv run python <script>            # 執行腳本
uv run uvicorn backend.main:app --reload --port 8001  # 啟動後端
```

> **ExFAT 注意**：`.venv` 是指向 SSD 上 venv 的 symlink。
> 若 symlink 失效，請重建：
>
> ```bash
> rm .venv && ln -s /Volumes/SSD/plan_a/visiumHD_pipeline_2_venv .venv
> find /Volumes/SSD/plan_a/visiumHD_pipeline_2_venv -name '._*' -delete
> UV_LINK_MODE=copy uv sync
> ```

### Node.js 前端

```bash
cd frontend
npm install
npm run dev    # port 3000
npm run build
```

---

## 3. 目錄結構規範

```
MSseg/
├── CLAUDE.md           # 本文件
├── config/
│   ├── pipeline.yaml          # 所有參數來源（禁止硬編碼）
│   └── profiles/              # Tissue Profile
│       ├── crc.yaml           # CRC：8 TME panels，MCseg v2 CRC 覆寫值
│       └── luad.yaml          # LUAD：4 TME panels，MCseg v2 LUAD 覆寫值
├── backend/
│   ├── main.py         # FastAPI 入口（port 8001）
│   └── src/
│       ├── api/
│       │   ├── data.py            # /api/data
│       │   ├── roi.py             # /api/roi
│       │   ├── segmentation.py    # /api/segmentation（MCseg v2 參數）
│       │   ├── cellpose_count.py  # /api/count
│       │   ├── analysis.py        # /api/analysis
│       │   └── export.py          # /api/export
│       ├── cellpose_counter/      # RNA 計數模組
│       │   └── counter.py
│       ├── roi/                   # Stage 0: ROI 裁切
│       ├── segmentation/          # Stage 1: MCseg v2
│       │   ├── cellpose_runner.py # MCseg v2 核心（run_mcseg_v2）
│       │   └── macenko.py         # 保留備用（MCseg v2 未使用）
│       ├── analysis/              # Stage 3: Scanpy
│       ├── export/                # Stage 4: Xenium/Loupe
│       └── utils/
│           └── constants.py       # 物理常數集中於此
└── frontend/src/
    ├── api/client.ts              # API 客戶端
    ├── stores/pipelineStore.ts    # Zustand
    └── pages/
        ├── DataSetup.tsx
        ├── Stage0_ROI.tsx
        ├── Stage1_Segmentation.tsx  # MCseg v2 參數 UI
        ├── Stage2_Count.tsx
        ├── Stage3_Analysis.tsx
        └── Stage4_Export.tsx
```

---

## 4. 物理常數（集中管理）

所有座標換算常數必須從 `backend/src/utils/constants.py` 引用，**嚴禁在腳本內硬編碼**：

```python
XENIUM_UM_PX   = 0.2125     # Xenium morphology 像素尺寸（µm/px）
XENIUM_NM_PX   = 212.5
PROSEG_UM_PX   = 0.2645833  # Proseg/Xenium Explorer 座標系
VISIUM_UM_PX   = 0.2737     # Visium HD fullres（µm/px）
```

---

## 5. ROI 裁切規範（Stage 0）

### 座標格式

```yaml
rois:
  - name: "test"
    tissue: "CRC"
    x: 48128
    y: 12657
    width_px: 1491
    height_px: 1210
    pixel_size_um: 0.2737
```

### 座標轉換公式

| 轉換 | 公式 |
|------|------|
| Visium fullres px → µm | `× 0.2737` |
| Visium fullres px → hires px | `× scalef` |

### BTF tile-based 讀取

大型 BTF/TIFF 必須使用 tile-based 讀取，**嚴禁全圖載入**。

---

## 6. MCseg v2 核心演算法（Stage 1）

`backend/src/segmentation/cellpose_runner.py` 核心函數：

### 主要函數

| 函數 | 說明 |
|------|------|
| `apply_clahe` | LAB 色彩空間局部對比增強 |
| `color_deconvolution_he` | Ruifrok & Johnston Hematoxylin 通道提取 |
| `voronoi_expand` | Voronoi 擴張（無重疊，distance capped） |
| `merge_masks_fast` | IoU 去重多遮罩集成 |
| `find_transcript_seeds` | 轉錄本密度補救（選配，需 vhd_pseudo_transcripts.csv） |
| `run_mcseg_v2` | 單張影像完整 MCseg v2 流程 |
| `run_segmentation_rois` | 多 ROI 批次執行 |
| `run_preview_patch` | 快速預覽（自動停用 cpsam 和 transcript rescue） |

### MCseg v2 參數鍵值（`_ROI_OVERRIDE_FIELDS`）

```python
use_gpu, batch_size,
dia_small, dia_mid, dia_large,
use_hematoxylin, use_cpsam,
voronoi_distance,
flow_threshold, cellprob_threshold,
min_size, max_size,
clahe_clip_limit, use_transcript_rescue
```

**每 ROI 覆寫**：`segmentation.py` API 的 `roi_overrides` 欄位可覆寫上述所有參數。

---

## 7. Stage 2 RNA 計數核心演算法

`backend/src/cellpose_counter/counter.py` 核心流程：

```python
# 1. 載入 MCseg v2 遮罩（Voronoi 擴張後）
mask = np.load("segmentation_masks.npy")

# 2. 等距擴張（填補 Voronoi 間隙）
from skimage.segmentation import expand_labels
if dilation_px > 0:
    mask = expand_labels(mask, distance=dilation_px)  # 預設 6px = 1.64µm

# 3. Bin 座標映射
coords = adata.obsm['spatial']
col = (coords[:, 0] - roi_x_px).astype(int)
row = (coords[:, 1] - roi_y_px).astype(int)

# 4. 稀疏矩陣聚合
A = csr_matrix(...)
count_matrix = A @ adata.X   # (n_cells, n_genes)
```

---

## 8. Browser 匯出規範（Stage 4）

### Xenium Explorer 匯出

```python
from skimage.measure import find_contours
# 從 MCseg v2 遮罩提取輪廓 → GeoJSON
```

### Loupe Browser 匯出

```python
# 10x_whitelist.txt 存於 backend/src/export/
# 分類欄位限制：唯一值 > 32,000 須轉整數
```

---

## 9. FastAPI 後端規範

### 長任務模式

```python
from fastapi import BackgroundTasks
@app.post("/api/segmentation/run")
async def run_seg(background_tasks: BackgroundTasks):
    background_tasks.add_task(run_mcseg_v2_pipeline)
    return {"status": "started"}
```

### API 回應格式

```python
{"status": "ok",      "data": {...}}
{"status": "error",   "message": "...", "detail": "..."}
{"status": "running", "progress": 0.68, "message": "..."}
```

---

## 10. React 前端規範

### MCseg v2 分割參數（Stage1_Segmentation.tsx）

介面分為左右兩欄：
- **左欄**：GPU 設定 + cyto3 多直徑（dia_small/mid/large + use_hematoxylin + use_cpsam）
- **右欄**：Voronoi/後處理 + Cellpose QC + 快速預設值

ROI override 表格欄位：name、cells、dia_mid（主徑）、voronoi_distance、overridden。

### 狀態管理（Zustand）

```typescript
interface StageState {
  status: "idle" | "running" | "done" | "error";
  progress: number;
  logs: string[];
}
```

---

## 11. Code Review 衍生規則

### DRY 原則

跨函數重複使用的邏輯必須提升至模組層級，禁止在函數內部重複定義。

### subprocess 安全使用

```python
# ✅ 正確
result = subprocess.run(cmd, check=False, capture_output=True, text=True)
if result.returncode != 0:
    logger.warning(f"失敗：{result.stderr.strip()}")
```

### ROI 座標偏移驗證

```python
coords[:, 0] -= roi_crop["x0"]
neg_mask = coords[:, 0] < 0
if neg_mask.any():
    logger.warning(f"⚠️ ROI 偏移後出現 {neg_mask.sum()} 個負座標")
```

### API 容錯：entry 層級 try/except

```python
for entry in entries:
    try:
        ...
    except (PermissionError, OSError) as e:
        logger.warning(f"跳過 {entry.name}：{e}")
        continue
```

---

## 12. 測試規範

```bash
uv run pytest backend/tests/ -v
cd frontend && npm run test
```

---

## 13. 效能注意事項

- Visium HD 2µm 全圖：>100 萬 bins，需 backed mode 或先裁切
- BTF 影像：10-80 GB，**必須** tile-based 讀取
- MCseg v2 batch_size ≤ 8（VRAM 限制）
- 每個 Stage 結束後顯式清理 GPU 快取：`torch.cuda.empty_cache()`

---

## 14. Tissue Profile 系統

換組織只需修改 `pipeline.yaml` 一行：

```yaml
global:
  tissue_profile: luad   # 改這一行
```

合併順序（後者優先）：

```text
config/profiles/{tissue}.yaml  ← 組織基底
config/pipeline.yaml            ← 專案設定
results/state.json              ← 執行期動態狀態
```

### Profile MCseg v2 可覆寫欄位

```yaml
segmentation:
  mcseg_v2:
    dia_small: 13.0
    dia_mid: 17.0
    dia_large: 22.0
    voronoi_distance: 9
    clahe_clip_limit: 3.0
    flow_threshold: 0.4
    cellprob_threshold: -2.0
    min_size: 20
    max_size: 6000
```

---

## 15. 資料安全規範

- **禁止執行 `cat` 大型生信檔案**（h5ad、BTF）
- **禁止使用 `pip`**，所有 Python 套件透過 `uv add`
- **所有參數必須從 `config/pipeline.yaml` 讀取**，禁止硬編碼
- **macOS 清理**：`find . -name "._*" -delete && find . -name ".DS_Store" -delete`

---

## 16. 前端 QC 視覺化規範

### 數值精度防護

```typescript
// 依據資料跨度決定捨入精度
const range = xMax - xMin
const rounded = range <= 10 ? Number(clamped.toFixed(3)) : Math.round(clamped)
```

---

## 17. 參考來源速查

| 需要移植的功能 | 來源路徑 |
|--------------|---------|
| MCseg v2 原始實作 | `autoresearch_seg/segment_best.py` |
| ROI AnnData 裁切 | `xenium_visiumhd_comparison/scripts/02_baseline/cluster_8um.py` |
| BTF tile-based 讀取 | `xenium_visiumhd_comparison/scripts/03_pipeline/segment_cellpose.py` |
| Xenium Explorer 匯出 | `Proseg-Zarr-Integration/scripts/export_to_xenium_full.py` |
| Loupe Browser 匯出 | `Proseg-Zarr-Integration/scripts/export_to_loupe_merged.py` |
| Scanpy 分析 | `visiumhd_pipeline/src/analysis/` |
