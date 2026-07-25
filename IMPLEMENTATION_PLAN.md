# MSseg 實作計畫 — 全圖分割串接 / 區域分析 / 對位配準 / NDPI 支援

> 建立日期：2026-07-25
> 前置：程式碼探勘已完成（見下方「現況事實」）。本計畫零佔位符，所有設計決策已定案。
> 測試指令一律使用 `.venv/bin/python -m pytest`（**禁用** `uv run pytest`，會摧毀 `.venv` symlink）。

---

## ⛔ P-PRE — 開工前的阻斷性前置（必須先做）

**本地 master 與 origin/master 已分歧：本地領先 1、落後 5。** Windows 那台於 2026-07-25 推了 5 個 commit（`347a964` Windows clone 修檔名、`cdcf247` QC 指標收斂、`cce57e1`+`70b6352` 匯出下沉、`04be6ab` roi_overrides 驗證），本地完全沒有。

- [x] **P-PRE-1 同步遠端並重放本地 commit** ✅ 2026-07-25
  - 預期行為：`git fetch origin` 後把本地唯一的 `1d2fa3f`（xenium_outs 測試對齊）重放到 `origin/master` 之上。
  - 驗證：`git rev-list --left-right --count master...origin/master` → `1  0`；`.venv/bin/python -m pytest backend/tests/ -q` 全綠
  - **結果**：rebase 後本地 commit 成為 `7f9163d`，疊在 `04be6ab` 之上，無衝突。全套 **106 passed**（本地 62 + 遠端新增 44），4m04s。
  - **坑**：`git rebase origin/master` 首次執行回報 `Aborting / Could not execute the todo command`（git 試圖開編輯器，此環境不支援互動）。解法：`GIT_EDITOR=true GIT_SEQUENCE_EDITOR=true git rebase --continue`。
  - commit：無（rebase）

- [x] **P-PRE-2 重新校準本計畫的行號引用** ✅ 2026-07-25
  - 預期行為：`04be6ab` 動過 `api/segmentation.py`（移除 `_ROI_OVERRIDE_FIELDS` import）、`cce57e1`/`70b6352` 把 `api/export.py` 由 679 砍到 200 行。**本計畫「現況事實」中的行號是 rebase 前的**，P0-8/P0-10/P0-11 動工前需重新 grep 確認。
  - 驗證：`grep -n "max_load_gb\|estimated_gb\|use_cpsam.*False\|full_image_segmentation_masks" backend/src/api/segmentation.py`
  - **結果**：本計畫引用的行號**全部未變** —— `estimated_gb > 6.0` 仍在 `:258`、`np.array(arr)` 在 `:264`、`use_cpsam=False` 在 `:275`、輸出路徑在 `:287`；`cli/segment.py` 的 `step_bin_attribution:167` / `_aggregate_cells_raw:199` / `VISIUM_UM_PX:535-536` 亦未變。`04be6ab` 的改動落在 `:315` 附近（本計畫關注區段之後），故無須改寫。
  - commit：`docs(plan): 校準行號引用至 origin/master`

> **測試檔編號衝突（已修正）**：遠端 5 個 commit 已佔用 `test_06_qc_metrics.py`、`test_07_export_inputs.py`、`test_08_export_jobs.py`、`test_09_roi_overrides.py`（共 44 項新測試）。本計畫的新測試檔已全部改為 **test_10 ~ test_13**。

---

## 目標摘要

| # | 目標 | 現況 | 完成後 |
|---|------|------|--------|
| P0 | 全圖分割接上 pipeline | UI 按鈕存在但輸出是死路，Stage 2 讀不到 | 全圖遮罩可直接產出 `cells.h5ad`，UI 可傳裁切座標 |
| P1 | 對位問題可視化 | 只有「>30% bins 超出範圍」的粗警告 | QC 疊圖 + 次像素位移估計數值 |
| P2 | 區域選取做圖分析 | 區域選取只存在於 Stage 0（分割前） | 全圖遮罩上可無限次框選區域做圖，不需重跑分割 |
| **P0.5** | **對位 JSON 支援** | **不讀 Loupe/CytAssist JSON，只用近似縮放** | **組合 homography 精確對位（實測多命中 13.4% bins）** |
| P3 | Registration 模組 | 無，attribution 假設座標系相同 | 殘餘位移估計（**降為驗證用**，主對位交給 P0.5） |
| P4 | NDPI 讀取 | 零支援，strip 解析邏輯綁死 BigTIFF | `SlideReader` 抽象層，支援 .ndpi/.svs/.mrxs |
| P5 | 串流 tiled 分割 | 整圖進 RAM，>6 GB 直接拒絕 | tile 串流 + memmap label map，不受 RAM 限制 |

### 現況事實（探勘結論，作為計畫依據）

- `/api/segmentation/run_full` 與 `/full_seg_status` **已存在**（`segmentation.py:207-305`），前端 `client.ts:37-38` + `Stage1_Segmentation.tsx:1004-1020` 已接線。
- 致命缺口：`segmentation.py:264` `img = np.array(arr)` 整圖進 RAM；`:258` 硬擋 6 GB；`:287` 輸出 `full_image_segmentation_masks.npy`，但 Stage 2 只找 `{output_dir}/roi/{name}/segmentation_masks.npy`（`counter.py:241`）→ **無法接續**。
- `segmentation.py:275` 強制 `use_cpsam=False`，UI 只能 4-pass，CLI 可 7-pass。
- CLI（`cli/segment.py`）已有完整全片流程且每步 `[SKIP]` 續跑：`step_bin_attribution`（:167）、`_aggregate_cells_raw`（:199）、`step_count_cells`（:255）。**P0 以 DRY 方式提升這些函式，不重寫。**
- `tile_server.py:8-14` 已記載 Space Ranger fullres 座標系 ≠ raw TIFF，但**只修了檢視器縮圖，未修 bin attribution**（`segment.py:188-189` 仍直接 `pxl_row_in_fullres - crop_y0`）。
- `tifffile 2026.2.24` 內含 `_ndpi_load_pages` / `_series_ndpi` → **原生支援 NDPI**；`tiffslide`、`openslide` 均未安裝。
- Stage 2 需要每 ROI 目錄同時有 `adata_002um.h5ad` 與 `segmentation_masks.npy`，並在 config `rois` 有 `x`/`y`/`pixel_size_um`。
- `results/state.json` 目前為空 `{}`，ROI 來自 `config/pipeline.yaml`。

---

## 文件架構圖

```text
MSseg/
├── IMPLEMENTATION_PLAN.md                    ← 本文件（新增）
├── config/pipeline.yaml                      ← 修改：新增 full_seg / alignment 區塊
├── pyproject.toml                            ← 修改：uv add tiffslide
├── backend/
│   ├── src/
│   │   ├── fullslide/                        ← 【新增目錄】P0
│   │   │   ├── __init__.py
│   │   │   └── pipeline.py                   ← 由 cli/segment.py 提升的共用函式
│   │   ├── registration/                     ← 【新增目錄】P1+P3
│   │   │   ├── __init__.py
│   │   │   ├── align.py                      ← 位移估計 + AffineAlignment
│   │   │   └── qc.py                         ← QC 疊圖產生
│   │   ├── utils/
│   │   │   └── slide_reader.py               ← 【新增】P4 SlideReader 抽象層
│   │   ├── api/
│   │   │   ├── segmentation.py               ← 修改：裁切座標、解除硬編碼、P5 串流
│   │   │   ├── cellpose_count.py             ← 修改：新增 /run_full、/full_status
│   │   │   ├── spatial.py                    ← 修改：region/polygon 過濾
│   │   │   └── registration.py               ← 【新增】P1 /api/registration/*
│   │   ├── cli/segment.py                    ← 修改：改 import fullslide.pipeline（DRY）
│   │   ├── roi/tile_server.py                ← 修改：P4 改用 SlideReader
│   │   └── segmentation/cellpose_runner.py   ← 修改：P5 tile_reader + memmap
│   └── tests/
│       ├── test_10_fullslide.py              ← 【新增】P0
│       ├── test_11_registration.py           ← 【新增】P1+P3
│       ├── test_12_region.py                 ← 【新增】P2
│       └── test_13_slide_reader.py           ← 【新增】P4
└── frontend/src/
    ├── api/client.ts                         ← 修改：新增端點
    ├── components/RegionSelector.tsx          ← 【新增】P2（由 Stage0 畫框邏輯提取）
    ├── i18n/translations.ts                  ← 修改：新增字串
    └── pages/
        ├── Stage1_Segmentation.tsx            ← 修改：裁切座標欄位
        ├── Stage2_Count.tsx                  ← 修改：全圖計數按鈕
        └── Stage35_SpatialExplorer.tsx        ← 修改：接 RegionSelector
```

---

## P0 — 全圖分割接上 pipeline（預估 4-6 h）

> 核心決策：**不**用 pseudo-ROI hack（會被迫為全片產生 `adata_002um.h5ad` 與 `he_crop.tif`，又撞 RAM 上限）。改為把 CLI 已驗證的 attribution 流程提升為共用模組，API 與 CLI 共用（CLAUDE.md §11 DRY）。

- [x] **P0-1 建立 fullslide 模組骨架**
  - 預期行為：`backend/src/fullslide/__init__.py`、`pipeline.py` 建立，`pipeline.py` 僅含模組 docstring。
  - 驗證：`.venv/bin/python -c "import backend.src.fullslide.pipeline"`
  - 檔案：`backend/src/fullslide/{__init__,pipeline}.py`
  - commit：`feat(fullslide): 建立全片流程模組骨架`

- [x] **P0-2 【紅燈】寫 bin_attribution 失敗測試**
  - 預期行為：`test_10_fullslide.py::test_bin_attribution_maps_bins_to_cells` 以合成資料（10×10 mask，label 1 佔 rows 0-4、label 2 佔 rows 5-9；4 個假 bin）斷言回傳 DataFrame 含 `barcode`/`cell_id` 且 cell_id 正確。此時 import 失敗。
  - 驗證：`.venv/bin/python -m pytest backend/tests/test_10_fullslide.py -q` → **FAILED（ImportError）**
  - 檔案：`backend/tests/test_10_fullslide.py`
  - commit：`test(fullslide): bin_attribution 紅燈測試`

- [x] **P0-3 【綠燈】搬移 bin_attribution 至 fullslide.pipeline**
  - 預期行為：把 `cli/segment.py:167-196` 的 `step_bin_attribution` 主體改寫為 `bin_attribution(mask, tp_path, crop_y0, crop_x0, out_path=None) -> pd.DataFrame`（快取邏輯改為 `out_path` 可選）。純函式、不含 log bar。
  - 驗證：同上指令 → **PASSED**
  - 檔案：`backend/src/fullslide/pipeline.py`
  - commit：`feat(fullslide): bin_attribution 純函式化`

- [x] **P0-4 【重構】CLI 改用共用函式**
  - 預期行為：`cli/segment.py` 的 `step_bin_attribution` 改為薄包裝（保留 `[SKIP]` log 與 parquet 快取），內部呼叫 `fullslide.pipeline.bin_attribution`。刪除重複邏輯。
  - 驗證：`.venv/bin/python -m pytest backend/tests/ -q` 全綠；`.venv/bin/python -m backend.src.cli.segment --help` 正常輸出
  - 檔案：`backend/src/cli/segment.py`
  - commit：`refactor(cli): bin_attribution 改用 fullslide 共用模組`

- [x] **P0-5 【紅燈】aggregate_cells 測試**
  - 預期行為：`test_aggregate_cells_sums_bins_per_cell` 用 3 bins→2 cells 的合成 h5（`anndata` 就地建立、寫臨時 h5ad 再讀）斷言輸出 `n_obs==2`、`obs['n_bins']==[2,1]`、counts 為加總。
  - 驗證：`.venv/bin/python -m pytest backend/tests/test_10_fullslide.py -q` → FAILED
  - 檔案：`backend/tests/test_10_fullslide.py`
  - commit：`test(fullslide): aggregate_cells 紅燈測試`

- [x] **P0-6 【綠燈】搬移 aggregate_cells**
  - 預期行為：`_aggregate_cells_raw`（`cli/segment.py:199-252`）搬為 `aggregate_cells(attribution, h5_path) -> AnnData`，簽章不變、移除 log。
  - 驗證：同上 → PASSED
  - 檔案：`backend/src/fullslide/pipeline.py`、`backend/src/cli/segment.py`
  - commit：`feat(fullslide): aggregate_cells 共用化`

- [x] **P0-7 加入 centroid 計算函式（含 fullres 座標）**
  - 預期行為：`add_centroids(cells, mask, pixel_size_um, origin_xy=(0, 0))` 就地補 `obs['centroid_x_px']`、`obs['centroid_y_px']`（**裁切局部**，維持 CLI 現行語意）、以及 `obs['centroid_x_fullres']`、`obs['centroid_y_fullres']`（局部 + `origin_xy`）、`obsm['spatial']`（µm，由局部座標換算）。邏輯取自 `cli/segment.py:276-286`。
  - **為何需要 fullres 欄位**：P2 的 RegionSelector 在全片座標系上框選，而 centroid 是裁切局部座標。若不在此處落地 fullres 欄位，P2 每次過濾都要回頭讀 meta sidecar 補償原點 —— 那是重複且易錯的轉換。
  - 驗證：`test_add_centroids_sets_local_and_fullres`：10×10 全 label-1 mask、`origin_xy=(100, 200)` → local ≈ (4.5, 4.5)、fullres ≈ (104.5, 204.5)
  - 檔案：`backend/src/fullslide/pipeline.py`
  - commit：`feat(fullslide): add_centroids 共用化`

- [x] **P0-8 解除 6 GB 與 cpsam 硬編碼**
  - 預期行為：`segmentation.py:258` 的 `6.0` 改讀 `config["full_seg"]["max_load_gb"]`（預設 6.0）；`:275` 的 `use_cpsam=False` 改讀 `config["full_seg"]["force_disable_cpsam"]`（預設 `true`，維持現行為）。`config/pipeline.yaml` 新增 `full_seg:` 區塊含 `tile_size: 1024`、`overlap: 128`、`max_load_gb: 6.0`、`force_disable_cpsam: true`。
  - 驗證：`.venv/bin/python -m pytest backend/tests/test_01_infra.py -q`；`.venv/bin/python -c "from backend.src.utils.config import load_config; c=load_config(); print(c['full_seg'])"`
  - 檔案：`backend/src/api/segmentation.py`、`config/pipeline.yaml`
  - commit：`refactor(segmentation): full_seg 參數移出硬編碼`

- [x] **P0-9 【紅燈】裁切座標驗證測試**
  - 預期行為：`test_full_seg_crop_validation` POST `/api/segmentation/run_full` 帶 `{"crop_x0":100,"crop_x1":50}` → 回 `status:error`、訊息含「crop_x1 必須大於 crop_x0」。
  - 驗證：`.venv/bin/python -m pytest backend/tests/test_10_fullslide.py -q` → FAILED
  - 檔案：`backend/tests/test_10_fullslide.py`
  - commit：`test(segmentation): 全圖裁切座標驗證紅燈測試`

- [x] **P0-10 【綠燈】run_full 接受裁切座標**
  - 預期行為：新增 `FullSegParams(BaseModel)`：`crop_y0/crop_y1/crop_x0/crop_x1: Optional[int] = None`、`use_cpsam: Optional[bool] = None`。`-1` 或 `None` = 全圖邊界。驗證 `0 <= v0 < v1 <= 影像邊界`，否則回 error。`_run_full_segmentation` 只讀取該窗格（`arr[y0:y1, x0:x1]`），並把 `(x0, y0)` 存入結果 metadata。
  - 驗證：同上 → PASSED
  - 檔案：`backend/src/api/segmentation.py`
  - commit：`feat(segmentation): run_full 支援裁切座標與 cpsam 開關`

- [x] **P0-11 全圖分割輸出加 metadata sidecar**
  - 預期行為：分割完成後除 `.npy` 另寫 `full_image_segmentation_meta.json`：`{crop_x0, crop_y0, width, height, n_cells, pixel_size_um, passes, created_at}`。下游據此還原座標。
  - 驗證：`test_full_seg_meta_schema` 斷言 json 含全部 8 個 key
  - 檔案：`backend/src/api/segmentation.py`
  - commit：`feat(segmentation): 全圖遮罩輸出 metadata sidecar`

- [x] **P0-12 【紅燈】/api/count/run_full 端點測試**
  - 預期行為：`test_count_run_full_requires_mask`：無遮罩時 POST `/api/count/run_full` → `status:error` 且訊息含「請先完成全圖分割」。
  - 驗證：`.venv/bin/python -m pytest backend/tests/test_10_fullslide.py -q` → FAILED
  - 檔案：`backend/tests/test_10_fullslide.py`
  - commit：`test(count): 全圖計數端點紅燈測試`

- [x] **P0-13 【綠燈】實作 /api/count/run_full + /full_status**
  - 預期行為：背景任務讀 `full_image_segmentation_masks.npy` + meta sidecar，依序呼叫 `bin_attribution` → `aggregate_cells` → `add_centroids`，輸出 `{output_dir}/fullslide/cells.h5ad`。`tp` 取 `paths.binned_002/spatial/tissue_positions.parquet`，`h5` 取 `paths.binned_002/filtered_feature_bc_matrix.h5`。進度回報沿用 `_full_status` 同型別 dict。
  - 驗證：同上 → PASSED
  - 檔案：`backend/src/api/cellpose_count.py`
  - commit：`feat(count): 新增全圖 RNA 計數端點`

- [x] **P0-14 per-sample mpp 覆寫常數**
  - 預期行為：新增 `fullslide.pipeline.resolve_pixel_size(config) -> float`：優先讀 `binned_002/spatial/scalefactors_json.json` 的 `microns_per_pixel`，缺失才回退 `VISIUM_UM_PX`。CLI `segment.py:535-536` 與 P0-13 改用它。**註記**：此值只影響 µm 換算與匯出比例尺，**不影響** attribution（純像素運算）。
  - 驗證：`test_resolve_pixel_size_prefers_scalefactors`（用臨時 json 斷言取到自訂值；檔案不存在時取 0.2737）
  - 檔案：`backend/src/fullslide/pipeline.py`、`backend/src/cli/segment.py`
  - commit：`fix(fullslide): 樣本 microns_per_pixel 覆寫預設常數`

- [x] **P0-15 前端：裁切座標與全圖計數 UI**
  - 預期行為：`client.ts` 新增 `runFullSegmentation(body)` 參數化、`runFullCount()`、`getFullCountStatus()`。`Stage1_Segmentation.tsx` 全圖區塊加 4 個座標輸入（空白 = 全圖）與 cpsam 勾選；`Stage2_Count.tsx` 加「全圖計數」按鈕 + 進度條。i18n 補 zh/en 字串。
  - 驗證：`cd frontend && npm run build`
  - 檔案：`frontend/src/api/client.ts`、`frontend/src/pages/Stage1_Segmentation.tsx`、`frontend/src/pages/Stage2_Count.tsx`、`frontend/src/i18n/translations.ts`
  - commit：`feat(ui): 全圖分割裁切座標與全圖計數介面`

---

### P0 完成紀錄（2026-07-25）

**全部 15 項完成**，7 個 commit（`05df83a` → `ea94214`）。backend 全套 **141 passed**（P0 新增 35 項，全為合成資料、Windows 可跑）；`npm run build`（含 `tsc`）通過。

實作過程與計畫的三處偏離：

1. **`-1` 哨兵語意收窄**（P0-9 測試抓出來的）：初稿把 `-1` 當成「取影像邊界」對四個座標一律適用，但這讓 `crop_x0=-1` 這種明顯錯誤變成合法輸入。改為 **`-1` 只對上界（`crop_x1`/`crop_y1`）有效**，下界收到負值即報錯 —— 與 CLI 的 `--crop-y0` 預設 0、`--crop-y1` 預設 -1 一致。
2. **記憶體上限改以裁切窗格計算**：原本 `estimated_gb` 算的是整張影像，加了裁切後若仍用全圖尺寸判斷，縮小範圍也會被擋。改成算實際窗格，超限訊息也從「請改用 ROI 模式」改為「請縮小裁切範圍或改用 ROI 模式」。
3. **多做了 `resolve_full_count_inputs`**：計畫的 P0-13 把輸入解析寫在端點內；實作時抽成 `fullslide.pipeline` 的純函式，才能在不啟動背景任務的情況下測錯誤路徑（4 項測試）。錯誤訊息一律不含絕對路徑。

sidecar 缺失的處理也在實作時定案為 **warning 而非 error**（原點退為 `(0,0)`），因為舊遮罩沒有 sidecar，而全圖模式的原點本來就是 `(0,0)`。

### P0 完成後補修：座標系縮放與等距擴張（commit `21fc26c`）

P0 收尾後查 EP 既有全片紀錄（`f981cda1`，dpcp01_vh_v114_02_hd_r004），發現 P0 交付的 `bin_attribution` 有**三個會靜默出錯**的缺陷，已修並用 EP 結果驗證：

1. **缺 SR fullres → 遮罩(TIFF) px 的縮放**。Space Ranger 的 fullres 校準到餵給它的影像，未必等於分割用的 raw TIFF。dpcp01：SR 0.5464 µm/px vs TIFF ~0.2737 → 差近 2 倍，原本每個 bin 都落在約一半位置。CRC 官方樣本因 mpp 恰為 0.2737（scale=1）巧合正確，所以一直沒被發現 —— **這就是「有時候需要 Space Ranger 重定位」的成因**。
2. **越界 bin 被 `.clip(0, h-1)` 夾到邊緣**，把界外 RNA 誤記到邊界細胞且不留痕跡 → 改為排除 + 越界 >30% 時警告。
3. **全圖路徑漏了 `expand_labels`**（ROI 路徑 `counter.py` 有）→ 兩者不可比，改為共用 `rna_counting.dilation_px`。

驗證鏈（EP 基準 1,545,019 bins / 35.3%）：

| 設定 | bins | 細胞 |
|------|------|------|
| `scale=1` 無擴張（修正前） | 746,020 (17.0%) | 9,653 |
| ＋自動推導 scale | 914,338 (20.9%) | 43,066 |
| EP 自己的 `vfr_CORRECTED` 遮罩 ＋ `scale=1` | 914,909 (20.9%) | 43,068 |
| ＋`dilation 6px` | **1,545,019 (35.3%)** | 43,212 |

中間兩列為**兩條獨立路徑交叉驗證**（差 2 個細胞為捨入）。注意 EP 摘要的「86,314 cells」是遮罩總分割細胞數（max label），非拿到 RNA 的細胞數。

---

## P0.5 — Loupe / CytAssist 對位 JSON 支援（預估 半天，**優先於 P1**）

> **本節取代 P3-3 / P3-4 的設計。** 那兩項假設「預設單位矩陣 + 估小位移」，實測不成立。

### 為何必須做（附量測證據）

使用者實際工作流程：**outs 內的高解析圖不夠清晰 → 另外輸出更高解析的圖 → 用 Loupe Browser 重新對位產生 JSON**。要讓 RNA 貼合那張新圖，必須讀該 JSON。

⚠️ **同時修正上面那條「自動推導 scale」的評價**：`resolve_bin_to_mask_scale`（hires 尺寸 ÷ scalef）在 **SR 畫布相對影像有 padding 時幾何上是錯的** —— 它把 bin 橫向壓縮進影像寬度。EP 的 `vfr_CORRECTED` resample 用同一慣例，所以兩者逐位元吻合卻**一起偏離真實幾何**，因此 **EP 的 35.3% 也偏低，不該當天花板**。

JSON 結構（`*_alignment_file.json` 與 `*-fiducials-image-registration.json` schema 相同），各帶兩個 3×3 homography：

| 欄位 | 映射 | dpcp01 實測 |
|------|------|------------|
| `transform` | slide µm → CytAssist px | scale 0.2158（CytAssist mpp 4.635） |
| `cytAssistInfo.transformImages` | 影像 px → CytAssist px | 0.1179（低解析）/ 0.0590（高解析），rot +90° 含鏡射 |

**已驗證的推導式**：`mpp_image = scale(transformImages) / scale(transform)`
dpcp01 無後綴 → 0.5465（scalefactors 記 0.5464 ✓）；`_0105` → 0.2732（＝ 47104×21504 高解析 TIFF）。
**dpcp01 的 `spatial/` 已存在一組 old + new 實例**，該流程等於做過。

正確變換：`高解析圖 px = inv(H_new) @ H_old @ (SR fullres px)`
dpcp01 組合結果：**等向 scale 2.00001、rot −0.0001°、平移 (−0.90, −0.85) px**。

量測比較（TIFF 空間遮罩，無 dilation）：

| 變換來源 | bins | 細胞 |
|---------|------|------|
| hires 推導 `(1.9088, 2.0)`（已 ship） | 914,347 (20.88%) | 43,066 |
| **JSON 組合 homography** | **1,036,566 (23.67%)** | **50,598** |

B 即使損失右緣落在影像外的 bin，仍多命中 **13.4% bins / 17.5% 細胞**。

`transformImages` 的 rot +90° 與鏡射由 Space Ranger 產生 `pxl_*_in_fullres` 時已套入，**不需再套**；但若有人拿方向不同的影像分割，scale-only 會靜默失敗，JSON 是唯一偵測依據。

### 任務

- [ ] **P0.5-1 【紅燈】load_alignment 測試**
  - 預期行為：`test_11_registration.py::test_load_alignment_derives_mpp`：以合成 JSON（`transform` scale 0.2、`transformImages` scale 0.1）斷言 `mpp == 0.5`、且回傳 `serial_number` / `area`。
  - 驗證：`.venv/bin/python -m pytest backend/tests/test_11_registration.py -q` → FAILED
  - commit：`test(registration): 對位 JSON 解析紅燈測試`

- [ ] **P0.5-2 【綠燈】實作 load_alignment**
  - 預期行為：`registration/alignment.py` 的 `load_alignment(path) -> Alignment`（dataclass：`transform` 3×3、`transform_images` 3×3、`mpp`、`serial_number`、`area`、`checksum`）。`mpp` 由上述推導式算出。損壞/缺鍵 raise `ValueError`（含檔名但不含完整路徑）。
  - 驗證：同上 → PASSED
  - 檔案：`backend/src/registration/alignment.py`
  - commit：`feat(registration): 解析 Loupe/CytAssist 對位 JSON`

- [ ] **P0.5-3 pick_source_alignment 自動選出 H_old**
  - 預期行為：`pick_source_alignment(paths, target_mpp, tol=0.01) -> (source, others)`：推導 mpp 與 `scalefactors.microns_per_pixel` 相對誤差 < tol 者即為 `H_old`。**同時解掉兩個陷阱**：EP 踩過的「套用不屬於本樣本的註冊檔」，以及 dpcp01 的「同 serial 兩版本 scale 差正好 2 倍」。無命中則 raise 並列出各候選的 mpp。
  - 驗證：`test_pick_source_alignment_matches_scalefactors_mpp`（兩份合成 JSON mpp 0.5465 / 0.2732，target 0.5464 → 選中前者）
  - commit：`feat(registration): 依 scalefactors mpp 自動選出來源對位檔`

- [ ] **P0.5-4 provenance 驗證**
  - 預期行為：`validate_pair(h_old, h_new)`：`serial_number` 或 `area` 不一致時 raise（訊息含兩邊的值）。防止把別片玻片的對位檔套進來。
  - 驗證：`test_validate_pair_rejects_serial_mismatch`
  - commit：`feat(registration): 對位檔 provenance 驗證`

- [ ] **P0.5-5 compose_bin_to_image**
  - 預期行為：`compose_bin_to_image(h_old, h_new) -> np.ndarray`（3×3，正規化 `M[2,2]=1`）＝ `inv(h_new.transform_images) @ h_old.transform_images`。
  - 驗證：`test_compose_recovers_isotropic_2x`：用 dpcp01 的真實兩份 JSON（`@pytest.mark.skipif` 檔案不存在時跳過），斷言 scale ≈ 2.0（±1e-3）、rot ≈ 0（±0.01°）
  - commit：`feat(registration): 組合 homography 得 bin→影像 變換`

- [ ] **P0.5-6 bin_attribution 改吃 3×3 homography**
  - 預期行為：新增 `transform: np.ndarray | None = None` 參數。給定時以 homography 映射 bin 座標（齊次除法），否則沿用現有 `scale` 對角特例 —— **不破壞既有呼叫端**。`scale` 與 `transform` 同時給時 `transform` 優先並記 warning。
  - 驗證：`test_bin_attribution_with_homography_equals_scale_for_diagonal`（對角 homography 結果須與 `scale` 路徑完全相同）
  - 檔案：`backend/src/fullslide/pipeline.py`
  - commit：`feat(fullslide): bin attribution 支援 3×3 homography`

- [ ] **P0.5-7 涵蓋率回報**
  - 預期行為：`bin_attribution` 回傳的 DataFrame 加 `.attrs["coverage"]`：`{n_total, n_in_bounds, n_assigned, frac_out_of_image}`。全圖計數完成訊息附「X% bins 落在影像範圍外」。
  - **為何需要**：dpcp01 正確變換下，SR 右緣約 980 TIFF px 寬的 bin 確實在高解析圖之外（SR 畫布 11266×2 = 22532 > TIFF 寬 21504；高度 23552×2 = 47104 完全吻合）。這是真實限制而非 bug，**必須明確回報而非靠壓縮硬塞**。
  - 驗證：`test_coverage_attrs_reports_out_of_image_fraction`
  - commit：`feat(fullslide): 回報 bin 影像涵蓋率`

- [ ] **P0.5-8 接進全圖計數流程並降級 fallback**
  - 預期行為：`resolve_full_count_inputs` 掃 `binned_002/spatial/*.json` 與 `alignment.extra_alignment_json`（新設定，指向 Loupe 新產生的檔），能組出 homography 就用；否則回退 `resolve_bin_to_mask_scale` 並 **log 標明為近似值**。
  - 驗證：`.venv/bin/python -m pytest backend/tests/ -q` 全綠
  - commit：`feat(count): 全圖計數優先採用對位 JSON`

> **P3 影響**：P3-3/P3-4 的 `alignment.matrix` 改為「JSON 缺失時的人工覆寫」；P1 的殘餘位移估計從**主要對位手段降為驗證手段**（範圍縮小到「JSON 正確但仍有殘餘偏移」的情況）。

---

## P1 — 對位可視化與位移估計（預估 3-4 h，最高性價比）

> 目的：先能**看見與量化**對位誤差，再談自動校正。20-50 px 的偏移不會觸發現有 30% 警告，卻會讓 RNA 落到隔壁細胞。

- [ ] **P1-1 建立 registration 模組骨架**
  - 預期行為：`backend/src/registration/{__init__,align,qc}.py` 建立。
  - 驗證：`.venv/bin/python -c "import backend.src.registration.align, backend.src.registration.qc"`
  - commit：`feat(registration): 建立配準模組骨架`

- [ ] **P1-2 【紅燈】bin 密度光柵化測試**
  - 預期行為：`test_render_bin_density_shape_and_mass`：4 個假 bin（含 1 個 `in_tissue=0` 應被排除）→ 斷言輸出 shape 正確、總和 == 3。
  - 驗證：`.venv/bin/python -m pytest backend/tests/test_11_registration.py -q` → FAILED
  - 檔案：`backend/tests/test_11_registration.py`
  - commit：`test(registration): bin 密度光柵化紅燈測試`

- [ ] **P1-3 【綠燈】實作 render_bin_density**
  - 預期行為：`render_bin_density(tp_path, full_shape, downsample) -> np.ndarray`：讀 parquet 的 `pxl_row/col_in_fullres` + `in_tissue`，只留 `in_tissue==1`，以 `np.add.at` 累加到 `(H//ds, W//ds)` float32 光柵。
  - 驗證：同上 → PASSED
  - 檔案：`backend/src/registration/align.py`
  - commit：`feat(registration): bin 密度光柵化`

- [ ] **P1-4 【紅燈】位移估計測試**
  - 預期行為：`test_estimate_shift_recovers_known_offset`：造一張 128×128 有結構的圖，`np.roll` 位移 (7, -5)，斷言 `estimate_shift` 回傳 dy/dx 誤差 < 0.5 px。
  - 驗證：`.venv/bin/python -m pytest backend/tests/test_11_registration.py -q` → FAILED
  - commit：`test(registration): 位移估計紅燈測試`

- [ ] **P1-5 【綠燈】實作 estimate_shift**
  - 預期行為：`estimate_shift(ref, mov, upsample_factor=10) -> tuple[float, float, float]` 回傳 `(dy, dx, error)`，單位為**輸入圖的像素**。用 `skimage.registration.phase_cross_correlation`；兩張圖先各自 z-score 標準化以抵抗亮度差（H&E 灰階 vs bin 密度量級差異極大，未標準化會失準）。
  - 驗證：同上 → PASSED
  - 檔案：`backend/src/registration/align.py`
  - commit：`feat(registration): 次像素位移估計`

- [ ] **P1-5b 兩階段（coarse-to-fine）估計**
  - 預期行為：`estimate_shift_fullres(btf_path, tp_path, full_shape) -> tuple[float, float, float]`：
    ① **粗估**：在 `THUMB_SCALE=32` 的全片縮圖上跑 `estimate_shift`，乘 32 得 fullres 初值；
    ② **精修**：依初值在 3 個高 bin 密度區域各取 1024×1024 窗格，於 `downsample=4` 重跑 `estimate_shift`，乘 4，取 **中位數**（抗離群）。
    回傳 fullres px 的 `(dy, dx, spread)`，`spread` = 三窗估計值的最大差，作為可信度指標。
  - **為何必須兩階段**（初稿缺陷）：`THUMB_SCALE=32` + `upsample_factor=10` 的理論極限是 0.1 縮圖 px = **3.2 fullres px**；而傷害最大的區間正是 20-50 px 偏移，在縮圖上僅 0.6-1.6 px，SNR 太低。且 ds=32 時主導訊號是組織輪廓（適合粗對位），細節紋理已被抹平（不適合精修）。單階段估計會給出「看起來有數字但不可信」的結果 —— 比沒有更危險。
  - 驗證：`test_estimate_shift_fullres_recovers_30px`：對合成全片圖施加已知 30 px 位移，斷言回推誤差 < 2 px 且 `spread < 3`
  - 檔案：`backend/src/registration/align.py`
  - commit：`feat(registration): 兩階段 coarse-to-fine 位移估計`

- [ ] **P1-5c 座標軸與單位轉換（防呆）**
  - 預期行為：`shift_to_matrix(dy, dx) -> AffineAlignment`（P3-2 後生效；P1 階段先只寫轉換函式與測試）。明確處理兩個轉換：**① 軸序**：`phase_cross_correlation` 回傳 `(row, col)` = `(dy, dx)`，而 affine/座標一律 `(x, y)` = `(col, row)` → 必須交換；**② 單位**：估計值在 downsample 圖上，需乘回 downsample 倍率。
  - **為何獨立成一個任務**：row/col ↔ x/y 交換與 downsample 倍率是此類配準程式最常見的靜默錯誤來源（結果不會報錯，只會位移錯方向或錯 32 倍）。用一個顯式測試把它釘死。
  - 驗證：`test_shift_to_matrix_axis_order`：`shift_to_matrix(dy=10, dx=-5)` 套用到點 `(0,0)` → 得 `(-5, 10)`（x 位移 -5、y 位移 10）
  - 檔案：`backend/src/registration/align.py`
  - commit：`feat(registration): 位移→仿射矩陣的軸序與單位轉換`

- [ ] **P1-6 QC 疊圖產生器**
  - 預期行為：`qc.render_overlay_patches(btf_path, tp_path, out_dir, n=3, size=512, seed=0) -> list[Path]`：從有 bin 的區域隨機取 n 個 patch，用 `read_btf_crop(btf_path, x0, y0, w, h)` 取 H&E（**注意回傳是 3-tuple `(img, actual_x0, actual_y0)`**，實際原點可能因 tile 對齊而異動，散點座標必須以回傳的 `actual_x0/y0` 為基準而非請求值），散點疊上該窗內 bin 質心，300 DPI PNG 存 `{results}/qc/alignment/patch_{i}.png`。`seed` 固定以可重現。
  - 驗證：`test_render_overlay_patches_writes_n_png`（以 CRC 資料，`@pytest.mark.skipif` 資料不存在時跳過）
  - 檔案：`backend/src/registration/qc.py`
  - commit：`feat(registration): 對位 QC 疊圖`

- [ ] **P1-7 API：/api/registration/estimate 與 /qc_patches**
  - 預期行為：`GET /api/registration/estimate?roi_name=` → `{dy, dx, error, downsample}`；`POST /api/registration/qc_patches` 產圖後 `GET /api/registration/qc_images` 回 base64 清單。縮圖沿用 `DZITileServer` 的 raw-TIFF 縮圖（`THUMB_SCALE=32`），避免重建。
  - 驗證：`.venv/bin/python -m pytest backend/tests/test_11_registration.py -q`
  - 檔案：`backend/src/api/registration.py`、`backend/main.py`（掛 router）
  - commit：`feat(api): 配準估計與 QC 疊圖端點`

- [ ] **P1-8 強化 counter 的偏移警告**
  - 預期行為：`counter.py:104-111` 的 30% 警告觸發時，額外對該 ROI 呼叫 `estimate_shift` 並 log `建議偏移 dy=.. dx=..`。失敗時只 warning 不中斷（依 CLAUDE.md §11 容錯規範）。
  - 驗證：`.venv/bin/python -m pytest backend/tests/ -q` 全綠
  - 檔案：`backend/src/cellpose_counter/counter.py`
  - commit：`feat(count): 偏移警告附建議位移值`

---

## P2 — 區域選取做圖分析（預估 1 天）

- [ ] **P2-1 【紅燈】區域過濾測試**
  - 預期行為：`test_filter_by_region_bbox` / `test_filter_by_polygon`：10 個已知 centroid 的假 AnnData，斷言 bbox 與三角 polygon 各自留下的細胞數與 index 正確；邊界點（正好落在邊上）計入。
  - 驗證：`.venv/bin/python -m pytest backend/tests/test_12_region.py -q` → FAILED
  - 檔案：`backend/tests/test_12_region.py`
  - commit：`test(spatial): 區域過濾紅燈測試`

- [ ] **P2-2 【綠燈】實作 filter_by_region**
  - 預期行為：`spatial.py` 新增 `_filter_by_region(adata, region, polygon, pixel_size_um) -> np.ndarray`（布林遮罩）。**座標系必須是全片 fullres px**，取用順序：① `obs['centroid_x_fullres']/['centroid_y_fullres']`（P0-7 產出，全圖流程）→ ② `obs['centroid_x_px'] + roi.x`（ROI 流程，由 config 補原點）→ ③ `obsm['spatial'] / pixel_size_um + roi.x`（舊資料回退）。polygon 用 `matplotlib.path.Path(...).contains_points(pts, radius=1e-9)`。`region` 與 `polygon` 同時給時 polygon 優先。
  - **⚠️ 這是計畫初稿的錯誤修正**：初稿直接拿 `centroid_x_px` 比對 RegionSelector 的全片座標，但前者是**裁切局部**座標 → 框選會系統性偏移一整個裁切原點（可達數萬 px），選到的是完全錯誤的區域。三層回退確保三種資料來源都落在同一座標系。
  - 驗證：同上 → PASSED
  - 檔案：`backend/src/api/spatial.py`
  - commit：`feat(spatial): 區域/多邊形細胞過濾`

- [ ] **P2-3 GenePlotRequest 加入區域欄位**
  - 預期行為：`GenePlotRequest` 新增 `region: Optional[dict] = None`（key: x0,y0,x1,y1，fullres px）與 `polygon: Optional[list[list[float]]] = None`。過濾後細胞數 < 10 時回 `status:error`「選取區域細胞數不足（<10）」。
  - 驗證：`test_gene_plot_rejects_tiny_region`
  - 檔案：`backend/src/api/spatial.py`
  - commit：`feat(spatial): gene_plot 支援區域選取`

- [ ] **P2-4 區域統計端點**
  - 預期行為：`POST /api/spatial/region_stats` 回 `{n_cells, median_counts, median_genes, cluster_counts}`（`cluster_counts` 取 `obs['leiden']` 或 `obs['celltypist_label']`，皆無則空 dict）。讓使用者框選後先看組成再決定做圖。
  - 驗證：`test_region_stats_returns_cluster_breakdown`
  - 檔案：`backend/src/api/spatial.py`
  - commit：`feat(spatial): 區域統計端點`

- [ ] **P2-5 前端：提取 RegionSelector 元件**
  - 預期行為：把 `Stage0_ROI.tsx` 的畫框邏輯提取為 `components/RegionSelector.tsx`，props：`{ imageUrl, fullWidth, fullHeight, mode: 'bbox'|'polygon', onChange }`。Stage0 改用該元件（行為不變，純重構）。
  - 驗證：`cd frontend && npm run build`；手動確認 Stage 0 仍可畫框
  - 檔案：`frontend/src/components/RegionSelector.tsx`、`frontend/src/pages/Stage0_ROI.tsx`
  - commit：`refactor(ui): 提取 RegionSelector 共用元件`

- [ ] **P2-6 前端：SpatialExplorer 接區域選取**
  - 預期行為：`Stage35_SpatialExplorer.tsx` 嵌入 `RegionSelector`（bbox + polygon 切換），選取後顯示 `region_stats`，做圖請求帶上 region/polygon。i18n 補字串。
  - 驗證：`cd frontend && npm run build`
  - 檔案：`frontend/src/pages/Stage35_SpatialExplorer.tsx`、`frontend/src/i18n/translations.ts`、`frontend/src/api/client.ts`
  - commit：`feat(ui): 空間探索器支援框選區域分析`

---

## P3 — Registration 模組正式化（預估 2-3 天）

- [ ] **P3-1 【紅燈】AffineAlignment 測試**
  - 預期行為：`test_affine_identity_is_noop`、`test_affine_translation_applies`、`test_affine_roundtrip_dict`：單位矩陣不動座標；平移矩陣正確位移；`to_dict`→`from_dict` 還原一致。
  - 驗證：`.venv/bin/python -m pytest backend/tests/test_11_registration.py -q` → FAILED
  - commit：`test(registration): AffineAlignment 紅燈測試`

- [ ] **P3-2 【綠燈】AffineAlignment dataclass**
  - 預期行為：`matrix: list[list[float]]`（2×3）、`source: str`、`target: str`、`estimated_error: float | None`。方法 `apply(coords_xy: np.ndarray) -> np.ndarray`、`to_dict()`、`from_dict()`、`identity()`。
  - 驗證：同上 → PASSED
  - 檔案：`backend/src/registration/align.py`
  - commit：`feat(registration): AffineAlignment 資料結構`

- [ ] **P3-3 pipeline.yaml 新增 alignment 區塊**
  - 預期行為：`alignment: {enabled: false, matrix: [[1,0,0],[0,1,0]], source: spaceranger_fullres, target: raw_btf, estimated_error: null}`。`enabled: false` 時全流程行為與現在**完全一致**（零回歸風險）。
  - 驗證：`.venv/bin/python -m pytest backend/tests/test_01_infra.py -q`
  - 檔案：`config/pipeline.yaml`
  - commit：`feat(config): 新增 alignment 設定區塊`

- [ ] **P3-4 bin_attribution 套用變換**
  - 預期行為：`bin_attribution` 新增 `alignment: AffineAlignment | None = None` 參數；非 None 且 enabled 時，先對 `(col, row)` 套 affine 再減裁切原點。預設 None → 行為不變。
  - 驗證：`test_bin_attribution_with_translation_shifts_assignment`（平移 5 px 後 cell_id 改變符合預期）
  - 檔案：`backend/src/fullslide/pipeline.py`
  - commit：`feat(registration): bin attribution 套用 affine 變換`

- [ ] **P3-5 仿射估計（含旋轉縮放）**
  - 預期行為：`estimate_affine(ref, mov) -> AffineAlignment`：先 `estimate_shift` 取平移初值，再以 `skimage.transform.estimate_transform('affine', src, dst)` 對多個分塊（4×4 grid，每塊各自 phase correlation）的對應點擬合；離群塊（誤差 > 3×median）剔除後重擬。
  - 驗證：`test_estimate_affine_recovers_rotation`（合成 2° 旋轉，斷言回推矩陣誤差 < 1%）
  - 檔案：`backend/src/registration/align.py`
  - commit：`feat(registration): 仿射變換估計（平移+旋轉+縮放）`

- [ ] **P3-6 API：估計→寫入設定**
  - 預期行為：`POST /api/registration/apply` 執行 `estimate_affine`，把結果寫入 `state.json` 的 `alignment`（`save_state`），並回傳矩陣與 residual 供前端顯示。**不自動啟用**——回傳後由使用者按「套用」才設 `enabled: true`（避免壞估計靜默污染資料）。
  - 驗證：`test_registration_apply_writes_state`
  - 檔案：`backend/src/api/registration.py`
  - commit：`feat(api): 配準結果寫入 state`

- [ ] **P3-7 前端：對位面板**
  - 預期行為：Stage 0 或 Stage 2 加「對位檢查」面板：顯示 QC 疊圖、估計出的 dy/dx/residual、「套用此變換」按鈕（含二次確認）。
  - 驗證：`cd frontend && npm run build`
  - 檔案：`frontend/src/pages/Stage0_ROI.tsx`、`frontend/src/i18n/translations.ts`、`frontend/src/api/client.ts`
  - commit：`feat(ui): 對位檢查與套用面板`

---

## P4 — NDPI 讀取支援（預估 1-2 天）

> 決策：用 **tiffslide**（純 Python、建於既有 tifffile+zarr 之上），**不用** openslide-python——後者需 `brew install openslide` / Windows DLL，與跨機（ExFAT SSD + macOS ⇄ Windows）情境衝突。

- [ ] **P4-1 安裝 tiffslide**
  - 預期行為：`UV_LINK_MODE=copy uv add tiffslide`（先 `find . -name '._*' -delete`）。
  - 驗證：`.venv/bin/python -c "import tiffslide; print(tiffslide.__version__)"`
  - 檔案：`pyproject.toml`、`uv.lock`
  - commit：`chore: 新增 tiffslide 依賴（NDPI/SVS 支援）`

- [ ] **P4-2 【紅燈】SlideReader 協定測試**
  - 預期行為：`test_open_slide_dispatches_by_suffix`：`.btf`/`.tif` → `TiffSlideReader`；`.ndpi`/`.svs`/`.mrxs` → `TiffslideReader`；未知副檔名 raise `ValueError`。用 monkeypatch 避免需要真實檔案。
  - 驗證：`.venv/bin/python -m pytest backend/tests/test_13_slide_reader.py -q` → FAILED
  - 檔案：`backend/tests/test_13_slide_reader.py`
  - commit：`test(slide): SlideReader 分派紅燈測試`

- [ ] **P4-3 【綠燈】SlideReader 抽象層**
  - 預期行為：`utils/slide_reader.py` 定義 `Protocol SlideReader`：`dimensions -> (W, H)`、`level_count -> int`、`read_region(x, y, w, h, level=0) -> np.ndarray`、`mpp -> float | None`。`TiffSlideReader` 包裝現有 `read_btf_crop`/`read_strip_crop`；`TiffslideReader` 包裝 `tiffslide.TiffSlide.read_region`（轉 RGB ndarray、丟棄 alpha）。`open_slide(path)` 依副檔名分派。
  - 驗證：同上 → PASSED
  - 檔案：`backend/src/utils/slide_reader.py`
  - commit：`feat(slide): SlideReader 抽象層支援 NDPI/SVS`

- [ ] **P4-4 tile_server 改用 SlideReader**
  - 預期行為：`DZITileServer.__init__` 改由 `open_slide()` 取得 reader；`full_width/full_height` 取自 `reader.dimensions`；縮圖建構改用 `read_region(level=最接近 THUMB_SCALE 的層)`——NDPI/SVS 自帶金字塔，可省下自建縮圖的整段成本。BTF 無金字塔時回退現行邏輯。
  - 驗證：`.venv/bin/python -m pytest backend/tests/test_04_roi.py -q`；手動確認 Stage 0 檢視器仍正常
  - 檔案：`backend/src/roi/tile_server.py`
  - commit：`refactor(roi): tile server 改用 SlideReader`

- [ ] **P4-5 discovery 認得 NDPI/SVS**
  - 預期行為：`utils/discovery.py` 掃描的影像副檔名清單加入 `.ndpi`、`.svs`、`.mrxs`，label 標明格式。
  - 驗證：`.venv/bin/python -m pytest backend/tests/test_02_data.py -q`
  - 檔案：`backend/src/utils/discovery.py`
  - commit：`feat(data): 資料掃描支援 NDPI/SVS/MRXS`

- [ ] **P4-6 文件：NDPI 需先配準的限制**
  - 預期行為：README.md / README_zh.md 新增小節說明：NDPI 通常是**另外掃描**的高解析影像，與 Visium 玻片無共同座標系（Hamamatsu 40x ≈ 0.226 µm/px vs Visium 0.2737），**必須先跑 P3 配準**才能貼合 RNA。
  - 驗證：`git diff --stat README.md README_zh.md`
  - 檔案：`README.md`、`README_zh.md`
  - commit：`docs: NDPI 使用前提與座標系限制`

---

## P5 — 串流 tiled 分割（預估 2-3 天，最後做）

> 放最後：等 P0-P4 確立流程正確性後再動核心演算法，避免同時改變太多變因。

- [ ] **P5-1 【紅燈】tile_reader 介面測試**
  - 預期行為：`test_run_tiled_accepts_tile_reader`：傳入一個從記憶體陣列切片的 callable，斷言結果與直接傳 ndarray 的舊路徑遮罩**完全相同**（`np.array_equal`）。
  - 驗證：`.venv/bin/python -m pytest backend/tests/test_10_fullslide.py -q` → FAILED
  - commit：`test(segmentation): tile_reader 介面紅燈測試`

- [ ] **P5-2 【綠燈】run_tiled_mcseg_v2 接受 tile_reader**
  - 預期行為：簽章加 `tile_reader: Callable[[int,int,int,int], np.ndarray] | None = None` 與 `full_shape: tuple[int,int] | None = None`。`img` 與 `tile_reader` 二者必有其一。Phase 1 逐塊改由 reader 取圖。
  - 驗證：同上 → PASSED
  - 檔案：`backend/src/segmentation/cellpose_runner.py`
  - commit：`feat(segmentation): tiled 分割支援 tile 串流讀取`

- [ ] **P5-3 label map 落地為 memmap**
  - 預期行為：Phase 1 的拼接標籤圖改用 `np.lib.format.open_memmap(path, mode='w+', dtype=np.int32, shape=(H,W))`（21504×47104 int32 ≈ 4 GB 落磁碟而非 RAM）。路徑取 `{output_dir}/tmp_labels_{config_hash}.npy`（`config_hash` 同 P5-6 定義），**僅在最終遮罩成功寫出後才刪除**。
  - **為何用 config_hash 命名且延後刪除**（初稿缺陷）：初稿寫固定檔名 `tmp_labels.npy` 且「完成後刪除」，與 P5-6 的續跑機制直接矛盾 —— 中斷後重跑時 label map 已被刪或被不同參數的執行覆寫，續跑會拼出混合兩組參數的錯誤遮罩。以 hash 命名確保不同參數互不干擾，延後刪除確保續跑有東西可接。
  - 驗證：`test_memmap_labels_not_in_ram`（以 psutil 斷言峰值 RSS 增幅 < 全圖 int32 大小的 50%）
  - 檔案：`backend/src/segmentation/cellpose_runner.py`
  - commit：`perf(segmentation): label map 改用 memmap 落地`

- [ ] **P5-4 Phase 2 Voronoi 分塊處理**
  - 預期行為：全圖 Voronoi 改為帶 overlap 的分塊處理（block 4096 + margin 256，margin 內結果丟棄），避免一次載入全圖距離場。驗證接縫：合成雙塊測試斷言跨塊細胞 label 連續。
  - 驗證：`test_blocked_voronoi_no_seam`
  - 檔案：`backend/src/segmentation/cellpose_runner.py`
  - commit：`perf(segmentation): Voronoi 分塊處理避免全圖距離場`

- [ ] **P5-5 run_full 改走串流路徑**
  - 預期行為：`_run_full_segmentation` 改以 `open_slide()` + `tile_reader` 呼叫，移除 `np.array(arr)` 與 `max_load_gb` 檢查（改為僅在 `tile_reader` 不可用時才檢查）。
  - 驗證：`.venv/bin/python -m pytest backend/tests/ -q` 全綠
  - 檔案：`backend/src/api/segmentation.py`
  - commit：`feat(segmentation): 全圖分割改走串流路徑`

- [ ] **P5-6 續跑機制**
  - 預期行為：Phase 1 每完成一列 tile 寫 `{output_dir}/full_seg_progress.json`（`{done_tiles: [[ty,tx],...], config_hash}`）；重啟時 `config_hash` 相同則跳過已完成 tile。`config_hash` 為 `seg_cfg` + 裁切座標的 sha256 前 12 碼。
  - 驗證：`test_resume_skips_done_tiles`（模擬中斷後重跑，斷言 reader 未被呼叫於已完成 tile）
  - 檔案：`backend/src/segmentation/cellpose_runner.py`
  - commit：`feat(segmentation): 全圖分割中斷續跑`

---

## 風險與已知取捨

| 風險 | 影響 | 緩解 |
|------|------|------|
| P3 配準估計不準卻被套用 | RNA 靜默錯位，比不配準更糟 | `alignment.enabled` 預設 false；P3-6 要求使用者手動確認才啟用；一律回報 residual |
| P0-13 全片 `aggregate_cells` 記憶體 | 全片 2µm bins > 100 萬，`sc.read_10x_h5` 一次載入 | 沿用 CLI 已驗證路徑（實測可行）；若 OOM 則改 backed mode，列為後續 |
| P5-4 分塊 Voronoi 產生接縫 | 跨塊細胞被切成兩個 | margin 256 px（> `voronoi_distance` 最大值 9 的 28 倍）+ 專門的接縫測試 |
| P4 tiffslide 對 NDPI 支援深度未實測 | 可能仍需 openslide | P4-2/P4-3 完成後**先用真實 NDPI 檔驗證**再往下做；失敗則退回 openslide 並記錄於 sb |
| 測試 fixture 硬編碼 CRC 路徑 | Windows / 換樣本時測試失敗 | 新測試一律 `@pytest.mark.skipif(not path.exists())`，合成資料優先 |
| P0 與 P5 都改 `segmentation.py` | 合併衝突 | P5 排在最後，且 P0-8 已把參數移出硬編碼，P5 只需改讀取路徑 |
| **P0-4/P0-6 搬移 CLI 內部函式時無回歸網** | CLI 全片流程目前無端到端測試，搬移可能靜默破壞 | P0-4 驗證含 `--help`（僅擋 import/語法錯誤，不擋邏輯錯誤）。**建議在 P0-4 前先補一個 CLI 小樣本煙霧測試**：用 512×512 合成 BTF + 20 bins 走完 `--he-crop` 路徑，斷言產出 `mcseg_mask.npy` 與 `cells.h5ad`。列為 P0-3b，非阻塞但強烈建議 |

---

## 計畫覆核紀錄（2026-07-25，撰寫後自我 review）

初稿寫完後對關鍵假設做了實地驗證，修正 4 項缺陷。**驗證通過**的假設：`filtered_feature_bc_matrix.h5` 位於 `binned_002` 根目錄、`spatial/tissue_positions.parquet` 存在、`psutil 7.2.2` 已可用（P5-3 測試可行）、`main.py` 以 `include_router(prefix=...)` 逐一掛載（P1-7 作法正確）。

| # | 缺陷 | 嚴重度 | 修正 |
|---|------|--------|------|
| 1 | P2 直接用 `centroid_x_px`（裁切局部）比對全片框選座標 | **高** — 框選會偏移一整個裁切原點（可達數萬 px），選到完全錯誤的區域 | P0-7 增設 `centroid_*_fullres` 欄位；P2-2 改為三層座標系回退 |
| 2 | P1 單階段估計解析度不足（ds=32 極限 3.2 fullres px，而目標區間 20-50 px 僅佔 0.6-1.6 縮圖 px） | **高** — 會產出「有數字但不可信」的估計，比沒有更危險 | 新增 P1-5b 兩階段 coarse-to-fine（ds=32 粗估 → ds=4 三窗精修取中位數 + spread 可信度） |
| 3 | P5-3 固定檔名 `tmp_labels.npy` + 完成即刪，與 P5-6 續跑機制矛盾 | 中 — 續跑會拼出混合兩組參數的錯誤遮罩 | 改 `tmp_labels_{config_hash}.npy`，延後至最終遮罩寫出後才刪 |
| 4 | 未處理 row/col ↔ x/y 軸序與 downsample 單位轉換 | 中 — 靜默錯方向或錯 32 倍，不會報錯 | 新增 P1-5c 專責轉換函式 + 顯式軸序測試 |

未修正但已記錄的弱點：P0-4/P0-6 的 CLI 重構缺乏端到端回歸測試（見風險表最後一列，建議補 P0-3b）。

## 不在此計畫範圍

- 全片 Xenium Explorer 匯出效能優化（`cli/segment.py:334-394` 的 `regionprops` 逐細胞輪廓，全片數十萬細胞會很慢）
- `CLAUDE.md` 因 gitignore 而漂移的結構性解法
- Stage 3 分析在全片規模下的效能（Leiden/UMAP 於 >50 萬細胞）

---

## 執行方式

每個任務完成後 `git commit`（訊息如各任務所列）。建議分支：

```bash
git checkout -b feat/fullslide-and-registration
```

階段性驗證（每個 P 結束時）：

```bash
find . -name '._*' -delete
.venv/bin/python -m pytest backend/tests/ -q     # 全綠
cd frontend && npm run build                      # 通過
```
