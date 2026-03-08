# Visium HD Pipeline 2 — 開發進度總結（更新：2026-03-08 Session 2）

---

## 已完成 Stage 總覽

| Stage | 功能 | 狀態 |
|-------|------|------|
| Setup | 資料自動掃描（discovery.py）| ✅ |
| Stage 0 | ROI 裁切（BTF tile-by-tile）| ✅ |
| Stage 1 | Cellpose + Logic A + 互動預覽 | ✅ |
| Stage 2 | Zarr 建構（macOS `._*` 防護）| ✅ |
| Stage 2.5 | Proseg 條件測試（Top 3 縮圖 + 排序表）| ✅ |
| Stage 3 | Proseg 完整執行 | ✅ |
| Stage 4 | Scanpy QC + UMAP + Leiden | ✅ |
| Stage 5 | Xenium Explorer + Loupe Browser 匯出 | ✅ |

## 本 Session 完成（2026-03-08 Session 1）

- Stage 1 預覽圖互動座標（hover 十字準線 + badge，點擊自動填快速預覽座標）
- NumberInput / RoiNumCell 鍵盤輸入修復（localStr + isEditing ref）
- ROI 覆寫表格新增 Eosin BG 欄位（per-ROI eosin_bg_threshold）
- params 變更時自動清除過時快速預覽
- 後端 /run 競態條件修正（status 在 lock 內設定）
- eosin_bg_threshold fallback 改從 postprocessing 讀取（之前誤讀 preprocessing）
- 統一 JPEG 品質常數 _PREVIEW_JPEG_QUALITY = 85
- Git commit: 2bc2052

## 本 Session 完成（2026-03-08 Session 2）

### Stage 4 程式碼審查與改善

**前端（Stage4_Analysis.tsx）**：
- 合併 `ChartView` / `UMAPChartView` 為統一元件（新增 `fullWidthKeys` prop）
- 修復 `useEffect` dependency array（加入 `chartDataMap`）
- 新增 `applyLabelError` state + try/catch in `handleApplyLabels`
- 新增 `.catch()` 至 `handleAnnotateResChange`
- `handleRunQC` 加入 `setClusterMeta({})` 清除舊標籤狀態

**後端（api/analysis.py）**：
- `_annot_suggestions` type 修正：`dict[str, Any]`（補 `Any` import）

**後端（analysis/pipeline.py）**：
- `get_cluster_ids()` return type annotation 修正：`-> tuple[list[str], dict[str, str]]`

### 多 ROI 合併匯出 Bug 修復（Xenium Explorer）

**問題**：合併 ROI 模式下 Xenium 匯出失敗：
```
[WARNING] 刪除 3349 個無有效多邊形的細胞
[ERROR] Reindexing only valid with uniquely valued Index objects
```

**根本原因**：
1. `merge_all_rois()` 呼叫 `ad.concat()` 合併多個 ROI 的 h5ad，每個 ROI 的 obs_names 均從 0 開始，造成重複 index
2. `_load_polygons_and_table()` 呼叫 `adata[valid_cell_ids].copy()` → pandas `.loc[]` 在重複 index 上失敗
3. `_run_xenium()` 只讀取第一個 ROI 的 GeoJSON，其他 ROI 的細胞找不到多邊形而被刪除

**修復 1（`backend/src/analysis/pipeline.py`，`merge_all_rois()`）**：
```python
# 加前綴避免多 ROI 合併後 obs_names 重複
adata.obs_names = [f"{roi_name}__{name}" for name in adata.obs_names]
```

**修復 2（`backend/src/api/export.py`，`_run_xenium()`）**：
- 偵測合併模式：`merged_h5ad = output_dir_base / "umap_computed.h5ad"` + `len(rois) > 1`
- 合併模式下從所有 ROI 的 GeoJSON 各自讀取，並在每個 feature 的 `full_id` 加上 `{roi_name}__` 前綴
- 合併後儲存至 `combined_all_rois.json`，確保 full_id 與 obs_names 一致

**注意**：需重新執行 Stage 4 QC（合併模式）才能套用 obs_names 前綴修復，之後再執行 UMAP → Xenium 匯出。

---

# 原始詳細紀錄（Detailed Historical Notes）

## 1. Eosin 細胞質遮罩 (Cyto Mask) 邏輯重構與修復

## 1. Eosin 細胞質遮罩 (Cyto Mask) 邏輯重構與修復

- **問題**：原本使用的 `R - B` (Red 減 Blue) 通道相減法來判斷 Eosin 組織分佈，對於判定組織與空腔背景的準確率極低。且在正式分割 `run_segmentation_rois` 中，Watershed 結果錯誤覆寫了原本 Cellpose 辛苦算出的細胞核/膜成果，導致偵測出的細胞數量從數百顆掉到 77 顆以內。
- **解法**：
    1. **改用「亮度」判斷 (Brightness Method)**：對於影像中所有像素取 RGB 最大值 $max(R, G, B)$，大於 `(255 - Eosin_BG_Threshold)` 的視為純空腔背景，其餘保留為組織細胞質。
    2. **分離最終遮罩**：取消對 `final_masks` (即 `segmentation_masks.npy`) 的覆寫與破壞。Cellpose / Logic A 產出的 `final_masks` 僅負責傳遞高品質細胞邊界，而 Eosin 亮度判斷另外生成二元矩陣並獨立儲存為 `cyto_mask.npy`，專供下游 Proseg 定位使用。

## 2. 預覽介面 (Preview UI) 大幅擴充

- **快速 Patch 預覽整合**：在 512x512 小塊預覽中，統一以 Tab 切換顯示：
  - H&E + 綠色邊界
  - Macenko 前處理影像 (Cellpose 真實輸入)
  - 光流方向圖 (Flows dP)
  - Cyto 遮罩
- **完整分割結果預覽支援多圖**：在正式跑完全圖/全 ROI 的任務後，除了能在「完整分割結果預覽」查看 H&E 疊圖外，現在新增：
  - Eosin 細胞質背景圖 (cyto_mask)
  - 小尺寸光流方向圖 (Flows dP)
- **前後端快取機制**：更新 Backend GET `/preview` 端點，支援自動偵測與讀取 `cyto_mask.npy` 與 `flows_preview.jpg` 轉為 base64 傳回給前端切換。

## 3. 面板參數 UI 體驗優化

- **新增滑鼠懸停提示 (Tooltips)**：對所有的 `NumberInput` 與 `Toggle` 元件加上說明支援。
- **翻譯與定義**：為每一個 segmentation 參數補齊繁體中文的實戰意義與調整建議（例如 Batch Size, Flow Threshold, Eosin BG Threshold, Logic A 雙尺寸設計等）。這在實驗與參數調整的過程中，幫助釐清各數值的增減邏輯。

## 4. Stage 3 Proseg 分割越界問題修復

- **問題**：在 Proseg 分子指派階段，即使啟用了 Watershed 隔離與 Cyto 遮罩空間防護，輸出的細胞多邊形（Polygon）依然像預設設定一樣盲目擴張，甚至筆直橫跨切斷了相鄰的細胞核。
- **解法**：盤查 `backend/src/proseg/pipeline.py` 後，發現 **Cyto 物理邊界過濾 (cyto_constraint)** 階段有一個極其隱蔽但致命的 NumPy 邏輯 Bug。
  - **錯誤寫法**：`np.where((dilated_mask > 0) & (cyto_constraint == dilated_mask), ...)`。因為 `cyto_constraint` 是剛產生的 0/1 二元背景遮罩，而 `dilated_mask` 是細胞的唯一編號 (2, 3... N)。這導致除了細胞 ID 1 以外，幾乎所有擴張防護在這一步被瞬間強制歸零。
  - **影響**：所有落入這個過濾盲區的細胞核周邊 RNA，全部被錯誤標記為 Background/Unassigned (0)。系統失去核周範圍後，Proseg 失去對各細胞核原有領地的認知，只能退回最無腦的 Gaussian Voronoi 全局擴張，產生切豆腐般的生硬幾何邊界。
  - **修復**：已將條件修正為 `(cyto_constraint > 0)`。現在 Watershed 的「楚河漢界」與 Eosin 背景邊界能完美無缺地保留給每個獨立細胞 ID。Proseg 將會老老實實順著我們給的防護邊緣進行分子歸屬，不再發生越界搶奪細胞的問題。

## 5. MCMC 演算法的細胞過度吞併問題終結

- **問題**：修復 Cyto 約束後，發現細胞總數仍然從原本 Cellpose 抓出的 308 顆掉至約 211 顆，且依然有細胞輪廓橫切過其他細胞核的情況，導致許多小細胞在運算過程中被抹除。
- **解法**：深入研究 Proseg （基於 MCMC 概率模型）的行為機制後發覺：
  - **根本原因**：即便我們透過 Pipeline 強制配發了精準的細胞核初始領地，Proseg 預設仍帶有 `--nuclear-reassignment-prob 0.2`（或 0.01）的功能。這意味著在預設高達數百次的迭代抽樣（Iteration Sampling）中，核內的 RNA 被允許依據「機率」重新歸屬給附近勢力更大或範圍更廣的強勢細胞。經過疊加，大量位於原核內的基因直接被周遭大細胞掠奪，被吞併的小細胞也就從最終地圖上消失了。
  - **強硬鎖死**：將 CLI 呼叫參數強制設置為 `--nuclear-reassignment-prob 0` 與 `--prior-seg-reassignment-prob 0`。在數學與演算法上**百分之百封印**這項「強盜行為」，徹底禁止 Proseg 將任何已被確認在核內的 RNA 分配給別人。現在它必須老實地遵守 Cellpose 取出的核心基礎數量與核位，不再有細胞被無故剝奪消滅或輪廓被橫切的情形。

## 6. Dask 與 SpatialData 底層套件衝突修復 (Dependency Hell)

- **問題**：在最新版的 `dask (>= 2025.x)` 之中，其內建強制移除了舊有的 Legacy DataFrame 機制 (`dask-expr`)，這導致了我們管線中依賴舊版 DataFrame 介面的核心套件 `spatialdata` 讀取崩潰。這在執行 Proseg 條件掃描時會拋出 `The legacy implementation is no longer supported` 的致命錯誤並中斷執行。另外 Zarr `LocalStore` 也因改版而出現 `LocalStore object has no attribute path`。
- **解法**：
  1. **Dask 引擎強制降級與限定**：手動更新了環境的依賴，強制安裝 `"dask[dataframe]<2025.1.0"`，回退至具有舊世代相容性的穩定版本，並確保 `dask.config.set({'dataframe.query-planning': False})` 參數生效。
  2. **Zarr 讀取容錯修正**：修改 `backend/src/proseg/zarr_handler.py`，增強了 `sdata._zarr_store.store` 的相容性與容錯語法，能正確讀取 `.root` 底下目錄，避開了 Zarr 改版帶來的封裝路徑存取問題。

## 下一步/注意事項

由於移除了破壞性 Watershed 邏輯與修復了 Proseg Cyto 約束 Bug，目前的細胞數量、細胞核邊界防護與空間對齊度理應處於最佳狀態。您可以透過 UI 的「Proseg (Stage 3)」確認修復後的結果。如果對局部效果仍有要求，可回到 Segmentation 微調 `Eosin BG Threshold` (控制物理擴張牆壁) 或嘗試不同的 Cellpose dP 容忍度 `Flow Threshold`。

---

## 21. 縮圖高畫質升級 + dilation=0 支援 (2026-03-07)

### 21-1. Stage 2.5 縮圖升級至 800×800

**需求**：原 `preview.jpg`（400×400）在 Finder 資料夾瀏覽時細胞輪廓太小，難以比較條件差異。

**修改**（`backend/src/proseg/condition_tester.py`）：

| 項目 | 舊值 | 新值 |
|------|------|------|
| `DISPLAY_W / DISPLAY_H` | 400 × 400 | **800 × 800** |
| resize interpolation | 預設（bilinear）| `cv2.INTER_CUBIC` |
| 輪廓線寬 | 1px | **2px** |
| 文字大小 / 粗細 | 0.45 / 1 | **0.8 / 2** |

`preview_hd.jpg`（200px 原圖裁切 → 800px 4× zoom，供前端 modal）維持不變。

### 21-2. 前端 Modal 點擊放大（HD 縮圖）

**新增元件**（`frontend/src/pages/Stage2b_ConditionTest.tsx`）：
- `ThumbnailModal`：全螢幕黑底 overlay，顯示 `preview_hd.jpg`（點擊外部關閉）
- `ThumbnailCard`：新增 🔍 badge + `cursor-zoom-in`，點擊呼叫 `/conditions/thumbnail_hd/{idx}` 取得 base64 後開啟 modal

**新增 API**（`backend/src/api/conditions.py`）：
- `GET /conditions/thumbnail_hd/{condition_idx}` → 回傳 `preview_hd.jpg` base64

**新增 client 函數**（`frontend/src/api/client.ts`）：
- `getConditionThumbnailHd(idx: number)`

### 21-3. dilation=0 選項

**需求**：Stage 2.5 前端原無 `dilation=0` 選項，無法測試純轉錄本驅動的細胞形狀。

**修改**：
- `config/pipeline.yaml`：`grid.dilation` 新增 `0`（→ `[0, 10, 20, 30]`）；`quick_grid.dilation` 新增 `0`（→ `[0, 10, 20]`）
- `frontend/.../Stage2b_ConditionTest.tsx`：dilation options 從 `[5, 10, 20, 30]` 改為 `[0, 5, 10, 20, 30]`

### 21-4. Polygon Smoothing 評估

- 測試新版 `simplify(1.2)+buffer(2.5, res=32)` 與 STRtree overlap removal
- 視覺對比腳本：`scripts/temp/test_smooth_polygon.py`（封存時已刪除）
- **結論**：用戶決定保留原始 `simplify(0.4)+buffer(1.0)` 版本

### 21-5. 已知待辦

- 重新執行 Stage 2.5，生成新 800×800 `preview.jpg`（舊檔案仍是 400×400）
- 考慮以 2×2 tiling 重跑 Stage 3，改善 tile 邊界縫隙問題（需 8-12 GB RAM）

---

## 7. ExFAT .venv 損毀修復與 start.sh 強化 (2026-03-05)

- **問題**：在 ExFAT 磁碟上執行 `uv sync` 時，macOS 的 `._*` resource fork 檔案會先於目錄被複製，導致 `fonttools` 等套件的快取目錄無法建立，進而 `.venv` 損毀無法啟動。
- **解法**：將 `.venv` 改為 symlink 指向 APFS 磁碟：`.venv → ~/.venvs/visiumHD_pipeline_2`，並新增：
  - `export UV_CACHE_DIR="$HOME/.cache/uv"` 避免 ExFAT cache 污染
  - `start.sh` 自動偵測 symlink 是否存在並重建
  - 啟動前執行 `lsof -ti:8000,3000 | xargs kill -9` 清除舊行程，避免 Port 衝突

## 8. 全專案程式碼品質稽核修復 (2026-03-05)

本次對整個 backend 進行全面 audit，修復以下 6 項問題：

### 8-1. `cellpose_runner.py`：33 個 print → 結構化 logging

- 新增 `logger = logging.getLogger("pipeline.segmentation")`
- 所有 33 個 `print()` 替換為 `logger.info()` / `logger.warning()`
- Stage 1 分割進度現在會正確輸出至 frontend terminal，不再消失

### 8-2. `cellpose_runner.py`：TiffFile 資源洩漏修復

- `run_segmentation()` 中使用 `tifffile.TiffFile(input_path)` 後原缺少 close 呼叫
- 修復：在 Macenko calibration 完成後立即執行 `tif.close()`，防止大型 BTF 檔案造成 fd 洩漏

### 8-3. `pipeline.py`：Proseg GeoJSON 靜默失敗防護

- 原 L943 `features = []  # TODO: Handle map?` 在 Proseg 輸出格式不符預期時，靜默丟棄全部 polygon 而不報錯
- 修復：改為 `raise ValueError(...)` 並附上詳細的 Key 診斷訊息，讓問題立即可見

### 8-4. `api/proseg.py` + `api/segmentation.py`：TOCTOU race condition 修復

- 原 `/run` endpoint 中 "check status == running" 與 "add task" 之間存在 race window
- 修復：兩個 API 均加入 `_task_lock = asyncio.Lock()`，check-and-set 操作包入 `async with _task_lock:`（使用 asyncio.Lock 而非 threading.Lock，避免阻塞事件循環）

### 8-5. `config/pipeline.yaml`：補齊 `burnin_samples`

- `proseg.golden_params` 原缺少 `burnin_samples` 設定，僅依賴 Proseg 預設值 200
- 新增 `burnin_samples: 150`（samples=500 時 30% burnin，符合 MCMC 收斂最佳實踐）

## 9. Dask 自動升級導致 Zarr 建構失敗修復 (2026-03-05)

- **問題**：`uv` 自動將 dask 從 `2024.12.1` 升至 `2026.1.2`，該版本完全移除了 `query-planning: False` 相容性開關，`import dask.dataframe` 直接噴 `NotImplementedError: The legacy implementation is no longer supported`，導致「建構 Zarr（Stage 2）」無法啟動。
- **解法**：`uv add "dask[dataframe]<2025.1.0"` 鎖定至 `dask==2024.12.1`，並寫入 `pyproject.toml` 版本約束，防止日後再次自動升級。

## 目前產出物清單

| 路徑 | 說明 |
|------|------|
| `results/analysis/roi/text/segmentation_masks.npy` | Cellpose Logic-A 分割遮罩 |
| `results/analysis/roi/text/segmentation_masks.tif` | 同上 TIFF 格式（ZLIB 壓縮）|
| `results/analysis/roi/text/cyto_mask.npy` | Eosin 亮度法細胞質二元遮罩 |
| `results/analysis/roi/text/flows_preview.jpg` | Cellpose dP 光流預覽 |
| `results/analysis/roi/text/adata_002um.h5ad` | Visium HD 2µm AnnData |
| `results/analysis/roi/text/adata_008um.h5ad` | Visium HD 8µm AnnData |
| `proseg-output.zarr/` | Proseg 輸出（shapes + tables + points）|
| `results/zarr/text/` | 整合後 SpatialData Zarr |
| `results/figures/alignment_*.png` | 空間對齊品質圖（4張）|
| `results/figures/stage1_*.png` | Stage 1 分割驗證圖（2張）|
| `results/figures/clip_before_after.png` | Polygon clip 修正前後對比圖 |

## 10. Proseg 條件邊界互咬與預覽圖偏移修復 (2026-03-06)

- **問題一（視覺錯位大災難）**：即使有 Watershed 保護，細胞分割的「綠色多邊形」在 UI 顯示上仍然集體偏移，如同切斷了原本的細胞核地盤。
  - **根本原因 (Visual Offset Error)**：Pipeline 給予 Proseg 原生模型的 `coordinate_scale` 使用了全域寫死的常數 `0.2645833`，然而組織的 UI 繪圖與裁切所使用的是 `0.2737`。這一微小的倍率差距，在動輒好幾百 pixel 的長寬畫佈下，產生了高達十幾像素的「集體大平移」，導致完美的分割網格在圖上蓋歪了。
  - **解法**：在 `runner.py` 與 `condition_tester.py` 中，將送給 Proseg 的倍率由常數強制改為動態讀取該 ROI 的真實倍率 `rois[0].get("pixel_size_um")`，確保計算出的 GeoJSON 坐標在除回像素空間時 100% 吻合 H&E 紫色細胞核的核心。
- **問題二（梯田狀/鋸齒狀交界）**：即便座標對齊，依然觀察到細胞與細胞之間的交界完全沒有空隙，且呈現方塊梯田狀的「互咬」。
  - **根本原因**：在高密度的空間中，未分配的游離 RNA 掉在兩個細胞的極狹窄縫隙（1~2 micron）間。Proseg MCMC 的 `Space-Filling` 特性會逼迫模型盡可能「吞噬並瓜分」所有的公海 Voxel（即使距離很遠），導致細胞外框長出梯田狀的觸手。
  - **解法**：降低 MCMC 的爭奪空間！將 `config/pipeline.yaml` 的預設 `samples` 迭代次數從 200 大幅下調至 **50** (recorded_samples=20)，且提升 `compactness` (0.06 ~ 0.1)。結合啟用 `--enforce-connectivity`，讓 Proseg 在確認完核心領域與 Watershed 邊緣後提早結束運算，遏止了無限膨脹的細胞膜觸手。

## 11. 最終空腔幽靈細胞清除與幾何削切 (2026-03-06)

- **問題一（FastAPI 顯示崩潰 NaN Error）**：執行參數掃描時，後端爆出 `ValueError: Out of range float values are not JSON compliant: nan` 造成 500 Server Error。
  - **根本原因**：參數掃描計算面積變異係數或失敗時傳遞了 `float("nan")`。`nan` 在 Python 可行，但無法被 `json.dumps()` 轉換導致後端全毀。
  - **解法**：將 `condition_tester.py` 強塞預估值的 `float("nan")` 安全轉換為 `0.0`。
- **問題二（細胞質邊界被突圍，綠色網格瘋狂擴張至空白處）**：雖然修復了 Zarr 容錯讀取，讓 `eosin_cyto` 實體防護牆成功傳給管線濾掉外部游離點，但 Proseg 的原罪（Space-Filling Voronoi 特性）使得就算沒有基因，它依然會強制把算好的細胞多邊形**「無限往外長」**，填滿整個畫布，長出蜘蛛網般的觸手跑到完全沒有組織的背景上。
  - **終極解法（Shapely Geometric Clipping 削切法）**：在 `pipeline.py` (Stage 3) 放棄依賴 Proseg 去學習邊緣，而是在它運算完後**追加一道具毀滅性的後製裁剪工序** `_clip_polygons_with_cyto`。利用 OpenCV 將 `eosin_cyto` 轉化成 `shapely.geometry.Polygon`。當 Proseg 產出細胞 GeoJSON 後，將組織防護圖當成剪刀與核邊界做幾何交集 (`intersection`)。
  - **效果**：所有跨越雷池的過長細胞膜，立刻如同切蛋糕般完美被沿著 Eosin 確切輪廓截斷於組織內部；而所有完全掉在純背景上的「虛構幽靈細胞」則因交集為空而被整顆剔除。成功確保了最高的視覺與組織擬真效果！

## 12. 下游單細胞分析管線升級 (Scanpy) (2026-03-06)

在解決完 Segmentation (Proseg) 之後，我們成功進入了 Stage 4 (下游聚類分析)，將產出的切割地圖做生物學分群，期間修復並實裝了以下功能：

- **問題一（找不到分析輸入檔）**：原先的下游管線固定尋找舊版的單一 H5AD 路徑，但由於引進了「核心巨型圖像分塊 (Tiling)」機制，Proseg 會動態產出一個接合所有圖塊的 `proseg_cells.h5ad`，導致路徑脫鉤讀不到檔案。
  - **解法**：修改 `analysis/pipeline.py`，改讓管線智能擷取當下設定檔的 `roi_name`，並指向正確的分塊縫合最終檔案 `/results/analysis/roi/{roi_name}/proseg_cells.h5ad`。
- **功能擴增（單細胞分析儀表板化與 UMI 深度過濾）**：
  - 以往分析參數直接寫死在 `pipeline.yaml`，需要手動改檔重啟。
  - **實作進度**：
    1. 前端 UI (`Stage4_Analysis.tsx`) 加入完整的深色調「控制面板」。
    2. 後端 `/api/analysis/run` 新增 Pydantic 模型接收前端的參數設定。
    3. 加入 `Resolution (聚類解析度)`、`n_pcs (PCA 維度)`、`min_genes (最低基因數)`、`max_pct_mito (粒線體上限)`，並能自動即時寫回 YAML 設定檔。
    4. **實裝了額外的品質控制閘門 `min_counts` (最低 UMI 數)**：於底層預處理模組加入 `sc.pp.filter_cells(adata, min_counts)` 功能，避免僅有基因種類達標但讀值極低（深度不足）的背景微滴雜訊干擾聚類，大幅提高 Leiden 分群結果的純度！

## 13. Stage 4 三區塊 UI 重構與圖表升級 (2026-03-06)

### 13-1. 三區塊 UI 設計

原本的 Stage 4 分析頁面為單一整合區塊，本次完全重構為三個**依序解鎖**的獨立區塊：

| 區塊 | 標題 | 內容 |
|------|------|------|
| Block 1 | QC 前處理 | QC → normalize → HVG → PCA； 顯示小提琴圖、散佈圖、PCA Elbow |
| Block 2 | UMAP 解析 | KNN → UMAP → Leiden，支援多組解析度同時比較；個別圖縮小為 50% 寬，Grid 全寬 |
| Block 3 | 熱圖輸出 | 選擇解析度 → 同時產生 Heatmap + Dotplot |

- Block 2 在 Block 1 完成前鎖定 (`opacity-50`)
- Block 3 在 Block 2 完成前鎖定
- 每次重新執行自動清空舊圖、後端快取同步清除，確保即時更新

### 13-2. 後端 pipeline 新增三步驟函數

- **`run_qc_step(config)`**：QC metrics → violin/scatter/PCA elbow 三張圖 → 儲存 `qc_preprocessed.h5ad`
- **`run_umap_step(config, resolutions, n_pcs, n_neighbors, min_dist)`**：KNN → UMAP → Leiden（多解析度）→ 各解析度 PNG + Grid PNG → 儲存 `umap_computed.h5ad`
- **`run_heatmap_step(config, resolution, n_top_genes)`**：同時產生兩張圖，回傳 `dict[str, str]`：
  - **Heatmap**：`seaborn.clustermap`，顯示**全部 HVGs**，行（cluster）與列（基因）均有樹枝圖（`row_cluster` / `col_cluster`），cluster 數 < 2 時自動關閉對應樹枝圖
  - **Dotplot**：`sc.pl.dotplot`，每 cluster 取 `n_top_genes` 個 marker 基因（由 `n_top_genes` 控制），點大小 = 表達比例，顏色 = 平均表達量

### 13-3. 後端 API 新增 12 條路由

| 路由 | 說明 |
|------|------|
| `POST /run_qc` | 啟動 QC 步驟 |
| `GET /qc_status` | 查詢 QC 狀態（含磁碟恢復） |
| `GET /qc_images` | 取得三張 QC 圖（base64） |
| `POST /run_umap` | 啟動 UMAP |
| `GET /umap_status` | 查詢 UMAP 狀態 |
| `GET /umap_images` | 取得各解析度 UMAP 圖 + Grid |
| `POST /run_heatmap` | 啟動 Heatmap + Dotplot |
| `GET /heatmap_status` | 查詢狀態 |
| `GET /heatmap` | 取得 `{heatmap, dotplot}` dict |

所有 GET 端點均實作**磁碟 fallback**：後端 `--reload` 重啟後自動從磁碟恢復狀態，不需重跑。

### 13-4. Bug 修復記錄

| # | 問題 | 修復 |
|---|------|------|
| 1 | `n_components > min(n_samples, n_features)` | `run_pca()` 自動夾緊 `n_pcs` 上限 |
| 2 | `pca_variance_ratio() got unexpected keyword 'ax'` | 改為直接讀 `adata.uns["pca"]["variance_ratio"]` 用 matplotlib 繪圖 |
| 3 | 後端 `--reload` 重啟後圖表消失 | 磁碟 fallback 機制 |
| 4 | Heatmap PNG 空白 | 原 `sc.pl.matrixplot()` 未加 `return_fig=True`，改用 `seaborn.clustermap` |
| 5 | `clustermap` 單一 cluster 報錯 | cluster < 2 時 `row_cluster=False` |
| 6 | 圖片重跑不即時更新 | 點擊執行時清空 state，後端執行前清空快取 |

## 14. Stage 0 ROI 表單驗證 + Stage 1 細胞包細胞修復 (2026-03-06)

### 14-1. Stage 0 ROI 表單靜默失敗修復

- **問題**：按「新增 ROI」沒有反應——原因是 `handleAdd` 中 `!form.name || !form.tissue` 條件靜默 return，組織欄位未填時完全沒有任何提示。
- **修復**：
  - 新增 `formError` state
  - 個別欄位驗證並顯示對應錯誤訊息（「請填寫名稱」/ 「請填寫組織（CRC/LUAD）」）
  - 任何欄位變動時自動清除錯誤

### 14-2. 細胞包細胞（雙輪廓）問題修復

**問題根因 1：LOGIC_A 孤立像素 bug**

`_merge_masks_logic_a` 替換大細胞時只覆蓋大細胞範圍內的像素，若小細胞有任何像素延伸到大細胞邊界外，這些殘留像素仍保留在 `merged` 中形成孤立的外圈細輪廓。

修復：替換前先清除所有 `small_ids_in_region` 的整體像素：
```python
for sid in small_ids_in_region:
    merged[merged == sid] = 0
merged[region] = next_id
```

**問題根因 2：`merge_enclosed_cells` 判斷邏輯歷程**

| 版本 | 方法 | 問題 |
|------|------|------|
| v1 | 1px dilation，ring 中無 0 | Cellpose 核/胞質 1-2px 間隙導致幾乎都跳過 |
| v2 | 拓撲連通（外部背景集合） | 間隙通往外部 → 仍跳過 |
| v3（最終） | 覆蓋率 + 四象限雙重判斷 | 能容忍間隙，且防止誤合併相鄰細胞 |

**最終演算法（v3）**：
1. 對每個細胞做 `dilation_px=6` px dilation ring
2. 外層較大細胞在 ring 中佔比 ≥ `coverage_threshold=0.55` → 可能被包圍
3. 再確認該外層細胞出現在**所有 4 個象限**（上左、上右、下左、下右）→ 排除只是鄰居而非包圍的情況
4. 通過後合併

**新增位置**：
- `cellpose_runner.py`：`merge_enclosed_cells()` 函數
- `run_segmentation()`、`_run_single_roi_segmentation()`：正式分割後呼叫
- `segmentation.py`：`_run_preview_sync()` 預覽也同步呼叫，保持 Quick Preview 與正式結果一致

**config**：`postprocessing.enable_merge_enclosed: true`（可關閉）

---

## 15. Stage 2.5 條件測試改用中心子區域 + Stage 5 匯出 API 修復 (2026-03-06)

### 15-1. Stage 2.5：條件測試從 ROI 全圖 → 中心子區域

**問題**：`condition_tester.py` 中 `test_roi_um`（config 預設 1000µm）已定義但從未使用，每次條件測試都跑整張 ROI Zarr，若 ROI 較大（例如 1324×1324µm），計算量龐大。

**修復**：在 `ConditionTester` 新增 `_get_center_test_roi()` helper：
- 輕量讀取 `zarr.open()` 取得 `labels/cellpose_nuclei` 的 H×W shape（不載入全陣列）
- 將 `test_roi_um ÷ pixel_size_um` 換算為像素大小
- 裁切影像正中心，clamp 到合法邊界
- 失敗時 fallback 回傳 `None`（退回全圖模式，不中斷流程）

**呼叫位置**：`_run_proseg_minimal()` 在建立 `ProsegPipeline` 前呼叫，並傳入 `fixed_roi=center_roi`。

**效能提升**：以 1000µm² 子區域 vs 1324×1324µm 全圖，資料量縮小約 **17 倍**，條件測試速度大幅提升。

**調整方式**：`config/pipeline.yaml` 的 `condition_test.test_roi_um`（設 0 或負數則關閉，改用全圖）。

---

### 15-2. Stage 5：Xenium / Loupe 匯出 API 修復

**問題根因**：`export.py` 的 `_run_xenium` 與 `_run_loupe` 兩個函式把整個 `config` dict 當作 `zarr_path`/`poly_json_path` 傳入 exporter 的建構子：
```python
exporter = XeniumExporter(config)  # TypeError: argument should be str/PathLike, not 'dict'
exporter = LoupeExporter(config)   # 同上
```

**修復**：改從 `config` 正確解析各路徑：

| 參數 | 來源 |
|------|------|
| `zarr_path` | `paths.zarr_dir / {roi_name} / proseg_integrated.zarr` |
| `poly_json_path` | `paths.proseg_dir / {roi_name} / combined_proseg_results_qc.json` |
| `transcripts_csv_path` | `paths.proseg_dir / {roi_name} / combined_transcripts.csv` |
| `h5ad_path` | 自動依序搜尋 `clustered_final.h5ad → umap_computed.h5ad → qc_preprocessed.h5ad` |
| `out_dir` | `paths.export_dir / xenium`（或 `loupe`） |

路徑不存在時傳 `None`（exporter 有 graceful fallback），h5ad 找不到才拋出明確 `FileNotFoundError`。

---

## 16. Stage 5 Xenium 匯出深度修復：NoneType + 路徑找不到 + GeoJSON 合併 (2026-03-06)

### 問題描述

上一節（Section 15）修復了 `dict-as-path` 錯誤後，重新點擊「匯出 Xenium」仍然失敗：

```
錯誤：'NoneType' object has no attribute 'uns'
```

Terminal 顯示：
```
[ERROR] ERROR pipeline.api.export: Xenium 匯出失敗：'NoneType' object has no attribute 'uns'
```

### 根因分析

| # | 問題 | 說明 |
|---|------|------|
| 1 | `sd_table=None` 傳入 SpatialData | `_write_xenium_bundle` 不論 `sd_table` 是否為 `None` 都強制設定 `tables`/`table`，造成 `spatialdata_xenium_explorer.write()` 內部嘗試存取 `None.uns` |
| 2 | `combined_proseg_results_qc.json` 不存在 | `export.py` 沿用舊路徑 `results/proseg/{roi_name}/`，但實際資料存放在 `results/analysis/roi/test/`，因此永遠找不到多邊形 |
| 3 | `h5ad` 搜尋範圍不完整 | 僅搜尋 `results/analysis/`，未搜尋 `results/analysis/roi/test/`（`proseg_cells.h5ad` 所在地）|

### 修復內容

#### `backend/src/export/xenium_exporter.py`

**1. `_write_xenium_bundle`：`sd_table=None` 防護**

```python
# 修復前：不管 None 都強制設定
if "tables" in sd_init_sig.parameters:
    sdata_kwargs["tables"] = {"table": sd_table}
else:
    sdata_kwargs["table"] = sd_table

# 修復後：None 時不加入 sdata_kwargs
if sd_table is not None:
    if "tables" in sd_init_sig.parameters:
        sdata_kwargs["tables"] = {"table": sd_table}
    else:
        sdata_kwargs["table"] = sd_table
```

**2. 新增 `generate_combined_geojson()` 函數（模組層級）**

由於 Stage 3 (Proseg) 使用 tile-based 分塊，每個 tile 各自輸出 `proseg_results.json`（gzip GeoJSON，相對座標 µm），不存在預先合併的單一 GeoJSON。

此函數從 12 個 tile 自動合併並套用 offset：

```python
abs_µm = tile_rel_µm + (ix * tile_w_px * scale_um_px,
                         iy * tile_h_px * scale_um_px)
```

`full_id` 格式為 `tile_y{iy}_x{ix}_{cell_idx}`，與 `proseg_cells.h5ad` 的 `obs_names` 格式完全一致（3150 個細胞，12 tiles × ~260 cells/tile）。

#### `backend/src/api/export.py`

**`_run_xenium`：路徑邏輯全面修正**

| 項目 | 修復前 | 修復後 |
|------|--------|--------|
| 路徑基底 | `results/proseg/{roi_name}/` | `results/analysis/roi/{roi_name}/` |
| GeoJSON | 只讀取現有檔案 | 不存在時自動從 tile 重建並儲存 |
| h5ad 搜尋 | 僅 `results/analysis/` | 同時搜尋 `results/analysis/` 和 `results/analysis/roi/test/` 含 `proseg_cells.h5ad` |

### 測試驗證

```
Total features: 3150
first full_id: tile_y0_x0_0
last full_id: tile_y2_x3_322
last 2 coords (abs um): [[450.57, 376.83], [450.57, 378.83]]  ← 正確絕對座標
```

合併函數在系統 Python 3.13 環境下通過驗證，3150 個細胞多邊形正確對應至 obs_names。

---

## 17. Stage 5 Xenium 匯出持續除錯：AssertionError 座標系 (2026-03-06) ⚠️ 未完成

### 當前狀態

Section 16 的修復（NoneType + 路徑 + GeoJSON 合併）讓流程推進至 `spatialdata_xenium_explorer.write()` 內部，但仍有新錯誤。

### 錯誤訊息

```
AssertionError
  File "spatialdata_xenium_explorer/converter.py", line 116:
    df = utils.to_intrinsic(sdata, df, image_key)
  File "spatialdata_xenium_explorer/utils.py", line 98:
    return sdata.transform_element_to_coordinate_system(element, cs)
  File "spatialdata/transformations/operations.py":
    assert isinstance(source_coordinate_system, str)
```

### 相關 WARNING（轉錄點部分）

```
[WARNING] Trying to get an element key of `sdata.points`, but it contains multiple values and no key was provided. It will not be saved to the xenium explorer.
```

### 根因分析（進行中）

| 觀察 | 推斷 |
|------|------|
| 431 個細胞多邊形已成功寫入 | shapes + table 部分正常 |
| `sdata.points` 含多個值（WARNING） | 但我們只傳入 images/shapes/tables，points 應為空 |
| `to_intrinsic(sdata, None, image_key)` 失敗 | `get_element` 因 multiple values 回傳 `None`，轉換 `None` 時找不到 source coordinate system |

**待確認假設**：

1. `sdata.points` 為何非空？可能來源：
   - zarr 的 `points/transcripts` 被 spatialdata 自動讀入？
   - `TableModel.parse` 的某版本在 table 含 spatialdata_attrs 時自動建立 points？
2. `image_key=None`（`_load_image` 靜默失敗）是否導致 `to_intrinsic` 接收 None 作為 coordinate system？

### 本次新增修復（已生效，但未完全解決）

- `xenium_exporter.py` 所有 `parse()` 呼叫加入 `transformations={'global': Identity()}`（座標系宣告）
- 小型測試確認此修法本身正確，`to_intrinsic` 在 element 有明確 transformation 時可正常使用

### 下次繼續計畫

**步驟 1**：在 debug 腳本加 `print(dict(sdata.points))` 確認 sdata.points 是否真的有內容

**步驟 2**：確認 `_load_image` 是否靜默失敗（`tissue_hires_image` key 應存在於 zarr，但 log 中無「影像載入完成」）

```python
# 可能修法 A：在 write() 前加印
logger.info(f"sdata elements: images={list(sdata.images.keys())}, shapes={list(sdata.shapes.keys())}, points={list(sdata.points.keys())}")
```

**步驟 3a（若 points 非空）**：確認 spatialdata 版本，查是否有 `points` 自動讀入行為

**步驟 3b（若 image_key=None）**：修復 `_load_image` 使其正確讀取 dask array，或改用 `sd.read_zarr` + element selection

**步驟 4**：最終備選方案 — 傳入 `mode="cbom"` 跳過 transcripts 寫入（`t` flag），配合確認版本 API

### 已安裝版本

```
.venv/lib/python3.12/site-packages/spatialdata_xenium_explorer/
.venv/lib/python3.12/site-packages/spatialdata/
```

---

## 20. Stage 5 Xenium Explorer 座標對齊根本修復 + Tile 邊界去重 (2026-03-06)

### 20-1. Xenium Explorer 雙重縮放 Bug（根本原因）

**問題**：Xenium Explorer 開啟後，細胞遮罩（cell mask）與 H&E 背景嚴重錯位。

**根本原因**：`xenium_exporter.py` 在呼叫 `spatialdata_xenium_explorer.write()` 前，手動將座標從 µm 換算為 px（÷XENIUM_UM_PX），但 SpatialData 內部已自動做此轉換，導致 **雙重縮放**（1.288× 放大）。

**修復三處（`backend/src/export/xenium_exporter.py`）**：

| 位置 | 修復前 | 修復後 |
|------|--------|--------|
| `_load_image()` | `ndimage.zoom(×1.288)` + `Identity()` | 直接用原圖 + `Scale([pixel_size_um, pixel_size_um])` |
| `_load_polygons_and_table()` | `shapely_scale(÷XENIUM_UM_PX)` 換算至 px | 直接傳 µm 多邊形，不手動換算 |
| `_load_transcripts()` | `tx_dd["x"] /= XENIUM_UM_PX` | 移除手動換算，直接傳 µm |

---

### 20-2. Tile 邊界重複細胞 → GeoJSON 邊界剪裁

**問題**：各 tile 因 nucleus mask 有 200px padding，Proseg 在 padding 區偵測到相鄰 tile 的細胞，產生重複多邊形。舊的 combined JSON（3150 features）包含大量重複（實際 tile 總計 1878 cells）。

**修復**（`generate_combined_geojson` in `xenium_exporter.py`）：
- 每個 feature 計算重心 global µm
- 只保留重心落在該 tile 名義邊界 `[ix×tile_w, (ix+1)×tile_w)` 內的細胞
- 結果：1844 features，無重複，tile 間距 < 2µm

---

### 20-3. Tile Transcript Overlap（根本修復邊界細胞）

**問題**：轉錄本嚴格切在 tile 邊界，邊界細胞只有半側轉錄本 → Proseg 分割品質差 → 邊界出現空白帶。

**修復三層**：

| 層 | 檔案 | 修改 |
|----|------|------|
| 轉錄本 filter | `pipeline.py` | transcript 範圍延伸 `padding` px 進相鄰 tile（overlap = padding = 200px） |
| roi_offset 固定 | `pipeline.py` | 不從 transcript min 重算（避免雙重 padding），直接固定為 `max(0, roi_x - padding)` |
| h5ad 去重 | `runner.py merge_tiles` | 按 obsm["spatial"] 重心裁剪到名義邊界，與 GeoJSON 去重邏輯一致 |

**注意**：此修復需刪除舊 tile 結果並重跑 Stage 3 才能生效（Smart Resume 會跳過已有 tile）。

---

### 20-4. 診斷工具（已清除）

本次工作期間在 `scripts/temp/` 建立了以下診斷腳本（封存時已刪除）：
- `debug_geojson_overlay.py`：生成 GeoJSON 多邊形 vs H&E 疊圖、重心比較圖、tile 邊界色碼圖

**診斷結果**：
- GeoJSON 多邊形與 H&E 組織完美對齊（分割正確）
- 重心偏移由 -71.6px 縮至 -15.2px（舊 combined JSON 過時所致，換新 JSON 後改善）

---

## 19. Stage 4 H&E 疊圖視覺化 + QC 參數修正 (2026-03-06)

### 19-1. H&E 疊圖 QC 前後比較（新功能）

**需求**：在 Stage 4 QC 後可視覺化細胞重心疊在 H&E 上，並觀察 QC 刪除了哪些細胞。

**實作**（commit `80cdca9`）：

| 產出檔案 | DPI | 說明 |
|---------|-----|------|
| `overlay_pre_qc.png` | 150 | QC 前全部細胞（青色），前端預覽用 |
| `overlay_pre_qc_hd.png` | 300 | QC 前 HD 存檔，供下載 |
| `overlay_post_qc.png` | 150 | 保留（綠）+ 刪除（紅）比較，前端預覽用 |
| `overlay_post_qc_hd.png` | 300 | QC 後 HD 存檔，供下載 |

**技術細節**：
- 座標系：`obsm["spatial"]`（µm）÷ `pixel_size_um`（0.2737）→ HE 像素
- 新函數 `_generate_overlay_images()` 於 `backend/src/analysis/pipeline.py`，在 `run_qc_step()` 過濾前後各記錄 obs_names，自動比對
- 前端 Stage 4 新增 2 個 tab（H&E 疊圖 QC 前/後）+ HD 下載連結
- 新 API：`GET /analysis/overlay_hd/{pre_qc|post_qc}` → FileResponse

### 19-2. QC 預設參數修正（Visium HD 適配）

**問題**：預設 `min_counts=100`、`max_pct_mito=20%` 對 Visium HD 過於嚴格，導致 75%+ 細胞被刪除。

**根本原因**：
- Visium HD 每顆 proseg 細胞 UMI 中位數僅 57（非 scRNA-seq 的 1000+）
- 粒線體 % 中位數 34.5%（probe-based 空間資料本來就高）

**修正**（`config/pipeline.yaml` + `Stage4_Analysis.tsx`）：

| 參數 | 舊值 | 新值 | 原因 |
|------|------|------|------|
| `min_counts` | 100 | 5 | 中位數僅 57，舊值過嚴 |
| `min_genes` | 20 | 10 | 空間資料基因稀疏 |
| `max_pct_mito` | 20% | 80% | 空間資料 mt% 本來就高 |

### 19-3. Violin Plot Y 軸修正

**問題**：`max_genes=8000` 門檻線撐爆 y 軸，KDE 延伸到負值（-400）。

**修正**：`pipeline.py` violin 繪圖後加 `ax.set_ylim(bottom=0, top=data_99*1.5)`，下限鎖 0，上限夾至 99th percentile × 1.5。

---

## 18. Stage 2.5 縮圖升級 + 空腔幽靈細胞三個根本問題定位 (2026-03-06)

### 18-1. Stage 2.5 條件測試縮圖升級為 400×400

**需求**：原縮圖尺寸 686×398（scaled down）太小，難以觀察細胞輪廓細節。

**修改**（`condition_tester.py`，commit `090104e`）：

| 檔案 | 尺寸 | 方式 |
|------|------|------|
| `preview.jpg` | 400×400 | 全 ROI 等比縮放至 400×400 |
| `preview_hd.jpg` | 400×400 | 原始解析度正中心裁切 400×400 px，2px 綠色輪廓 |

`preview_hd.jpg` 生成邏輯：
1. 在原始 HE 影像上直接繪製 2px 綠色多邊形（不縮放）
2. 取影像正中心 400×400 crop（`y0 = h//2 - 200, x0 = w//2 - 200`）
3. 尺寸不足時以黑色 canvas padding
4. 右下角疊加條件標籤文字（`FONT_HERSHEY_SIMPLEX, 0.7`）
5. JPEG quality 95 輸出

---

### 18-2. 空腔區域幽靈細胞三個根本問題（已定位，尚未修復）

**問題現象**：`condition_comparison.png` 右下角白色空腔區域仍顯示綠色細胞輪廓，即使 eosin 遮罩理論上應覆蓋該區域。

**調查方法**：
- 讀取 `cyto_mask.npy`（1210×1491，87% 組織）並疊加 HE 影像（`he_cyto_overlay.png`）
- 確認 `eosin_cyto` zarr label 與 `cyto_mask.npy` 一致，遮罩本身正確

**定位出的三個 Bug（均在 `backend/src/proseg/pipeline.py`）**：

#### Bug 1：`_clip_polygons_with_cyto` 使用 `RETR_EXTERNAL`，無法感知孔洞

| 項目 | 說明 |
|------|------|
| 問題 | `cv2.findContours(..., cv2.RETR_EXTERNAL)` 只取最外層輪廓，vessel lumen（空腔）是「洞」而非「外輪廓」，完全被忽略 |
| 錯誤邏輯 | 建立 `cyto_union`（組織外輪廓） → 與細胞多邊形取 `intersection` → 空腔內細胞因在組織外輪廓「內部」而通過 |
| 正確修法 | 找背景（black=0）輪廓建立 `bg_union` → 用 `cell_poly.difference(bg_union)` 削去背景區域 → 加 1px erosion 避免切除 1-2px 邊緣雜訊 |

#### Bug 2：Smart Resume 略過削切步驟

| 項目 | 說明 |
|------|------|
| 問題 | `run_proseg()` 在所有輸出檔案已存在時提前 return（L657），跳過 L712-714 的 `_clip_polygons_with_cyto` 呼叫 |
| 正確修法 | 在 Smart Resume 早返回前補充呼叫削切函數 |

#### Bug 3：`proseg_results.json` 為 gzip 壓縮格式

| 項目 | 說明 |
|------|------|
| 問題 | 以 `open(..., 'r', encoding='utf-8')` 讀取 → `UnicodeDecodeError: 0x8b at position 1` |
| 偵測 | magic bytes `\x1f\x8b` 代表 gzip |
| 正確修法 | 偵測 magic bytes → 以 `gzip.open(..., 'rt', encoding='utf-8')` 讀取 → 寫出時同樣以 gzip 寫回 |

**最終決定**：三個 bug 均已確認，但因用戶決定以 `git restore` 返回先前版本，**所有 `pipeline.py` 修改均已還原**。`condition_tester.py` 的縮圖升級保留（commit `090104e`）。

**後續**：若需重新修復，三個問題的修法均已明確記錄於此，可直接套用。

---

## Section 21: Xenium Explorer 細胞輪廓錯位根因分析與修復（2026-03-07）

### 21-1. 問題描述
開啟 Xenium Explorer bundle 後，細胞輪廓位置與 H&E 不對齊（與 Section 20-2/20-3 同類問題再次出現）。

### 21-2. 根本原因：proseg_cells.h5ad 與 combined_proseg_results_qc.json 版本不一致

| 檔案 | 來自的 Stage 3 run |
|------|------------------|
| `combined_proseg_results_qc.json` | 08:14 最舊一次 tile run（tile_y0_x0 含 125 cells）|
| `proseg_cells.h5ad` | 08:43 中間一次 tile run（tile_y0_x0 含 ~154 cells）|
| 當前 tile files | 最新 run（tile_y0_x0 含 243 cells，IDs 0–242）|

三個版本錯開 → cell ID 比對率僅 35%（592/1693）→ 大量細胞沒有對應多邊形 → 位置異常。

### 21-3. 結構性修復

**一次性修復**：重跑 `merge_tiles` 與 `generate_combined_geojson`，強制從當前 tile files 同步兩個檔案。

**根本修復**（`backend/src/proseg/runner.py`）：  
在 `run_tiled_proseg` 的 `merge_tiles()` 呼叫後立即同步重建 `combined_proseg_results_qc.json`，確保兩者永遠來自同一次 tile run：

```python
# merge_tiles 完成後緊接著同步 GeoJSON
from backend.src.export.xenium_exporter import generate_combined_geojson
geojson = generate_combined_geojson(tile_proseg_dir=output_dir, ...)
with open(geojson_path, "w") as f:
    json.dump(geojson, f)
```

修復後比對率：98.9%（1480/1497 cells）。

---

## Section 22: 多 ROI 合併分析功能實作（2026-03-07）

### 22-1. 背景
當 `pipeline.yaml` 定義 2-3 個 ROI（來自同一 H&E + Visium），Stage 2/3 各自獨立處理。
新增合併功能讓 Stage 4 可以整合所有 ROI 一起做 QC + UMAP + Leiden。

### 22-2. 設計決策
- 同一 H&E + Visium → **無需 batch correction**
- 座標還原：各 ROI local µm + (roi.x × pixel_size_um, roi.y × pixel_size_um) → 全局 µm
- `obs["roi"]` 欄位記錄細胞來源，可在 UMAP 上著色

### 22-3. 實作位置

| 檔案 | 變更 |
|------|------|
| `config/pipeline.yaml` | 新增 `analysis.merge_rois: false`（預設關閉）|
| `backend/src/analysis/pipeline.py` | 新增 `merge_all_rois(config)` 函數；修改 `run_qc_step` 支援合併模式 |
| `backend/src/proseg/runner.py` | Stage 3 所有 ROI 完成後自動呼叫 `merge_all_rois`（若啟用）|

### 22-4. 使用方式
```yaml
# pipeline.yaml
analysis:
  merge_rois: true
```
Stage 3 完成後自動生成 `results/analysis/merged_all_rois.h5ad`；Stage 4 QC 自動從此檔載入。

---

## Section 23: Stage 4 UI 多 ROI 選擇器 + Heatmap 改善（2026-03-07）

### 23-1. Stage 4 來源選擇 UI

**需求**：Stage 4 讓使用者在執行 QC 前選擇：分析單一 ROI 還是合併所有 ROI。

**後端變更**（`backend/src/api/analysis.py`）：
- `QCParams` 新增 `merge_rois: Optional[bool]` 與 `roi_name: Optional[str]`
- 新增 `GET /analysis/available_rois` endpoint：列出所有有 `proseg_cells.h5ad` 的 ROI
- `run_qc` endpoint：in-memory 覆寫 config（不寫回 YAML），`roi_name` 指定時將該 ROI 移至列表第一位

**前端變更**（`frontend/src/pages/Stage4_Analysis.tsx`）：
- 新增 `analysisMode: 'single' | 'merge'` 狀態
- Radio group 選擇模式；單一模式下顯示 ROI 下拉選單
- 僅有多個 ROI 時才顯示選擇器
- `handleRunQC` 依模式傳入 `merge_rois` / `roi_name`

**前端 API**（`frontend/src/api/client.ts`）：
- 新增 `getAvailableRois()`

### 23-2. Heatmap 色域改善（Z-score + RdBu_r）

**問題**：min-max scaling 使每個基因的最高表達 cluster 都達到顏色上限（亮黃），導致所有區塊飽和。

**修復**（`backend/src/analysis/pipeline.py`，`run_heatmap_step`）：

| 項目 | 舊值 | 新值 |
|------|------|------|
| 縮放方式 | min-max | Z-score + clip(-2.5, 2.5) |
| 色圖 | `viridis` | `RdBu_r` |
| 色域固定 | 無 | `vmin=-2.5, vmax=2.5` |
| 色條標籤 | Scaled mean expr. | Z-score (mean expr.) |

Z-score：藍 = 低於平均，白 = 平均，紅 = 高於平均；個別 outlier 不再撐爆色域。

### 23-3. Heatmap 基因數控制 + 圖形比例修復

**問題 1**：前端 `n_top_genes` 只控制 dotplot，heatmap 永遠使用全部 HVG（最多 2000 基因），修改無效。

**問題 2**：大量基因時寬度封頂 60 吋、高度只跟 cluster 數相關，形成極扁熱圖。

**修復**：

| 層次 | 變更 |
|------|------|
| `HeatmapParams`（api/analysis.py）| 新增 `n_heatmap_genes: int = 50`（獨立控制熱圖基因數）|
| `run_heatmap_step`（pipeline.py）| 接受 `n_heatmap_genes`，從 HVG 中取方差最高的前 N 個基因 |
| 大小公式 | 寬度 `每基因 0.25 吋`，上限 80 吋；高度 `max(4, clusters×0.6+2, 寬度/4)`，確保高:寬 ≥ 1:4 |
| 基因標籤閾值 | ≤80 個才顯示（原 ≤150）|
| 前端 UI | 新增「熱圖基因數 (n_heatmap_genes)」slider；原「n_top_genes」改名為「Dotplot 每群基因數」|

**選基因邏輯**：計算 HVG 在全部細胞的方差（稀疏矩陣用 E[X²]-E[X]²），取 top N，確保最具鑑別力的基因優先顯示。
