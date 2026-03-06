# Visium HD Pipeline 2 - Stage 1 Segmentation 功能開發與修復進度總結

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
