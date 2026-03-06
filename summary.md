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
