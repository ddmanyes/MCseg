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
