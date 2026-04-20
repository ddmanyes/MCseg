# IMPLEMENTATION_PLAN.md
# MSseg Skill 實作計畫

## 目標摘要

建立兩個可移植 AI Skill（`msseg-segment`、`msseg-analyze`）的支援腳本層，
使 Skill 文件中的所有指令可以實際執行。

**核心發現（研究後）：**
- `run_mcseg_v2(img, cfg)` — 接收 numpy 陣列，Skill 需先讀取 BTF ROI crop
- `run_tiled_mcseg_v2(img, cfg, tile_size, overlap)` — 全圖 tiled 分割已存在
- `count_rna_per_cell(adata_path, mask_path, roi_x_px, roi_y_px)` — 計數已存在
- `XeniumExporter` — 設計給 Proseg，需建立 MCseg 適配器
- Skill 文件中的函數呼叫介面需修正以符合實際 API

---

## 文件架構圖

```
MSseg/
├── skills/
│   ├── msseg-segment.md        ← 修正 Step 3 函數呼叫介面
│   └── msseg-analyze.md        ← 修正 Step 7 export 呼叫介面
├── scripts/
│   ├── roi_sampler.py          ← 新建：組織感知 ROI 隨機抽樣
│   ├── seg_quality.py          ← 新建：NUC + MCseg 雙模式分割包裝器
│   ├── qc_metrics.py           ← 新建：FTC / NED / Co-exp 計算
│   ├── write_handoff.py        ← 新建：產出 handoff_report.json
│   ├── build_full_adata.py     ← 新建：全圖 AnnData 組裝（無預設 ROI）
│   ├── run_analysis.py         ← 新建：QC→normalize→PCA→UMAP→Leiden 一鍵執行
│   └── export_mcseg.py         ← 新建：MCseg → Xenium zarr + h5ad 匯出
└── backend/src/export/
    └── xenium_exporter.py      ← 新增 export_from_mcseg() 適配函數
```

---

## 任務列表

### Phase A：msseg-segment 支援腳本

---

#### A1｜建立 `scripts/roi_sampler.py`

**預期行為：**
讀取 `tissue_positions.parquet`，隨機抽取 3–5 個組織覆蓋率 ≥ 60% 的 ROI，
寫入 `results/qc_rois.json`。若 40 次嘗試後仍不足 2 個，降閾值至 40% 並警告。

**相關檔案：** `scripts/roi_sampler.py`（新建）

- [ ] **Red**：建立測試 `backend/tests/test_skill_scripts.py`，
  測試 `roi_sampler` 在假造的 `tissue_positions.parquet`（200 個 in_tissue bins）下
  能產出 ≥1 個 ROI，且每個 ROI 的 `coverage` 欄位存在。
  ```bash
  uv run pytest backend/tests/test_skill_scripts.py::test_roi_sampler_produces_output -v
  # 預期：FAILED（腳本不存在）
  ```

- [ ] **Green**：建立 `scripts/roi_sampler.py`（從 `skills/msseg-segment.md` Step 2 提取程式碼），
  加入以下改動：
  - 降採樣格子大小改為 `stride = max(8, ROI_PX // 64)`（自適應加速）
  - 覆蓋率計算使用 vectorized pandas，不用雙層迴圈
  - 輸出 JSON 包含 `timestamp`、`threshold_used`、`rois` 三個欄位
  ```bash
  uv run pytest backend/tests/test_skill_scripts.py::test_roi_sampler_produces_output -v
  # 預期：PASSED
  ```

- [ ] **Refactor**：確認無硬編碼路徑，全部從 `config/pipeline.yaml` 讀取。
  ```bash
  grep -n "hardcode\|/Volumes\|C:\\\\" scripts/roi_sampler.py
  # 預期：無輸出
  ```

- [ ] **Commit**：`git commit -m "feat(scripts): add tissue-aware ROI sampler"`

---

#### A2｜建立 `scripts/seg_quality.py`

**預期行為：**
讀取 `results/qc_rois.json`，對每個 ROI：
1. 用 `tifffile` tile-based 裁切 BTF 影像（禁止全圖載入）
2. 以 NUC 模式（dia=15, clahe=0, voronoi=0）呼叫 `run_mcseg_v2(img, cfg)`
3. 以完整 MCseg 模式呼叫 `run_mcseg_v2(img, cfg)`
4. 儲存 `results/qc/{roi_name}_nuc.npy` 和 `results/qc/{roi_name}_mcseg.npy`

**相關檔案：** `scripts/seg_quality.py`（新建）

- [ ] **Red**：在 `test_skill_scripts.py` 加入 `test_seg_quality_creates_masks`，
  使用 300×300 白色假影像，驗證產出兩個 `.npy` 檔。
  ```bash
  uv run pytest backend/tests/test_skill_scripts.py::test_seg_quality_creates_masks -v
  # 預期：FAILED
  ```

- [ ] **Green**：建立 `scripts/seg_quality.py`：
  ```python
  # BTF 裁切（tile-based，禁止全圖載入）
  import tifffile
  with tifffile.TiffFile(he_path) as tif:
      page = tif.pages[0]
      crop = page.asarray(out='memmap')[y0:y0+h, x0:x0+w]

  # NUC cfg
  nuc_cfg = {"dia_small":15,"dia_mid":15,"dia_large":15,
             "use_hematoxylin":False,"use_cpsam":False,
             "voronoi_distance":0,"clahe_clip_limit":0,
             "min_size":20,"max_size":6000,"use_gpu":False,
             "batch_size":4,"flow_threshold":0.4,"cellprob_threshold":-2.0}

  # MCseg cfg 從 pipeline.yaml 讀取
  ```
  ```bash
  uv run pytest backend/tests/test_skill_scripts.py::test_seg_quality_creates_masks -v
  # 預期：PASSED
  ```

- [ ] **Refactor**：加入每個 ROI 的 try/except，失敗時 logger.warning 並 continue（參考 CLAUDE.md §11）。

- [ ] **Commit**：`git commit -m "feat(scripts): add NUC+MCseg dual-mode seg quality runner"`

---

#### A3｜建立 `scripts/qc_metrics.py`

**預期行為：**
讀取 `results/qc/` 下的 `.npy` 遮罩，計算每個 ROI × 每個方法的
FTC、NED（抽樣 200 細胞）、Artificial Co-expression Rate，
輸出 `results/qc_metrics.csv`（欄位：roi, method, n_cells, ftc, ned, coexp_rate）。

**相關檔案：** `scripts/qc_metrics.py`（新建）

- [ ] **Red**：在 `test_skill_scripts.py` 加入 `test_qc_metrics_columns`，
  使用 2 個假遮罩（100×100, 5 個細胞）+ 假 AnnData（10 個 bins, 20 基因），
  驗證輸出 CSV 包含正確欄位且 NED 介於 0–1。
  ```bash
  uv run pytest backend/tests/test_skill_scripts.py::test_qc_metrics_columns -v
  # 預期：FAILED
  ```

- [ ] **Green**：建立 `scripts/qc_metrics.py`（從 `skills/msseg-segment.md` Step 4 提取），
  關鍵修正：
  - NED 計算使用 `skimage.morphology.dilation` 而非雙層迴圈
  - `count_mat` 使用 `csr_matrix` sparse 格式，避免 `.todense()` OOM
  - 若基因不在 `var_names` 中，靜默跳過該對（不拋出 KeyError）
  ```bash
  uv run pytest backend/tests/test_skill_scripts.py::test_qc_metrics_columns -v
  # 預期：PASSED
  ```

- [ ] **Refactor**：NED 的 Hellinger 計算提升為模組層級函數 `_hellinger(p, q)`，
  避免在迴圈內重複定義（參考 CLAUDE.md §11 DRY 原則）。

- [ ] **Commit**：`git commit -m "feat(scripts): add FTC/NED/coexp QC metrics calculator"`

---

#### A4｜建立 `scripts/write_handoff.py`

**預期行為：**
讀取 `results/qc_metrics.csv`，計算各指標平均值，
依 AI 判斷規則自動推薦 `min_genes`、`max_pct_mt`，
寫入 `results/handoff_report.json`。

**相關檔案：** `scripts/write_handoff.py`（新建）

- [ ] **Red**：在 `test_skill_scripts.py` 加入 `test_handoff_report_keys`，
  驗證 JSON 包含 `segmentation_complete`、`roi_qc`、`recommended_analysis_params`、
  `masks_dir`、`binned_dir`、`tissue_profile` 六個頂層鍵。
  ```bash
  uv run pytest backend/tests/test_skill_scripts.py::test_handoff_report_keys -v
  # 預期：FAILED
  ```

- [ ] **Green**：建立 `scripts/write_handoff.py`，
  自動推薦邏輯：
  - `ned_delta >= 0.03` → quality = "good"；`0.01–0.03` → "marginal"；`<0.01` → "poor"
  - `coexp_mean > 0.06` → 建議減少 `voronoi_distance`，寫入 `warnings` 欄位
  ```bash
  uv run pytest backend/tests/test_skill_scripts.py::test_handoff_report_keys -v
  # 預期：PASSED
  ```

- [ ] **Commit**：`git commit -m "feat(scripts): add handoff report generator"`

---

#### A5｜修正 `skills/msseg-segment.md` Step 3

**預期行為：**
Step 3 原本錯誤呼叫 `run_mcseg_v2(roi, use_gpu=False, ...)` 需改為
呼叫 `scripts/seg_quality.py`，符合實際 API（`run_mcseg_v2(img, cfg)`）。

**相關檔案：** `skills/msseg-segment.md`

- [ ] 將 STEP 3「ROI 品質評估」的 Python 嵌入程式碼，
  改為呼叫已建立的腳本：
  ```bash
  uv run python scripts/seg_quality.py
  ```
  並移除 skill 內的 `run_mcseg_v2` 直接呼叫。

- [ ] 將 STEP 8 的手動計算程式碼改為：
  ```bash
  uv run python scripts/write_handoff.py
  ```

- [ ] **Commit**：`git commit -m "fix(skills): align msseg-segment step3/8 with actual API"`

---

### Phase B：msseg-analyze 支援腳本

---

#### B1｜建立 `scripts/build_full_adata.py`

**預期行為：**
讀取 `handoff_report.json` 取得 `binned_dir` 與 `masks_dir`，
呼叫 `count_rna_per_cell` 對所有遮罩進行計數，
合併多個 ROI 的 AnnData，寫入 `results/analysis/cellpose_cells.h5ad`。

**相關檔案：** `scripts/build_full_adata.py`（新建）

- [ ] **Red**：在 `test_skill_scripts.py` 加入 `test_build_adata_shape`，
  使用假遮罩（50×50, 3 cells）+ 假 bin matrix（9 bins, 5 genes），
  驗證輸出 h5ad 的 `n_obs = 3`、`n_vars = 5`。
  ```bash
  uv run pytest backend/tests/test_skill_scripts.py::test_build_adata_shape -v
  # 預期：FAILED
  ```

- [ ] **Green**：建立 `scripts/build_full_adata.py`：
  - 使用 `count_rna_per_cell(adata_path, mask_path, roi_x_px, roi_y_px, dilation_px=6)`
  - 多 ROI 時用 `anndata.concat(adatas, merge="same")`
  - 在 `adata.obs` 加入 `roi_name` 欄位，保留空間座標
  ```bash
  uv run pytest backend/tests/test_skill_scripts.py::test_build_adata_shape -v
  # 預期：PASSED
  ```

- [ ] **Refactor**：加入 ROI 座標偏移後負座標驗證（參考 CLAUDE.md §11）。

- [ ] **Commit**：`git commit -m "feat(scripts): add full AnnData builder from MCseg masks"`

---

#### B2｜建立 `scripts/run_analysis.py`

**預期行為：**
接收 QC 閾值參數，執行完整 Scanpy 流程：
`filter_cells → normalize → log1p → HVG → scale → PCA → neighbors → UMAP → leiden`，
中間每階段寫入 checkpoint h5ad，最終輸出 `cellpose_cells_clustered.h5ad`。

**相關檔案：** `scripts/run_analysis.py`（新建）

- [ ] **Red**：在 `test_skill_scripts.py` 加入 `test_run_analysis_umap_exists`，
  使用 100 細胞假 AnnData，驗證輸出包含 `X_umap`、`leiden` 欄位。
  ```bash
  uv run pytest backend/tests/test_skill_scripts.py::test_run_analysis_umap_exists -v
  # 預期：FAILED
  ```

- [ ] **Green**：建立 `scripts/run_analysis.py`，
  接受 CLI 參數：
  ```
  --min-genes INT    --max-pct-mt FLOAT    --min-counts INT
  --resolution FLOAT --n-pcs INT           --input PATH
  --output-dir PATH
  ```
  自動 resolution 邏輯（依細胞數）封裝為 `_auto_resolution(n_cells)`。
  ```bash
  uv run pytest backend/tests/test_skill_scripts.py::test_run_analysis_umap_exists -v
  # 預期：PASSED
  ```

- [ ] **Commit**：`git commit -m "feat(scripts): add scanpy analysis pipeline runner"`

---

#### B3｜在 `xenium_exporter.py` 新增 `export_from_mcseg()` 適配函數

**預期行為：**
新增 module-level 函數 `export_from_mcseg(adata, masks_dir, output_dir, pixel_size_um)`，
從 MCseg 遮罩提取輪廓（`skimage.measure.find_contours`），轉換為 GeoJSON，
並呼叫 `XeniumExporter` 建立 zarr bundle。

**相關檔案：** `backend/src/export/xenium_exporter.py`

- [ ] **Red**：在 `backend/tests/test_05_xenium.py` 加入 `test_export_from_mcseg_creates_files`，
  使用 50×50 假遮罩（3 cells）+ 假 AnnData，
  驗證輸出目錄包含 `cell_boundaries.geojson` 或 `cell_boundaries.parquet`。
  ```bash
  uv run pytest backend/tests/test_05_xenium.py::test_export_from_mcseg_creates_files -v
  # 預期：FAILED
  ```

- [ ] **Green**：在 `xenium_exporter.py` 末端新增：
  ```python
  def export_from_mcseg(
      adata: "anndata.AnnData",
      masks_dir: str | Path,
      output_dir: str | Path,
      pixel_size_um: float = VISIUM_UM_PX,
  ) -> Path:
      """MCseg 遮罩 → Xenium Explorer zarr + GeoJSON"""
      from skimage.measure import find_contours
      import json
      # 從遮罩提取多邊形輪廓（µm 座標）
      # 寫入 cell_boundaries.geojson
      # 呼叫 XeniumExporter.export(h5ad_path, output_dir)
  ```
  ```bash
  uv run pytest backend/tests/test_05_xenium.py::test_export_from_mcseg_creates_files -v
  # 預期：PASSED
  ```

- [ ] **Refactor**：輪廓提取邏輯提升為私有函數 `_masks_to_geojson(masks_dir, pixel_size_um)`，
  不在 `export_from_mcseg` 主體內重複定義。

- [ ] **Commit**：`git commit -m "feat(export): add export_from_mcseg adapter for MCseg masks"`

---

#### B4｜建立 `scripts/export_mcseg.py`

**預期行為：**
CLI 包裝器，呼叫 `export_from_mcseg`，同時儲存乾淨的 `msseg_final.h5ad`。
接受 `--format xenium|h5ad|both`，預設 `both`。

**相關檔案：** `scripts/export_mcseg.py`（新建）

- [ ] **Red**：在 `test_skill_scripts.py` 加入 `test_export_cli_both_format`，
  驗證 `--format both` 時同時產出 `msseg_final.h5ad` 與 `xenium/` 目錄。
  ```bash
  uv run pytest backend/tests/test_skill_scripts.py::test_export_cli_both_format -v
  # 預期：FAILED
  ```

- [ ] **Green**：建立 `scripts/export_mcseg.py`，
  使用 `argparse`，呼叫 `export_from_mcseg` 與 `adata.write`。
  ```bash
  uv run pytest backend/tests/test_skill_scripts.py::test_export_cli_both_format -v
  # 預期：PASSED
  ```

- [ ] **Commit**：`git commit -m "feat(scripts): add export CLI for Xenium + h5ad output"`

---

#### B5｜修正 `skills/msseg-analyze.md` Step 1、Step 7

**預期行為：**
- Step 1 改為呼叫 `scripts/build_full_adata.py`
- Step 7 改為呼叫 `scripts/export_mcseg.py --format xenium`
- 移除 skill 內嵌的 `export_to_xenium(...)` 呼叫（該函數已改名）

**相關檔案：** `skills/msseg-analyze.md`

- [ ] 將 STEP 1 的計數程式碼改為：
  ```bash
  uv run python scripts/build_full_adata.py
  ```
- [ ] 將 STEP 7 改為：
  ```bash
  uv run python scripts/export_mcseg.py --format xenium --output results/export
  ```
- [ ] 將 STEP 8 改為：
  ```bash
  uv run python scripts/export_mcseg.py --format h5ad --output results/export
  ```

- [ ] **Commit**：`git commit -m "fix(skills): align msseg-analyze steps with actual script API"`

---

### Phase C：整合測試

---

#### C1｜端對端 Smoke Test（使用現有 CRC ROI）

**預期行為：**
使用 `config/pipeline.yaml` 中已定義的 `test` ROI（1491×1210 px），
完整跑一遍 msseg-segment 的 Step 2–8，不報錯，產出 `handoff_report.json`。

- [ ] 執行：
  ```bash
  uv run python scripts/roi_sampler.py
  uv run python scripts/seg_quality.py
  uv run python scripts/qc_metrics.py
  uv run python scripts/write_handoff.py
  cat results/handoff_report.json
  ```
  驗證：JSON 有效，`segmentation_complete: true`，`n_rois_evaluated >= 1`

- [ ] **Commit**：`git commit -m "test(e2e): validate msseg-segment pipeline on CRC test ROI"`

---

#### C2｜更新 `summary.md`

- [ ] 在 `summary.md` 新增「Skill 支援腳本」章節，列出所有新建腳本與用途。

- [ ] **Commit**：`git commit -m "docs: update summary with skill support scripts"`

---

## 執行順序

```
A1 → A2 → A3 → A4 → A5（Phase A 序列執行）
B1 → B2 → B3 → B4 → B5（Phase B 序列執行）
C1 → C2（整合測試最後執行）
```

Phase A 與 Phase B 可平行開發，但 C1 需等 A 全部完成。

---

這個計畫看起來沒問題嗎？準備好執行 `/sp-executing-plans` 指令了嗎？
