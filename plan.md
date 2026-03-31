# MSseg 實作計畫

> 版本：v1.0 (2026-03-31)
> 基於 visiumHD_pipeline_3，整合 MCseg v2 多模型 Voronoi 集成分割引擎。

---

## 一、專案目標與方向

在 `visiumHD_pipeline_3` 基礎上，以三項核心升級建構 MSseg：

| # | 升級項目 | 實作策略 |
|---|----------|----------|
| 1 | **MCseg v2 分割引擎** | 以 cyto3 三直徑（13/17/22px）+ Hematoxylin pass + Voronoi 擴張取代 LOGIC_A 雙尺寸策略。LUAD PQ 從 0.432 提升至 0.554（+28%）。 |
| 2 | **移除 Proseg** | 刪除 Stage 2.5 全部程式碼（backend/src/proseg/、api/proseg_rna.py），去除對外部 Rust binary 的依賴。 |
| 3 | **per-ROI 參數覆寫** | 保留每個 ROI 獨立覆寫 MCseg v2 參數的能力（dia_mid、voronoi_distance 等），透過 `roi_overrides` API 欄位實現。 |

---

## 二、系統架構（MSseg v1）

```
MSseg/
├── CLAUDE.md                    # 開發維護規範
├── plan.md                      # 本文件
├── README.md                    # 使用說明
├── pyproject.toml               # Python 依賴（msseg 1.0.0）
│
├── config/
│   ├── pipeline.yaml            # 專案參數（segmentation.mcseg_v2 子段）
│   └── profiles/                # Tissue Profile（CRC, LUAD）
│
├── backend/
│   └── src/
│       ├── api/
│       │   └── segmentation.py  # MCseg v2 Pydantic 參數模型
│       ├── segmentation/
│       │   └── cellpose_runner.py  # MCseg v2 核心實作
│       ├── cellpose_counter/    # Stage 2: RNA 計數
│       ├── roi/                 # Stage 0: ROI 裁切
│       ├── analysis/            # Stage 3: Scanpy
│       └── export/              # Stage 4: Browser 匯出
│
└── frontend/
    └── src/
        └── pages/
            └── Stage1_Segmentation.tsx  # MCseg v2 參數 UI
```

---

## 三、MCseg v2 演算法詳規

### 核心流程（`run_mcseg_v2`）

```text
輸入：he_image（H×W×3 RGB）+ params dict

Step 1：CLAHE 前處理
  apply_clahe(image, clip_limit=clahe_clip_limit)
  → LAB 色彩空間均衡化，增強局部對比

Step 2：多 pass cyto3 推論
  pass 1：dia_small=13px  → small_mask
  pass 2：dia_mid=17px    → mid_mask（主要基底）
  pass 3：dia_large=22px  → large_mask
  pass 4（選配）：cyto3 on Ruifrok Hematoxylin 通道 → hema_mask
  pass 5-7（選配）：cpsam × 3 → cpsam_masks

Step 3：merge_masks_fast 集成
  IoU 閾值去重，保留最大非重疊細胞集合
  → merged_mask

Step 4：Voronoi 擴張（取代 expand_labels）
  voronoi_expand(merged_mask, max_distance=voronoi_distance)
  每個背景像素 → Voronoi 最近細胞（distance capped）
  ⚠️ 無重疊保證（expand_labels 在複雜拓撲下可重疊）
  → voronoi_mask

Step 5（選配）：轉錄本密度補救
  find_transcript_seeds(voronoi_mask, vhd_pseudo_transcripts.csv)
  高密度無細胞區域 → 植入漏偵測細胞

Step 6：尺寸過濾
  clean_mask(mask, min_size, max_size)
  → final_mask（H×W int32，0=背景）
```

### 與 LOGIC_A 對比

| 特性 | LOGIC_A（舊）| MCseg v2（新）|
|------|-------------|--------------|
| 模型 | cyto2 雙直徑 | cyto3 三直徑 + Hema pass |
| 合併策略 | 覆蓋率投票 | IoU 去重集成 |
| 擴張方式 | expand_labels | Voronoi（無重疊）|
| 轉錄本補救 | 無 | 選配（vhd_pseudo_transcripts.csv）|
| LUAD PQ@0.5 | 0.432 | **0.554**（+28%）|

---

## 四、各 Stage 詳細規格

### Stage 0：ROI 定義與全自動裁切

BTF OME-TIFF tile-by-tile 記憶體串流。ROI 定義格式：

```yaml
rois:
  - name: "test"
    tissue: CRC
    x: 48128
    y: 12657
    width_px: 1491
    height_px: 1210
    pixel_size_um: 0.2737
```

### Stage 1：MCseg v2 分割

**前端 UI（Stage1_Segmentation.tsx）**：

- 左欄：GPU 設定 + cyto3 多直徑（dia_small/mid/large）+ use_hematoxylin + use_cpsam
- 右欄：Voronoi/後處理 + Cellpose QC + 快速預設值
- ROI override 表：name / cells / dia_mid / voronoi_distance / overridden

**快速預覽（`run_preview_patch`）**：
- 自動停用 cpsam 和 transcript rescue（加速預覽）
- 回傳 raw、CLAHE enhanced、Hematoxylin 三張影像

**前處理預覽（`preview_preproc`）**：
- 回傳 `raw_b64`、`clahe_b64`、`hema_b64` 三張 base64 影像

### Stage 2：RNA 計數（稀疏矩陣法）

`counter.py`：

1. MCseg v2 遮罩（已含 Voronoi 擴張）
2. `expand_labels(distance=dilation_px=6)`：補填 Voronoi 間隙 bins
3. 稀疏矩陣聚合 `A @ adata_002um.X`

輸出：`cellpose_cells.h5ad`（cells × 18K genes）

### Stage 3：Scanpy 下游分析

Scanpy pipeline：QC → Normalize → log1p → HVG → PCA → kNN → UMAP → Leiden → Marker Gene → TME Panel

TME panels 由 `config/profiles/{tissue}.yaml` 動態讀取（`_build_tme_config(config)`）。

### Stage 4：Browser 匯出

- Xenium Explorer：`skimage.find_contours` → GeoJSON
- Loupe Browser：barcode whitelist + cluster assignment

---

## 五、config/pipeline.yaml segmentation 結構

```yaml
segmentation:
  mcseg_v2:
    use_gpu: true
    batch_size: 4
    dia_small: 13.0
    dia_mid: 17.0
    dia_large: 22.0
    use_hematoxylin: true
    use_cpsam: false
    voronoi_distance: 9
    min_size: 20
    max_size: 6000
    flow_threshold: 0.4
    cellprob_threshold: -2.0
    clahe_clip_limit: 3.0
    use_transcript_rescue: true
  output:
    mask_filename: segmentation_masks.npy
    mask_tif_filename: segmentation_masks.tif
```

---

## 六、Tissue Profile MCseg v2 覆寫欄位

Profile YAML 可覆寫 `segmentation.mcseg_v2` 下的任何參數。例如 LUAD 可降低 `dia_mid` 配合較小的肺泡細胞：

```yaml
# config/profiles/luad.yaml
segmentation:
  mcseg_v2:
    dia_small: 10.0
    dia_mid: 14.0
    dia_large: 18.0
    clahe_clip_limit: 2.5
    cellprob_threshold: -1.5
```

---

## 七、已完成項目（2026-03-31）

- [x] `backend/src/segmentation/cellpose_runner.py`：MCseg v2 全實作
- [x] `backend/src/api/segmentation.py`：MCseg v2 Pydantic 參數模型
- [x] `backend/main.py`：移除 Proseg router
- [x] `backend/src/proseg/`：已刪除
- [x] `backend/src/api/proseg_rna.py`：已刪除
- [x] `frontend/src/App.tsx`：移除 Proseg 路由
- [x] `frontend/src/components/layout/TopNav.tsx`：移除 Proseg Stage，更新品牌
- [x] `frontend/src/components/layout/Sidebar.tsx`：更新標題 MSseg / MCseg v2
- [x] `frontend/src/pages/Stage1_Segmentation.tsx`：MCseg v2 參數 UI
- [x] `frontend/src/pages/Stage25_ProsegRNA.tsx`：已刪除
- [x] `config/pipeline.yaml`：更新 segmentation 段為 mcseg_v2
- [x] `pyproject.toml`：name=msseg, version=1.0.0

---

## 八、未來優化項目

- [ ] LUAD / CRC 組織別 MCseg v2 dia 預設值細調（根據 PQ benchmark）
- [ ] cpsam 整合完整測試（目前可選，預設 off）
- [ ] 轉錄本密度補救 end-to-end 驗證
- [ ] Batch Mode：多 ROI 一鍵排程
- [ ] Stage 3 PAGA 軌跡推論（選配）
