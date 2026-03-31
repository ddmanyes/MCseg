# MSseg — 前端建置與驗證計畫（歷史紀錄）

> 原始建立：2026-03-03（visiumHD_pipeline_2 基礎架構）
> MSseg v1.0 遷移：2026-03-31（MCseg v2 整合，Proseg 移除）
> 狀態：✅ 完成

---

## 專案現況

後端模組（FastAPI + 8 個 API 路由 + WebSocket log 串流）已 **全部完成** 。前端原始碼也已完整撰寫（19 個檔案），結構如下：

| 類別 | 檔案 |
| ---- | ---- |
| 入口 | `main.tsx`, `App.tsx`, `index.html`, `index.css` |
| 佈局 | `Sidebar.tsx`, `Header.tsx` |
| 共用元件 | `StageCard.tsx`, `Terminal.tsx` |
| 頁面（8 個） | `DataSetup.tsx` + `Stage0_ROI.tsx` ~ `Stage5_Export.tsx` |
| 狀態管理 | `pipelineStore.ts` (Zustand) |
| API 層 | `client.ts` (Axios) |
| WebSocket | `useStageLog.ts` (custom hook) |
| 型別 | `pipeline.ts` |
| 配置 | `vite.config.ts`, `tailwind.config.js`, `postcss.config.js`, `tsconfig.json` |

---

## 執行清單

### Phase 1：安裝與編譯

- [x] 在 `frontend/` 目錄執行 `npm install`（197 packages）
- [x] 執行 `npx tsc --noEmit` 檢查 TypeScript 編譯 → **零錯誤**
- [x] 無需修復，所有型別正確

### Phase 2：啟動驗證

- [x] 執行 `npm run dev` 啟動 Vite dev server（port 3000）→ 698ms 就緒
- [x] 使用瀏覽器訪問 `http://localhost:3000` → 正常載入
- [x] 驗證所有頁面正常：DataSetup + Stage 0 ~ Stage 5 均可導航
- [x] 深色主題正確渲染
- [x] `npm run build` 生產建構成功（2376 modules, 2.78s）

### Phase 3：Git 初始化

- [x] `git init` + 初始提交（64 files, commit `00b7771`）

### Phase 4：資料自動發現功能

- [x] 後端 `discovery.py` — 掃描邏輯
- [x] 後端 `data.py` — API 路由（scan / apply / status）
- [x] 前端 `DataSetup.tsx` — 資料設定 UI 頁面
- [x] 修改路由、Sidebar、Header、client.ts
- [x] Git commit: `8955455`

---

## 驗證結果

| 項目 | 結果 |
| ---- | ---- |
| `npm install` | ✅ 197 packages |
| `tsc --noEmit` | ✅ 零錯誤 |
| `npm run dev` | ✅ port 3000 |
| 瀏覽器載入 | ✅ 所有 8 頁正常 |
| `npm run build` | ✅ 2376 modules, 2.78s |
| Git 提交 | ✅ `00b7771` + `8955455` |

---

## 變更紀錄

| 日期 | 變更內容 |
| ---- | ---- |
| 2026-03-03 | 初始建立，完成前端安裝、編譯、驗證與 Git 初始化 |
| 2026-03-03 | 新增資料自動發現功能（DataSetup 頁面 + 後端掃描邏輯） |
| 2026-03-04 | 完成全域整合架構：導入 Zarr Tiling 分塊合併策略防爆記憶體、套用 Dask-Expr 查詢防護（Fix `legacy implementation`）與 React 掛載狀態優化 |
| 2026-03-19 | **Tissue Profile 系統**：新增 `config/profiles/`（crc.yaml, luad.yaml），`config.py` 支援三層 merge，換組織只需改 `tissue_profile` 一行 |
| 2026-03-19 | **RNA 計數優化**：`counter.py` 新增 `expand_labels` dilation（6px = 1.64 µm），PQ 0.397 → 0.432（+9%） |
| 2026-03-19 | **TME Panel 動態化**：`analysis/pipeline.py` 新增 `_build_tme_config(config)`，TME panels 改由 profile YAML 驅動，不再硬編碼 |
| 2026-03-19 | **Stage 2.5 Method E + Hungarian 回填**：`runner.py` 新增 `_combine_refine_rna()`（Proseg centroid → Cellpose lookup，全向量化）；`pipeline.yaml` 更新最佳參數（md10_c004_sp20，roi1 PQ=0.502）；`hungarian_refine: true` 開關控制 |
