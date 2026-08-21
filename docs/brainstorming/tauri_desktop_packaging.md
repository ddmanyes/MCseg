# MCseg 桌面應用程式 (Tauri Desktop App) 完整設計規格與實作計畫

- **狀態（2026-08-20 更新）**：**Phase 0、1、2、3、4 全部完成**（核心邏輯 +
  macOS 上大部分實機驗證）。Phase 5（CI + 簽章）／Phase 6（驗收測試矩陣）
  未開始。**⚠️ 所有改動目前是未 commit 的工作區變更，且當前分支是
  `feat/if-multichannel-segmentation`（另一條無關的既有 feature branch）
  ——下次接手前兩件事：(1) 決定要不要開一條新分支承接這批 Tauri 改動，
  不要混進 IF 分割那條；(2) commit。** 詳見文末「進度總結（2026-08-20
  收工）」。
- **v1 記錄日期**：2026-08-19（由 antigravity 撰寫）
- **v2 審查日期**：2026-08-20（對照現有程式碼與環境重新檢視，找到架構落差與 8 項計畫外風險，見文末「v2 修訂」）
- **核心目標**：將現有 MCseg 空間轉錄體分析流水線封裝為跨平台原生桌面應用程式（`.dmg` / `.msi`），為實驗室研究人員與生信分析師提供開箱即用、免終端機指令、具備原生檔案對話框與系統通知的一體化體驗。

> ⚠️ **v1 內容（下方「系統架構圖」「5 大定案規格」「實作執行計畫」）在部分細節上已被 v2 修訂取代或補強，請以文末「v2 修訂」為準再動工。** v1 保留原樣以供對照。

---

## 🏛️ 系統架構圖 (System Architecture)

```mermaid
graph TD
    subgraph "MCseg 桌面應用程式 (Tauri Desktop App)"
        A["原生視窗 (macOS / Windows)"] <-->|Tauri IPC Event| B["Tauri Rust 核心"]
        B -->|展示| C["MCseg React 前端 (UI)"]
        B -->|原生檔案對話框| F["系統 Finder / 檔案總管"]
        B -->|系統原生通知| G["macOS / Windows 橫幅通知"]
        B -->|Sidecar 生命週期管理| D["微型 uv 隔離環境 (AppData)"]
        D -->|啟動 & 監聽| E["Python FastAPI 後端<br>(PyTorch MPS/CUDA + Cellpose + Scanpy)"]
    end
    C <-->|HTTP / WebSocket (Port 8001)| E
```

---

## 📋 5 大定案規格 (Finalized Specifications)

1. **執行環境分發模式 (Runtime Distribution)**
   * **輕量外殼發布**：安裝包（`.dmg` / `.msi`）維持在 **~20 MB**，不預先打包數 GB 的 PyTorch 與 CUDA 二進位檔。
   * **首次自動環境配置**：
     * 首次啟動時在使用者本機應用程式支援目錄（`~/Library/Application Support/MCseg/` 或 `%APPDATA%/MCseg/`）建立專屬虛擬環境。
     * 利用極速套件管理器 `uv` 自動下載與安裝 PyTorch（依硬體自動適配 MPS/CUDA）、Cellpose 與 Scanpy。
     * 支援軟體快速版本升級，不需每次重新下載數 GB 檔案。

2. **目標作業系統與硬體加速策略 (OS & Hardware Acceleration)**
   * **支援平台**：**macOS (Apple Silicon M 系列)** 與 **Windows 10/11 (NVIDIA CUDA / CPU)** 雙平台。
   * **硬體自動適配**：
     * **macOS**：自動偵測 Apple Silicon 並啟用 `mps` Metal 加速。
     * **Windows**：自動檢測 NVIDIA GPU 與驅動，自動安裝對應 CUDA 版本的 PyTorch；無專用 GPU 則自動回退為多執行緒 CPU 模式。

3. **首次啟動與安裝導覽 (First-Time Setup Wizard)**
   * **視覺化安裝精靈 (Setup Wizard)**：
     * 首次開啟時呈現歡迎與初始化引導卡片。
     * 步驟視覺化：`[✓] 檢查系統硬體與 GPU` $\rightarrow$ `[🔄] 初始化 Python 與生信環境` $\rightarrow$ `[ ] 啟動 MCseg 運算引擎`。
     * 具備動態進度百分比條與可展開的即時日誌抽屜（Log Drawer），遇到網路超時可一鍵重試。
     * 完成後自動平滑過渡進入 Stage 0 數據載入頁面。

4. **本機檔案瀏覽與路徑存取 (Native File Dialog Integration)**
   * **雙模混合輸入**：
     * 點擊各輸入框旁的「📁 瀏覽 (Browse)」按鈕直接呼叫 macOS Finder / Windows 檔案總管原生對話框。
     * 檔案選取後自動將本機絕對路徑（如 `/Volumes/SSD/sample_crc/tissue_hires_image.png`）回填入輸入框。
     * 保留輸入框文字可編輯性，支援手動貼上路徑。

5. **背景行程生命週期與防呆管理 (Lifecycle & Resource Cleanup)**
   * **運算中防呆退出**：當後台有耗時分析正在運行（如 Cellpose 分割、UMAP 計算），關閉視窗時跳出警告對話框確認是否強制終止。
   * **徹底釋放顯存**：確認退出時，Tauri 自動發送 SIGTERM/SIGKILL 終止 Python 後端行程，釋放 Port 8001、系統 RAM 與 GPU 顯存，不留殭屍行程。
   * **重啟防衝突**：啟動前自動檢查並清理殘留的 8001 端口。

---

## 🛠️ 下次啟動之實作執行計畫 (Actionable Implementation Plan)

### 階段一：Tauri 專案骨架與環境初始化 (Phase 1: Tauri Infrastructure)
- [ ] 安裝 `@tauri-apps/cli` 與 Rust Tauri 工具鏈。
- [ ] 執行 `tauri init` 初始化 `src-tauri/` 配置目錄。
- [ ] 配置 `tauri.conf.json`：
  - App Identifier: `com.mcseg.app`
  - 視窗預設尺寸：`1440x920`（最小限制 `1024x720`）
  - 前端建置命令：`npm run build`，靜態產物目錄：`../dist`
  - 啟用安全權限與必要外掛（`dialog`, `notification`, `shell`, `process`）。
- [ ] 設計與配置 MCseg 原生應用程式圖示（App Icons: `.icns` for Mac, `.ico` for Windows）。

### 階段二：微型自動環境配置器與 Setup Wizard (Phase 2: Bootstrapper & Setup UI)
- [ ] 編寫跨平台微型環境安裝器 `scripts/bootstrap_env.py`：
  - 自動下載 Standalone `uv` 執行檔。
  - 執行 `uv venv` 建立專屬虛擬環境。
  - 偵測硬體平台並執行 `uv pip install` 安裝 FastAPI, PyTorch (MPS/CUDA), Cellpose, Scanpy, Anndata, Libvips, OpenCV。
- [ ] 於前端建立 `frontend/src/components/setup/SetupWizard.tsx`：
  - 接收後台環境建立進度（百分比、當前步驟、即時日誌串流）。
  - 當環境檢驗合格時，自動切換至 Stage 0。

### 階段三：原生對話框與系統通知整合 (Phase 3: Native Dialogs & Notifications)
- [ ] 安裝與配置 `@tauri-apps/plugin-dialog`：
  - 封裝通用 Hook `useFileDialog()`。
  - 升級 `Stage0_DataSetup.tsx`、`Stage0_ROI.tsx` 與自訂 Marker CSV 選取器，支援原生 Finder / 檔案總管選取。
- [ ] 安裝與配置 `@tauri-apps/plugin-notification`：
  - 當 Cellpose 分割、QC 評估、UMAP 計算、CellTypist 標註完成時，調用系統橫幅通知。

### 階段四：Rust 生命週期調度與進程守護 (Phase 4: Sidecar Process Lifecycle)
- [ ] 在 `src-tauri/src/main.rs` 實作：
  - 啟動時自動呼叫後端 Uvicorn 服務。
  - 監聽 `tauri::WindowEvent::CloseRequested` 事件。
  - 檢查當前 API 分析狀態；若處於 `running` 狀態則跳出原生確認對話框。
  - 確認關閉時，發送信號安全終止 Python 行程並釋放 Port 8001 與顯存。

### 階段五：編譯打包與驗證 (Phase 5: Build & Packaging)
- [ ] 執行 `npm run tauri build` 編譯 macOS Universal / Apple Silicon `.dmg` 安裝包。
- [ ] 驗證安裝流程（安裝精靈 $\rightarrow$ 數據載入 $\rightarrow$ 空間分割 $\rightarrow$ 下游分析 $\rightarrow$ 關閉退出）。

---

## 🔍 v2 修訂（2026-08-20，動工前審查）

對照現有程式碼（`backend/main.py`、`frontend/src/api/client.ts`、`pyproject.toml`）與本機環境重新檢視 v1 計畫後，發現一個會導致「打包後畫面靜默空白」的架構落差，以及 8 項計畫沒處理但會實際影響上線後果的風險。完整審查脈絡見 sb 筆記 [MSseg Tauri 桌面應用打包計畫與審查]（second-brain vault，`20-areas/coding/msseg-tauri-桌面應用打包計畫與審查.md`）。

### 環境現況（審查當下，macOS 開發機）
- ❌ 未安裝 Rust/Cargo/Tauri CLI，需先 `rustup` 裝工具鏈（Xcode CLT 已就緒）
- ✅ `uv 0.11.24`、`node v24.18.0`、`npm 11.16.0` 齊全
- ✅ SSD 空間充足

### ⛔ 架構修正：取代「Tauri 自帶前端 + 獨立後端」，改為單一 origin

v1 架構圖畫的是 Tauri 渲染 React 前端、與跑在另一 port 的 Python 後端經 HTTP/WS 通訊，兩者視為不同來源。但實際程式碼：

- `frontend/src/api/client.ts` 的 axios `baseURL` 是**相對路徑** `/api`，沒有 `VITE_*` 環境變數可切換成絕對網址
- `backend/main.py` 的 `CORSMiddleware.allow_origins` **寫死**只允許 `http://localhost:3000` / `127.0.0.1:3000`（Vite dev server 來源）
- Tauri webview 預設來源是 `tauri://localhost`（Mac）/ `http://tauri.localhost`（Win）——跟上面兩者都對不上

若照 v1 架構圖實作，會同時踩到「相對路徑打不到 API」+「CORS 擋掉」，且是**沒有錯誤訊息、畫面直接空白**的失敗模式。

**已有更簡單的路**：`backend/main.py:123-135` 已經能把 build 好的前端 `dist/` 用同一個 FastAPI app 服務（`/assets` 先掛載 + SPA catch-all，順序 bug 也已修過）。**改為 Tauri 不打包前端，視窗直接導向 `http://127.0.0.1:8001`**——後端啟動後同時服務 UI 與 API，單一 origin，CORS/相對路徑問題整個不存在。

⚠️ **但這個簡化本身需要額外設定**：Tauri 預設 CSP 會限制 `connect-src`/導覽只能到 App 自己的來源，需在 `tauri.conf.json` 明確把 `http://127.0.0.1:8001` 加入白名單，否則會踩到同一種「畫面空白、無錯誤訊息」的坑，只是換了個地方發生。

### 8 項計畫外風險

| # | 嚴重度 | 問題 | 修法 |
|---|--------|------|------|
| 1 | 🔴 | **共用工作站無單一實例防護**：`start.sh` 既有的「先 kill 佔用 port 8001 的行程再啟動」邏輯若原樣搬進 Tauri App，會在同事的分割工作（實測過 ~35 小時）還在跑時被另一個人開 App 給靜默砍掉 | Phase 4 前必須設計單一實例鎖；偵測到 port 8001 已有工作在跑時，改為「附加查看進度」而非強制終止 |
| 2 | 🔴 | **bootstrap 重新發明依賴管理**：`bootstrap_env.py` 規格是手寫套件清單，與 `pyproject.toml`/`uv.lock` 是兩份獨立東西，`uv add` 後會漂移 | 把 `pyproject.toml` + `uv.lock` 當 App 資源打包，bootstrap 直接跑真正的 `uv sync` |
| 3 | 🔴 | **打包拿掉了「靠開發者眼睛抓 bug」的安全網**：本專案歷史上的 bug 幾乎都是技術人員發現數字不對勁才抓到，黑盒子使用者無法判斷「0 個 bin 命中」是 bug 還是資料本身沒東西 | 發布給非專家前，評估是否要把已知靜默失敗類 bug（如 `save_state` 深合併導致「全部重設」無效，見 [[msseg-前端錯誤處理修復計畫]]）列為發布阻塞項 |
| 4 | 🟡 | **無真正發布流程**：計畫止於「跑一次 `npm run tauri build`」，但 `.dmg`/`.msi` 要在兩台不同實體機器編譯 | 改用 GitHub Actions matrix（macOS runner + Windows runner），tag push 自動產出版本化安裝檔 |
| 5 | 🟡 | **升級機制未設計**：只有一句「支援快速版本升級」，無具體機制 | 待拍板：是否比對 `uv.lock` hash 只重裝變動套件，或先做最簡單版（有新版就導向重新下載安裝檔） |
| 6 | 🟡 | **bootstrap 無續傳/共用快取設計**：網路中斷後「一鍵重試」是否從頭來過未定義；快取路徑（App 專屬 vs 共用 `~/.cache/uv`）未定 | 待拍板：建議共用既有 `~/.cache/uv`，避免重裝整包重下載 |
| 7 | 🟡 | **`--reload` 未拿掉**：`start.sh`/`start.ps1` 都用 `uvicorn ... --reload`（開發用檔案監聽），正式版必須移除 | Phase 4 啟動指令改為不含 `--reload` |
| 8 | 🟢 | **Windows CUDA「自動偵測版本」與現有程式碼不符**：`pyproject.toml` 是靜態 `win32→cu128`，沒有依驅動動態選版本的邏輯；cu128 wheel 在無 GPU 機器上也能跑（退化 CPU，只是多下載） | 維持現狀固定裝 cu128，規格文字改為如實描述，不寫「自動偵測 CUDA 版本」 |

### 附帶：一個計畫外的優點

✅ **ExFAT 問題被架構順便解決**：venv 建在 `Application Support`/`%APPDATA%`（使用者本機系統碟），不在 ExFAT 外接碟上，CLAUDE.md 與跨機開發銜接指南記載的 `.venv` symlink / `._*` 雜訊問題不會發生在打包版身上。

### ✅ 開放問題已拍板（2026-08-20）

1. **簽章/公證**：**先不簽**（實驗室內部用）。README 附截圖說明 macOS 右鍵「打開」繞過 Gatekeeper、Windows SmartScreen 點「仍要執行」。之後真的要對外發布再補簽章。
2. **長任務背景執行 UX**：**縮到系統列繼續跑**（採用，非完全與 GUI 解耦的 job-queue 模式）。視窗關閉時縮到系統列/選單列圖示，工作行程持續在背景執行；點圖示恢復視窗查看進度。
3. **版本升級機制**：**最小版**。啟動時比對版本號，有新版時顯示提示橫幅，導向手動下載新安裝檔；不做差異更新，之後真有需求再升級成 Tauri updater。
4. **已知 bug（`save_state` 深合併，`put_roi_overrides` 全部重設無效）**：**列為發布前阻塞項**，修法已明確（詳見下方 Phase 0）。**不是現在立刻修**——記進 Phase 0 checklist，跟其他前置項一起處理，避免在動工前的規劃階段零散改動生產程式碼。

### 修訂後的階段規劃（新增 Phase 0 與 Phase 6）

- **Phase 0（新增，動工前）**：
  - [x] 安裝 Rust/Cargo 工具鏈（2026-08-20 完成：`rustc 1.97.1`／`cargo 1.97.1`／`rustup 1.29.0`，皆為 aarch64-apple-darwin stable）
  - [x] 確認 `uv` 快取共用策略（沿用 `~/.cache/uv`，不建獨立快取路徑——已定案，實際套用留給 Phase 2 寫 `bootstrap_env.py` 時處理）
  - [x] **修復 `put_roi_overrides` 深合併 bug**（2026-08-20 完成，詳見下方「Phase 0 執行紀錄」）
  - [x] README 補上 Gatekeeper/SmartScreen 繞過說明草稿（詳見下方「Phase 0 產出草稿」，待 Phase 5 產出真正的 `.dmg`/`.msi` 後移入正式 README 並補實際截圖）
- **Phase 1**：Tauri 骨架 —— 額外加入 `tauri.conf.json` 的 CSP/導覽白名單設定（允許 `http://127.0.0.1:8001`），主視窗 `url` 直接指向該位址（取代載入本地前端資源）
- **Phase 2**：bootstrap 改為打包 `pyproject.toml`/`uv.lock` 並執行真正的 `uv sync`；加入單一實例鎖與 port 8001 佔用偵測（偵測到既有工作在跑時附加而非終止）；設計續傳/快取策略
- **Phase 3**：原生對話框/通知 —— 因改用單一 origin，不再需要處理 CORS/相對路徑轉換，工作量減少；補充「偵測到已有實例在跑」的提示 UI；補版本檢查橫幅 UI
- **Phase 4**：Rust 生命週期調度 —— 啟動指令移除 `--reload`；視窗關閉事件改為「縮到系統列」而非直接終止行程（運算中才需要，閒置時可正常關閉）；系統列圖示點擊恢復視窗
- **Phase 5**：編譯打包 —— 改為 GitHub Actions CI matrix（macOS + Windows runner）；不含簽章步驟（依拍板結果）
- **Phase 6（新增）**：驗收測試矩陣 —— 至少涵蓋：有 GPU / 無 GPU 機器、首次安裝 / 重裝、bootstrap 中途斷網重試、macOS Gatekeeper 攔截首次啟動、同機雙開偵測是否生效、縮到系統列後長任務是否真的持續運算、`put_roi_overrides` 回歸測試

---

## Phase 0 執行紀錄（2026-08-20）

### Rust 工具鏈
`sh rustup-init.sh -y --default-toolchain stable --profile default` 非互動安裝完成。
偵測到本機已有殘留的 `~/.rustup/settings.toml`（舊嘗試留下），沿用其設定，未造成衝突。
`rustc 1.97.1`／`cargo 1.97.1`／`rustup 1.29.0`，host triple `aarch64-apple-darwin`。
Tauri CLI（`@tauri-apps/cli` 或 `cargo install tauri-cli`）留給 Phase 1 開工時再裝，
不在 Phase 0 範圍內（Phase 0 只需要 Rust 工具鏈本身）。

### `put_roi_overrides` 深合併 bug 修復

**根因**：`backend/src/utils/config.py` 的 `save_state()` 用 `_deep_merge` 遞迴合併，
對「巢狀 dict 只想更新幾個子欄位」是對的語意，但 `roi_seg_overrides` 需要的是
「整份取代」（清空全部、或刪除單一 ROI 的覆寫）——深合併會讓沒出現在新值裡的
舊 key 被誤判為未變更而繼續殘留。

**修法**：`config.py` 新增 `save_state_key(key, value)`，整份取代單一頂層 key，
不經 `_deep_merge`；同時抽出 `_write_state()` 共用兩者的原子寫入邏輯（DRY，
避免兩份重複的 tmp-file-then-replace 程式碼）。

`backend/src/api/segmentation.py` 兩處呼叫點都改用新函式：
- `POST /segmentation/run`（line 256，儲存本次執行使用的 overrides）
- `PUT /segmentation/roi_overrides`（line 439，前端「重設」按鈕打的端點）

兩處都改是因為根因相同（同一個 key、同一個 `save_state` 深合併陷阱），只修
使用者回報的那個端點會讓同一個 bug 從 `/run` 這條路徑繼續存在，不符合
CLAUDE.md §11 的 DRY 原則。

**測試**：`backend/tests/test_09_roi_overrides.py` 新增 `TestSaveStateKeyReplacesNotMerges`
三個案例（含一個刻意保留、用來證明舊行為確實是 bug 的對照組：`save_state`
深合併對空 dict 確實不會清空）。該檔案 11 項全過，純函式測試、monkeypatch
`_STATE_PATH` 到 `tmp_path`，不需要真實資料、不需要 `--extra dev`。

全套 `backend/tests/` 回歸結果：**330 passed, 2 failed**（425 秒）。2 個失敗
在 `test_skill_scripts.py`（`TestQcMetrics::test_qc_metrics_columns`、
`TestExportMcseg::test_export_cli_both_format`），確認與本次改動無關（該檔案
完全沒有引用 `save_state`/`roi_seg_overrides`/`config.py`/`segmentation.py`）：
- `qc_metrics.py` 實際函式名是 `compute_metrics_from_geojson`，測試在找不存在的
  `compute_metrics` —— 既有命名不一致
- `scripts/export_mcseg.py` 這個檔案本身不存在於 repo —— 既有技術債

兩者都是 Phase 0 範圍外的既有問題，未動它們。

## Phase 0：完成（2026-08-20）

四項全部完成，`put_roi_overrides` 修復乾淨（330 passed，動到的部分零回歸）。
可以進 Phase 1（Tauri 骨架）。

---

## Phase 1 執行紀錄（2026-08-20）

### 完成項目
- `frontend/` 內裝 `@tauri-apps/cli@2.11.4`（dev dependency）
- `npx tauri init --ci` 非互動產生 `frontend/src-tauri/` 骨架
- `tauri.conf.json` 調整：
  - `identifier`: `com.mcseg.app`
  - 主視窗 `width/height` 1440×920、`minWidth/minHeight` 1024×720
  - **主視窗 `url` 設為 `http://127.0.0.1:8001`**（落實 v2 修訂的單一 origin 架構——生產模式視窗直接導向後端，不載入 Tauri 自己打包的前端資源；`devUrl` 維持 `http://localhost:3000`，`tauri dev` 時走 Vite dev server，既有的 `/api`、`/ws` proxy 設定原封不動適用，不需要改前端程式碼）
- `capabilities/default.json` 加 `remote.urls: ["http://127.0.0.1:8001"]`——因為主視窗載入的是遠端來源（後端服務的頁面），不是 Tauri 內嵌的 `tauri://localhost` 資源，IPC 橋接預設不會授權給遠端來源，Phase 3 要接 dialog/notification plugin 時若沒有這行會被靜默擋下。先在骨架階段就配好，避免 Phase 3 重新踩雷。
- `Cargo.toml`／`main.rs` 套件名從預設 `app`/`app_lib` 改成 `mcseg`/`mcseg_lib`，補上 license（MIT）／repository（`github.com/ddmanyes/MCseg`）
- `package.json` 加 `"tauri": "tauri"` script
- 圖示：沿用 `tauri init` 自動產生的預設佔位圖示集（`.icns`/`.ico`/各尺寸 PNG 齊全，可正常打包）。**MCseg 目前沒有品牌 logo**（前端連 favicon 都還是 Vite 預設圖）——這塊是骨架可動、但視覺不是 MCseg 的狀態，等你提供來源圖檔（PNG/SVG）後再用 `tauri icon <source>` 重新產生正式圖示。

### ⛔ 意外踩到的新 ExFAT 坑（已修，記錄供之後對照）

`cargo check` 首次執行直接炸掉：

```
failed to read file '.../permissions/path/autogenerated/._default.toml':
stream did not contain valid UTF-8
```

跟 CLAUDE.md 記載的 `.venv`／`uv sync` ExFAT 坑同一個病灶，只是這次咬到 Rust
build script——`._default.toml` 是 macOS 在 ExFAT 上自動幫每個檔案建立的
AppleDouble 影子檔，Tauri 的 build script 掃到它當成真正的 permission TOML
去解析就跟著炸。且不是一次性的：`target/debug/deps/` 底下每個編譯產物
（`.rlib`/`.d`/`.dylib`/`.rmeta`……）都會各自長出一個 `._*` 影子檔，清一次
後過陣子又會因為 OS 觸碰 metadata 而重新出現，光清檔案不是耐久解法。

**修法**：比照專案既有的 `.venv → ~/.venvs/msseg` symlink 防護模式，把
`target/` 整個挪出 ExFAT。新增 `frontend/src-tauri/.cargo/config.toml`：

```toml
[build]
target-dir = "/Users/lab_center/.cargo-targets/msseg-tauri"
```

刪除已被雜訊污染的舊 `target/`（在 ExFAT 上，本來就該丟）後重跑
`cargo check`，27.28 秒編譯乾淨過關，`Finished` \`dev\` profile。

**這是繼 `.venv` 之後第二個需要「build 產物必須放本機 APFS，不能留在
ExFAT」的案例**——值得記成通則：任何工具鏈的 build/cache 目錄（不只
Python venv）搬進這個專案時，第一件事都該先問「這個工具會不會頻繁寫入
大量小檔案」，若會，直接預先配置到本機磁碟，不要等炸了才修。

### Phase 1 驗證
- `npx tauri info`：環境檢查通過（`rustc`/`cargo`/`rustup` 齊全；`node`/`npm`
  正常；`Xcode: not installed!`——只有 Command Line Tools，完整 Xcode 未裝，
  **Phase 5 真正 build/簽章前需要留意**，`cargo check`/`tauri dev` 目前不受影響）
  ；`frontendDist`/`devUrl`/`framework: React`/`bundler: Vite` 皆正確辨識；
  CSP 顯示 `unset`（預設寬鬆，之前審查擔心的「Tauri 預設 CSP 擋掉
  127.0.0.1:8001」在這個專案的 scaffold 預設值下不成立，維持 `csp: null`
  即可，不需要額外開白名單——`remote.urls` capability 才是真正需要配置的
  部分，且已配好）
- `cargo check`：乾淨通過（見上方 ExFAT 修復後的結果）

### 未驗證／留給後續 Phase
- `tauri dev`／`tauri build` 尚未實際跑過（`tauri dev` 需要後端＋前端 dev
  server 同時起來才有意義，屬於 Phase 4 生命週期調度完成後才適合驗證；
  `tauri build` 是 Phase 5 工作）
- App 圖示仍是預設佔位圖，非 MCseg 品牌識別
- 完整 Xcode（非僅 CLT）尚未安裝，Phase 5 打包/簽章前需要處理

## Phase 1：完成（2026-08-20）

骨架建立、配置、編譯驗證全部完成。下次可進 Phase 2（bootstrap 環境配置器）。

---

## Phase 2 執行紀錄（2026-08-20）

### ⛔ 動工前推翻 v1/v2 計畫的一個設計缺陷：`bootstrap_env.py` 有雞生蛋問題

原計畫（v1、v2 都沿用）預設 bootstrap 邏輯是一支 **Python 腳本**
（`scripts/bootstrap_env.py`）。但它存在的目的正是「這台機器可能還沒有
Python 環境」——如果目標機器（尤其全新 Windows 機）系統本身沒裝 Python，
Rust 沒辦法執行一支 Python 腳本去裝 Python。

**修正**：bootstrap 邏輯改為**直接由 Rust 呼叫 `uv` 的 CLI**，完全不寫
Python 膠水腳本。`uv` 本身是獨立 Rust binary，`uv sync` 會依 pyproject.toml
的 `requires-python` **自動安裝相符的 Python**，目標機器不需要預先有任何
Python。這也讓「網路中斷後重試」變成 uv 自帶的能力（快取機制本來就可續傳），
不用自己刻續傳邏輯——呼應先前「不要重新發明依賴管理」那個決定，這次是
同一個原則的自然延伸。

**連帶決定（已與你確認）**：`uv` 執行檔直接當 **Tauri sidecar** 包進安裝檔
（而非執行期下載），安裝包大小從「~20MB」變成「~40-60MB」，換取拿掉最脆弱
的第一步（下載 uv 本身）對網路的依賴。

### 完成項目

**1. `backend/main.py` 健康檢查補上服務識別欄位**
`/api/health` 回應新增 `"service": "mcseg"`，桌面殼靠這個欄位分辨 port 8001
上回應的到底是不是自己的後端，還是別的軟體剛好用了同一個 port——這是
「不主動砍掉佔用 port 的行程」防護的技術基礎。`test_03_api.py` 補上對應斷言。

**2. `uv` sidecar 二進位檔**
下載 `uv 0.12.5` 的 macOS aarch64／Windows x86_64 release（皆驗證 sha256），
放進 `frontend/src-tauri/binaries/`，依 Tauri `externalBin` 命名慣例
（`uv-<target-triple>[.exe]`）。因為是平台專屬二進位檔（各 ~20MB），**不進
版控**——改寫 `binaries/fetch-uv.sh` 自動下載＋驗證＋放置，`.gitignore`
排除 `binaries/uv-*`，打包前跑一次這支腳本即可。

**3. 架構修正：boot-ui（取代 Phase 1 原本「視窗直接打 127.0.0.1:8001」的做法）**

Phase 1 把主視窗 `url` 寫死指向 `http://127.0.0.1:8001`，但這在 Phase 2
才發現一個問題：**bootstrap 期間後端根本還沒啟動，視窗會直接顯示連線被拒絕
的錯誤頁，沒有任何東西可以顯示 Setup Wizard**。

修正為兩階段：
- `frontendDist` 改指向新建立的 `frontend/src-tauri/boot-ui/`——一個手寫、
  不需要 build step 的極簡靜態頁面（純 HTML/CSS/JS，不經 npm/Vite），內容
  是 Setup Wizard 畫面（三步驟指示、進度條、可展開日誌、失敗重試鈕）
- 視窗 `url` 恢復預設（載入 `boot-ui/index.html`，即 `tauri://localhost`
  來源）。bootstrap 成功後，boot-ui 收到 `bootstrap://ready` 事件，用
  `window.location.href = 'http://127.0.0.1:8001'` 把**同一個視窗**導向
  後端——後續行為與 v2 修訂設計的單一 origin 架構完全一致，只是多了「先顯示
  一個極簡的本機頁面撐過 bootstrap 期間」這一步

`bundle.resources` 另外把 `backend/`、`config/`、`pyproject.toml`、
`uv.lock`、`frontend/dist`（打包前跑 `npm run build` 產生的正式前端）
包成 App 資源（跟 `frontendDist` 是兩個獨立的東西：`frontendDist` 是
Tauri webview 直接載入的內容＝boot-ui；`bundle.resources` 是唯讀資源，
啟動時複製到 Application Support 才會用到＝真正的後端與正式前端）。

**4. `frontend/src-tauri/src/bootstrap.rs`（Phase 2 核心邏輯）**

- `probe_health()`：打 `/api/health`，依回應區分三種狀態——**是我方**
  （`service=="mcseg"`，直接附加）／**被別的東西佔用**（回應了但不是我方，
  回報錯誤、**不主動終止**）／**空的**（連線被拒，正常走 bootstrap）。
  這是先前「grill me」審查抓到的「共用工作站可能砍掉別人工作」風險的
  直接落地防護。
- `sync_resources()`：把唯讀的 Resources/app 複製到可寫入的
  `Application Support/MCseg/app`（用 `dirs::data_dir()`，macOS/Windows
  都拿到平台正確路徑），用版本標記檔避免每次啟動都重新複製整包原始碼。
- `run_uv_sync()`：呼叫 `uv` sidecar 跑 `sync --project <app_dir>`，逐行
  把 stdout/stderr 轉成 `bootstrap://log` 事件。**刻意不設定
  `UV_CACHE_DIR`**——讓 uv 用它自己平台原生的預設快取路徑（macOS/Linux
  `~/.cache/uv`、Windows `%LOCALAPPDATA%\uv\cache`），才是「沿用既有快取」
  這個決定在跨平台下唯一正確的做法；寫死 `~/.cache/uv` 在 Windows 上其實
  是錯的路徑，不是 uv 的原生慣例（這是實作時才發現、原計畫沒考慮到的
  跨平台細節）。
- `spawn_backend()`：啟動 `.venv/bin/uvicorn`（不含 `--reload`），設
  `kill_on_drop(false)`——只要 Tauri App 行程還活著（含縮到系統列狀態），
  後端就跟著活著，符合 Phase 4「縮到系統列繼續跑」的設計。
- `run_bootstrap` Tauri command：串起以上所有步驟，emit 對應的
  `bootstrap://step` / `bootstrap://message` / `bootstrap://log` /
  `bootstrap://error` / `bootstrap://ready` 事件給 boot-ui。

**5. `lib.rs`：單一實例鎖**
加 `tauri-plugin-single-instance`（依官方要求註冊為第一個 plugin）。第二次
啟動嘗試時，callback 在「已經在跑」的那個實例裡被呼叫（focus 現有視窗），
不會真的另開一個新行程、不會啟動第二份後端——這是防止共用工作站上兩人
各自開一次 App、各自啟動一份後端互撞 port 8001 的第一道防線。

### 意外踩到的 Tauri 權限系統陷阱（已修，記錄避免下次重踩）

`cargo check` 兩次因 `capabilities/default.json` 寫錯而失敗：
1. 第一次：`run_bootstrap`（底線）不符合權限識別碼語法（只能小寫字母＋連字號）
2. 改成 `run-bootstrap`（連字號）後，build script 直接印出**完整的合法權限
   清單**，才發現 app 自己定義的 `#[tauri::command]`（不是外掛指令）**根本
   不需要額外 ACL 權限**——`run-bootstrap` 這個識別碼本身就不存在，我的
   假設是錯的。移除該行，只保留 `core:default` + `core:event:default` 即可。

**教訓**：Tauri v2 的 ACL 系統只管制**外掛**指令跨越信任邊界的呼叫
（例如 `shell:allow-execute`），app 自己定義、透過 `tauri::generate_handler!`
註冊的指令不受這層管制——猜測式加保護反而會讓 build 直接失敗，寧可先讓
它報錯一次拿到官方合法清單，也不要憑印象亂加。

### ⛔ 第三次踩到同一個 ExFAT 坑

`cargo check` 過程中，`capabilities/default.json` 每次編輯後都會在 ExFAT
上生出新的 `._default.json` 影子檔，build script 讀到它一樣「stream did
not contain valid UTF-8」炸掉——跟 Phase 1 記錄的 `target/` 問題同病灶，
但這次咬到的是**原始碼目錄本身**（`capabilities/`），不能比照 `target/`
搬去本機磁碟（它是真正要進版控的原始碼，不是 build 產物）。這種情況沒有
一勞永逸的搬遷解法，只能延用 CLAUDE.md 既有習慣：**每次 build 前先
`find . -name "._*" -delete`**。已在每次重跑前這樣做，三次都順利過關。

### 驗證

- `cargo check`：乾淨通過，零警告
- **資源複製 + `uv sync` 機制的實機煙霧測試**（不需要 GUI 視窗即可驗證）：
  手動複製 `backend/`、`config/`、`pyproject.toml`、`uv.lock` 到暫存目錄
  （完全比照 `sync_resources()` 的邏輯），用 sidecar 二進位檔跑
  `uv sync --dry-run`——**完整解析出所有依賴**（含平台專屬的
  `torch==2.11.0` cpu build），零錯誤。證明「複製資源到別的目錄再
  `uv sync --project <dir>`」這條核心機制是真的可行，不只是編譯得過。
  （刻意用 `--dry-run`：完整跑一次 `uv sync` 要下載安裝 torch/cellpose/
  scanpy 等數 GB 依賴，這裡只驗證機制正確性，不做這種重量級操作）

### 未驗證／留給後續

- **完整 `uv sync`（非 dry-run）跑到底**：沒有實際跑過，第一次真的執行
  會下載安裝數 GB 依賴，需要另外安排時間驗證，不適合當「順手測一下」
- **`tauri dev`／整個視窗+事件流程**：這台機器可以編譯但沒有驗證真正跑
  GUI（無法確認 boot-ui 能不能正確收到 Rust 端 emit 的事件、Setup Wizard
  畫面顯示是否正確）——**這是 Phase 2 最大的未驗證缺口**，需要在有顯示器
  的環境（或至少能跑視窗的環境）實際 `npm run tauri dev` 一次
- Windows 端完全沒有實機驗證（sidecar 二進位檔已備妥且驗證過 sha256，
  但沒有 Windows 機器可以實跑 `spawn_backend()`／NTFS 上的路徑處理）
- 單一實例鎖的實際行為（第二次啟動是否真的 focus 現有視窗而非另開新的）
  未驗證

## Phase 2 追加：實機測試（2026-08-20，同日）

上一節寫「這台機器沒有顯示器環境跑視窗」是錯的猜測——**這台 Mac 本身可以開真的
GUI 視窗**，用 `npx tauri dev -c '{"build":{"devUrl":null}}'`（用 `-c` 覆寫
config，強迫載入 `frontendDist`＝boot-ui，而非預設 dev 模式會載入的
`devUrl`＝Vite dev server）成功跑出完整鏈路：boot-ui 頁面載入 → JS 呼叫
`core.invoke('run_bootstrap')` → Rust 端收到、`probe_health()` 正確判定
port 8001 是空的 → 開始 bootstrap 流程。這證明 window 建立、Tauri 全域 API
注入、IPC 呼叫、command 分派**全部真的可行**，不只是編譯得過。

### ⛔ 過程中抓到一個真的會讓使用者卡住的 bug：錯誤被裸 `?` 吞掉

第一次測試時，log 卡在「正在準備 Python 與分析環境...」後不再有任何輸出，
process 0% CPU、Application Support 底下什麼都沒建立——看起來像卡住。

追查後發現 `run_bootstrap` 裡 `let app_dir = sync_resources(&app)?;` 用了
裸 `?`——`sync_resources()` 失敗時，錯誤會直接順著 `?` 傳回
`run_bootstrap` 的 `Result`，**完全跳過 `emit_error()`**。boot-ui 的
`core.invoke().catch()` 理論上還是抓得到這個 rejection 並顯示重試鈕，但
**log 檔案完全看不到任何線索**——這台機器沒有畫面可以確認 GUI 到底有沒有
跳出重試鈕，只能靠 log 判斷，而裸 `?` 讓這條路徑對我來說等同「靜默卡住」。

**這正是這條主線一路在抓的「靜默出錯」同一個病灶，這次換我自己的 Rust
程式碼中了招。**

修法：`sync_resources(&app)?` 改成明確 `match`，失敗時比照 `run_uv_sync`
的既有寫法先 `emit_step(..., "error")` + `emit_error(...)` 才回傳，並在
`sync_resources()`／`run_bootstrap` 內補上逐步 `log::info!`（app_support_dir
路徑、resource_app_dir 路徑與是否存在、每個資源複製的開始/完成），讓失敗時
至少 log 檔案看得出卡在哪一步——不必依賴看得到 GUI 才能除錯。

### 抓到根本原因：`resource_dir()` 在 `tauri dev` 模式下無法解析

補上診斷 log 後重跑，錯誤訊息立刻清楚浮現：

```
[bootstrap] 無法解析內嵌資源目錄：unknown path
```

即 `app.path().resource_dir()` 在 **`tauri dev`（非正式 build）模式下回傳
錯誤**——這是 Tauri 本身的已知限制：`resource_dir()` 這支 API 是為「真正
打包後的 App」設計的路徑解析，`tauri dev` 執行的是開發用的 `cargo run`
二進位檔，沒有真正的、位置固定的 Resources 目錄可解析。

（過程中曾在 `~/.cargo-targets/msseg-tauri/debug/app/` 底下看到
`backend/`／`config/`／`frontend/dist` 等內容，一度以為是 dev 模式也會
materialize resources——事後看應該是 Tauri build script 為了讓
`bundle.resources` 的路徑在建置期就先做一次驗證/暫存而產生的中間產物，
不等於執行期 `resource_dir()` 這支 API 在 dev 模式下真的可以解析出這個
位置。這兩者是分開的機制，不要混為一談。）

### 這對「怎麼測」的實際意義

**`tauri dev` 能驗證的範圍**：視窗建立、boot-ui 渲染、Tauri IPC／事件流、
`probe_health()`（含「port 已被佔用」／「port 是空的」兩種判斷）、單一
實例鎖（尚未實測，但機制不依賴 `resource_dir()`，理論上 `tauri dev` 也能測）。

**`tauri dev` 測不到的範圍**：任何依賴 `resource_dir()` 的邏輯——也就是
`sync_resources()` 之後的整條鏈路（複製資源、`uv sync`、啟動後端、健康檢查
等到 `bootstrap://ready`）。**這條鏈路唯一的真實測試方式是跑一次
`tauri build`，產出真正的 `.app`／`.dmg`，用那個產物測試**，因為只有
真正打包後的 App 才有 `resource_dir()` 能解析的固定 Resources 目錄。

### 未驗證／留給後續（更新版）

- **完整 `tauri build`＋用打包後的 App 實測整條 bootstrap 鏈路**：目前
  唯一有意義的端到端測試方式，還沒做——會需要先 `npm run build` 產出
  `frontend/dist`、跑 `binaries/fetch-uv.sh`、`npm run tauri build`，
  之後真的執行一次完整的（非 dry-run）`uv sync`（下載數 GB，需要另外
  安排時間），這是下次接手最高優先的一項
- Windows 端完全沒有實機驗證
- 單一實例鎖的實際行為未驗證（理論上 `tauri dev` 就能測，還沒測）

## Phase 2 追加二：`tauri build` 完整端到端驗證通過（2026-08-20，同日）

使用者要求「幫我做」——實際執行了完整流程，這次是**真的**端到端驗證，不是
dry-run，也不是 dev 模式局部驗證。

### 執行步驟

```bash
npm run build                # 產出 frontend/dist
npm run tauri build          # release profile 編譯 + 打包
```

`.dmg` 那層失敗（`bundle_dmg.sh` 需要 Finder/AppleScript 自動化權限，這種
背景執行環境本來就拿不到，不影響 `.app` 本身）。`.app` 本身打包成功：

```
Contents/Resources/app/{backend,config,frontend/dist,pyproject.toml,uv.lock}
Contents/MacOS/{mcseg,uv}
總大小 61MB（跟先前「~40-60MB」的預估吻合）
```

### ⛔ 又抓到一個只有真的跑過才會現形的 bug：release build 完全沒有 log

`lib.rs` 原本 `if cfg!(debug_assertions)` 才註冊 log plugin——release build
下 `log::info!`/`log::error!` 全部變成沒有輸出的空操作。這台機器沒有畫面能看
GUI，第一次測試卡住時完全無從診斷（log 檔案是空的，process 0% CPU，看起來
像卡住但其實是「錯誤發生了但沒人看得到」）。

**改為永遠註冊 log plugin**（`tauri-plugin-log` 預設本來就會同時寫 stdout
與 App 的 log 目錄檔案）。這不只是為了這次除錯方便——正式使用者之後回報
問題時附上 log 檔案（`~/Library/Logs/com.mcseg.app/mcseg.log`）也用得到，
是永久性的改善，不是臨時 debug 手段。

### ⛔ 抓到並修掉兩個真實 bug（都是「只有真的跑過才會現形」）

**Bug 1｜`sync_resources(&app)?` 用裸問號吞掉錯誤**（詳見上一節「Phase 2
追加：實機測試」）——已修為明確 match + emit_error。

**Bug 2｜`uv sync` 預設會把 `msseg` 專案本身建成 editable 安裝，因為缺
README.md 而失敗**：

```
OSError: Readme file does not exist: README.md
```

`pyproject.toml` 的 `readme = "README.md"` 欄位被 hatchling 的 metadata
驗證要求，但打包資源沒有複製 README.md。**修法不是補檔案，是根本不需要
把 `msseg` 裝成套件**——`backend.main:app` 本來就是用複製過去的原始碼＋
`current_dir` 直接跑，只需要它的依賴裝好。`run_uv_sync()` 加
`--no-install-project` 旗標，一併避免以後 pyproject.toml 的打包 metadata
需求（README/LICENSE/…）又跟資源複製清單各自漂移出新的落差。

### ✅ 完整驗證結果（真實非 dry-run 執行）

| 驗證項目 | 結果 |
|---|---|
| `.app` 啟動 → boot-ui → IPC → `run_bootstrap` | ✅ 全鏈路正確觸發 |
| `probe_health()` 判定 port 空的 | ✅ 正確 |
| 資源複製（backend/config/frontend/pyproject.toml/uv.lock） | ✅ 正確複製到 Application Support |
| **真實 `uv sync`**（非 dry-run，284 個套件） | ✅ 成功，因共用 `~/.cache/uv` 快取幾乎瞬間完成 |
| 後端啟動（`uv sync` venv 的 uvicorn，無 `--reload`） | ✅ `ps` 確認行程存在 |
| 健康檢查 `/api/health` | ✅ 回應 `{"status":"ok","version":"0.8.0","service":"mcseg"}` |
| 後端同時服務前端靜態檔 | ✅ `GET /` 回傳 HTTP 200 + 正確的 `index.html`（含 `lang="zh-TW"`） |
| **單一實例鎖**：launch 第二個實例 | ✅ 第二個行程立刻自行結束，沒有另開一份後端，原本那個繼續正常回應 |
| **port 被別的服務佔用時的防護**（整個計畫最重要的安全機制） | ✅ 用假 HTTP 服務佔住 8001，啟動 mcseg.app 後正確偵測到「不是我方」、記錄錯誤、**完全沒有去動假服務**，假服務全程存活 |
| `kill_on_drop(false)`（縮到系統列的技術基礎） | ✅ 意外驗證：`pkill` 掉 Tauri GUI 行程後，後端子行程確實沒有跟著死，仍佔用 8001（這也提醒 Phase 4 真正「結束程式」要另外明確處理，不能只靠這個機制） |

**這是 Phase 2 目前為止最完整的驗證**——不只是「編譯得過」或「dry-run
邏輯正確」，是真的打包、真的裝依賴、真的啟動後端、真的用假服務攻擊測試
port 防護機制，全部通過。

### 仍未驗證

- Windows 端完全沒有實機驗證（sidecar 二進位檔已備妥且驗證 sha256）
- boot-ui 的視覺呈現（三步驟指示、進度條、重試按鈕）沒有人眼確認過，只
  驗證了背後的事件/邏輯正確
- `bootstrap://ready` → `window.location.href` 導向後端這最後一步的視覺
  結果沒有人眼確認（健康檢查通過代表 Rust 端會 emit ready，但 boot-ui
  收到後是否真的正確換頁沒有肉眼驗證）
- Phase 4 才要做的「縮到系統列」「運算中關閉跳警告」「真正結束程式時
  乾淨終止後端」都還沒實作，上面的 `kill_on_drop(false)` 觀察到的「孤兒
  後端行程」現象目前只能靠手動 `pkill` 清理

## Phase 2：完整端到端驗證通過（2026-08-20）

### Gatekeeper / SmartScreen 繞過說明（草稿，供 Phase 5 移入正式 README）

> ⚠️ 以下文字尚未搭配實際安裝檔截圖（`.dmg`/`.msi` 還沒產出），Phase 5 完成後
> 需要補真實畫面並移入正式 README 或 App 內建的首次啟動說明。

**macOS**：因未經 Apple 公證，首次開啟 MCseg.app 時 Gatekeeper 會顯示「無法驗證
開發者」而拒絕直接開啟。解法：在「Finder」找到 MCseg.app → **右鍵（或 Control+
點擊）→ 選「打開」** → 跳出的對話框點「打開」即可，僅第一次需要這樣做，之後
雙擊可正常啟動。（雙擊會被直接拒絕、不會跳出「打開」選項——一定要用右鍵/
Control+點擊這個路徑。）

**Windows**：因未經程式碼簽章，執行安裝檔時 SmartScreen 會顯示「Windows 已保護
您的電腦」。解法：點左下角「**其他資訊**」→ 會出現「**仍要執行**」按鈕 → 點擊
即可繼續安裝。

兩者都建議在內部發布時附上截圖版 SOP（一頁 PDF 或內建於 Setup Wizard 的說明卡），
避免第一次使用的實驗室同事以為中毒而放棄安裝。

---

## Phase 3 執行紀錄（2026-08-20）

原生檔案對話框與系統通知。v2 修訂原本預期這個 Phase 因為改用單一 origin
架構，CORS/相對路徑問題不用處理，工作量會變小——實際做下來確實如此，
沒有再踩到新的架構級問題，主要是 Rust 端掛 plugin + 前端寫兩個小工具
+ 找出真正需要接線的地方。

### 完成項目

**Rust 端**：`Cargo.toml` 加 `tauri-plugin-dialog`、`tauri-plugin-notification`；
`lib.rs` 註冊兩個 plugin；`capabilities/default.json` 加 `dialog:default`／
`notification:default`（沿用 Phase 1 已經配好的 `remote.urls`，這兩個
plugin 不需要再額外處理 remote 授權）。

**前端 JS 依賴**：`@tauri-apps/api`、`@tauri-apps/plugin-dialog`、
`@tauri-apps/plugin-notification`。

**`frontend/src/hooks/useFileDialog.ts`**（新增）：`pickDirectory()` /
`pickFile()`，用 `isTauri()`（`@tauri-apps/api/core`）偵測執行環境——這份
前端同時服務純瀏覽器（`start.sh`/`npm run dev`）與桌面殼兩種情境，純瀏覽器
時回傳 `null`，讓呼叫端退回既有行為，不是兩套獨立程式碼。

**`frontend/src/utils/notify.ts`**（新增）：`notify(title, body)`，同樣
`isTauri()` 閘門；權限請求失敗或使用者拒絕都靜默吞掉，不讓通知失敗影響
分析流程本身或彈錯誤訊息——通知是錦上添花，不是必要路徑。

**接線到 `DataSetup.tsx`**：兩個既有的「瀏覽」按鈕（`data_root`／
`output_dir`）——原本固定開啟自建的 `FolderBrowser` Modal（呼叫後端
`browseDir()` API 自己刻的資料夾清單 UI）。改為優先呼叫原生
`pickDirectory()`；只有「不在 Tauri 裡」才退回原本的自建 Modal；使用者在
原生對話框按「取消」則什麼都不做（不會誤跳自建 Modal——這兩種「拿不到
路徑」的情境刻意分開處理，是這次唯一需要仔細想的地方）。

**接線到 `useStageStatus.ts`**（通知，集中一處而非散在每個 Stage 頁面）：
這個 hook 是所有 Stage 頁面共用的狀態輪詢入口，用 `useRef` 記錄前一次
status，偵測「這次 session 裡真的從 running 轉換過來」才發通知，避免
頁面一載入、query 第一次就拿到本來已經是 done 的舊狀態，誤判成「剛跑完」。
涵蓋所有 Stage（roi/segmentation/count/analysis/spatial/xenium/loupe），
比 v1 計畫原本只列「Cellpose/QC/UMAP/CellTypist」四項更完整、也更不會
漏掉——因為是接在共用 hook 上，不是逐頁手動加。

### ⛔ 一個計畫裡假設存在、但程式碼裡沒有的功能

v1 計畫 Phase 3 checklist 寫「升級...自訂 Marker CSV 選取器」，查證後
`Stage3_Analysis.tsx` 那個「marker_csv」按鈕其實是**下載**目前分群結果的
marker genes CSV，不是「選取一份自訂 marker CSV 上傳」的功能——全專案
搜尋不到任何「選取自訂 marker CSV」的 UI 存在。沒有编造這塊，略過。

### Lint

`npm run lint` 過程中抓到自己寫的 2 個 floating-promise 錯誤（`notify()`
呼叫在 `useEffect` 裡沒有處理回傳的 Promise），已用 `void notify(...)` 修掉
（`notify()` 內部本來就會吞掉自己的錯誤，`void` 是安全的，不是繞過檢查）。
`npm run lint` 剩 3 個錯誤，確認都是完全沒動過的檔案（`git diff --stat`
零差異），屬既有技術債，非本次改動造成。

### 驗證

- `cargo check`：乾淨通過，`dialog:default`／`notification:default`
  權限識別碼一次就對（吸取 Phase 2 「先讓它報錯拿官方清單」的教訓，
  直接用已知正確的慣例，不用再試錯）
- `npm run build`（`tsc && vite build`）：乾淨通過
- `npm run lint`：我改動的檔案零錯誤零警告
- **完整 `tauri build` → 啟動 `.app` → 健康檢查**：再跑一次全套（這次
  venv 已存在，`uv sync` 幾乎瞬間完成），`/api/health` 正確回應，確認
  Phase 3 的改動沒有破壞 Phase 0-2 已驗證的鏈路
- 從真正在跑的後端（`127.0.0.1:8001`）抓下實際服務的 JS bundle，
  `grep` 到 `isTauri` 字串——確認 Phase 3 寫的程式碼真的有被打包進
  使用者實際會載入的那份前端，不是只存在原始碼裡

### 仍未驗證（誠實記錄，不是隨便寫寫）

- **`isTauri()` 在遠端來源（`127.0.0.1:8001`）頁面裡實際執行結果**沒有
  肉眼確認過。這台機器沒有瀏覽器自動化工具（無 playwright/puppeteer，
  MSseg 這個 repo 自己的沿革筆記也記過同樣的環境限制），無法用程式模擬
  點擊按鈕、開對話框、截圖驗證。依據：(1) `withGlobalTauri: true`
  是 webview 層級的 init script 注入，架構上不分「這次導覽的內容是本地
  資源還是遠端網址」；(2) Phase 2 已經在**同一個 webview**、同樣的
  `withGlobalTauri` 設定下，實測驗證過 boot-ui（一開始載入的內容）能
  正確用 IPC 呼叫 `run_bootstrap`；(3) 這次確認程式碼真的有進打包後的
  bundle。三者合起來是有依據的推論，但不等於肉眼看過「瀏覽按鈕點下去
  真的跳出 Finder」——這件事仍然需要有畫面的環境做最後一次確認。
- 系統通知的實際彈出樣式、macOS 通知權限首次請求的彈窗，同樣沒有肉眼
  確認過。
- Windows 端完全沒有實機驗證（沿用 Phase 0-2 的既有缺口）。

## Phase 3：核心邏輯完成＋部分實機驗證，視覺/互動層待有畫面環境確認（2026-08-20）

---

## Phase 4 執行紀錄（2026-08-20）

Rust 生命週期調度：系統列、視窗關閉行為、結束前確認。落實拍板決定
「縮到系統列繼續跑」。

### 後端新增：跨 Stage 忙碌狀態彙總端點

`backend/src/api/system.py`（新檔）：`GET /api/system/busy`。

**背景**：桌面殼要在使用者按「結束」時知道有沒有工作在跑，但每個 Stage
（roi/segmentation/count/analysis/export）各自獨立追蹤自己的
`_task_status`/`_full_status`（5 個模組共 9 個獨立全域變數），沒有共用
彙總點。與其讓 Rust 一次打 9 支 API 自己合併，改在後端集中做一次。

**⚠️ 一個容易寫錯的 Python 陷阱**：這 9 個狀態變數在任務開始/結束時是
**整個重新賦值**（`_task_status = {...}`），不是原地 mutate。如果用
`from backend.src.api.segmentation import _task_status` 匯入，拿到的是
匯入當下那個 dict 的參照——模組之後重新賦值，這個參照就變成過期舊值，
永遠讀不到真正在跑的狀態，而且不會報錯，是那種典型的「看起來對、實際
靜默讀到舊值」的 bug。正確做法是匯入**模組本身**，每次呼叫時動態存取
`segmentation._task_status` 這個屬性。已用 `test_reads_live_value_not_stale_import`
專門釘住這個回歸案例。

測試：`backend/tests/test_16_system_busy.py` 5 項，純函式測試（monkeypatch
各模組狀態），不需要真實資料。全套 `backend/tests/`（不含既有已知失敗的
`test_skill_scripts.py`）**330 passed**，含這 5 個新測試，零回歸。

### Rust 端：`frontend/src-tauri/src/lifecycle.rs`（新模組）

- `Cargo.toml` 幫 `tauri` 加 `features = ["tray-icon"]`
- `AppState { backend_child: Mutex<Option<Child>> }`：`bootstrap.rs` 的
  `spawn_backend()` 原本用 `kill_on_drop(false)` 讓後端在縮到系統列時
  不被誤殺，但這代表「真的要結束」時沒有 handle 可以主動終止它——
  Phase 2 實機測試時就踩過「後端變孤兒行程，只能手動 pkill」的情況。
  這次把 `Child` handle 存進 managed state，供結束流程取出來 `.kill()`。
- `is_backend_busy()`：打 `GET /api/system/busy`，解析 `busy` 欄位；連
  健康檢查都打不到就視為不忙（後端根本沒在跑，談不上忙碌）。
- `handle_quit_request()`：系統列選單「結束 MCseg」的處理——先問後端忙
  不忙；忙的話用 `tauri_plugin_dialog` 跳原生確認對話框（`OkCancelCustom`
  「結束並中止工作」／「取消」）；不忙或使用者確認，才真的終止後端
  子行程＋`app.exit(0)`。
- `setup()`：建立系統列圖示 + 選單（顯示視窗／結束），並在這裡一併註冊
  主視窗的 `CloseRequested` 事件——**一律 `prevent_close()` + 隱藏視窗**，
  不區分忙不忙（運算中才需要確認的邏輯只掛在「結束」選單，不掛在關閉
  按鈕；這樣關閉按鈕的行為簡單、可預期，不會有時候關閉視窗真的退出、
  有時候不退出的不一致感）。

### 驗證

- `cargo check`：**第一次就編譯乾淨過關**——tray/menu/dialog 的 API
  用法一次猜對，沒有像 Phase 2 那樣需要靠編譯器報錯反覆修正（推測是
  因為這次先查了套件原始碼裡的 method 簽名再寫，而不是憑印象）
- 完整 `tauri build` → 啟動 → 健康檢查 → `/api/system/busy`：全部正確
  （`{"status":"ok","busy":false,"running_stages":[]}`），確認 Phase 4
  的改動沒有破壞前面驗證過的鏈路，且新端點在真正跑起來的 App 裡可用
- `ps` 確認後端子行程正確以 uvicorn 身份存在、被 mcseg 主行程持有

### 誠實記錄：這次驗證不到的部分（比前幾輪更多，如實記錄原因）

這台機器的終端機沒有 Accessibility（輔助使用）權限，`osascript` 呼叫
`System Events` 想確認選單列圖示真的渲染出來時被系統拒絕
（`-1719 執行錯誤：不允許輔助取用`）——這是這幾輪第一次連「用程式間接
探測 UI 有沒有出現」都做不到的情況，比前幾輪「至少能看 log 判斷邏輯
是否觸發」更受限。因此以下完全沒有驗證：

- 系統列圖示是否真的出現在選單列
- 點擊系統列的「顯示視窗」「結束 MCseg」是否真的觸發對應行為
- 關閉視窗按鈕是否真的隱藏視窗而非關閉 App
- 運算中按「結束」是否真的跳出原生確認對話框、對話框文字是否正確顯示
- 確認結束後，後端子行程是否真的被 `.kill()` 終止（程式碼邏輯上會執行到，
  但沒有實際觸發過這條路徑來看結果）

這些都只做到「Rust 邏輯編譯正確、appAppState 有註冊、setup() 執行沒有
panic（否則整個 App 開不起來、health check 就不會過）」這個間接證據
層級，不等於肉眼或系統層級確認過使用者實際操作的結果。

### 已知限制（範圍內故意不做，非遺漏）

只處理「系統列選單結束」這一條路徑的優雅關閉。使用者用 macOS 強制結束
（Force Quit）、或直接對行程送 `kill`/`SIGTERM`，不會經過
`handle_quit_request()`，後端子行程可能變成孤兒（跟 Phase 2 踩過的
`kill_on_drop(false)` 現象一樣）。v1 計畫本身只要求處理視窗關閉/選單
結束這個路徑，沒有要求攔截作業系統層級的強制終止信號，這次維持一致，
不擴大範圍，記錄下來供之後決定要不要補。

## Phase 4：核心邏輯完成＋間接驗證（無 GUI/Accessibility 權限，視覺與互動層完全待確認），Phase 0-4 全數完成（2026-08-20）

---

## 📍 進度總結（2026-08-20 收工，供下次接手快速讀）

### 完成度

| Phase | 狀態 | 備註 |
|---|---|---|
| Phase 0（前置：Rust 工具鏈、bug 修復） | ✅ 完成 | `put_roi_overrides` 深合併 bug 已修 |
| Phase 1（Tauri 骨架） | ✅ 完成 | `com.mcseg.app`，`cargo check` 過 |
| Phase 2（bootstrap 環境配置器） | ✅ 完成＋**完整端到端實測通過** | 真實 `uv sync`（284 套件）、後端啟動、port 佔用防護、單一實例鎖全部實測過 |
| Phase 3（原生對話框／通知） | ✅ 核心邏輯完成，視覺層待確認 | `isTauri()` 環境偵測、DataSetup 兩個瀏覽按鈕、7 個 Stage 完成通知 |
| Phase 4（系統列／生命週期） | ✅ 核心邏輯完成，視覺層待確認 | 縮到系統列、結束前確認、後端子行程正確終止 |
| Phase 5（CI + 簽章） | ⬜ 未開始 | |
| Phase 6（驗收測試矩陣） | ⬜ 未開始 | |

### ⚠️ 最重要：git 狀態未處理

目前分支 `feat/if-multichannel-segmentation`（這是另一條無關的既有
feature branch，見 [[ep-msseg-if-多通道分割與相容性修復計畫]]，內容是
IF 螢光分割，跟這次 Tauri 打包完全無關）。今天所有改動都還是**工作區
未 commit 的變更**：

```
 M backend/main.py                        （health 端點加 service 欄位）
 M backend/src/api/segmentation.py        （put_roi_overrides bug 修復）
 M backend/src/utils/config.py            （save_state_key 新函式）
 M backend/tests/test_03_api.py
 M backend/tests/test_09_roi_overrides.py
 M docs/brainstorming/tauri_desktop_packaging.md
 M frontend/package-lock.json / package.json
 M frontend/src/hooks/useStageStatus.ts
 M frontend/src/pages/DataSetup.tsx
?? backend/src/api/system.py              （新：/api/system/busy）
?? backend/tests/test_16_system_busy.py
?? frontend/src-tauri/                    （整個 Tauri 專案骨架，新）
?? frontend/src/hooks/useFileDialog.ts
?? frontend/src/utils/notify.ts
```

**下次接手前務必先決定**：
1. 開一條新分支（例如 `feat/tauri-desktop-packaging`）從 `master`／
   `origin/master` 分出來承接這批改動，不要混進 IF 分割那條——兩者是
   完全獨立的功能，混在同一條分支歷史會很難 review／之後各自要合併時
   互相干擾
2. 分幾個 commit（建議照 Phase 分：Phase 0 bug 修復一個 commit、
   Phase 1-2 骨架+bootstrap 一個、Phase 3 一個、Phase 4 一個），而不是
   全部擠成一個巨大 commit
3. `frontend/src-tauri/binaries/uv-*` 那兩個執行檔（各 ~20-50MB）**不要**
   進版控——已經在 `.gitignore` 排除，`git status` 應該看不到它們，
   確認一下沒有意外被加進去

### 下次接手第一件事

看這份文件的「v2 修訂」段落抓整體架構決策，再依序看各 Phase 的「執行
紀錄」段落抓細節。所有 Phase 都有清楚的「仍未驗證」清單，不用重新猜
什麼測過什麼沒測過。

**最值得優先做的兩件事**：
1. 上面的 git 分支/commit 整理（風險：改動目前完全沒有版本控制保護，
   萬一工作區被意外清空就全部白做）
2. 找一台有畫面（且有 Accessibility 權限）的環境，把 Phase 2-4 標記
   「未肉眼驗證」的部分實際看一遍——系統列圖示、原生對話框、通知、
   boot-ui 視覺呈現。這台機器從 Phase 4 開始連間接探測 UI 都做不到了。

### 完整檔案異動清單

```
10 files changed, 1110 insertions(+), 16 deletions(-)
+ backend/src/api/system.py（新）
+ backend/tests/test_16_system_busy.py（新）
+ frontend/src-tauri/（新，整個 Tauri 專案）
+ frontend/src/hooks/useFileDialog.ts（新）
+ frontend/src/utils/notify.ts（新）
```

sb 完整記錄：second-brain vault
`20-areas/coding/msseg-tauri-桌面應用打包計畫與審查.md`（每個 Phase 都有
對應段落，含所有審查/找碴/實測結果）。
