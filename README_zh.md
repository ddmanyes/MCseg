# MCseg

[English](README.md) | **繁體中文**

### AI 代理引導的工作流程搜尋，用於免寫程式的細胞分割與空間轉錄體轉錄本歸屬

[![bioRxiv preprint](https://img.shields.io/badge/bioRxiv-2026.09.20.752837-b31b1b)](https://doi.org/10.64898/2026.09.20.752837)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](pyproject.toml)

**MCseg（Multiple Cellpose Segmentation）** 協助空間轉錄體研究者透過本機免寫程式介面，將 **Visium HD 的 H&E 影像與 2 µm 空間表現量網格（bins）轉成細胞層級的空間轉錄體資料**。平台以同一個網頁介面串接感興趣區域（ROI）選取、細胞分割、轉錄本歸屬、品質控制、分群、細胞型別註解、空間視覺化與結果匯出；也提供 CLI 進行全切片批次處理。使用者可檢視分割結果，依組織特性調整參數。方法設計、效能評估與適用限制請參閱[bioRxiv 預印本](https://doi.org/10.64898/2026.09.20.752837)。

AI 代理協助方法開發階段的候選流程搜尋。**日常分析在本機執行開發後保留的流程，不需要重新執行 AI 代理搜尋、使用外部語言模型 API，或提供 Xenium 參考資料。** 首次安裝與模型下載需要網路。

**[閱讀論文](https://doi.org/10.64898/2026.09.20.752837) · [下載與安裝](#桌面安裝windows-與-macos) · [操作指南（英文）](docs/usage.md) · [引用 MCseg](#引用)**

[快速開始](#快速開始) · [工作流程](#工作流程) · [命令列](#命令列操作) · [重現性](#重現性)

**首次分析前請注意：**Windows 0.2.2 已包含 H&E 修正，但尚未在真實組織上完成端到端驗證。macOS 0.2.1 安裝包不含此修正，macOS 使用者請[從目前原始碼安裝](#從原始碼安裝)。選擇下載版本前，請先閱讀[發行狀態](#桌面安裝windows-與-macos)。

## 工作流程

<p align="center">
  <img src="docs/fig1_development_deployment.png" width="1000" alt="MCseg 開發與部署：研究者定義操作元件，由 AI 代理透過 Xenium 參考評分搜尋並保留分割流程；使用者在本機介面匯入影像、分割、RNA 計數、分析與匯出。右側顯示 H&E、MCseg 分割遮罩與 Xenium 參考邊界。">
</p>

**方法開發與日常部署分為兩個階段。** (a) 研究者定義操作元件、參考資料、目標與執行限制，由代理提出並評估候選流程，再由研究者審查保留的流程。(b) 使用者透過本機平台分析自己的影像與空間表現量資料。(c) 代表性的 H&E 影像、MCseg 分割遮罩與 Xenium 參考邊界。Xenium 邊界是計算得到的參考，不是人工描繪的完整細胞真值。

<details>
<summary>圖片來源與模型名稱</summary>

此圖由文章原始圖稿直接裁切，未重新取樣或更改面板內容。原圖中的「spsam」為 `cpsam` 的拼字誤植，模型說明見[模型實作](#模型實作)。裁切來源與範圍記錄於 [figure-source.md](docs/figure-source.md)。

</details>

| 步驟 | 操作 | 主要產出 |
| --- | --- | --- |
| 設定資料 | 選取 H&E 與對應的 Space Ranger 輸出 | 已驗證的輸入路徑 |
| 選取區域 | 繪製 ROI，或使用全切片 CLI | 影像裁切與空間座標 |
| 分割細胞 | 執行多輪整合分割與受約束擴張 | 細胞分割遮罩（masks） |
| 歸屬表現量 | 將空間 bins 分配至遮罩並加總計數 | 細胞 × 基因 AnnData 矩陣 |
| 分析 | 篩選細胞、計算 PCA／UMAP／Leiden，並以 CellTypist 註解 | 細胞特徵與型別標籤 |
| 探索與匯出 | 檢視空間表現量並匯出結果 | AnnData、Xenium Explorer 或 Loupe Browser 輸出 |

詳細操作請見[介面導覽](docs/usage.md#interface-tour)、[逐步操作](docs/usage.md#usage-guide)與[輸出結構](docs/usage.md#output-structure)（英文）。

## 快速開始

### 輸入資料與系統需求

請準備 H&E 影像及**與影像相符**的 Space Ranger Visium HD 輸出：

- H&E 影像，通常為分塊 BigTIFF（`.btf`、`.tif` 或 `.tiff`）。
- **2 µm bins** 的 `tissue_positions.parquet` 與 `filtered_feature_bc_matrix.h5`。
- 所選流程需要的空間中繼資料。另行掃描的影像必須先完成對位，再進行轉錄本歸屬；見[影像格式與對位](docs/usage.md#supported-image-formats)。

從原始碼安裝需要 **Python ≥3.10**、**uv**，以及網頁介面使用的 **Node.js/npm**；只使用 CLI 則不需要 Node.js。建議準備至少 16 GB RAM，大型影像需要更多記憶體；實際記憶體與磁碟需求取決於資料及安裝的依賴。文章分析使用 Apple Silicon 的 MPS 加速。系統支援 CPU，但未系統性測試 CPU 執行時間；文章中使用 MPS 的幾何評估約需每個 ROI 20–40 分鐘。

### 桌面安裝：Windows 與 macOS

桌面安裝包內含 MCseg 介面、後端程式與 `uv` 環境管理工具，**不需要自行安裝 Python、Node.js、Rust 或 uv**。首次啟動會下載並建立 Python／PyTorch／Cellpose 環境。這不是離線安裝包：請保持網路連線，預留至少 15 GB 可用磁碟空間供初始化使用，並額外保留資料與分析結果的空間。

| 平台 | 安裝方式 | 狀態 |
| --- | --- | --- |
| Windows 10/11，Intel／AMD x64 | [下載安裝檔 0.2.2](https://github.com/ddmanyes/MCseg/releases/download/desktop-v0.2.2/mcseg_0.2.2_x64-setup.exe) | 包含 H&E 修正；已在 Windows 11 驗證安裝與首次啟動 |
| macOS，Apple Silicon | [從目前原始碼安裝](#從原始碼安裝) | 已發布的 0.2.1 DMG 不含 H&E 修正，目前沒有含修正的 macOS 安裝包 |

**驗證範圍：**Windows 0.2.2 已通過建置、安裝包內容、安裝與首次啟動檢查，但此次發行驗證未執行真實組織分析。請先以小 ROI 檢查分割遮罩，再進行較大範圍的分析。詳見[發行資訊](https://github.com/ddmanyes/MCseg/releases/tag/desktop-v0.2.2)、[SHA-256 校驗檔](https://github.com/ddmanyes/MCseg/releases/download/desktop-v0.2.2/mcseg_0.2.2_x64-setup.exe.sha256)與[建置來源紀錄](https://github.com/ddmanyes/MCseg/releases/download/desktop-v0.2.2/mcseg_0.2.2_x64-setup.exe.build.json)。

H&E 修正會改變分割輸入；上述發行檢查並未重新驗證既有論文結果，比較結果前請閱讀[重現性](#重現性)。技術細節見[解卷積修正紀錄](docs/color_deconvolution_fix.md)。

桌面版與原始碼／Python 套件（`0.8.0`）使用不同版本編號。GitHub 的「Source code」壓縮檔不是桌面安裝包。發行版本資訊核對日期：2026 年 10 月 2 日。

#### Windows

1. 雙擊 **`mcseg_0.2.2_x64-setup.exe`**，依安裝精靈完成安裝。
2. 若 Microsoft Defender SmartScreen 顯示無法辨識的應用程式，先確認安裝檔來自 MCseg 維護者，再於系統政策允許時選擇 **其他資訊 → 仍要執行**。
3. 從開始功能表啟動 **MCseg**。環境建立與分析引擎啟動期間，請保持初始化視窗開啟。
4. 初始化完成後會開啟主介面。選取自己的資料，再依[操作指南](docs/usage.md#usage-guide)開始分析。

#### macOS（Apple Silicon）

請依下方[原始碼安裝步驟](#從原始碼安裝)使用含修正的流程。

<details>
<summary>歷史桌面安裝包</summary>

[macOS 0.2.1 DMG](https://github.com/ddmanyes/MCseg/releases/download/desktop-v0.2.1/mcseg_0.2.1_aarch64.dmg) 是由 `f73665f` 建置的 Apple Silicon 預發行版，**不含 H&E 修正**。此處保留連結供版本追溯，不建議用於新分析。未列出 Intel Mac 安裝包。

[v0.8.0](https://github.com/ddmanyes/MCseg/releases/tag/v0.8.0) 底下的 0.2.0 安裝包也早於修正，且包含開發分析狀態，已被取代。

</details>

#### 首次啟動疑難排解

| 狀況 | 檢查方式 |
| --- | --- |
| 初始化看似停住 | 展開初始化紀錄，檢查網路、磁碟空間與套件是否仍在下載。先排除紀錄中的錯誤，再選擇 **Retry**。 |
| 連接埠 8001 已被占用 | 關閉自己啟動的其他 MCseg／後端工作階段；若為其他程式占用，先確認程式身分再處理。 |
| 安裝包不適用於電腦 | 確認 Windows x64 與 Apple Silicon macOS 的差別；其他環境可從原始碼安裝。 |

### 從原始碼安裝

適合開發者、CLI 使用者，或沒有對應桌面安裝包的系統。

安裝 [uv](https://docs.astral.sh/uv/getting-started/installation/) 與 [Node.js](https://nodejs.org/)，將專案複製到原生檔案系統（macOS 使用 APFS、Windows 使用 NTFS，或 Linux 原生檔案系統）：

```bash
git clone https://github.com/ddmanyes/MCseg.git
cd MCseg
uv sync
npm --prefix frontend install
```

在專案根目錄啟動後端：

```bash
uv run uvicorn backend.main:app --host 127.0.0.1 --port 8001
```

於第二個終端機切換到同一專案根目錄，啟動前端：

```bash
npm --prefix frontend run dev
```

開啟 **[http://localhost:3000](http://localhost:3000)**，選取資料並依[操作指南](docs/usage.md#usage-guide)開始分析。上述指令也適用於 PowerShell。首次下載依賴與模型可能需要較長時間。

<details>
<summary>啟動腳本與外接硬碟</summary>

專案另提供 `start.sh` 與 `start.ps1`。使用前請先檢查內容：目前的 macOS shell 啟動腳本會將專案環境重建為指向 `~/.venvs/msseg` 的符號連結，並終止占用 8001／3000 連接埠的程序。在原生檔案系統上使用上述兩個終端機的啟動方式即可，不需要這些額外處理。

若資料位於 ExFAT 外接硬碟，建議將程式與 Python 環境放在系統磁碟，再讀取外接硬碟的資料。Windows 的 ExFAT 安裝可能需要設定 `UV_LINK_MODE=copy`；詳見[疑難排解](docs/usage.md#troubleshooting)。

</details>

## 命令列操作

執行指令保留既有名稱 **`msseg-segment`**。

```bash
uv run msseg-segment \
  --btf /path/to/image.btf \
  --tp /path/to/tissue_positions.parquet \
  --h5 /path/to/filtered_feature_bc_matrix.h5 \
  --out /path/to/output \
  --tissue crc \
  --cpsam
```

此指令會執行分割、bin 歸屬、細胞矩陣彙整與 CellTypist 註解。加入 `--export-xenium` 可匯出 Explorer 套件，`--skip-celltypist` 可略過註解，`--no-gpu` 則使用 CPU。若只提供 `--btf` 與 `--out`，會僅執行分割。`--cpsam` 開啟額外三輪推論，組成七輪配置；主要與額外輪次均使用 `cpsam`。

```bash
uv run msseg-segment --help
```

詳見[完整選項與 PowerShell 範例](docs/usage.md#cli-no-ui-whole-slide-pipeline)。

## 方法與設定

保留的流程結合 CLAHE 前處理、不同影像表示與直徑設定下的多輪 Cellpose `cpsam` 推論、依優先順序整合遮罩、可選的轉錄本密度補救，以及 Voronoi 約束邊界擴張。應用程式提供較精簡的配置，也能開啟額外輪次，使用七輪流程。

轉錄本歸屬會將 bin 中心點轉換到影像／遮罩座標，再將計數加總為稀疏的細胞 × 基因矩陣。正確的影像對位、像素尺度與 ROI 位移十分重要。RNA 計數階段還可透過 `rna_counting.dilation_px` 額外擴張遮罩；這會改變歸屬範圍，比較結果或重現基準時必須記錄。

設定檔位於 [`config/pipeline.yaml`](config/pipeline.yaml) 與 [`config/profiles/`](config/profiles/)。組織設定檔提供起始值，pipeline 與執行時設定可能覆寫這些數值。每份資料都應檢查最終遮罩與實際生效參數。

### 模型實作

MCseg 的多輪整合分割使用 **Cellpose `cpsam`**。各輪差異在影像表示、直徑與細胞機率閾值（cell-probability threshold），並非混合不同模型家族。七輪配置包含三輪不同直徑、一輪 hematoxylin 輸入，以及三輪額外 `cpsam` 推論。[`_load_primary_model`](backend/src/segmentation/cellpose_runner.py) 載入器會在執行紀錄中寫入實際模型權重路徑。

為便於重現，請保留分析版本、環境、實際模型權重與生效設定。

## 重現性

**論文分析與目前應用程式必須區分。** 封存的 [LUAD benchmark 腳本](analysis/scripts/analysis/08_luad_benchmark.py) 仍保留舊版 H&E 逆矩陣轉置公式，[目前應用程式](backend/src/segmentation/cellpose_runner.py) 則使用修正後的公式。[修正紀錄](docs/color_deconvolution_fix.md) 明確記載未重跑歷史分割或模型推論，因此已發表的基準數值不能視為修正版應用程式的驗證。重現論文結果需使用對應的分析程式與設定；評估修正版則需重新執行並記錄版本。

| 資源 | 內容 |
| --- | --- |
| [`analysis/scripts/`](analysis/scripts/) | 基準分析與製圖腳本 |
| [`analysis/data/`](analysis/data/) | 已納入版本控制的指標與摘要表格 |
| [`analysis/supplementary/`](analysis/supplementary/) | 補充說明與表格；使用前請核對與稿件的版本一致性 |
| [`docs/autoResearch/`](docs/autoResearch/) | 開發提示詞、執行器與起始範本 |
| [`backend/src/`](backend/src/) | 部署使用的分析實作 |

代理引導的開發迴圈參考 AutoResearch 方法：研究者選定候選操作、參考資料、評分方式、提示詞與執行限制，由代理提出並評估可執行流程，再由研究者審查保留的配置。開發範本與日常分析分開，執行範本需要另行設定 API；提供範本不代表已完整封存每次歷史搜尋紀錄。

重現分析時，請保留 Git 版本、解析後的依賴、模型權重識別資訊、輸入資料與座標、生效的分割／計數參數，以及執行紀錄。專案套件沿用歷史名稱 `msseg`。以目前原始碼、Windows 0.2.2 安裝包或清理後的 macOS 0.2.1 安裝包全新安裝時，預設輸入路徑與 ROI 選取為空，請先設定自己的資料；既有安裝可能保留先前儲存的設定。

### 資料取得

- **LUAD：**10x Genomics 的配對 [post-Xenium Visium HD](https://www.10xgenomics.com/datasets/visium-hd-cytassist-gene-expression-human-lung-cancer-post-xenium-expt) 與 [Xenium Prime 5K](https://www.10xgenomics.com/cn/datasets/xenium-human-lung-cancer-post-xenium-technote) 資料（Experiment 2）；使用六個開發 ROI。
- **CRC：**[GEO GSE280318](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE280318)；包含 15 個轉錄本評估 ROI，以及同一切片上另一個經專家審閱的 ENACT 區域。
- **乳癌：**10x Genomics 的公開[新鮮冷凍人類乳癌 Visium HD 資料](https://www.10xgenomics.com/datasets/visium-hd-cytassist-gene-expression-human-breast-cancer-fresh-frozen)，用於固定流程的跨組織應用。

文章規劃將處理後的 AnnData 與分割遮罩存入 Zenodo。目前此處尚未提供公開的 Zenodo 典藏識別碼；預印本 DOI 請見[引用](#引用)。

## 引用

方法設計、效能評估與適用限制請參閱以下 bioRxiv 預印本；若使用 MCseg，請引用此文章：

> Chan, C.-R., Chang, N.-W., Wang, C.-Y., Tan, H.-Y., and Lin, S.-J. (2026). **MCseg: AI agent-guided workflow search for no-code cell segmentation and transcript attribution in spatial transcriptomics.** *bioRxiv* [preprint]. [https://doi.org/10.64898/2026.09.20.752837](https://doi.org/10.64898/2026.09.20.752837).

Chan 與 Chang 為共同第一作者。預印本於 2026 年 9 月 25 日公開，尚未經同儕審查認證。

## 支援與授權

操作細節請見[操作指南（英文）](docs/usage.md)。如遇問題，請至 [GitHub Issues](https://github.com/ddmanyes/MCseg/issues) 提供作業系統、Git 版本、套件版本、指令／設定及相關紀錄。

MCseg 採用 [MIT License](LICENSE)。第三方元件資訊請見 [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md)，桌面建置與驗證方式見[桌面建置指南（英文）](docs/desktop-build.md)。
