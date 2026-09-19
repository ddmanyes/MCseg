# MCseg

### AI 代理引導的工作流程搜尋，用於免寫程式的細胞分割與空間轉錄體轉錄本歸屬

[English](README.md) · [操作指南](docs/usage.md) · [桌面建置](docs/desktop-build.md) · [MIT 授權](LICENSE)

**MCseg（Multiple Cellpose Segmentation）** 將 Visium HD 的 H&E 影像與 2 µm 表現量 bins 轉成細胞層級資料，串接 ROI 選取、分割、RNA 計數、QC、分群、細胞型別註解、空間探索與匯出。網頁介面可免寫程式操作，CLI 支援全切片批次處理。

AI 代理只用於方法開發時的候選流程搜尋。**日常分析執行固定流程，不需要外部語言模型 API、重新搜尋架構或提供 Xenium 參考資料。** 首次安裝與模型下載需要網路。

## 工作流程

![MCseg 開發與部署流程](docs/fig1_development_deployment.png)

研究者設定操作元件、參考資料、評分與執行限制，AI 代理提出並評估候選流程，最後由研究者審查保留的配置。使用者透過本機介面匯入資料、選取 ROI、分割、RNA 計數、分析與匯出。Xenium 邊界是計算得到的參考，不是人工描繪的完整細胞真值。原圖「spsam」為 `cpsam` 的拼字誤植；圖片依原像素裁切保留。

## 研究結果

| 評估 | 文章結果 | 解讀範圍 |
| --- | --- | --- |
| LUAD 固定參數，6 個開發 ROI | PQ **0.472 ± 0.072**，2Cseg **0.432 ± 0.037** | 絕對增加 0.040，約 9% 相對提升；屬開發集結果 |
| LUAD 參考引導校準 | PQ **0.554 ± 0.063** | 使用 Xenium 逐 ROI 選擇擴張策略與距離的上限分析，不代表日常部署效能 |
| CRC ENACT 人工審閱參考 | micro-F1 **0.805 vs 0.723** | 僅比較兩方法共同涵蓋的 **10,275** 個細胞；MCseg 涵蓋全部 20,991 個參考 centroid 的 65.1% |
| CRC，15 個 ROI | NED **0.727 vs 0.712**；互斥 lineage 共表現率 **0.49% vs 0.67%** | 相較 Space Ranger，sign-flip p 分別為 0.008、0.010 |
| 新鮮冷凍乳癌固定流程 | **96,876** 個細胞、FTC **0.514**、UMI 中位數 **1,502**、NED **0.519** | 未重新搜尋架構；不代表已證實優於乳癌比較方法 |

CRC 中兩者 UMI 密度接近（11.6 vs 11.7 UMIs/µm²），但 FTC 不同（0.737 vs 0.934）。**密度接近不等於捕獲量相同。** NED 反映鄰近細胞的表現分離程度，不是絕對幾何準確度，也不能單獨作為方法排名。

LUAD 參與了方法開發。ENACT 區域與 15 個 CRC ROI 不重疊，但來自同一組織切片；其 centroid 原先由 StarDist 偵測，再經人工審閱。跨組織驗證仍有限，擴張設定需要配合影像尺度、組織形態與結果檢查。

## 安裝與開始使用

### 桌面版：Windows 與 macOS

最新可下載版本與平台狀態請見[英文首頁安裝表](README.md#desktop-installation-windows-and-macos)及 [GitHub Releases](https://github.com/ddmanyes/MCseg/releases)。請使用與 Windows x64 或 Apple Silicon Mac 相符的安裝包；原始碼 ZIP 不是桌面安裝檔。

- **Windows**：執行 `mcseg_<version>_x64-setup.exe`，依精靈安裝後從開始功能表啟動 MCseg。
- **macOS Apple Silicon**：開啟 `mcseg_<version>_aarch64.dmg`，將 MCseg 拖入 Applications，再啟動。
- 首次啟動會建立 Python 環境並下載依賴與模型；保持連線，預留至少 15 GB 加上資料與結果所需空間。無須自行安裝 Python、Node.js、Rust 或 uv。
- 若系統提示未簽章／未公證，先確認安裝檔來源與 SHA-256，再依系統提供的允許開啟流程處理。

### 原始碼與 CLI

安裝 uv、Node.js/npm 後，在原生檔案系統上執行：

```bash
git clone https://github.com/ddmanyes/MCseg.git
cd MCseg
uv sync
npm --prefix frontend install
uv run uvicorn backend.main:app --host 127.0.0.1 --port 8001
```

第二個終端機在同一專案目錄執行：

```bash
npm --prefix frontend run dev
```

開啟 [http://localhost:3000](http://localhost:3000)，從 Data Setup 設定 H&E、同一切片的 Space Ranger 2 µm outputs 與輸出目錄。另行掃描的影像需要先對位。**發行預設沒有資料路徑或 ROI**，請自行選取資料。

CLI 名稱保留為 `msseg-segment`：

```bash
uv run msseg-segment --btf /path/to/image.btf \
  --tp /path/to/tissue_positions.parquet \
  --h5 /path/to/filtered_feature_bc_matrix.h5 \
  --out /path/to/output --tissue crc --cpsam
```

`--cpsam` 開啟額外三輪推論，組成七輪流程；主模型與額外輪次均使用 cpsam。可加 `--no-gpu`、`--skip-celltypist` 或 `--export-xenium`。詳見[操作指南](docs/usage.md)。

## 模型與重現性

所有推論輪次使用 **Cellpose cpsam**，差異在影像輸入、直徑、cellprob 閾值及推論設定。流程結合 CLAHE、hematoxylin 輸入、mask 整合、可選轉錄本密度補救與 Voronoi 約束擴張。詳見[實作對照說明](analysis/supplementary/Supplementary_Note_1.md)。

記錄 Git commit、依賴版本、實際權重、輸入資料／座標與有效參數。RNA 計數的 `dilation_px` 會改變歸屬幾何，也需記錄。研究數值取自文章，本次文件更新沒有重新執行科學分析。文章分析腳本與指標在 [analysis](analysis/)，開發範例在 [docs/autoResearch](docs/autoResearch/)。

LUAD 與乳癌資料來源為 10x Genomics dataset portal；CRC 為 [GSE280318](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE280318)。文章規劃將處理後 AnnData 與 masks 存入 Zenodo，目前此處尚無正式 deposit identifier。

## 引用

> Chan, C.-R., Chang, N.-W., Wang, C.-Y., Tan, H.-Y., and Lin, S.-J. (2026). **MCseg: AI agent-guided workflow search for no-code cell segmentation and transcript attribution in spatial transcriptomics.** Manuscript.

Chan 與 Chang 為共同第一作者。規劃先發表於 bioRxiv；正式上線並取得 DOI 後再更新為預印本引用。尚未將文章標示為已發表或期刊審稿中。

## 支援與授權

問題請附版本、作業系統、參數與相關 log 至 [GitHub Issues](https://github.com/ddmanyes/MCseg/issues)。程式採 [MIT License](LICENSE)；第三方元件資訊見 [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md)。
