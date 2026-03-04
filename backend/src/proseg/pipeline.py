"""
Proseg Pipeline Module

Integrated workflow for Proseg cell segmentation refinement using SpatialData Zarr format.
"""

import os
import subprocess
import shutil
import json
from pathlib import Path
from typing import Optional, Dict, Any, Tuple
import numpy.typing as npt
import warnings
import logging
import cv2  # Add OpenCV for fast dilation

# Configure logger
logger = logging.getLogger("pipeline.proseg")

import dask
dask.config.set({"dataframe.query-planning": True})
import dask.dataframe as dd

import numpy as np
import pandas as pd
import spatialdata as sd
from anndata import AnnData

from backend.src.proseg.zarr_handler import (
    load_zarr,
    extract_transcripts,
    extract_masks,
    extract_image,
    verify_alignment,
    get_scale_factors,
    calculate_roi,
    add_shapes_to_zarr,
)


class ProsegPipeline:
    """
    整合化的 Proseg 工作流程

    此類別封裝了從 Zarr 讀取、執行 Proseg、組裝 AnnData 到回寫 Zarr 的完整流程。

    Parameters
    ----------
    zarr_path : str
        輸入的 SpatialData Zarr 檔案路徑
    output_dir : str
        輸出目錄路徑
    proseg_bin : str, optional
        Proseg 執行檔路徑，預設為 ~/.cargo/bin/proseg
    max_dist : float, default=20
        Proseg 最大轉錄點至核心距離參數
    compactness : float, default=0.1
        Proseg 細胞緊密度參數
    dilation_radius : int, default=5
        核遮罩擴張半徑(像素)
    cyto_mask_path : str, optional
        外部細胞質遮罩檔案路徑(.npy格式)。如果提供,將用於細胞ID分配,
        取代核遮罩+擴張的方式
    use_cyto_mask_from_zarr : bool, default=False
        是否從 Zarr 檔案載入細胞質遮罩(labels/cellpose_cyto)。
        此選項優先於 cyto_mask_path
    nucleus_label_name : str, default="cellpose_nuclei"
        Zarr 檔案中細胞核標籤的名稱
    cyto_label_name : str, default="cellpose_cyto"
        Zarr 檔案中細胞質標籤的名稱

    Attributes
    ----------
    sdata : SpatialData
        載入的 SpatialData 物件
    transcripts_df : pd.DataFrame
        提取的轉錄點位
    mask : np.ndarray
        提取的細胞核遮罩
    cyto_mask : Optional[np.ndarray]
        提取的細胞質遮罩 (如果啟用)
    proseg_results : dict
        Proseg 執行結果

    Examples
    --------
    >>> pipeline = ProsegPipeline(
    ...     zarr_path="data/proseg_integrated.zarr",
    ...     output_dir="results/proseg_output"
    ... )
    >>> pipeline.run_full_pipeline()
    """

    def __init__(
        self,
        zarr_path: str,
        output_dir: str,
        proseg_bin: Optional[str] = None,
        max_dist: float = 20.0,
        compactness: float = 0.1,
        dilation_radius: int = 5,
        samples: int = 200,          # New: Control MCMC iterations
        burnin_samples: int = 200,   # New: Control burn-in phase
        coordinate_scale: float = 1.0, # New: Control pixel-to-micron scale
        padding: int = 100,            # New: Configurable padding
        force_scale: Optional[float] = None, # New: Manual override for Zarr scale
        cyto_mask_path: Optional[str] = None, # New: Cytoplasm mask path for spatial constraints
        use_cyto_mask_from_zarr: bool = False, # New: Load cyto mask from Zarr labels/cellpose_cyto
        nucleus_label_name: str = "cellpose_nuclei", # New: Configurable nucleus label name
        cyto_label_name: str = "cellpose_cyto", # New: Configurable cyto label name
        use_watershed: bool = True,           # New: Use Watershed seeding
        recorded_samples: int = 150,           # New: MCMC recorded samples
        enforce_connectivity: bool = True,      # New: Force contiguous polygons
        fixed_roi: Optional[Tuple[float, float, float, float]] = None # New: (x, y, w, h) in Global Units
    ):
        self.zarr_path = zarr_path
        self.output_dir = Path(output_dir)
        self.max_dist = max_dist
        self.compactness = compactness
        self.dilation_radius = dilation_radius
        self.samples = samples
        self.burnin_samples = burnin_samples
        self.coordinate_scale = coordinate_scale
        self.padding = padding
        self.force_scale = force_scale
        
        # 強制修正邏輯 (防禦 um/px 對齊災難)
        if self.force_scale == 1.0:
            import traceback
            caller = "".join(traceback.format_stack()[-2:-1])
            logger.warning(f"  🔍 偵測到強制 Scale 1.0 (呼叫源：{caller.strip()})")
            logger.warning("  ⚠️  自動修復：將 force_scale 設回 None 以防止對齊漂移！ (如果要保持 1.0 請改設定為 1.000001)")
            self.force_scale = None
            
        self.cyto_mask_path = cyto_mask_path
        self.use_cyto_mask_from_zarr = use_cyto_mask_from_zarr
        self.nucleus_label_name = nucleus_label_name
        self.cyto_label_name = cyto_label_name
        self.use_watershed = use_watershed
        self.recorded_samples = recorded_samples
        self.enforce_connectivity = enforce_connectivity
        self.fixed_roi = fixed_roi

        # Proseg 執行檔路徑
        if proseg_bin is None:
            default_path = os.path.expanduser("~/.cargo/bin/proseg")
            if os.path.exists(default_path):
                self.proseg_bin = default_path
            elif shutil.which("proseg"):
                self.proseg_bin = "proseg"
            else:
                raise FileNotFoundError(
                    "找不到 Proseg 執行檔。請安裝 Proseg 或指定執行檔路徑。\n"
                    "安裝指令: cargo install proseg"
                )
        else:
            self.proseg_bin = proseg_bin

        # 初始化屬性
        if self.coordinate_scale == 1.0:
            logger.warning("⚠️ coordinate_scale is set to 1.0. For most Proseg-Zarr tasks, "
                           "please ensure this matches your nm/px scale (e.g., 0.2645833).")

        self.sdata: Optional[sd.SpatialData] = None
        self.transcripts_df: Optional[pd.DataFrame] = None
        self.mask: Optional[npt.NDArray[np.int32]] = None
        self.cyto_mask: Optional[npt.NDArray[np.int32]] = None  # Store cyto mask if loaded from Zarr
        self.proseg_results: Dict[str, Any] = {}

        # ROI 資訊
        self.roi_offset = (0, 0) # (min_x, min_y) pixel coordinates
        self.roi_shape = None # (h, w)
        self.scale_factors = (None, None) # (y, x)


        # 建立輸出目錄
        self.output_dir.mkdir(parents=True, exist_ok=True)

        logger.info(f"🚀 Proseg Pipeline 初始化完成")
        logger.info(f"  - 輸入 Zarr: {self.zarr_path}")
        logger.info(f"  - 輸出目錄: {self.output_dir}")
        logger.info(f"  - Proseg 執行檔: {self.proseg_bin}")
        logger.info(f"  - 參數: max_dist={self.max_dist}, compactness={self.compactness}, dilation={self.dilation_radius}")

    def load_data(self) -> None:
        """
        載入 Zarr 資料並提取轉錄點位與遮罩
        """
        logger.info("=" * 60)
        logger.info("步驟 1: 載入 Zarr 資料")
        logger.info("=" * 60)

        # 載入 Zarr
        self.sdata = load_zarr(self.zarr_path)

        # 驗證對齊
        # 驗證對齊 (移除全域驗證以節省記憶體，改用 ROI 驗證)
        # logger.info("\n驗證座標系對齊...")
        # is_aligned, stats = verify_alignment(self.sdata)
        # if not is_aligned:
        #     warnings.warn(
        #         f"座標系可能未對齊！重疊率: {stats['overlap_rate']:.1%}",
        #         UserWarning
        #     )

        # 提取轉錄點位
        logger.info("提取轉錄點位...")
        self.transcripts_df = extract_transcripts(self.sdata, points_name="transcripts")

        # 計算 ROI (Bounding Box + Padding)
        # 計算 ROI (Bounding Box + Padding)
        logger.info("計算 ROI...")
        # 取得 Scale Factor
        # 取得 Scale Factor
        if self.force_scale is not None:
             logger.info(f"  ⚠️  使用強制 Scale Factor: {self.force_scale}")
             self.scale_factors = (self.force_scale, self.force_scale)
        else:
             self.scale_factors = get_scale_factors(self.zarr_path, "images/tissue_hires_image")

        scale_y, scale_x = self.scale_factors
        logger.debug(f"  - Coordinate Scale Factors: y={scale_y:.4f}, x={scale_x:.4f}")

        # 使用 centralized ROI 計算函式
        if self.fixed_roi:
            roi_x_px, roi_y_px, roi_w_px, roi_h_px = self.fixed_roi
            logger.info(f"  🔍 使用固定 ROI (Pixels): x={roi_x_px}, y={roi_y_px}, w={roi_w_px}, h={roi_h_px}")

            # Convert fixed_roi (pixels) to Global Units (nm) for transcript filtering
            roi_min_x_nm = roi_x_px * scale_x
            roi_max_x_nm = (roi_x_px + roi_w_px) * scale_x
            roi_min_y_nm = roi_y_px * scale_y
            roi_max_y_nm = (roi_y_px + roi_h_px) * scale_y

            self.transcripts_df = self.transcripts_df[
                (self.transcripts_df['x'] >= roi_min_x_nm) & (self.transcripts_df['x'] <= roi_max_x_nm) &
                (self.transcripts_df['y'] >= roi_min_y_nm) & (self.transcripts_df['y'] <= roi_max_y_nm)
            ].copy()
            logger.info(f"  ✅ 篩選轉錄點位完成: 剩餘 {len(self.transcripts_df)} 個點位")

            # Re-calculate or Set fixed ROI parameters
            # roi format needs to be in label space (pixels)
            self.roi_offset, roi = calculate_roi(
                self.transcripts_df,
                self.scale_factors,
                padding=self.padding
            )
        else:
            self.roi_offset, roi = calculate_roi(
                self.transcripts_df,
                self.scale_factors,
                padding=self.padding
            )
        # roi is (min_y, max_y, min_x, max_x)

        # 提取遮罩 (使用 ROI)
        logger.info(f"提取細胞核遮罩 ({self.nucleus_label_name})...")
        self.mask = extract_masks(self.sdata, label_name=self.nucleus_label_name, scale=0, roi=roi)
        self.roi_shape = self.mask.shape

        # 提取細胞質遮罩 (如果啟用)
        if self.use_cyto_mask_from_zarr:
            logger.info(f"提取細胞質遮罩 ({self.cyto_label_name}) 從 Zarr...")
            try:
                self.cyto_mask = extract_masks(self.sdata, label_name=self.cyto_label_name, scale=0, roi=roi)
                if self.cyto_mask.shape != self.mask.shape:
                    logger.warning(f"  ⚠️  細胞質遮罩尺寸 {self.cyto_mask.shape} 與核遮罩 {self.mask.shape} 不符！")
                    self.cyto_mask = None
                else:
                    logger.info(f"  ✅ 細胞質遮罩載入成功，尺寸: {self.cyto_mask.shape}")
                    logger.info(f"  - 細胞質遮罩包含 {len(np.unique(self.cyto_mask)) - 1} 個細胞")
            except (KeyError, ValueError, IOError, OSError) as e:
                logger.error(f"  ❌ 載入細胞質遮罩失敗: {e}")
                logger.warning(f"  ⚠️  將回退至核遮罩 + 擴張模式")
                self.cyto_mask = None
            except Exception as e:
                logger.critical(f"  ❌ 嚴重錯誤: {type(e).__name__}: {e}")
                logger.warning(f"  ⚠️  將回退至核遮罩 + 擴張模式")
                self.cyto_mask = None

        # 驗證 ROI 對齊
        logger.info("內部驗證 ROI 對齊狀態...")
        h_mask, w_mask = self.mask.shape
        roi_min_x, roi_min_y = self.roi_offset

        # 採樣檢查
        if len(self.transcripts_df) == 0:
            logger.warning("  ⚠️  無轉錄點位，跳過對齊驗證")
        else:
            sample_df = self.transcripts_df.sample(min(1000, len(self.transcripts_df)))
            # Global -> Local Pixel
            scale_y, scale_x = self.scale_factors
            local_x = (sample_df['x'].values / scale_x) - roi_min_x
            local_y = (sample_df['y'].values / scale_y) - roi_min_y

            # 檢查邊界
            in_bounds = (local_x >= 0) & (local_x < w_mask) & (local_y >= 0) & (local_y < h_mask)
            valid_x = local_x[in_bounds].astype(int)
            valid_y = local_y[in_bounds].astype(int)

            # 檢查重疊 (Cell ID > 0)
            hits = 0
            if len(valid_x) > 0:
                cell_ids = self.mask[valid_y, valid_x]
                hits = (cell_ids > 0).sum()

            hit_rate = hits / len(sample_df)
            logger.info(f"  - 採樣點落入細胞率: {hit_rate:.1%} ({hits}/{len(sample_df)})")
            if hit_rate < 0.1:
                logger.warning(f"ROI 區域內細胞重疊率過低 ({hit_rate:.1%})，請檢查座標對齊!")

        logger.info(f"✅ 資料載入完成")

    def prepare_proseg_input(self) -> Path:
        """
        準備 Proseg 輸入 CSV 檔案

        Returns
        -------
        Path
            CSV 檔案路徑
        """
        logger.info("=" * 60)
        logger.info("步驟 2: 準備 Proseg 輸入")
        logger.info("=" * 60)

        csv_path = self.output_dir / "transcripts_for_proseg.csv"

        # 檢查必要欄位
        required_cols = ['x', 'y', 'gene']
        missing_cols = set(required_cols) - set(self.transcripts_df.columns)
        if missing_cols:
            raise ValueError(f"轉錄點位缺少必要欄位: {missing_cols}")

        # 準備資料
        df = self.transcripts_df.copy()



        # 加入 cell_id（從遮罩查找）
        logger.info("從遮罩查找細胞 ID...")
        h_mask, w_mask = self.mask.shape

        # ============================================================
        # 準備 Lookup Mask 的正確邏輯
        # ============================================================
        # 邏輯: 核遮罩 -> 擴張 -> Cyto 約束
        # 1. 從核遮罩開始
        # 2. 如果有 dilation,對核遮罩進行擴張
        # 3. 如果有 cyto_mask,使用它約束擴張邊界
        # ============================================================

        lookup_mask = None

        # 步驟 1: 取得初始 Lookup Mask (擴張或分水嶺)
        if self.use_watershed:
            logger.info(f"使用 Watershed 分水嶺種子分配 (半徑: {self.dilation_radius}px)...")
            try:
                from scipy.ndimage import distance_transform_edt
                # 核遮罩擴張作為有效邊界
                kernel_size = 2 * self.dilation_radius + 1
                kernel = cv2.getStructuringElement(cv2.MORPH_RECT, (kernel_size, kernel_size))
                dilated_area = cv2.dilate((self.mask > 0).astype(np.float32), kernel).astype(bool)

                # 計算最近核 ID (反向距離轉換)
                inv_seeds = (self.mask == 0)
                _, (iy, ix) = distance_transform_edt(inv_seeds, return_indices=True)
                nearest_labels = self.mask[iy, ix]

                # 在有效擴張範圍內保留最近 ID
                dilated_mask = np.where(dilated_area, nearest_labels, 0).astype(np.int32)
                logger.info(f"  ✅ Watershed 分配完成")
            except Exception as e:
                logger.warning(f"  ❌ Watershed 失敗 ({e}), 回退至傳統擴張")
                kernel_size = 2 * self.dilation_radius + 1
                kernel = cv2.getStructuringElement(cv2.MORPH_RECT, (kernel_size, kernel_size))
                dilated_mask = cv2.dilate(self.mask.astype(np.float32), kernel).astype(np.int32)
        elif self.dilation_radius > 0:
            logger.info(f"對核遮罩進行傳統擴張 (半徑: {self.dilation_radius}px)...")
            kernel_size = 2 * self.dilation_radius + 1
            kernel = cv2.getStructuringElement(cv2.MORPH_RECT, (kernel_size, kernel_size))
            dilated_mask = cv2.dilate(self.mask.astype(np.float32), kernel).astype(np.int32)
        else:
            dilated_mask = self.mask
            logger.info("  - 不進行擴張,使用原始核遮罩")

        # 步驟 2: 取得 Cyto 約束 (如果有)
        cyto_constraint = None

        # 2a. 從 Zarr 載入的 cyto mask
        if self.cyto_mask is not None:
            logger.info(f"使用 Zarr 細胞質遮罩 ({self.cyto_label_name}) 作為約束...")
            cyto_constraint = self.cyto_mask

        # 2b. 從外部 .npy 檔案載入 cyto mask
        elif self.cyto_mask_path:
            logger.info(f"載入外部細胞質遮罩: {self.cyto_mask_path}")
            try:
                cyto_mask_full = np.load(self.cyto_mask_path, mmap_mode='r')
                roi_min_x, roi_min_y = self.roi_offset
                h_mask, w_mask = self.roi_shape
                cyto_mask_roi = cyto_mask_full[roi_min_y:roi_min_y+h_mask,
                                                roi_min_x:roi_min_x+w_mask].copy()

                if cyto_mask_roi.shape != (h_mask, w_mask):
                    raise ValueError(f"Cyto mask ROI 尺寸不符: {cyto_mask_roi.shape} != {(h_mask, w_mask)}")

                cyto_constraint = cyto_mask_roi
                logger.info(f"  ✅ 外部細胞質遮罩載入成功")

            except (FileNotFoundError, ValueError, IOError) as e:
                logger.error(f"  ❌ 載入失敗: {e}")
                logger.warning(f"  ⚠️  將不使用 cyto 約束")
                cyto_constraint = None
            except Exception as e:
                logger.critical(f"  ❌ 嚴重錯誤: {type(e).__name__}: {e}")
                cyto_constraint = None


        # 步驟 3: 應用 Cyto 約束 (如果有)
        if cyto_constraint is not None:
            logger.info("應用細胞質遮罩約束...")

            # 約束邏輯: 只保留在對應 cyto 區域內的擴張
            # 對於每個像素,如果擴張後有值且在對應的 cyto 內,則保留
            lookup_mask = np.where(
                (dilated_mask > 0) & (cyto_constraint == dilated_mask),
                dilated_mask,
                0
            ).astype(np.int32)

            # 統計約束效果
            dilated_cells = len(np.unique(dilated_mask)) - 1
            constrained_cells = len(np.unique(lookup_mask)) - 1
            logger.info(f"  - 擴張後細胞數: {dilated_cells}")
            logger.info(f"  - Cyto 約束後細胞數: {constrained_cells}")

            if constrained_cells < dilated_cells:
                logger.info(f"  ✅ Cyto 約束生效: {dilated_cells - constrained_cells} 個細胞被限制")
            else:
                logger.info(f"  ℹ️  擴張未超出 cyto 邊界")
        else:
            # 無 cyto 約束,使用擴張結果
            lookup_mask = dilated_mask
            logger.info("  - 無 cyto 約束,使用完整擴張結果")

        # 統計最終結果
        original_cells = len(np.unique(self.mask)) - 1
        final_cells = len(np.unique(lookup_mask)) - 1
        logger.info(f"  - 原始核遮罩細胞數: {original_cells}")
        logger.info(f"  - 最終 lookup mask 細胞數: {final_cells}")


        # 使用 NumPy 向量化索引加速查找 (比 apply 快 100x 以上)
        # 1. 取得座標並轉換為 Mask 座標 (Pixel)
        scale_y, scale_x = self.scale_factors
        # 保持 float 以進行減法，最後轉 int
        x_pixels = (df['x'].values / scale_x)
        y_pixels = (df['y'].values / scale_y)


        # 2. 建立有效索引遮罩 (確保座標在影像範圍內 - 這裡指全域影像)
        # 但實際上我們只關心 ROI 內
        roi_min_x, roi_min_y = self.roi_offset

        local_x = x_pixels - roi_min_x
        local_y = y_pixels - roi_min_y

        # 更新 valid_mask: 必須在裁剪後的 Mask 範圍內
        valid_mask = (local_x >= 0) & (local_x < w_mask) & \
                     (local_y >= 0) & (local_y < h_mask)

        # 3. 初始化 cell_id 陣列 (預設為 0/背景)
        cell_ids = np.zeros(len(df), dtype=self.mask.dtype)

        # 4. 直接索引查找
        valid_loc_x = local_x[valid_mask].astype(int)
        valid_loc_y = local_y[valid_mask].astype(int)

        # 先用可能已擴張的 lookup_mask (包含 Watershed / Cyto 邊界) 進行查找填寫
        cell_ids[valid_mask] = lookup_mask[valid_loc_y, valid_loc_x]

        # 🎯【核心突破】：強制恢復「不可被侵犯」的核心核區域
        # 為什麼？因為 lookup_mask 經過擴張/分水嶺後，可能會誤把 B 細胞的某些轉錄點「分配」給 A。
        # 如果一開始餵錯，Proseg 就會順理成章吃掉 B 的核。
        # 我們在這裡做「硬重置」：只要轉錄點落在最原始、沒有擴張過的細胞核 (self.mask) 裡，
        # 它就絕對必須屬於這個核！
        pure_nuc_ids = self.mask[valid_loc_y, valid_loc_x]
        pure_nuc_mask = (pure_nuc_ids > 0)
        
        # 只取有落在純核中的 df indices，並強制覆寫回最真實的 cell ID
        if pure_nuc_mask.sum() > 0:
            # np.where(valid_mask)[0] 能把 valid 的子集 mapping 回原本 df 的 index
            global_indices = np.where(valid_mask)[0][pure_nuc_mask]
            cell_ids[global_indices] = pure_nuc_ids[pure_nuc_mask]



        # 5. Smart Filtering (記憶體優化)
        # -----------------------------
        # 過濾掉距離任何細胞過遠的背景點，減小 CSV 大小
        # 原理: Proseg 只會考慮 max_dist 內的點，其他的點對分配無貢獻(除背景估計外)
        # 我們保留 max_dist * 1.5 範圍內的背景點以確保背景估計足夠

        logger.info("執行 Smart Filtering (過濾遠離細胞的背景點)...")

        # 計算擴張半徑 (像素)
        max_dist_um = self.max_dist
        # 若被 force_scale 覆蓋為 1.0（處理 pixel space 時），應改用真實的 self.coordinate_scale 計算物理半徑
        if self.force_scale is not None and self.force_scale == 1.0:
            scale_in_um = self.coordinate_scale
            logger.info(f"  🔍 偵測到 force_scale=1.0，改採自定義物理比例計算濾波: {scale_in_um} um/px")
        else:
            scale_avg = (scale_x + scale_y) / 2
            scale_in_um = scale_avg

        # Scale 單位判斷 (Heuristic)
        # 如果 scale > 50, 假設是 nm/px，轉為 um/px
        if scale_in_um > 50:
            logger.info(f"  🔍 偵測到 Scale 為 nm 單位 ({scale_in_um:.2f})，轉換為 um ({scale_in_um/1000:.4f})")
            scale_in_um = scale_in_um / 1000.0

        filter_radius_px = int((max_dist_um / scale_in_um) * 1.5)

        logger.info(f"  - 濾波半徑: {max_dist_um} um * 1.5 = {filter_radius_px} px")

        if filter_radius_px > 0:
            # 建立 "有效背景區域" Mask
            # 1. 取出所有 > 0 的區域 (細胞)
            cell_mask_binary = (lookup_mask > 0).astype(np.uint8)

            # 2. 擴張
            kernel_size = 2 * filter_radius_px + 1
            kernel = cv2.getStructuringElement(cv2.MORPH_RECT, (kernel_size, kernel_size))
            valid_bg_mask = cv2.dilate(cell_mask_binary, kernel, iterations=1)

            # 3. 篩選 DataFrame
            # 使用 valid_loc 檢查每個點是否在 valid_bg_mask 內
            # 注意: valid_loc 是已經在 ROI 內的局部坐標

            # 為了快速索引，我們可以用 boolean indexing
            # valid_mask (前面定義的) 標記了哪些點在 ROI 內
            # 現在我們進一步: 哪些點在 valid_bg_mask 內?

            # a. 標記在 valid_bg_mask 內的點
            # Record total before filtering
            total_before = len(df)

            # Note: Actual df filtering happens at line 493 using final_keep_mask
            # to prevent shape mismatch with valid_mask.
            # cell_ids 是對應 df 的。

            # 修正邏輯:
            # 前面 L424: cell_ids[valid_mask] = lookup_mask[...]
            # 這裡我們需要一個能篩選 df 的 boolean array

            final_keep_mask = np.zeros(len(df), dtype=bool)

            # 只有在 valid_mask (ROI內) 且 在 valid_bg_mask (擴張範圍內) 的才保留
            # 我們需要構建一個與 valid_mask 長度相同的 array 來存放 "in_bg_mask"
            subset_in_bg = valid_bg_mask[valid_loc_y, valid_loc_x] > 0

            # 將 subset 結果映射回 full mask
            # valid_mask 是 True 的位置，放入 subset_in_bg 的值
            final_keep_mask[valid_mask] = subset_in_bg

            df = df[final_keep_mask].copy()

            # 也要更新 local_x, local_y (它們是 numpy array，長度等於原 df)
            # 但我們已經把 df filter 了，所以需要重新 assign x, y
            # 或者先 filter local_x, local_y

            # 比較簡單的方法：直接用過濾後的 df 的 x, y 重算 local (或保留先前的計算)
            # 由於 Proseg 需要 local x, y 作為 columns
            # 我們在 L430 賦值，所以我們得確保 L430 賦的值是對的

            # 重新整理:
            # 1. 前面計算了 local_x, local_y (全長)
            # 2. 前面計算了 cell_ids (全長)
            # 3. 過濾

            df['cell_id'] = cell_ids[final_keep_mask]
            df['x'] = local_x[final_keep_mask]
            df['y'] = local_y[final_keep_mask]

            logger.info(f"  - 過濾後: {len(df):,} / {total_before:,} (減少 {(1 - len(df)/total_before)*100:.1f}%)")

        else:
            # 若無濾波 (半徑0)，則只做基本的 cell_id 賦值 & ROI transform
            df['cell_id'] = cell_ids
            df['x'] = local_x
            df['y'] = local_y




        # 加入其他必要欄位
        df['qv'] = 40  # 品質分數（Proseg 需要，但不影響結果）
        df['z'] = 0.0  # Z 座標（2D 資料）

        # 儲存
        df[['x', 'y', 'gene', 'qv', 'cell_id', 'z']].to_csv(csv_path, index=False)

        # 統計
        n_in_cells = (df['cell_id'] > 0).sum()
        logger.info(f"✅ Proseg 輸入準備完成")
        logger.info(f"  - 檔案: {csv_path}")
        logger.info(f"  - 總轉錄點數: {len(df):,}")
        logger.info(f"  - 落在細胞內: {n_in_cells:,} ({n_in_cells/len(df)*100:.1f}%)")
        logger.info(f"  - 背景點數: {(df['cell_id'] == 0).sum():,}")

        return csv_path

    def run_proseg(self, csv_path: Path) -> Dict[str, Path]:
        """
        執行 Proseg CLI

        Parameters
        ----------
        csv_path : Path
            輸入 CSV 路徑

        Returns
        -------
        dict
            輸出檔案路徑字典
        """
        logger.info("=" * 60)
        logger.info("步驟 3: 執行 Proseg")
        logger.info("=" * 60)

        # 定義輸出路徑
        outputs = {
            'polygons': self.output_dir / "proseg_results.json",
            'counts': self.output_dir / "counts.csv.gz",  # Changed to .gz
            'cells': self.output_dir / "cells.csv",
            'genes': self.output_dir / "genes.csv",
        }

        # Smart Resume Check
        all_exist = all(path.exists() for path in outputs.values())
        if all_exist:
             logger.info("✅ 偵測到 Proseg 輸出檔案已存在，跳過執行 (Smart Resume)")
             for name, path in outputs.items():
                 logger.info(f"  - {name}: {path}")
             self.proseg_results = outputs
             return outputs

        # 建立 Proseg 指令
        cmd = [
            self.proseg_bin,
            "--overwrite",
            "--output-cell-polygons", str(outputs['polygons']),
            "--output-counts", str(outputs['counts']), # Already has .gz
            "--output-counts-fmt", "csv-gz",
            "--output-cell-metadata", str(outputs['cells']),
            "--output-cell-metadata-fmt", "csv",
            "--output-gene-metadata", str(outputs['genes']),
            "--output-gene-metadata-fmt", "csv",
            "--coordinate-scale", str(self.coordinate_scale),
            "--gene-column", "gene",
            "--x-column", "x",
            "--y-column", "y",
            "--z-column", "z",
            "--cell-id-column", "cell_id",
            "--cell-id-unassigned", "0",
            "--ignore-z-coord",
            "--min-qv", "0",
            "--max-transcript-nucleus-distance", str(self.max_dist),
            "--cell-compactness", str(self.compactness),
            "--samples", str(self.samples),
            "--burnin-samples", str(self.burnin_samples),
            "--recorded-samples", str(self.recorded_samples),
            "--nuclear-reassignment-prob", "0.01",  # 新增：嚴格保護核區域，防止越界擴張
        ]

        if self.enforce_connectivity:
            cmd.append("--enforce-connectivity")

        cmd.append(str(csv_path))

        logger.info(f"執行指令:")
        logger.info(" ".join(cmd))
        logger.info("⏳ Proseg 執行中...")

        try:
            result = subprocess.run(
                cmd,
                check=True,
                capture_output=True,
                text=True
            )
            logger.debug(result.stdout)
            logger.info("✅ Proseg 執行完成")

            # 驗證輸出
            for name, path in outputs.items():
                if not path.exists():
                    raise FileNotFoundError(f"Proseg 輸出缺少: {name} ({path})")
                logger.debug(f"  - {name}: {path}")

            self.proseg_results = outputs
            return outputs

        except subprocess.CalledProcessError as e:
            logger.error(f"❌ Proseg 執行失敗！")
            logger.error(f"錯誤訊息: {e.stderr}")
            raise

    def assemble_anndata(self) -> AnnData:
        """
        組裝 Proseg 結果為 AnnData

        Returns
        -------
        AnnData
            組裝的 AnnData 物件
        """
        logger.info("=" * 60)
        logger.info("步驟 4: 組裝 AnnData")
        logger.info("=" * 60)

        # 讀取 Proseg 輸出
        logger.info("讀取 Proseg 輸出檔案...")

        def read_csv_safe(path):
            """Helper to read CSV that might be gzipped"""
            try:
                # Try reading as plain CSV first
                return pd.read_csv(path)
            except UnicodeDecodeError:
                logger.warning(f"  ⚠️  偵測到 gzip 壓縮: {path.name}")
                # Retry with gzip compression
                return pd.read_csv(path, compression='gzip')

        counts_df = None # 不再使用 DataFrame 讀取 counts
        cells_df = read_csv_safe(self.proseg_results['cells'])
        genes_df = read_csv_safe(self.proseg_results['genes'])

        logger.debug(f"  - Cells columns: {list(cells_df.columns)}")
        logger.debug(f"  - Genes columns: {list(genes_df.columns)}")

        # 讀取計數矩陣 (Matrix Market format)
        # Note: Proseg seems to force Matrix Market format even if --output-counts-fmt csv is requested,
        # or at least the file header indicates MatrixMarket. We handle this robustly.
        logger.info("讀取計數矩陣 (Matrix Market / Gzip)...")
        from scipy.io import mmread
        from scipy.sparse import csr_matrix
        import gzip

        def read_mtx_safe(path):
            """Helper to read MTX that might be gzipped"""
            path = Path(path)
            # Check for gzip magic number
            try:
                with open(path, 'rb') as f:
                    header = f.read(2)
            except Exception:
                # Could vary if path is string
                pass

            is_gzipped = (header == b'\x1f\x8b')

            if is_gzipped:
                logger.debug(f"  - 偵測到 gzip 壓縮: {path.name}")
                with gzip.open(path, 'rt') as f:
                    return mmread(f)
            else:
                return mmread(str(path))

        try:
             X = read_mtx_safe(self.proseg_results['counts'])
        except Exception as e:
             logger.error(f"  ❌ MTX 讀取失敗: {e}")
             raise ValueError(f"無法讀取計數矩陣: {e}")


        # 轉換為 CSR 格式以便快速操作
        X = X.tocsr()

        logger.debug(f"  - Matrix shape: {X.shape}")

        # 驗證維度
        # MTX header: rows(cells) cols(genes)
        # cells_df rows should match X.shape[0]
        # genes_df rows should match X.shape[1]

        if X.shape[0] != len(cells_df):
             # 嘗試轉置 (有時候是 gene x cell)
             if X.shape[1] == len(cells_df) and X.shape[0] == len(genes_df):
                 logger.warning("  ⚠️  Matrix dimensions transposed (Genes x Cells), transposing to (Cells x Genes)...")
                 X = X.T
             else:
                 raise ValueError(
                     f"Matrix dimensions {X.shape} do not match metadata: "
                     f"Cells ({len(cells_df)}), Genes ({len(genes_df)})"
                 )

        # 再次確認
        if X.shape[0] != len(cells_df):
            raise ValueError(f"Cell count mismatch: Matrix {X.shape[0]} != Metadata {len(cells_df)}")
        if X.shape[1] != len(genes_df):
            raise ValueError(f"Gene count mismatch: Matrix {X.shape[1]} != Metadata {len(genes_df)}")


        # 建立 AnnData
        logger.info("建立 AnnData...")

        # 設定 metadata 索引
        # Proseg 通常保證 metadata 順序與矩陣索引一致 (1-based index in MTX maps to 0-based index in file)

        # 處理 Cells metadata
        if 'cell' in cells_df.columns:
            cells_df = cells_df.set_index('cell')
        elif 'cell_id' in cells_df.columns:
            cells_df = cells_df.set_index('cell_id')

        # 處理 Genes metadata
        if 'gene' in genes_df.columns:
            genes_df = genes_df.set_index('gene')
        elif 'gene_id' in genes_df.columns:
            genes_df = genes_df.set_index('gene_id')

        adata = AnnData(X=X)
        adata.obs = cells_df
        adata.var = genes_df

        # 加入空間座標
        # columns might be 'centroid_x', 'centroid_y' based on previous output or 'cx', 'cy'
        # Previous debug output: ['cell', 'original_cell_id', 'centroid_x', 'centroid_y', 'centroid_z', ...]

        spatial_cols = []
        if 'centroid_x' in cells_df.columns and 'centroid_y' in cells_df.columns:
             spatial_cols = ['centroid_x', 'centroid_y']
        elif 'cx' in cells_df.columns and 'cy' in cells_df.columns:
             spatial_cols = ['cx', 'cy']

        if spatial_cols:
            logger.info(f"  - 加入空間座標 (自動校正回全域): {spatial_cols}")
            # 取得局部座標
            local_coords = cells_df[spatial_cols].values

            # 校正回全域座標
            # global = ( (local_um / coordinate_scale) + roi_offset_px ) * scale_factors_px_to_nm
            roi_min_x, roi_min_y = self.roi_offset
            scale_y, scale_x = self.scale_factors

            global_coords = local_coords.copy()
            # 1. 將 Proseg 輸出的微米單位轉回局部像素
            local_px_x = global_coords[:, 0] / self.coordinate_scale
            local_px_y = global_coords[:, 1] / self.coordinate_scale

            # 2. 加上 ROI 偏移量得到全域像素，再乘上 Zarr 的 Scale 得到全域奈米
            global_coords[:, 0] = (local_px_x + roi_min_x) * scale_x # x (Global NM)
            global_coords[:, 1] = (local_px_y + roi_min_y) * scale_y # y (Global NM)

            adata.obsm['spatial'] = global_coords

        logger.info(f"✅ AnnData 組裝完成")
        logger.info(f"  - 形狀: {adata.shape}")
        logger.info(f"  - 細胞數: {adata.n_obs}")
        logger.info(f"  - 基因數: {adata.n_vars}")

        return adata

    def save_outputs(self, adata: AnnData) -> None:
        """
        儲存輸出（H5AD 和 Zarr）

        Parameters
        ----------
        adata : AnnData
            要儲存的 AnnData 物件
        """
        logger.info("=" * 60)
        logger.info("步驟 5: 儲存輸出")
        logger.info("=" * 60)

        # 儲存 H5AD
        h5ad_path = self.output_dir / "proseg_integrated.h5ad"
        logger.info(f"儲存 H5AD: {h5ad_path}")
        adata.write_h5ad(h5ad_path)

        # 回寫至 Zarr
        zarr_out_path = self.output_dir / "proseg-output.zarr"
        logger.info(f"儲存 Zarr: {zarr_out_path}")

        # Proseg Polygons (Shapes) Integration
        # ------------------------------------
        polygons_json_path = self.proseg_results.get('polygons')
        if polygons_json_path and polygons_json_path.exists():
            logger.info(f"處理 Proseg 多邊形遮罩: {polygons_json_path}")
            try:
                import geopandas as gpd
                from shapely.geometry import Polygon, MultiPolygon
                import gzip

                def read_json_safe(path):
                    with open(path, 'rb') as f:
                        header = f.read(2)
                    is_gzipped = (header == b'\x1f\x8b')

                    if is_gzipped:
                        logger.info(f"  🔍 偵測到 Gzip 壓縮: {path.name}")
                        with gzip.open(path, 'rt', encoding='utf-8') as f:
                            return json.load(f)
                    else:
                        with open(path, 'r', encoding='utf-8') as f:
                            return json.load(f)

                poly_data = read_json_safe(polygons_json_path)

                # Proseg JSON Format check (it might be list of lists or dict)
                # Usually: {"cell_id": [[x1,y1], [x2,y2], ...]} or list of dicts
                # Assuming Proseg outputs GeoJSON or similar structure

                # Check structure
                geometries = []
                ids = []

                # Global Transform Parameters
                roi_min_x, roi_min_y = self.roi_offset
                scale_y, scale_x = self.scale_factors

                # Proseg Output is GeoJSON FeatureCollection
                if isinstance(poly_data, dict) and 'features' in poly_data:
                    logger.info("  🔍 偵測到 GeoJSON FeatureCollection")
                    features = poly_data['features']
                else:
                    # Fallback or older format
                    logger.warning(f"  ⚠️  未知 JSON 結構，嘗試簡易解析... (Keys: {list(poly_data.keys())[:3]})")
                    features = [] # TODO: Handle map?

                count = 0
                for feat in features:
                    props = feat.get('properties', {})
                    # Proseg GeoJSON uses 'cell' as the ID key based on inspection
                    cell_id = props.get('cell', props.get('id', props.get('cell_id')))

                    if cell_id is None: continue

                    # Geometry
                    geom_data = feat.get('geometry', {})
                    geom_type = geom_data.get('type')
                    coords = geom_data.get('coordinates')
                    if not coords: continue

                    # Robust coordinate parsing
                    try:
                        ring = None

                        if geom_type == 'Polygon':
                            # Polygon: [ [x,y]... ] (rings)
                            # coords[0] is exterior ring
                            if len(coords) > 0:
                                ring = coords[0]

                        elif geom_type == 'MultiPolygon':
                            # MultiPolygon: [ [ [x,y]... ] ] (polygons -> rings)
                            # coords[0] is first polygon
                            # coords[0][0] is exterior ring of first polygon
                            if len(coords) > 0 and len(coords[0]) > 0:
                                ring = coords[0][0]

                        if ring is None or len(ring) < 3:
                            continue

                        local_pts = np.array(ring)

                        # Validate shape (N, 2)
                        if local_pts.ndim != 2 or local_pts.shape[1] != 2:
                             logger.warning(f"  ⚠️  細胞 {cell_id} 多邊形維度錯誤: {local_pts.shape}")
                             continue

                        # Transform to Global Physical
                        global_pts = local_pts.copy()
                        global_pts[:, 0] = (global_pts[:, 0] + roi_min_x) * scale_x # x
                        global_pts[:, 1] = (global_pts[:, 1] + roi_min_y) * scale_y # y

                        poly = Polygon(global_pts)
                        geometries.append(poly)
                        ids.append(int(cell_id))
                        count += 1

                    except Exception as e:
                         # logger.debug(f"Parsing geometry for cell {cell_id} failed: {e}")
                         continue

                if geometries:
                    gdf = gpd.GeoDataFrame(
                        {'cell_id': ids},
                        geometry=geometries
                    )

                    # Write to Zarr
                    add_shapes_to_zarr(
                        zarr_path=str(self.zarr_path), # Use input Zarr path (modify in-place)
                        shapes_name="proseg_polygons",
                        gdf=gdf,
                        overwrite=True
                    )
                    logger.info(f"✅ Proseg 多邊形 ({count} 個) 已寫入 Zarr: shapes/proseg_polygons")
                else:
                    logger.warning("  ⚠️  無有效多邊形可轉換")

            except Exception as e:
                logger.error(f"  ❌ 多邊形處理失敗: {e}")
                # Log traceback for debugging
                import traceback
                logger.debug(traceback.format_exc())

        # 檢查是否支援寫入 (SimplifiedSpatialData 不支援)
        if not hasattr(self.sdata, 'write') or not hasattr(self.sdata, 'copy'):
            logger.warning(f"  ⚠️  警告: 目前的 SpatialData 物件 (容錯模式) 不支援寫回 Zarr。")
            logger.warning(f"  ⚠️  已跳過 Zarr 輸出。請使用 {h5ad_path} 進行後續分析。")
            logger.info(f"✅ 輸出儲存完成 (僅 H5AD)")
            logger.info(f"  - H5AD: {h5ad_path}")
            return

        # 複製原始 Zarr
        # if zarr_out_path.exists():
        #    shutil.rmtree(zarr_out_path)
        #
        # try:
        #     # 讀取原始 SpatialData 並加入新 table
        #     sdata_out = self.sdata.copy()
        #     sdata_out.tables['proseg_cells'] = adata
        #     sdata_out.write(zarr_out_path)

        #     logger.info(f"✅ 輸出儲存完成")
        #     logger.info(f"  - H5AD: {h5ad_path}")
        #     logger.info(f"  - Zarr: {zarr_out_path}")
        # except Exception as e:
        #     logger.error(f"  ❌ Zarr 儲存失敗: {e}")
        #     logger.error(f"  (這不影響 H5AD 的使用)")


    def run_full_pipeline(self) -> AnnData:
        """
        執行完整 Pipeline

        Returns
        -------
        AnnData
            最終的 AnnData 物件
        """
        logger.info("=" * 80)
        logger.info(" " * 20 + "🚀 Proseg Pipeline 開始執行")
        logger.info("=" * 80)

        try:
            # 1. 載入資料
            self.load_data()

            # 2. 準備 Proseg 輸入
            csv_path = self.prepare_proseg_input()

            # 3. 執行 Proseg
            self.run_proseg(csv_path)

            # 4. 組裝 AnnData
            adata = self.assemble_anndata()

            # 5. 儲存輸出
            self.save_outputs(adata)

            logger.info("=" * 80)
            logger.info(" " * 20 + "✅ Pipeline 執行成功！")
            logger.info("=" * 80)

            return adata

        except Exception as e:
            logger.error("=" * 80)
            logger.error(" " * 20 + "❌ Pipeline 執行失敗")
            logger.error("=" * 80)
            logger.error(f"錯誤: {e}")
            raise

    def run_full(self) -> None:
        """執行完整 Proseg 流程（API 呼叫入口）"""
        self.run_full_pipeline()
