"""
Stage 5 — Xenium Explorer 匯出模組

將 Proseg 分割結果（h5ad + GeoJSON 多邊形 + 轉錄點 CSV）組裝成
spatialdata_xenium_explorer 可讀取的 Xenium Explorer bundle。

座標系設計
----------
遵循 SpatialData 慣例：shapes/points 以物理 µm 傳入，影像以 pixel 傳入並
附加 Scale 轉換標記其解析度。write() 與 Scale 使用相同的 pixel_size_um，
避免 library 重採樣不一致（與 Proseg-Zarr-Integration 做法相同）。

1. 影像：原始 pixel + Scale([pixel_size_um, pixel_size_um]) 轉換
2. 多邊形：GeoJSON µm → ShapesModel（Identity，物理 µm）
3. 轉錄點：CSV/zarr µm → PointsModel（Identity，物理 µm）
4. experiment.xenium：pixel_size = pixel_size_um（修補 spatialdata_xenium_explorer bug）

移植自：Proseg-Zarr-Integration/scripts/export_to_xenium_full.py
"""

from __future__ import annotations

import inspect
import json
import logging
from pathlib import Path
from typing import Optional

import anndata
import geopandas as gpd
import numpy as np
import pandas as pd
import shapely.geometry
import spatialdata as sd
from spatialdata.models import (
    Image2DModel,
    PointsModel,
    ShapesModel,
    TableModel,
)
from spatialdata.transformations import Identity, Scale

from backend.src.utils.constants import VISIUM_UM_PX

logger = logging.getLogger("pipeline.export.xenium")


def _apply_xenium_explorer_patches() -> None:
    """
    Monkey-patch spatialdata_xenium_explorer.utils.to_intrinsic 使其能
    處理 element=None 的情況（無轉錄點時 get_element 回傳 None，
    原始實作直接呼叫 transform_element_to_coordinate_system(None, cs)
    觸發 AssertionError）。
    """
    try:
        import spatialdata_xenium_explorer.utils as _sx_utils

        if getattr(_sx_utils, "_to_intrinsic_patched", False):
            return  # 已 patch，不重複

        _orig = _sx_utils.to_intrinsic

        def _safe_to_intrinsic(sdata, element, element_cs):
            if element is None:
                return None
            return _orig(sdata, element, element_cs)

        _sx_utils.to_intrinsic = _safe_to_intrinsic
        _sx_utils._to_intrinsic_patched = True
        logger.debug("spatialdata_xenium_explorer.utils.to_intrinsic patched (None-safe)")
    except Exception as exc:
        logger.warning(f"xenium_explorer patch 失敗（不影響主流程）：{exc}")


class XeniumExporter:
    """
    將 Proseg 分割結果匯出為 Xenium Explorer bundle。

    Parameters
    ----------
    zarr_path:
        原始 Zarr 資料集路徑（可選，提供時會包含影像層）。
    poly_json_path:
        combined_proseg_results_qc.json 路徑，包含多邊形（µm 座標）。
    transcripts_csv_path:
        combined_transcripts.csv 路徑（x/y 欄位為 µm）。
    ram_threshold_gb:
        spatialdata_xenium_explorer.write() 的 RAM 閾值。
    he_image_path:
        H&E 影像路徑（BTF/TIFF，可選）。合併模式下無 zarr 時使用。
    he_crop_bounds:
        (x0, y0, x1, y1) 影像 pixel 座標，指定從 BTF 裁切的區域。
        None = 讀取整張（小圖使用）。
    """

    def __init__(
        self,
        zarr_path: Optional[str | Path] = None,
        poly_json_path: Optional[str | Path] = None,
        transcripts_csv_path: Optional[str | Path] = None,
        ram_threshold_gb: float = 4.0,
        pixel_size_um: float = VISIUM_UM_PX,
        he_image_path: Optional[str | Path] = None,
        he_crop_bounds: Optional[tuple] = None,
    ) -> None:
        self.zarr_path = Path(zarr_path) if zarr_path else None
        self.poly_json_path = Path(poly_json_path) if poly_json_path else None
        self.transcripts_csv_path = (
            Path(transcripts_csv_path) if transcripts_csv_path else None
        )
        self.ram_threshold_gb = ram_threshold_gb
        self.pixel_size_um = pixel_size_um  # µm/px of the morphology image
        self.he_image_path = Path(he_image_path) if he_image_path else None
        self.he_crop_bounds = he_crop_bounds  # (x0, y0, x1, y1) in image pixels

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def export(self, h5ad_path: str | Path, output_dir: str | Path) -> Path:
        """
        主要匯出入口。

        Parameters
        ----------
        h5ad_path:
            已過 QC 的 AnnData h5ad 檔案路徑。
        output_dir:
            Xenium Explorer bundle 輸出目錄（會自動建立）。

        Returns
        -------
        Path
            實際寫出的 Xenium Explorer 目錄路徑。
        """
        h5ad_path = Path(h5ad_path)
        out_xenium_dir = Path(output_dir)
        out_xenium_dir.mkdir(parents=True, exist_ok=True)

        logger.info("=== Stage 5 Xenium Explorer 匯出開始 ===")
        logger.info(f"h5ad: {h5ad_path}")
        logger.info(f"輸出目錄: {out_xenium_dir}")

        # 1. 載入影像（可選）
        sd_image = self._load_image()

        # 2. 載入 AnnData
        adata = self._load_adata(h5ad_path)
        if adata is None:
            raise RuntimeError(f"無法載入 AnnData: {h5ad_path}")

        # 3. 載入多邊形並建立 ID 重映射
        sd_shapes, sd_table, id_remap = self._load_polygons_and_table(adata)

        # 4. 載入轉錄點（優先 CSV，其次 zarr parquet）
        sd_points = self._load_transcripts(id_remap)
        if sd_points is None:
            sd_points = self._load_transcripts_from_zarr()

        # 5. 組裝 SpatialData 並寫出
        self._write_xenium_bundle(
            out_xenium_dir=out_xenium_dir,
            sd_image=sd_image,
            sd_shapes=sd_shapes,
            sd_table=sd_table,
            sd_points=sd_points,
        )

        logger.info(f"=== Xenium Explorer bundle 完成：{out_xenium_dir} ===")
        return out_xenium_dir

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    def _load_image(self) -> Optional[sd.models.Image2DModel]:
        """
        從 Zarr 惰性載入影像，附加 Scale 轉換（pixel_size_um µm/px）。

        spatialdata_xenium_explorer.write() 遵循 SpatialData 慣例：
        影像座標系為 pixel，透過 Scale 轉換告知 library 每 pixel 的物理尺寸（µm）。
        Library 在 write() 時自行根據 pixel_size 將影像重採樣到 Xenium 輸出解析度。
        不在此手動縮放影像，避免座標空間錯亂。
        """
        if self.zarr_path is None or not self.zarr_path.exists():
            logger.info("未提供 zarr_path，跳過影像層。")
            # 嘗試供給 he_image_path（BTF/TIFF）
            if self.he_image_path and self.he_image_path.exists():
                return self._load_he_image()
            return None

        try:
            import dask.array as da
            import numpy as np
            import zarr
            from spatialdata.transformations import Scale

            logger.info(f"載入 Zarr 影像：{self.zarr_path}")
            z = zarr.open(str(self.zarr_path), mode="r")
            img_array = z["images"]["tissue_hires_image"]["0"]

            # 統一為 (c, y, x) 格式
            raw = np.array(img_array)
            if raw.shape[0] not in (1, 3):
                raw = raw.transpose((2, 0, 1))

            dask_img = da.from_array(raw, chunks=(1, 2048, 2048))

            # Scale 轉換：告知 spatialdata 每 pixel = self.pixel_size_um µm
            # library 根據此資訊與 write() 的 pixel_size 自行處理重採樣
            sd_image = Image2DModel.parse(
                dask_img, dims=("c", "y", "x"),
                transformations={"global": Scale(
                    [self.pixel_size_um, self.pixel_size_um], axes=("y", "x")
                )},
            )
            logger.info(f"影像載入完成，shape={dask_img.shape}，pixel_size={self.pixel_size_um} µm/px")
            return sd_image

        except Exception as exc:
            logger.warning(f"Zarr 影像載入失敗（繼續執行）：{exc}")
            return None

    def _load_he_image(self) -> Optional["sd.models.Image2DModel"]:
        """
        從 BTF/TIFF 以 tiled 讀取方式惰性載入 H&E 影像。

        避開 zarr v2/v3 相容性問題：對無壓縮 tiled TIFF，直接以 file seek +
        dask.delayed 讀取 crop 區域的各 tile，不依賴 zarr store。

        Transform：Sequence([Scale, Translation])
        - Scale:       image pixel → µm
        - Translation: crop 左上角偏移（全域 µm）
        """
        try:
            import tifffile
            import dask
            import dask.array as da
            from spatialdata.transformations import Scale

            if self.he_image_path is None or not self.he_image_path.exists():
                return None

            logger.info(f"載入 H&E 影像（BTF）：{self.he_image_path}")

            with tifffile.TiffFile(str(self.he_image_path)) as tif:
                page = tif.series[0].pages[0]
                img_h = page.imagelength
                img_w = page.imagewidth
                tile_h = int(page.tilelength) if page.is_tiled else img_h
                tile_w = int(page.tilewidth)  if page.is_tiled else img_w
                n_channels = page.shape[2] if len(page.shape) > 2 else 1
                dtype = page.dtype
                offsets = list(page.dataoffsets)
                bytecounts = list(page.databytecounts)

            # 裁切範圍（image pixels）
            if self.he_crop_bounds is not None:
                x0, y0, x1, y1 = (int(v) for v in self.he_crop_bounds)
            else:
                x0, y0, x1, y1 = 0, 0, img_w, img_h

            x0 = max(0, x0);  y0 = max(0, y0)
            x1 = min(img_w, x1);  y1 = min(img_h, y1)
            crop_w = x1 - x0;  crop_h = y1 - y0

            logger.info(f"  擷取區域：x=[{x0},{x1}], y=[{y0},{y1}]，大小 {crop_w}×{crop_h} px")

            path_str = str(self.he_image_path)

            @dask.delayed
            def _read_tiled_crop(
                path, y0_, y1_, x0_, x1_,
                img_w_, img_h_, tile_h_, tile_w_, nc_, offsets_, bytecounts_
            ):
                """
                讀取 uncompressed tiled TIFF 的 crop 區域，只讀必要 tiles。
                對壓縮 tile 自動降回 tifffile.imread 整張讀取後裁切。
                """
                import numpy as np
                import tifffile as _tf

                crop_h_ = y1_ - y0_
                crop_w_ = x1_ - x0_
                result = np.zeros((crop_h_, crop_w_, nc_), dtype=np.uint8)

                ntiles_x_ = (img_w_ + tile_w_ - 1) // tile_w_
                ntiles_y_ = (img_h_ + tile_h_ - 1) // tile_h_

                ty0_ = y0_ // tile_h_
                ty1_ = min((y1_ - 1) // tile_h_ + 1, ntiles_y_)
                tx0_ = x0_ // tile_w_
                tx1_ = min((x1_ - 1) // tile_w_ + 1, ntiles_x_)

                try:
                    with open(path, "rb") as fh:
                        for ty in range(ty0_, ty1_):
                            for tx in range(tx0_, tx1_):
                                tidx = ty * ntiles_x_ + tx
                                if tidx >= len(offsets_) or offsets_[tidx] == 0:
                                    continue
                                fh.seek(offsets_[tidx])
                                raw_bytes = fh.read(bytecounts_[tidx])
                                tile = np.frombuffer(raw_bytes, dtype=np.uint8).reshape(
                                    tile_h_, tile_w_, nc_
                                )
                                gy0 = ty * tile_h_;  gy1 = min(gy0 + tile_h_, img_h_)
                                gx0 = tx * tile_w_;  gx1 = min(gx0 + tile_w_, img_w_)
                                oy0 = max(y0_, gy0);  oy1 = min(y1_, gy1)
                                ox0 = max(x0_, gx0);  ox1 = min(x1_, gx1)
                                if oy1 <= oy0 or ox1 <= ox0:
                                    continue
                                result[oy0-y0_:oy1-y0_, ox0-x0_:ox1-x0_] = \
                                    tile[oy0-gy0:oy1-gy0, ox0-gx0:ox1-gx0]
                except Exception:
                    # 壓縮 tile 無法直接讀取 raw bytes，退回 tifffile.imread 整張裁切
                    full = _tf.imread(path, series=0, level=0)
                    result = full[y0_:y1_, x0_:x1_]

                return result

            raw = da.from_delayed(
                _read_tiled_crop(
                    path_str, y0, y1, x0, x1,
                    img_w, img_h, tile_h, tile_w, n_channels,
                    offsets, bytecounts,
                ),
                shape=(crop_h, crop_w, n_channels),
                dtype=dtype,
            )

            # (y, x, c) → (c, y, x)，只取 RGB
            raw = da.transpose(raw, (2, 0, 1))[:3]

            # Transform：Sequence([Scale, Translation])
            # Scale: image pixel → µm
            # Translation: 加上 crop 左上角的全域偏移（µm），確保影像原點對齊多邊形全域座標
            ps = self.pixel_size_um
            from spatialdata.transformations import Sequence, Translation
            if x0 != 0 or y0 != 0:
                transform = Sequence([
                    Scale([ps, ps], axes=("y", "x")),
                    Translation([y0 * ps, x0 * ps], axes=("y", "x")),
                ])
            else:
                transform = Scale([ps, ps], axes=("y", "x"))

            sd_image = Image2DModel.parse(
                raw,
                dims=("c", "y", "x"),
                transformations={"global": transform},
            )
            logger.info(
                f"  H&E 影像載入完成，crop={crop_w}×{crop_h} px，"
                f"pixel_size={ps:.4f} µm/px，"
                f"全域偏移=({x0 * ps:.1f}, {y0 * ps:.1f}) µm"
            )
            return sd_image

        except Exception as exc:
            logger.warning(f"H&E 影像載入失敗（繼續執行）：{exc}")
            return None

    def _load_adata(self, h5ad_path: Path) -> Optional[anndata.AnnData]:
        """載入並初始化 AnnData。"""
        try:
            logger.info(f"載入 AnnData：{h5ad_path}")
            adata = anndata.read_h5ad(h5ad_path)
            adata.obs["region"] = "cell_boundaries"
            logger.info(f"AnnData 載入完成，{len(adata)} 個細胞，{adata.n_vars} 個基因。")
            return adata
        except Exception as exc:
            logger.error(f"AnnData 載入失敗：{exc}")
            return None

    def _load_polygons_and_table(
        self, adata: anndata.AnnData
    ) -> tuple[
        Optional[gpd.GeoDataFrame],
        Optional[sd.models.TableModel],
        dict,
    ]:
        """
        讀取 GeoJSON 多邊形，進行座標轉換（µm → px）與平滑化，
        並同步修剪 AnnData（只保留有有效多邊形的細胞），
        最後建立密集 instance_id 重映射。

        Returns
        -------
        sd_shapes, sd_table, id_remap
        """
        if self.poly_json_path is None or not self.poly_json_path.exists():
            logger.warning(f"多邊形 JSON 不存在：{self.poly_json_path}，跳過 shapes/table。")
            return None, None, {}

        logger.info(f"載入多邊形 JSON：{self.poly_json_path}")

        try:
            with open(self.poly_json_path, "r", encoding="utf-8") as f:
                geo_data = json.load(f)
        except Exception as exc:
            logger.error(f"讀取多邊形 JSON 失敗：{exc}")
            return None, None, {}

        # 建立 ID 查找表，支援多種 ID 格式（原始 ID, cell_id 數值, 'cell_N' 格式）
        id_to_obs = {str(name): name for name in adata.obs_names}
        if "cell_id" in adata.obs.columns:
            for obs_name, cid in zip(adata.obs_names, adata.obs["cell_id"]):
                id_to_obs[str(int(cid))] = obs_name

        # 處理 'cell_N' 格式的對應（單 ROI）
        for name in adata.obs_names:
            if name.startswith("cell_"):
                try:
                    num_id = name.split("_")[1]
                    id_to_obs[num_id] = name
                except Exception:
                    pass

        # Merge 模式：obs_names 為 "{roi}__cell_{N}"，GeoJSON full_id 為 "{roi}__{N}"
        # 補上 "{roi}__{N}" → "{roi}__cell_{N}" 的映射
        for name in adata.obs_names:
            if "__cell_" in name:
                parts = name.split("__cell_", 1)
                if len(parts) == 2:
                    id_to_obs[f"{parts[0]}__{parts[1]}"] = name

        features = geo_data.get("features", [])
        polygons: list = []
        valid_cell_ids: list[str] = []

        obs_set = set(adata.obs_names)

        for feat in features:
            props = feat.get("properties", {})
            
            # 優先從 full_id 找，找不到則從 cell/cell_id 找
            raw_id = props.get("full_id") or props.get("cell_id") or props.get("cell")
            if raw_id is None:
                continue
            
            str_id = str(raw_id)
            if str_id not in id_to_obs:
                # 針對 LUAD 這種 'cell_N' 與 'N' 的對應
                alt_id = f"cell_{str_id}"
                if alt_id in obs_set:
                    full_id = alt_id
                else:
                    # 嘗試將 '1.0' 轉為 '1'
                    try:
                        f_id = str(int(float(str_id)))
                        if f_id in id_to_obs: full_id = id_to_obs[f_id]
                        elif f"cell_{f_id}" in obs_set: full_id = f"cell_{f_id}"
                        else: continue
                    except Exception: continue
            else:
                full_id = id_to_obs[str_id]

            geom_type = feat["geometry"]["type"]
            coords = feat["geometry"]["coordinates"]

            try:
                if geom_type == "Polygon":
                    poly_um = shapely.geometry.Polygon(coords[0])
                elif geom_type == "MultiPolygon":
                    poly_um = shapely.geometry.MultiPolygon(
                        [shapely.geometry.Polygon(r[0]) for r in coords]
                    )
                else:
                    continue

                if not poly_um.is_valid or poly_um.is_empty:
                    continue

                # 保持 µm 座標傳入 SpatialData（SpatialData 慣例：全域座標為物理 µm）
                # spatialdata_xenium_explorer.write() 根據 pixel_size 自行轉換為 Xenium px
                poly_out = self._smooth_polygon(poly_um)

                # Xenium Explorer 要求全部為 MultiPolygon
                if poly_out.geom_type == "Polygon":
                    poly_out = shapely.geometry.MultiPolygon([poly_out])

                polygons.append(poly_out)
                valid_cell_ids.append(full_id)

            except Exception:
                continue

        logger.info(f"有效多邊形數量：{len(polygons)}（原始 AnnData：{len(adata)} 個細胞）")

        # 同步修剪 AnnData：只保留有有效多邊形的細胞
        # spatialdata_xenium_explorer 要求 instance_id 嚴格連續 0..N-1
        if len(valid_cell_ids) < len(adata):
            dropped = len(adata) - len(valid_cell_ids)
            logger.warning(f"刪除 {dropped} 個無有效多邊形的細胞。")
            adata = adata[valid_cell_ids].copy()

        # 建立密集 instance_id
        adata.obs["instance_id"] = pd.array(range(len(adata)), dtype="int64")
        id_remap: dict[str, int] = {
            cid: idx for idx, cid in enumerate(adata.obs_names)
        }
        logger.info(f"ID 重映射建立完成：{len(id_remap)} 個細胞。")

        # 建立 SpatialData 物件
        sd_table = TableModel.parse(
            adata,
            region="cell_boundaries",
            region_key="region",
            instance_key="instance_id",
        )

        gdf = gpd.GeoDataFrame(
            {"geometry": polygons},
            index=pd.Index(range(len(polygons)), dtype=int),
        )
        sd_shapes = ShapesModel.parse(gdf, transformations={"global": Identity()})
        logger.info(f"ShapesModel 建立完成：{len(gdf)} 個多邊形。")

        return sd_shapes, sd_table, id_remap

    @staticmethod
    def _smooth_polygon(poly: shapely.geometry.base.BaseGeometry) -> shapely.geometry.base.BaseGeometry:
        """
        消除 MCseg v2 分割產生的鋸齒邊緣，同時嚴格控制頂點數量。

        步驟：
        1. simplify(0.4)       — 移除微小變異（去掉像素鋸齒）
        2. buffer(+1.0, res=4) — 對外膨脹，去除尖角；resolution=4 避免插入過多曲線點
        3. buffer(-1.0, res=4) — 收縮回原大小，整平邊緣
        4. simplify(0.5)       — 最終壓縮，確保頂點數可被 Xenium Explorer 快速載入

        注意：預設 resolution=16 會讓頂點數爆增（20 頂點 → 100+），
        導致 cells.zarr.zip 極大，Xenium Explorer 「Loading cells...」卡住。
        resolution=4 + 最終 simplify 可將頂點控制在 ~30 個以內。
        """
        try:
            smooth = poly.simplify(0.4, preserve_topology=True)
            # resolution=4：每 1/4 圓弧只插 4 個點（vs 預設 16 個），大幅減少頂點數
            smooth = (
                smooth
                .buffer(1.0, join_style=1, resolution=4)
                .buffer(-1.0, join_style=1, resolution=4)
            )
            # 最終再 simplify：壓縮殘留冗餘頂點，容差 0.5µm 對 Visium HD 影響可忽略
            smooth = smooth.simplify(0.5, preserve_topology=True)
            if (
                smooth.is_valid
                and not smooth.is_empty
                and smooth.geom_type in ("Polygon", "MultiPolygon")
            ):
                return smooth
        except Exception:
            pass
        return poly

    def _load_transcripts(self, id_remap: dict) -> Optional[sd.models.PointsModel]:
        """
        惰性載入轉錄點 CSV（dask），將 x/y µm → px，
        並將字串 cell_id 重映射為連續整數。
        """
        if self.transcripts_csv_path is None or not self.transcripts_csv_path.exists():
            logger.info("未提供轉錄點 CSV，跳過 points 層。")
            return None

        try:
            import dask.dataframe as dd

            logger.info(f"載入轉錄點 CSV：{self.transcripts_csv_path}")
            tx_dd = dd.read_csv(str(self.transcripts_csv_path))

            # 保持 µm 座標（SpatialData 慣例），不手動 ÷ XENIUM_UM_PX
            # spatialdata_xenium_explorer.write() 根據 pixel_size 自行轉換

            # 預覽欄位名稱
            sample = pd.read_csv(self.transcripts_csv_path, nrows=3)
            logger.info(f"轉錄點欄位：{sample.columns.tolist()}")

            parse_kwargs: dict = {"coordinates": {"x": "x", "y": "y"}}

            if "cell_id" in sample.columns and id_remap:
                logger.info("映射 cell_id 字串 → 連續整數（透過 Dask map）...")
                tx_dd["cell_id"] = (
                    tx_dd["cell_id"].map(id_remap).fillna(-1).astype("int32")
                )
                parse_kwargs["instance_key"] = "cell_id"

            if "gene" in sample.columns:
                parse_kwargs["feature_key"] = "gene"

            sd_points = PointsModel.parse(
                tx_dd,
                sort=False,   # sort=True 會對大型 CSV 做全排序，速度極慢
                transformations={"global": Identity()},
                **parse_kwargs
            )
            logger.info("轉錄點 PointsModel 建立完成。")
            return sd_points

        except Exception as exc:
            logger.error(f"轉錄點載入失敗（繼續執行）：{exc}")
            return None

    def _load_transcripts_from_zarr(self) -> Optional[object]:
        """
        從 Zarr 的 points/transcripts parquet 直接載入轉錄點（µm 座標）。
        轉換為 Xenium px 座標後建立 PointsModel。
        """
        if self.zarr_path is None or not self.zarr_path.exists():
            return None

        parquet_path = (
            self.zarr_path / "points" / "transcripts"
            / "points.parquet" / "part.0.parquet"
        )
        if not parquet_path.exists():
            logger.info("zarr 中無轉錄點 parquet，跳過。")
            return None

        try:
            import dask.dataframe as dd

            logger.info(f"從 zarr 載入轉錄點：{parquet_path}")
            tx_dd = dd.read_parquet(str(parquet_path))

            # 保持 µm 座標（SpatialData 慣例），不手動 ÷ XENIUM_UM_PX

            sd_points = PointsModel.parse(
                tx_dd,
                sort=False,   # sort=True 對大型 parquet 全排序，速度極慢
                coordinates={"x": "x", "y": "y"},
                feature_key="gene",
                transformations={"global": Identity()},
            )
            logger.info("zarr 轉錄點 PointsModel 建立完成。")
            return sd_points

        except Exception as exc:
            logger.warning(f"zarr 轉錄點載入失敗（繼續執行）：{exc}")
            return None

    def _create_dummy_image(
        self, sd_shapes: Optional[gpd.GeoDataFrame]
    ) -> "sd.models.Image2DModel":
        """
        當無 zarr_path 時，依多邊形邊界框建立最小灰階佔位影像。

        spatialdata_xenium_explorer.write() 硬性要求 sdata.images 中至少有一張影像，
        此方法用空白 numpy array 滿足該需求，不影響下游分析。

        影像以 pixel 座標傳入，並附加 Scale 轉換標記物理解析度（self.pixel_size_um µm/px）。
        """
        import dask.array as da
        from spatialdata.transformations import Scale

        # 從多邊形邊界框推算影像範圍（µm → pixel）
        ps = self.pixel_size_um
        margin_um = 50.0
        if sd_shapes is not None and len(sd_shapes) > 0:
            bounds = sd_shapes.total_bounds  # (minx, miny, maxx, maxy) in µm
            x0_um = bounds[0] - margin_um
            y0_um = bounds[1] - margin_um
            width_px = max(1, int((bounds[2] - bounds[0] + 2 * margin_um) / ps))
            height_px = max(1, int((bounds[3] - bounds[1] + 2 * margin_um) / ps))
        else:
            x0_um, y0_um = 0.0, 0.0
            width_px, height_px = 512, 512

        # 限制最大大小避免 OOM（超大影像區域用稀疏表示即可）
        MAX_PX = 32768
        width_px = min(width_px, MAX_PX)
        height_px = min(height_px, MAX_PX)

        logger.info(
            f"建立佔位影像（{height_px}×{width_px} px，pixel_size={ps} µm/px，"
            f"全域偏移=({x0_um:.1f}, {y0_um:.1f}) µm）"
        )
        dummy = da.zeros((1, height_px, width_px), dtype=np.uint8, chunks=(1, 2048, 2048))

        # 加入 Translation 讓佔位影像與多邊形全域座標對齊
        from spatialdata.transformations import Sequence, Translation
        if x0_um != 0.0 or y0_um != 0.0:
            transform = Sequence([
                Scale([ps, ps], axes=("y", "x")),
                Translation([y0_um, x0_um], axes=("y", "x")),
            ])
        else:
            transform = Scale([ps, ps], axes=("y", "x"))

        return Image2DModel.parse(
            dummy,
            dims=("c", "y", "x"),
            transformations={"global": transform},
        )

    def _patch_experiment_xenium(self, out_xenium_dir: Path, pixel_size: float) -> None:
        """
        修補 spatialdata_xenium_explorer 的已知 bug：
        write() 可能在 experiment.xenium 寫入錯誤的 pixel_size，
        導致 Xenium Explorer 座標顯示錯誤。

        強制將 pixel_size 覆寫為實際使用的值，確保：
        - 細胞邊界座標（µm / pixel_size）與影像像素一一對應
        - Xenium Explorer 的比例尺正確
        """
        import json as _json

        exp_file = out_xenium_dir / "experiment.xenium"
        if not exp_file.exists():
            logger.warning(f"experiment.xenium 不存在，跳過 pixel_size 修補：{exp_file}")
            return

        try:
            with open(exp_file, "r", encoding="utf-8") as f:
                exp_data = _json.load(f)

            old_ps = exp_data.get("pixel_size", "未知")
            exp_data["pixel_size"] = pixel_size

            with open(exp_file, "w", encoding="utf-8") as f:
                _json.dump(exp_data, f, indent=2)

            logger.info(
                f"experiment.xenium pixel_size 修補完成：{old_ps} → {pixel_size} µm/px"
            )

            # 補齊 Xenium Explorer 必要的影像檔案
            # spatialdata_xenium_explorer 只寫出 morphology.ome.tif；
            # Xenium Explorer 還需要 morphology_mip.ome.tif 與 morphology_focus.ome.tif，
            # 缺少時會卡在 "Loading cells..." 等待這兩個檔案。
            # 用 morphology.ome.tif 複製補齊，H&E 無 z-stack，三檔等同。
            self._ensure_morphology_files(out_xenium_dir)

        except Exception as exc:
            logger.error(f"修補 experiment.xenium 失敗：{exc}")

    @staticmethod
    def _ensure_morphology_files(out_xenium_dir: Path) -> None:
        """
        Xenium Explorer morphology 衍生檔處理。

        Xenium Explorer 在 experiment.xenium 中引用 morphology_mip_filepath 與
        morphology_focus_filepath，但找不到時會直接跳過（graceful skip）。

        ⚠️ 不複製 morphology.ome.tif 作為 MIP/Focus：
        morphology.ome.tif 是 3 通道 RGB H&E 影像；
        Xenium Explorer 期望 MIP/Focus 為單通道 DAPI 灰階影像。
        格式不符時 Xenium Explorer 會在載入階段卡住（Loading cells... 數分鐘）。

        正確做法：讓 MIP/Focus 路徑不存在，Xenium Explorer 會自動跳過。
        """
        src = out_xenium_dir / "morphology.ome.tif"
        if not src.exists():
            logger.warning("morphology.ome.tif 不存在。")
            return

        # 若有舊版錯誤複製的 MIP/Focus 檔案，主動刪除以免卡住載入
        for fname in ("morphology_mip.ome.tif", "morphology_focus.ome.tif"):
            dst = out_xenium_dir / fname
            if dst.exists():
                try:
                    dst.unlink()
                    logger.info(f"移除格式錯誤的衍生影像：{fname}")
                except Exception as exc:
                    logger.warning(f"移除 {fname} 失敗（不影響主流程）：{exc}")

    @staticmethod
    def _rebuild_cells_zarr_v4(
        out_xenium_dir: Path,
        sd_shapes: "gpd.GeoDataFrame",
        pixel_size_um: float,
    ) -> None:
        """
        以 Xenium Explorer 4.x 新格式完整重建 cells.zarr.zip。

        spatialdata_xenium_explorer v0.1.7 寫出的是舊格式（polygon_vertices），
        Xenium Explorer 4.1.1 需要新格式：
          - polygon_sets/{0,1}/vertices   shape [n_cells, 50] float32，µm 座標
          - polygon_sets/{0,1}/num_vertices   shape [n_cells] int32
          - polygon_sets/{0,1}/cell_index     shape [n_cells] int32（0-based）
          - polygon_sets/{0,1}/method         shape [n_cells] uint32
          - masks/{0,1}                   shape [H, W] uint32，tile-chunked
          - masks/homogeneous_transform   4×4 float32，scale=1/pixel_size_um
          - cell_summary                  [n_cells, 8] float64
          - cell_id                       [n_cells, 2] uint32

        多邊形座標以 µm 傳入（ROI-local，與 morphology.ome.tif 像素一一對應）。
        mask 從 polygon 光柵化產生；chunks 為 256×256，避免單一大 chunk 阻塞載入。
        """
        import io
        import os
        import zipfile
        import zarr
        import numpy as np
        import numcodecs
        from zarr.storage import MemoryStore

        N_VERTS_MAX = 25   # 每個多邊形最多頂點數（real Xenium 亦為 25）
        MASK_CHUNK  = 256  # tile-based chunk size（避免單一 9MB+ chunk 阻塞 Xenium Explorer）
        BLOSC_LZ4   = numcodecs.Blosc(cname="lz4",  clevel=5, shuffle=numcodecs.Blosc.SHUFFLE)
        BLOSC_ZSTD  = numcodecs.Blosc(cname="zstd", clevel=1, shuffle=numcodecs.Blosc.SHUFFLE)

        cells_zarr_path = out_xenium_dir / "cells.zarr.zip"
        if not cells_zarr_path.exists():
            logger.warning("找不到 cells.zarr.zip（spatialdata_xenium_explorer 未寫出），跳過格式轉換。")
            return

        # ── 取得 morphology 影像尺寸（決定 mask 大小）────────────────────────
        morph_path = out_xenium_dir / "morphology.ome.tif"
        H, W = 512, 512  # 預設
        if morph_path.exists():
            try:
                import tifffile
                tf = tifffile.TiffFile(str(morph_path))
                page0 = tf.pages[0]
                H, W = page0.shape[0], page0.shape[1]
                tf.close()
            except Exception as exc:
                logger.warning(f"讀取 morphology.ome.tif 尺寸失敗：{exc}，使用預設 {H}×{W}")

        n_cells = len(sd_shapes)
        logger.info(f"重建 cells.zarr.zip：{n_cells} 個細胞，mask {H}×{W} px，tile={MASK_CHUNK}×{MASK_CHUNK}")

        # ── 多邊形 → 頂點陣列 ────────────────────────────────────────────────
        def _poly_to_verts(poly, n_max: int):
            """將 shapely Polygon/MultiPolygon 化簡為 ≤n_max 頂點，回傳 flat float32 與計數。"""
            import shapely
            if poly.geom_type == "MultiPolygon":
                poly = max(poly.geoms, key=lambda p: p.area)
            coords = list(poly.exterior.coords)
            if len(coords) > 1 and coords[0] == coords[-1]:
                coords = coords[:-1]
            # 化簡直到頂點數 ≤ n_max
            tol = 0.1
            while len(coords) > n_max and tol < 100:
                s = poly.simplify(tol, preserve_topology=True)
                if s.geom_type == "MultiPolygon":
                    s = max(s.geoms, key=lambda p: p.area)
                new_coords = list(s.exterior.coords)
                if new_coords[0] == new_coords[-1]:
                    new_coords = new_coords[:-1]
                coords = new_coords
                tol *= 2.0
            coords = coords[:n_max]
            flat = np.zeros(n_max * 2, dtype=np.float32)
            for i, (x, y) in enumerate(coords):
                flat[i * 2]     = float(x)
                flat[i * 2 + 1] = float(y)
            return flat, len(coords)

        vertices   = np.zeros((n_cells, N_VERTS_MAX * 2), dtype=np.float32)
        num_verts  = np.zeros(n_cells, dtype=np.int32)
        cell_index = np.arange(n_cells, dtype=np.int32)
        cx = np.zeros(n_cells, dtype=np.float64)
        cy = np.zeros(n_cells, dtype=np.float64)
        cell_area  = np.zeros(n_cells, dtype=np.float64)

        for i, poly in enumerate(sd_shapes.geometry):
            if poly is None or (hasattr(poly, "is_empty") and poly.is_empty):
                continue
            try:
                flat, nv    = _poly_to_verts(poly, N_VERTS_MAX)
                vertices[i] = flat
                num_verts[i] = nv
                cx[i]        = poly.centroid.x
                cy[i]        = poly.centroid.y
                cell_area[i] = poly.area
            except Exception as exc:
                logger.debug(f"細胞 {i} 頂點提取失敗：{exc}")

        # ── 光柵化多邊形 → mask ────────────────────────────────────────────────
        from skimage.draw import polygon as sk_polygon

        mask = np.zeros((H, W), dtype=np.uint32)
        scale = 1.0 / pixel_size_um   # µm → pixel
        n_rasterized = 0
        for i in range(n_cells):
            if num_verts[i] == 0:
                continue
            xs = vertices[i, 0::2][:num_verts[i]] * scale   # column (x)
            ys = vertices[i, 1::2][:num_verts[i]] * scale   # row    (y)
            # 若 polygon 超出 mask 邊界則跳過
            if xs.max() < 0 or xs.min() >= W or ys.max() < 0 or ys.min() >= H:
                continue
            try:
                rr, cc = sk_polygon(ys, xs, shape=(H, W))
                mask[rr, cc] = i + 1   # 1-indexed，0 = background
                n_rasterized += 1
            except Exception:
                pass
        logger.info(f"光柵化完成：{n_rasterized}/{n_cells} 個細胞")

        # ── cell_summary（8 欄）────────────────────────────────────────────────
        cell_summary = np.zeros((n_cells, 8), dtype=np.float64)
        cell_summary[:, 0] = cx          # cell_centroid_x (µm)
        cell_summary[:, 1] = cy          # cell_centroid_y (µm)
        cell_summary[:, 2] = cell_area   # cell_area (µm²)
        cell_summary[:, 3] = cx          # nucleus_centroid_x
        cell_summary[:, 4] = cy          # nucleus_centroid_y
        # [:, 5] nucleus_area = 0
        # [:, 6] z_level = 0
        cell_summary[:, 7] = 1.0         # nucleus_count

        # ── homogeneous_transform：µm → mask pixel ─────────────────────────
        transform = np.eye(4, dtype=np.float32)
        transform[0, 0] = scale
        transform[1, 1] = scale

        # ── cell_id：[n_cells, 2]，第 2 欄為 z-level（1 = in-focus） ────────
        # ⚠️ 必須 0-indexed（從 0 開始），與 cell_feature_matrix.zarr.zip 的 library 輸出一致。
        # 若用 1-indexed，Xenium Explorer 點選 cell 時會無限旋轉（ID 對不上）。
        cell_id_arr = np.zeros((n_cells, 2), dtype=np.uint32)
        cell_id_arr[:, 0] = np.arange(0, n_cells, dtype=np.uint32)
        cell_id_arr[:, 1] = 1

        # ── 寫入 MemoryStore，再序列化成 ZIP ──────────────────────────────────
        mem  = MemoryStore()
        root = zarr.open(mem, mode="w")

        # Root .zattrs — Xenium Explorer 4.x 必要 schema
        # polygon_sets/0 = cell boundaries，polygon_sets/1 = cell boundaries（無獨立核分割）
        root.attrs["major_version"]             = 6
        root.attrs["minor_version"]             = 2
        root.attrs["name"]                      = "CellSegmentationDataset"
        root.attrs["number_cells"]              = n_cells
        root.attrs["polygon_set_descriptions"]  = [
            "H&E cell segmentation by MCseg v2",
            "H&E cell segmentation by MCseg v2",
        ]
        root.attrs["polygon_set_display_names"] = ["Cell boundaries", "Cell boundaries"]
        root.attrs["polygon_set_names"]         = ["cell", "cell"]
        root.attrs["segmentation_methods"]      = ["MCseg v2 H&E cell segmentation"]
        root.attrs["spatial_units"]             = "microns"

        # cell_id
        root.array(
            "cell_id", cell_id_arr,
            chunks=(n_cells, 1), dtype=np.uint32,
            compressor=BLOSC_LZ4, overwrite=True,
        )

        # cell_summary + attrs
        cs_arr = root.array(
            "cell_summary", cell_summary,
            chunks=(n_cells, 1), dtype=np.float64,
            compressor=BLOSC_LZ4, overwrite=True,
        )
        cs_arr.attrs["column_names"] = [
            "cell_centroid_x", "cell_centroid_y", "cell_area",
            "nucleus_centroid_x", "nucleus_centroid_y",
            "nucleus_area", "z_level", "nucleus_count",
        ]
        cs_arr.attrs["column_descriptions"] = [
            "Cell centroid in X", "Cell centroid in Y", "Cell area",
            "Nucleus centroid in X", "Nucleus centroid in Y",
            "Nucleus area", "z_level", "Nucleus count",
        ]

        # masks
        masks_grp = root.require_group("masks")
        masks_grp.array(
            "homogeneous_transform", transform,
            chunks=(4, 4), dtype=np.float32,
            compressor=BLOSC_LZ4, overwrite=True,
        )
        for set_idx in (0, 1):
            masks_grp.array(
                str(set_idx), mask,
                chunks=(MASK_CHUNK, MASK_CHUNK), dtype=np.uint32,
                compressor=BLOSC_ZSTD,
                overwrite=True,
            )

        # polygon_sets
        ps_grp = root.require_group("polygon_sets")
        for set_idx in (0, 1):
            grp = ps_grp.require_group(str(set_idx))
            grp.array(
                "vertices", vertices,
                chunks=(n_cells, N_VERTS_MAX), dtype=np.float32,
                compressor=BLOSC_LZ4, overwrite=True,
            )
            grp.array(
                "num_vertices", num_verts,
                chunks=(n_cells,), dtype=np.int32,
                compressor=BLOSC_LZ4, overwrite=True,
            )
            grp.array(
                "cell_index", cell_index,
                chunks=(n_cells,), dtype=np.int32,
                compressor=BLOSC_LZ4, overwrite=True,
            )
            method_val = 3 if set_idx == 0 else 0  # 3=cell boundary, 0=nucleus boundary
            grp.array(
                "method",
                np.full(n_cells, method_val, dtype=np.uint32),
                chunks=(n_cells,), dtype=np.uint32,
                compressor=BLOSC_LZ4, overwrite=True,
            )

        # 序列化 MemoryStore → ZIP bytes
        buf = io.BytesIO()
        with zipfile.ZipFile(buf, "w", compression=zipfile.ZIP_STORED) as zw:
            for key, data in mem.items():
                zw.writestr(key, bytes(data))

        # 原地取代（同一 filesystem → os.replace 接近原子性）
        tmp_path = cells_zarr_path.parent / "_cells_v4_rebuild.zarr.zip"
        tmp_path.write_bytes(buf.getvalue())
        os.replace(str(tmp_path), str(cells_zarr_path))
        logger.info(
            f"cells.zarr.zip 已重建為 Xenium Explorer 4.x 格式："
            f"{n_cells} 個細胞，mask {H}×{W}，"
            f"centroid_x=[{cx.min():.1f}, {cx.max():.1f}] µm"
        )

    # ── 保留舊方法作為備用（僅處理舊格式 polygon_vertices） ────────────────
    @staticmethod
    def _patch_cell_summary_centroids(out_xenium_dir: Path) -> None:
        """
        （已棄用）修補舊格式 cells.zarr.zip 的 cell_summary 中心點。
        新流程改用 _rebuild_cells_zarr_v4 直接寫出 Xenium Explorer 4.x 格式。
        保留此方法供 legacy bundle 診斷使用。
        """
        import os
        import json
        import zipfile
        import zarr
        import numpy as np
        from zarr.storage import MemoryStore

        cells_zarr_path = out_xenium_dir / "cells.zarr.zip"
        if not cells_zarr_path.exists():
            logger.warning("找不到 cells.zarr.zip，跳過中心點修補。")
            return

        tmp_path = cells_zarr_path.parent / "_cells_centroid_patch.zarr.zip"

        try:
            with zipfile.ZipFile(cells_zarr_path, "r") as zr:
                all_entries: dict[str, bytes] = {n: zr.read(n) for n in zr.namelist()}

            # 舊格式才有 polygon_vertices
            if "polygon_vertices/.zarray" not in all_entries:
                logger.info("cells.zarr.zip 已為新格式（polygon_sets），跳過舊格式中心點修補。")
                return

            store_r = zarr.storage.ZipStore(str(cells_zarr_path), mode="r")
            root_r  = zarr.open(store_r, mode="r")
            verts     = root_r["polygon_vertices"][:]
            num_verts = root_r["polygon_num_vertices"][:]
            cs        = root_r["cell_summary"][:]
            store_r.close()

            coords = verts[0]
            n_arr  = num_verts[0]
            x_vals = coords[:, 0::2]
            y_vals = coords[:, 1::2]
            n_cells = len(n_arr)
            x_c = np.array(
                [x_vals[i, :int(n_arr[i])].mean() if n_arr[i] > 0 else 0.0
                 for i in range(n_cells)], dtype=np.float64,
            )
            y_c = np.array(
                [y_vals[i, :int(n_arr[i])].mean() if n_arr[i] > 0 else 0.0
                 for i in range(n_cells)], dtype=np.float64,
            )
            cs8 = np.zeros((n_cells, 8), dtype=np.float64)
            cs8[:, :cs.shape[1]] = cs
            cs8[:, 0] = x_c
            cs8[:, 1] = y_c
            cs8[:, 3] = x_c
            cs8[:, 4] = y_c
            cs8[:, 7] = 1.0

            cs_zarray_key = "cell_summary/.zarray"
            cs_zattrs_key = "cell_summary/.zattrs"
            old_cs_meta   = json.loads(all_entries[cs_zarray_key])
            new_cs_zarray = dict(old_cs_meta)
            new_cs_zarray["shape"]  = [n_cells, 8]
            new_cs_zarray["chunks"] = [n_cells, 1]
            new_cs_zattrs = {
                "column_names": [
                    "cell_centroid_x", "cell_centroid_y", "cell_area",
                    "nucleus_centroid_x", "nucleus_centroid_y",
                    "nucleus_area", "z_level", "nucleus_count",
                ],
                "column_descriptions": [
                    "Cell centroid in X", "Cell centroid in Y", "Cell area",
                    "Nucleus centroid in X", "Nucleus centroid in Y",
                    "Nucleus area", "z_level", "Nucleus count",
                ],
            }

            import numcodecs as _nc
            mem_store = MemoryStore()
            comp = _nc.get_codec(new_cs_zarray["compressor"]) if new_cs_zarray.get("compressor") else None
            cs_arr_mem = zarr.open_array(
                mem_store, path="cell_summary", mode="w",
                shape=(n_cells, 8), dtype=np.float64,
                chunks=(n_cells, 1), compressor=comp, order="C",
            )
            cs_arr_mem[:] = cs8

            cs_prefix = "cell_summary/"
            new_chunks: dict[str, bytes] = {
                k: bytes(v) for k, v in mem_store.items()
                if k.startswith(cs_prefix) and not k.split("/")[-1].startswith(".")
            }
            old_cs_chunk_keys = {
                k for k in all_entries
                if k.startswith(cs_prefix) and not k.split("/")[-1].startswith(".")
            }

            with zipfile.ZipFile(tmp_path, "w", compression=zipfile.ZIP_STORED) as zw:
                for name, data in all_entries.items():
                    if name == cs_zarray_key:
                        zw.writestr(name, json.dumps(new_cs_zarray).encode())
                    elif name == cs_zattrs_key:
                        zw.writestr(name, json.dumps(new_cs_zattrs).encode())
                    elif name in old_cs_chunk_keys:
                        pass
                    else:
                        zw.writestr(name, data)
                for k, v in new_chunks.items():
                    zw.writestr(k, v)

            os.replace(str(tmp_path), str(cells_zarr_path))
            logger.info(
                f"cell_summary 修補完成（8 欄）："
                f"x=[{x_c.min():.1f}, {x_c.max():.1f}] µm，"
                f"y=[{y_c.min():.1f}, {y_c.max():.1f}] µm"
            )

        except Exception as exc:
            logger.warning(f"cell_summary centroid patch 失敗（繼續執行）：{exc}")
            if tmp_path.exists():
                try:
                    tmp_path.unlink()
                except OSError:
                    pass

    @staticmethod
    def _patch_seg_mask_value(out_xenium_dir: Path) -> None:
        """
        修補 cells.zarr.zip 的 seg_mask_value 為 0-indexed [0, 1, ..., N-1]。

        spatialdata_xenium_explorer bug：寫出 1-indexed [1, 2, ..., N]，
        而 cell_id 為 0-indexed [0, 1, ..., N-1]。兩者不一致會導致
        Xenium Explorer "Loading cells..." 永遠卡住。
        """
        import os, json, zipfile
        import numpy as np
        import numcodecs

        cells_path = out_xenium_dir / "cells.zarr.zip"
        if not cells_path.exists():
            logger.warning("找不到 cells.zarr.zip，跳過 seg_mask_value 修補。")
            return

        tmp_path = cells_path.parent / "_cells_seg_patch.zarr.zip"
        try:
            with zipfile.ZipFile(str(cells_path), "r") as zr:
                all_entries: dict[str, bytes] = {n: zr.read(n) for n in zr.namelist()}

            meta = json.loads(all_entries["seg_mask_value/.zarray"])
            n_cells = meta["shape"][0]
            chunk_size = meta["chunks"][0]
            comp_info = meta.get("compressor")
            comp = numcodecs.get_codec(comp_info) if comp_info else None

            # 確認是否已是 0-indexed（避免重複修補）
            first_chunk_raw = all_entries.get("seg_mask_value/0", b"")
            if first_chunk_raw:
                first_val = np.frombuffer(
                    comp.decode(first_chunk_raw) if comp else first_chunk_raw, dtype=np.uint32
                )[0]
                if first_val == 0:
                    logger.info("seg_mask_value 已為 0-indexed，跳過修補。")
                    return

            new_arr = np.arange(0, n_cells, dtype=np.uint32)
            new_chunks: dict[str, bytes] = {}
            for i in range(0, n_cells, chunk_size):
                chunk_data = new_arr[i : i + chunk_size]
                encoded = comp.encode(chunk_data.tobytes()) if comp else chunk_data.tobytes()
                new_chunks[f"seg_mask_value/{i // chunk_size}"] = encoded

            old_keys = {
                k for k in all_entries
                if k.startswith("seg_mask_value/") and not k.split("/")[-1].startswith(".")
            }

            with zipfile.ZipFile(str(tmp_path), "w", compression=zipfile.ZIP_STORED) as zw:
                for name, data in all_entries.items():
                    if name not in old_keys:
                        zw.writestr(name, data)
                for k, v in new_chunks.items():
                    zw.writestr(k, v)

            os.replace(str(tmp_path), str(cells_path))
            logger.info(f"seg_mask_value 修補完成：[0..{n_cells - 1}]（0-indexed）")

        except Exception as exc:
            logger.warning(f"seg_mask_value patch 失敗（繼續執行）：{exc}")
            if tmp_path.exists():
                try:
                    tmp_path.unlink()
                except OSError:
                    pass

    @staticmethod
    def _patch_analysis_indptr(out_xenium_dir: Path) -> None:
        """
        修補 analysis.zarr.zip 中每個 cell_groups 的 indptr（CSR 格式缺少結尾 entry）。

        spatialdata_xenium_explorer bug：每個 grouping 的 indptr 只寫 N 個值，
        標準 CSR 格式需要 N+1 個（最後一個值 = len(indices)），
        否則 Xenium Explorer 無法解析最後一個群組，annotation 完全不顯示。

        修補策略：重建 analysis.zarr.zip，對每個 grouping 追加 indptr 結尾值。
        使用同一目錄下的臨時檔案，以 os.replace 原地取代，避免跨 filesystem 複製風險。
        """
        import os
        import zarr
        import numpy as np

        analysis_path = out_xenium_dir / "analysis.zarr.zip"
        if not analysis_path.exists():
            logger.warning("找不到 analysis.zarr.zip，跳過 indptr 修補。")
            return

        # 臨時檔案放在與 analysis.zarr.zip 相同的目錄（同一 filesystem），
        # 讓 os.replace 可在同一 ExFAT volume 上原地替換，避免跨 filesystem 複製。
        tmp_path = analysis_path.parent / "_analysis_indptr_patch.zarr.zip"
        store_r = store_w = store_v = None
        try:
            # 1. 讀取現有資料
            store_r = zarr.storage.ZipStore(str(analysis_path), mode="r")
            root_r  = zarr.open(store_r, mode="r")
            cg_attrs    = dict(root_r["cell_groups"].attrs)
            n_groupings = cg_attrs["number_groupings"]
            group_names = cg_attrs["group_names"]   # List[List[str]]

            groupings_data: list[dict] = []
            for i in range(n_groupings):
                grp     = root_r["cell_groups"][str(i)]
                indices = grp["indices"][:]
                indptr  = grp["indptr"][:]
                # 若 indptr 少了結尾 entry（N 個而非 N+1 個），補上 len(indices)
                # ⚠️ 必須保留原始 dtype（uint32），np.append 預設升型為 int64。
                # Xenium Explorer 預期 uint32，讀到 int64 會永遠卡住 "Loading cells..."。
                if len(indptr) == len(group_names[i]):
                    indptr = np.append(indptr, np.array([len(indices)], dtype=indptr.dtype))
                groupings_data.append({"indices": indices, "indptr": indptr})

        finally:
            if store_r is not None:
                store_r.close()

        try:
            # 2. 重建 zarr 至臨時檔案（同目錄，同 filesystem）
            store_w = zarr.storage.ZipStore(str(tmp_path), mode="w")
            root_w  = zarr.open(store_w, mode="w")
            root_w.require_group("cell_groups")
            root_w["cell_groups"].attrs.update(cg_attrs)
            for i, gd in enumerate(groupings_data):
                grp_w = root_w["cell_groups"].require_group(str(i))
                grp_w.array("indices", gd["indices"], dtype=gd["indices"].dtype, overwrite=True)
                grp_w.array("indptr",  gd["indptr"],  dtype=gd["indptr"].dtype,  overwrite=True)
        finally:
            if store_w is not None:
                store_w.close()

        # 3. 原地取代（ExFAT 同目錄，盡可能原子性）
        try:
            os.replace(str(tmp_path), str(analysis_path))
        except OSError as exc:
            logger.warning(f"analysis.zarr.zip 取代原始檔案失敗：{exc}")
            raise
        finally:
            # 若 os.replace 失敗，tmp_path 尚存；成功則已消耗
            if tmp_path.exists():
                try:
                    tmp_path.unlink()
                except OSError:
                    pass

        # 4. 驗證
        store_v = None
        try:
            store_v = zarr.storage.ZipStore(str(analysis_path), mode="r")
            root_v  = zarr.open(store_v, mode="r")
            gn_key  = "grouping_names"
            fixed = [
                f"{cg_attrs.get(gn_key, ['?'] * n_groupings)[i]}({len(group_names[i])}g,"
                f"indptr={len(root_v['cell_groups'][str(i)]['indptr'])})"
                for i in range(n_groupings)
            ]
            logger.info(f"analysis.zarr.zip indptr 修補完成：{', '.join(fixed)}")
        except Exception as exc:
            logger.warning(f"indptr 修補後驗證失敗（不影響主要結果）：{exc}")
        finally:
            if store_v is not None:
                store_v.close()

    def _write_xenium_bundle(
        self,
        out_xenium_dir: Path,
        sd_image: Optional[object],
        sd_shapes: Optional[gpd.GeoDataFrame],
        sd_table: Optional[object],
        sd_points: Optional[object],
    ) -> None:
        """組裝 SpatialData 並呼叫 spatialdata_xenium_explorer.write()。"""
        import spatialdata_xenium_explorer

        # 修補 to_intrinsic：讓 element=None 時直接回傳 None，
        # 避免無轉錄點時觸發 AssertionError（library bug）。
        _apply_xenium_explorer_patches()

        # spatialdata_xenium_explorer.write() 硬性要求至少一張影像；
        # 未提供 zarr_path 時以佔位空白影像滿足需求。
        if sd_image is None:
            logger.info("未偵測到影像層，自動建立佔位影像以滿足 Xenium Explorer 需求。")
            sd_image = self._create_dummy_image(sd_shapes)

        # 確保影像至少有一個名為 'global' 的座標系轉換
        # 如果 sd_image 已經有轉換，則維持原樣

        # 相容 spatialdata 新舊 API（tables vs table 參數）
        sdata_kwargs: dict = {}
        sd_init_sig = inspect.signature(sd.SpatialData.__init__)
        if sd_table is not None:
            if "tables" in sd_init_sig.parameters:
                sdata_kwargs["tables"] = {"table": sd_table}
            else:
                sdata_kwargs["table"] = sd_table

        if sd_image is not None:
            sdata_kwargs["images"] = {"tissue_image": sd_image}
        if sd_shapes is not None:
            sdata_kwargs["shapes"] = {"cell_boundaries": sd_shapes}
        if sd_points is not None:
            sdata_kwargs["points"] = {"transcripts": sd_points}

        logger.info(f"組裝 SpatialData 物件 (可用座標系: {['global']})")
        sdata = sd.SpatialData(**sdata_kwargs)

        # 使用與影像 Scale 相同的 pixel_size，避免 library 重採樣不一致。
        # 原則：write(pixel_size=X) 與影像 Scale([X, X]) 必須相同，
        # 對應 Proseg-Zarr-Integration 的做法（pixel_size=PROSEG_SCALE_UM_PX）。
        out_pixel_size = self.pixel_size_um
        logger.info(f"寫出 Xenium Explorer bundle 至：{out_xenium_dir}（pixel_size={out_pixel_size}，大型資料集耗時較長）")
        spatialdata_xenium_explorer.write(
            path=str(out_xenium_dir),
            sdata=sdata,
            image_key="tissue_image" if sd_image is not None else None,
            shapes_key="cell_boundaries" if sd_shapes is not None else None,
            points_key="transcripts" if sd_points is not None else None,
            gene_column="gene" if sd_points is not None else None,
            pixel_size=out_pixel_size,
            lazy=True,
            ram_threshold_gb=self.ram_threshold_gb,
        )
        logger.info("spatialdata_xenium_explorer.write() 完成。")

        # 6. 修補 experiment.xenium pixel_size Bug
        # spatialdata_xenium_explorer 有已知 bug：write() 寫出的 experiment.xenium
        # pixel_size 可能不正確（不等於傳入的 pixel_size 參數）。
        # 強制覆寫確保 Xenium Explorer 使用正確比例。
        self._patch_experiment_xenium(out_xenium_dir, out_pixel_size)

        # 7. 修補 cells.zarr.zip 的 cell_summary 中心點 + seg_mask_value（0-indexed）
        # spatialdata_xenium_explorer 寫出的是舊格式（major_version=5, polygon_vertices）。
        # 此格式與 Xenium Explorer 4.1.1 相容（比 v6 新格式的 2D mask 載入快得多）。
        # ⚠️ 不呼叫 _rebuild_cells_zarr_v4——新格式需載入整張 2D 遮罩，在 Xenium Explorer
        #    4.1.1 中反而導致 "Loading cells..." 卡頓數分鐘。
        # ⚠️ seg_mask_value 必須 0-indexed（與 cell_id 一致），否則 Xenium Explorer 卡住。
        #    library 寫出 1-indexed [1..N]，此修補改為 [0..N-1]。
        self._patch_cell_summary_centroids(out_xenium_dir)
        self._patch_seg_mask_value(out_xenium_dir)

        # 8. analysis.zarr.zip 的 indptr 格式說明（不需修補）
        # Xenium Explorer 期待每個 grouping 的 indptr 有 N_groups 個 start positions：
        #   indptr[i] = group i 在 indices 中的起始位置；最後一個 group 隱含到 len(indices)。
        # spatialdata_xenium_explorer 的輸出格式正確（N_groups entries, dtype=uint32）。
        # ⚠️ 不要呼叫 _patch_analysis_indptr——它原本是錯誤假設（以為需要 N+1 entries），
        #    會多加一個 trailing entry 導致 annotation 完全不顯示。
        # self._patch_analysis_indptr(out_xenium_dir)  # 已停用

        # 9. 建立 analysis_summary.html（Xenium Explorer 4.x 要求此檔存在）
        # 缺少時 Xenium Explorer 會持續等待，卡在 "Loading cells..." 畫面。
        summary_html = out_xenium_dir / "analysis_summary.html"
        if not summary_html.exists():
            summary_html.write_text(
                "<!DOCTYPE html><html><body><p>MCseg v2 Analysis</p></body></html>",
                encoding="utf-8",
            )
            logger.info("analysis_summary.html 建立完成。")


def generate_combined_geojson(
    tile_proseg_dir: "Path",
    zarr_path: "Path" = None,
    config: dict = None,
) -> dict:
    """
    將各 tile 的 proseg_results.json（gzip GeoJSON）合併為 ROI 絕對座標的 GeoJSON。

    座標換算說明
    -----------
    Proseg 每個 tile 的 GeoJSON 座標為**局部座標**（相對於 padded tile 起點），
    需加上正確的全域偏移才能得到 ROI µm 座標。

    正確偏移公式（考慮 padding）：
        x_start_px = max(0, ix × tile_w_px - padding_px)
        x_offset_µm = x_start_px × coordinate_scale_um_px

    原始錯誤公式的兩個 bug：
        1. 使用 PROSEG_UM_PX (0.2645833) 而非 ROI coordinate_scale (0.2737)
        2. 未減去 padding，導致偏移量虛高且座標超出影像範圍

    Parameters
    ----------
    tile_proseg_dir:
        ``results/analysis/roi/{roi_name}/proseg_tiles`` 目錄路徑。
    zarr_path:
        Zarr 路徑，用於讀取 label 尺寸（可選；若 None 則從 tile 目錄推斷）。
    config:
        Pipeline config dict，讀取 proseg.tiling 與 ROI pixel_size_um。

    Returns
    -------
    dict
        GeoJSON FeatureCollection（可直接 json.dump）。
    """
    import gzip
    import json as _json
    import re
    import zarr as _zarr

    # 從 config 取得分塊參數
    cfg = config or {}
    tiling = cfg.get("proseg", {}).get("tiling", {})
    grid_nx: int = int(tiling.get("grid_nx", 4))
    grid_ny: int = int(tiling.get("grid_ny", 3))
    padding: int = int(tiling.get("padding", 200))

    # coordinate_scale：優先用 ROI pixel_size_um，退回 VISIUM_UM_PX
    rois = cfg.get("rois", [{}])
    coordinate_scale: float = (
        rois[0].get("pixel_size_um", VISIUM_UM_PX) if rois else VISIUM_UM_PX
    )

    # 從 Zarr 取得完整 label 尺寸（決定 tile 邊界）
    w_full: Optional[int] = None
    h_full: Optional[int] = None
    if zarr_path is not None:
        try:
            z = _zarr.open(str(zarr_path), mode="r")
            label_arr = z["labels"]["cellpose_nuclei"]["0"]
            h_full, w_full = int(label_arr.shape[-2]), int(label_arr.shape[-1])
        except Exception as exc:
            logger.warning(f"無法讀取 zarr label 尺寸（{exc}），使用 tile 目錄命名推斷。")

    _pattern = re.compile(r"tile_y(\d+)_x(\d+)$")
    all_features: list = []

    tile_dirs = sorted(
        [d for d in Path(tile_proseg_dir).iterdir()
         if d.is_dir() and _pattern.match(d.name)],
        key=lambda d: d.name,
    )

    for tile_dir in tile_dirs:
        json_path = tile_dir / "proseg_results.json"
        if not json_path.exists():
            logger.warning(f"缺少 proseg_results.json：{tile_dir.name}，略過")
            continue

        m = _pattern.match(tile_dir.name)
        iy, ix = int(m.group(1)), int(m.group(2))

        # 計算 tile 像素邊界（使用與 runner.py 相同的分塊邏輯）
        tile_w = (w_full // grid_nx) if w_full else 372
        tile_h = (h_full // grid_ny) if h_full else 403
        # padded_start = 去掉 padding 的全域像素起點
        x_start_px = max(0, ix * tile_w - padding)
        y_start_px = max(0, iy * tile_h - padding)
        x_offset_um = x_start_px * coordinate_scale
        y_offset_um = y_start_px * coordinate_scale

        # 名義 tile 邊界（µm）—— 用於剪裁超出 tile 範圍的重疊多邊形
        # Proseg 因 mask padding 會讓多邊形延伸進相鄰 tile，造成重複細胞
        x_tile_min_um = ix * tile_w * coordinate_scale
        x_tile_max_um = min((ix + 1) * tile_w, w_full if w_full else (ix + 1) * tile_w) * coordinate_scale
        y_tile_min_um = iy * tile_h * coordinate_scale
        y_tile_max_um = min((iy + 1) * tile_h, h_full if h_full else (iy + 1) * tile_h) * coordinate_scale

        try:
            with gzip.open(str(json_path), "rt") as f:
                data = _json.load(f)
        except Exception as exc:
            logger.warning(f"讀取 {tile_dir.name} GeoJSON 失敗：{exc}，略過")
            continue

        skipped = 0
        for feat in data.get("features", []):
            cell_idx = feat["properties"].get("cell", 0)
            full_id = f"{tile_dir.name}_{cell_idx}"
            geom = feat["geometry"]
            geom_type = geom["type"]

            if geom_type == "MultiPolygon":
                new_coords = [
                    [
                        [[c[0] + x_offset_um, c[1] + y_offset_um] for c in ring]
                        for ring in poly
                    ]
                    for poly in geom["coordinates"]
                ]
                new_geom = {"type": "MultiPolygon", "coordinates": new_coords}
                # 用第一個 polygon 的外環計算重心
                flat = new_coords[0][0]
            elif geom_type == "Polygon":
                new_coords = [
                    [[c[0] + x_offset_um, c[1] + y_offset_um] for c in ring]
                    for ring in geom["coordinates"]
                ]
                new_geom = {"type": "Polygon", "coordinates": new_coords}
                flat = new_coords[0]
            else:
                continue

            # 邊界剪裁：重心必須落在此 tile 的名義範圍內
            # 避免因 mask padding 導致相鄰 tile 出現重複多邊形
            cx = sum(c[0] for c in flat) / len(flat)
            cy = sum(c[1] for c in flat) / len(flat)
            if not (x_tile_min_um <= cx < x_tile_max_um and
                    y_tile_min_um <= cy < y_tile_max_um):
                skipped += 1
                continue

            all_features.append({
                "type": "Feature",
                "properties": {"full_id": full_id},
                "geometry": new_geom,
            })

        if skipped:
            logger.debug(f"  {tile_dir.name}: 剪裁 {skipped} 個超出名義邊界的多邊形")

    logger.info(
        f"合併 {len(tile_dirs)} 個 tile，共 {len(all_features)} 個多邊形。"
        f"（coordinate_scale={coordinate_scale}, padding={padding}px）"
    )
    return {"type": "FeatureCollection", "features": all_features}
