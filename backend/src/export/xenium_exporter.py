"""
Stage 5 — Xenium Explorer 匯出模組

將 Proseg 分割結果（h5ad + GeoJSON 多邊形 + 轉錄點 CSV）組裝成
spatialdata_xenium_explorer 可讀取的 Xenium Explorer bundle。

關鍵 Bug 修復
--------------
spatialdata_xenium_explorer.write() 內部硬編碼 pixel_size=0.2125，
每次寫出後必須修補 experiment.xenium 的 pixel_size 欄位為 PROSEG_UM_PX。

座標系對應
----------
- GeoJSON 多邊形：µm（Proseg 輸出）
- 轉錄點 CSV：µm（combined_transcripts.csv）
- Xenium Explorer 輸入：像素（÷ PROSEG_UM_PX = ÷ 0.2645833）
- 影像（Zarr）：像素，直接使用

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
from shapely.affinity import scale as shapely_scale
from spatialdata.models import (
    Image2DModel,
    PointsModel,
    ShapesModel,
    TableModel,
)
from spatialdata.transformations import Identity

from backend.src.utils.constants import PROSEG_UM_PX, PROSEG_NM_PX

logger = logging.getLogger("pipeline.export.xenium")


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
    """

    def __init__(
        self,
        zarr_path: Optional[str | Path] = None,
        poly_json_path: Optional[str | Path] = None,
        transcripts_csv_path: Optional[str | Path] = None,
        ram_threshold_gb: float = 4.0,
    ) -> None:
        self.zarr_path = Path(zarr_path) if zarr_path else None
        self.poly_json_path = Path(poly_json_path) if poly_json_path else None
        self.transcripts_csv_path = (
            Path(transcripts_csv_path) if transcripts_csv_path else None
        )
        self.ram_threshold_gb = ram_threshold_gb

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

        # 4. 載入轉錄點
        sd_points = self._load_transcripts(id_remap)

        # 5. 組裝 SpatialData 並寫出
        self._write_xenium_bundle(
            out_xenium_dir=out_xenium_dir,
            sd_image=sd_image,
            sd_shapes=sd_shapes,
            sd_table=sd_table,
            sd_points=sd_points,
        )

        # 6. 修補 experiment.xenium pixel_size Bug
        self._patch_experiment_xenium(out_xenium_dir)

        logger.info(f"=== Xenium Explorer bundle 完成：{out_xenium_dir} ===")
        return out_xenium_dir

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    def _load_image(self) -> Optional[sd.models.Image2DModel]:
        """從 Zarr 惰性載入影像（失敗時回傳 None，不中斷流程）。"""
        if self.zarr_path is None or not self.zarr_path.exists():
            logger.info("未提供 zarr_path，跳過影像層。")
            return None

        try:
            import dask.array as da
            import zarr

            logger.info(f"載入 Zarr 影像：{self.zarr_path}")
            z = zarr.open(str(self.zarr_path), mode="r")
            img_array = z["images"]["tissue_hires_image"]["0"]

            # 統一為 (c, y, x) 格式
            if img_array.shape[0] in (1, 3):
                dask_img = da.from_zarr(img_array)
            else:
                dask_img = da.from_zarr(img_array).transpose((2, 0, 1))

            sd_image = Image2DModel.parse(
                dask_img, dims=("c", "y", "x"),
                transformations={"global": Identity()},
            )
            logger.info(f"影像載入完成，shape={dask_img.shape}")
            return sd_image

        except Exception as exc:
            logger.warning(f"Zarr 影像載入失敗（繼續執行）：{exc}")
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

        features = geo_data.get("features", [])
        polygons: list = []
        valid_cell_ids: list[str] = []

        obs_set = set(adata.obs_names)

        for feat in features:
            props = feat.get("properties", {})
            full_id = props.get("full_id")

            # 只保留通過 QC 的細胞
            if full_id not in obs_set:
                continue

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

                # µm → px
                poly_px = shapely_scale(
                    poly_um,
                    xfact=1.0 / PROSEG_UM_PX,
                    yfact=1.0 / PROSEG_UM_PX,
                    origin=(0, 0),
                )

                # 平滑化（消除 watershed 鋸齒）
                poly_px = self._smooth_polygon(poly_px)

                # Xenium Explorer 要求全部為 MultiPolygon
                if poly_px.geom_type == "Polygon":
                    poly_px = shapely.geometry.MultiPolygon([poly_px])

                polygons.append(poly_px)
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
        消除 Proseg watershed 分割產生的鋸齒邊緣。

        步驟：
        1. simplify(0.4)  — 移除微小變異
        2. buffer(+1.0)   — 對外膨脹，去除尖角
        3. buffer(-1.0)   — 收縮回原大小，整平邊緣
        """
        try:
            smooth = poly.simplify(0.4, preserve_topology=True)
            smooth = smooth.buffer(1.0, join_style=1).buffer(-1.0, join_style=1)
            if (
                smooth.is_valid
                and not smooth.is_empty
                and smooth.geom_type in ("Polygon", "MultiPolygon")
            ):
                return smooth
        except Exception:
            pass
        return poly  # 失敗時回退到原始多邊形

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

            # µm → px
            tx_dd["x"] = tx_dd["x"] / PROSEG_UM_PX
            tx_dd["y"] = tx_dd["y"] / PROSEG_UM_PX

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

            sd_points = PointsModel.parse(tx_dd, sort=True, transformations={"global": Identity()}, **parse_kwargs)
            logger.info("轉錄點 PointsModel 建立完成。")
            return sd_points

        except Exception as exc:
            logger.error(f"轉錄點載入失敗（繼續執行）：{exc}")
            return None

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

        logger.info("組裝 SpatialData 物件...")
        sdata = sd.SpatialData(**sdata_kwargs)

        logger.info(f"寫出 Xenium Explorer bundle 至：{out_xenium_dir}（大型資料集耗時較長）")
        spatialdata_xenium_explorer.write(
            path=str(out_xenium_dir),
            sdata=sdata,
            image_key="tissue_image" if sd_image is not None else None,
            shapes_key="cell_boundaries" if sd_shapes is not None else None,
            points_key="transcripts" if sd_points is not None else None,
            gene_column="gene" if sd_points is not None else None,
            pixel_size=PROSEG_UM_PX,
            lazy=True,
            ram_threshold_gb=self.ram_threshold_gb,
        )
        logger.info("spatialdata_xenium_explorer.write() 完成。")

    @staticmethod
    def _patch_experiment_xenium(out_xenium_dir: Path) -> None:
        """
        修補 experiment.xenium 的 pixel_size Bug。

        spatialdata_xenium_explorer 內部硬編碼 pixel_size=0.2125（Xenium 原生值），
        即使呼叫時傳入正確的 PROSEG_UM_PX，寫出後仍會被覆蓋。
        此方法在寫出後強制將 pixel_size 改寫為正確的 PROSEG_UM_PX。
        """
        exp_file = out_xenium_dir / "experiment.xenium"
        if not exp_file.exists():
            logger.warning(f"experiment.xenium 不存在，跳過修補：{exp_file}")
            return

        try:
            with open(exp_file, "r") as f:
                exp_data = json.load(f)

            old_val = exp_data.get("pixel_size", "N/A")
            exp_data["pixel_size"] = PROSEG_UM_PX  # 0.2645833（修補硬編碼 0.2125）

            with open(exp_file, "w") as f:
                json.dump(exp_data, f, indent=2)

            logger.info(
                f"experiment.xenium pixel_size 修補完成：{old_val} → {PROSEG_UM_PX}"
            )
        except Exception as exc:
            logger.error(f"experiment.xenium 修補失敗：{exc}")


def generate_combined_geojson(
    tile_proseg_dir: "Path",
    zarr_path: "Path",
    config: dict,
) -> dict:
    """
    將各 tile 的 proseg_results.json（gzip GeoJSON）合併為 ROI 絕對座標的 GeoJSON。

    座標轉換：
        abs_µm = tile_rel_µm + (ix * tile_w_px * scale_um_px,
                                 iy * tile_h_px * scale_um_px)

    每個 feature 的 full_id 設為 ``tile_y{iy}_x{ix}_{cell_idx}``，
    與 proseg_cells.h5ad 的 obs_names 格式一致。

    Parameters
    ----------
    tile_proseg_dir:
        ``results/analysis/roi/{roi_name}/proseg_tiles`` 目錄路徑。
    zarr_path:
        ``proseg_integrated.zarr`` 路徑，用于讀取 label 尺寸。
    config:
        Pipeline config dict，需含 proseg.tiling 與 proseg.constants。

    Returns
    -------
    dict
        GeoJSON FeatureCollection（可直接 json.dump）。
    """
    import gzip
    import json as _json
    import re
    import zarr as _zarr

    tiling = config.get("proseg", {}).get("tiling", {})
    scale_um_px: float = (
        config.get("proseg", {})
        .get("constants", {})
        .get("scale_um_px", PROSEG_UM_PX)
    )
    grid_nx: int = int(tiling.get("grid_nx", 4))
    grid_ny: int = int(tiling.get("grid_ny", 3))

    # 從 Zarr 取得完整 label 尺寸
    z = _zarr.open(str(zarr_path), mode="r")
    label_arr = z["labels"]["cellpose_nuclei"]["0"]
    h_full, w_full = label_arr.shape[-2], label_arr.shape[-1]
    tile_w_px = w_full // grid_nx
    tile_h_px = h_full // grid_ny

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
        x_offset_um = ix * tile_w_px * scale_um_px
        y_offset_um = iy * tile_h_px * scale_um_px

        try:
            with gzip.open(str(json_path), "rt") as f:
                data = _json.load(f)
        except Exception as exc:
            logger.warning(f"讀取 {tile_dir.name} GeoJSON 失敗：{exc}，略過")
            continue

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
            elif geom_type == "Polygon":
                new_coords = [
                    [[c[0] + x_offset_um, c[1] + y_offset_um] for c in ring]
                    for ring in geom["coordinates"]
                ]
                new_geom = {"type": "Polygon", "coordinates": new_coords}
            else:
                continue

            all_features.append({
                "type": "Feature",
                "properties": {"full_id": full_id},
                "geometry": new_geom,
            })

    logger.info(
        f"合併 {len(tile_dirs)} 個 tile，共 {len(all_features)} 個多邊形。"
    )
    return {"type": "FeatureCollection", "features": all_features}
