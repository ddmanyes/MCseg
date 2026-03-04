"""
Stage 2: Zarr 建構模組
整合 H&E 影像、細胞核遮罩、Visium HD binned matrix 為 SpatialData OME-Zarr 格式

移植自 visiumhd_pipeline/scripts/02_build_zarr/create_zarr.py
並整合 visiumhd_pipeline/src/proseg/zarr_handler.py 的 dask monkey patch 守則。
"""

import gc
import json
import logging
import os
import shutil
import subprocess
import warnings
from pathlib import Path
from typing import Any, Optional, Tuple

import dask
dask.config.set({"dataframe.query-planning": True})

import dask.array as da
import numpy as np

from backend.src.utils.config import resolve_path
from backend.src.utils.constants import VISIUM_UM_PX

logger = logging.getLogger("pipeline.zarr")

CHUNK_SIZE = 1024


# ---------------------------------------------------------------------------
# macOS ExFAT 相容性：Monkey patch zarr DirectoryStore
# ---------------------------------------------------------------------------
# ExFAT 上每次寫入目錄都會產生 macOS AppleDouble 附屬檔案（._*），
# zarr.convenience.consolidate_metadata 在遍歷 store 時會透過 _keys_fast
# 取得所有檔案 key，並嘗試將 is_zarr_key(.zarray/.zgroup/.zattrs) 的內容
# 以 json_loads 解碼。若 ._* 附屬檔案的名稱恰好符合這些後綴，
# 就會嘗試解碼含有 0xB0 等非 UTF-8 位元組的 AppleDouble 二進位格式，
# 引發 UnicodeDecodeError。
# 最安全的解法是在 _keys_fast 層面直接過濾掉這些檔案。

def _patch_zarr_for_exfat() -> None:
    """Monkey patch zarr.storage.DirectoryStore._keys_fast 以過濾 macOS 垃圾檔案。"""
    try:
        import os as _os
        import zarr.storage as _zs

        @staticmethod  # type: ignore[misc]
        def _filtered_keys_fast(path, walker=_os.walk):
            for dirpath, _, filenames in walker(path):
                dirpath = _os.path.relpath(dirpath, path)
                clean = [
                    f for f in filenames
                    if not f.startswith("._") and f != ".DS_Store"
                ]
                if dirpath == _os.curdir:
                    yield from clean
                else:
                    dirpath = dirpath.replace("\\", "/")
                    for f in clean:
                        yield "/".join((dirpath, f))

        _zs.DirectoryStore._keys_fast = _filtered_keys_fast
        logger.info("✅ zarr DirectoryStore._keys_fast monkey patch 已套用（ExFAT 相容）")
    except Exception as _e:
        logger.warning(f"⚠️  zarr monkey patch 失敗（無影響，但 ExFAT 上可能出現 UTF-8 錯誤）：{_e}")


_patch_zarr_for_exfat()


# ---------------------------------------------------------------------------
# macOS junk cleanup
# ---------------------------------------------------------------------------

def _clean_mac_junk(path: Path) -> None:
    """清理 macOS 系統垃圾檔案（._* 和 .DS_Store）。"""
    count = 0
    for f in path.rglob("._*"):
        f.unlink(missing_ok=True)
        count += 1
    for f in path.rglob(".DS_Store"):
        f.unlink(missing_ok=True)
        count += 1
    if count:
        logger.info(f"清理 {count} 個 macOS 垃圾檔案於 {path}")


# ---------------------------------------------------------------------------
# Zarr structure creation
# ---------------------------------------------------------------------------

def _create_zarr_structure(store_path: Path):
    """
    建立（或重建）SpatialData OME-Zarr 根目錄結構。

    Returns
    -------
    zarr.Group
        根群組物件
    """
    import zarr

    store_str = str(store_path)
    if store_path.exists():
        logger.info(f"移除既有 Zarr：{store_path}")
        try:
            shutil.rmtree(store_str)
        except Exception as e:
            logger.warning(f"shutil.rmtree 失敗（{e}），改用 rm -rf...")
            # macOS 臨時隱藏檔案可能導致 exit 1，忽略即可
            subprocess.run(["rm", "-rf", store_str])

    root = zarr.open_group(store=store_str, mode="w")
    root.create_group("images")
    root.create_group("labels")
    root.create_group("points")
    root.create_group("shapes")
    root.create_group("tables")
    root.attrs["spatialdata_attrs"] = {"version": "0.1"}
    return root


# ---------------------------------------------------------------------------
# H&E image pyramid
# ---------------------------------------------------------------------------

def _write_image_pyramid(root, image_path: Path, pixel_size_override: Optional[float] = None):
    """
    讀取 H&E 影像並寫出多尺度 OME-Zarr 金字塔（儲存於 root/images/tissue_hires_image）。

    Returns
    -------
    (num_levels, pixel_size_dict)
    """
    import cv2
    import tifffile

    image_path_str = str(image_path)
    logger.info(f"處理 H&E 影像：{image_path_str}")

    pixel_size = {"unit": "micrometer", "scale": 1.0}
    if pixel_size_override is not None:
        pixel_size = {"unit": "micrometer", "scale": float(pixel_size_override)}
        logger.info(f"使用設定的 pixel size：{pixel_size['scale']:.4f} µm/px")

    # --- 讀取影像 ---
    img_data = None
    if image_path_str.lower().endswith((".tiff", ".tif", ".btf")):
        logger.info("以 tifffile memmap 讀取 TIFF/BTF...")
        try:
            with tifffile.TiffFile(image_path_str) as tif:
                try:
                    tags = tif.pages[0].tags
                    unit_tag = tags.get(296)
                    x_res_tag = tags.get(282)
                    if unit_tag and x_res_tag and pixel_size_override is None:
                        unit_val = unit_tag.value
                        x_res_val = x_res_tag.value
                        if isinstance(x_res_val, tuple):
                            res_val = x_res_val[0] / x_res_val[1]
                        else:
                            res_val = float(x_res_val)
                        if unit_val == 2:  # Inch
                            scale = 25400.0 / res_val
                            pixel_size = {"unit": "micrometer", "scale": scale}
                            logger.info(f"自動偵測 pixel size（Inch）：{scale:.4f} µm/px")
                        elif unit_val == 3:  # Centimeter
                            scale = 10000.0 / res_val
                            pixel_size = {"unit": "micrometer", "scale": scale}
                            logger.info(f"自動偵測 pixel size（cm）：{scale:.4f} µm/px")
                        else:
                            logger.warning(f"未知解析度單位 tag：{unit_val}，使用預設")
                except Exception as meta_e:
                    logger.warning(f"TIFF metadata 提取失敗：{meta_e}")

            img_data = tifffile.memmap(image_path_str)
        except Exception as e:
            logger.warning(f"tifffile memmap 失敗：{e}，改用 tifffile pages[0]...")
            try:
                with tifffile.TiffFile(image_path_str) as tif:
                    img_data = tif.pages[0].asarray()
                logger.info("tifffile pages[0] 讀取成功")
            except Exception as e2:
                logger.warning(f"tifffile pages[0] 失敗：{e2}，改用 cv2...")
                img_data = cv2.imread(image_path_str)
                if img_data is not None:
                    img_data = cv2.cvtColor(img_data, cv2.COLOR_BGR2RGB)
    else:
        img_data = cv2.imread(image_path_str)
        if img_data is not None:
            img_data = cv2.cvtColor(img_data, cv2.COLOR_BGR2RGB)

    if img_data is None:
        raise ValueError(f"無法讀取影像：{image_path_str}")

    logger.info(f"原始 shape：{img_data.shape}，dtype：{img_data.dtype}")

    # 標準化為 (C, H, W)
    if img_data.ndim == 2:
        img_data = img_data[:, :, np.newaxis]
    if img_data.shape[0] == 3 and img_data.shape[-1] != 3:
        pass  # 已是 (C, H, W)
    else:
        img_data = np.transpose(img_data, (2, 0, 1))

    logger.info(f"標準化後 shape (C, H, W)：{img_data.shape}")

    c_dim = img_data.shape[0]
    darr = da.from_array(img_data, chunks=(1, CHUNK_SIZE, CHUNK_SIZE))

    images_grp = root["images"]
    img_grp = images_grp.create_group("tissue_hires_image", overwrite=True)

    logger.info("寫出 Level 0...")
    z0 = img_grp.create_dataset(
        name="0",
        shape=darr.shape,
        dtype=darr.dtype,
        chunks=(1, CHUNK_SIZE, CHUNK_SIZE),
        overwrite=True,
    )
    da.to_zarr(darr, z0)

    # 金字塔降採樣
    current_data = darr
    level = 1
    while True:
        if current_data.shape[1] < CHUNK_SIZE and current_data.shape[2] < CHUNK_SIZE:
            break
        logger.info(f"生成 Level {level}...")
        current_data = da.coarsen(
            np.mean, current_data, {0: 1, 1: 2, 2: 2}, trim_excess=True
        )
        current_data = current_data.astype(img_data.dtype)
        z_next = img_grp.create_dataset(
            name=str(level),
            shape=current_data.shape,
            dtype=current_data.dtype,
            chunks=(c_dim, CHUNK_SIZE, CHUNK_SIZE),
            overwrite=True,
        )
        da.to_zarr(current_data, z_next)
        level += 1
        if level > 6:
            break

    del darr, current_data
    gc.collect()

    return level, pixel_size


# ---------------------------------------------------------------------------
# OME-NGFF metadata
# ---------------------------------------------------------------------------

def _write_ome_metadata(root, num_levels: int, pixel_size: Optional[dict] = None) -> None:
    """將 OME-NGFF multiscales metadata 寫入 images/tissue_hires_image 群組。"""
    if pixel_size is None:
        pixel_size = {"unit": "micrometer", "scale": 1.0}

    base_scale = pixel_size["scale"]
    unit = pixel_size["unit"]
    logger.info(f"寫出 OME-NGFF metadata（base_scale={base_scale} {unit}）")

    final_zattrs = {
        "multiscales": [
            {
                "version": "0.4",
                "axes": [
                    {"name": "c", "type": "channel"},
                    {"name": "y", "type": "space", "unit": unit},
                    {"name": "x", "type": "space", "unit": unit},
                ],
                "coordinateTransformations": [
                    {
                        "type": "scale",
                        "scale": [1.0, 1.0, 1.0],
                    }
                ],
                "datasets": [
                    {
                        "path": str(i),
                        "coordinateTransformations": [
                            {
                                "type": "scale",
                                "scale": [1.0, base_scale * (2.0 ** i), base_scale * (2.0 ** i)],
                            }
                        ],
                    }
                    for i in range(num_levels)
                ],
            }
        ]
    }

    if "images" in root and "tissue_hires_image" in root["images"]:
        root["images"]["tissue_hires_image"].attrs.update(final_zattrs)
    else:
        logger.warning("images/tissue_hires_image 群組不存在，metadata 無法寫入")


# ---------------------------------------------------------------------------
# Expression table
# ---------------------------------------------------------------------------

def _write_table(root, table_name: str, config: dict, store_path: str, roi_crop=None) -> None:
    """
    載入 Visium HD H5 matrix + tissue_positions，對齊座標後寫入 tables/{table_name}。
    roi_crop: {"x0": int, "y0": int, "x1": int, "y1": int} fullres px，若設定則篩選 ROI 範圍。
    """
    import pandas as pd
    import scanpy as sc
    import zarr

    logger.info(f"處理 Table：{table_name}")

    matrix_path = str(config["matrix"])
    pos_path = str(config["positions"])

    logger.info(f"  載入 H5：{matrix_path}")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        adata = sc.read_10x_h5(matrix_path)

    logger.info(f"  載入 Parquet：{pos_path}")
    df_pos = pd.read_parquet(pos_path)

    if "barcode" in df_pos.columns:
        df_pos = df_pos.set_index("barcode")

    # ROI 篩選（僅保留 crop 範圍內的 spots）
    if roi_crop is not None and "pxl_col_in_fullres" in df_pos.columns:
        x0, y0 = roi_crop["x0"], roi_crop["y0"]
        x1, y1 = roi_crop["x1"], roi_crop["y1"]
        in_roi = (
            (df_pos["pxl_col_in_fullres"] >= x0) & (df_pos["pxl_col_in_fullres"] < x1) &
            (df_pos["pxl_row_in_fullres"] >= y0) & (df_pos["pxl_row_in_fullres"] < y1)
        )
        df_pos = df_pos[in_roi]
        logger.info(f"  ROI 篩選後：{len(df_pos)} spots in x[{x0},{x1}) y[{y0},{y1})")

    adata.var_names_make_unique()
    adata_original = adata.copy()
    common = adata.obs_names.intersection(df_pos.index)

    if len(adata_original) > 0:
        overlap_rate = len(common) / len(adata_original)
        logger.info(f"  Barcode 覆蓋率：{overlap_rate:.2%}")
        if overlap_rate < 0.8:
            logger.warning(f"  低覆蓋率（{overlap_rate:.2%}），請確認資料對齊！")

    adata = adata[common].copy()
    df_pos = df_pos.loc[common]
    logger.info(f"  Cells: {adata.shape[0]}, Genes: {adata.shape[1]}")

    sf = config["scale_factor"]
    logger.info(f"  套用 scale factor：{sf}")

    if "pxl_col_in_fullres" in df_pos.columns:
        coords = df_pos[["pxl_col_in_fullres", "pxl_row_in_fullres"]].values.astype(float)
    else:
        logger.error("  座標欄位不存在，跳過 Table 寫出")
        return

    # ROI 座標偏移（轉為 crop-local 座標）
    if roi_crop is not None:
        coords[:, 0] -= roi_crop["x0"]
        coords[:, 1] -= roi_crop["y0"]

    adata.obsm["spatial"] = coords * sf

    for col in df_pos.columns:
        if col not in ["pxl_col_in_fullres", "pxl_row_in_fullres"]:
            adata.obs[col] = df_pos[col]

    table_path = os.path.join(store_path, "tables", table_name)
    logger.info(f"  寫出 AnnData：{table_path}")

    # Pre-delete existing table dir to avoid ExFAT shutil.rmtree fd-unlink failures
    if os.path.exists(table_path):
        subprocess.run(["rm", "-rf", table_path], check=False)
        logger.info(f"  已清除舊 table：{table_path}")

    adata.uns["spatialdata_attrs"] = {
        "version": "0.1",
        "region": "tissue_hires_image",
        "region_key": "region",
        "instance_key": "instance_id",
    }
    adata.obs["region"] = "tissue_hires_image"
    adata.obs["region"] = adata.obs["region"].astype("category")
    adata.obs["instance_id"] = np.arange(adata.shape[0])

    # Sanitize bytes→str to prevent zarr UTF-8 decode errors
    # (10x H5 files sometimes store gene names/barcodes as raw bytes)
    def _decode_bytes(v):
        if isinstance(v, bytes):
            try:
                return v.decode("utf-8")
            except UnicodeDecodeError:
                return v.decode("latin-1")
        return v

    adata.var_names = pd.Index([_decode_bytes(v) for v in adata.var_names])
    adata.obs_names = pd.Index([_decode_bytes(v) for v in adata.obs_names])
    for _col in list(adata.var.columns):
        if adata.var[_col].dtype == object:
            adata.var[_col] = adata.var[_col].map(_decode_bytes)
    for _col in list(adata.obs.columns):
        if adata.obs[_col].dtype == object:
            adata.obs[_col] = adata.obs[_col].map(_decode_bytes)

    adata.write_zarr(table_path)

    try:
        z_grp = zarr.open_group(table_path, mode="r+")
        z_grp.attrs.update(
            {
                "version": "0.1",
                "region": "tissue_hires_image",
                "region_key": "region",
                "instance_key": "instance_id",
                "spatialdata-encoding-type": "ngff:regions_table",
            }
        )
        logger.info(f"  spatialdata attrs 已寫入 {table_path}")
    except Exception as e:
        logger.warning(f"  無法更新 table group attrs：{e}")

    del adata, df_pos
    gc.collect()


# ---------------------------------------------------------------------------
# Segmentation masks (labels)
# ---------------------------------------------------------------------------

def _add_masks_to_zarr(
    root,
    mask_path: Path,
    physical_scale: float = 1.0,
    label_name: str = "cellpose",
) -> None:
    """
    將 .npy 細胞分割遮罩寫入 Zarr labels/{label_name}（含多尺度金字塔）。
    """
    from numcodecs import Blosc

    mask_path_str = str(mask_path)
    logger.info(f"處理 Mask：{mask_path_str}")

    if not mask_path.exists():
        logger.warning(f"  找不到 Mask 檔案：{mask_path}，跳過")
        return

    if "labels" not in root:
        labels_grp = root.create_group("labels")
    else:
        labels_grp = root["labels"]

    logger.info(f"  載入 Mask：{mask_path_str}")
    try:
        mask_data = np.load(mask_path_str, mmap_mode="r")
    except Exception as e:
        logger.error(f"  載入 Mask 失敗：{e}")
        return

    logger.info(f"  Mask shape：{mask_data.shape}")
    label_data_grp = labels_grp.create_group(label_name, overwrite=True)

    logger.info("  寫出 Mask Level 0...")
    label_data_grp.create_dataset(
        name="0",
        shape=mask_data.shape,
        data=mask_data,
        chunks=(CHUNK_SIZE, CHUNK_SIZE),
        dtype="i4",
        compressor=Blosc(cname="zstd", clevel=5),
        overwrite=True,
    )

    darr = da.from_array(mask_data, chunks=(CHUNK_SIZE, CHUNK_SIZE))
    current_mask = darr
    level = 1

    while True:
        if current_mask.shape[0] < CHUNK_SIZE and current_mask.shape[1] < CHUNK_SIZE:
            break
        logger.info(f"  生成 Mask Level {level}...")
        current_mask = current_mask[::2, ::2]
        z_next = label_data_grp.create_dataset(
            name=str(level),
            shape=current_mask.shape,
            dtype=mask_data.dtype,
            chunks=(CHUNK_SIZE, CHUNK_SIZE),
            compressor=Blosc(cname="zstd", clevel=5),
            overwrite=True,
        )
        da.to_zarr(current_mask, z_next)
        level += 1
        if level > 6:
            break

    datasets_meta = [
        {
            "path": str(i),
            "coordinateTransformations": [
                {
                    "type": "scale",
                    "scale": [
                        physical_scale * (2.0 ** i),
                        physical_scale * (2.0 ** i),
                    ],
                }
            ],
        }
        for i in range(level)
    ]

    label_data_grp.attrs["multiscales"] = [
        {
            "version": "0.4",
            "name": label_name,
            "axes": [
                {"name": "y", "type": "space", "unit": "micrometer"},
                {"name": "x", "type": "space", "unit": "micrometer"},
            ],
            "coordinateTransformations": [
                {"type": "scale", "scale": [1.0, 1.0]}
            ],
            "datasets": datasets_meta,
        }
    ]

    logger.info(f"  Mask 寫入完成：labels/{label_name}")


# ---------------------------------------------------------------------------
# Transcript points
# ---------------------------------------------------------------------------

def _write_points(root, config: dict, store_path: str, physical_scale: float = 1.0, roi_crop=None) -> None:
    """
    從 Visium HD binned matrix 爆炸轉錄點位（含隨機 jitter），
    以 Dask Parquet 格式寫入 points/transcripts。
    roi_crop: {"x0", "y0", "x1", "y1"} fullres px，若設定則篩選 ROI 範圍。
    """
    import pandas as pd
    import scanpy as sc
    import dask.dataframe as dd

    logger.info("處理 Points（Gene Transcripts）")

    pos_path = str(config["positions"])
    matrix_path = str(config["matrix"])

    if not os.path.exists(pos_path):
        logger.error(f"  Positions 不存在：{pos_path}")
        return
    if not os.path.exists(matrix_path):
        logger.error(f"  Matrix 不存在：{matrix_path}")
        return

    df_pos = pd.read_parquet(pos_path)
    if "pxl_col_in_fullres" not in df_pos.columns or "pxl_row_in_fullres" not in df_pos.columns:
        logger.error("  Positions 缺少座標欄位")
        return
    if "barcode" in df_pos.columns:
        df_pos = df_pos.set_index("barcode")

    # ROI 篩選
    if roi_crop is not None:
        x0, y0 = roi_crop["x0"], roi_crop["y0"]
        x1, y1 = roi_crop["x1"], roi_crop["y1"]
        in_roi = (
            (df_pos["pxl_col_in_fullres"] >= x0) & (df_pos["pxl_col_in_fullres"] < x1) &
            (df_pos["pxl_row_in_fullres"] >= y0) & (df_pos["pxl_row_in_fullres"] < y1)
        )
        df_pos = df_pos[in_roi]
        logger.info(f"  ROI 篩選後：{len(df_pos)} spots")

    logger.info(f"  載入 H5 Matrix：{matrix_path}")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        adata = sc.read_10x_h5(matrix_path)
    adata.var_names_make_unique()
    # Sanitize bytes→str (same as _write_table)
    def _dec(v):
        if isinstance(v, bytes):
            try:
                return v.decode("utf-8")
            except UnicodeDecodeError:
                return v.decode("latin-1")
        return v
    adata.var_names = pd.Index([_dec(v) for v in adata.var_names])
    adata.obs_names = pd.Index([_dec(v) for v in adata.obs_names])

    common = adata.obs_names.intersection(df_pos.index)
    if len(common) == 0:
        logger.error("  barcode 無交集！")
        return
    logger.info(f"  對齊 {len(common):,} bins")
    adata = adata[common].copy()
    df_pos = df_pos.loc[common]

    sf = config.get("scale_factor", 2.0)
    logger.info(f"  scale factor：{sf}")

    coords_unique = df_pos["pxl_col_in_fullres"].sort_values().unique()
    if len(coords_unique) > 1:
        stride = np.median(np.diff(coords_unique))
    else:
        stride = 1.0
        logger.warning("  無法計算 stride，使用預設 1.0")

    bin_size = stride * sf
    half_bin = bin_size / 2.0
    logger.info(f"  bin_size={bin_size:.4f}, jitter=±{half_bin:.4f}")

    X_coo = adata.X.tocoo()
    counts = X_coo.data.astype(int)
    total_points = int(counts.sum())
    logger.info(f"  總轉錄點位數：{total_points:,}")

    raw_coords = df_pos[["pxl_col_in_fullres", "pxl_row_in_fullres"]].values.astype(float)
    if roi_crop is not None:
        raw_coords[:, 0] -= roi_crop["x0"]
        raw_coords[:, 1] -= roi_crop["y0"]
    coords_base = raw_coords * sf
    row_repeats = np.repeat(X_coo.row, counts)
    col_repeats = np.repeat(X_coo.col, counts)

    x_base = coords_base[row_repeats, 0]
    y_base = coords_base[row_repeats, 1]

    rng = np.random.default_rng(42)
    jitter_x = rng.uniform(-half_bin, half_bin, size=len(row_repeats))
    jitter_y = rng.uniform(-half_bin, half_bin, size=len(row_repeats))

    x_final = (x_base + jitter_x) * physical_scale
    y_final = (y_base + jitter_y) * physical_scale
    gene_names = adata.var_names[col_repeats].values

    del X_coo, row_repeats, col_repeats, x_base, y_base, jitter_x, jitter_y
    gc.collect()

    logger.info("  建構 Points DataFrame...")
    points_df = pd.DataFrame({"x": x_final, "y": y_final, "gene": gene_names})
    del x_final, y_final, gene_names
    gc.collect()

    element_name = "transcripts"
    n_partitions = max(1, total_points // 5_000_000)
    logger.info(f"  使用 {n_partitions} 分區寫出 Parquet...")
    ddf = dd.from_pandas(points_df, npartitions=n_partitions)

    points_store_path = os.path.join(store_path, "points", element_name)
    points_data_path = os.path.join(points_store_path, "points.parquet")

    if os.path.exists(points_data_path):
        shutil.rmtree(points_data_path)

    ddf.to_parquet(points_data_path)

    # macOS 清理
    _clean_mac_junk(Path(points_data_path))

    # .zattrs
    zattrs_path = os.path.join(points_store_path, ".zattrs")
    points_attrs = {
        "axes": ["x", "y"],
        "coordinateTransformations": [
            {
                "type": "identity",
                "input": {
                    "name": "xy",
                    "axes": [
                        {"name": "x", "type": "space", "unit": "micrometer"},
                        {"name": "y", "type": "space", "unit": "micrometer"},
                    ],
                },
                "output": {
                    "name": "global",
                    "axes": [
                        {"name": "x", "type": "space", "unit": "micrometer"},
                        {"name": "y", "type": "space", "unit": "micrometer"},
                    ],
                },
            }
        ],
        "encoding-type": "ngff:points",
        "spatialdata_attrs": {"version": "0.1"},
    }
    with open(zattrs_path, "w") as f:
        json.dump(points_attrs, f)

    zgroup_path = os.path.join(points_store_path, ".zgroup")
    with open(zgroup_path, "w") as f:
        json.dump({"zarr_format": 2}, f)

    logger.info(f"  Points 寫入完成：points/{element_name}")


# ---------------------------------------------------------------------------
# 8µm bin shapes
# ---------------------------------------------------------------------------

def _write_shapes(root, config: dict, store_path: str, physical_scale: float = 1.0, roi_crop=None) -> None:
    """
    將 8µm 正方形 bins 以 GeoParquet Shapes 格式寫入 shapes/grid_008um。
    roi_crop: {"x0", "y0", "x1", "y1"} fullres px，若設定則篩選 ROI 範圍。
    """
    import pandas as pd
    import geopandas as gpd
    from shapely.geometry import box

    logger.info("處理 Shapes：8µm Grid")

    pos_path = str(config["positions"])
    if not os.path.exists(pos_path):
        logger.error(f"  Positions 不存在：{pos_path}")
        return

    df = pd.read_parquet(pos_path)

    # ROI 篩選
    if roi_crop is not None and "pxl_col_in_fullres" in df.columns:
        x0, y0 = roi_crop["x0"], roi_crop["y0"]
        x1, y1 = roi_crop["x1"], roi_crop["y1"]
        in_roi = (
            (df["pxl_col_in_fullres"] >= x0) & (df["pxl_col_in_fullres"] < x1) &
            (df["pxl_row_in_fullres"] >= y0) & (df["pxl_row_in_fullres"] < y1)
        )
        df = df[in_roi]
        logger.info(f"  ROI 篩選後：{len(df)} shapes")

    sf = config.get("scale_factor", 2.0)
    logger.info(f"  scale factor：{sf}")

    coords_unique = df["pxl_col_in_fullres"].sort_values().unique()
    if len(coords_unique) > 1:
        stride = np.median(np.diff(coords_unique))
    else:
        stride = 14.64  # 近似值 fallback
    side_length = stride * sf
    half_side = side_length / 2.0
    logger.info(f"  stride={stride:.4f}, side={side_length:.4f}")

    # ROI 座標偏移
    col_vals = df["pxl_col_in_fullres"].values.astype(float)
    row_vals = df["pxl_row_in_fullres"].values.astype(float)
    if roi_crop is not None:
        col_vals -= roi_crop["x0"]
        row_vals -= roi_crop["y0"]

    cx = (col_vals * sf) * physical_scale
    cy = (row_vals * sf) * physical_scale

    geometries = [
        box(
            x - (half_side * physical_scale),
            y - (half_side * physical_scale),
            x + (half_side * physical_scale),
            y + (half_side * physical_scale),
        )
        for x, y in zip(cx, cy)
    ]

    cols_to_keep = [c for c in ["barcode", "array_row", "array_col"] if c in df.columns]
    gdf = gpd.GeoDataFrame(df[cols_to_keep], geometry=geometries)

    element_name = "grid_008um"
    shapes_store_path = os.path.join(store_path, "shapes", element_name)
    os.makedirs(shapes_store_path, exist_ok=True)

    parquet_path = os.path.join(shapes_store_path, "shapes.parquet")
    logger.info(f"  寫出 Shapes Parquet：{parquet_path}")
    gdf.to_parquet(parquet_path)

    zattrs_path = os.path.join(shapes_store_path, ".zattrs")
    shapes_attrs = {
        "axes": ["x", "y"],
        "coordinateTransformations": [
            {
                "type": "identity",
                "input": {
                    "name": "xy",
                    "axes": [
                        {"name": "x", "type": "space", "unit": "unit"},
                        {"name": "y", "type": "space", "unit": "unit"},
                    ],
                },
                "output": {
                    "name": "global",
                    "axes": [
                        {"name": "x", "type": "space", "unit": "unit"},
                        {"name": "y", "type": "space", "unit": "unit"},
                    ],
                },
            }
        ],
        "encoding-type": "ngff:shapes",
        "spatialdata_attrs": {"version": "0.2"},
    }
    with open(zattrs_path, "w") as f:
        json.dump(shapes_attrs, f)

    zgroup_path = os.path.join(shapes_store_path, ".zgroup")
    with open(zgroup_path, "w") as f:
        json.dump({"zarr_format": 2}, f)

    _clean_mac_junk(Path(shapes_store_path))
    logger.info(f"  Shapes 寫入完成：shapes/{element_name}")


# ---------------------------------------------------------------------------
# H&E image resolution helper
# ---------------------------------------------------------------------------

_LARGE_BTF_THRESHOLD = 10000  # pixels — above this, prefer hires PNG


def _resolve_he_image(
    he_path: Path,
    binned_002_dir: Path,
    pixel_size_um: float,
) -> tuple[Path, float]:
    """
    若 H&E 影像為超大 BTF（> LARGE_BTF_THRESHOLD px），
    改用 binned_002/spatial/tissue_hires_image.png，並根據
    scalefactors_json.json 自動換算 pixel size。
    """
    import tifffile, json

    he_str = str(he_path)
    if not he_str.lower().endswith((".btf", ".tiff", ".tif")):
        return he_path, pixel_size_um

    try:
        with tifffile.TiffFile(he_str) as tif:
            p0 = tif.pages[0]
            h, w = p0.imagelength, p0.imagewidth
    except Exception:
        return he_path, pixel_size_um

    if max(h, w) <= _LARGE_BTF_THRESHOLD:
        return he_path, pixel_size_um

    logger.info(
        f"BTF 解析度過大（{w}×{h}），改用 tissue_hires_image.png"
    )
    hires_png = binned_002_dir / "spatial" / "tissue_hires_image.png"
    sf_json   = binned_002_dir / "spatial" / "scalefactors_json.json"

    if not hires_png.exists():
        logger.warning("tissue_hires_image.png 不存在，仍使用原始 BTF（可能很慢）")
        return he_path, pixel_size_um

    hires_scalef = 1.0
    if sf_json.exists():
        try:
            with open(sf_json) as f:
                sf = json.load(f)
            hires_scalef = float(sf.get("tissue_hires_scalef", 1.0))
        except Exception as e:
            logger.warning(f"無法讀取 scalefactors_json：{e}")

    adjusted_pixel_size = pixel_size_um / hires_scalef
    logger.info(
        f"使用 hires PNG：{hires_png}  "
        f"scalef={hires_scalef:.6f}  "
        f"pixel_size={adjusted_pixel_size:.4f} µm/px"
    )
    return hires_png, adjusted_pixel_size


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def build_zarr(config: dict[str, Any]) -> Path:
    """
    建構 SpatialData OME-Zarr。

    整合流程：
    1. 建立 Zarr 根結構
    2. 寫出 H&E 影像金字塔
    3. 寫出 OME-NGFF metadata
    4. 寫出 2µm binned matrix 表格
    5. 寫出轉錄點位（Points）
    6. 寫出細胞核分割遮罩（Labels）
    7. 寫出 8µm bin Shapes
    8. 清理 macOS 垃圾

    Parameters
    ----------
    config : dict
        pipeline.yaml 配置字典

    Returns
    -------
    Path
        輸出 Zarr 路徑
    """
    paths = config["paths"]
    zarr_cfg = config.get("zarr_builder", {})
    files_cfg = zarr_cfg.get("files", {})
    params = zarr_cfg.get("parameters", {})
    seg_cfg = config.get("segmentation", {})

    he_image_path = resolve_path(paths["he_image"])
    binned_002_dir = resolve_path(paths["binned_002"])
    binned_008_dir = resolve_path(paths.get("binned_008", paths["binned_002"]))
    masks_dir = resolve_path(paths["masks_dir"])
    zarr_dir = resolve_path(paths["zarr_dir"])
    zarr_dir.mkdir(parents=True, exist_ok=True)

    pixel_size_um = params.get("pixel_size_um", VISIUM_UM_PX)
    out_filename = params.get("output_filename", "proseg_integrated.zarr")
    gene_scale_factor = params.get("gene_scale_factor", 1.0)
    mask_filename = seg_cfg.get("output", {}).get("mask_filename", "segmentation_masks.npy")

    h5_matrix      = files_cfg.get("h5_matrix",       "filtered_feature_bc_matrix.h5")
    spatial_dir    = files_cfg.get("spatial_dir",      "spatial")
    tissue_positions = files_cfg.get("tissue_positions", "tissue_positions.parquet")

    datasets = {
        "square_002um": {
            "matrix":       binned_002_dir / h5_matrix,
            "positions":    binned_002_dir / spatial_dir / tissue_positions,
            "scale_factor": gene_scale_factor,
        },
        "square_008um": {
            "matrix":       binned_008_dir / h5_matrix,
            "positions":    binned_008_dir / spatial_dir / tissue_positions,
            "scale_factor": gene_scale_factor,
        },
    }

    # ── 偵測 ROI 模式 ────────────────────────────────────────────────────────
    output_dir = resolve_path(paths.get("output_dir", "results/analysis"))
    roi_base   = output_dir / "roi"
    roi_list   = config.get("rois", [])

    roi_he_crops: list[tuple[str, Path, dict]] = []  # (roi_name, he_crop_path, roi_def)
    for roi in roi_list:
        roi_name = roi.get("name", "")
        he_crop  = roi_base / roi_name / "he_crop.tif"
        if he_crop.exists():
            roi_he_crops.append((roi_name, he_crop, roi))

    # fallback：掃描目錄
    if not roi_he_crops and roi_base.exists():
        for d in sorted(roi_base.iterdir()):
            if d.is_dir() and not d.name.startswith("."):
                he_crop = d / "he_crop.tif"
                if he_crop.exists():
                    roi_he_crops.append((d.name, he_crop, {}))

    if roi_he_crops:
        # ── Per-ROI 模式：每個 ROI 建立獨立 Zarr ────────────────────────────
        logger.info("=" * 60)
        logger.info(f"Stage 2: Zarr 建構（per-ROI 模式，{len(roi_he_crops)} 個 ROI）")
        logger.info("=" * 60)

        last_zarr = None
        for roi_name, he_crop_path, roi_def in roi_he_crops:
            roi_zarr_dir = zarr_dir / roi_name
            roi_zarr_dir.mkdir(parents=True, exist_ok=True)
            out_zarr = roi_zarr_dir / out_filename

            # ROI crop 範圍（fullres px）
            roi_crop = None
            if roi_def.get("x") is not None:
                rx  = int(roi_def["x"])
                ry  = int(roi_def["y"])
                rw  = int(roi_def["width_px"])
                rh  = int(roi_def["height_px"])
                roi_crop = {"x0": rx, "y0": ry, "x1": rx + rw, "y1": ry + rh}

            roi_mask_path = roi_base / roi_name / mask_filename

            logger.info(f"\n  ROI: {roi_name}")
            logger.info(f"  H&E crop:   {he_crop_path}")
            logger.info(f"  Mask:       {roi_mask_path}")
            logger.info(f"  ROI crop:   {roi_crop}")
            logger.info(f"  輸出 Zarr:  {out_zarr}")

            root = _create_zarr_structure(out_zarr)
            _clean_mac_junk(out_zarr)

            num_levels, pixel_size_obj = _write_image_pyramid(
                root, he_crop_path, pixel_size_override=pixel_size_um
            )
            physical_scale = pixel_size_obj["scale"]
            _write_ome_metadata(root, num_levels, pixel_size_obj)

            store_path_str = str(out_zarr)
            _write_table(root, "table", datasets["square_002um"], store_path_str, roi_crop=roi_crop)
            _write_points(root, datasets["square_002um"], store_path_str,
                          physical_scale=physical_scale, roi_crop=roi_crop)
            _add_masks_to_zarr(root, roi_mask_path, physical_scale=physical_scale,
                               label_name="cellpose_nuclei")
            _write_shapes(root, datasets["square_008um"], store_path_str,
                          physical_scale=physical_scale, roi_crop=roi_crop)
            _clean_mac_junk(out_zarr)

            logger.info(f"  ✅ ROI {roi_name} Zarr 完成：{out_zarr}")
            last_zarr = out_zarr

        logger.info("=" * 60)
        logger.info("Zarr 建構完成（per-ROI）")
        logger.info("=" * 60)
        return last_zarr

    # ── 全圖模式（fallback）───────────────────────────────────────────────────
    logger.info("=" * 60)
    logger.info("Stage 2: Zarr 建構（全圖模式）")
    logger.warning("未找到 he_crop.tif，改用全圖 BTF（可能耗時很長）")
    logger.info("=" * 60)

    out_zarr = zarr_dir / out_filename
    mask_path = masks_dir / mask_filename

    root = _create_zarr_structure(out_zarr)
    _clean_mac_junk(out_zarr)

    # 全圖模式：BTF 太大時使用 hires PNG
    effective_he, effective_px = _resolve_he_image(he_image_path, binned_002_dir, pixel_size_um)

    num_levels, pixel_size_obj = _write_image_pyramid(
        root, effective_he, pixel_size_override=effective_px
    )
    physical_scale = pixel_size_obj["scale"]
    _write_ome_metadata(root, num_levels, pixel_size_obj)

    store_path_str = str(out_zarr)
    _write_table(root, "table", datasets["square_002um"], store_path_str)
    _write_points(root, datasets["square_002um"], store_path_str, physical_scale=physical_scale)
    _add_masks_to_zarr(root, mask_path, physical_scale=physical_scale, label_name="cellpose_nuclei")
    _write_shapes(root, datasets["square_008um"], store_path_str, physical_scale=physical_scale)
    _clean_mac_junk(out_zarr)

    logger.info("=" * 60)
    logger.info("Zarr 建構完成")
    logger.info(f"  輸出：{out_zarr}")
    logger.info("=" * 60)
    return out_zarr
