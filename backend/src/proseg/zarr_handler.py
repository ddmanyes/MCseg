"""
Zarr Handler Module

Core module for loading and extracting data from SpatialData Zarr files.
Provides functions to read transcripts, masks, and images with proper
coordinate system alignment verification.
"""

import os
from typing import Optional, Tuple, Dict, List
from pathlib import Path
import warnings
import json
import logging

# Configure logger
logger = logging.getLogger(__name__)

import dask
try:
    dask.config.set({'dataframe.query-planning': False})
except Exception:
    pass # Older dask might not have this key, which is fine as it defaults to legacy
import dask.dataframe as dd

import numpy as np
import pandas as pd
import spatialdata as sd
from spatialdata import SpatialData

# Import centralized constants
try:
    from src.constants import SCALE_NM_PX, SCALE_UM_PX
except ImportError:
    # Fallback if not installed as regular package
    SCALE_NM_PX = 264.5833333333333
    SCALE_UM_PX = 0.2645833333333333
from spatial_image import SpatialImage
import dask.array as da


def load_zarr(zarr_path: str) -> SpatialData:
    """
    載入 SpatialData Zarr 檔案
    
    Parameters
    ----------
    zarr_path : str
        Zarr 檔案路徑
        
    Returns
    -------
    SpatialData
        載入的 SpatialData 物件
        
    Raises
    ------
    FileNotFoundError
        如果 Zarr 檔案不存在
    ValueError
        如果 Zarr 檔案格式不正確
    
    Examples
    --------
    >>> sdata = load_zarr("data/proseg_integrated.zarr")
    >>> print(sdata)
    """
    if not os.path.exists(zarr_path):
        raise FileNotFoundError(f"Zarr file not found: {zarr_path}")
    
    try:
        # Suppress Zarr warnings about hidden files
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=UserWarning, module="zarr")
            sdata = sd.read_zarr(zarr_path)
        
        logger.info(f"✅ 成功載入 Zarr: {zarr_path}")
        logger.debug(f"  - Images: {list(sdata.images.keys())}")
        logger.debug(f"  - Labels: {list(sdata.labels.keys())}")
        logger.debug(f"  - Points: {list(sdata.points.keys())}")
        logger.debug(f"  - Shapes: {list(sdata.shapes.keys())}")
        logger.debug(f"  - Tables: {list(sdata.tables.keys())}")
        
        return sdata
    
    except Exception as e:
        # 如果 spatialdata.read_zarr 失敗（例如 ValueError 或 metadata assertion failure）
        # 一旦失敗，我們就嘗試手動容錯載入，因為本專案的萃取函數支援 _load_zarr_fallback 的資料結構
        logger.warning(f"⚠️  標準載入失敗 ({type(e).__name__}: {e})，嘗試容錯模式...")
        try:
            return _load_zarr_fallback(zarr_path)
        except Exception as fallback_e:
            raise ValueError(f"無法讀取 Zarr 檔案: {fallback_e}")


def _load_zarr_fallback(zarr_path: str) -> SpatialData:
    """
    Fallback: 手動重建 SpatialData 物件
    
    當 spatialdata.read_zarr 因為 metadata 問題失敗時使用。
    """
    import zarr
    import xarray as xr
    import dask.array as da
    
    logger.info(f"  📦 使用容錯模式載入 Zarr: {zarr_path}")
    
    z = zarr.open(zarr_path, mode='r')
    elements = {}
    
    # 簡化的 SpatialData 模擬類別
    class SimplifiedSpatialData:
        def __init__(self, zarr_path, elements, zarr_store):
            self.zarr_path = zarr_path
            self._elements = elements
            self._zarr_store = zarr_store
            
        @property
        def images(self):
            return {k: v for k, v in self._elements.items() if k.startswith('tissue')}
        
        @property
        def labels(self):
            return {k: v for k, v in self._elements.items() if k.startswith('cellpose')}
        
        @property
        def points(self):
            # Points 需要特殊處理，直接從 Zarr 讀取
            if 'points' in self._zarr_store and 'transcripts' in self._zarr_store['points']:
                 return {'transcripts': self._zarr_store['points']['transcripts']}
            return {}
        
        @property
        def coordinate_systems(self):
             return ["global"] # Dummy

    # 1. 載入 Images (tissue_hires_image)
    if 'images' in z and 'tissue_hires_image' in z['images']:
        try:
            img_group = z['images']['tissue_hires_image']
            # 建立 DataArray (簡化版，取 scale 0)
            if '0' in img_group:
                scale0_path = os.path.join(zarr_path, 'images/tissue_hires_image/0')
                # Load as dask array
                darr = da.from_zarr(scale0_path)
                # Assuming (c, y, x)
                arr = xr.DataArray(darr, dims=['c', 'y', 'x'], name='tissue_hires_image')
                arr = xr.DataArray(darr, dims=['c', 'y', 'x'], name='tissue_hires_image')
                elements['tissue_hires_image'] = arr
                logger.debug(f"  - 載入 Images: tissue_hires_image (Scale 0)")
        except Exception as e:
            logger.warning(f"    ⚠️  Images 載入失敗: {e}")
    
    # 2. 載入 Labels (cellpose_nuclei, cellpose_cyto)
    if 'labels' in z:
        for label_name in ['cellpose_nuclei', 'cellpose_cyto', 'cellpose']:  # Try all possible names
            if label_name in z['labels']:
                try:
                    label_group = z['labels'][label_name]
                    if '0' in label_group:
                        scale0_path = os.path.join(zarr_path, f'labels/{label_name}/0')
                        darr = da.from_zarr(scale0_path)
                        # Assuming (y, x)
                        arr = xr.DataArray(darr, dims=['y', 'x'], name=label_name)
                        elements[label_name] = arr
                        logger.debug(f"  - 載入 Labels: {label_name} (Scale 0)")
                except Exception as e:
                    logger.warning(f"    ⚠️  Labels {label_name} 載入失敗: {e}")

    sdata = SimplifiedSpatialData(zarr_path, elements, z)
    logger.info(f"  ✅ 容錯載入完成")
    return sdata


def extract_transcripts(
    sdata: SpatialData,
    points_name: str = "transcripts",
    required_columns: Optional[list] = None
) -> pd.DataFrame:
    """
    提取轉錄點位資料
    
    Parameters
    ----------
    sdata : SpatialData
        SpatialData 物件
    points_name : str, default="transcripts"
        Points 元素名稱
    required_columns : list, optional
        必須包含的欄位，預設為 ['x', 'y', 'gene']
        
    Returns
    -------
    pd.DataFrame
        轉錄點位 DataFrame，包含 x, y, gene 等欄位
        
    Raises
    ------
    KeyError
        如果指定的 Points 元素不存在
    ValueError
        如果缺少必要欄位
    
    Examples
    --------
    >>> sdata = load_zarr("data/proseg_integrated.zarr")
    >>> df = extract_transcripts(sdata)
    >>> print(df.head())
    """
    if required_columns is None:
        required_columns = ['x', 'y', 'gene']
    
    # Check if this is our simplified fallback object (which might store store in _zarr_store)
    # Standard SpatialData object usually has .points
    try:
        points_dict = sdata.points
    except AttributeError:
        # Should not happen with our SimplifiedSpatialData since we defined .points
        points_dict = {}

    if points_name not in points_dict:
        available = list(points_dict.keys())
        # Check if fallback storage is available
        # Check if fallback storage is available
        if hasattr(sdata, '_zarr_store'):
             logger.warning(f"  ⚠️  使用容錯模式讀取 Points (from store)")
             return _extract_transcripts_from_zarr(sdata._zarr_store, points_name, required_columns)
             
        raise KeyError(
            f"Points 元素 '{points_name}' 不存在。"
            f"可用的元素: {available}"
        )
    
    # 取得 Points（可能是 Dask DataFrame 或 Zarr Group if fallback）
    points_data = points_dict[points_name]
    
    # Handle Fallback Zarr Group directly
    import zarr
    if isinstance(points_data, (zarr.Array, zarr.Group)):
         return _extract_transcripts_from_zarr(sdata._zarr_store, points_name, required_columns)

    # 如果是 Dask DataFrame，計算成 Pandas DataFrame
    # 如果是 Dask DataFrame，計算成 Pandas DataFrame
    if hasattr(points_data, 'compute'):
        logger.info(f"⏳ 正在載入 {points_name} 轉錄點位（Dask 延遲載入）...")
        df = points_data.compute()
    else:
        df = points_data
    
    # 驗證必要欄位
    missing_cols = set(required_columns) - set(df.columns)
    if missing_cols:
        raise ValueError(
            f"轉錄點位缺少必要欄位: {missing_cols}。"
            f"可用欄位: {list(df.columns)}"
        )
    
    logger.info(f"✅ 提取轉錄點位: {len(df):,} 個點")
    logger.debug(f"  - 欄位: {list(df.columns)}")
    if 'gene' in df.columns:
        logger.debug(f"  - 基因數: {df['gene'].nunique():,}")
    
    return df


def _extract_transcripts_from_zarr(zarr_store, points_name: str, required_columns: list) -> pd.DataFrame:
    """
    從 Zarr store 直接讀取轉錄點位（容錯模式）
    """
    import pyarrow.parquet as pq
    
    try:
        # Points 通常儲存為 Parquet 格式在 .zarr/points/NAME/points.parquet
        # 但 SpatialData Zarr 的 points 結構通常是: /points/NAME 且包含 parquet 檔案
        # 或者 /points/NAME 是一個 zarr group, 且 .zattrs 指向 parquet
        
        # Let's try to find the parquet file.
        # Construct path manually assuming standard structure
        base_path = zarr_store.store.path # This might be the root of the .zarr directory
        
        # Path might be .../points/transcripts/points.parquet or .../points/transcripts/parquet
        # Let's check common locations
        import glob
        
        # Assuming zarr_store is the root group
        # If it's a DirectoryStore
        if hasattr(zarr_store.store, 'path'):
             root_path = zarr_store.store.path
             search_path = os.path.join(root_path, 'points', points_name, '*.parquet')
             files = glob.glob(search_path)
             if not files:
                 # Check subdirectories
                 search_path = os.path.join(root_path, 'points', points_name, '**', '*.parquet')
                 files = glob.glob(search_path, recursive=True)
                 
             if files:
                 parquet_file = files[0] # Take the first one
                 logger.debug(f"    - 讀取 Parquet: {parquet_file}")
                 table = pq.read_table(parquet_file)
                 df = table.to_pandas()
                 return df
        
        # If we can't find parquet via file system (e.g. S3), or no files found
        raise NotImplementedError("無法找到 Points 對應的 Parquet 檔案")

    except Exception as e:
        # Last resort: Try reading as Zarr array if it's not parquet? 
        # But points are usually tables.
        raise ValueError(f"無法從 Zarr 讀取 Points: {e}")


def extract_masks(
    sdata: SpatialData,
    label_name: str = "cellpose",
    scale: int = 0,
    roi: Optional[Tuple[int, int, int, int]] = None
) -> np.ndarray:
    """
    提取細胞分割遮罩
    
    Parameters
    ----------
    sdata : SpatialData
        SpatialData 物件
    label_name : str, default="cellpose"
        Labels 元素名稱
    scale : int, default=0
        多尺度影像的尺度層級（0 為最高解析度）
    roi : tuple of int, optional
        感興趣區域 (min_y, max_y, min_x, max_x)。
        如果提供，將只載入並提取此區域的遮罩。
        
    Returns
    -------
    np.ndarray
        2D 遮罩陣列，值為細胞 ID（0 為背景）
        
    Raises
    ------
    KeyError
        如果指定的 Labels 元素不存在
    IndexError
        如果指定的 scale 不存在
    
    Examples
    --------
    >>> sdata = load_zarr("data/proseg_integrated.zarr")
    >>> mask = extract_masks(sdata, label_name='cellpose', scale=0)
    >>> print(f"Mask shape: {mask.shape}")
    >>> # With ROI: (min_y, max_y, min_x, max_x)
    >>> roi_mask = extract_masks(sdata, scale=0, roi=(1000, 2000, 1000, 2000))
    >>> print(f"ROI Mask shape: {roi_mask.shape}")
    """
    if label_name not in sdata.labels:
        available = list(sdata.labels.keys())
        raise KeyError(
            f"Labels 元素 '{label_name}' 不存在。"
            f"可用的元素: {available}"
        )
    
    labels_data = sdata.labels[label_name]
    
    # 處理多尺度 DataTree 結構
    try:
        import xarray as xr
        # 如果是 DataArray (單尺度)，直接使用
        if isinstance(labels_data, xr.DataArray):
            mask_array = labels_data.data # Use .data to keep it lazy (dask)
        # 如果是 DataTree 或 dict-like (多尺度)，取指定 scale
        elif hasattr(labels_data, '__getitem__') and hasattr(labels_data, 'keys'):
            if scale in labels_data:
                scale_data = labels_data[scale]
            elif str(scale) in labels_data:
                scale_data = labels_data[str(scale)]
            else:
                 keys = list(labels_data.keys())
                 if scale < len(keys):
                     scale_data = labels_data[keys[scale]]
                 else:
                     raise KeyError(f"Scale {scale} out of range for keys {keys}")
            
            if isinstance(scale_data, xr.DataArray):
                 mask_array = scale_data.data
            elif hasattr(scale_data, 'data_vars'):
                 # DataTree or Dataset: take first variable
                 var_names = list(scale_data.data_vars)
                 if var_names:
                     mask_array = scale_data[var_names[0]].data
                 else:
                     raise ValueError(f"No data variables found in scale {scale}")
            else:
                 mask_array = scale_data.data if hasattr(scale_data, 'data') else scale_data
        else:
             mask_array = labels_data.data if hasattr(labels_data, 'data') else labels_data
        
        # 如果是 Dask array，計算成 NumPy
        if isinstance(mask_array, da.Array):
            if roi:
                min_y, max_y, min_x, max_x = roi
                logger.info(f"✂️  應用 ROI 裁剪: y[{min_y}:{max_y}], x[{min_x}:{max_x}]")
                # 確保 ROI 在範圍內
                h, w = mask_array.shape[-2:] # Handle (C, Y, X) or (Y, X)
                min_y = max(0, min_y)
                min_x = max(0, min_x)
                max_y = min(h, max_y)
                max_x = min(w, max_x)
                
                if mask_array.ndim == 2:
                    mask_array = mask_array[min_y:max_y, min_x:max_x]
                else:
                    # (C, Y, X) or similar, slice spatial dims
                    mask_array = mask_array[..., min_y:max_y, min_x:max_x]
 
            logger.info(f"⏳ 正在載入 {label_name} 遮罩（Dask 延遲載入）...")
            mask_array = mask_array.compute()
        
        # 確保是 2D
        if mask_array.ndim > 2:
            # 如果有多個 channel，取第一個
            mask_array = mask_array[0] if mask_array.shape[0] < mask_array.shape[-1] else mask_array.squeeze()
        
        n_cells = np.unique(mask_array).size - 1  # 扣除背景 (0)
        logger.info(f"✅ 提取細胞遮罩: {mask_array.shape}")
        logger.debug(f"  - Scale: {scale}")
        logger.debug(f"  - 細胞數: {n_cells:,}")
        logger.debug(f"  - 背景像素比例: {(mask_array == 0).sum() / mask_array.size * 100:.1f}%")
        
        return mask_array
    
    except (IndexError, KeyError) as e:
        raise IndexError(f"無法取得 scale {scale} 的遮罩: {e}")


def extract_image(
    sdata: SpatialData,
    image_name: str = "tissue_hires_image",
    scale: int = 0,
    as_rgb: bool = True,
    roi: Optional[Tuple[int, int, int, int]] = None
) -> np.ndarray:
    """
    提取影像資料
    
    Parameters
    ----------
    sdata : SpatialData
        SpatialData 物件
    image_name : str, default="tissue_hires_image"
        Images 元素名稱
    scale : int, default=0
        多尺度影像的尺度層級（0 為最高解析度）
    as_rgb : bool, default=True
        是否轉換為 RGB 格式 (H, W, 3)
    roi : tuple, optional
        裁剪區域 (min_y, max_y, min_x, max_x)
        
    Returns
    -------
    np.ndarray
        影像陣列，形狀為 (H, W, 3) 或 (C, H, W)
        
    Raises
    ------
    KeyError
        如果指定的 Images 元素不存在
    
    Examples
    --------
    >>> sdata = load_zarr("data/proseg_integrated.zarr")
    >>> img = extract_image(sdata, scale=2)  # 較低解析度以節省記憶體
    >>> print(f"Image shape: {img.shape}")
    """
    if image_name not in sdata.images:
        available = list(sdata.images.keys())
        raise KeyError(
            f"Images 元素 '{image_name}' 不存在。"
            f"可用的元素: {available}"
        )
    
    image_data = sdata.images[image_name]
    logger.debug(f"DEBUG: image_data type: {type(image_data)}")
    if hasattr(image_data, 'keys'):
        logger.debug(f"DEBUG: keys: {list(image_data.keys())}")
    if hasattr(image_data, 'values'):
         logger.debug(f"DEBUG: image_data.values type: {type(image_data.values)}")
    
    # 處理多尺度 DataTree 結構
    try:
        import xarray as xr
        # 如果是 DataArray (單尺度)，直接使用
        if isinstance(image_data, xr.DataArray):
            img_array = image_data.data
        # 如果是 DataTree 或 dict-like (多尺度)，取指定 scale
        elif hasattr(image_data, '__getitem__') and hasattr(image_data, 'keys'):
            # 嘗試取得 scale
            if scale in image_data:
                scale_data = image_data[scale]
            elif str(scale) in image_data:
                scale_data = image_data[str(scale)]
            elif f"scale{scale}" in image_data:
                scale_data = image_data[f"scale{scale}"]
            else:
                 keys = list(image_data.keys())
                 if scale < len(keys):
                     scale_data = image_data[keys[scale]]
                 else:
                     raise KeyError(f"Scale {scale} out of range for keys {keys}")
                     
            if isinstance(scale_data, xr.DataArray):
                 img_array = scale_data.data
            elif "image" in scale_data:
                 img_array = scale_data["image"].data
            else:
                 keys = list(scale_data.keys())
                 if keys:
                     img_array = scale_data[keys[0]].data
                 else:
                     img_array = scale_data.data
        else:
             # Fallback
             img_array = getattr(image_data, 'data', image_data)

        # Lazy array (Dask) → NumPy
        if hasattr(img_array, 'compute'):
            if roi:
                min_y, max_y, min_x, max_x = roi
                logger.info(f"✂️  應用影像 ROI 裁剪: y[{min_y}:{max_y}], x[{min_x}:{max_x}]")
                h, w = img_array.shape[-2:]
                min_y = max(0, min_y)
                min_x = max(0, min_x)
                max_y = min(h, max_y)
                max_x = min(w, max_x)
                
                if img_array.ndim == 2:
                    img_array = img_array[min_y:max_y, min_x:max_x]
                else:
                    img_array = img_array[..., min_y:max_y, min_x:max_x]

            logger.info(f"⏳ 正在載入 {image_name} 影像（Scale {scale}）...")
            img_array = img_array.compute()
        
        # 轉換為 RGB 格式 (H, W, 3)
        if as_rgb and img_array.ndim == 3:
            if img_array.shape[0] == 3:  # (C, H, W) → (H, W, C)
                img_array = np.transpose(img_array, (1, 2, 0))
        
        logger.info(f"✅ 提取影像: {img_array.shape}")
        logger.debug(f"  - Scale: {scale}")
        logger.debug(f"  - Dtype: {img_array.dtype}")
        
        return img_array
    
    except (IndexError, KeyError) as e:
        raise IndexError(f"無法取得 scale {scale} 的影像: {e}")


def verify_alignment(
    sdata: SpatialData,
    points_name: str = "transcripts",
    label_name: str = "cellpose",
    sample_size: int = 100
) -> Tuple[bool, dict]:
    """
    驗證 Points 和 Labels 的座標系對齊
    
    Parameters
    ----------
    sdata : SpatialData
        SpatialData 物件
    points_name : str, default="transcripts"
        Points 元素名稱
    label_name : str, default="cellpose"
        Labels 元素名稱
    sample_size : int, default=100
        採樣點數用於驗證
        
    Returns
    -------
    is_aligned : bool
        是否對齊
    stats : dict
        統計資訊（重疊率、背景比例等）
        
    Examples
    --------
    >>> sdata = load_zarr("data/proseg_integrated.zarr")
    >>> is_aligned, stats = verify_alignment(sdata)
    >>> if is_aligned:
    ...     print(f"✅ 座標對齊驗證通過！重疊率: {stats['overlap_rate']:.1%}")
    """
    logger.info("🔍 驗證座標系對齊...")
    
    # 1. 檢查 coordinate systems
    try:
        coord_systems = sdata.coordinate_systems
        logger.debug(f"  - Coordinate systems: {coord_systems}")
    except AttributeError:
        logger.warning("  ⚠️  無法取得 coordinate systems 資訊")
    
    # 2. 提取資料
    transcripts_df = extract_transcripts(sdata, points_name=points_name)
    mask = extract_masks(sdata, label_name=label_name, scale=0)
    
    # 3. 隨機採樣點位
    if len(transcripts_df) > sample_size:
        samples = transcripts_df.sample(sample_size, random_state=42)
    else:
        samples = transcripts_df
    
    # 4. 檢查點位是否落在遮罩範圍內
    h_mask, w_mask = mask.shape
    in_bounds = 0
    in_cells = 0
    
    for _, row in samples.iterrows():
        x, y = int(row['x']), int(row['y'])
        
        # 檢查是否在影像範圍內
        if 0 <= x < w_mask and 0 <= y < h_mask:
            in_bounds += 1
            cell_id = mask[y, x]
            if cell_id > 0:  # 不是背景
                in_cells += 1
    
    # 5. 計算統計
    in_bounds_rate = in_bounds / len(samples)
    overlap_rate = in_cells / len(samples)
    
    stats = {
        'total_samples': len(samples),
        'in_bounds': in_bounds,
        'in_cells': in_cells,
        'in_bounds_rate': in_bounds_rate,
        'overlap_rate': overlap_rate,
        'mask_shape': mask.shape,
        'transcripts_count': len(transcripts_df),
    }
    
    logger.info(f"\n📊 對齊驗證結果:")
    logger.info(f"  - 採樣點數: {stats['total_samples']}")
    logger.info(f"  - 在影像範圍內: {stats['in_bounds']} ({in_bounds_rate:.1%})")
    logger.info(f"  - 落在細胞內: {stats['in_cells']} ({overlap_rate:.1%})")
    
    # 6. 判斷是否對齊
    is_aligned = in_bounds_rate > 0.95 and overlap_rate > 0.5
    
    
    if is_aligned:
        logger.info("  ✅ 座標系對齊驗證通過！")
    else:
        logger.warning("  ⚠️  座標系可能未對齊，請檢查！")
        if in_bounds_rate <= 0.95:
            logger.warning(f"     原因: 太多點位超出影像範圍 ({in_bounds_rate:.1%})")
        if overlap_rate <= 0.5:
            logger.warning(f"     原因: 太少點位落在細胞內 ({overlap_rate:.1%})")
    
    return is_aligned, stats

def get_scale_factors(zarr_path: str, element_path: str = "images/tissue_hires_image") -> Tuple[float, float]:
    """
    從 Zarr Metadata 讀取 Scale Factors (y, x)
    
    Parameters
    ----------
    zarr_path : str
        Zarr 根目錄路徑
    element_path : str
        元素相對路徑 (例如 images/tissue_hires_image)
        
    Returns
    -------
    (scale_y, scale_x)
    """
    try:
        attrs_path = os.path.join(zarr_path, element_path, ".zattrs")
        if not os.path.exists(attrs_path):
             # Try assuming it's not nested? (unlikely for NGFF)
             logger.warning(f"  ⚠️  Metadata 不存在: {attrs_path}，使用預設高精度 Scale {SCALE_NM_PX}")
             return (SCALE_NM_PX, SCALE_NM_PX)
             
        with open(attrs_path, 'r') as f:
            attrs = json.load(f)
            
        # Parse multiscales
        if "multiscales" in attrs and len(attrs["multiscales"]) > 0:
            multiscales = attrs["multiscales"][0]
            
            # Identify axes to find y and x specifically
            # Axes list: [{"name": "c", "type": "channel"}, {"name": "y", "type": "space"}, ...]
            axes = multiscales.get("axes", [])
            y_idx, x_idx = -2, -1 # Default fallback
            
            if axes:
                for idx, ax in enumerate(axes):
                     name = ax.get("name", "").lower()
                     if name == "y": y_idx = idx
                     if name == "x": x_idx = idx
            
            if "datasets" in multiscales and len(multiscales["datasets"]) > 0:
                # Get level 0 (base resolution)
                dataset0 = multiscales["datasets"][0]
                
                # Check coordinateTransformations
                if "coordinateTransformations" in dataset0:
                    for t in dataset0["coordinateTransformations"]:
                        if t["type"] == "scale":
                            scale = t["scale"]
                            # Use identified indices
                            try:
                                sy = float(scale[y_idx])
                                sx = float(scale[x_idx])
                                return (sy, sx)
                            except IndexError:
                                logger.warning(f"  ⚠️  Scale array ({scale}) dimension mismatch with axes indices y={y_idx}, x={x_idx}")
                                # Fallback to last two
                                if len(scale) >= 2:
                                    return (float(scale[-2]), float(scale[-1]))
                                
        logger.warning(f"  ⚠️  無法解析 Scale Metadata，使用預設 1.0")
        return (1.0, 1.0)
        
    except Exception as e:
        logger.warning(f"  ⚠️  讀取 Scale 失敗 ({e})，使用預設 1.0")
        return (1.0, 1.0)


def calculate_roi(
    transcripts_df: pd.DataFrame,
    scale_factors: Tuple[float, float],
    padding: int = 100
) -> Tuple[Tuple[int, int], Tuple[int, int, int, int]]:
    """
    計算感興趣區域 (ROI)，包含 Padding 與 Scale 轉換
    
    Parameters
    ----------
    transcripts_df : pd.DataFrame
        轉錄點位 DataFrame (欄位: x, y)
    scale_factors : tuple
        (scale_y, scale_x)
    padding : int, default=100
        外擴像素 (Pixels)
        
    Returns
    -------
    roi_offset : tuple (min_x, min_y)
        ROI 左上角偏移量 (Pixels)
    roi_box : tuple (min_y, max_y, min_x, max_x)
        ROI 邊界盒 (Pixels)，可用於 mask slicing
    """
    scale_y, scale_x = scale_factors
    
    # Global Extent (Physical units)
    min_x = int(transcripts_df['x'].min())
    max_x = int(transcripts_df['x'].max())
    min_y = int(transcripts_df['y'].min())
    max_y = int(transcripts_df['y'].max())
    
    logger.debug(f"  - Transcripts extent (Global): x[{min_x}:{max_x}], y[{min_y}:{max_y}]")
    
    # Physical -> Pixels
    pixel_min_x = int(min_x / scale_x)
    pixel_max_x = int(max_x / scale_x)
    pixel_min_y = int(min_y / scale_y)
    pixel_max_y = int(max_y / scale_y)
    
    # Apply Padding
    roi_min_x = max(0, pixel_min_x - padding)
    roi_min_y = max(0, pixel_min_y - padding)
    roi_max_x = pixel_max_x + padding
    roi_max_y = pixel_max_y + padding
    
    logger.info(f"  - ROI (Pixels): x[{roi_min_x}:{roi_max_x}], y[{roi_min_y}:{roi_max_y}] (Padding: {padding})")
    
    roi_offset = (roi_min_x, roi_min_y)
    roi_box = (roi_min_y, roi_max_y, roi_min_x, roi_max_x)
    
    return roi_offset, roi_box


def add_shapes_to_zarr(
    zarr_path: str,
    shapes_name: str,
    gdf: pd.DataFrame,
    overwrite: bool = True
) -> None:
    """
    將 GeoDataFrame 形狀資料寫入 Zarr 的 shapes 群組
    
    Parameters
    ----------
    zarr_path : str
        Zarr 檔案路徑
    shapes_name : str
        Shapes 元素名稱 (例如 "proseg_polygons")
    gdf : geopandas.GeoDataFrame
        包含 geometry 欄位的 GeoDataFrame
    overwrite : bool
        是否覆蓋已存在的 shapes
    """
    import zarr
    import shutil
    from spatialdata.models import ShapesModel
    from spatialdata import SpatialData
    
    logger.info(f"💾 正在將 {shapes_name} 寫入 Zarr {zarr_path}...")
    
    # 1. 轉換為 SpatialData ShapesModel
    # 確保有 geometry 欄位
    if 'geometry' not in gdf.columns:
        raise ValueError("GeoDataFrame 必須包含 'geometry' 欄位")
        
    # 設定變換矩陣 - 假設已經是 Global 座標
    transform = {"global": sd.transformations.Identity()}
    
    parsed_shapes = ShapesModel.parse(gdf, transformations=transform)
    
    # 2. 寫入 Zarr - 使用 Temp Zarr + Copy 策略避免鎖定與覆蓋問題
    import tempfile
    
    try:
        # Create a temporary directory for the isolated Zarr
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_zarr_path = os.path.join(temp_dir, "temp.zarr")
            
            # Create a temporary SpatialData object with ONLY the new shapes
            # This ensures clean metadata generation
            temp_sdata = SpatialData(shapes={shapes_name: parsed_shapes})
            
            # Write to temporary Zarr
            temp_sdata.write(temp_zarr_path)
            
            # Now copy the specific shape group to the real Zarr
            # Source: temp.zarr/shapes/shapes_name
            # Target: real.zarr/shapes/shapes_name
            
            src_store = zarr.DirectoryStore(temp_zarr_path)
            dst_store = zarr.DirectoryStore(zarr_path)
            
            src_group = zarr.open_group(store=src_store, mode='r')
            dst_group = zarr.open_group(store=dst_store, mode='r+') # Read-write
            
            # Ensure 'shapes' group exists in destination
            if 'shapes' not in dst_group:
                dst_group.create_group('shapes')
                
            dst_shapes_group = dst_group['shapes']
            src_shapes_group = src_group['shapes']
            
            if shapes_name in src_shapes_group:
                # Check if target exists
                if shapes_name in dst_shapes_group:
                    if overwrite:
                        logger.warning(f"  ⚠️  Shapes {shapes_name} 已存在，正在覆蓋...")
                        # Delete existing using shutil to handle macOS hidden files
                        import shutil
                        shapes_dir = Path(zarr_path) / "shapes" / shapes_name
                        if shapes_dir.exists():
                            shutil.rmtree(shapes_dir, ignore_errors=True)
                    else:
                        logger.warning(f"  ⚠️  Shapes {shapes_name} 已存在，跳過寫入。")
                        return

                # Copy
                logger.info(f"  📋 複製 Shapes 資料到目標 Zarr...")
                zarr.copy(
                    source=src_shapes_group[shapes_name],
                    dest=dst_shapes_group,
                    name=shapes_name,
                    if_exists='replace' if overwrite else 'raise'
                )
                
                # We also need to update the root .zattrs if necessary? 
                # SpatialData 0.0.12+ puts checking in individual elements usually.
                # But sometimes 'shapes' group has valid info.
                
                # Consolidate metadata if needed (optional)
                # zarr.consolidate_metadata(dst_store) 
                
                logger.info(f"✅ 成功寫入 Shapes: {shapes_name}")
            else:
                 logger.error(f"❌ 暫存 Zarr 中未找到 Shapes: {shapes_name}")

    except Exception as e:
        logger.error(f"❌ 寫入 Shapes 失敗: {e}")
        # Detailed debugging
        import traceback
        logger.debug(traceback.format_exc())
        raise
