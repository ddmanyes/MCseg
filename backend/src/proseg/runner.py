import logging
import gc
from pathlib import Path
import anndata as ad
import zarr

from backend.src.proseg.pipeline import ProsegPipeline

logger = logging.getLogger("pipeline.proseg.tiling")

def merge_tiles(results_path: Path, output_file: Path):
    logger.info(f"🔍 正在搜尋分塊結果於：{results_path}")
    tile_dirs = sorted([d for d in results_path.iterdir() if d.is_dir() and d.name.startswith("tile_")])
    logger.info(f"  - 找到 {len(tile_dirs)} 個分塊目錄")
    
    adatas = []
    for tile_dir in tile_dirs:
        h5ad_file = tile_dir / "proseg_integrated.h5ad"
        if not h5ad_file.exists():
            logger.warning(f"  ⚠️  分塊 {tile_dir.name} 缺少 H5AD，跳過")
            continue
            
        logger.info(f"  📖 載入 {tile_dir.name}...")
        try:
            adata = ad.read_h5ad(h5ad_file)
            adata.obs['tile_id'] = tile_dir.name
            adata.obs_names = [f"{tile_dir.name}_{name}" for name in adata.obs_names]
            adatas.append(adata)
        except Exception as e:
            logger.error(f"  ❌ 載入失敗 {tile_dir.name}: {e}")
            
    if not adatas:
        raise ValueError("❌ 未找到任何有效的分塊 AnnData 進行合併！")
        
    logger.info(f"🔗 合併 {len(adatas)} 個分塊...")
    merged = ad.concat(adatas, join='outer', fill_value=0)
    
    logger.info(f"💾 儲存合併結果至：{output_file}")
    merged.write_h5ad(output_file)
    logger.info(f"✅ 合併成功！總細胞數：{merged.n_obs:,}，總基因數：{merged.n_vars:,}")

def run_tiled_proseg(config: dict) -> None:
    golden = config.get("proseg", {}).get("golden_params", {})
    tiling = config.get("proseg", {}).get("tiling", {})
    paths = config.get("paths", {})
    scale_um_px = config.get("proseg", {}).get("constants", {}).get("scale_um_px", 0.2645833)
    
    zarr_base = Path(paths.get("zarr_dir", "results/zarr"))
    out_base = Path(paths.get("output_dir", "results/analysis")) / "roi"
    
    rois = config.get("rois", [])
    if not rois:
        raise ValueError("未定義 ROI")
        
    for roi in rois:
        roi_name = roi["name"]
        
        # Override the global scale with the ROI-specific scale if available.
        # This prevents coordinates shifting between um and px if they differ.
        roi_scale_um_px = roi.get("pixel_size_um", scale_um_px)
        
        logger.info(f"[{roi_name}] 開始進行分塊 Proseg 分析")
        
        zarr_path = zarr_base / roi_name / "proseg_integrated.zarr"
        output_dir = out_base / roi_name / "proseg_tiles"
        output_dir.mkdir(parents=True, exist_ok=True)
        final_h5ad = out_base / roi_name / "proseg_cells.h5ad"
        
        if not zarr_path.exists():
            logger.error(f"找不到 Zarr 儲存庫：{zarr_path}，請先執行 Stage 2")
            continue
            
        # 取得 Label 大小
        label_path = zarr_path / "labels" / "cellpose_nuclei" / "0"
        if not label_path.exists():
            logger.error(f"Zarr 中無細胞核遮罩：{label_path}")
            continue
            
        z = zarr.open(str(label_path), "r")
        h_full, w_full = z.shape[-2:]
        logger.info(f"  - ROI 影像大小: {w_full} x {h_full}")
        
        grid_nx = tiling.get("grid_nx", 4)
        grid_ny = tiling.get("grid_ny", 3)
        padding = tiling.get("padding", 200)
        
        tile_w = w_full // grid_nx
        tile_h = h_full // grid_ny
        
        tasks = []
        for iy in range(grid_ny):
            for ix in range(grid_nx):
                x_min = ix * tile_w
                y_min = iy * tile_h
                x_max = min((ix + 1) * tile_w, w_full)
                y_max = min((iy + 1) * tile_h, h_full)
                w = x_max - x_min
                h = y_max - y_min
                if w <= 0 or h <= 0: continue
                tasks.append({
                    "id": f"tile_y{iy}_x{ix}",
                    "roi": (x_min, y_min, w, h)
                })
                
        logger.info(f"  - 劃分為 {len(tasks)} 個分塊")
        
        for task in tasks:
            tile_id = task["id"]
            roi_px = task["roi"]
            tile_out = output_dir / tile_id
            
            if (tile_out / "counts.csv.gz").exists():
                logger.info(f"⏩ 略過 {tile_id} (已完成)")
                continue
                
            # 自動尋找 cyto_mask.npy 如果 Zarr 沒有重新建構
            cyto_npy = None
            cyto_npy_path = roi_out_dir / "cyto_mask.npy"
            if cyto_npy_path.exists():
                cyto_npy = str(cyto_npy_path)

            logger.info(f"\n📦 處理 {tile_id} | 範圍: {roi_px}")
            pipeline = ProsegPipeline(
                zarr_path=str(zarr_path),
                output_dir=str(tile_out),
                max_dist=golden.get("max_dist", 40.0),
                compactness=golden.get("compactness", 0.06),
                dilation_radius=golden.get("dilation", 20),
                samples=golden.get("samples", 500),
                burnin_samples=golden.get("burnin_samples", int(golden.get("samples", 500) * 0.3)),
                recorded_samples=golden.get("recorded_samples", 150),
                coordinate_scale=roi_scale_um_px,
                padding=padding,
                nucleus_label_name="cellpose_nuclei",
                use_cyto_mask_from_zarr=True,       # 優先從 Zarr，失敗則退回 cyto_mask_path
                cyto_mask_path=cyto_npy,            # [新增] 自動探測的外部遮罩
                cyto_label_name="eosin_cyto",
                use_watershed=golden.get("use_watershed", True),
                enforce_connectivity=golden.get("enforce_connectivity", True),
                fixed_roi=roi_px
            )
            pipeline.run_full_pipeline()
            # 每個 Tile 後主動清理記憶體
            gc.collect()
            
        # 合併 Tile
        merge_tiles(output_dir, final_h5ad)
        logger.info(f"[{roi_name}] Tiled Proseg 處理完成！")
