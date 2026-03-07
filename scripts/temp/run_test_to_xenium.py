"""
一次性測試腳本：在 test ROI 上跑 1 個 Proseg 條件（golden params）→ QC → 匯出 Xenium

執行方式（從 visiumHD_pipeline_2/ 根目錄）：
    uv run python scripts/temp/run_test_to_xenium.py

流程：
    1. ConditionTester - 1 個條件（max_dist=40, compactness=0.06, dilation=20）
    2. Scanpy QC（過濾低品質細胞 / 基因 → 正規化）
    3. XeniumExporter → results/export/test/xenium_from_cond_test/
"""

import logging
import sys
from pathlib import Path

# Working directory = visiumHD_pipeline_2/
PROJECT_ROOT = Path(__file__).parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

import scanpy as sc

from backend.src.utils.config import load_config, resolve_path
from backend.src.proseg.condition_tester import ConditionTester
from backend.src.analysis.preprocessing import Preprocessor
from backend.src.export.xenium_exporter import XeniumExporter

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger("run_test_to_xenium")


def main() -> None:
    config = load_config()
    paths = config["paths"]
    rois = config.get("rois", [{}])
    roi = rois[0]
    roi_name = roi.get("name", "test")
    pixel_size_um = roi.get("pixel_size_um", 0.2737)

    output_dir = resolve_path(paths["output_dir"])
    roi_out_dir = output_dir / "roi" / roi_name
    conditions_dir = resolve_path(paths["conditions_dir"])
    zarr_dir = resolve_path(paths["zarr_dir"])
    export_dir = resolve_path(paths.get("export_dir", "results/export"))

    # ── Stage 2.5: 條件測試（單一 golden params 條件） ──────────────────
    logger.info("=== Stage 2.5: 條件測試（1 個 golden params 條件）===")
    tester = ConditionTester(config)
    grid = {
        "max_dist": [40.0],
        "compactness": [0.06],
        "dilation": [20],
    }
    results = tester.run_grid(grid, roi_name=roi_name)

    best = max(
        results,
        key=lambda r: r.get("n_cells", 0) * r.get("median_genes", 0),
    )
    best_idx = best["condition_idx"]
    logger.info(
        f"最佳條件：idx={best_idx}, n_cells={best['n_cells']}, "
        f"median_genes={best.get('median_genes', 0):.1f}, label={best.get('label', '')}"
    )

    cond_dir = conditions_dir / f"cond_{best_idx:02d}"
    cond_transcripts = cond_dir / "transcripts_for_proseg.csv"

    # Stage 4/5 使用 Stage 3 的完整 proseg 輸出（cell IDs 與 polygon full_id 對齊）
    # Stage 2.5 只負責確認 golden params，QC/匯出必須用完整跑的 proseg_cells.h5ad
    main_h5ad = roi_out_dir / "proseg_cells.h5ad"
    if not main_h5ad.exists():
        raise FileNotFoundError(
            f"找不到 Stage 3 proseg 輸出（請先完成 Stage 3）：{main_h5ad}"
        )

    # ── Stage 4: QC ────────────────────────────────────────────────────
    logger.info("=== Stage 4: Scanpy QC ===")
    adata = sc.read_h5ad(str(main_h5ad))
    logger.info(f"  載入：{adata.n_obs:,} 細胞  {adata.n_vars:,} 基因")

    analysis_cfg = config.get("analysis", {})
    pre_cfg = analysis_cfg.get("preprocessing", {})
    qc_params = pre_cfg.get("cellular", {})

    preprocessor = Preprocessor(analysis_cfg)
    adata = preprocessor.calculate_qc_metrics(adata, qc_params)
    adata = preprocessor.filter_cells(adata, qc_params)
    adata = preprocessor.filter_genes(adata, qc_params)
    adata = preprocessor.normalize(adata)

    logger.info(f"  QC 後：{adata.n_obs:,} 細胞  {adata.n_vars:,} 基因")

    qc_h5ad = roi_out_dir / "proseg_cells_qc.h5ad"
    roi_out_dir.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(str(qc_h5ad))
    logger.info(f"  QC h5ad 已儲存：{qc_h5ad}")

    # ── Stage 5: Xenium 匯出 ───────────────────────────────────────────
    logger.info("=== Stage 5: Xenium Explorer 匯出 ===")

    # 多邊形來源：優先使用 ROI 已有的 combined_proseg_results_qc.json
    poly_path = roi_out_dir / "combined_proseg_results_qc.json"
    if not poly_path.exists():
        if cond_poly.exists():
            poly_path = cond_poly
            logger.info(f"  使用條件結果 proseg_results.json：{poly_path}")
        else:
            poly_path = None
            logger.warning("  找不到任何多邊形 JSON，匯出將不含 shapes")

    # 轉錄點 CSV
    transcript_path = cond_transcripts if cond_transcripts.exists() else None
    if transcript_path is None:
        logger.warning("  找不到 transcripts_for_proseg.csv，匯出將不含 points")

    # Zarr 影像（可選）
    zarr_full = zarr_dir / roi_name / "proseg_integrated.zarr"
    zarr_arg = str(zarr_full) if zarr_full.exists() else None

    xenium_out_dir = export_dir / roi_name / "xenium_from_cond_test"
    exporter = XeniumExporter(
        zarr_path=zarr_arg,
        poly_json_path=str(poly_path) if poly_path else None,
        transcripts_csv_path=str(transcript_path) if transcript_path else None,
        pixel_size_um=pixel_size_um,
    )
    out = exporter.export(h5ad_path=str(qc_h5ad), output_dir=str(xenium_out_dir))
    logger.info(f"Xenium Explorer bundle 已輸出至：{out}")


if __name__ == "__main__":
    main()
