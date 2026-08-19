"""
品質控制與形態學指標計算入口 (QC Metrics)
架構深化版：透過 backend.src.segmentation.results.SegmentationResult 提供統一計算。
"""

import json
import logging
import sys
from pathlib import Path

# 加入專案路徑
sys.path.insert(0, str(Path(__file__).parent.parent))

from backend.src.segmentation.results import SegmentationResult

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger("mcseg.qc_metrics")


def compute_metrics_from_geojson(geojson_path: str | Path, summary_json_path: str | Path | None = None) -> dict:
    """從 GeoJSON 檔案載入並計算形態學 QC 指標。"""
    geojson_path = Path(geojson_path)
    with open(geojson_path, "r", encoding="utf-8") as f:
        data = json.load(f)

    cells = []
    for feat in data.get("features", []):
        props = feat.get("properties", {})
        geom = feat.get("geometry", {})
        cells.append({
            "cell_id": props.get("cell_id", len(cells)),
            "polygon": geom.get("coordinates", []),
            "area": props.get("area_px", 0.0),
            "solidity": props.get("solidity", 0.0),
        })

    result = SegmentationResult(
        sample_id=geojson_path.stem.replace("_mcseg_if_fullslide", ""),
        cells=cells,
        image_shape=(0, 0),
    )

    qc = result.compute_qc_metrics()
    logger.info(f"Computed QC for {geojson_path.name}: {qc[total_cells]} cells, mean area: {qc[area_mean_px]:.1f} px")

    if summary_json_path:
        out_p = Path(summary_json_path)
        out_p.parent.mkdir(parents=True, exist_ok=True)
        with open(out_p, "w", encoding="utf-8") as f:
            json.dump(qc, f, indent=2)

    return qc


if __name__ == "__main__":
    if len(sys.argv) > 1:
        compute_metrics_from_geojson(sys.argv[1])
