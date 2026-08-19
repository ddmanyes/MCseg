import pytest
import json
import numpy as np
from pathlib import Path
from backend.src.segmentation.results import SegmentationResult

def test_segmentation_result_qc_and_geojson(tmp_path):
    cells = [
        {
            "cell_id": 1,
            "polygon": [[[10.0, 10.0], [20.0, 10.0], [20.0, 20.0], [10.0, 20.0]]],
            "centroid": [15.0, 15.0],
            "area": 100.0,
            "solidity": 0.95
        },
        {
            "cell_id": 2,
            "polygon": [[[50.0, 50.0], [70.0, 50.0], [70.0, 70.0], [50.0, 70.0]]],
            "centroid": [60.0, 60.0],
            "area": 400.0,
            "solidity": 0.90
        }
    ]
    
    result = SegmentationResult(
        sample_id="test_sample",
        cells=cells,
        image_shape=(1000, 1000),
        metadata={"raw_instances": 3, "suppressed_duplicates": 1},
        runtime_seconds=12.5
    )
    
    assert len(result) == 2
    assert result.total_cells == 2
    
    # Test QC
    qc = result.compute_qc_metrics()
    assert qc["total_cells"] == 2
    assert qc["area_mean_px"] == 250.0
    assert qc["area_median_px"] == 250.0
    assert qc["solidity_mean"] == pytest.approx(0.925)
    assert qc["dedup_suppression_rate"] == pytest.approx(1/3)
    
    # Test GeoJSON export
    geojson_path = tmp_path / "test.geojson"
    result.to_geojson(geojson_path)
    assert geojson_path.exists()
    
    with open(geojson_path) as f:
        data = json.load(f)
    assert data["type"] == "FeatureCollection"
    assert len(data["features"]) == 2
    assert data["features"][0]["properties"]["area_px"] == 100.0
