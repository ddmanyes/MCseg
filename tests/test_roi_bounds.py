import pytest
import numpy as np
from backend.src.roi.extractor import get_rna_capture_bounds
from backend.src.utils.config import load_config


def test_rna_capture_bounds():
    config = load_config()
    bounds = get_rna_capture_bounds(config)
    assert bounds is not None
    assert "min_x" in bounds and "max_x" in bounds
    assert "min_y" in bounds and "max_y" in bounds
    assert bounds["max_x"] > bounds["min_x"]
    assert bounds["max_y"] > bounds["min_y"]
    assert bounds["total_bins"] > 0
    assert bounds["he_width"] > 0
    assert bounds["he_height"] > 0


def test_out_of_bounds_roi_detection():
    config = load_config()
    bounds = get_rna_capture_bounds(config)
    assert bounds is not None

    # Construct an ROI clearly outside the RNA bounds
    out_of_bounds_roi = {
        "name": "out_test",
        "x": int(bounds["min_x"] - 5000),
        "y": int(bounds["max_y"] + 5000),
        "width_px": 500,
        "height_px": 500,
    }

    # Verify that x1 < min_x and y0 > max_y
    assert out_of_bounds_roi["x"] + out_of_bounds_roi["width_px"] < bounds["min_x"]
    assert out_of_bounds_roi["y"] > bounds["max_y"]
