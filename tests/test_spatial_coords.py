import pytest
import numpy as np
from backend.src.registration.spatial_coords import SpatialCoordinateTransformer

def test_spatial_coords_scaling():
    transformer = SpatialCoordinateTransformer(
        vfr_width=3000,
        vfr_height=3000,
        raw_width=30000,
        raw_height=30000,
        hires_scalef=0.1,
    )
    
    # 1. Point scaling
    raw_pt = np.array([15000.0, 15000.0])
    vfr_pt = transformer.transform_points(raw_pt, from_space="raw_tiff", to_space="virtual_fullres")
    assert np.allclose(vfr_pt, [1500.0, 1500.0])
    
    # 2. Reverse scaling
    back_to_raw = transformer.transform_points(vfr_pt, from_space="virtual_fullres", to_space="raw_tiff")
    assert np.allclose(back_to_raw, raw_pt)
    
    # 3. Hires scaling
    hires_pt = transformer.transform_points(vfr_pt, from_space="virtual_fullres", to_space="tissue_hires")
    assert np.allclose(hires_pt, [150.0, 150.0])

def test_spatial_coords_rotation():
    transformer = SpatialCoordinateTransformer(
        vfr_width=1000,
        vfr_height=2000,
        raw_width=1000,
        raw_height=2000,
        rotation_deg=90
    )
    
    pt = np.array([[100.0, 200.0]])
    # 90 deg clockwise with H=2000: (x, y) -> (H - y, x) = (2000 - 200, 100) = (1800, 100)
    rotated = transformer.transform_points(pt, from_space="raw_tiff", to_space="raw_tiff", apply_rotation=True)
    assert np.allclose(rotated, [[1800.0, 100.0]])

def test_transform_polygons():
    transformer = SpatialCoordinateTransformer(
        vfr_width=1000,
        vfr_height=1000,
        raw_width=10000,
        raw_height=10000,
    )
    
    poly = [[[1000.0, 1000.0], [2000.0, 1000.0], [2000.0, 2000.0], [1000.0, 2000.0]]]
    t_poly = transformer.transform_polygons([poly], from_space="raw_tiff", to_space="virtual_fullres")
    expected = [[[[100.0, 100.0], [200.0, 100.0], [200.0, 200.0], [100.0, 200.0]]]]
    assert np.allclose(t_poly, expected)
