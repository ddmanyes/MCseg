import pytest
import numpy as np
from PIL import Image
from backend.src.fullslide.streaming_segmenter import StreamingSlideSegmenter

def test_thumbnail_mask_generation(tmp_path):
    # 創建合成測試影像 (1000x1000)，中間有 200x200 亮區 (組織)，其餘全黑 (背景)
    img_arr = np.zeros((1000, 1000, 3), dtype=np.uint8)
    img_arr[400:600, 400:600] = 200  # 組織區
    
    test_img_path = tmp_path / "test_slide.png"
    Image.fromarray(img_arr).save(test_img_path)
    
    segmenter = StreamingSlideSegmenter(
        tile_size=200,
        stride=180,
        downsample_factor=10,
        device="cpu"
    )
    
    fg_mask, W, H = segmenter.generate_thumbnail_tissue_mask(test_img_path)
    assert W == 1000
    assert H == 1000
    assert fg_mask.shape == (100, 100)
    
    # 驗證背景區被標為 False，組織區被標為 True
    assert not segmenter._is_tissue_tile(fg_mask, x=0, y=0)
    assert segmenter._is_tissue_tile(fg_mask, x=400, y=400)
