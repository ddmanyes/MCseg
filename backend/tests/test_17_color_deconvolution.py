"""Regression: recover known H concentrations from row-oriented stain vectors."""

import cv2
import numpy as np
import pytest

from backend.src.segmentation.cellpose_runner import color_deconvolution_he


@pytest.mark.parametrize("quantized", [False, True])
def test_recovers_hematoxylin_from_stain_mixtures(monkeypatch, quantized):
    # Forward Beer-Lambert model; H varies independently of E and residual.
    matrix = np.array([
        [0.6500286, 0.7041680, 0.2860126],
        [0.0728940, 0.9904310, 0.1155140],
        [0.2688350, 0.5706770, 0.7768750],
    ])
    matrix /= np.linalg.norm(matrix, axis=1, keepdims=True)
    h, e, residual = np.meshgrid(
        [0.0, 0.25, 0.5, 1.0], [0.0, 0.5, 1.0], [0.0, 0.5],
        indexing="ij",
    )
    concentrations = np.stack([h, e, residual], axis=-1).reshape(4, 6, 3)
    rgb = 256.0 * np.exp(-(concentrations @ matrix)) - 1.0
    if quantized:
        rgb = np.rint(rgb).astype(np.uint8)

    # Inspect the channel before local contrast can obscure cross-talk.
    class IdentityClahe:
        def apply(self, image):
            return image

    monkeypatch.setattr(cv2, "createCLAHE", lambda **kwargs: IdentityClahe())
    actual = color_deconvolution_he(rgb)
    expected = concentrations[..., 0] * 255
    np.testing.assert_allclose(actual, expected, atol=5 if quantized else 1)


def test_real_clahe_output_contract_and_rgba():
    rgb = np.random.default_rng(17).integers(0, 256, (32, 48, 3), dtype=np.uint8)
    actual = color_deconvolution_he(rgb)
    assert actual.shape == rgb.shape[:2]
    assert actual.dtype == np.uint8
    rgba = np.concatenate([rgb, np.zeros((32, 48, 1), dtype=np.uint8)], axis=-1)
    np.testing.assert_array_equal(color_deconvolution_he(rgba), actual)
