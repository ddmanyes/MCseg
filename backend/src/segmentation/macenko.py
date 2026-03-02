# ==========================================================
# Provenance: cell_segmentation/src/utils.py
# Original Project: cellpose-he-segmentation
# Extracted: 2026-02-20
# ==========================================================
"""
Macenko 色彩標準化與影像前處理工具

提供 H&E 染色組織切片的色彩解構功能：
- MacenkoNormalizer: Macenko 方法的染色分離與 Hematoxylin 提取
- apply_clahe: 局部對比度增強
"""

import numpy as np
import cv2


class MacenkoNormalizer:
    """
    Macenko method for stain normalization and separation.
    Separates Hematoxylin (H) and Eosin (E) stains.
    """
    def __init__(self):
        self.HERef = np.array([[0.5626, 0.2159],
                               [0.7201, 0.8012],
                               [0.4062, 0.5581]])
        self.maxCRef = np.array([1.9705, 1.0308])
        self.stain_matrix = None

    def _convert_rgb_to_od(self, I: np.ndarray, Io: int = 240, beta: float = 0.15):
        """Convert RGB to Optical Density (OD)."""
        OD = -np.log((I.reshape((-1, 3)).astype(np.float64) + 1) / (Io + 1))
        mask = np.any(OD < beta, axis=1)
        ODhat = OD[~mask]
        return OD, ODhat

    def fit(self, I: np.ndarray, Io: int = 240, beta: float = 0.15) -> bool:
        """
        Fit the stain vectors (stain_matrix) based on the input image.
        Returns True if successful, False otherwise.
        """
        if I.ndim == 2:
            return False

        if I.ndim == 3 and I.shape[-1] == 4:
            I = I[..., :3]

        if I.ndim == 3 and I.shape[-1] != 3:
            return False

        I_reshaped = I.reshape((-1, 3))
        OD, ODhat = self._convert_rgb_to_od(I_reshaped, Io, beta)

        if ODhat.shape[0] < 100:
            print(f"  -> Warning: Only {ODhat.shape[0]} valid pixels at beta={beta}. Retrying with beta=0.05...")
            OD, ODhat = self._convert_rgb_to_od(I_reshaped, Io, beta=0.05)

        if ODhat.shape[0] < 100:
            print("  -> Debug: Image is still too empty/white for calibration (Valid Pixels < 100).")
            return False

        try:
            _, V = np.linalg.eigh(np.cov(ODhat, rowvar=False))
            V = V[:, [2, 1]]
            Phi = np.arctan2(np.dot(ODhat, V[:, 1]), np.dot(ODhat, V[:, 0]))

            minPhi = np.percentile(Phi, 1)
            maxPhi = np.percentile(Phi, 99)

            vMin = np.dot(V, np.array([np.cos(minPhi), np.sin(minPhi)]))
            vMax = np.dot(V, np.array([np.cos(maxPhi), np.sin(maxPhi)]))

            if vMin[0] > vMax[0]:
                HE = np.array([vMin, vMax])
            else:
                HE = np.array([vMax, vMin])

            self.stain_matrix = HE / np.linalg.norm(HE, axis=1)[:, None]
            return True

        except Exception as e:
            print(f"❌ Error calculating stain vectors: {e}")
            return False

    def extract_hematoxylin(self, I: np.ndarray, Io: int = 240, beta: float = 0.15) -> np.ndarray:
        """
        Extract Hematoxylin channel from RGB image using fitted stain matrix.
        Fallback to fit on the fly if not fitted, or grayscale if that fails.
        """
        # If already grayscale, just return as is (normalized to 8-bit if needed)
        if I.ndim == 2:
            if I.dtype == np.uint8: return I
            return cv2.normalize(I, None, 0, 255, cv2.NORM_MINMAX).astype(np.uint8)

        h, w = I.shape[:2]

        if self.stain_matrix is None:
            success = self.fit(I, Io, beta)
            if not success:
               print("⚠️ Macenko normalization not fitted and fit failed. Falling back to grayscale conversion.")
               if I.ndim == 3:
                   if I.shape[-1] == 4: I = I[..., :3]
                   if I.shape[-1] == 3: return cv2.cvtColor(I, cv2.COLOR_RGB2GRAY)
                   return I[..., 0]
               return I

        if I.ndim == 3 and I.shape[-1] == 4:
             I = I[..., :3]

        try:
            OD = -np.log((I.reshape((-1, 3)).astype(np.float64) + 1) / (Io + 1))
            C = np.linalg.lstsq(self.stain_matrix.T, OD.T, rcond=None)[0].T
            H_conc = C[:, 0].reshape(h, w)
            H_norm = cv2.normalize(H_conc, None, 0, 255, cv2.NORM_MINMAX)
            return H_norm.astype(np.uint8)

        except Exception as e:
             print(f"❌ Error extracting concentration: {e}. Falling back to grayscale.")
             if I.ndim == 3:
                if I.shape[-1] == 3: return cv2.cvtColor(I, cv2.COLOR_RGB2GRAY)
                return I[..., 0]
             return I


def apply_clahe(img: np.ndarray, clip_limit: float = 2.0, grid_size: tuple = (8, 8)) -> np.ndarray:
    """Apply CLAHE (Contrast Limited Adaptive Histogram Equalization) to an image."""
    clahe = cv2.createCLAHE(clipLimit=clip_limit, tileGridSize=grid_size)
    return clahe.apply(img)
