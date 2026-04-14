"""
segment_template.py — AutoResearch Sandbox Starter Script

This is the file the AI agent modifies freely during autonomous discovery.
The only constraint: preserve the `build_and_predict` function signature.

Adapt `prepare.py` in your own project to load your image and ground truth,
then replace this starter implementation with your baseline.
"""

import cv2
import numpy as np
from cellpose import models
from scipy.ndimage import label


# ── Helper functions (agent may add, remove, or replace these) ────────────────

def apply_clahe(img: np.ndarray, clip_limit: float = 3.0, tile_size: int = 8) -> np.ndarray:
    """CLAHE contrast enhancement in LAB colour space."""
    lab = cv2.cvtColor(img, cv2.COLOR_RGB2LAB)
    l, a, b = cv2.split(lab)
    clahe = cv2.createCLAHE(clipLimit=clip_limit, tileGridSize=(tile_size, tile_size))
    cl = clahe.apply(l)
    return cv2.cvtColor(cv2.merge((cl, a, b)), cv2.COLOR_LAB2RGB)


def extract_hematoxylin(img: np.ndarray) -> np.ndarray:
    """HED colour deconvolution — returns hematoxylin channel as uint8."""
    img_float = img.astype(np.float64) + 1.0
    od = -np.log(img_float / 256.0)
    he_matrix = np.array([
        [0.6500286, 0.7041680, 0.2860126],
        [0.0728940, 0.9904310, 0.1155140],
        [0.2688350, 0.5706770, 0.7768750],
    ])
    for i in range(3):
        norm = np.linalg.norm(he_matrix[i])
        if norm > 0:
            he_matrix[i] /= norm
    stains = od.reshape(-1, 3) @ np.linalg.inv(he_matrix).T
    hema = np.clip(stains.reshape(img.shape)[:, :, 0], 0, None)
    h_max = np.percentile(hema, 99.5)
    if h_max > 0:
        hema = np.clip(hema / h_max, 0, 1)
    return (hema * 255).astype(np.uint8)


def voronoi_expand(mask: np.ndarray, max_distance: int = 8) -> np.ndarray:
    """Voronoi-constrained boundary expansion: unclaimed pixels assigned to
    nearest cell centroid, up to max_distance pixels."""
    from scipy.ndimage import distance_transform_edt
    binary = mask > 0
    if not binary.any():
        return mask.copy()
    dist_from_cell, nearest_idx = distance_transform_edt(~binary, return_indices=True)
    expanded = mask.copy()
    in_range = (dist_from_cell > 0) & (dist_from_cell <= max_distance)
    expanded[in_range] = mask[nearest_idx[0][in_range], nearest_idx[1][in_range]]
    return expanded


def merge_masks(masks: list[np.ndarray], overlap_threshold: float = 0.15) -> np.ndarray:
    """Greedy non-maximum suppression: merge multiple instance masks,
    resolving overlaps by retaining whichever mask was added first."""
    if not masks:
        return np.zeros((1, 1), dtype=np.int32)
    h, w = masks[0].shape
    merged = np.zeros((h, w), dtype=np.int32)
    next_id = 1
    for m in masks:
        for cell_id in np.unique(m):
            if cell_id == 0:
                continue
            cell_pixels = m == cell_id
            overlap = np.sum(cell_pixels & (merged > 0)) / np.sum(cell_pixels)
            if overlap < overlap_threshold:
                merged[cell_pixels & (merged == 0)] = next_id
                next_id += 1
    return merged


# ── Main segmentation function (preserve this signature exactly) ──────────────

def build_and_predict(img: np.ndarray, vhd_csv: str, gt_mask=None) -> np.ndarray:
    """
    Perform cell instance segmentation on an H&E image patch.

    Parameters
    ----------
    img : np.ndarray
        RGB image array, shape (H, W, 3), dtype uint8.
    vhd_csv : str
        Path to Visium HD pseudo-transcript CSV (columns: x, y, gene).
        May be empty/unused depending on the approach.
    gt_mask : np.ndarray, optional
        Ground-truth instance mask for intermediate scoring during development.

    Returns
    -------
    np.ndarray
        2D int32 array: 0 = background, >0 = cell instance ID.
    """
    # ── Step 1: Preprocessing ─────────────────────────────────────────────────
    enhanced = apply_clahe(img, clip_limit=3.0, tile_size=8)
    hema = extract_hematoxylin(img)

    # ── Step 2: Multi-pass detection ──────────────────────────────────────────
    # Agent is free to change models, diameters, and number of passes.
    model_cyto3 = models.CellposeModel(pretrained_model='cyto3', gpu=True)
    model_cpsam = models.CellposeModel(pretrained_model='cpsam', gpu=True)

    all_masks = []
    for diameter in [13, 17, 22]:
        m, _, _ = model_cyto3.eval(enhanced, diameter=diameter,
                                   flow_threshold=0.4, cellprob_threshold=0.0)
        all_masks.append(m)

    for diameter in [0, 16]:  # 0 = auto
        m, _, _ = model_cpsam.eval(enhanced, diameter=diameter,
                                   flow_threshold=0.4, cellprob_threshold=0.0,
                                   augment=False, resample=False)
        m = m[:img.shape[0], :img.shape[1]]  # guard against resample size mismatch
        all_masks.append(m)

    m_hema, _, _ = model_cpsam.eval(hema, diameter=0,
                                    flow_threshold=0.4, cellprob_threshold=0.0,
                                    augment=False, resample=False)
    m_hema = m_hema[:img.shape[0], :img.shape[1]]
    all_masks.append(m_hema)

    # ── Step 3: Ensemble merging ──────────────────────────────────────────────
    merged = merge_masks(all_masks, overlap_threshold=0.15)

    # ── Step 4: Voronoi boundary expansion ───────────────────────────────────
    final_mask = voronoi_expand(merged, max_distance=8)

    # ── Step 5: Quality filtering ─────────────────────────────────────────────
    for cell_id in np.unique(final_mask):
        if cell_id == 0:
            continue
        area = np.sum(final_mask == cell_id)
        if area < 20 or area > 6000:
            final_mask[final_mask == cell_id] = 0

    return final_mask.astype(np.int32)


# ── Entry point for scoring ───────────────────────────────────────────────────
if __name__ == "__main__":
    import prepare  # your project's data loading + evaluation script
    img, gt_mask, vhd_csv = prepare.load_data()
    pred_mask = build_and_predict(img, vhd_csv, gt_mask=gt_mask)
    score = prepare.evaluate_iou(pred_mask, gt_mask)
