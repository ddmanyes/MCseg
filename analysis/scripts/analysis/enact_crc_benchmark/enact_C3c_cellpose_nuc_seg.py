"""
enact_C3c_cellpose_nuc_seg.py
------------------------------
Cellpose nuclei segmentation on the ENACT CRC H&E crop.
Generates cellpose_nuc_mask.npy — pixel-level nuclear segmentation
to replace the Space Ranger bin-grid nuc_mask.npy.

CP4 only has cpsam; nucleitorch_0 (CP3) is incompatible.
Uses cpsam on the hematoxylin channel (Ruifrok colour deconvolution)
which isolates nuclear stain for clean nucleus detection.

Parameters tuned for CRC H&E at 0.2737 µm/px:
  - model: cpsam (CP4)
  - input: hematoxylin channel from colour deconvolution
  - diameter: 20px (~5.5µm, typical CRC nucleus)
  - tiled: 2048px tiles with 256px overlap

Usage:
    cd /Volumes/SSD/plan_a
    /Users/zhanqiru/.venvs/msseg/bin/python \
        submission_bioinformatics/scripts/analysis/enact_crc_benchmark/enact_C3c_cellpose_nuc_seg.py
"""

import logging
import time
from pathlib import Path

import numpy as np
import tifffile
from cellpose import models

# ── Hematoxylin colour deconvolution (Ruifrok & Johnston) ────────────────────
# Same stain matrix as MSseg/backend/src/segmentation/cellpose_runner.py
_HE_STAIN = np.array([
    [0.65, 0.70, 0.29],   # Hematoxylin
    [0.07, 0.99, 0.11],   # Eosin
    [0.27, 0.57, 0.78],   # DAB (residual)
], dtype=np.float64)
_HE_INV = np.linalg.inv(_HE_STAIN)


def extract_hematoxylin(rgb: np.ndarray) -> np.ndarray:
    """Return hematoxylin channel as uint8 (bright = high concentration)."""
    rgb_f = rgb.astype(np.float64) / 255.0 + 1e-6
    od = -np.log(np.clip(rgb_f, 1e-6, 1.0))
    h_od = (od @ _HE_INV[:, 0:1]).squeeze()
    h_od = np.clip(h_od, 0, None)
    h_norm = (h_od / (h_od.max() + 1e-8) * 255).astype(np.uint8)
    return h_norm

logging.basicConfig(level=logging.INFO, format="%(asctime)s  %(message)s",
                    datefmt="%H:%M:%S")
log = logging.getLogger(__name__)

# ─── Paths ────────────────────────────────────────────────────────────────────
PLAN_A     = Path("/Volumes/SSD/plan_a")
F1_DIR     = PLAN_A / "submission_bioinformatics" / "results" / "enact_crc_f1"
RESULT_DIR = PLAN_A / "submission_bioinformatics" / "results" / "enact_method_comparison"
HE_CROP    = F1_DIR / "he_crop.tif"
OUT_MASK   = RESULT_DIR / "cellpose_nuc_mask.npy"

# ─── Segmentation params ──────────────────────────────────────────────────────
NUC_DIAMETER    = 20     # px ≈ 5.5µm at 0.2737µm/px
TILE_SIZE       = 2048
OVERLAP         = 256
USE_GPU         = True
FLOW_THRESH     = 0.4
CELLPROB_THRESH = 0.0    # nuclei model default


def run_tiled_nuclei(img: np.ndarray) -> np.ndarray:
    H, W = img.shape[:2]
    log.info(f"Image: {W}×{H}px  tile={TILE_SIZE}px  overlap={OVERLAP}px")

    # CP4: only cpsam is available; nucleitorch_0 is CP3-only
    nuc_model = models.CellposeModel(model_type="cpsam", gpu=USE_GPU)

    global_mask = np.zeros((H, W), dtype=np.int32)
    cell_offset = 0

    steps_y = list(range(0, H, TILE_SIZE - OVERLAP))
    steps_x = list(range(0, W, TILE_SIZE - OVERLAP))
    total_tiles = len(steps_y) * len(steps_x)
    tile_idx = 0

    for y0 in steps_y:
        y1 = min(H, y0 + TILE_SIZE)
        for x0 in steps_x:
            x1 = min(W, x0 + TILE_SIZE)
            tile_idx += 1
            t0 = time.time()

            tile_rgb = img[y0:y1, x0:x1]
            # Use hematoxylin channel: isolates nuclear stain, suppresses eosin
            tile_hema = extract_hematoxylin(tile_rgb)

            masks_pred, _, _ = nuc_model.eval(
                tile_hema,
                diameter=NUC_DIAMETER,
                flow_threshold=FLOW_THRESH,
                cellprob_threshold=CELLPROB_THRESH,
            )

            n_cells = int(masks_pred.max())
            elapsed = time.time() - t0
            log.info(f"  Tile {tile_idx}/{total_tiles} ({x0},{y0})  cells={n_cells}  {elapsed:.0f}s")

            if n_cells == 0:
                continue

            # Write core region (exclude overlap border except at image edges)
            oy0 = OVERLAP // 2 if y0 > 0 else 0
            oy1 = (y1 - y0) - OVERLAP // 2 if y1 < H else (y1 - y0)
            ox0 = OVERLAP // 2 if x0 > 0 else 0
            ox1 = (x1 - x0) - OVERLAP // 2 if x1 < W else (x1 - x0)

            core = masks_pred[oy0:oy1, ox0:ox1]
            nonzero = core > 0
            global_mask[y0+oy0:y0+oy1, x0+ox0:x0+ox1][nonzero] = core[nonzero] + cell_offset
            cell_offset += n_cells

    # Remap accumulated IDs to contiguous 1..N WITHOUT merging touching cells.
    # label(>0) would merge adjacent cells from different tiles into one component.
    log.info("Remapping to contiguous IDs...")
    unique_ids = np.unique(global_mask)
    unique_ids = unique_ids[unique_ids > 0]
    remap = np.zeros(int(global_mask.max()) + 1, dtype=np.int32)
    for new_id, old_id in enumerate(unique_ids, start=1):
        remap[old_id] = new_id
    global_mask = remap[global_mask]
    log.info(f"Final cell count: {int(global_mask.max()):,}")
    return global_mask


def main() -> None:
    if OUT_MASK.exists():
        log.info(f"Already exists: {OUT_MASK.name} — skipping")
        m = np.load(str(OUT_MASK))
        log.info(f"  shape={m.shape}  cells={int(m.max()):,}")
        return

    log.info(f"Loading H&E crop: {HE_CROP.name}")
    img = tifffile.imread(str(HE_CROP))
    if img.ndim == 3 and img.shape[-1] == 4:
        img = img[..., :3]
    log.info(f"  shape={img.shape}")

    t0 = time.time()
    mask = run_tiled_nuclei(img)
    elapsed = time.time() - t0

    np.save(str(OUT_MASK), mask)
    log.info(f"\nDone in {elapsed/60:.1f} min")
    log.info(f"Saved: {OUT_MASK.name}  shape={mask.shape}  cells={int(mask.max()):,}")
    log.info(f"Pixel coverage: {float((mask>0).mean()):.1%}")


if __name__ == "__main__":
    main()
