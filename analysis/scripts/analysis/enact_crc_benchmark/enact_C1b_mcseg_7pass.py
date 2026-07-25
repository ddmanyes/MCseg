"""
enact_C1b_mcseg_7pass.py
------------------------
Re-run MCseg v2 with full 7-pass configuration (use_cpsam=True) on the
ENACT CRC H&E crop and save to mcseg_mask_7pass.npy.

7-pass = cyto3 × 3 diameters + hematoxylin + cpsam_auto + cpsam_dia16 + cpsam_hema
(4-pass deployment mode omits the 3 cpsam passes)

Usage:
    cd /Volumes/SSD/plan_a
    uv run python submission_bioinformatics/scripts/analysis/enact_crc_benchmark/enact_C1b_mcseg_7pass.py
"""

import logging
import sys
import time
from pathlib import Path

import numpy as np
import tifffile

logging.basicConfig(level=logging.INFO, format="%(asctime)s  %(message)s",
                    datefmt="%H:%M:%S")
log = logging.getLogger(__name__)

# ─── Paths ────────────────────────────────────────────────────────────────────
PLAN_A      = Path("/Volumes/SSD/plan_a")
MSSEG_ROOT  = PLAN_A / "MSseg"
F1_DIR      = PLAN_A / "submission_bioinformatics" / "results" / "enact_crc_f1"
HE_CROP     = F1_DIR / "he_crop.tif"
OUT_MASK    = F1_DIR / "mcseg_mask_7pass.npy"


def main() -> None:
    if not HE_CROP.exists():
        log.error(f"H&E crop not found: {HE_CROP}")
        log.error("Run enact_C1_crc_f1_comparison.py first to generate he_crop.tif")
        sys.exit(1)

    if OUT_MASK.exists():
        log.info(f"7-pass mask already exists: {OUT_MASK.name} — skipping")
        mask = np.load(str(OUT_MASK))
        log.info(f"  shape={mask.shape}  cells={int(mask.max()):,}")
        return

    log.info(f"Loading H&E crop: {HE_CROP.name}")
    img = tifffile.imread(str(HE_CROP))
    if img.ndim == 3 and img.shape[-1] == 4:
        img = img[..., :3]
    log.info(f"  shape={img.shape}")

    sys.path.insert(0, str(MSSEG_ROOT / "backend"))
    from src.segmentation.cellpose_runner import run_tiled_mcseg_v2

    cfg = {
        "use_gpu":               True,
        "batch_size":            2,
        "dia_small":             13.0,
        "dia_mid":               17.0,
        "dia_large":             22.0,
        "use_hematoxylin":       True,
        "use_cpsam":             True,    # ← enables 3 extra cpsam passes → 7-pass total
        "voronoi_distance":      8,
        "flow_threshold":        0.4,
        "cellprob_threshold":    -2.0,
        "min_size":              20,
        "max_size":              6000,
        "clahe_clip_limit":      3.0,
        "use_transcript_rescue": False,
    }

    log.info("Running 7-pass MCseg v2 (use_cpsam=True) — expect ~2–3× longer than 4-pass")

    def _progress(p: float, msg: str) -> None:
        log.info(f"  [{p*100:.0f}%] {msg}")

    t0 = time.time()
    mask = run_tiled_mcseg_v2(
        img,
        cfg,
        tile_size=1024,
        overlap=128,
        progress_callback=_progress,
    )
    elapsed = time.time() - t0

    np.save(str(OUT_MASK), mask)
    log.info(f"\nDone in {elapsed/60:.1f} min")
    log.info(f"Saved: {OUT_MASK.name}  shape={mask.shape}  cells={int(mask.max()):,}")
    coverage = float((mask > 0).mean())
    log.info(f"Pixel coverage: {coverage:.1%}")


if __name__ == "__main__":
    main()
