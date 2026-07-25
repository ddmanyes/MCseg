"""
enact_C3b_proseg_seg.py
-----------------------
ProSeg segmentation on the ENACT CRC crop region (Visium HD 2µm bins).

Generates proseg_mask.npy (shape: CROP_H × CROP_W) in the ENACT crop local
coordinate system, compatible with enact_C3_method_comparison.py.

Method:
  1. Load filtered_feature_bc_matrix.h5 + tissue_positions.parquet
  2. Filter barcodes to ENACT crop (BTF_COL_MIN–MAX, CROP_Y0–Y1)
  3. Expand count matrix → transcript CSV (each count = 1 row)
     Each transcript placed at bin centroid ± jitter within one bin width
  4. Run ProSeg with --visiumhd --coordinate-scale 0.2737 (px → µm)
  5. Rasterize output cell polygons → proseg_mask.npy

Note: Visium HD 2µm bins lack sub-bin spatial resolution, so transcript
positions are approximated by bin centroids with random jitter. This is a
known limitation described in the manuscript Methods.

Usage:
    cd /Volumes/SSD/plan_a
    uv run python submission_bioinformatics/scripts/analysis/enact_crc_benchmark/enact_C3b_proseg_seg.py
"""

import logging
import subprocess
import sys
import time
from pathlib import Path

import cv2
import h5py
import numpy as np
import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(asctime)s  %(message)s",
                    datefmt="%H:%M:%S")
log = logging.getLogger(__name__)

# ─── Paths ────────────────────────────────────────────────────────────────────
PLAN_A     = Path(__file__).resolve().parents[4]
ENACT_CRC  = PLAN_A / "tissue sample" / "ENACT_supporting_files" / "public_data" / "human_colorectal"
H5_PATH    = ENACT_CRC / "input_files" / "filtered_feature_bc_matrix.h5"
TP_PATH    = ENACT_CRC / "input_files" / "tissue_positions.parquet"
RESULT_DIR = PLAN_A / "submission_bioinformatics" / "results" / "enact_method_comparison"
PROSEG_OUT = RESULT_DIR / "proseg_out"
TRANSCRIPT_CSV     = PROSEG_OUT / "transcripts_enact_crop.csv"
TRANSCRIPT_CSV_XFM = PROSEG_OUT / "transcripts_xenium_fmt.csv"  # Xenium-format version
POLYGON_CSV        = PROSEG_OUT / "proseg_cell_polygons.csv"
MASK_PATH      = RESULT_DIR / "proseg_mask.npy"
NUC_MASK_PATH  = RESULT_DIR / "nuc_mask.npy"   # Cellpose init for ProSeg
PROSEG_BIN     = Path.home() / ".cargo" / "bin" / "proseg"

# ─── ENACT crop constants (identical to enact_C3_method_comparison.py) ────────
CROP_X0, CROP_X1 = 5154, 15242
CROP_Y0, CROP_Y1 = 4635, 18599
COL_OFFSET  = 40598
BTF_COL_MIN = CROP_X0 + COL_OFFSET   # 45752
BTF_COL_MAX = CROP_X1 + COL_OFFSET   # 55840
CROP_W = CROP_X1 - CROP_X0           # 10088
CROP_H = CROP_Y1 - CROP_Y0           # 13964

# Visium HD fullres pixel size — used by --coordinate-scale to convert px → µm
VHD_PIXEL_UM = 0.2737
# Jitter half-width: ±half of 2µm bin spacing in fullres pixels (~3.65 px)
BIN_HALF_PX = 3.5

# ─── ProSeg CLI parameters ───────────────────────────────────────────────────
PROSEG_FLAGS = [
    "--coordinate-scale", str(VHD_PIXEL_UM),
    "--voxel-size",        "2",   # match 2µm bin resolution
    "--burnin-voxel-size", "4",
    "--cell-compactness",  "0.06",
    "--samples",           "100",
    "--recorded-samples",  "50",
    "--enforce-connectivity",
    "--overwrite",
]


# ─── Step 1: Build transcript CSV ────────────────────────────────────────────

def build_transcript_csv(rng: np.random.Generator) -> None:
    if TRANSCRIPT_CSV.exists():
        log.info(f"Transcript CSV already cached ({TRANSCRIPT_CSV.name}) — skipping.")
        return

    log.info("Loading tissue positions...")
    tp = pd.read_parquet(str(TP_PATH))
    in_crop = (
        (tp["pxl_col_in_fullres"] >= BTF_COL_MIN) & (tp["pxl_col_in_fullres"] < BTF_COL_MAX) &
        (tp["pxl_row_in_fullres"] >= CROP_Y0)      & (tp["pxl_row_in_fullres"] < CROP_Y1)
    )
    tp_crop = tp[in_crop].reset_index(drop=True)
    log.info(f"  {len(tp_crop):,} barcodes in ENACT crop")

    log.info("Loading h5 count matrix...")
    with h5py.File(str(H5_PATH), "r") as f:
        all_barcodes = [b.decode() if isinstance(b, bytes) else b
                        for b in f["matrix/barcodes"][:]]
        gene_names   = [g.decode() if isinstance(g, bytes) else g
                        for g in f["matrix/features/name"][:]]
        indptr  = f["matrix/indptr"][:]
        indices = f["matrix/indices"][:]
        data    = f["matrix/data"][:]

    barcode_to_idx = {b: i for i, b in enumerate(all_barcodes)}
    crop_barcodes  = tp_crop["barcode"].tolist()
    crop_col_px    = tp_crop["pxl_col_in_fullres"].values
    crop_row_px    = tp_crop["pxl_row_in_fullres"].values

    log.info("Expanding count matrix → transcript CSV (~3–5 min)...")
    PROSEG_OUT.mkdir(parents=True, exist_ok=True)

    t0 = time.time()
    CHUNK = 100_000
    total_written = 0

    with open(str(TRANSCRIPT_CSV), "w") as fout:
        fout.write("x_loc,y_loc,gene\n")

        for chunk_start in range(0, len(crop_barcodes), CHUNK):
            chunk_end = min(chunk_start + CHUNK, len(crop_barcodes))
            rows_x: list[float] = []
            rows_y: list[float] = []
            rows_gene: list[str] = []

            for k in range(chunk_start, chunk_end):
                bc     = crop_barcodes[k]
                bc_idx = barcode_to_idx.get(bc)
                if bc_idx is None:
                    continue
                x_base = crop_col_px[k] - BTF_COL_MIN
                y_base = crop_row_px[k] - CROP_Y0

                start, end = int(indptr[bc_idx]), int(indptr[bc_idx + 1])
                if start == end:
                    continue
                gene_idx_arr = indices[start:end]
                count_arr    = data[start:end].astype(np.int32)
                total_umi    = int(count_arr.sum())
                if total_umi == 0:
                    continue

                jx = rng.uniform(-BIN_HALF_PX, BIN_HALF_PX, size=total_umi)
                jy = rng.uniform(-BIN_HALF_PX, BIN_HALF_PX, size=total_umi)

                t_offset = 0
                for gi, cnt in zip(gene_idx_arr, count_arr):
                    gene_str = gene_names[gi]
                    for _ in range(int(cnt)):
                        rows_x.append(round(x_base + jx[t_offset], 3))
                        rows_y.append(round(y_base + jy[t_offset], 3))
                        rows_gene.append(gene_str)
                        t_offset += 1

            if rows_x:
                chunk_df = pd.DataFrame({
                    "x_loc": rows_x,
                    "y_loc": rows_y,
                    "gene":  rows_gene,
                })
                fout.write(chunk_df.to_csv(index=False, header=False))
                total_written += len(rows_x)

            if (chunk_start // CHUNK) % 10 == 0:
                pct     = 100 * chunk_end / len(crop_barcodes)
                elapsed = time.time() - t0
                log.info(f"  {chunk_end:,}/{len(crop_barcodes):,} barcodes ({pct:.0f}%)  "
                         f"{total_written:,} rows  {elapsed:.0f}s")

    elapsed = time.time() - t0
    log.info(f"  Transcript CSV done: {total_written:,} rows in {elapsed:.0f}s")


# ─── Step 2: Run ProSeg ──────────────────────────────────────────────────────

def run_proseg() -> bool:
    if not PROSEG_BIN.exists():
        log.error(f"ProSeg binary not found at {PROSEG_BIN}")
        return False

    if not NUC_MASK_PATH.exists():
        log.error(f"Cellpose init mask not found: {NUC_MASK_PATH}")
        return False

    # ProSeg requires uint32 — convert if needed
    nuc_mask_u32 = PROSEG_OUT / "nuc_mask_uint32.npy"
    if not nuc_mask_u32.exists():
        log.info("Converting nuc_mask int32 → uint32 for ProSeg...")
        m = np.load(str(NUC_MASK_PATH))
        np.save(str(nuc_mask_u32), m.astype(np.uint32))
        del m
        log.info(f"  Saved: {nuc_mask_u32}")

    cmd = (
        [str(PROSEG_BIN), str(TRANSCRIPT_CSV_XFM)]
        + ["--gene-column",    "feature_name"]
        + ["--x-column",       "x_location"]
        + ["--y-column",       "y_location"]
        + ["--z-column",       "z_location"]
        + ["--qv-column",      "qv"]
        + ["--cell-id-column", "cell_id"]
        + ["--cell-id-unassigned", "0"]
        + ["--cellpose-masks", str(nuc_mask_u32)]
        + ["--cellpose-scale", str(VHD_PIXEL_UM)]
        + ["--output-path",         str(PROSEG_OUT) + "/"]
        + ["--output-cell-polygons", str(POLYGON_CSV)]
        + ["--output-cell-metadata", str(PROSEG_OUT / "cells.csv")]
        + ["--output-counts",        str(PROSEG_OUT / "counts.csv")]
        + ["--output-counts-fmt",    "csv"]
        + PROSEG_FLAGS
    )

    log.info("Running ProSeg (~40–60 min for 74M transcripts)...")
    log.info(f"  {' '.join(cmd[:6])} ...")
    t0 = time.time()
    result = subprocess.run(cmd)
    elapsed = time.time() - t0

    if result.returncode != 0:
        log.error(f"ProSeg failed (exit {result.returncode}) after {elapsed/60:.1f} min")
        return False

    log.info(f"ProSeg completed in {elapsed/60:.1f} min")
    return True


# ─── Step 3: Rasterize polygons → mask ───────────────────────────────────────

def rasterize_polygons() -> np.ndarray:
    log.info(f"Rasterizing ProSeg polygons → {CROP_H}×{CROP_W} mask...")

    if not POLYGON_CSV.exists():
        raise FileNotFoundError(
            f"Expected polygon output not found: {POLYGON_CSV}\n"
            f"ProSeg output dir contents: {[f.name for f in PROSEG_OUT.iterdir() if not f.name.startswith('._')]}"
        )

    import gzip, json

    with gzip.open(str(POLYGON_CSV), "rt") as fh:
        geojson = json.load(fh)

    features = geojson["features"]
    log.info(f"  GeoJSON FeatureCollection: {len(features):,} features")

    mask = np.zeros((CROP_H, CROP_W), dtype=np.int32)
    n_ok = 0

    for feat in features:
        cell_id = feat["properties"]["cell"] + 1  # 0-indexed → 1-indexed
        geom    = feat["geometry"]
        gtype   = geom["type"]

        # coords are in µm; convert back to pixels
        if gtype == "Polygon":
            rings = geom["coordinates"]
        elif gtype == "MultiPolygon":
            rings = [ring for poly in geom["coordinates"] for ring in poly]
        else:
            continue

        for ring in rings:
            pts = (np.array(ring, dtype=np.float64)[:, :2] / VHD_PIXEL_UM).astype(np.int32)
            pts = np.clip(pts, [0, 0], [CROP_W - 1, CROP_H - 1])
            pts = pts.reshape(-1, 1, 2)
            cv2.fillPoly(mask, [pts], color=cell_id)
        n_ok += 1

    log.info(f"  Rasterized {n_ok:,} cells  (mask max cell_id={mask.max():,})")
    return mask


# ─── Main ────────────────────────────────────────────────────────────────────

def main():
    RESULT_DIR.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(42)

    build_transcript_csv(rng)
    if not run_proseg():
        sys.exit(1)

    mask = rasterize_polygons()
    np.save(str(MASK_PATH), mask)
    log.info(f"Saved: {MASK_PATH}  shape={mask.shape}  non-zero px={int((mask > 0).sum()):,}")


if __name__ == "__main__":
    main()
