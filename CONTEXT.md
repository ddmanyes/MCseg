# MCseg Spatial Transcriptomics Segmentation & Analysis Pipeline

High-resolution spatial transcriptomics pipeline integrating multi-pass deep learning cell segmentation (Cellpose / CPSAM), spatial bin-to-cell attribution, and single-cell downstream analysis for 10x Visium HD and multimodal assays.

## Coordinate Spaces

**RawTiffPixel**:
Pixel coordinate in the full-resolution H&E or IF slide image (typically scanned at ~0.2737 µm/px or higher). Used for image cropping, DZI tile rendering, and segmentation input.
_Avoid_: Generic pixel, image coord, fullres

**SpaceRangerPixel**:
Coordinate system output by Space Ranger in `tissue_positions.parquet` (`pxl_col_in_fullres`, `pxl_row_in_fullres`). May have a different scale or translation offset relative to the raw slide scan.
_Avoid_: Fullres px, array coordinate

**MicronCoordinate**:
Physical distance measurement in micrometers (µm) based on slide metadata `microns_per_pixel`. Used for cross-modality comparison (e.g. Visium HD vs Xenium).
_Avoid_: Physical pixel, physical coord

**CaptureAreaBoundingBox**:
The physical bounding rectangle of active Visium HD sequencing probes (typically 6.5 mm × 6.5 mm) mapped to RawTiffPixel space. Regions outside this box contain only histological tissue with zero sequencing data.
_Avoid_: Slide boundary, image bounds

## Spatial & Biological Entities

**Bin**:
A regular spatial grid unit (2 µm, 8 µm, or 16 µm square) defined by Visium HD sequencing barcodes.
_Avoid_: Spot, pixel, barcode point

**CellMask**:
An integer-labeled 2D array or vector polygon representing the segmented boundary of an individual whole cell.
_Avoid_: Segmentation boundary, outline, contour

**NuclearMask**:
An integer-labeled mask representing only the cell nucleus (derived from StarDist, DAPI, or Hematoxylin extraction).
_Avoid_: Cell mask, NUC label

**CellAttribution**:
The spatial mapping process that assigns sequencing Bins and their transcript counts to overlapping CellMask labels.
_Avoid_: Counting, bin sum, simple allocation

## Pipeline Lifecycle

**Stage0_RoiCrop**:
Interactive region-of-interest selection on whole-slide images, extracting paired H&E image crops and AnnData bin subsets with unified coordinate projection.
_Avoid_: Crop step, slicing

**Stage1_EnsembleSegmentation**:
Multi-pass Cellpose and CPSAM segmentation across multiple diameters and color spaces, merged with Voronoi cell expansion.
_Avoid_: Cellpose runner, segmentation pass

**Stage2_SpatialCounting**:
CellAttribution aggregation generating the cell-by-gene AnnData matrix (`cellpose_cells.h5ad`).
_Avoid_: Matrix aggregation, bin counter

**Stage3_DownstreamAnalysis**:
Quality control (QC) filtering, UMAP dimensionality reduction, Leiden clustering, and CellTypist cell-type annotation.
_Avoid_: Scanpy runner, clustering

**TierClassification**:
Hierarchical cell-type labeling structure: Tier 1 (Broad tissue type), Tier 2 (Specific cell lineage), Tier 3 (Functional state/subtype).
_Avoid_: Cluster name, cell label
