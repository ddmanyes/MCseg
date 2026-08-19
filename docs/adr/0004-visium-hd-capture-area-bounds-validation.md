# 4. Visium HD Capture Area Bounds Validation at Stage 0

Visium HD capture chips (6.5 mm × 6.5 mm) only cover a subset of whole-slide histological scans, meaning large areas of tissue have zero sequencing probes. We decided to compute the active RNA capture bounding box during Stage 0 and proactively reject or warn on ROIs with zero overlapping bins. This prevents downstream CellAttribution from generating empty cells that cause silent or confusing single-cell QC crashes.
