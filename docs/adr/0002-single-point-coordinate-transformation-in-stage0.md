# 2. Single-Point Coordinate Transformation in Stage 0

Space Ranger pixel coordinates (`virtual_fullres`) and raw slide image pixels (`RawTiffPixel`) often differ in resolution and offset. We decided that coordinate projection must happen exclusively in Stage 0 (`subset_anndata_roi`) and be persisted directly to `adata.obsm['spatial']`. Downstream stages (Stage 2 CellAttribution, Stage 3 Analysis, Stage 4 Export) must consume `obsm['spatial']` as-is and never re-apply transformations.
