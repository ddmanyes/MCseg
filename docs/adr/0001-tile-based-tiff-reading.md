# 1. Tile-based TIFF/BTF Reading for Slide Image Access

Full-slide histological images (WSI) range from 10 GB to over 80 GB with dimensions exceeding 50,000 × 50,000 pixels. We decided to strictly require tile-based or strip-based random access reading (`read_btf_crop` / `PyramidSlideReader`) and forbid whole-image loading (`page.asarray()`) anywhere in the pipeline to guarantee bounded memory usage and prevent out-of-memory (OOM) crashes.
