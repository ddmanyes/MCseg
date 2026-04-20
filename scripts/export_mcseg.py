"""
MCseg 匯出 CLI
從 MCseg 分割遮罩 + AnnData 產出：
  --format xenium  → results/export/xenium/ (GeoJSON + zarr)
  --format h5ad    → results/export/msseg_final.h5ad
  --format both    → 兩者皆產出（預設）

不修改任何現有後端檔案。
"""
from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path

import anndata as ad
import numpy as np

logger = logging.getLogger(__name__)


# ── GeoJSON 輪廓提取 ──────────────────────────────────────────────────────────

def _masks_to_geojson(
    masks_dir:     Path,
    pixel_size_um: float = 0.2737,
) -> dict:
    """從所有 _mcseg.npy 遮罩提取細胞輪廓，轉為 GeoJSON（µm 座標）。"""
    from skimage.measure import find_contours

    features = []
    cell_global_id = 1

    for mask_path in sorted(masks_dir.glob("*_mcseg.npy")):
        mask = np.load(mask_path)
        n_cells = int(mask.max())

        for cell_id in range(1, n_cells + 1):
            cell_binary = (mask == cell_id).astype(np.uint8)
            contours = find_contours(cell_binary, level=0.5)
            if not contours:
                continue
            # 取最大輪廓
            contour = max(contours, key=len)
            # row,col → x,y (µm)
            coords = [[float(c[1]) * pixel_size_um, float(c[0]) * pixel_size_um]
                      for c in contour]
            if len(coords) < 3:
                continue
            coords.append(coords[0])  # 閉合多邊形

            features.append({
                "type": "Feature",
                "id": cell_global_id,
                "geometry": {
                    "type": "Polygon",
                    "coordinates": [coords],
                },
                "properties": {"cell_id": cell_global_id},
            })
            cell_global_id += 1

    return {"type": "FeatureCollection", "features": features}


# ── Xenium Explorer 匯出 ──────────────────────────────────────────────────────

def export_xenium(
    adata:         ad.AnnData,
    masks_dir:     Path,
    output_dir:    Path,
    pixel_size_um: float = 0.2737,
) -> Path:
    """
    產出 Xenium Explorer 相容的 GeoJSON + cell_feature_matrix.zarr。
    不依賴也不修改 XeniumExporter（設計給 Proseg）。
    """
    import zarr
    import scipy.sparse as sp

    output_dir.mkdir(parents=True, exist_ok=True)

    # 1. GeoJSON 輪廓
    geojson = _masks_to_geojson(masks_dir, pixel_size_um)
    geojson_path = output_dir / "cell_boundaries.geojson"
    geojson_path.write_text(
        json.dumps(geojson, indent=None, ensure_ascii=False),
        encoding="utf-8",
    )
    print(f"  GeoJSON：{len(geojson['features'])} cells → {geojson_path}")

    # 2. cell_feature_matrix.zarr（稀疏矩陣）
    zarr_path = output_dir / "cell_feature_matrix.zarr"
    store     = zarr.open(str(zarr_path), mode="w")

    X = adata.X
    if sp.issparse(X):
        X_csr = X.tocsr()
    else:
        X_csr = sp.csr_matrix(X)

    store.create_dataset("data",    data=X_csr.data,    overwrite=True)
    store.create_dataset("indices", data=X_csr.indices, overwrite=True)
    store.create_dataset("indptr",  data=X_csr.indptr,  overwrite=True)
    store.attrs["shape"]     = list(X_csr.shape)
    store.attrs["barcodes"]  = list(adata.obs_names)
    store.attrs["features"]  = list(adata.var_names)
    store.attrs["pixel_size_um"] = pixel_size_um
    print(f"  zarr：{X_csr.shape} → {zarr_path}")

    # 3. experiment.xenium 元數據
    meta = {
        "cell_count":    int(adata.n_obs),
        "gene_count":    int(adata.n_vars),
        "pixel_size_um": pixel_size_um,
        "software_version": "msseg-skill-1.0",
    }
    (output_dir / "experiment.xenium").write_text(
        json.dumps(meta, indent=2), encoding="utf-8"
    )

    print(f"[OK] Xenium export done: {output_dir}")
    return output_dir


# ── h5ad 匯出 ─────────────────────────────────────────────────────────────────

def export_h5ad(adata: ad.AnnData, output_dir: Path) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    out_path = output_dir / "msseg_final.h5ad"
    adata.write(out_path)
    print(f"[OK] h5ad saved: {out_path}  ({adata.n_obs} cells x {adata.n_vars} genes)")
    return out_path


# ── CLI ───────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(description="MSseg 匯出工具")
    parser.add_argument("--input",      required=True, help="輸入 h5ad 路徑")
    parser.add_argument("--masks-dir",  required=True, help="MCseg 遮罩目錄（含 *_mcseg.npy）")
    parser.add_argument("--output",     required=True, help="輸出根目錄")
    parser.add_argument("--format",     default="both",
                        choices=["xenium", "h5ad", "both"],
                        help="輸出格式（預設 both）")
    parser.add_argument("--pixel-size", type=float, default=0.2737,
                        help="µm/px（預設 0.2737 Visium HD）")
    args = parser.parse_args()

    adata      = ad.read_h5ad(args.input)
    masks_dir  = Path(args.masks_dir)
    output_dir = Path(args.output)

    if args.format in ("xenium", "both"):
        export_xenium(adata, masks_dir, output_dir / "xenium", args.pixel_size)

    if args.format in ("h5ad", "both"):
        export_h5ad(adata, output_dir)


if __name__ == "__main__":
    main()
