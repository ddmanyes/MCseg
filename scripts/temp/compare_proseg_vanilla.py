"""
比較 proseg 兩種執行方式在 tile_y0_x0 上的結果：
  A) 我們的管線：Python Watershed 注入 cell_id + nuclear-reassignment-prob=0
  B) 原版 proseg：--cellpose-masks + expand-initialized-cells + 預設 reassignment-prob=0.2

輸出：figures/proseg_vanilla_comparison.png（2×2 grid）
      figures/proseg_vanilla_comparison_zoom.png（放大中心區塊）
"""
import gzip
import json
import os
import subprocess
import sys
import time
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import tifffile

# ── 路徑設定 ──────────────────────────────────────────────────────────────
ROOT       = Path("/Volumes/SSD/plan_a/visiumHD_pipeline_2")
TILE_DIR   = ROOT / "results/analysis/roi/test/proseg_tiles/tile_y0_x0"
MASK_NPY   = ROOT / "results/analysis/roi/test/segmentation_masks.npy"
HE_TIF     = ROOT / "results/analysis/roi/test/he_crop.tif"
OUT_DIR    = ROOT / "results/analysis/roi/test/proseg_vanilla"
FIG_DIR    = ROOT / "figures"
PROSEG_BIN = Path("~/.cargo/bin/proseg").expanduser()

# proseg 參數
COORD_SCALE   = 0.2737   # µm/px（ROI pixel size，和管線一致）
MAX_DIST      = 40.0
COMPACTNESS   = 0.06
SAMPLES       = 200
BURNIN        = 150
RECORDED      = 50
EXPAND_PX     = 20       # 取代 Python dilation=20px

# tile_y0_x0 pixel 邊界（grid 4×3, padding=200）
TILE_X0, TILE_X1 = 0, 572
TILE_Y0, TILE_Y1 = 0, 603

# ── 建立輸出目錄 ───────────────────────────────────────────────────────────
OUT_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)


# ── Step 1: 裁切 Cellpose mask 至 tile 範圍 ────────────────────────────────
print("Step 1: 裁切 Cellpose mask …")
full_mask = np.load(str(MASK_NPY))
tile_mask = full_mask[TILE_Y0:TILE_Y1, TILE_X0:TILE_X1].copy()

# 重新編號：讓 cell id 從 1 連續（proseg 需要）
uniq = np.unique(tile_mask)
uniq = uniq[uniq > 0]
remap = np.zeros(int(full_mask.max()) + 1, dtype=np.int32)
for new_id, old_id in enumerate(uniq, start=1):
    remap[old_id] = new_id
tile_mask_renum = remap[tile_mask]

mask_path = OUT_DIR / "tile_mask.npy"
np.save(str(mask_path), tile_mask_renum.astype(np.uint32))
print(f"  mask shape: {tile_mask_renum.shape}, cells: {len(uniq)}")


# ── Step 2: 建立「無 cell_id」transcript CSV ───────────────────────────────
print("Step 2: 準備 vanilla transcript CSV …")
df = pd.read_csv(TILE_DIR / "transcripts_for_proseg.csv")
df_vanilla = df.copy()
df_vanilla["cell_id"] = 0   # 全部設為 unassigned

vanilla_csv = OUT_DIR / "transcripts_vanilla.csv"
df_vanilla.to_csv(vanilla_csv, index=False)
print(f"  transcripts: {len(df_vanilla):,} rows")


# ── Step 3: 執行 proseg（原版參數） ────────────────────────────────────────
poly_out  = OUT_DIR / "polygons_vanilla.json"
count_out = OUT_DIR / "counts_vanilla.csv.gz"
meta_out  = OUT_DIR / "cells_vanilla.csv"
gene_out  = OUT_DIR / "genes_vanilla.csv"

cmd = [
    str(PROSEG_BIN),
    "--overwrite",
    "--cellpose-masks", str(mask_path),
    "--cellpose-scale", str(COORD_SCALE),
    "--expand-initialized-cells", str(EXPAND_PX),
    "--output-cell-polygons", str(poly_out),
    "--output-counts", str(count_out),
    "--output-counts-fmt", "csv-gz",
    "--output-cell-metadata", str(meta_out),
    "--output-cell-metadata-fmt", "csv",
    "--output-gene-metadata", str(gene_out),
    "--output-gene-metadata-fmt", "csv",
    "--coordinate-scale", str(COORD_SCALE),
    "--gene-column", "gene",
    "--x-column", "x",
    "--y-column", "y",
    "--z-column", "z",
    "--cell-id-column", "cell_id",
    "--cell-id-unassigned", "0",
    "--ignore-z-coord",
    "--min-qv", "0",
    "--max-transcript-nucleus-distance", str(MAX_DIST),
    "--cell-compactness", str(COMPACTNESS),
    "--samples", str(SAMPLES),
    "--burnin-samples", str(BURNIN),
    "--recorded-samples", str(RECORDED),
    "--enforce-connectivity",
    # 原版預設：不加 --nuclear-reassignment-prob 0（保留 0.2）
    # diffusion 模型也保留預設（不加 --no-diffusion）
    str(vanilla_csv),
]

print("Step 3: 執行 proseg（原版參數）…")
print("  指令：", " ".join(cmd))
t0 = time.time()
result = subprocess.run(cmd, capture_output=True, text=True)
elapsed = time.time() - t0
print(f"  耗時：{elapsed:.0f}s  returncode={result.returncode}")
if result.returncode != 0:
    print("STDERR:", result.stderr[-3000:])
    sys.exit(1)
print("  proseg 完成")


# ── Step 4: 讀取兩份 GeoJSON ──────────────────────────────────────────────
def load_geojson(path: Path):
    path = Path(path)
    magic = path.read_bytes()[:2]
    if magic == b'\x1f\x8b':
        with gzip.open(path, 'rt') as f:
            return json.load(f)
    with open(path) as f:
        return json.load(f)

print("Step 4: 讀取 GeoJSON …")
geo_ours    = load_geojson(TILE_DIR / "proseg_results.json")
geo_vanilla = load_geojson(poly_out)
print(f"  我們的管線: {len(geo_ours['features'])} cells")
print(f"  原版 proseg: {len(geo_vanilla['features'])} cells")


# ── Step 5: 讀取 H&E ───────────────────────────────────────────────────────
print("Step 5: 讀取 H&E …")
he = tifffile.imread(str(HE_TIF))
if he.ndim == 3 and he.shape[0] in (3, 4):
    he = np.moveaxis(he, 0, -1)
he_rgb = he[:, :, :3] if he.ndim == 3 and he.shape[-1] >= 3 else he

# 裁切到 tile 範圍
he_tile = he_rgb[TILE_Y0:TILE_Y1, TILE_X0:TILE_X1]
print(f"  H&E tile: {he_tile.shape}")


# ── Step 6: 轉換 GeoJSON 座標為 H&E 像素 ────────────────────────────────────
UM_PX = COORD_SCALE  # µm → px 換算（proseg 輸出 µm，÷ COORD_SCALE 得 px）

def geojson_to_px_polys(geo):
    """從 GeoJSON features 取出所有 polygon 外環，轉換為 H&E px 座標 list。"""
    polys = []
    for feat in geo["features"]:
        geom = feat["geometry"]
        if geom["type"] == "Polygon":
            rings = [geom["coordinates"][0]]
        elif geom["type"] == "MultiPolygon":
            rings = [p[0] for p in geom["coordinates"]]
        else:
            continue
        for ring in rings:
            arr = np.array(ring)
            xs = arr[:, 0] / UM_PX
            ys = arr[:, 1] / UM_PX
            polys.append((xs, ys))
    return polys


polys_ours    = geojson_to_px_polys(geo_ours)
polys_vanilla = geojson_to_px_polys(geo_vanilla)


# ── Step 7: 繪圖（全圖 + zoom） ───────────────────────────────────────────
# zoom 區域：中央 200×200 px
CX = (TILE_X1 - TILE_X0) // 2
CY = (TILE_Y1 - TILE_Y0) // 2
ZW = ZH = 200
zx0 = max(0, CX - ZW // 2); zx1 = zx0 + ZW
zy0 = max(0, CY - ZH // 2); zy1 = zy0 + ZH

def draw_polys(ax, polys, color, lw=0.8, alpha=0.85):
    for xs, ys in polys:
        ax.plot(xs, ys, color=color, lw=lw, alpha=alpha)

labels_info = [
    ("我們的管線\n(Watershed + reassign-prob=0)",
     polys_ours,    len(geo_ours["features"]),    "cyan"),
    ("原版 proseg\n(cellpose-masks + reassign-prob=0.2)",
     polys_vanilla, len(geo_vanilla["features"]), "lime"),
]

# --- 圖 1：全 tile 比較（1×2）---
fig1, axes1 = plt.subplots(1, 2, figsize=(18, 8), dpi=150)
fig1.patch.set_facecolor("#1a1a1a")

for ax, (label, polys, n_cells, color) in zip(axes1, labels_info):
    ax.imshow(he_tile, aspect="equal")
    draw_polys(ax, polys, color=color, lw=0.7)
    ax.set_title(f"{label}\n({n_cells} cells)", fontsize=10, color="white", pad=6)
    ax.axis("off")
    # zoom 矩形
    rect = plt.Rectangle((zx0, zy0), ZW, ZH,
                          linewidth=1.5, edgecolor="yellow", facecolor="none")
    ax.add_patch(rect)
    ax.set_facecolor("#1a1a1a")

plt.tight_layout(pad=1.5)
out1 = FIG_DIR / "proseg_vanilla_comparison.png"
plt.savefig(str(out1), dpi=150, bbox_inches="tight", facecolor=fig1.get_facecolor())
print(f"儲存：{out1}")

# --- 圖 2：zoom 比較（1×2）---
he_zoom = he_tile[zy0:zy1, zx0:zx1]
fig2, axes2 = plt.subplots(1, 2, figsize=(16, 8), dpi=200)
fig2.patch.set_facecolor("#1a1a1a")

for ax, (label, polys, n_cells, color) in zip(axes2, labels_info):
    ax.imshow(he_zoom, aspect="equal", extent=[zx0, zx1, zy1, zy0])
    # 只畫 zoom 範圍內的 poly（用 centroid 過濾）
    for xs, ys in polys:
        cx, cy = xs.mean(), ys.mean()
        if zx0 <= cx <= zx1 and zy0 <= cy <= zy1:
            ax.plot(xs, ys, color="white", lw=1.5, alpha=0.9)
            ax.plot(xs, ys, color=color,   lw=0.8, alpha=0.8)
    ax.set_xlim(zx0, zx1); ax.set_ylim(zy1, zy0)
    ax.set_title(f"Zoom — {label}\n({n_cells} cells)", fontsize=9, color="white", pad=5)
    ax.axis("off")
    ax.set_facecolor("#1a1a1a")

plt.tight_layout(pad=1.5)
out2 = FIG_DIR / "proseg_vanilla_comparison_zoom.png"
plt.savefig(str(out2), dpi=200, bbox_inches="tight", facecolor=fig2.get_facecolor())
print(f"儲存：{out2}")

# ── Step 8: 統計摘要 ──────────────────────────────────────────────────────
print()
print("=" * 60)
print("比較摘要")
print("=" * 60)
for label, geo, _ in [
    ("我們的管線", geo_ours, None),
    ("原版 proseg", geo_vanilla, None),
]:
    n = len(geo["features"])
    # 面積（um²）→ 換算面積
    areas = []
    for feat in geo["features"]:
        coords = feat["geometry"]["coordinates"]
        if feat["geometry"]["type"] == "Polygon":
            ring = np.array(coords[0])
        else:
            ring = np.array(coords[0][0])
        # Shoelace formula
        x_, y_ = ring[:, 0], ring[:, 1]
        area = 0.5 * abs(np.dot(x_, np.roll(y_, 1)) - np.dot(y_, np.roll(x_, 1)))
        areas.append(area)
    areas = np.array(areas)
    print(f"{label}: {n} cells  "
          f"面積 median={np.median(areas):.1f} µm²  "
          f"mean={np.mean(areas):.1f}  "
          f"p5={np.percentile(areas,5):.1f}  p95={np.percentile(areas,95):.1f}")

print()
print(f"完成！圖表已儲存至 {FIG_DIR}/proseg_vanilla_comparison*.png")
