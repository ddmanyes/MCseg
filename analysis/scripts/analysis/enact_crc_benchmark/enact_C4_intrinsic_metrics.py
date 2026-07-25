"""
enact_C4_intrinsic_metrics.py
------------------------------
Intrinsic segmentation quality metrics for the ENACT CRC crop region.
Does NOT use StarDist GT cell-type labels — all metrics are label-free.

Metrics per method:
  1. n_cells        — number of detected cells
  2. FTC            — fraction of transcripts captured inside cell masks
  3. umi_density    — median UMI / cell area (µm²)  [area-normalised]
  4. NED            — Neighbour Expression Divergence (mean Hellinger dist, HVG=1000)
  5. doublet_rate   — mean co-expression rate of 4 mutually-exclusive marker pairs
  6. gt_coverage    — fraction of ENACT GT cells whose centroid falls inside a mask

Usage:
    cd /Volumes/SSD/plan_a
    uv run python submission_bioinformatics/scripts/analysis/enact_crc_benchmark/enact_C4_intrinsic_metrics.py
"""

import logging
import time
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.ndimage import grey_dilation

logging.basicConfig(level=logging.INFO, format="%(asctime)s  %(message)s",
                    datefmt="%H:%M:%S")
log = logging.getLogger(__name__)

# ─── Paths ────────────────────────────────────────────────────────────────────
PLAN_A      = Path(__file__).resolve().parents[4]
RESULT_DIR  = PLAN_A / "submission_bioinformatics" / "results" / "enact_method_comparison"
ENACT_F1    = PLAN_A / "submission_bioinformatics" / "results" / "enact_crc_f1"
TRANSCRIPT_CSV = RESULT_DIR / "proseg_out" / "transcripts_enact_crop.csv"
OUT_CSV    = RESULT_DIR / "intrinsic_metrics.csv"
OUT_FIG    = RESULT_DIR / "fig_intrinsic_metrics.png"

# ─── ENACT crop constants ─────────────────────────────────────────────────────
CROP_W, CROP_H = 10088, 13964
VHD_PIXEL_UM   = 0.2737
PIXEL_AREA_UM2 = VHD_PIXEL_UM ** 2

# ─── Method definitions ───────────────────────────────────────────────────────
METHODS = {
    "StarDist": {
        "mask": RESULT_DIR / "stardist_mask.npy",
        "gt_coverage": 1.000,
    },
    "MCseg": {
        "mask": ENACT_F1 / "mcseg_mask_7pass.npy",
        "gt_coverage": 0.674,
    },
    "SR": {
        "mask": RESULT_DIR / "sr_mask.npy",
        "gt_coverage": 0.938,
    },
    "NUC": {
        "mask": RESULT_DIR / "cellpose_nuc_mask.npy",
        "gt_coverage": 0.772,
    },
    "ProSeg": {
        "mask": RESULT_DIR / "proseg_mask.npy",
        "gt_coverage": 0.997,
    },
}

# Mutually exclusive marker pairs (C1 / doublet rate)
# From CRC transcript attribution study — manuscript Result 3
IMPOSSIBLE_PAIRS = [
    ("EPCAM",  "CD3E"),
    ("MUC2",   "NKG7"),
    ("ACTA2",  "CD3E"),
    ("PECAM1", "EPCAM"),
]

N_HVGS = 1000


# ─── Load transcripts ─────────────────────────────────────────────────────────

def load_transcripts() -> pd.DataFrame:
    log.info(f"Loading transcripts ({TRANSCRIPT_CSV.stat().st_size / 1e9:.1f} GB)...")
    t0 = time.time()
    df = pd.read_csv(str(TRANSCRIPT_CSV),
                     dtype={"x_loc": np.float32, "y_loc": np.float32, "gene": "category"})
    df["xi"] = df["x_loc"].round().astype(np.int32).clip(0, CROP_W - 1)
    df["yi"] = df["y_loc"].round().astype(np.int32).clip(0, CROP_H - 1)
    log.info(f"  {len(df):,} rows  {df['gene'].nunique():,} genes  {time.time()-t0:.0f}s")
    return df


# ─── Assign transcripts → cell IDs ───────────────────────────────────────────

def assign_transcripts(df: pd.DataFrame, mask: np.ndarray) -> np.ndarray:
    return mask[df["yi"].values, df["xi"].values]


# ─── Build expression matrix ──────────────────────────────────────────────────

def build_expr_matrix(df: pd.DataFrame, cell_ids: np.ndarray,
                      max_cells: int) -> tuple[np.ndarray, np.ndarray, list[str]]:
    gene_cats = df["gene"].cat.categories.tolist()
    n_genes   = len(gene_cats)

    assigned_mask = cell_ids > 0
    c_ids   = cell_ids[assigned_mask] - 1            # 0-indexed rows
    g_codes = df["gene"].cat.codes.values[assigned_mask]

    X_sp = sp.csr_matrix(
        (np.ones(c_ids.shape[0], dtype=np.float32), (c_ids, g_codes)),
        shape=(max_cells, n_genes), dtype=np.float32,
    )

    cell_umis  = np.asarray(X_sp.sum(axis=1)).ravel()
    valid      = cell_umis > 0
    cell_id_arr = np.where(valid)[0] + 1             # 1-indexed
    X_dense = np.asarray(X_sp[valid].todense(), dtype=np.float32)

    return X_dense, cell_id_arr, gene_cats


# ─── Metric functions ─────────────────────────────────────────────────────────

def compute_ftc(cell_ids: np.ndarray) -> float:
    return float((cell_ids > 0).mean())


def compute_umi_density(mask: np.ndarray, X: np.ndarray,
                        cell_id_arr: np.ndarray) -> float:
    pixel_counts = np.bincount(mask.ravel(), minlength=int(mask.max()) + 1)
    umi_per_cell = X.sum(axis=1)
    area_um2     = pixel_counts[cell_id_arr] * PIXEL_AREA_UM2
    valid        = area_um2 > 0
    density      = umi_per_cell[valid] / area_um2[valid]
    return float(np.median(density))


def compute_ned(mask: np.ndarray, X: np.ndarray,
                cell_id_arr: np.ndarray) -> float:
    if X.shape[0] < 5:
        return np.nan

    Xf = X
    if Xf.shape[1] > N_HVGS:
        hvg_idx = np.argpartition(Xf.var(axis=0), -N_HVGS)[-N_HVGS:]
        Xf = Xf[:, hvg_idx]

    row_sums = np.maximum(Xf.sum(axis=1, keepdims=True), 1e-10)
    X_prob   = (Xf / row_sums).astype(np.float32)
    cid_to_row = {int(c): i for i, c in enumerate(cell_id_arr)}

    struct  = np.ones((3, 3), dtype=np.int32)
    dilated = grey_dilation(mask.astype(np.int32), footprint=struct)
    bnd     = (mask > 0) & (dilated != mask)
    ci = mask[bnd].astype(np.int32)
    cj = dilated[bnd].astype(np.int32)
    ok = (cj > 0) & (cj != ci)
    ci, cj = ci[ok], cj[ok]

    if len(ci) == 0:
        return np.nan

    pairs = np.unique(np.sort(np.stack([ci, cj], axis=1), axis=1), axis=0)
    known = set(cid_to_row.keys())
    filt  = np.array([(int(a) in known) and (int(b) in known) for a, b in pairs])
    pairs = pairs[filt]

    if len(pairs) < 5:
        return np.nan

    if len(pairs) > 3000:
        pairs = pairs[np.random.default_rng(42).choice(len(pairs), 3000, replace=False)]

    i_idx = np.array([cid_to_row[int(a)] for a in pairs[:, 0]])
    j_idx = np.array([cid_to_row[int(b)] for b in pairs[:, 1]])

    sqrt_i = np.sqrt(np.maximum(X_prob[i_idx], 0))
    sqrt_j = np.sqrt(np.maximum(X_prob[j_idx], 0))
    hell   = np.sqrt(np.sum((sqrt_i - sqrt_j) ** 2, axis=1) / 2)

    return float(np.clip(np.mean(hell), 0, 1))


def compute_doublet_rate(X: np.ndarray, gene_list: list[str]) -> float:
    g2i   = {g: i for i, g in enumerate(gene_list)}
    rates = []
    for ga, gb in IMPOSSIBLE_PAIRS:
        if ga not in g2i or gb not in g2i:
            continue
        col_a = X[:, g2i[ga]]
        col_b = X[:, g2i[gb]]
        rates.append(float(((col_a > 0) & (col_b > 0)).mean()))
    return float(np.mean(rates)) if rates else np.nan


# ─── Main loop ────────────────────────────────────────────────────────────────

def run() -> pd.DataFrame:
    df_tx   = load_transcripts()
    records = []

    for name, info in METHODS.items():
        log.info(f"\n{'─'*50}")
        log.info(f"[{name}]")
        t0 = time.time()

        mask      = np.load(str(info["mask"]))
        max_cells = int(mask.max())

        cell_ids = assign_transcripts(df_tx, mask)
        ftc      = compute_ftc(cell_ids)
        log.info(f"  n_cells={max_cells:,}  FTC={ftc:.3f}")

        log.info("  Building expression matrix...")
        X, cell_id_arr, gene_list = build_expr_matrix(df_tx, cell_ids, max_cells)

        umi_dens = compute_umi_density(mask, X, cell_id_arr)
        log.info(f"  UMI density={umi_dens:.4f} UMI/µm²")

        log.info("  Computing NED...")
        ned = compute_ned(mask, X, cell_id_arr)
        log.info(f"  NED={ned:.3f}")

        dbl = compute_doublet_rate(X, gene_list)
        log.info(f"  Doublet rate={dbl*100:.2f}%")

        log.info(f"  GT coverage={info['gt_coverage']:.1%}  [{time.time()-t0:.0f}s]")

        records.append({
            "method":       name,
            "n_cells":      max_cells,
            "ftc":          round(ftc, 4),
            "umi_density":  round(umi_dens, 4),
            "ned":          round(ned, 4),
            "doublet_rate": round(dbl, 4),
            "gt_coverage":  round(info["gt_coverage"], 3),
        })
        del mask, X

    return pd.DataFrame(records)


# ─── Figure ───────────────────────────────────────────────────────────────────

def plot_results(df: pd.DataFrame) -> None:
    methods = df["method"].tolist()
    colors  = ["#4C8BE2", "#E05C5C", "#6DBF7E", "#F7A23B", "#A67FBF"]

    PANELS = [
        ("n_cells",      "# Detected cells",                    False),
        ("ftc",          "FTC (fraction transcripts in cells)", False),
        ("umi_density",  "UMI density (UMI / µm²)",            False),
        ("ned",          "NED (Hellinger distance)",            False),
        ("doublet_rate", "Doublet rate (C1)\nco-expression",   True),
        ("gt_coverage",  "GT cell coverage\n(ENACT-defined)",  False),
    ]

    fig, axes = plt.subplots(2, 3, figsize=(15, 9))

    for ax, (col, ylabel, lower_better) in zip(axes.flat, PANELS):
        vals     = df[col].values
        best_idx = int(np.argmin(vals) if lower_better else np.argmax(vals))
        bars     = ax.bar(methods, vals, color=colors, edgecolor="white", linewidth=0.5)
        bars[best_idx].set_edgecolor("black")
        bars[best_idx].set_linewidth(2.5)

        ax.set_title(ylabel, fontsize=10, fontweight="bold", pad=8)
        note = "lower = better" if lower_better else "higher = better"
        ax.set_ylabel(note, fontsize=8, color="gray")
        ax.tick_params(axis="x", labelsize=9, rotation=25)

        for bar, val in zip(bars, vals):
            if col == "n_cells":
                label = f"{int(val):,}"
            elif col in ("ftc", "doublet_rate", "gt_coverage"):
                label = f"{val*100:.1f}%"
            else:
                label = f"{val:.3f}"
            ax.text(bar.get_x() + bar.get_width() / 2,
                    bar.get_height() * 1.02,
                    label, ha="center", va="bottom", fontsize=8)

    fig.suptitle(
        "ENACT CRC Crop — Intrinsic Segmentation Quality\n"
        "(label-free; no StarDist GT cell-type labels used)",
        fontsize=12, fontweight="bold",
    )
    plt.tight_layout()
    fig.savefig(str(OUT_FIG), dpi=300, bbox_inches="tight")
    log.info(f"Figure saved: {OUT_FIG.name}")


# ─── Entry point ──────────────────────────────────────────────────────────────

def main() -> None:
    df = run()
    df.to_csv(str(OUT_CSV), index=False)
    log.info(f"\nSaved: {OUT_CSV.name}")
    log.info("\n" + df.to_string(index=False))
    plot_results(df)


if __name__ == "__main__":
    main()
