"""
27_fig3_revised.py
==================
Fig. 3 (fig4a_capture.png) — 6-panel combined, 2 rows × 3 cols:
  a) FTC | b) UMI/cell | c) Genes/cell
  d) UMI density | e) NED | f) Doublet rate (all boxplots, no panel titles)

Supp. Fig. S8 (fig4d_roi_heatmap.png):
  Per-ROI heatmap

Output:
  crc_transcript_attribution/results/figures/
  manuscript/figures/04_crc_tas/
"""

from __future__ import annotations

import json
import shutil
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.ticker as mticker
import yaml
from pathlib import Path
from scipy import stats

warnings.filterwarnings("ignore")

# ── 路徑 ──────────────────────────────────────────────────────────────────

ROOT       = Path(__file__).parent.parent
cfg        = yaml.safe_load((ROOT / "config.yaml").read_text())
PATHS      = cfg["paths"]
DATA       = cfg["data"]
PLT_CFG    = cfg["plotting"]

METRICS_DIR    = ROOT / PATHS["metrics_dir"]
FIGURES_DIR    = ROOT / PATHS["figures_dir"]
MANUSCRIPT_FIG_DIR  = Path("/Volumes/SSD/plan_a/manuscript/figures/04_crc_tas")
MANUSCRIPT_SUPP_DIR = Path("/Volumes/SSD/plan_a/manuscript/supplementary")

FIGURES_DIR.mkdir(parents=True, exist_ok=True)

METHODS       = DATA["methods"]
METHOD_COLORS = PLT_CFG["method_colors"]
DPI           = PLT_CFG["dpi"]

METHOD_LABELS = {
    "v12": "MCseg v2",
    "sr":  "SR",
    "p3":  "MCseg v1",
    "nuc": "NUC",
}

METHOD_ORDER = [m for m in ["sr", "v12", "p3", "nuc"] if m in METHODS]

MM_TO_IN   = 1 / 25.4
DOUBLE_COL = 183 * MM_TO_IN

with open(PATHS["roi_info"]) as f:
    ROI_INFO = json.load(f)

# ── 全域樣式 ───────────────────────────────────────────────────────────────

plt.rcParams.update({
    "font.family":        ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size":          7,
    "axes.labelsize":     8,
    "axes.titlesize":     8,
    "xtick.labelsize":    7,
    "ytick.labelsize":    7,
    "legend.fontsize":    7,
    "figure.dpi":         150,
    "savefig.dpi":        300,
    "savefig.bbox":       "tight",
    "savefig.facecolor":  "white",
    "axes.spines.top":    False,
    "axes.spines.right":  False,
})


def clean_ax(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


_SUPP_NAMES = {"fig4d_roi_heatmap.png", "SuppFigS8.png"}

def save_fig(fig, name: str):
    dst = FIGURES_DIR / name
    fig.savefig(dst, dpi=DPI, bbox_inches="tight", facecolor="white")
    if name in _SUPP_NAMES:
        # 補充圖 → supplementary/
        out_name = "SuppFigS8.png"
        if MANUSCRIPT_SUPP_DIR.exists():
            shutil.copy2(dst, MANUSCRIPT_SUPP_DIR / out_name)
            print(f"  ✓ {name}  →  manuscript/supplementary/{out_name}")
        else:
            print(f"  ✓ {name}  (supp dir not found, skipped copy)")
    else:
        # 主圖 → figures/04_crc_tas/
        if MANUSCRIPT_FIG_DIR.exists():
            shutil.copy2(dst, MANUSCRIPT_FIG_DIR / name)
            print(f"  ✓ {name}  →  manuscript/figures/04_crc_tas/")
        else:
            print(f"  ✓ {name}  (fig dir not found, skipped copy)")
    plt.close(fig)


def _p_stars(p: float) -> str:
    if p < 0.001: return "***"
    if p < 0.01:  return "**"
    if p < 0.05:  return "*"
    return "ns"


def _add_box_bracket(ax, x1: int, x2: int, y: float, p: float,
                     val_range: float, fontsize: float = 6.5):
    """兩箱型圖之間畫 Wilcoxon 括號。"""
    tick_h = val_range * 0.025
    ax.plot([x1, x1, x2, x2], [y - tick_h, y, y, y - tick_h],
            color="#333", lw=0.75, solid_capstyle="round")
    ax.text((x1 + x2) / 2, y + val_range * 0.015, _p_stars(p),
            ha="center", va="bottom", fontsize=fontsize, color="#333")


def _boxplot_panel(ax, data_list, colors, seed=0, s=12, alpha=0.75):
    """畫箱型圖 + jitter，回傳 bp 物件。"""
    n = len(data_list)
    bp = ax.boxplot(
        data_list, positions=range(n), widths=0.5,
        patch_artist=True,
        medianprops=dict(color="black", linewidth=1.8),
        whiskerprops=dict(linewidth=0.75),
        capprops=dict(linewidth=0.75),
        flierprops=dict(marker="", markersize=0),
        zorder=3,
    )
    for patch, color in zip(bp["boxes"], colors):
        patch.set_facecolor(color); patch.set_alpha(0.60); patch.set_linewidth(0.75)
    for w, c in zip(bp["whiskers"], [c for c in colors for _ in range(2)]):
        w.set_color(c)
    for cap, c in zip(bp["caps"], [c for c in colors for _ in range(2)]):
        cap.set_color(c)
    rng = np.random.default_rng(seed=seed)
    for i, (vals, color) in enumerate(zip(data_list, colors)):
        jit = rng.uniform(-0.13, 0.13, len(vals))
        ax.scatter(np.full(len(vals), i) + jit, vals,
                   color=color, s=s, alpha=alpha, zorder=6,
                   linewidths=0.4, edgecolors="white")
    return bp


def _whisker_tops(data_list):
    """各資料集的 whisker 上緣（用於括號定位）。"""
    tops = []
    for vals in data_list:
        if not vals: tops.append(0.0); continue
        q3  = float(np.percentile(vals, 75))
        iqr = float(np.percentile(vals, 75) - np.percentile(vals, 25))
        tops.append(min(q3 + 1.5 * iqr, max(vals)))
    return tops


def _wilcoxon_p(vals_i, vals_j):
    """配對 Wilcoxon（跳過 NaN），回傳 p 值。"""
    paired = [(a, b) for a, b in zip(vals_i, vals_j)
              if not (np.isnan(float(a)) or np.isnan(float(b)))]
    if len(paired) >= 3:
        try:
            _, p = stats.wilcoxon(*zip(*paired))
            return p
        except Exception:
            pass
    return 1.0


# ═══════════════════════════════════════════════════════════════════════════
# Fig. 3  (fig4a_capture.png)
# 2 rows × 3 cols:
#   Row 1:  a) FTC  b) UMI/cell  c) Genes/cell
#   Row 2:  d) UMI density  e) NED  f) Doublet rate
# ═══════════════════════════════════════════════════════════════════════════

def fig3_combined(df_ac: pd.DataFrame):
    mord   = METHOD_ORDER
    labels = [METHOD_LABELS.get(m, m) for m in mord]
    colors = [METHOD_COLORS.get(m, "#888") for m in mord]
    n_meth = len(mord)
    rois   = list(ROI_INFO.keys())

    idx_sr  = mord.index("sr")  if "sr"  in mord else None
    idx_v12 = mord.index("v12") if "v12" in mord else None
    idx_p3  = mord.index("p3")  if "p3"  in mord else None

    fig = plt.figure(figsize=(DOUBLE_COL * 1.30, 170 * MM_TO_IN))
    gs  = fig.add_gridspec(2, 3, wspace=0.42, hspace=0.38)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[0, 2])
    ax_d = fig.add_subplot(gs[1, 0])
    ax_e = fig.add_subplot(gs[1, 1])
    ax_f = fig.add_subplot(gs[1, 2])

    # ── a: FTC ────────────────────────────────────────────────────────────
    bp_a_data = [df_ac[df_ac["method"] == m]["a1_capture"].dropna().values.tolist()
                 for m in mord]
    _boxplot_panel(ax_a, bp_a_data, colors, seed=0, s=12, alpha=0.75)

    all_ftc   = [v for d in bp_a_data for v in d]
    ftc_range = max(all_ftc) - min(all_ftc) if all_ftc else 0.1
    wt_a      = _whisker_tops(bp_a_data)
    b_base_a  = max(wt_a) + ftc_range * 0.12
    gap_a     = ftc_range * 0.14

    for k, (ii, jj) in enumerate([(idx_sr, idx_v12), (idx_sr, idx_p3)]):
        if ii is None or jj is None: continue
        p = _wilcoxon_p(bp_a_data[ii], bp_a_data[jj])
        _add_box_bracket(ax_a, ii, jj, b_base_a + k * gap_a, p, ftc_range)

    ax_a.set_xticks(range(n_meth)); ax_a.set_xticklabels([""] * n_meth)
    ax_a.set_ylabel("Transcript capture rate (FTC)", fontsize=8)
    ax_a.set_ylim(0, b_base_a + 2 * gap_a + ftc_range * 0.15)
    ax_a.yaxis.set_major_locator(mticker.MultipleLocator(0.2))
    ax_a.axhline(1.0, color="#bbb", lw=0.5, ls="--", zorder=1)
    ax_a.annotate("↑ better", xy=(0.97, 0.05), xycoords="axes fraction",
                  ha="right", va="bottom", fontsize=6.5, color="#555")
    clean_ax(ax_a)
    ax_a.text(-0.15, 1.06, "b", transform=ax_a.transAxes,
              fontsize=9, fontweight="bold", va="top", ha="right")

    # ── b: UMI/cell ───────────────────────────────────────────────────────
    bp_b_data = [df_ac[df_ac["method"] == m]["a2_median_umi"].dropna().values.tolist()
                 for m in mord]
    _boxplot_panel(ax_b, bp_b_data, colors, seed=10, s=14, alpha=0.8)

    all_umi   = [v for d in bp_b_data for v in d]
    umi_range = max(all_umi) - min(all_umi) if all_umi else 200
    wt_b      = _whisker_tops(bp_b_data)
    b_base_b  = max(wt_b) + umi_range * 0.10
    gap_b     = umi_range * 0.16

    for k, (ii, jj) in enumerate([
        (idx_sr, idx_v12),
        (idx_v12, idx_p3),
        (idx_sr, idx_p3),
    ]):
        if ii is None or jj is None: continue
        p = _wilcoxon_p(bp_b_data[ii], bp_b_data[jj])
        _add_box_bracket(ax_b, ii, jj, b_base_b + k * gap_b, p, umi_range)

    ax_b.axhline(50, color="#e74c3c", lw=0.9, ls="--", zorder=2)
    ax_b.set_xticks(range(n_meth)); ax_b.set_xticklabels([""] * n_meth)
    ax_b.set_ylabel("Median UMI per cell", fontsize=8)
    ax_b.set_ylim(0, b_base_b + 3 * gap_b + umi_range * 0.10)
    ax_b.yaxis.set_major_locator(mticker.MultipleLocator(200))
    clean_ax(ax_b)
    ax_b.text(-0.15, 1.06, "c", transform=ax_b.transAxes,
              fontsize=9, fontweight="bold", va="top", ha="right")

    # ── c: Genes/cell ─────────────────────────────────────────────────────
    bp_c_data = [df_ac[df_ac["method"] == m]["a3_median_genes"].dropna().values.tolist()
                 for m in mord]
    _boxplot_panel(ax_c, bp_c_data, colors, seed=77, s=14, alpha=0.8)

    all_genes   = [v for d in bp_c_data for v in d]
    gene_range  = max(all_genes) - min(all_genes) if all_genes else 200
    wt_c        = _whisker_tops(bp_c_data)
    b_base_c    = max(wt_c) + gene_range * 0.10
    gap_c       = gene_range * 0.16

    for k, (ii, jj) in enumerate([(idx_sr, idx_p3), (idx_v12, idx_p3)]):
        if ii is None or jj is None: continue
        p = _wilcoxon_p(bp_c_data[ii], bp_c_data[jj])
        _add_box_bracket(ax_c, ii, jj, b_base_c + k * gap_c, p, gene_range)

    gene_ref = 300
    ax_c.axhline(gene_ref, color="#e74c3c", lw=0.9, ls="--", zorder=2)
    ax_c.set_xticks(range(n_meth)); ax_c.set_xticklabels([""] * n_meth)
    ax_c.set_ylabel("Median genes per cell", fontsize=8)
    ax_c.set_ylim(0, b_base_c + 2 * gap_c + gene_range * 0.15)
    ax_c.yaxis.set_major_locator(mticker.MultipleLocator(200))
    clean_ax(ax_c)
    ax_c.text(-0.15, 1.06, "d", transform=ax_c.transAxes,
              fontsize=9, fontweight="bold", va="top", ha="right")

    # ── d: UMI density ────────────────────────────────────────────────────
    bp_d_data = [df_ac[df_ac["method"] == m]["a1_umi_density"].dropna().values.tolist()
                 for m in mord]
    _boxplot_panel(ax_d, bp_d_data, colors, seed=99, s=14, alpha=0.8)

    all_dens   = [v for d in bp_d_data for v in d]
    dens_range = max(all_dens) - min(all_dens) if all_dens else 5
    wt_d       = _whisker_tops(bp_d_data)

    if idx_sr is not None and idx_v12 is not None:
        p = _wilcoxon_p(bp_d_data[idx_sr], bp_d_data[idx_v12])
        _add_box_bracket(ax_d, idx_sr, idx_v12,
                         max(wt_d) + dens_range * 0.18, p, dens_range)

    ax_d.set_xticks(range(n_meth)); ax_d.set_xticklabels([""] * n_meth)
    ax_d.set_ylabel("UMI density (UMI/µm²)", fontsize=8)
    ax_d.set_ylim(0, max(all_dens) * 1.55 if all_dens else 20)
    ax_d.annotate("FTC Paradox →", xy=(0.97, 0.95), xycoords="axes fraction",
                  ha="right", va="top", fontsize=6, color="#e74c3c", style="italic")
    clean_ax(ax_d)
    ax_d.text(-0.15, 1.06, "e", transform=ax_d.transAxes,
              fontsize=9, fontweight="bold", va="top", ha="right")

    # ── e: NED ────────────────────────────────────────────────────────────
    method_ned = {}
    for method in mord:
        sub = df_ac[df_ac["method"] == method].set_index("roi")
        method_ned[method] = [
            sub.loc[r, "ned"] if r in sub.index else np.nan for r in rois
        ]
    bp_e_data = [[v for v in method_ned[m] if not np.isnan(v)] for m in mord]
    _boxplot_panel(ax_e, bp_e_data, colors, seed=42, s=14, alpha=0.85)

    mean_ned  = {m: float(np.nanmean(method_ned[m])) for m in mord}
    wt_e      = _whisker_tops(bp_e_data)
    all_ned   = [v for d in bp_e_data for v in d]
    ned_range = max(all_ned) - min(all_ned) if len(all_ned) > 1 else 0.1

    b_base_e = max(wt_e) + ned_range * 0.16
    gap_e    = ned_range * 0.20

    for k, (ii, jj) in enumerate([(0, 1), (1, 2), (2, 3)]):
        if jj >= n_meth: continue
        p = _wilcoxon_p(method_ned[mord[ii]], method_ned[mord[jj]])
        h      = b_base_e + k * gap_e
        tick_h = ned_range * 0.025
        ax_e.plot([ii, ii, jj, jj], [h - tick_h, h, h, h - tick_h],
                  color="#333", lw=0.75, solid_capstyle="round")
        ax_e.text((ii + jj) / 2, h + ned_range * 0.015, _p_stars(p),
                  ha="center", va="bottom", fontsize=6.5, color="#333")

    y_ceil_e  = b_base_e + 3 * gap_e + ned_range * 0.08
    y_floor_e = min(all_ned) - ned_range * 0.05 if all_ned else 0.0
    ax_e.set_xlim(-0.65, n_meth - 0.35)
    ax_e.set_ylim(y_floor_e, y_ceil_e)
    ax_e.set_xticks(range(n_meth)); ax_e.set_xticklabels([""] * n_meth)
    ax_e.set_ylabel("NED (Hellinger distance)", fontsize=8)
    ax_e.annotate("↑ better", xy=(0.97, 0.95), xycoords="axes fraction",
                  ha="right", va="top", fontsize=6.5, color="#555")
    clean_ax(ax_e)
    ax_e.text(-0.15, 1.06, "f", transform=ax_e.transAxes,
              fontsize=9, fontweight="bold", va="top", ha="right")

    # ── f: Doublet rate (boxplot) ─────────────────────────────────────────
    bp_f_data = [(df_ac[df_ac["method"] == m]["c1_coexpr"].dropna() * 100).values.tolist()
                 for m in mord]
    _boxplot_panel(ax_f, bp_f_data, colors, seed=55, s=14, alpha=0.8)

    all_dr   = [v for d in bp_f_data for v in d]
    dr_range = max(all_dr) - min(all_dr) if all_dr else 1.0
    wt_f     = _whisker_tops(bp_f_data)

    if idx_sr is not None and idx_v12 is not None:
        p = _wilcoxon_p(bp_f_data[idx_sr], bp_f_data[idx_v12])
        _add_box_bracket(ax_f, idx_sr, idx_v12,
                         max(wt_f) + dr_range * 0.18, p, dr_range)

    for i, method in enumerate(mord):
        mean_f = np.mean(bp_f_data[i]) if bp_f_data[i] else np.nan
        if not np.isnan(mean_f):
            ax_f.text(i, max(all_dr) * 1.25, f"{mean_f:.2f}%",
                      ha="center", va="top", fontsize=6,
                      color="#444", fontweight="bold",
                      bbox=dict(boxstyle="round,pad=0.1", fc="white", ec="none", alpha=0.85))

    ax_f.set_xticks(range(n_meth)); ax_f.set_xticklabels([""] * n_meth)
    ax_f.set_ylabel("Doublet rate (%)", fontsize=8)
    ax_f.set_ylim(0, max(all_dr) * 1.65 if all_dr else 5)
    ax_f.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{v:.2f}%"))
    ax_f.annotate("↓ better", xy=(0.97, 0.95), xycoords="axes fraction",
                  ha="right", va="top", fontsize=6.5, color="#555")
    clean_ax(ax_f)
    ax_f.text(-0.15, 1.06, "g", transform=ax_f.transAxes,
              fontsize=9, fontweight="bold", va="top", ha="right")

    # 共用圖例
    handles = [mpatches.Patch(color=c, label=l, alpha=0.88)
               for c, l in zip(colors, labels)]
    fig.legend(handles=handles, loc="lower center",
               ncol=n_meth, bbox_to_anchor=(0.5, -0.03),
               frameon=False, fontsize=7.5)

    save_fig(fig, "fig4a_capture.png")


# ═══════════════════════════════════════════════════════════════════════════
# Supp Fig S8  (fig4d_roi_heatmap.png)
# FTC / UMI Density / NED / Doublet Rate  ×  ROIs × methods
# ═══════════════════════════════════════════════════════════════════════════

def suppfig_s8_heatmap(df_ac: pd.DataFrame):
    """
    Supplementary Fig S8 — Per-ROI metric heatmap
    Rows: method × ROI
    Cols: FTC | UMI Density | NED | Doublet Rate
    """
    mord = METHOD_ORDER
    rois = list(ROI_INFO.keys())

    cols_raw = {
        "a1_capture":     "FTC",
        "a1_umi_density": "UMI Density\n(UMI/µm²)",
        "ned":            "NED",
        "c1_coexpr":      "Doublet Rate",
    }

    records = []
    for method in mord:
        sub = df_ac[df_ac["method"] == method].copy()
        for roi in rois:
            row = sub[sub["roi"] == roi]
            if row.empty:
                continue
            rec = {
                "Method": METHOD_LABELS.get(method, method),
                "ROI": roi.upper(),
            }
            for col in cols_raw:
                rec[col] = float(row.iloc[0][col]) if col in row.columns else np.nan
            records.append(rec)

    df_heat = pd.DataFrame(records)
    df_heat["row_label"] = (
        df_heat["Method"].str.replace("(Ours)", "", regex=False).str.strip()
        + " – " + df_heat["ROI"]
    )
    df_heat = df_heat.set_index("row_label")

    mat_raw = df_heat[list(cols_raw.keys())].astype(float)
    mat_raw.columns = list(cols_raw.values())

    mat_norm = mat_raw.copy()
    for col in mat_norm.columns:
        col_min, col_max = mat_norm[col].min(), mat_norm[col].max()
        if col_max > col_min:
            mat_norm[col] = (mat_norm[col] - col_min) / (col_max - col_min)
        if "Doublet" in col:
            mat_norm[col] = 1 - mat_norm[col]

    n_row = len(df_heat)
    n_col = len(cols_raw)
    fig_h = max(140, n_row * 3.2) * MM_TO_IN
    fig_w = (50 + n_col * 22) * MM_TO_IN
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    im = ax.imshow(mat_norm.values, aspect="auto", cmap="RdYlGn",
                   vmin=0.0, vmax=1.0)

    ax.set_xticks(range(n_col))
    ax.set_xticklabels(mat_raw.columns, fontsize=7.5, rotation=30, ha="right")
    ax.set_yticks(range(n_row))
    ax.set_yticklabels(mat_raw.index, fontsize=5.8)
    ax.tick_params(left=False, bottom=False)

    for k in range(1, len(mord)):
        ax.axhline(k * len(rois) - 0.5, color="white", lw=2.0)

    for i in range(n_row):
        for j, col in enumerate(mat_raw.columns):
            val = mat_raw.values[i, j]
            if np.isnan(val):
                continue
            if "UMI Density" in col:
                txt = f"{val:.1f}"
            elif "Doublet" in col:
                txt = f"{val*100:.2f}%"
            else:
                txt = f"{val:.3f}"
            ax.text(j, i, txt, ha="center", va="center",
                    fontsize=5.0, color="black")

    cbar = plt.colorbar(im, ax=ax, fraction=0.020, pad=0.02)
    cbar.ax.tick_params(labelsize=6)
    cbar.set_label("Relative performance\n(per column min–max)", fontsize=6.5)

    ax.set_title(
        f"Supp. Fig. S8 | Per-ROI Transcript Attribution Metrics "
        f"({len(mord)} Methods × {len(rois)} ROIs)",
        fontsize=8, pad=6
    )
    for spine in ax.spines.values():
        spine.set_visible(False)

    fig.text(0.01, -0.01,
             "Colour: column-wise min–max normalisation; Doublet Rate colour inverted (lower = better = green).",
             fontsize=5.5, ha="left", va="top", color="#555")

    save_fig(fig, "fig4d_roi_heatmap.png")


if __name__ == "__main__":
    print("[27_fig3_revised] 讀取 metrics_ac.csv ...")
    df_ac = pd.read_csv(METRICS_DIR / "metrics_ac.csv")

    print("\n[Fig 3] Combined 6-panel (a–f) ...")
    fig3_combined(df_ac)

    print("\n[Supp Fig S8] Per-ROI Heatmap ...")
    suppfig_s8_heatmap(df_ac)

    print("\n✅ 完成。輸出至：")
    print(f"   {FIGURES_DIR}")
    print(f"   主圖：{MANUSCRIPT_FIG_DIR}")
    print(f"   補充：{MANUSCRIPT_SUPP_DIR}")
