"""
fig3a_transcript_metrics.py
============================
Fig. 3d–i — 6-panel boxplot (2 rows × 3 cols):
  d) FTC   e) UMI/cell   f) Genes/cell
  g) UMI density   h) NED   i) Doublet rate

Methods (6): SR, MCseg, 2Cseg, NUC, Proseg, StarDist+WBA

Output:
  submission_bioinformatics/figures/fig3/fig3d.png
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

_CRC_ROOT  = Path("/Volumes/SSD/plan_a/crc_transcript_attribution")
cfg        = yaml.safe_load((_CRC_ROOT / "config.yaml").read_text())
PATHS      = cfg["paths"]
PLT_CFG    = cfg["plotting"]

METRICS_DIR    = _CRC_ROOT / PATHS["metrics_dir"]
FIGURES_DIR    = Path(__file__).parents[2] / "figures" / "fig3"
MANUSCRIPT_FIG_DIR  = Path("/Volumes/SSD/plan_a/manuscript/figures/04_crc_tas")

FIGURES_DIR.mkdir(parents=True, exist_ok=True)

DPI = PLT_CFG["dpi"]

# 6-method order and display labels
METHOD_ORDER = ["sr", "v12", "p3", "nuc", "proseg", "stardist_wba"]
METHOD_LABELS = {
    "v12":          "MCseg",
    "sr":           "SR",
    "p3":           "2Cseg",
    "nuc":          "NUC",
    "proseg":       "Proseg",
    "stardist_wba": "SD+WBA",
}
METHOD_COLORS = {
    "v12":          "#2196F3",   # blue
    "sr":           "#FF5722",   # orange-red
    "p3":           "#4CAF50",   # green
    "nuc":          "#FF9800",   # amber
    "proseg":       "#2CA02C",   # darker green (supplementary)
    "stardist_wba": "#17BECF",   # cyan (ENACT)
}

MM_TO_IN   = 1 / 25.4
DOUBLE_COL = 183 * MM_TO_IN

with open(PATHS["roi_info"]) as f:
    ROI_INFO = json.load(f)

# ── 全域樣式 ───────────────────────────────────────────────────────────────

plt.rcParams.update({
    "font.family":        ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size":          11,
    "axes.labelsize":     13,
    "axes.titlesize":     13,
    "xtick.labelsize":    12,
    "ytick.labelsize":    12,
    "legend.fontsize":    11,
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


def save_fig(fig, name: str):
    dst = FIGURES_DIR / name
    fig.savefig(dst, dpi=DPI, bbox_inches="tight", facecolor="white")
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
                     val_range: float, fontsize: float = 12):
    tick_h = val_range * 0.025
    ax.plot([x1, x1, x2, x2], [y - tick_h, y, y, y - tick_h],
            color="#333", lw=0.75, solid_capstyle="round")
    ax.text((x1 + x2) / 2, y + val_range * 0.015, _p_stars(p),
            ha="center", va="bottom", fontsize=fontsize, color="#333")


def _sig_brackets(ax, pairs_idx, data_list, val_range, base, gap, color="#555"):
    """Draw brackets only for significant pairs (p < 0.05); skip ns."""
    k = 0
    for ii, jj in pairs_idx:
        p = _wilcoxon_p(data_list[ii], data_list[jj])
        if p >= 0.05:
            continue
        y = base + k * gap
        tick_h = val_range * 0.025
        ax.plot([ii, ii, jj, jj], [y - tick_h, y, y, y - tick_h],
                color=color, lw=0.75, solid_capstyle="round", linestyle="--")
        ax.text((ii + jj) / 2, y + val_range * 0.015, _p_stars(p),
                ha="center", va="bottom", fontsize=12, color=color)
        k += 1
    return k  # number of significant brackets drawn


def _boxplot_panel(ax, data_list, colors, seed=0, s=10, alpha=0.75):
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
    tops = []
    for vals in data_list:
        if not vals: tops.append(0.0); continue
        q3  = float(np.percentile(vals, 75))
        iqr = float(np.percentile(vals, 75) - np.percentile(vals, 25))
        tops.append(min(q3 + 1.5 * iqr, max(vals)))
    return tops


def _wilcoxon_p(vals_i, vals_j):
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
# Fig. 3c — 6-panel boxplot (6 methods)
# ═══════════════════════════════════════════════════════════════════════════

def fig3c_combined(df_naive: pd.DataFrame, df_wba: pd.DataFrame):
    rois   = list(ROI_INFO.keys())
    mord   = METHOD_ORDER
    labels = [METHOD_LABELS[m] for m in mord]
    colors = [METHOD_COLORS[m] for m in mord]
    n_meth = len(mord)

    def _naive_col(method_key: str, col: str) -> list:
        sub = df_naive[df_naive["method"] == method_key]
        return sub[col].dropna().values.tolist()

    def _sd_wba_col(col: str) -> list:
        sub = df_wba[(df_wba["method"] == "stardist") & (df_wba["attribution"] == "wba")]
        if col == "a1_umi_density" and "a1_umi_density_corrected" in df_wba.columns:
            col = "a1_umi_density_corrected"
        return sub[col].dropna().values.tolist()

    def _get_col(m: str, col: str) -> list:
        return _sd_wba_col(col) if m == "stardist_wba" else _naive_col(m, col)

    def _ned_series(m: str) -> list:
        if m == "stardist_wba":
            sub = df_wba[(df_wba["method"] == "stardist") & (df_wba["attribution"] == "wba")].set_index("roi")
        else:
            sub = df_naive[df_naive["method"] == m].set_index("roi")
        return [sub.loc[r, "ned"] if r in sub.index else np.nan for r in rois]

    ftc_data  = [_get_col(m, "a1_capture")      for m in mord]
    umi_data  = [_get_col(m, "a2_median_umi")   for m in mord]
    gene_data = [_get_col(m, "a3_median_genes")  for m in mord]
    dens_data = [_get_col(m, "a1_umi_density")   for m in mord]
    ned_data  = [_ned_series(m) for m in mord]
    ned_clean = [[v for v in d if not np.isnan(v)] for d in ned_data]

    def _dr(m: str) -> list:
        if m == "stardist_wba":
            sub = df_wba[(df_wba["method"] == "stardist") & (df_wba["attribution"] == "wba")]
            return (sub["c1_coexpr"].dropna() * 100).values.tolist()
        return (df_naive[df_naive["method"] == m]["c1_coexpr"].dropna() * 100).values.tolist()

    dr_data = [_dr(m) for m in mord]

    idx_sr   = mord.index("sr")
    idx_v12  = mord.index("v12")
    idx_p3   = mord.index("p3")
    idx_nuc  = mord.index("nuc")
    idx_pros = mord.index("proseg")
    idx_sdw  = mord.index("stardist_wba")

    fig = plt.figure(figsize=(DOUBLE_COL * 1.55, 170 * MM_TO_IN))
    gs  = fig.add_gridspec(2, 3, wspace=0.46, hspace=0.55)
    ax_d = fig.add_subplot(gs[0, 0])
    ax_e = fig.add_subplot(gs[0, 1])
    ax_f = fig.add_subplot(gs[0, 2])
    ax_g = fig.add_subplot(gs[1, 0])
    ax_h = fig.add_subplot(gs[1, 1])
    ax_i = fig.add_subplot(gs[1, 2])

    # ── Panel d: FTC ─────────────────────────────────────────────────────
    _boxplot_panel(ax_d, ftc_data, colors, seed=0)
    all_v = [v for d in ftc_data for v in d]
    rng_v = max(all_v) - min(all_v) if all_v else 0.1
    wt    = _whisker_tops(ftc_data)
    base  = max(wt) + rng_v * 0.12
    gap   = rng_v * 0.14
    for k, (ii, jj) in enumerate([(idx_sr, idx_v12), (idx_sr, idx_p3)]):
        _add_box_bracket(ax_d, ii, jj, base + k * gap,
                         _wilcoxon_p(ftc_data[ii], ftc_data[jj]), rng_v)
    n_supp = _sig_brackets(ax_d, [(idx_v12, idx_pros), (idx_v12, idx_sdw)],
                            ftc_data, rng_v, base + 2 * gap, gap, color="#777")
    ax_d.set_xticks(range(n_meth)); ax_d.set_xticklabels(labels, fontsize=12, rotation=30, ha="right")
    ax_d.set_ylabel("Transcript capture rate (FTC)", fontsize=13)
    ax_d.set_ylim(0, base + (2 + n_supp) * gap + rng_v * 0.15)
    ax_d.yaxis.set_major_locator(mticker.MultipleLocator(0.2))
    ax_d.axhline(1.0, color="#bbb", lw=0.5, ls="--", zorder=1)
    ax_d.annotate("↑ better", xy=(0.97, 0.05), xycoords="axes fraction",
                  ha="right", va="bottom", fontsize=11, color="#555")
    clean_ax(ax_d)
    ax_d.text(-0.30, 1.10, "d", transform=ax_d.transAxes,
              fontsize=13, fontweight="bold", va="top", ha="right")

    # ── Panel e: UMI/cell ─────────────────────────────────────────────────
    _boxplot_panel(ax_e, umi_data, colors, seed=10)
    all_v = [v for d in umi_data for v in d]
    rng_v = max(all_v) - min(all_v) if all_v else 200
    wt    = _whisker_tops(umi_data)
    base  = max(wt) + rng_v * 0.10
    gap   = rng_v * 0.16
    # Only key comparisons: MCseg vs SR, MCseg vs SD+WBA (most informative)
    for k, (ii, jj) in enumerate([(idx_sr, idx_v12)]):
        _add_box_bracket(ax_e, ii, jj, base + k * gap,
                         _wilcoxon_p(umi_data[ii], umi_data[jj]), rng_v)
    n_supp = _sig_brackets(ax_e, [(idx_v12, idx_sdw)],
                            umi_data, rng_v, base + gap, gap, color="#777")
    ax_e.axhline(50, color="#e74c3c", lw=0.9, ls="--", zorder=2)
    ax_e.text(n_meth - 0.5, 52, "QC threshold", fontsize=10,
              color="#e74c3c", ha="right", va="bottom")
    ax_e.set_xticks(range(n_meth)); ax_e.set_xticklabels(labels, fontsize=12, rotation=30, ha="right")
    ax_e.set_ylabel("Median UMI per cell", fontsize=13)
    ax_e.set_ylim(0, base + (1 + n_supp) * gap + rng_v * 0.10)
    ax_e.yaxis.set_major_locator(mticker.MultipleLocator(200))
    clean_ax(ax_e)
    ax_e.text(-0.30, 1.10, "e", transform=ax_e.transAxes,
              fontsize=13, fontweight="bold", va="top", ha="right")

    # ── Panel f: Genes/cell ───────────────────────────────────────────────
    _boxplot_panel(ax_f, gene_data, colors, seed=77)
    all_v = [v for d in gene_data for v in d]
    rng_v = max(all_v) - min(all_v) if all_v else 200
    wt    = _whisker_tops(gene_data)
    base  = max(wt) + rng_v * 0.10
    gap   = rng_v * 0.16
    for k, (ii, jj) in enumerate([(idx_sr, idx_v12)]):
        _add_box_bracket(ax_f, ii, jj, base + k * gap,
                         _wilcoxon_p(gene_data[ii], gene_data[jj]), rng_v)
    n_supp = _sig_brackets(ax_f, [(idx_v12, idx_sdw)],
                            gene_data, rng_v, base + gap, gap, color="#777")
    ax_f.axhline(300, color="#e74c3c", lw=0.9, ls="--", zorder=2)
    ax_f.text(n_meth - 0.5, 308, "QC threshold", fontsize=10,
              color="#e74c3c", ha="right", va="bottom")
    ax_f.set_xticks(range(n_meth)); ax_f.set_xticklabels(labels, fontsize=12, rotation=30, ha="right")
    ax_f.set_ylabel("Median genes per cell", fontsize=13)
    ax_f.set_ylim(0, base + (1 + n_supp) * gap + rng_v * 0.15)
    ax_f.yaxis.set_major_locator(mticker.MultipleLocator(200))
    clean_ax(ax_f)
    ax_f.text(-0.30, 1.10, "f", transform=ax_f.transAxes,
              fontsize=13, fontweight="bold", va="top", ha="right")

    # ── Panel g: UMI density ──────────────────────────────────────────────
    _boxplot_panel(ax_g, dens_data, colors, seed=99)
    all_v = [v for d in dens_data for v in d]
    rng_v = max(all_v) - min(all_v) if all_v else 5
    wt    = _whisker_tops(dens_data)
    base  = max(wt) + rng_v * 0.18
    gap   = rng_v * 0.18
    for k, (ii, jj) in enumerate([(idx_sr, idx_v12), (idx_sr, idx_p3)]):
        _add_box_bracket(ax_g, ii, jj, base + k * gap,
                         _wilcoxon_p(dens_data[ii], dens_data[jj]), rng_v)
    n_supp = _sig_brackets(ax_g,
                            [(idx_v12, idx_sdw), (idx_nuc, idx_v12)],
                            dens_data, rng_v, base + 2 * gap, gap, color="#777")
    ax_g.set_xticks(range(n_meth)); ax_g.set_xticklabels(labels, fontsize=12, rotation=30, ha="right")
    ax_g.set_ylabel("Size-corrected UMI density (UMI/µm²)", fontsize=12.5)
    ax_g.set_ylim(0, base + (2 + n_supp) * gap + rng_v * 0.15 if all_v else 20)
    clean_ax(ax_g)
    ax_g.text(-0.30, 1.10, "g", transform=ax_g.transAxes,
              fontsize=13, fontweight="bold", va="top", ha="right")

    # Legend — anchored to panel g's right edge in figure coordinates
    leg_handles = [mpatches.Patch(color=METHOD_COLORS[m], label=METHOD_LABELS[m], alpha=0.88)
                   for m in mord]
    # Get panel g bounding box in figure coordinates, then place legend just to its right
    fig.canvas.draw()
    bbox_i = ax_i.get_position()   # Bbox in figure fraction: x0, y0, width, height
    leg_x  = bbox_i.x1 + 0.01     # just right of panel i
    leg_y  = (bbox_i.y0 + bbox_i.y1) / 2  # vertically centred on panel i
    fig.legend(leg_handles, [METHOD_LABELS[m] for m in mord],
               fontsize=11, frameon=True, framealpha=0.95, edgecolor="#ccc",
               loc="center left", bbox_to_anchor=(leg_x, leg_y),
               bbox_transform=fig.transFigure,
               ncol=1, handlelength=1.0, handletextpad=0.5,
               borderpad=0.6, labelspacing=0.4)

    # ── Panel h: NED ──────────────────────────────────────────────────────
    _boxplot_panel(ax_h, ned_clean, colors, seed=42)
    all_v = [v for d in ned_clean for v in d]
    rng_v = max(all_v) - min(all_v) if len(all_v) > 1 else 0.1
    wt    = _whisker_tops(ned_clean)
    base  = max(wt) + rng_v * 0.16
    gap   = rng_v * 0.20
    for k, (ii, jj) in enumerate([(idx_sr, idx_v12), (idx_v12, idx_p3)]):
        h_y    = base + k * gap
        tick_h = rng_v * 0.025
        ax_h.plot([ii, ii, jj, jj], [h_y - tick_h, h_y, h_y, h_y - tick_h],
                  color="#333", lw=0.75, solid_capstyle="round")
        ax_h.text((ii + jj) / 2, h_y + rng_v * 0.015,
                  _p_stars(_wilcoxon_p(ned_data[ii], ned_data[jj])),
                  ha="center", va="bottom", fontsize=12, color="#333")
    n_supp = _sig_brackets(ax_h, [(idx_v12, idx_sdw)],
                            ned_data, rng_v, base + 2 * gap, gap, color="#777")
    y_ceil  = base + (2 + n_supp) * gap + rng_v * 0.08
    y_floor = min(all_v) - rng_v * 0.05 if all_v else 0.0
    ax_h.set_xlim(-0.65, n_meth - 0.35)
    ax_h.set_ylim(y_floor, y_ceil)
    ax_h.set_xticks(range(n_meth)); ax_h.set_xticklabels(labels, fontsize=12, rotation=30, ha="right")
    ax_h.set_ylabel("NED (Hellinger distance)", fontsize=13)
    ax_h.annotate("↑ better", xy=(0.97, 0.95), xycoords="axes fraction",
                  ha="right", va="top", fontsize=11, color="#555")
    clean_ax(ax_h)
    ax_h.text(-0.30, 1.10, "h", transform=ax_h.transAxes,
              fontsize=13, fontweight="bold", va="top", ha="right")

    # ── Panel i: Doublet rate ─────────────────────────────────────────────
    _boxplot_panel(ax_i, dr_data, colors, seed=55)
    all_v = [v for d in dr_data for v in d]
    rng_v = max(all_v) - min(all_v) if all_v else 1.0
    wt    = _whisker_tops(dr_data)
    base_i = max(wt) + rng_v * 0.18
    gap_i  = rng_v * 0.20
    _add_box_bracket(ax_i, idx_sr, idx_v12, base_i,
                     _wilcoxon_p(dr_data[idx_sr], dr_data[idx_v12]), rng_v)
    n_supp = _sig_brackets(ax_i, [(idx_v12, idx_sdw)],
                            dr_data, rng_v, base_i + gap_i, gap_i, color="#777")
    ax_i.set_xticks(range(n_meth)); ax_i.set_xticklabels(labels, fontsize=12, rotation=30, ha="right")
    ax_i.set_ylabel("Doublet rate (%)", fontsize=13)
    ax_i.set_ylim(-0.05, base_i + (1 + n_supp) * gap_i + rng_v * 0.15 if all_v else 5)
    ax_i.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{v:.1f}%"))
    ax_i.annotate("↓ better", xy=(0.97, 0.95), xycoords="axes fraction",
                  ha="right", va="top", fontsize=11, color="#555")
    clean_ax(ax_i)
    ax_i.text(-0.30, 1.10, "i", transform=ax_i.transAxes,
              fontsize=13, fontweight="bold", va="top", ha="right")

    save_fig(fig, "fig3efghij.png")   # Figure 3 panels e–j（六格箱型圖）


if __name__ == "__main__":
    print("[fig3a_transcript_metrics] Loading data ...")
    df_naive = pd.read_csv(METRICS_DIR / "metrics_ac.csv")
    df_wba   = pd.read_csv(METRICS_DIR / "metrics_ac_wba.csv")

    print("\n[Fig 3c] 6-method boxplot (panels i–vi) ...")
    fig3c_combined(df_naive, df_wba)

    print("\n✅ Done. Output:", FIGURES_DIR)
