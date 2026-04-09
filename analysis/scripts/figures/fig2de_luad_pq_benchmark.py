"""
32_fig2_luad_pq_plots.py
========================
生成 LUAD 幾何標竿的兩張新圖：

  (A) fig_pq_boxplot.png  — PQ 箱型圖（MCseg v2 vs MCseg v1, n=6 ROI, 含個別 ROI 點）
  (B) fig_pq_metrics.png  — PQ / SQ / RQ 三指標比較長條圖（帶 SD 誤差棒）

輸出至 manuscript/figures/03_luad_benchmark/
"""

from __future__ import annotations
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path

OUT_DIR = Path("/Volumes/SSD/plan_a/manuscript/figures/03_luad_benchmark")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ── Nature Communications 全域樣式 ─────────────────────────────────────────
plt.rcParams.update({
    "font.family":      "Arial",
    "font.size":        8,
    "axes.titlesize":   9,
    "axes.labelsize":   8,
    "xtick.labelsize":  7,
    "ytick.labelsize":  7,
    "legend.fontsize":  7,
    "axes.linewidth":   0.8,
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
    "xtick.major.size": 3,
    "ytick.major.size": 3,
    "pdf.fonttype":     42,
    "ps.fonttype":      42,
})

# ── 數據（逐 ROI，依 MEMORY.md 與 luad_benchmark_results.csv 確認）──────────
#  ROI 命名對齊：roi1, roi2, roi3, roi4, roi5, roi6
#  組織類型標籤
ROI_LABELS = ["roi1\n(Tumor\nBoundary)",
              "roi2\n(Stroma)",
              "roi3\n(Mixed)",
              "roi4\n(Normal-\nTumor)",
              "roi5\n(Alveolar)",
              "roi6\n(Tumor\nCore)"]
ROI_SHORT  = ["ROI1", "ROI2", "ROI3", "ROI4", "ROI5", "ROI6"]

# MCseg v2 per-ROI PQ / SQ / RQ（從 MEMORY.md 確認）
MCSEG_V2 = {
    "pq": np.array([0.6615, 0.5534, 0.5229, 0.5101, 0.4683, 0.6097]),
    "sq": np.array([0.831,  0.777,  0.761,  0.755,  0.734,  0.805]),
    "rq": np.array([0.796,  0.713,  0.687,  0.676,  0.638,  0.758]),
}

# MCseg v1 (cellpose_dilate) per-ROI（從 xenium_he_seg benchmark_results.csv）
MCSEG_V1 = {
    "pq": np.array([0.502222, 0.437619, 0.401194, 0.409299, 0.394689, 0.449939]),
    "sq": np.array([0.642549, 0.596133, 0.550167, 0.555841, 0.550859, 0.613643]),
    "rq": np.array([0.781609, 0.734096, 0.729221, 0.736361, 0.716496, 0.733226]),
}

COLOR_V2 = "#E63946"   # 鮮紅（MCseg v2 主色）
COLOR_V1  = "#457B9D"   # 深藍（MCseg v1 基準）


# ══════════════════════════════════════════════════════════════════════════════
# 圖 A：PQ 箱型圖（帶個別 ROI 資料點 + 配對連線）
# ══════════════════════════════════════════════════════════════════════════════

def plot_pq_boxplot():
    fig, ax = plt.subplots(figsize=(88/25.4, 70/25.4))   # 88mm wide, ~70mm tall

    data    = [MCSEG_V1["pq"], MCSEG_V2["pq"]]
    labels  = ["MCseg v1\n(Baseline)", "MCseg v2\n(Ours)"]
    colors  = [COLOR_V1, COLOR_V2]
    x_pos   = [1, 2]

    # 箱型圖
    bp = ax.boxplot(data, positions=x_pos, widths=0.32, patch_artist=True,
                    medianprops=dict(color="white", linewidth=1.5),
                    whiskerprops=dict(linewidth=0.8),
                    capprops=dict(linewidth=0.8),
                    flierprops=dict(marker="", markersize=0))

    for patch, color in zip(bp["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.55)
        patch.set_linewidth(0.8)

    # 個別 ROI 點 + 配對連線
    jitter_p3  = np.linspace(-0.08, 0.08, 6)
    jitter_v12 = np.linspace(-0.08, 0.08, 6)
    np.random.seed(42)
    np.random.shuffle(jitter_p3)
    np.random.shuffle(jitter_v12)

    xs_p3  = np.ones(6) + jitter_p3
    xs_v12 = np.ones(6) * 2 + jitter_v12

    for i in range(6):
        ax.plot([xs_p3[i], xs_v12[i]],
                [MCSEG_V1["pq"][i], MCSEG_V2["pq"][i]],
                color="grey", linewidth=0.5, alpha=0.5, zorder=2)

    ax.scatter(xs_p3,  MCSEG_V1["pq"],  color=COLOR_V1,  s=20, zorder=3, linewidths=0.4,
               edgecolors="white")
    ax.scatter(xs_v12, MCSEG_V2["pq"], color=COLOR_V2, s=20, zorder=3, linewidths=0.4,
               edgecolors="white")

    # ROI 標籤（標在線段中點稍右）
    for i in range(6):
        mid_x = (xs_p3[i] + xs_v12[i]) / 2
        mid_y = (MCSEG_V1["pq"][i] + MCSEG_V2["pq"][i]) / 2
        ax.annotate(ROI_SHORT[i],
                    xy=(mid_x + 0.04, mid_y),
                    fontsize=5.5, color="dimgray",
                    va="center", ha="left")

    # 顯著性標記（Wilcoxon，n=6，配對）
    from scipy.stats import wilcoxon
    stat, p = wilcoxon(MCSEG_V1["pq"], MCSEG_V2["pq"], alternative="less")
    star = "***" if p < 0.001 else ("**" if p < 0.01 else ("*" if p < 0.05 else "ns"))
    y_sig = max(MCSEG_V2["pq"]) + 0.03
    ax.plot([1, 2], [y_sig, y_sig], color="k", linewidth=0.8)
    ax.text(1.5, y_sig + 0.01, f"{star}\np = {p:.3f}", ha="center", fontsize=6.5)

    # 均值標記
    for x, vals, c in zip(x_pos, data, colors):
        ax.hlines(np.mean(vals), x - 0.20, x + 0.20,
                  colors=c, linewidths=1.5, linestyles="--", zorder=4)

    ax.set_xlim(0.5, 2.7)
    ax.set_xticks(x_pos)
    ax.set_xticklabels(labels)
    ax.set_ylabel("Panoptic Quality (PQ)", fontsize=8)
    ax.set_ylim(0.32, max(MCSEG_V2["pq"]) + 0.10)
    ax.set_title("LUAD Geometric Benchmark\n(n = 6 ROIs, Xenium GT)", fontsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", alpha=0.25, linewidth=0.5)

    # 圖例（均值虛線）
    handles = [
        mpatches.Patch(color=COLOR_V1,  alpha=0.6, label=f"MCseg v1  (mean={np.mean(MCSEG_V1['pq']):.3f})"),
        mpatches.Patch(color=COLOR_V2, alpha=0.6, label=f"MCseg v2 (mean={np.mean(MCSEG_V2['pq']):.3f})"),
    ]
    ax.legend(handles=handles, fontsize=6, loc="upper left",
              frameon=True, framealpha=0.7, edgecolor="none")

    out = OUT_DIR / "fig_pq_boxplot.png"
    plt.tight_layout(pad=0.4)
    plt.savefig(out, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"✅ Saved: {out}")


# ══════════════════════════════════════════════════════════════════════════════
# 圖 B：PQ / SQ / RQ 三指標比較（grouped bar + SD 誤差棒）
# ══════════════════════════════════════════════════════════════════════════════

def plot_pq_metrics():
    fig, ax = plt.subplots(figsize=(120/25.4, 72/25.4))   # 120mm wide

    metrics = ["PQ", "SQ", "RQ"]
    x       = np.arange(len(metrics))
    w       = 0.32

    means_v12 = [np.mean(MCSEG_V2[k.lower()]) for k in metrics]
    stds_v12  = [np.std(MCSEG_V2[k.lower()])  for k in metrics]
    means_p3  = [np.mean(MCSEG_V1[k.lower()])  for k in metrics]
    stds_p3   = [np.std(MCSEG_V1[k.lower()])   for k in metrics]

    bars_p3 = ax.bar(x - w/2, means_p3, w,
                     yerr=stds_p3, capsize=3,
                     color=COLOR_V1, alpha=0.7, linewidth=0.6,
                     error_kw=dict(elinewidth=0.7, capthick=0.7),
                     label="MCseg v1 (Baseline)", zorder=3)

    bars_v12 = ax.bar(x + w/2, means_v12, w,
                      yerr=stds_v12, capsize=3,
                      color=COLOR_V2, alpha=0.7, linewidth=0.6,
                      error_kw=dict(elinewidth=0.7, capthick=0.7),
                      label="MCseg v2 (Ours)", zorder=3)

    # 個別 ROI 散點（與 bar 疊加）
    rng = np.random.default_rng(0)
    for i, key in enumerate(["pq", "sq", "rq"]):
        j3  = rng.uniform(-0.08, 0.08, 6)
        j12 = rng.uniform(-0.08, 0.08, 6)
        ax.scatter(i - w/2 + j3,  MCSEG_V1[key],  s=14, color=COLOR_V1,
                   zorder=4, linewidths=0.3, edgecolors="white", alpha=0.8)
        ax.scatter(i + w/2 + j12, MCSEG_V2[key], s=14, color=COLOR_V2,
                   zorder=4, linewidths=0.3, edgecolors="white", alpha=0.8)

    # Δ 標記（在每組正上方標出提升量）
    for i, key in enumerate(["pq", "sq", "rq"]):
        delta = np.mean(MCSEG_V2[key]) - np.mean(MCSEG_V1[key])
        y_top = max(np.mean(MCSEG_V2[key]) + np.std(MCSEG_V2[key]),
                    np.mean(MCSEG_V1[key])  + np.std(MCSEG_V1[key])) + 0.02
        ax.text(i, y_top, f"+{delta:.3f}", ha="center", fontsize=6.5,
                color="#C0392B", fontweight="bold")

    ax.set_xticks(x)
    ax.set_xticklabels(["PQ\n(Panoptic Quality)",
                         "SQ\n(Segmentation Quality)",
                         "RQ\n(Recognition Quality)"], fontsize=7.5)
    ax.set_ylabel("Score (mean ± SD, n = 6 ROIs)", fontsize=8)
    ax.set_ylim(0, 1.05)
    ax.set_title("MCseg v2 vs MCseg v1: Three-Component Segmentation Metrics\n(LUAD, Xenium Ground Truth)",
                 fontsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", alpha=0.25, linewidth=0.5)
    ax.legend(fontsize=7, loc="upper left",
              frameon=True, framealpha=0.7, edgecolor="none")

    # 參考線（0.5）
    ax.axhline(0.5, color="grey", linewidth=0.5, linestyle=":", alpha=0.6)

    out = OUT_DIR / "fig_pq_metrics.png"
    plt.tight_layout(pad=0.4)
    plt.savefig(out, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"✅ Saved: {out}")


# ══════════════════════════════════════════════════════════════════════════════
# 圖 C：逐 ROI PQ 條形圖（取代原 fig_pq_bar.png — 並排條形 + 組織類型標示）
# ══════════════════════════════════════════════════════════════════════════════

def plot_pq_bar_by_roi():
    """原版 fig_pq_bar 的改良版：逐 ROI 並排，以組織類型著色背景。"""
    fig, ax = plt.subplots(figsize=(140/25.4, 70/25.4))

    x   = np.arange(6)
    w   = 0.36
    rois_ordered = ["roi1", "roi2", "roi3", "roi4", "roi5", "roi6"]
    tissue_types = ["Tumor\nBoundary", "Tumor\nStroma", "Mixed\nTumor-Stroma",
                    "Normal-Tumor\nInterface", "Alveolar", "Tumor\nCore"]

    # 背景色帶（區分組織類型）
    bg_colors = ["#fff5f5", "#f5f0ff", "#f0f5ff", "#f5fff0", "#fffff0", "#fff5f5"]
    for i, bc in enumerate(bg_colors):
        ax.axvspan(i - 0.5, i + 0.5, color=bc, alpha=0.4, zorder=0)

    bars_p3 = ax.bar(x - w/2, MCSEG_V1["pq"], w,
                     color=COLOR_V1, alpha=0.80, linewidth=0.5,
                     label="MCseg v1 (Baseline)", zorder=3)
    bars_v12 = ax.bar(x + w/2, MCSEG_V2["pq"], w,
                      color=COLOR_V2, alpha=0.80, linewidth=0.5,
                      label="MCseg v2 (Ours)", zorder=3)

    # 數值標籤（MCseg v2 bar 頂端）
    for i, (p, v) in enumerate(zip(MCSEG_V1["pq"], MCSEG_V2["pq"])):
        ax.text(i - w/2, p + 0.008, f"{p:.3f}", ha="center", va="bottom",
                fontsize=5.5, color=COLOR_V1)
        ax.text(i + w/2, v + 0.008, f"{v:.3f}", ha="center", va="bottom",
                fontsize=5.5, color=COLOR_V2)

    # 提升量 Δ
    for i in range(6):
        delta = MCSEG_V2["pq"][i] - MCSEG_V1["pq"][i]
        y_mid = (MCSEG_V1["pq"][i] + MCSEG_V2["pq"][i]) / 2
        ax.text(i + w/2 + 0.22, y_mid, f"Δ+{delta:.3f}",
                ha="left", va="center", fontsize=5, color="#C0392B")

    # mean 水平線
    ax.axhline(np.mean(MCSEG_V2["pq"]), color=COLOR_V2, linewidth=0.8,
               linestyle="--", alpha=0.8,
               label=f"MCseg v2 mean = {np.mean(MCSEG_V2['pq']):.3f}")
    ax.axhline(np.mean(MCSEG_V1["pq"]),  color=COLOR_V1,  linewidth=0.8,
               linestyle="--", alpha=0.8,
               label=f"MCseg v1  mean = {np.mean(MCSEG_V1['pq']):.3f}")

    ax.set_xticks(x)
    ax.set_xticklabels([f"{r.upper()}\n{t}" for r, t in zip(rois_ordered, tissue_types)],
                       fontsize=6.5)
    ax.set_ylabel("Panoptic Quality (PQ)", fontsize=8)
    ax.set_ylim(0, max(MCSEG_V2["pq"]) + 0.12)
    ax.set_title("LUAD Geometric Benchmark: PQ per ROI (MCseg v2 vs MCseg v1, Xenium GT)",
                 fontsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", alpha=0.25, linewidth=0.5)
    ax.legend(fontsize=6.5, loc="upper right",
              frameon=True, framealpha=0.7, edgecolor="none", ncol=2)

    out = OUT_DIR / "fig_pq_bar_revised.png"
    plt.tight_layout(pad=0.4)
    plt.savefig(out, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"✅ Saved: {out}")


# ══════════════════════════════════════════════════════════════════════════════
# 圖 D：主圖 Fig. 2（3-row layout）
# ══════════════════════════════════════════════════════════════════════════════

def plot_pq_3boxplots():
    """
    Main Fig. 2: 3-row layout
      Row 0 (a/b/c):  PQ | SQ | RQ boxplots — geometric benchmark (n=6 ROIs)
      Row 1 (d):      ROI10 seg image  |  AT1/AT2 dotplot
                      Visium HD 18K unbiased panel enables AT2 ID via SFTPC/SFTPB/SFTPA1
                      (absent from Xenium 5K targeted panel)
      Row 2 (e):      ROI9 seg image   |  SPP1⁺ TAM bar chart
                      SPP1 absent from Xenium 5K targeted panel; captured by Visium HD 18K

    Output: fig_pq_metrics.png
    """
    import pandas as pd
    from scipy.stats import wilcoxon
    import matplotlib.gridspec as gridspec
    import matplotlib.patches as mpatches
    import matplotlib.colors as mcolors
    from matplotlib.lines import Line2D
    from PIL import Image

    DATA_DIR = Path("/Volumes/SSD/plan_a/manuscript/data")
    ARCHIVE  = Path("/Volumes/SSD/plan_a/manuscript/figures/_archive/root")

    df_spp1 = pd.read_csv(DATA_DIR / "luad_roi9_spp1_macrophage_markers.csv")

    # ── Dotplot data (precomputed from roi10_v12.h5ad) ────────────────────────
    # Cells: AT2=488 (45.0%), AT1 excl.=23 (2.1%), Other=574 (52.9%)
    DOT_GENES  = ['SFTPC', 'SFTPB', 'SFTPA1', 'LAMP3', 'AGER', 'CAV1', 'HOPX']
    DOT_LABELS = ['AT2\n(n=488, 45%)', 'AT1\n(n=23, 2%)', 'Other\n(n=574, 53%)']

    # pct_positive [group × gene]: AT2 row / AT1-excl row / Other row
    pct_arr = np.array([
        [75.8, 49.4, 26.4,  1.4,  2.0,  2.3,  2.9],
        [ 0.0,  0.0,  0.0,  0.0, 30.4, 43.5, 21.7],
        [ 0.0,  0.0,  0.0,  1.2,  0.0,  0.0,  0.0],
    ])
    # mean log(1+UMI) [group × gene]
    mean_arr = np.array([
        [0.809, 0.471, 0.221, 0.010, 0.014, 0.016, 0.021],
        [0.000, 0.000, 0.000, 0.000, 0.246, 0.301, 0.151],
        [0.000, 0.000, 0.000, 0.010, 0.000, 0.000, 0.000],
    ])

    # ── Figure layout ─────────────────────────────────────────────────────────
    fig = plt.figure(figsize=(183 / 25.4, 224 / 25.4))
    gs  = gridspec.GridSpec(
        3, 12,
        height_ratios=[1.0, 1.30, 1.00],
        hspace=0.62, wspace=0.45,
        left=0.065, right=0.975,
        top=0.955, bottom=0.042,
    )

    ax_a     = fig.add_subplot(gs[0, 0:4])   # PQ
    ax_b     = fig.add_subplot(gs[0, 4:8])   # SQ
    ax_c     = fig.add_subplot(gs[0, 8:12])  # RQ
    ax_d_img = fig.add_subplot(gs[1, 0:4])   # ROI10 seg image
    ax_d_dot = fig.add_subplot(gs[1, 4:12])  # AT1/AT2 dotplot
    ax_e_img = fig.add_subplot(gs[2, 0:4])   # ROI9 seg image
    ax_e_bar = fig.add_subplot(gs[2, 4:12])  # SPP1+ bar chart

    rng = np.random.default_rng(42)

    # ── Row 0: PQ / SQ / RQ boxplots ─────────────────────────────────────────
    for ax, (key, label) in zip([ax_a, ax_b, ax_c],
                                 [("pq", "PQ"), ("sq", "SQ"), ("rq", "RQ")]):
        v1 = MCSEG_V1[key]
        v2 = MCSEG_V2[key]

        bp = ax.boxplot([v1, v2], positions=[1, 2], widths=0.36,
                        patch_artist=True,
                        medianprops=dict(color="white", linewidth=1.5),
                        whiskerprops=dict(linewidth=0.8),
                        capprops=dict(linewidth=0.8),
                        flierprops=dict(marker=""))
        bp["boxes"][0].set(facecolor=COLOR_V1, alpha=0.55, linewidth=0.8)
        bp["boxes"][1].set(facecolor=COLOR_V2, alpha=0.55, linewidth=0.8)

        j1 = rng.uniform(-0.07, 0.07, 6)
        j2 = rng.uniform(-0.07, 0.07, 6)
        xs1 = np.ones(6) + j1
        xs2 = np.full(6, 2.0) + j2
        for i in range(6):
            ax.plot([xs1[i], xs2[i]], [v1[i], v2[i]],
                    color="grey", lw=0.5, alpha=0.45, zorder=2)
        ax.scatter(xs1, v1, color=COLOR_V1, s=22, zorder=3, lw=0.4, edgecolors="white")
        ax.scatter(xs2, v2, color=COLOR_V2, s=22, zorder=3, lw=0.4, edgecolors="white")

        _, p = wilcoxon(v1, v2, alternative="two-sided")
        star = ("***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns")
        dr = max(max(v1), max(v2)) - min(min(v1), min(v2))
        y_bar  = max(max(v1), max(v2)) + dr * 0.06
        y_star = y_bar + dr * 0.025
        y_delt = y_bar + dr * 0.10
        ax.plot([1, 2], [y_bar, y_bar], color="k", lw=0.8)
        ax.text(1.5, y_star, star, ha="center", fontsize=7.5)
        delta = np.mean(v2) - np.mean(v1)
        ax.text(1.5, y_delt, f"+{delta:.3f}", ha="center", fontsize=8,
                color="#C0392B", fontweight="bold")

        ax.set_xlim(0.5, 2.75)
        ax.set_xticks([1, 2])
        ax.set_xticklabels(["MCseg v1", "MCseg v2"], fontsize=7)
        ax.set_ylabel("Score (n = 6 ROIs)" if key == "pq" else "", fontsize=7.5)
        ax.set_title(label, fontsize=9, fontweight="bold")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.grid(axis="y", alpha=0.25, lw=0.5)
        if key != "rq":
            ax.axhline(0.5, color="grey", lw=0.5, ls=":", alpha=0.6)
        ax.set_ylim(min(min(v1), min(v2)) - dr * 0.05, y_delt + dr * 0.08)

    # ── Panel d-img: ROI10 segmentation image ────────────────────────────────
    img10 = np.array(Image.open(ARCHIVE / "roi10_panel_A_winner_v12.png").convert("RGB"))
    crop_top10 = int(img10.shape[0] * 0.08)
    ax_d_img.imshow(img10[crop_top10:], aspect="auto")
    ax_d_img.axis("off")
    ax_d_img.set_title("MCseg v2 — ROI10\nNormal Alveolar (n=1,085 cells)",
                       fontsize=6.5, fontweight="bold", pad=2)

    # ── Panel d-dot: AT1/AT2 dotplot ──────────────────────────────────────────
    n_g  = len(DOT_GENES)
    n_gr = len(DOT_LABELS)

    cmap_dot = plt.cm.YlOrRd
    norm_dot = mcolors.Normalize(vmin=0, vmax=0.85)

    xs_all, ys_all, sz_all, ce_all = [], [], [], []
    for gi in range(n_gr):
        y = n_gr - 1 - gi   # AT2 on top (y=2), Other on bottom (y=0)
        for gei in range(n_g):
            xs_all.append(gei)
            ys_all.append(y)
            sz_all.append(max((pct_arr[gi, gei] / 100) * 400, 4))
            ce_all.append(mean_arr[gi, gei])

    sc = ax_d_dot.scatter(xs_all, ys_all, s=sz_all, c=ce_all,
                          cmap=cmap_dot, norm=norm_dot,
                          edgecolors="grey", linewidths=0.3, zorder=3)

    # Colorbar (inset, upper right)
    cbar_ax = ax_d_dot.inset_axes([0.875, 0.08, 0.028, 0.52])
    cb = fig.colorbar(sc, cax=cbar_ax, orientation="vertical")
    cb.set_label("Mean log(1+UMI)", fontsize=5.5, labelpad=3)
    cb.ax.tick_params(labelsize=5.5, pad=1)
    cb.set_ticks([0, 0.2, 0.4, 0.6, 0.8])

    # Size legend (lower left)
    leg_handles = [
        Line2D([0], [0], marker='o', color='w',
               markerfacecolor='#aaa', markeredgecolor='grey',
               markeredgewidth=0.4,
               markersize=np.sqrt(max((p / 100) * 400, 4) / np.pi) * 2,
               label=f'{p}%')
        for p in [25, 50, 75]
    ]
    ax_d_dot.legend(handles=leg_handles, title="% cells\npositive",
                    loc="lower left", fontsize=5.5, title_fontsize=5.5,
                    handletextpad=0.3, labelspacing=0.4,
                    framealpha=0.85, edgecolor="#ccc")

    ax_d_dot.set_xlim(-0.6, n_g - 0.4)
    ax_d_dot.set_ylim(-0.6, n_gr - 0.4)
    ax_d_dot.set_xticks(range(n_g))
    ax_d_dot.set_xticklabels(DOT_GENES, fontsize=7.5, rotation=45, ha="right")
    ax_d_dot.set_yticks(range(n_gr))
    ax_d_dot.set_yticklabels(DOT_LABELS[::-1], fontsize=7.5)
    ax_d_dot.grid(alpha=0.18, lw=0.4)
    ax_d_dot.spines["top"].set_visible(False)
    ax_d_dot.spines["right"].set_visible(False)
    ax_d_dot.set_title(
        "AT2/AT1 Markers in ROI10 (Normal Alveolar) — Visium HD, MCseg v2 (n=1,085 cells)\n"
        "SFTPC/SFTPB/SFTPA1 absent from Xenium 5K targeted panel  |  LAMP3: Xenium AT2 marker",
        fontsize=7, fontweight="bold", pad=4
    )

    # Color x-tick labels by gene category
    # red = Visium HD only (AT2), purple = Xenium AT2 marker, blue = AT1
    xtick_colors = ["#C0392B", "#C0392B", "#C0392B",   # SFTPC SFTPB SFTPA1
                    "#7B2D8B",                           # LAMP3
                    "#2166AC", "#2166AC", "#2166AC"]     # AGER CAV1 HOPX
    for tick, col in zip(ax_d_dot.get_xticklabels(), xtick_colors):
        tick.set_color(col)

    # Dashed separators between gene groups
    ax_d_dot.axvline(2.5, color="grey", lw=0.6, ls="--", alpha=0.4)
    ax_d_dot.axvline(3.5, color="grey", lw=0.6, ls="--", alpha=0.4)

    # ── Panel e-img: ROI9 segmentation image ─────────────────────────────────
    img9 = np.array(Image.open(ARCHIVE / "roi9_panel_A_winner.png").convert("RGB"))
    crop_top9 = int(img9.shape[0] * 0.08)
    ax_e_img.imshow(img9[crop_top9:], aspect="auto")
    ax_e_img.axis("off")
    ax_e_img.set_title("MCseg v2 — ROI9\nPigmented Macrophage Zone (n=8,393 cells)",
                       fontsize=6.5, fontweight="bold", pad=2)

    # ── Panel e-bar: SPP1+ macrophage markers ────────────────────────────────
    genes  = df_spp1["gene"].tolist()
    pcts   = df_spp1["pct_positive"].tolist()
    colors_bar = ["#C0392B" if g == "SPP1" else "#E8A0A0" for g in genes]
    x_pos  = np.arange(len(genes))
    bars = ax_e_bar.bar(x_pos, pcts, color=colors_bar, alpha=0.85,
                        edgecolor="white", lw=0.5)
    for bar, pct in zip(bars, pcts):
        ax_e_bar.text(bar.get_x() + bar.get_width() / 2, pct + 0.03,
                      f"{pct:.2f}%", ha="center", va="bottom",
                      fontsize=6.5, fontweight="bold")

    ax_e_bar.set_xticks(x_pos)
    ax_e_bar.set_xticklabels(genes, fontsize=8)
    ax_e_bar.set_ylabel("% Cells Positive\n(MCseg v2, ROI9, n=8,393 cells)", fontsize=7)
    ax_e_bar.set_ylim(0, max(pcts) * 1.55)
    ax_e_bar.set_title(
        "SPP1⁺ TAM Markers — ROI9 Pigmented Macrophage Zone\n"
        "(Visium HD 18K unbiased panel; SPP1 absent from Xenium 5K targeted panel)",
        fontsize=7, fontweight="bold"
    )
    ax_e_bar.spines["top"].set_visible(False)
    ax_e_bar.spines["right"].set_visible(False)
    ax_e_bar.grid(axis="y", alpha=0.25, lw=0.5)

    # Annotation: Xenium gene coverage note
    ax_e_bar.text(
        0.98, 0.97,
        "CD68 / TREM2 / MMP9 / CD163:\n  also in Xenium panel\nSPP1: absent from Xenium panel",
        transform=ax_e_bar.transAxes, fontsize=5.5, ha="right", va="top",
        color="#555", style="italic",
        bbox=dict(boxstyle="round,pad=0.3", fc="#fffff8", ec="#ddd", lw=0.5)
    )

    # ── Panel labels (a–e) ────────────────────────────────────────────────────
    label_axes = [ax_a, ax_b, ax_c, ax_d_img, ax_e_img]
    for ax, letter in zip(label_axes, ["a", "b", "c", "d", "e"]):
        ax.text(-0.10, 1.06, letter, transform=ax.transAxes,
                fontsize=11, fontweight="bold", va="top", ha="right")

    out = OUT_DIR / "fig_pq_metrics.png"
    fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"✅ Saved: {out}")


def plot_fig2abc():
    """
    Fig. 2a,b,c only — PQ / SQ / RQ boxplots (1-row, 3-panel).
    Output: fig2abc.png
    """
    from scipy.stats import wilcoxon

    fig = plt.figure(figsize=(183 / 25.4, 72 / 25.4))
    import matplotlib.gridspec as gridspec
    gs = gridspec.GridSpec(
        1, 3,
        hspace=0.0, wspace=0.40,
        left=0.08, right=0.97,
        top=0.88, bottom=0.14,
    )
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[0, 2])

    rng = np.random.default_rng(42)

    for ax, (key, label) in zip([ax_a, ax_b, ax_c],
                                  [("pq", "PQ"), ("sq", "SQ"), ("rq", "RQ")]):
        v1 = MCSEG_V1[key]
        v2 = MCSEG_V2[key]

        bp = ax.boxplot([v1, v2], positions=[1, 2], widths=0.36,
                        patch_artist=True,
                        medianprops=dict(color="white", linewidth=1.5),
                        whiskerprops=dict(linewidth=0.8),
                        capprops=dict(linewidth=0.8),
                        flierprops=dict(marker=""))
        bp["boxes"][0].set(facecolor=COLOR_V1, alpha=0.55, linewidth=0.8)
        bp["boxes"][1].set(facecolor=COLOR_V2, alpha=0.55, linewidth=0.8)

        j1 = rng.uniform(-0.07, 0.07, 6)
        j2 = rng.uniform(-0.07, 0.07, 6)
        xs1 = np.ones(6) + j1
        xs2 = np.full(6, 2.0) + j2
        for i in range(6):
            ax.plot([xs1[i], xs2[i]], [v1[i], v2[i]],
                    color="grey", lw=0.5, alpha=0.45, zorder=2)
        ax.scatter(xs1, v1, color=COLOR_V1, s=22, zorder=3, lw=0.4, edgecolors="white")
        ax.scatter(xs2, v2, color=COLOR_V2, s=22, zorder=3, lw=0.4, edgecolors="white")

        _, p = wilcoxon(v1, v2, alternative="two-sided")
        star = ("***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns")
        dr = max(max(v1), max(v2)) - min(min(v1), min(v2))
        y_bar  = max(max(v1), max(v2)) + dr * 0.06
        y_star = y_bar + dr * 0.025
        y_delt = y_bar + dr * 0.10
        ax.plot([1, 2], [y_bar, y_bar], color="k", lw=0.8)
        ax.text(1.5, y_star, star, ha="center", fontsize=7.5)
        delta = np.mean(v2) - np.mean(v1)
        ax.text(1.5, y_delt, f"+{delta:.3f}", ha="center", fontsize=8,
                color="#C0392B", fontweight="bold")

        ax.set_xlim(0.5, 2.75)
        ax.set_xticks([1, 2])
        ax.set_xticklabels(["MCseg v1", "MCseg v2"], fontsize=7)
        ax.set_ylabel("Score (n = 6 ROIs)" if key == "pq" else "", fontsize=7.5)
        ax.set_title(label, fontsize=9, fontweight="bold")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.grid(axis="y", alpha=0.25, lw=0.5)
        if key != "rq":
            ax.axhline(0.5, color="grey", lw=0.5, ls=":", alpha=0.6)
        ax.set_ylim(min(min(v1), min(v2)) - dr * 0.05, y_delt + dr * 0.08)

    # Panel labels
    for ax, letter in zip([ax_a, ax_b, ax_c], ["a", "b", "c"]):
        ax.text(-0.10, 1.06, letter, transform=ax.transAxes,
                fontsize=11, fontweight="bold", va="top", ha="right")

    out = OUT_DIR / "fig2abc.png"
    fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"✅ Saved: {out}")


if __name__ == "__main__":
    print("Generating LUAD PQ figures...")
    plot_fig2abc()             # Fig. 2a,b,c → fig2abc.png
    plot_pq_3boxplots()        # Main Fig. 2  → fig_pq_metrics.png
    plot_pq_boxplot()          # Supp Fig.    → fig_pq_boxplot.png
    plot_pq_bar_by_roi()       # Supp Fig.    → fig_pq_bar_revised.png
    print("\nDone. Output:")
    print("  fig_pq_metrics.png     — 主圖 Fig. 2: PQ/SQ/RQ 三欄箱型圖")
    print("  fig_pq_boxplot.png     — 補充圖: PQ 箱型圖（含 ROI 標籤，詳細版）")
    print("  fig_pq_bar_revised.png — 補充圖: 逐 ROI PQ 條形圖")
