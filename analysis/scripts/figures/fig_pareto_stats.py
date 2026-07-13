"""
fig_pareto_stats.py
====================
Capture–Purity Pareto 散點圖 + 全指標統計檢定

指標：
  - FTC (a1_capture): Fraction of Transcripts in Cells ↑
  - UMI density (corrected): UMI/µm² ↑
  - NED: Normalized Expression Distance ↑
  - Coexpr (c1_coexpr): Impossible co-expression rate ↓

方法：
  主要組  — SR, MCseg (v12), 2Cseg (p3), NUC (nuc)  [naive attribution]
  補充對照 — StarDist+lookup, StarDist+WBA           [from metrics_ac_combined & wba]
           — Proseg                                  [from metrics_ac.csv, naive]

統計：
  所有配對以 Wilcoxon signed-rank test（n=15 配對 ROI），
  MCseg 為參考，比較每個對手；另對四個主要方法跑 Friedman test。
  多重比較用 Bonferroni 校正（所有配對）。

執行：
  cd /Volumes/SSD/plan_a/crc_transcript_attribution
  uv run python scripts/analysis/fig_pareto_stats.py

輸出：
  results/figures/fig_pareto_capture_purity.png
  results/figures/stats_all_metrics.csv
"""

from __future__ import annotations

import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore", category=FutureWarning)

ROOT     = Path(__file__).resolve().parent.parent.parent
FIG_DIR  = ROOT / "results" / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)
STATS_OUT = FIG_DIR / "stats_all_metrics.csv"

# ─── 載入資料 ──────────────────────────────────────────────────────────────────

df_naive = pd.read_csv(ROOT / "results/metrics/metrics_ac.csv")
df_comb  = pd.read_csv(ROOT / "results/metrics/metrics_ac_combined.csv")
df_wba   = pd.read_csv(ROOT / "results/metrics/metrics_ac_wba.csv")


def get_per_roi(method_key: str, label: str, df: pd.DataFrame,
                attr: str | None = None) -> pd.DataFrame:
    sub = df[df["method"] == method_key].copy()
    if attr is not None and "attribution" in df.columns:
        sub = sub[sub["attribution"] == attr]
    sub = sub.sort_values("roi").reset_index(drop=True)
    assert len(sub) == 15, f"{label}: expected 15 ROIs, got {len(sub)}"
    umi_col = ("a1_umi_density_corrected"
               if "a1_umi_density_corrected" in sub.columns
               and sub["a1_umi_density_corrected"].notna().all()
               else "a1_umi_density")
    return pd.DataFrame({
        "roi":    sub["roi"].values,
        "label":  label,
        "ftc":    sub["a1_capture"].values,
        "umi":    sub[umi_col].values,
        "ned":    sub["ned"].values,
        "coexpr": sub["c1_coexpr"].values * 100,
    })


methods_data: dict[str, pd.DataFrame] = {
    "SR":              get_per_roi("sr",       "SR",              df_naive),
    "MCseg":           get_per_roi("v12",      "MCseg",           df_naive),
    "2Cseg":           get_per_roi("p3",       "2Cseg",           df_naive),
    "NUC":             get_per_roi("nuc",      "NUC",             df_naive),
    "Proseg":          get_per_roi("proseg",   "Proseg",          df_naive),
    "StarDist+lookup": get_per_roi("stardist", "StarDist+lookup", df_comb, attr="naive"),
    "StarDist+WBA":    get_per_roi("stardist", "StarDist+WBA",    df_wba),
}

# ─── 均值摘要 ─────────────────────────────────────────────────────────────────

means = {
    label: {"FTC": df["ftc"].mean(), "UMI": df["umi"].mean(),
            "NED": df["ned"].mean(), "Coexpr": df["coexpr"].mean()}
    for label, df in methods_data.items()
}
means_df = pd.DataFrame(means).T
print("\n=== 各方法均值 ===")
print(means_df.round(3).to_string())

# ─── 統計檢定 ─────────────────────────────────────────────────────────────────

metrics = {"FTC": "ftc", "UMI": "umi", "NED": "ned", "Coexpr": "coexpr"}

# MCseg vs 每個其他方法 + 額外關注配對
ref = "MCseg"
pairs = [(ref, k) for k in methods_data if k != ref]
pairs += [("SR", "StarDist+WBA"), ("SR", "StarDist+lookup")]
n_tests = len(pairs) * len(metrics)

rows = []
for m_name, m_col in metrics.items():
    for a, b in pairs:
        va = methods_data[a][m_col].values
        vb = methods_data[b][m_col].values
        w_stat, p_raw = stats.wilcoxon(va, vb, alternative="two-sided")
        p_bonf = min(p_raw * n_tests, 1.0)
        # rank-biserial effect size
        n = len(va)
        r_rb = 1 - (2 * w_stat) / (n * (n + 1) / 2)
        rows.append({
            "Metric":     m_name,
            "Comparison": f"{a} vs {b}",
            "W":          w_stat,
            "p_raw":      p_raw,
            "p_bonf":     p_bonf,
            "r_rb":       r_rb,
            "mean_A":     va.mean(),
            "mean_B":     vb.mean(),
            "direction":  "A>B" if va.mean() > vb.mean() else "A<B",
        })

def sig_stars(p: float) -> str:
    if p < 0.001: return "***"
    if p < 0.01:  return "**"
    if p < 0.05:  return "*"
    return "ns"

stats_df = pd.DataFrame(rows)
stats_df["sig_raw"]  = stats_df["p_raw"].apply(sig_stars)
stats_df["sig_bonf"] = stats_df["p_bonf"].apply(sig_stars)
stats_df.to_csv(STATS_OUT, index=False, float_format="%.4f")
print(f"\n統計結果已儲存：{STATS_OUT}")

print("\n=== MCseg 配對統計 ===")
mask = stats_df["Comparison"].str.startswith("MCseg")
cols = ["Metric", "Comparison", "mean_A", "mean_B", "direction",
        "p_raw", "sig_raw", "p_bonf", "sig_bonf", "r_rb"]
print(stats_df[mask][cols].to_string(index=False))

print("\n=== Friedman test（SR / MCseg / 2Cseg / NUC） ===")
for m_name, m_col in metrics.items():
    groups = [methods_data[k][m_col].values for k in ["SR", "MCseg", "2Cseg", "NUC"]]
    chi2, p_f = stats.friedmanchisquare(*groups)
    print(f"  {m_name}: χ²={chi2:.2f}, p={p_f:.4f} {sig_stars(p_f)}")

# ─── Pareto 散點圖（Bioinformatics 發表風格）────────────────────────────────

# Colorbrewer-safe palette，色盲友善
COLORS = {
    "SR":              "#D62728",   # red
    "MCseg":           "#1F77B4",   # blue  (主角)
    "2Cseg":           "#9467BD",   # purple
    "NUC":             "#FF7F0E",   # orange
    "Proseg":          "#2CA02C",   # green  (supplementary)
    "StarDist+lookup": "#7F7F7F",   # grey   (supplementary)
    "StarDist+WBA":    "#17BECF",   # cyan   (supplementary)
}
MARKERS = {
    "SR":              "o",
    "MCseg":           "o",
    "2Cseg":           "s",
    "NUC":             "^",
    "Proseg":          "D",
    "StarDist+lookup": "v",
    "StarDist+WBA":    "P",
}
HOLLOW  = {"Proseg", "StarDist+lookup", "StarDist+WBA"}
MAIN    = {"SR", "MCseg", "2Cseg", "NUC"}

LABEL_OFFSETS = {
    "SR":              ( 0.013,  0.002),
    "MCseg":           ( 0.013,  0.002),
    "2Cseg":           (-0.095,  0.002),
    "NUC":             ( 0.010, -0.011),
    "Proseg":          ( 0.010, -0.011),
    "StarDist+lookup": (-0.010, -0.012),
    "StarDist+WBA":    ( 0.010,  0.002),
}

# Bioinformatics 標準：Arial/Helvetica，7–9 pt，去掉頂部和右側 spine
plt.rcParams.update({
    "font.family":      "sans-serif",
    "font.sans-serif":  ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size":        12,
    "axes.linewidth":   0.8,
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
    "xtick.major.size":  3,
    "ytick.major.size":  3,
    "pdf.fonttype":     42,   # embeds fonts as TrueType
    "ps.fonttype":      42,
})

# 單欄寬 8.5 cm ≈ 3.35 in；雙欄 17.6 cm ≈ 6.93 in
# 此圖設計為雙欄寬
fig, ax = plt.subplots(figsize=(6.93, 4.80))

for label, df in methods_data.items():
    x = df["ftc"].mean()
    y = df["ned"].mean()
    c = COLORS[label]
    m = MARKERS[label]
    # 主要方法較大，補充對照較小；MCseg 再大一些
    sz = 90 if label == "MCseg" else (70 if label in MAIN else 55)

    if label in HOLLOW:
        ax.scatter(x, y, s=sz, marker=m, facecolor="white",
                   edgecolors=c, linewidths=1.4, zorder=5, clip_on=False)
    else:
        ax.scatter(x, y, s=sz, marker=m, color=c,
                   edgecolors="white", linewidths=0.6, zorder=5, clip_on=False)

    dx, dy = LABEL_OFFSETS.get(label, (0.013, 0.002))
    weight = "bold" if label == "MCseg" else "normal"
    ax.text(x + dx, y + dy, label,
            fontsize=12, color=c, fontweight=weight, va="center")

# ── MCseg vs SR 統計標注：bracket 連線 ──
x_mc = means["MCseg"]["FTC"];  y_mc = means["MCseg"]["NED"]
x_sr = means["SR"]["FTC"];     y_sr = means["SR"]["NED"]

ned_row = stats_df[(stats_df["Comparison"]=="MCseg vs SR") &
                   (stats_df["Metric"]=="NED")].iloc[0]
p_val = ned_row["p_raw"]
stars = ned_row["sig_raw"]

# 在兩點之間畫一條淡線，中點放 p 值
mid_x = (x_mc + x_sr) / 2
mid_y = (y_mc + y_sr) / 2 + 0.012
ax.annotate("", xy=(x_sr, y_sr + 0.006), xytext=(x_mc, y_mc + 0.006),
            arrowprops=dict(arrowstyle="-", color="#555555",
                            lw=0.8, linestyle="dashed"))
ax.text(mid_x, mid_y, f"NED: p = {p_val:.3f} ({stars})",
        fontsize=11, color="#444444", ha="center", va="bottom")

# ── 象限說明（右下 & 左上），灰色斜體 ──
# 左上：靠左對齊並貼到頂，避免壓到 NUC 標記（0.245, 0.824）
ax.text(0.183, 0.847,
        "Low capture,\nhigh purity\n(nuclear methods)",
        fontsize=10.5, color="#999999", ha="left", style="italic", va="top",
        linespacing=1.4)
ax.text(0.935, 0.703,
        "High capture,\nlow purity\n(over-expansion)",
        fontsize=10.5, color="#999999", ha="center", style="italic", va="top",
        linespacing=1.4)

# ── 軸標籤（Bioinformatics 不用 title，改放 panel 字母）──
ax.set_xlabel("FTC (fraction of transcripts in cells) ↑", fontsize=13)
ax.set_ylabel("NED (normalized expression distance) ↑", fontsize=13)

# ── 圖例：主要 vs 補充 ──
from matplotlib.lines import Line2D
leg_elements = [
    Line2D([0], [0], marker="o", color="w", markerfacecolor="#555555",
           markeredgecolor="w", markersize=6, label="Primary methods (filled)"),
    Line2D([0], [0], marker="o", color="w", markerfacecolor="white",
           markeredgecolor="#555555", markeredgewidth=1.2,
           markersize=6, label="Supplementary (open)"),
]
# 圖例移到右上（該區為空白）：原本的 lower right 會蓋住
# "High capture, low purity" 象限說明。
ax.legend(handles=leg_elements, fontsize=11, frameon=True,
          framealpha=0.9, edgecolor="#cccccc", loc="upper right",
          handletextpad=0.4, borderpad=0.5)

ax.set_xlim(0.17, 1.06)
# 下緣留到 0.680：放大字級後右下角象限說明（3 行）需要更多垂直空間
ax.set_ylim(0.680, 0.848)
ax.tick_params(labelsize=11)

# 只保留左側和底部 spine
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

# 淡格線，不搶眼
ax.grid(True, axis="both", alpha=0.18, linestyle=":", linewidth=0.6, color="#000000")

plt.tight_layout(pad=0.5)
# 同時輸出到投稿圖檔目錄：此圖即 Figure 3 panel d
SUBMISSION_FIG3 = Path("/Volumes/SSD/plan_a/submission_bioinformatics/figures/fig3/fig3d.png")
SUBMISSION_FIG3.parent.mkdir(parents=True, exist_ok=True)
for fig_path in (FIG_DIR / "fig_pareto_capture_purity.png", SUBMISSION_FIG3):
    plt.savefig(fig_path, dpi=300, bbox_inches="tight")
    print(f"\n圖片已儲存：{fig_path}")
plt.close()
