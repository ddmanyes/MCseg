"""
08_stats.py
===========
非參數統計顯著性分析：
  - Friedman 檢定（4 methods × 7 ROIs repeated measures）
  - Pairwise Wilcoxon signed-rank 檢定 + Bonferroni 校正
  - 針對 TAS 及各子指標

輸出：
  results/metrics/stats_friedman.csv   — 每指標 Friedman 結果
  results/metrics/stats_pairwise.csv   — 每指標 × 每方法對 pairwise 結果
  results/metrics/stats_summary.txt    — 人類可讀摘要
"""

from __future__ import annotations

import sys
import itertools
import numpy as np
import pandas as pd
from pathlib import Path
from scipy import stats

ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT))

import yaml
cfg = yaml.safe_load((ROOT / "config.yaml").read_text())
PATHS = cfg["paths"]

METRICS_DIR = ROOT / PATHS["metrics_dir"]


def friedman_test(df: pd.DataFrame, metric: str) -> dict:
    """
    Friedman 檢定：4 methods × 7 ROIs（ROI 為 block）。
    """
    pivot = df.pivot(index="roi", columns="method", values=metric).dropna()
    if pivot.shape[1] < 3 or pivot.shape[0] < 3:
        return {"metric": metric, "statistic": np.nan, "p_value": np.nan, "n_blocks": 0}

    groups = [pivot[col].values for col in pivot.columns]
    stat, p = stats.friedmanchisquare(*groups)
    return {
        "metric": metric,
        "statistic": round(float(stat), 4),
        "p_value": float(p),
        "n_blocks": int(pivot.shape[0]),
        "methods": list(pivot.columns),
    }


def pairwise_wilcoxon(df: pd.DataFrame, metric: str) -> list[dict]:
    """
    Pairwise Wilcoxon signed-rank（配對，ROI 為 pair）+ Bonferroni 校正。
    """
    pivot = df.pivot(index="roi", columns="method", values=metric).dropna()
    methods = list(pivot.columns)
    n_pairs = len(list(itertools.combinations(methods, 2)))
    rows = []
    for m1, m2 in itertools.combinations(methods, 2):
        x = pivot[m1].values
        y = pivot[m2].values
        try:
            stat, p_raw = stats.wilcoxon(x, y, alternative="two-sided")
        except ValueError:
            stat, p_raw = np.nan, np.nan
        p_bonf = min(float(p_raw) * n_pairs, 1.0) if not np.isnan(p_raw) else np.nan
        rows.append({
            "metric": metric,
            "method_a": m1,
            "method_b": m2,
            "wilcoxon_stat": round(float(stat), 4) if not np.isnan(stat) else np.nan,
            "p_raw": float(p_raw),
            "p_bonferroni": p_bonf,
            "significant_0.05": p_bonf < 0.05 if not np.isnan(p_bonf) else False,
            "mean_diff": round(float(pivot[m1].mean() - pivot[m2].mean()), 4),
        })
    return rows


def run_stats():
    tas_path = METRICS_DIR / "tas.csv"
    if not tas_path.exists():
        print(f"⚠️  找不到 {tas_path}，請先執行 05_tas_score.py")
        return

    df = pd.read_csv(tas_path)

    # 分析的指標
    metrics = [
        "tas", "core_tas", "capture_score", "purity_score",
        "a1_capture", "a1_umi_density", "a2_median_umi",
        "a3_median_genes", "c1_coexpr", "ned",
        "d1_norm", "e1_norm",
    ]
    metrics = [m for m in metrics if m in df.columns]

    # ── Friedman 檢定 ─────────────────────────────────────
    friedman_rows = []
    for metric in metrics:
        result = friedman_test(df, metric)
        friedman_rows.append(result)
        p = result["p_value"]
        sig = "***" if p < 0.001 else ("**" if p < 0.01 else ("*" if p < 0.05 else "ns"))
        print(f"  Friedman {metric:25s}: χ²={result['statistic']:.3f}  p={p:.4f}  {sig}")

    df_friedman = pd.DataFrame(friedman_rows)[["metric", "statistic", "p_value", "n_blocks"]]
    df_friedman.to_csv(METRICS_DIR / "stats_friedman.csv", index=False)

    # ── Pairwise Wilcoxon ─────────────────────────────────
    pairwise_rows = []
    for metric in metrics:
        pairwise_rows.extend(pairwise_wilcoxon(df, metric))

    df_pw = pd.DataFrame(pairwise_rows)
    df_pw.to_csv(METRICS_DIR / "stats_pairwise.csv", index=False)

    # ── 人類可讀摘要 ──────────────────────────────────────
    summary_lines = [
        "統計顯著性摘要（Friedman + Pairwise Wilcoxon）",
        "=" * 60,
        f"設計：4 methods × 7 ROIs repeated-measures",
        f"多重比較校正：Bonferroni（n_pairs=6 per metric）",
        "",
        "── Friedman 檢定 ────────────────────────────────────────",
    ]
    for _, row in df_friedman.iterrows():
        p = row["p_value"]
        sig = "***" if p < 0.001 else ("**" if p < 0.01 else ("*" if p < 0.05 else "ns"))
        summary_lines.append(
            f"  {row['metric']:25s}: χ²={row['statistic']:.3f}  p={p:.4f}  {sig}"
        )

    summary_lines += ["", "── Pairwise Wilcoxon（TAS, Bonferroni 校正）──────────────"]
    tas_pw = df_pw[df_pw["metric"] == "tas"].sort_values("p_bonferroni")
    for _, row in tas_pw.iterrows():
        sig = "*" if row["significant_0.05"] else "ns"
        summary_lines.append(
            f"  {row['method_a']:6s} vs {row['method_b']:6s}: "
            f"Δmean={row['mean_diff']:+.4f}  p_bonf={row['p_bonferroni']:.4f}  {sig}"
        )

    summary_lines += ["", "── Pairwise Wilcoxon（所有顯著對）──────────────────────"]
    sig_pw = df_pw[df_pw["significant_0.05"]].sort_values(["metric", "p_bonferroni"])
    if len(sig_pw) == 0:
        summary_lines.append("  （無顯著差異對）")
    for _, row in sig_pw.iterrows():
        summary_lines.append(
            f"  [{row['metric']:20s}] {row['method_a']:6s} vs {row['method_b']:6s}: "
            f"Δ={row['mean_diff']:+.4f}  p_bonf={row['p_bonferroni']:.4f}  *"
        )

    summary_txt = "\n".join(summary_lines)
    out_txt = METRICS_DIR / "stats_summary.txt"
    out_txt.write_text(summary_txt, encoding="utf-8")

    print(f"\n✅ 08_stats.py 完成")
    print(f"   {METRICS_DIR / 'stats_friedman.csv'}")
    print(f"   {METRICS_DIR / 'stats_pairwise.csv'}")
    print(f"   {out_txt}")
    print("\n" + summary_txt)


if __name__ == "__main__":
    run_stats()
