"""
09_bootstrap_permutation.py
============================
對 CRC 15 ROI 的 NED 與共表達率（Doublet_rate）進行：
  1. Bootstrap 95% CI（MCseg - SR 差值，10,000 次重抽樣）
  2. ROI-level permutation test（打亂方法標籤，10,000 次）

資料來源：
  crc_transcript_attribution/results/metrics/SuppTable1_crc_metrics_per_roi.csv
  （方法欄位：sr / v12 / p3 / nuc）

輸出：
  submission_bioinformatics/results/stats_bootstrap_permutation.csv
  submission_bioinformatics/results/stats_bootstrap_permutation.txt
"""

from __future__ import annotations
import numpy as np
import pandas as pd
from pathlib import Path

RNG_SEED = 42
N_BOOT   = 10_000
N_PERM   = 10_000

DATA_PATH = Path("/Volumes/SSD/plan_a/crc_transcript_attribution/results/metrics/SuppTable1_crc_metrics_per_roi.csv")
OUT_CSV   = Path("/Volumes/SSD/plan_a/submission_bioinformatics/results/stats_bootstrap_permutation.csv")
OUT_TXT   = Path("/Volumes/SSD/plan_a/submission_bioinformatics/results/stats_bootstrap_permutation.txt")

# (metric, method_a, method_b, higher_is_better_for_a)
COMPARISONS = [
    ("NED",          "v12", "sr",  True),   # MCseg NED > SR NED
    ("Doublet_rate", "v12", "sr",  False),  # MCseg doublet < SR doublet
    ("NED",          "v12", "p3",  True),   # MCseg NED vs Proseg
]

METHOD_LABEL = {"v12": "MCseg", "sr": "SR", "p3": "Proseg", "nuc": "NUC"}


def bootstrap_ci(diff: np.ndarray, n_boot: int, rng: np.random.Generator):
    boot_means = np.array([
        rng.choice(diff, size=len(diff), replace=True).mean()
        for _ in range(n_boot)
    ])
    lo, hi = np.percentile(boot_means, [2.5, 97.5])
    return float(lo), float(hi), float(boot_means.mean())


def permutation_test(a: np.ndarray, b: np.ndarray, n_perm: int,
                     rng: np.random.Generator, higher_is_better: bool):
    """
    Sign-flip permutation test on paired differences.
    H0: sign of each ROI's (a - b) is exchangeable.
    """
    diff     = a - b
    obs_stat = diff.mean()
    signs    = rng.choice([-1, 1], size=(n_perm, len(diff)))
    perm_stats = (signs * diff).mean(axis=1)

    if higher_is_better:
        p_val = float((perm_stats >= obs_stat).mean())
    else:
        p_val = float((perm_stats <= obs_stat).mean())
    return p_val, float(obs_stat)


def main():
    rng = np.random.default_rng(RNG_SEED)
    df  = pd.read_csv(DATA_PATH)

    records = []
    lines   = []

    for metric, ma, mb, hib in COMPARISONS:
        la, lb = METHOD_LABEL[ma], METHOD_LABEL[mb]
        sub    = df[df["method"].isin([ma, mb])].copy()
        pivot  = sub.pivot(index="roi", columns="method", values=metric).dropna()

        if ma not in pivot.columns or mb not in pivot.columns:
            print(f"[WARN] Missing data for {metric} {ma} vs {mb}")
            continue

        a    = pivot[ma].values
        b    = pivot[mb].values
        diff = a - b
        n    = len(diff)

        ci_lo, ci_hi, boot_mean = bootstrap_ci(diff, N_BOOT, rng)
        p_perm, obs_mean        = permutation_test(a, b, N_PERM, rng, hib)

        ci_excludes_zero = (ci_lo > 0) if hib else (ci_hi < 0)

        records.append({
            "metric":           metric,
            "comparison":       f"{la} vs {lb}",
            "n_rois":           n,
            "obs_mean_diff":    round(obs_mean, 5),
            "boot_mean_diff":   round(boot_mean, 5),
            "ci_95_lo":         round(ci_lo, 5),
            "ci_95_hi":         round(ci_hi, 5),
            "ci_excludes_zero": ci_excludes_zero,
            "perm_p":           round(p_perm, 4),
            "direction":        "higher_better" if hib else "lower_better",
        })

        direction_str = f"{la} > {lb}" if hib else f"{la} < {lb}"
        sig_str = "✓ significant" if p_perm < 0.05 else "✗ not significant"
        ci_str  = "excludes 0" if ci_excludes_zero else "includes 0"

        block = (
            f"\n{'='*60}\n"
            f"  {metric}  |  {direction_str}  |  n = {n} ROIs\n"
            f"{'='*60}\n"
            f"  Observed mean diff : {obs_mean:+.5f}\n"
            f"  Bootstrap 95% CI   : [{ci_lo:+.5f}, {ci_hi:+.5f}]  ({ci_str})\n"
            f"  Permutation p-value: {p_perm:.4f}  ({sig_str})\n"
        )
        lines.append(block)
        print(block)

    pd.DataFrame(records).to_csv(OUT_CSV, index=False)

    summary = "".join(lines)
    summary += (
        f"\n\nParameters: n_boot={N_BOOT}, n_perm={N_PERM}, seed={RNG_SEED}\n"
        f"Data: {DATA_PATH}\n"
    )
    OUT_TXT.write_text(summary)
    print(f"\n✓ {OUT_CSV}")
    print(f"✓ {OUT_TXT}")


if __name__ == "__main__":
    main()
