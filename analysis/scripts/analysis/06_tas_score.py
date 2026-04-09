"""
05_tas_score.py
===============
讀取 metrics_ac.csv 和 metrics_de.csv，計算 TAS：

  Capture_score = mean(A1_norm, A2_norm, A3_norm)
  Purity_score  = mean(C1_norm, [C3_norm])   # C3 可選
  Core_TAS      = sqrt(Capture_score × Purity_score)
  TAS           = 0.65×Core_TAS + 0.20×D1_norm + 0.15×E1_norm

輸出：
  results/metrics/tas.csv（per-ROI）
  results/metrics/tas_summary.csv（per-method 平均）
"""

from __future__ import annotations

import sys
import numpy as np
import pandas as pd
from pathlib import Path

ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT))

import yaml
cfg = yaml.safe_load((ROOT / "config.yaml").read_text())
PATHS   = cfg["paths"]
DATA    = cfg["data"]
TAS_CFG = cfg["tas"]

METRICS_DIR = ROOT / PATHS["metrics_dir"]

UMI_REF  = TAS_CFG["umi_ref"]    # 1000
GENE_REF = TAS_CFG["genes_ref"]  # 300
W_CORE   = TAS_CFG["weight_core_tas"]   # 0.65
W_BIO    = TAS_CFG["weight_biology"]    # 0.20
W_IMM    = TAS_CFG["weight_immune"]     # 0.15
USE_C3   = TAS_CFG["use_c3_admixture"] # True / False


def normalize(df: pd.DataFrame) -> pd.DataFrame:
    """
    歸一化（1 = 最好）：

    A1_norm = A1（已是 0–1）
    A2_norm = clip(log(A2) / log(UMI_REF), 0, 1)
    A3_norm = clip(A3 / GENE_REF, 0, 1)
    C1_norm = 1 – C1
    C3_norm = 1 – C3_norm_raw（若有）
    D1_norm = (D1 + 1) / 2
    E1_norm = E1（已是 0–1）
    """
    out = df.copy()

    out["a1_norm"] = out["a1_capture"].clip(0, 1)

    # A1d: UMI density（max 歸一化，消除大細胞 FTC 偏差）
    if "a1_umi_density" in out.columns and out["a1_umi_density"].notna().any():
        max_density = out["a1_umi_density"].max()
        if max_density > 0:
            out["a1d_norm"] = (out["a1_umi_density"] / max_density).clip(0, 1)
        else:
            out["a1d_norm"] = np.nan
    else:
        out["a1d_norm"] = np.nan

    # A2: log normalization（UMI=0 → 0）
    out["a2_norm"] = np.where(
        out["a2_median_umi"] > 0,
        (np.log(out["a2_median_umi"].clip(lower=1)) / np.log(UMI_REF)).clip(0, 1),
        0.0
    )

    out["a3_norm"] = (out["a3_median_genes"] / GENE_REF).clip(0, 1)

    out["c1_norm"] = (1 - out["c1_coexpr"]).clip(0, 1)

    # NED_norm：Neighbor Expression Divergence（已是 [0,1]，1 = 最好）
    if "ned" in out.columns:
        out["ned_norm"] = out["ned"].clip(0, 1)
    else:
        out["ned_norm"] = np.nan

    out["d1_norm"] = ((out["d1_silhouette"] + 1) / 2).clip(0, 1)

    out["e1_norm"] = out["e1_immune_survival"].clip(0, 1)

    return out


def compute_tas(df: pd.DataFrame) -> pd.DataFrame:
    """計算 TAS 及各子分數。"""
    out = df.copy()

    # capture_score：優先使用 UMI density（消除大細胞偏差），若無則用 FTC
    if "a1d_norm" in out.columns and out["a1d_norm"].notna().any():
        a1_for_capture = out["a1d_norm"].fillna(out["a1_norm"])
        print("  [A1d] 使用 UMI density 取代 FTC 作為 capture score")
    else:
        a1_for_capture = out["a1_norm"]
        print("  [A1] 使用 FTC（UMI density 未計算）")
    # A1d 雙加權（2:1:1）：A1d 為唯一大小無關指標，A2/A3 受細胞大小影響
    out["capture_score"] = (
        2 * a1_for_capture + out["a2_norm"] + out["a3_norm"]
    ) / 4

    # Purity 組合：C3 > NED+C1 組合 > C1 alone（降級順序）
    if USE_C3 and "c3_norm" in out.columns and out["c3_norm"].notna().any():
        out["purity_score"] = out[["c1_norm", "c3_norm"]].mean(axis=1)
        print("  [C3] 使用 cellAdmix admixture score")
    elif "ned_norm" in out.columns and out["ned_norm"].notna().any():
        # NED + C1 幾何平均：同時考慮全域轉錄本相似度（NED）與基因對共現（C1）
        ned_vals = out["ned_norm"].fillna(out["c1_norm"])  # NED missing → fall back to C1
        out["purity_score"] = np.sqrt(out["c1_norm"] * ned_vals)
        print("  [NED+C1] 使用 Neighbor Expression Divergence + Artificial co-expression 組合")
    else:
        out["purity_score"] = out["c1_norm"]
        print("  [C1] 使用 Artificial co-expression only（NED 未計算）")

    # 幾何平均（防止 sqrt(0×...)）
    cs = out["capture_score"].values.clip(0, 1)
    ps = out["purity_score"].values.clip(0, 1)
    out["core_tas"] = np.sqrt(cs * ps)

    e1_contrib = W_IMM * out["e1_norm"] if W_IMM > 0 else 0.0
    out["tas"] = (
        W_CORE * out["core_tas"] +
        W_BIO  * out["d1_norm"] +
        e1_contrib
    ).clip(0, 1)

    return out


def run_tas():
    # 讀取 A/C 指標
    ac_path = METRICS_DIR / "metrics_ac.csv"
    de_path = METRICS_DIR / "metrics_de.csv"

    if not ac_path.exists():
        print(f"⚠️  找不到 {ac_path}，請先執行 03_metrics_ac.py")
        return
    if not de_path.exists():
        print(f"⚠️  找不到 {de_path}，請先執行 04_metrics_de.py")
        return

    df_ac = pd.read_csv(ac_path)
    df_de = pd.read_csv(de_path)

    # 合併（D1 是 per-method，E1 是 per-ROI）
    # df_de 有 method, roi, e1_immune_survival, d1_silhouette
    df = df_ac.merge(df_de, on=["method", "roi"], how="left")

    print(f"合併後：{len(df)} 筆記錄（{df['method'].nunique()} methods × {df['roi'].nunique()} ROIs）")

    # 歸一化
    df = normalize(df)

    # 計算 TAS
    df = compute_tas(df)

    # 輸出 per-ROI
    out_roi = METRICS_DIR / "tas.csv"
    df.to_csv(out_roi, index=False)

    # per-method 平均
    agg_cols = [
        "a1_capture", "a1_umi_density", "a2_median_umi", "a3_median_genes",
        "c1_coexpr", "ned", "d1_silhouette", "e1_immune_survival",
        "a1_norm", "a1d_norm", "a2_norm", "a3_norm", "c1_norm", "ned_norm",
        "capture_score", "purity_score", "core_tas",
        "d1_norm", "e1_norm", "tas"
    ]
    agg_cols = [c for c in agg_cols if c in df.columns]
    df_summary = df.groupby("method")[agg_cols].mean().round(4)

    out_summary = METRICS_DIR / "tas_summary.csv"
    df_summary.to_csv(out_summary)

    print(f"\n✅ 05_tas_score.py 完成")
    print(f"   per-ROI  → {out_roi}")
    print(f"   summary  → {out_summary}")
    print("\n=== TAS Summary (per-method mean) ===")
    display_cols = ["capture_score", "purity_score", "core_tas", "d1_norm", "e1_norm", "tas"]
    display_cols = [c for c in display_cols if c in df_summary.columns]
    print(df_summary[display_cols].sort_values("tas", ascending=False).to_string())


if __name__ == "__main__":
    run_tas()
