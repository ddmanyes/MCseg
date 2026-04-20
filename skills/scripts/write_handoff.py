"""
分割品質交接報告產生器
讀取 qc_metrics.csv，計算摘要統計，依規則推薦分析參數，
寫入 results/handoff_report.json。
"""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import yaml

ROOT = Path(__file__).parent.parent


def _quality_label(ned_delta: float) -> str:
    if ned_delta >= 0.03:
        return "good"
    if ned_delta >= 0.01:
        return "marginal"
    return "poor"


def write_handoff(
    metrics_csv:    Path | str | None = None,
    masks_dir:      str | None = None,
    binned_dir:     str | None = None,
    tissue_profile: str | None = None,
    out_path:       Path | str | None = None,
) -> dict:
    if metrics_csv is None:
        metrics_csv = ROOT / "results" / "qc_metrics.csv"
    if out_path is None:
        out_path = ROOT / "results" / "handoff_report.json"

    cfg_raw = yaml.safe_load((ROOT / "config" / "pipeline.yaml").read_text(encoding="utf-8"))
    if masks_dir is None:
        masks_dir = cfg_raw["paths"].get("masks_dir", "results/masks")
    if binned_dir is None:
        binned_dir = cfg_raw["paths"]["binned_002"]
    if tissue_profile is None:
        tissue_profile = cfg_raw["global"].get("tissue_profile", "crc")

    df = pd.read_csv(metrics_csv)
    mcseg = df[df["method"] == "mcseg"]
    nuc   = df[df["method"] == "nuc"]

    ned_mcseg  = float(mcseg["ned"].mean())
    ned_nuc    = float(nuc["ned"].mean()) if len(nuc) > 0 else 0.0
    ned_delta  = ned_mcseg - ned_nuc
    ftc_mean   = float(mcseg["ftc"].mean())
    coexp_mean = float(mcseg["coexp_rate"].mean())

    warnings = []
    if _quality_label(ned_delta) == "poor":
        warnings.append("NED 提升 < 0.01：建議調整 dia_mid (±2px) 或 clahe_clip_limit")
    if coexp_mean > 0.06:
        warnings.append("Co-expression rate > 6%：邊界過度擴張，建議減少 voronoi_distance")
    if ftc_mean < 0.60:
        warnings.append("FTC < 60%：建議增加 voronoi_distance 或 expand_labels 距離")

    # 建議 QC 閾值（保守估計，分析階段可再調整）
    n_cells_mean = float(mcseg["n_cells"].mean()) if len(mcseg) > 0 else 0
    recommended = {
        "min_genes":    150,
        "max_pct_mt":   12.0,
        "min_counts":   80,
        "leiden_resolution": (
            0.8 if n_cells_mean < 5000 else
            0.5 if n_cells_mean < 50000 else
            0.3
        ),
    }

    report = {
        "segmentation_complete": True,
        "n_rois_evaluated":      int(len(mcseg)),
        "roi_qc": {
            "ned_mcseg":  round(ned_mcseg, 4),
            "ned_nuc":    round(ned_nuc, 4),
            "ned_delta":  round(ned_delta, 4),
            "quality":    _quality_label(ned_delta),
            "ftc_mean":   round(ftc_mean, 4),
            "coexp_mean": round(coexp_mean, 4),
        },
        "warnings":                  warnings,
        "recommended_analysis_params": recommended,
        "masks_dir":      str(masks_dir),
        "binned_dir":     str(binned_dir),
        "tissue_profile": tissue_profile,
    }

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")

    print(f"✅ handoff_report.json 已寫入 {out_path}")
    print(f"   NED delta={ned_delta:+.4f} ({report['roi_qc']['quality']})")
    if warnings:
        for w in warnings:
            print(f"   ⚠️  {w}")

    return report


if __name__ == "__main__":
    write_handoff()
