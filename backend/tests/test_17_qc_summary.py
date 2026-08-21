"""Test 17: QC 直方圖統計邏輯（`analysis/qc_summary.py`）

架構深化 P7：把 `api/analysis.py::_hist_metric` 的 MAD/percentile 數學搬到這裡，
使其可獨立單元測試，不需要透過 FastAPI route 才能驗證。
"""
import numpy as np
import pytest


class TestComputeQcHistogram:
    def test_normal_distribution_bounds_and_bins(self):
        """常態分布：mad_min < p50 < mad_max，bin_edges 長度為 n_bins+1"""
        from backend.src.analysis.qc_summary import compute_qc_histogram

        rng = np.random.default_rng(0)
        arr = rng.normal(loc=500, scale=50, size=1000).clip(min=0)

        result = compute_qc_histogram(arr, "Transcripts Per Cell", n_bins=60)

        assert result is not None
        assert len(result.bin_edges) == 61
        assert len(result.counts) == 60
        assert result.mad_min < result.p50 < result.mad_max

    def test_skewed_distribution_uses_log_space_mad(self):
        """右偏分布：log1p 空間 MAD 算出的上界要明顯大於線性空間直接算的 median+3*MAD

        驗證真的是在 log1p 空間做 MAD 轉換，而不是線性空間的 MAD（右偏分布下
        線性 MAD 會嚴重低估合理上界，這正是原本 _hist_metric 選擇 log 空間的理由）。
        """
        from backend.src.analysis.qc_summary import compute_qc_histogram

        rng = np.random.default_rng(1)
        arr = rng.lognormal(mean=5, sigma=1.5, size=2000)

        result = compute_qc_histogram(arr, "Skewed Metric", n_bins=60)

        assert result is not None
        linear_median = float(np.median(arr))
        linear_mad = float(np.median(np.abs(arr - linear_median)))
        linear_naive_max = linear_median + 3 * linear_mad
        assert result.mad_max > linear_naive_max

    def test_all_zero_array_returns_none(self):
        """全 0 陣列：np.histogram 會產生無效 bin_edges，必須提早回傳 None"""
        from backend.src.analysis.qc_summary import compute_qc_histogram

        result = compute_qc_histogram(np.zeros(50), "All Zero")
        assert result is None

    def test_empty_array_returns_none(self):
        from backend.src.analysis.qc_summary import compute_qc_histogram

        result = compute_qc_histogram(np.array([]), "Empty")
        assert result is None

    def test_single_value_array_mad_min_equals_mad_max(self):
        """單一值（MAD=0）：mad_min == mad_max == 該值"""
        from backend.src.analysis.qc_summary import compute_qc_histogram

        result = compute_qc_histogram(np.full(10, 5.0), "Constant")

        assert result is not None
        assert result.mad_min == pytest.approx(5.0, abs=0.01)
        assert result.mad_max == pytest.approx(5.0, abs=0.01)

    def test_to_dict_matches_expected_keys(self):
        """to_dict() 的 key 集合須與原本 _hist_metric 回傳的 dict 一致（JSON 契約）"""
        from backend.src.analysis.qc_summary import compute_qc_histogram

        rng = np.random.default_rng(2)
        result = compute_qc_histogram(rng.normal(100, 10, 200).clip(min=0), "X", unit="µm")

        assert result is not None
        d = result.to_dict()
        assert set(d.keys()) == {
            "label", "unit", "bin_edges", "counts",
            "mad_min", "mad_max", "p5", "p50", "p95", "p99", "mean",
        }
        assert d["label"] == "X"
        assert d["unit"] == "µm"
