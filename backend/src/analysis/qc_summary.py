"""QC 指標的直方圖摘要與範圍建議（QcRangeSuggestion，見 CONTEXT.md）。

與 `analysis/preprocessing.py::compute_qc_metrics` 分工：那邊算「每個細胞的
QC 指標值」，這裡把一組已知指標值摘要成直方圖，並從分布推導建議的篩選範圍，
供 Stage3 QC 門檻 UI 預填、使用者仍可手動覆寫。不依賴 FastAPI，可獨立測試。
"""
from dataclasses import asdict, dataclass


@dataclass
class QcMetricHistogram:
    """單一 QC 指標的直方圖 + 建議篩選範圍。"""

    label: str
    unit: str
    bin_edges: list
    counts: list
    mad_min: float
    mad_max: float
    p5: float
    p50: float
    p95: float
    p99: float
    mean: float

    def to_dict(self) -> dict:
        return asdict(self)


def compute_qc_histogram(arr, label: str, unit: str = "", n_bins: int = 60) -> "QcMetricHistogram | None":
    """把一組指標值摘要成直方圖，並推導建議篩選範圍（QcRangeSuggestion）。

    退化情況（濾完後空陣列、或全部為 0）回傳 `None`——全 0 時 `np.histogram`
    會產生 `[-0.5, 0.5]` 這種無意義的 bin_edges，寧可讓呼叫端自行決定如何處理，
    也不要回傳一份看似有效實則誤導的直方圖。
    """
    import numpy as np

    arr = arr[np.isfinite(arr) & (arr >= 0)]
    if len(arr) == 0:
        return None
    if arr.max() == 0:
        return None

    counts, bin_edges = np.histogram(arr, bins=n_bins)
    median = float(np.median(arr))
    # MAD 在 log1p 空間計算後轉回線性空間，適合 count data（高度右偏分布）。
    # 不乘 1.4826 比例因子，用原始 MAD 提供更直覺且緊湊的建議範圍，避免 exp
    # 轉回後上界被過度放大（Exponential Inflation）。
    log_arr = np.log1p(arr)
    log_median = float(np.median(log_arr))
    log_mad = float(np.median(np.abs(log_arr - log_median)))
    mad_min = float(max(0.0, np.expm1(log_median - 3 * log_mad)))
    mad_max = float(np.expm1(log_median + 3 * log_mad))

    return QcMetricHistogram(
        label=label,
        unit=unit,
        bin_edges=bin_edges.tolist(),
        counts=counts.tolist(),
        mad_min=round(mad_min, 2),
        mad_max=round(mad_max, 2),
        p5=round(float(np.percentile(arr, 5)), 2),
        p50=round(median, 2),
        p95=round(float(np.percentile(arr, 95)), 2),
        p99=round(float(np.percentile(arr, 99)), 2),
        mean=round(float(arr.mean()), 2),
    )
