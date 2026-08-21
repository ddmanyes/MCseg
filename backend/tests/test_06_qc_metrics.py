"""Test 6: QC 指標唯一擁有者 — compute_qc_metrics（不需真實資料）"""
import numpy as np
import pytest


@pytest.fixture
def toy_adata():
    """3 細胞 x 3 基因；第 3 顆為空細胞（total_counts=0）"""
    import anndata as ad
    import pandas as pd

    X = np.array([
        [10.0, 0.0, 5.0],    # 2 genes, 15 counts, 5 mito
        [0.0, 4.0, 0.0],     # 1 gene,  4 counts,  0 mito
        [0.0, 0.0, 0.0],     # 空細胞
    ])
    var = pd.DataFrame(index=["GENE1", "GENE2", "MT-CO1"])
    return ad.AnnData(X=X, var=var)


class TestComputeQCMetrics:
    """指標定義只有一份，且 API / pipeline 都穿過它"""

    def test_mito_prefix_from_params(self, toy_adata):
        """mito 前綴取自 QC 參數（大小寫不敏感）"""
        from backend.src.analysis.preprocessing import compute_qc_metrics

        adata = compute_qc_metrics(toy_adata, {"mito_prefix": "MT-"})
        assert adata.var["mt"].tolist() == [False, False, True]
        assert adata.obs["pct_counts_mt"].iloc[0] == pytest.approx(100 * 5 / 15)

    def test_empty_cell_mito_is_zero_not_nan(self, toy_adata):
        """空細胞的 pct_counts_mt 填 0，不留 NaN（NaN 會污染直方圖與過濾）"""
        from backend.src.analysis.preprocessing import compute_qc_metrics

        adata = compute_qc_metrics(toy_adata, {"mito_prefix": "MT-"})
        assert adata.obs["pct_counts_mt"].isna().sum() == 0
        assert adata.obs["pct_counts_mt"].iloc[2] == 0.0

    def test_complexity_finite(self, toy_adata):
        """complexity = log10(genes)/log10(counts)，極端值收斂為 0"""
        from backend.src.analysis.preprocessing import compute_qc_metrics

        adata = compute_qc_metrics(toy_adata, {"mito_prefix": "MT-"})
        complexity = adata.obs["complexity"].values
        assert np.isfinite(complexity).all()
        assert complexity[0] == pytest.approx(np.log10(2) / np.log10(15))
        assert complexity[2] == 0.0

    def test_preprocessor_delegates(self, toy_adata):
        """Preprocessor.calculate_qc_metrics 與模組函式結果一致"""
        from backend.src.analysis.preprocessing import Preprocessor, compute_qc_metrics

        pre = Preprocessor({"preprocessing": {"cellular": {"mito_prefix": "MT-"}}})
        via_method = pre.calculate_qc_metrics(toy_adata.copy())
        via_func = compute_qc_metrics(toy_adata.copy(), {"mito_prefix": "MT-"})

        assert via_method.var["mt"].tolist() == via_func.var["mt"].tolist()
        for col in ("total_counts", "n_genes_by_counts", "pct_counts_mt", "complexity"):
            np.testing.assert_allclose(via_method.obs[col].values, via_func.obs[col].values)

    def test_api_histogram_uses_shared_entry(self):
        """raw_histogram 不得自行重算 QC 指標（避免定義分歧）"""
        from pathlib import Path

        src = Path(__file__).resolve().parents[1] / "src" / "api" / "analysis.py"
        body = src.read_text(encoding="utf-8")
        start = body.index("async def get_raw_histogram")
        end = body.index("async def get_qc_status")
        histogram_src = body[start:end]

        assert "compute_qc_metrics" in histogram_src
        assert "sc.pp.calculate_qc_metrics" not in histogram_src
        assert "np.log10" not in histogram_src

    def test_api_histogram_delegates_to_qc_summary(self):
        """raw_histogram 的直方圖/MAD 建議範圍統計，必須委派給 analysis.qc_summary

        （不得自己在 API route 裡重算 log1p 空間 MAD——那段數學屬於領域層，
        搬到 qc_summary.py 後才能獨立單元測試，見架構深化 P7）
        """
        from pathlib import Path

        src = Path(__file__).resolve().parents[1] / "src" / "api" / "analysis.py"
        body = src.read_text(encoding="utf-8")
        start = body.index("async def get_raw_histogram")
        end = body.index("async def get_qc_status")
        histogram_src = body[start:end]

        assert "compute_qc_histogram" in histogram_src
