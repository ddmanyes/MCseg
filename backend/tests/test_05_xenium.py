"""Test 5: Xenium Data — Xenium outs 資料讀取驗證"""
import pytest
from pathlib import Path

from .conftest import CRC_XENIUM


class TestXeniumTranscripts:
    """Xenium transcripts.parquet 讀取"""

    def test_read_transcripts_schema(self):
        """transcripts.parquet 包含必要欄位"""
        import pandas as pd
        tx_path = CRC_XENIUM / "transcripts.parquet"
        # Only read first 100 rows to verify schema
        df = pd.read_parquet(tx_path, columns=None)
        cols = set(df.columns)
        # Standard Xenium transcript columns
        for c in ["x_location", "y_location", "feature_name"]:
            assert c in cols, f"Missing column: {c}"
        print(f"  Transcripts: {len(df):,} rows, {len(cols)} columns")

    def test_transcript_coords_range(self):
        """transcript 座標在合理範圍內"""
        import pandas as pd
        tx_path = CRC_XENIUM / "transcripts.parquet"
        df = pd.read_parquet(tx_path, columns=["x_location", "y_location"])
        assert df["x_location"].min() >= 0
        assert df["y_location"].min() >= 0
        # CRC Xenium: x: 6~10688, y: 6~4352 µm (from rois.yaml)
        assert df["x_location"].max() < 20000  # safety margin
        assert df["y_location"].max() < 20000


class TestXeniumCells:
    """Xenium cells.parquet 讀取"""

    def test_read_cells_schema(self):
        """cells.parquet 包含 cell_id 和座標"""
        import pandas as pd
        cells_path = CRC_XENIUM / "cells.parquet"
        df = pd.read_parquet(cells_path)
        assert "cell_id" in df.columns
        assert "x_centroid" in df.columns
        assert "y_centroid" in df.columns
        print(f"  Cells: {len(df):,}")


class TestXeniumExperiment:
    """experiment.xenium JSON"""

    def test_read_experiment(self):
        """experiment.xenium 是有效 JSON"""
        import json
        exp_path = CRC_XENIUM / "experiment.xenium"
        with open(exp_path) as f:
            data = json.load(f)
        assert isinstance(data, dict)
        # Should have pixel_size or similar field
        assert "pixel_size" in data or "analysis_sw_version" in data


class TestXeniumNucleusBoundaries:
    """nucleus_boundaries.parquet"""

    def test_read_nucleus_boundaries(self):
        """nucleus_boundaries.parquet 可讀取"""
        import pandas as pd
        nb_path = CRC_XENIUM / "nucleus_boundaries.parquet"
        df = pd.read_parquet(nb_path)
        assert len(df) > 0
        assert "cell_id" in df.columns
        print(f"  Nucleus boundaries: {len(df):,} rows")
