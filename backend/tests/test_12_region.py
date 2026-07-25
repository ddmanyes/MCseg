"""Test 12: 全片座標系上的區域／多邊形細胞過濾（backend/src/api/spatial.py）

全部合成資料。重點在**座標系**：RegionSelector 在全片 fullres px 上框選，
而細胞的座標可能來自三種不同來源（全圖流程、ROI 流程、舊資料），
必須全部落回同一個座標系才不會選到錯誤區域。
"""
import numpy as np
import pytest


def _make_adata(coords_fullres, *, kind="fullres", roi_x=0, roi_y=0, pixel_size_um=0.2737):
    """造一份 AnnData，把已知的全片座標以指定來源型態放進去。

    kind:
      - "fullres" : obs['centroid_x_fullres'] / _y_fullres（P0-7 全圖流程）
      - "local"   : obs['centroid_x_px'] / _y_px（ROI 流程，需加回 ROI 原點）
      - "spatial" : 只有 obsm['spatial']（µm，舊資料）
    """
    import anndata as ad

    coords = np.asarray(coords_fullres, dtype=float)
    n = len(coords)
    adata = ad.AnnData(X=np.ones((n, 2), dtype=np.float32))
    adata.obs_names = [f"cell_{i}" for i in range(n)]
    adata.var_names = ["GENE_A", "GENE_B"]

    if kind == "fullres":
        adata.obs["centroid_x_fullres"] = coords[:, 0]
        adata.obs["centroid_y_fullres"] = coords[:, 1]
        adata.obsm["spatial"] = coords * pixel_size_um
    elif kind == "local":
        adata.obs["centroid_x_px"] = coords[:, 0] - roi_x
        adata.obs["centroid_y_px"] = coords[:, 1] - roi_y
        adata.obsm["spatial"] = (coords - [roi_x, roi_y]) * pixel_size_um
    elif kind == "spatial":
        adata.obsm["spatial"] = (coords - [roi_x, roi_y]) * pixel_size_um
    else:                                   # pragma: no cover - 測試自身錯誤
        raise ValueError(kind)
    return adata


GRID = [(x, y) for y in (100, 200, 300) for x in (100, 200, 300)]   # 9 個細胞
BOX = {"x0": 150, "y0": 150, "x1": 350, "y1": 250}                  # 涵蓋 y=200 那列的 x=200,300


# ── bbox 過濾 ───────────────────────────────────────────────────────────────

class TestFilterByRegion:
    """矩形框選"""

    def test_filter_by_region_bbox(self):
        from backend.src.api.spatial import _filter_by_region

        adata = _make_adata(GRID)
        mask = _filter_by_region(adata, BOX, None, {}, None)

        assert mask.sum() == 2
        assert adata.obs_names[mask].tolist() == ["cell_4", "cell_5"]

    def test_boundary_points_are_included(self):
        """正好落在邊上的細胞須計入 —— 否則使用者框住的細胞會莫名消失。"""
        from backend.src.api.spatial import _filter_by_region

        adata = _make_adata([(100, 100), (200, 200)])
        mask = _filter_by_region(adata, {"x0": 100, "y0": 100, "x1": 150, "y1": 150}, None, {}, None)

        assert mask.tolist() == [True, False]

    def test_region_accepts_reversed_corners(self):
        """使用者從右下往左上拉框時，x0 > x1 —— 須自動正規化而非回傳空集合。"""
        from backend.src.api.spatial import _filter_by_region

        adata = _make_adata(GRID)
        rev = {"x0": 350, "y0": 250, "x1": 150, "y1": 150}

        assert _filter_by_region(adata, rev, None, {}, None).sum() == 2


# ── polygon 過濾 ────────────────────────────────────────────────────────────

class TestFilterByPolygon:
    """多邊形框選"""

    def test_filter_by_polygon(self):
        from backend.src.api.spatial import _filter_by_region

        adata = _make_adata(GRID)
        triangle = [[90, 90], [310, 90], [90, 310]]     # 左上直角三角形
        mask = _filter_by_region(adata, None, triangle, {}, None)

        kept = adata.obs_names[mask].tolist()
        # (100,100) (200,100) (300,100) (100,200) (200,200) (100,300) 在三角形內
        assert kept == ["cell_0", "cell_1", "cell_2", "cell_3", "cell_4", "cell_6"]

    def test_polygon_wins_over_region(self):
        """兩者同時給定時以 polygon 為準（較精確的那個）。"""
        from backend.src.api.spatial import _filter_by_region

        adata = _make_adata(GRID)
        tiny = [[95, 95], [105, 95], [105, 105], [95, 105]]
        mask = _filter_by_region(adata, BOX, tiny, {}, None)

        assert adata.obs_names[mask].tolist() == ["cell_0"]

    def test_degenerate_polygon_selects_nothing(self):
        """點數不足的多邊形視為未框選任何細胞，不可拋錯。"""
        from backend.src.api.spatial import _filter_by_region

        adata = _make_adata(GRID)
        assert _filter_by_region(adata, None, [[1, 1], [2, 2]], {}, None).sum() == 0


# ── 三層座標系回退 ───────────────────────────────────────────────────────────

class TestCoordinateFallback:
    """框選座標一律是全片 fullres px；三種資料來源都必須落回同一系統"""

    def test_uses_fullres_columns_when_present(self):
        from backend.src.api.spatial import _filter_by_region

        adata = _make_adata(GRID, kind="fullres")
        assert _filter_by_region(adata, BOX, None, {}, None).sum() == 2

    def test_local_centroid_gets_roi_origin_added(self):
        """ROI 流程的 centroid 是**裁切局部**座標，必須加回 ROI 原點。

        不加就會整批偏移一整個裁切原點（真實 ROI 可達數萬 px）——
        框選到的會是完全錯誤的區域，而且不會報錯。
        """
        from backend.src.api.spatial import _filter_by_region

        config = {"rois": [{"name": "r1", "x": 48128, "y": 12657, "pixel_size_um": 0.2737}]}
        shifted = [(x + 48128, y + 12657) for x, y in GRID]
        box = {"x0": 150 + 48128, "y0": 150 + 12657, "x1": 350 + 48128, "y1": 250 + 12657}

        adata = _make_adata(shifted, kind="local", roi_x=48128, roi_y=12657)
        mask = _filter_by_region(adata, box, None, config, "r1")

        assert adata.obs_names[mask].tolist() == ["cell_4", "cell_5"]

    def test_spatial_um_fallback_converts_and_offsets(self):
        """只有 obsm['spatial']（µm）的舊資料：先 ÷ pixel_size 再加 ROI 原點。"""
        from backend.src.api.spatial import _filter_by_region

        config = {"rois": [{"name": "r1", "x": 48128, "y": 12657, "pixel_size_um": 0.2737}]}
        shifted = [(x + 48128, y + 12657) for x, y in GRID]
        box = {"x0": 150 + 48128, "y0": 150 + 12657, "x1": 350 + 48128, "y1": 250 + 12657}

        adata = _make_adata(shifted, kind="spatial", roi_x=48128, roi_y=12657)
        mask = _filter_by_region(adata, box, None, config, "r1")

        assert adata.obs_names[mask].tolist() == ["cell_4", "cell_5"]

    def test_local_without_roi_config_stays_local(self):
        """找不到 ROI 設定時退為局部座標，並如實反映（不可假裝成功）。"""
        from backend.src.api.spatial import _filter_by_region

        adata = _make_adata(GRID, kind="local", roi_x=0, roi_y=0)
        assert _filter_by_region(adata, BOX, None, {}, "unknown").sum() == 2


# ── 未框選 ──────────────────────────────────────────────────────────────────

def test_no_region_selects_all():
    from backend.src.api.spatial import _filter_by_region

    adata = _make_adata(GRID)
    assert _filter_by_region(adata, None, None, {}, None).all()


# ── API：gene_plot 區域欄位與 region_stats ─────────────────────────────────

def _write_h5ad(tmp_path, adata, clusters=None):
    """把合成 AnnData 寫成 cellpose_cells.h5ad，並回傳可用的 config。"""
    if clusters is not None:
        adata.obs["leiden"] = [str(c) for c in clusters]
    out = tmp_path / "roi" / "r1"
    out.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(str(out / "cellpose_cells.h5ad"))
    return {"paths": {"output_dir": str(tmp_path)}, "rois": [{"name": "r1"}], "analysis": {}}


class TestRegionAPI:
    """/api/spatial/region_stats 與 gene_plot 的 region 欄位"""

    @pytest.mark.asyncio
    async def test_gene_plot_rejects_tiny_region(self, monkeypatch, tmp_path):
        """框到不足 10 個細胞時回 400 —— 少數細胞的圖只會誤導。"""
        from httpx import ASGITransport, AsyncClient

        import backend.src.api.spatial as sp
        from backend.main import app

        adata = _make_adata(GRID)
        config = _write_h5ad(tmp_path, adata)
        monkeypatch.setattr(sp, "load_config", lambda: config)
        monkeypatch.setattr(sp, "_get_roi_h5ad", lambda c, n: tmp_path / "roi/r1/cellpose_cells.h5ad")

        transport = ASGITransport(app=app)
        async with AsyncClient(transport=transport, base_url="http://test") as c:
            r = await c.post("/api/spatial/gene_plot", json={
                "roi_name": "r1", "genes": ["GENE_A"],
                "region": {"x0": 90, "y0": 90, "x1": 110, "y1": 110},
            })

        assert r.status_code == 400
        assert "細胞數不足" in r.json()["detail"]

    @pytest.mark.asyncio
    async def test_region_stats_returns_cluster_breakdown(self, monkeypatch, tmp_path):
        from httpx import ASGITransport, AsyncClient

        import backend.src.api.spatial as sp
        from backend.main import app

        adata = _make_adata(GRID)
        config = _write_h5ad(tmp_path, adata, clusters=[0, 0, 0, 1, 1, 1, 2, 2, 2])
        monkeypatch.setattr(sp, "load_config", lambda: config)
        monkeypatch.setattr(sp, "_get_roi_h5ad", lambda c, n: tmp_path / "roi/r1/cellpose_cells.h5ad")

        transport = ASGITransport(app=app)
        async with AsyncClient(transport=transport, base_url="http://test") as c:
            r = await c.post("/api/spatial/region_stats", json={
                "roi_name": "r1",
                "region": {"x0": 150, "y0": 150, "x1": 350, "y1": 250},
            })

        data = r.json()["data"]
        assert data["n_cells"] == 2
        assert data["cluster_counts"] == {"1": 2}
        assert data["median_counts"] == pytest.approx(2.0)   # 每個細胞 2 個基因各 1 count
        assert data["median_genes"] == pytest.approx(2.0)

    @pytest.mark.asyncio
    async def test_region_stats_empty_selection_is_ok(self, monkeypatch, tmp_path):
        """框到空白區不是錯誤 —— 前端據此提示「此區無細胞」。"""
        from httpx import ASGITransport, AsyncClient

        import backend.src.api.spatial as sp
        from backend.main import app

        adata = _make_adata(GRID)
        config = _write_h5ad(tmp_path, adata)
        monkeypatch.setattr(sp, "load_config", lambda: config)
        monkeypatch.setattr(sp, "_get_roi_h5ad", lambda c, n: tmp_path / "roi/r1/cellpose_cells.h5ad")

        transport = ASGITransport(app=app)
        async with AsyncClient(transport=transport, base_url="http://test") as c:
            r = await c.post("/api/spatial/region_stats", json={
                "region": {"x0": 9000, "y0": 9000, "x1": 9100, "y1": 9100},
            })

        assert r.json()["data"] == {
            "n_cells": 0, "median_counts": 0.0, "median_genes": 0.0,
            "cluster_counts": {}, "roi_name": None,
        }
