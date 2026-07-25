"""Test 3: API Endpoints — FastAPI 端點測試"""
import pytest
import pytest_asyncio
from httpx import AsyncClient, ASGITransport
from backend.main import app

pytestmark = pytest.mark.asyncio


@pytest_asyncio.fixture
async def client():
    """Async test client using httpx"""
    transport = ASGITransport(app=app)
    async with AsyncClient(transport=transport, base_url="http://test") as c:
        yield c


class TestHealthEndpoints:
    """基礎端點"""

    async def test_health(self, client):
        r = await client.get("/api/health")
        assert r.status_code == 200
        data = r.json()
        assert data["status"] == "ok"

    async def test_config(self, client):
        r = await client.get("/api/config")
        assert r.status_code == 200


class TestDataApi:
    """資料設定 API"""

    async def test_data_status(self, client):
        """GET /api/data/status 回傳 3 個必填欄位"""
        r = await client.get("/api/data/status")
        assert r.status_code == 200
        data = r.json()
        assert data["status"] == "ok"
        status = data["data"]
        for key in ["he_image", "binned_002", "binned_008"]:
            assert key in status

    async def test_data_status_xenium_optional(self, client):
        """xenium_outs 是 add-on 對照資料，不列入必填；若出現須具備標準欄位"""
        r = await client.get("/api/data/status")
        status = r.json()["data"]
        if "xenium_outs" in status:
            assert set(status["xenium_outs"]) >= {"path", "configured"}

    async def test_data_status_configured(self, client):
        """CRC 資料已設定 → configured = True"""
        r = await client.get("/api/data/status")
        status = r.json()["data"]
        assert status["he_image"]["configured"] is True
        assert status["binned_002"]["configured"] is True

    async def test_browse_home(self, client):
        """GET /api/data/browse?path=~ 回傳目錄列表"""
        r = await client.get("/api/data/browse", params={"path": "~"})
        assert r.status_code == 200
        data = r.json()
        assert data["status"] == "ok"
        assert "current" in data["data"]
        assert "items" in data["data"]
        assert len(data["data"]["items"]) > 0

    async def test_browse_invalid_path(self, client):
        """Browse 不存在路徑回傳 error"""
        r = await client.get("/api/data/browse", params={"path": "/nonexistent/path/abc"})
        assert r.status_code == 200
        data = r.json()
        assert data["status"] == "error"

    async def test_scan_crc(self, client):
        """POST /api/data/scan 掃描 CRC 根目錄"""
        r = await client.post("/api/data/scan", json={
            "data_root": "/Volumes/SSD/plan_a/tissue sample/CRC"
        })
        assert r.status_code == 200
        data = r.json()
        assert data["status"] == "ok"
        result = data["data"]
        assert result["he_image"] is not None


class TestRoiApi:
    """ROI API"""

    async def test_roi_list(self, client):
        """GET /api/roi/list 回傳 ROI 列表"""
        r = await client.get("/api/roi/list")
        assert r.status_code == 200

    async def test_roi_status(self, client):
        """GET /api/roi/status 回傳狀態"""
        r = await client.get("/api/roi/status")
        assert r.status_code == 200


class TestStageStatusApis:
    """各 Stage 的 /status 端點"""

    @pytest.mark.parametrize("endpoint", [
        "/api/segmentation/status",
        "/api/zarr/status",
        "/api/conditions/status",
        "/api/proseg/status",
        "/api/analysis/status",
    ])
    async def test_stage_status(self, client, endpoint):
        """Stage 狀態端點回傳 200"""
        r = await client.get(endpoint)
        assert r.status_code == 200


class TestStaticAssetServing:
    """生產模式（`/` 直接吃 frontend/dist）的靜態資源

    SPA catch-all 若註冊在 `/assets` 掛載之前，會把 `/assets/*.js` 一起吃掉並
    回傳 index.html —— 瀏覽器把 HTML 當 ES module 載入就**整頁空白且無錯誤訊息**。
    掛載順序是有意義的，用測試釘死。
    """

    async def test_js_assets_are_not_swallowed_by_spa_fallback(self, client):
        from pathlib import Path

        dist = Path(__file__).parents[2] / "frontend" / "dist" / "assets"
        js = sorted(dist.glob("*.js")) if dist.exists() else []
        if not js:
            pytest.skip("frontend 尚未 build，無 dist/assets")

        r = await client.get(f"/assets/{js[0].name}")

        assert r.status_code == 200
        assert "javascript" in r.headers["content-type"], (
            f"/assets/*.js 回傳 {r.headers['content-type']}，"
            "表示被 SPA catch-all 攔截（掛載順序錯誤）"
        )

    async def test_deep_link_still_returns_index_html(self, client):
        from pathlib import Path

        if not (Path(__file__).parents[2] / "frontend" / "dist" / "index.html").exists():
            pytest.skip("frontend 尚未 build")

        r = await client.get("/roi")

        assert r.status_code == 200
        assert "text/html" in r.headers["content-type"]
