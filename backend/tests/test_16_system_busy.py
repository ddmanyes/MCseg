"""Test 16: /api/system/busy —— 跨 Stage 彙總「有沒有任何工作在跑」（不需真實資料）"""
import pytest

pytestmark = pytest.mark.asyncio


@pytest.fixture(autouse=True)
def _reset_all_status():
    """每個測試前後都把各模組狀態重設回 idle，避免測試互相污染。"""
    from backend.src.api import analysis, cellpose_count, export, roi, segmentation

    def _idle():
        return {"status": "idle", "progress": 0.0, "message": ""}

    yield

    roi._task_status = _idle()
    segmentation._task_status = _idle()
    segmentation._full_status = _idle()
    cellpose_count._status = _idle()
    cellpose_count._full_status = _idle()
    analysis._task_status = _idle()
    export._xenium_status = _idle()
    export._loupe_status = _idle()
    export._result_status = _idle()


class TestGetBusy:
    async def test_all_idle_reports_not_busy(self):
        from backend.src.api.system import get_busy

        result = await get_busy()
        assert result["busy"] is False
        assert result["running_stages"] == []

    async def test_one_stage_running_reports_busy(self):
        from backend.src.api import segmentation
        from backend.src.api.system import get_busy

        segmentation._task_status = {"status": "running", "progress": 0.3, "message": "..."}
        result = await get_busy()
        assert result["busy"] is True
        assert result["running_stages"] == ["segmentation"]

    async def test_reads_live_value_not_stale_import(self):
        """回歸：必須用『匯入模組再存取屬性』，不能用『from module import name』。

        後者會在匯入當下就把 dict 參照複製走——之後模組重新賦值
        `_task_status = {...}`（這是專案內所有 Stage 更新狀態的既有寫法，
        不是原地 mutate），用 `from module import name` 拿到的參照就會
        變成過期的舊值，永遠讀不到真正在跑的狀態。
        """
        from backend.src.api import roi
        from backend.src.api.system import get_busy

        # 先確認初始狀態是 idle
        result = await get_busy()
        assert "roi" not in result["running_stages"]

        # 模組重新賦值（整個換一個新 dict，模擬 roi.py 實際的更新方式）
        roi._task_status = {"status": "running", "progress": 0.1, "message": "開始裁切..."}

        result = await get_busy()
        assert "roi" in result["running_stages"]

    async def test_multiple_stages_running_all_listed(self):
        from backend.src.api import analysis, export
        from backend.src.api.system import get_busy

        analysis._task_status = {"status": "running", "progress": 0.5, "message": "..."}
        export._xenium_status = {"status": "running", "progress": 0.2, "message": "..."}

        result = await get_busy()
        assert result["busy"] is True
        assert set(result["running_stages"]) == {"analysis", "export_xenium"}

    async def test_done_and_error_not_counted_as_busy(self):
        from backend.src.api import cellpose_count
        from backend.src.api.system import get_busy

        cellpose_count._status = {"status": "done", "progress": 1.0, "message": "..."}
        cellpose_count._full_status = {"status": "error", "progress": 0.0, "message": "..."}

        result = await get_busy()
        assert result["busy"] is False
