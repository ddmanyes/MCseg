"""
系統層級端點：跨 Stage 彙總查詢，目前只有「有沒有任何工作在跑」一項。

背景：桌面殼（Tauri）決定「使用者按結束/Quit 時要不要跳確認對話框」需要
知道**任何一個** Stage 是不是還在跑，但每個 Stage 各自獨立追蹤自己的
`_task_status`/`_full_status`（roi.py、segmentation.py、cellpose_count.py、
analysis.py、export.py 五個模組共 9 個獨立的狀態全域變數），沒有共用的
彙總點。與其讓桌面殼一次打 9 支 API 再自己合併，這裡集中做一次——之後
新增 Stage 只需要在下面補一行，桌面殼完全不用跟著改。
"""
from fastapi import APIRouter

from backend.src.api import analysis, cellpose_count, export, roi, segmentation

router = APIRouter()


@router.get("/busy")
async def get_busy():
    """回傳目前是否有任何 Stage 正在執行中，供桌面殼關閉/結束前確認用。

    ⚠️ 這裡刻意 import 模組本身（`from backend.src.api import roi` 這種），
    不是 `from backend.src.api.roi import _task_status`——每個 Stage 的
    狀態變數在任務開始/結束時是整個重新賦值（`_task_status = {...}`），
    不是原地 mutate。`from module import name` 只會抓到匯入當下那個 dict
    的參照，模組之後重新賦值就變成讀到過期舊值的靜默 bug。必須透過模組
    物件動態存取屬性（`roi._task_status`），每次呼叫才讀得到即時狀態。
    """
    checks = {
        "roi": roi._task_status,
        "segmentation": segmentation._task_status,
        "segmentation_fullslide": segmentation._full_status,
        "count": cellpose_count._status,
        "count_fullslide": cellpose_count._full_status,
        "analysis": analysis._task_status,
        "export_xenium": export._xenium_status,
        "export_loupe": export._loupe_status,
        "export_result": export._result_status,
    }
    running = [name for name, status in checks.items() if status.get("status") == "running"]
    return {"status": "ok", "busy": len(running) > 0, "running_stages": running}
