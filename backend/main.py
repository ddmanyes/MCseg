"""
VisiumHD Pipeline 3 — FastAPI 主應用程式

啟動方式：
    uv run uvicorn backend.main:app --reload --port 8000
"""
import asyncio
import json
import logging
from contextlib import asynccontextmanager
from pathlib import Path

from fastapi import FastAPI, WebSocket, WebSocketDisconnect
from fastapi.middleware.cors import CORSMiddleware
from fastapi.staticfiles import StaticFiles

from backend.src.utils.config import load_config
from backend.src.utils.logging import (
    register_ws_queue,
    setup_logging,
    unregister_ws_queue,
)
from backend.src.api import (
    analysis,
    cellpose_count,
    data,
    export,
    roi,
    segmentation,
)

logger = logging.getLogger("pipeline.main")


@asynccontextmanager
async def lifespan(app: FastAPI):
    """啟動/關閉鉤子"""
    config = load_config()
    setup_logging(config.get("global", {}).get("log_level", "INFO"))
    logger.info("VisiumHD Pipeline 3 started")
    yield
    logger.info("VisiumHD Pipeline 3 shutting down")


app = FastAPI(
    title="VisiumHD Pipeline 3",
    version="3.0.0",
    description="空間轉錄體分析流水線 API（直接 Cellpose 分析，無 Proseg）",
    lifespan=lifespan,
)

# CORS（允許 React dev server 跨域）
app.add_middleware(
    CORSMiddleware,
    allow_origins=["http://localhost:3000", "http://127.0.0.1:3000"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# ── REST API 路由 ────────────────────────────────────────────
app.include_router(data.router,          prefix="/api/data",        tags=["Data Setup"])
app.include_router(roi.router,           prefix="/api/roi",         tags=["Stage 0: ROI"])
app.include_router(segmentation.router,  prefix="/api/segmentation",tags=["Stage 1: Segmentation"])
app.include_router(cellpose_count.router,prefix="/api/count",       tags=["Stage 2: Count"])
app.include_router(analysis.router,      prefix="/api/analysis",    tags=["Stage 3: Analysis"])
app.include_router(export.router,        prefix="/api/export",      tags=["Stage 4: Export"])


@app.get("/api/health")
async def health():
    return {"status": "ok", "version": "3.0.0"}


@app.get("/api/config")
async def get_config():
    try:
        return {"status": "ok", "data": load_config()}
    except Exception as e:
        return {"status": "error", "message": str(e)}


# ── WebSocket Log 串流 ───────────────────────────────────────
@app.websocket("/ws/log/{stage}")
async def websocket_log(websocket: WebSocket, stage: str):
    """
    前端連接此端點以接收指定 stage 的即時 log。
    stage 可為：global | roi | segmentation | count | analysis | export
    """
    await websocket.accept()
    queue: asyncio.Queue = asyncio.Queue(maxsize=500)
    register_ws_queue(stage, queue)

    try:
        while True:
            try:
                msg = await asyncio.wait_for(queue.get(), timeout=30.0)
                await websocket.send_text(json.dumps(msg))
            except asyncio.TimeoutError:
                # 保活心跳
                await websocket.send_text(json.dumps({"type": "ping"}))
    except WebSocketDisconnect:
        pass
    finally:
        unregister_ws_queue(stage, queue)


# ── 生產模式：服務前端靜態資源 ───────────────────────────────
_frontend_dist = Path(__file__).parent.parent / "frontend" / "dist"
if _frontend_dist.exists():
    app.mount("/", StaticFiles(directory=str(_frontend_dist), html=True), name="frontend")
