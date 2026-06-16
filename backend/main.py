"""
MSseg — FastAPI 主應用程式

啟動方式：
    uv run uvicorn backend.main:app --reload --port 8001
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
    spatial,
)

logger = logging.getLogger("pipeline.main")

VALID_STAGES = frozenset({"global", "roi", "segmentation", "count", "analysis", "export"})


@asynccontextmanager
async def lifespan(app: FastAPI):
    """啟動/關閉鉤子"""
    config = load_config()
    setup_logging(config.get("global", {}).get("log_level", "INFO"))
    logger.info("MSseg started")
    yield
    logger.info("MSseg shutting down")


app = FastAPI(
    title="MSseg",
    version="0.8.0",
    description="MCseg v2 Visium HD 空間轉錄體分析流水線 API",
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
app.include_router(spatial.router,       prefix="/api/spatial",     tags=["Stage 3.5: Spatial Explorer"])
app.include_router(export.router,        prefix="/api/export",      tags=["Stage 4: Export"])


@app.get("/api/health")
async def health():
    return {"status": "ok", "version": "0.8.0"}


@app.get("/api/config")
async def get_config():
    try:
        cfg = load_config()
        # Sanitize: strip internal filesystem paths before returning
        paths = cfg.get("paths", {})
        safe_paths = {k: str(v) for k, v in paths.items()} if paths else {}
        return {"status": "ok", "data": {**cfg, "paths": safe_paths}}
    except Exception as e:
        logger.error(f"Failed to load config: {e}")
        return {"status": "error", "message": "設定載入失敗"}


# ── WebSocket Log 串流 ───────────────────────────────────────
@app.websocket("/ws/log/{stage}")
async def websocket_log(websocket: WebSocket, stage: str):
    """
    前端連接此端點以接收指定 stage 的即時 log。
    stage 可為：global | roi | segmentation | count | analysis | export
    """
    if stage not in VALID_STAGES:
        await websocket.close(code=1008)  # Policy violation
        return
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
    from fastapi.responses import FileResponse

    # SPA catch-all：所有非 /api 路徑都回傳 index.html，讓前端 React Router 處理
    @app.get("/{full_path:path}", include_in_schema=False)
    async def spa_fallback(full_path: str):
        index = _frontend_dist / "index.html"
        return FileResponse(str(index))

    app.mount("/assets", StaticFiles(directory=str(_frontend_dist / "assets")), name="assets")
