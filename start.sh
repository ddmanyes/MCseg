#!/bin/bash
# VisiumHD Pipeline 2 — 一鍵啟動（開發模式）
# 使用方式：bash start.sh

set -e
ROOT="$(cd "$(dirname "$0")" && pwd)"

# 顏色
GREEN='\033[0;32m'; BLUE='\033[0;34m'; NC='\033[0m'

echo -e "${BLUE}VisiumHD Pipeline 2${NC}"
echo "Root: $ROOT"

# 檢查 uv
if ! command -v uv &>/dev/null; then
    echo "錯誤：找不到 uv，請先安裝：curl -LsSf https://astral.sh/uv/install.sh | sh"
    exit 1
fi

# 啟動後端（背景）
echo -e "${GREEN}[1/2] 啟動後端（port 8000）...${NC}"
cd "$ROOT"
uv run uvicorn backend.main:app --reload --port 8000 &
BACKEND_PID=$!

sleep 2

# 啟動前端（背景）
echo -e "${GREEN}[2/2] 啟動前端（port 3000）...${NC}"
cd "$ROOT/frontend"

if [ ! -d "node_modules" ]; then
    echo "安裝前端依賴..."
    npm install
fi

npm run dev &
FRONTEND_PID=$!

echo ""
echo -e "${GREEN}========================================${NC}"
echo -e "${GREEN}  後端 API:  http://localhost:8000${NC}"
echo -e "${GREEN}  前端 UI:   http://localhost:3000${NC}"
echo -e "${GREEN}  API Docs:  http://localhost:8000/docs${NC}"
echo -e "${GREEN}========================================${NC}"
echo ""
echo "按 Ctrl+C 停止所有服務"

# 等待並清理
trap "kill $BACKEND_PID $FRONTEND_PID 2>/dev/null; echo '已停止'" INT TERM
wait
