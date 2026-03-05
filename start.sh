#!/bin/bash
# VisiumHD Pipeline 2 — 一鍵啟動（開發模式）
# 使用方式：bash start.sh

ROOT="$(cd "$(dirname "$0")" && pwd)"

# 顏色
RED='\033[0;31m'; GREEN='\033[0;32m'; BLUE='\033[0;34m'; YELLOW='\033[1;33m'; NC='\033[0m'

echo -e "${BLUE}VisiumHD Pipeline 2${NC}"
echo "Root: $ROOT"

# 檢查 uv
if ! command -v uv &>/dev/null; then
    echo -e "${RED}錯誤：找不到 uv，請先安裝：curl -LsSf https://astral.sh/uv/install.sh | sh${NC}"
    exit 1
fi

# ExFAT 防護：把 .venv 實體放在本機 APFS（避免 ._* resource fork 破壞安裝）
LOCAL_VENV="$HOME/.venvs/visiumHD_pipeline_2"
VENV_DIR="$ROOT/.venv"
export UV_CACHE_DIR="$HOME/.cache/uv"

if [ -L "$VENV_DIR" ] && [ -d "$(readlink "$VENV_DIR")" ]; then
    : # symlink 正常
elif [ -d "$VENV_DIR" ] && [ ! -L "$VENV_DIR" ]; then
    echo -e "${YELLOW}清理 ExFAT .venv，改用本機 symlink...${NC}"
    find "$VENV_DIR" -name "._*" -delete 2>/dev/null || true
    rm -rf "$VENV_DIR"
    mkdir -p "$LOCAL_VENV"
    ln -s "$LOCAL_VENV" "$VENV_DIR"
else
    echo -e "${YELLOW}建立本機 .venv symlink...${NC}"
    mkdir -p "$LOCAL_VENV"
    ln -sf "$LOCAL_VENV" "$VENV_DIR"
fi

# 清除舊行程（避免 Address already in use）
echo -e "${YELLOW}清除舊行程（port 8000/3000）...${NC}"
lsof -ti:8000,3000 | xargs kill -9 2>/dev/null || true
sleep 1

# 啟動後端（背景）
echo -e "${GREEN}[1/2] 啟動後端（port 8000）...${NC}"
cd "$ROOT"
uv run uvicorn backend.main:app --reload --port 8000 &
BACKEND_PID=$!

# 等待後端健康確認（最多 10 秒）
echo -n "等待後端就緒"
BACKEND_OK=false
for i in $(seq 1 10); do
    sleep 1
    echo -n "."
    if curl -sf http://localhost:8000/api/health &>/dev/null; then
        BACKEND_OK=true
        break
    fi
    # 如果進程已結束代表啟動失敗
    if ! kill -0 "$BACKEND_PID" 2>/dev/null; then
        break
    fi
done
echo ""

if [ "$BACKEND_OK" = false ]; then
    echo -e "${RED}❌ 後端啟動失敗！請嘗試手動修復：${NC}"
    echo -e "${YELLOW}  rm -rf .venv${NC}"
    echo -e "${YELLOW}  mkdir -p ~/.venvs/visiumHD_pipeline_2 && ln -s ~/.venvs/visiumHD_pipeline_2 .venv${NC}"
    echo -e "${YELLOW}  UV_CACHE_DIR=~/.cache/uv uv sync${NC}"
    kill "$BACKEND_PID" 2>/dev/null || true
    exit 1
fi
echo -e "${GREEN}✅ 後端已就緒${NC}"

# 啟動前端（背景）
echo -e "${GREEN}[2/2] 啟動前端（port 3000）...${NC}"
cd "$ROOT/frontend"

if [ ! -d "node_modules" ]; then
    echo "安裝前端依賴..."
    npm install
fi

npm run dev &
FRONTEND_PID=$!

# 等待前端就緒（最多 8 秒，Vite 通常 2-3 秒）
echo -n "等待前端就緒"
for i in $(seq 1 8); do
    sleep 1
    echo -n "."
    if curl -sf http://localhost:3000 &>/dev/null; then
        echo ""
        echo -e "${GREEN}✅ 前端已就緒${NC}"
        break
    fi
    if ! kill -0 "$FRONTEND_PID" 2>/dev/null; then
        echo ""
        echo -e "${RED}❌ 前端啟動失敗！請確認 node_modules 與 vite 設定${NC}"
        break
    fi
done

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
