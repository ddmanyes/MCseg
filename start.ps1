# MSseg — 一鍵啟動（開發模式，Windows / PowerShell）
# 使用方式：  powershell -ExecutionPolicy Bypass -File start.ps1
#     或於 PowerShell 內：  .\start.ps1
#
# 對應 macOS 的 start.sh。差異：
#   - 不需要 .venv symlink / ExFAT ._* 防護（NTFS 無此問題）
#   - 若專案位於 ExFAT 外接磁碟，會自動設定 UV_LINK_MODE=copy

$ErrorActionPreference = "Stop"
$Root = Split-Path -Parent $MyInvocation.MyCommand.Path

Write-Host "MSseg — MCseg v2 Visium HD 分析平台" -ForegroundColor Blue
Write-Host "Root: $Root"

# 檢查 uv
if (-not (Get-Command uv -ErrorAction SilentlyContinue)) {
    Write-Host "錯誤：找不到 uv。請先安裝（PowerShell）：" -ForegroundColor Red
    Write-Host '  powershell -c "irm https://astral.sh/uv/install.ps1 | iex"' -ForegroundColor Yellow
    exit 1
}

# ExFAT 外接磁碟防護：非 NTFS 時改用 copy link mode，避免硬連結失敗
try {
    $drive = (Get-Item $Root).PSDrive.Name
    $fs = (Get-Volume -DriveLetter $drive -ErrorAction SilentlyContinue).FileSystem
    if ($fs -and $fs -ne "NTFS") {
        Write-Host "偵測到非 NTFS 磁碟（$fs），啟用 UV_LINK_MODE=copy..." -ForegroundColor Yellow
        $env:UV_LINK_MODE = "copy"
    }
} catch { }

# 確保 Python 環境已同步
Write-Host "同步 Python 環境（uv sync）..." -ForegroundColor Yellow
uv sync
if ($LASTEXITCODE -ne 0) {
    Write-Host "❌ uv sync 失敗。若磁碟為 ExFAT，請手動執行： `$env:UV_LINK_MODE='copy'; uv sync" -ForegroundColor Red
    exit 1
}

# 清除佔用 port 8001 / 3000 的舊行程
Write-Host "清除舊行程（port 8001/3000）..." -ForegroundColor Yellow
foreach ($port in 8001, 3000) {
    try {
        Get-NetTCPConnection -LocalPort $port -State Listen -ErrorAction SilentlyContinue |
            Select-Object -ExpandProperty OwningProcess -Unique |
            ForEach-Object { Stop-Process -Id $_ -Force -ErrorAction SilentlyContinue }
    } catch { }
}
Start-Sleep -Seconds 1

# 啟動後端（背景）
Write-Host "[1/2] 啟動後端（port 8001）..." -ForegroundColor Green
$backend = Start-Process -PassThru -FilePath "uv" `
    -ArgumentList "run", "uvicorn", "backend.main:app", "--reload", "--port", "8001" `
    -WorkingDirectory $Root

# 等待後端健康確認（最多 15 秒；uvicorn --reload 首次啟動較慢）
Write-Host -NoNewline "等待後端就緒"
$backendOk = $false
for ($i = 0; $i -lt 15; $i++) {
    Start-Sleep -Seconds 1
    Write-Host -NoNewline "."
    try {
        $r = Invoke-WebRequest -Uri "http://localhost:8001/api/health" -UseBasicParsing -TimeoutSec 2
        if ($r.StatusCode -eq 200) { $backendOk = $true; break }
    } catch { }
    if ($backend.HasExited) { break }
}
Write-Host ""

if (-not $backendOk) {
    Write-Host "❌ 後端啟動失敗！請嘗試手動修復：" -ForegroundColor Red
    Write-Host '  $env:UV_LINK_MODE="copy"; uv sync' -ForegroundColor Yellow
    Write-Host '  uv run uvicorn backend.main:app --port 8001' -ForegroundColor Yellow
    if (-not $backend.HasExited) { Stop-Process -Id $backend.Id -Force -ErrorAction SilentlyContinue }
    exit 1
}
Write-Host "✅ 後端已就緒" -ForegroundColor Green

# 啟動前端（背景）
Write-Host "[2/2] 啟動前端（port 3000）..." -ForegroundColor Green
$frontendDir = Join-Path $Root "frontend"
if (-not (Test-Path (Join-Path $frontendDir "node_modules"))) {
    Write-Host "安裝前端依賴..."
    Push-Location $frontendDir; npm install; Pop-Location
}
$frontend = Start-Process -PassThru -FilePath "npm" -ArgumentList "run", "dev" `
    -WorkingDirectory $frontendDir

Write-Host -NoNewline "等待前端就緒"
for ($i = 0; $i -lt 10; $i++) {
    Start-Sleep -Seconds 1
    Write-Host -NoNewline "."
    try {
        $r = Invoke-WebRequest -Uri "http://localhost:3000" -UseBasicParsing -TimeoutSec 2
        if ($r.StatusCode -eq 200) { Write-Host ""; Write-Host "✅ 前端已就緒" -ForegroundColor Green; break }
    } catch { }
    if ($frontend.HasExited) { Write-Host ""; Write-Host "❌ 前端啟動失敗！請確認 node_modules 與 vite 設定" -ForegroundColor Red; break }
}

Write-Host ""
Write-Host "========================================" -ForegroundColor Green
Write-Host "  後端 API:  http://localhost:8001" -ForegroundColor Green
Write-Host "  前端 UI:   http://localhost:3000" -ForegroundColor Green
Write-Host "  API Docs:  http://localhost:8001/docs" -ForegroundColor Green
Write-Host "========================================" -ForegroundColor Green
Write-Host ""
Write-Host "按 Ctrl+C 停止所有服務"

# 等待並在結束時清理子行程
try {
    while (-not $backend.HasExited -and -not $frontend.HasExited) {
        Start-Sleep -Seconds 1
    }
} finally {
    foreach ($p in $backend, $frontend) {
        if ($p -and -not $p.HasExited) { Stop-Process -Id $p.Id -Force -ErrorAction SilentlyContinue }
    }
    Write-Host "已停止"
}
