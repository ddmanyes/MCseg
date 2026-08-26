#!/usr/bin/env bash
# tauri build 的產物（dmg/app）本身留在 CARGO_TARGET_DIR（見
# ../.cargo/config.toml：刻意挪去本機 APFS，避開 ExFAT 上 AppleDouble
# ._* 影子檔讓 build script 解析炸掉的問題，不能改回專案內）。
#
# 這支腳本只負責把打包好的 dmg 額外複製一份進專案的 release/，方便
# 直接從 MSseg 資料夾內取用安裝檔（例如要拿去另一台機器測試）。
# target dir 用 `cargo metadata` 動態查，不在這裡重複寫死路徑，避免
# 以後 .cargo/config.toml 改了、這裡忘記跟著改。
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRC_TAURI_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
RELEASE_DIR="$(cd "$SCRIPT_DIR/../../.." && pwd)/release"

TARGET_DIR="$(cd "$SRC_TAURI_DIR" && cargo metadata --format-version=1 --no-deps 2>/dev/null | jq -r '.target_directory')"
if [ -z "$TARGET_DIR" ] || [ "$TARGET_DIR" = "null" ]; then
    echo "[copy-release] ⚠️ 抓不到 cargo target dir，略過複製" >&2
    exit 0
fi

DMG="$(find "$TARGET_DIR/release/bundle/dmg" -maxdepth 1 -name '*.dmg' 2>/dev/null | head -1)"
if [ -z "$DMG" ]; then
    echo "[copy-release] ⚠️ 在 $TARGET_DIR/release/bundle/dmg 找不到 .dmg，略過複製" >&2
    exit 0
fi

mkdir -p "$RELEASE_DIR"
cp -f "$DMG" "$RELEASE_DIR/"
echo "[copy-release] 已複製 $(basename "$DMG") → $RELEASE_DIR/"
