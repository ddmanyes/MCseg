#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRC_TAURI_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
FRONTEND_DIR="$(cd "$SRC_TAURI_DIR/.." && pwd)"
TARGET_DIR="$(cd "$SRC_TAURI_DIR" && cargo metadata --format-version=1 --no-deps | node -e 'let s="";process.stdin.on("data",x=>s+=x);process.stdin.on("end",()=>console.log(JSON.parse(s).target_directory))')"
VERSION="$(node -p "require('$SRC_TAURI_DIR/tauri.conf.json').version")"
case "$(uname -m)" in
  arm64) ARCH=aarch64 ;;
  x86_64) ARCH=x64 ;;
  *) echo "Unsupported macOS architecture" >&2; exit 1 ;;
esac
DMG="$TARGET_DIR/release/bundle/dmg/mcseg_${VERSION}_${ARCH}.dmg"
[ -f "$DMG" ] || { echo "Expected new DMG missing: $DMG" >&2; exit 1; }
node "$FRONTEND_DIR/scripts/finalize-release.mjs" "$DMG"
