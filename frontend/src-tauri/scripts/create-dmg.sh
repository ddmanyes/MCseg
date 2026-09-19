#!/usr/bin/env bash
# Package the compiled app without Finder automation (works in headless builds).
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRC_TAURI_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
TARGET_DIR="$(cd "$SRC_TAURI_DIR" && cargo metadata --format-version=1 --no-deps | node -e 'let s="";process.stdin.on("data",x=>s+=x);process.stdin.on("end",()=>console.log(JSON.parse(s).target_directory))')"
VERSION="$(node -p "require('$SRC_TAURI_DIR/tauri.conf.json').version")"
case "$(uname -m)" in
  arm64) ARCH=aarch64 ;;
  x86_64) ARCH=x64 ;;
  *) echo 'Unsupported architecture' >&2; exit 1 ;;
esac
APP="$TARGET_DIR/release/bundle/macos/mcseg.app"
[ -d "$APP" ] || { echo 'Build the macOS app before creating a DMG' >&2; exit 1; }
LAYOUT=$(mktemp -d)
trap 'rm -rf "$LAYOUT"' EXIT
ditto "$APP" "$LAYOUT/mcseg.app"
ln -s /Applications "$LAYOUT/Applications"
cp "$SRC_TAURI_DIR/../../LICENSE" "$LAYOUT/LICENSE"
mkdir -p "$TARGET_DIR/release/bundle/dmg"
DMG="$TARGET_DIR/release/bundle/dmg/mcseg_${VERSION}_${ARCH}.dmg"
hdiutil create -ov -format UDZO -fs HFS+ -volname MCseg -srcfolder "$LAYOUT" "$DMG"
hdiutil verify "$DMG"
bash "$SCRIPT_DIR/copy-release.sh"
