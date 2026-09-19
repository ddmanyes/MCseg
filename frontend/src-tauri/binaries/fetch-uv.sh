#!/usr/bin/env bash
# Download and verify the pinned native macOS sidecar. Windows uses fetch-uv.ps1.
set -euo pipefail
cd "$(dirname "$0")"
UV_VERSION="${1:-0.12.5}"
case "$(uname -m)" in
  arm64) TARGET=aarch64-apple-darwin ;;
  x86_64) TARGET=x86_64-apple-darwin ;;
  *) echo 'Unsupported macOS architecture' >&2; exit 1 ;;
esac
ARCHIVE="uv-${TARGET}.tar.gz"
TEMP_DIR=$(mktemp -d)
trap 'rm -rf "$TEMP_DIR"' EXIT
URL="https://github.com/astral-sh/uv/releases/download/${UV_VERSION}/${ARCHIVE}"
curl -fLsS "$URL" -o "$TEMP_DIR/$ARCHIVE"
curl -fLsS "$URL.sha256" -o "$TEMP_DIR/$ARCHIVE.sha256"
(cd "$TEMP_DIR" && shasum -a 256 -c "$ARCHIVE.sha256")
tar xzf "$TEMP_DIR/$ARCHIVE" -C "$TEMP_DIR"
cp "$TEMP_DIR/uv-${TARGET}/uv" "uv-${TARGET}"
chmod +x "uv-${TARGET}"
"./uv-${TARGET}" --version
