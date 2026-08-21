#!/bin/bash
# 下載 uv 的 standalone 執行檔，當 Tauri sidecar 用（見 tauri.conf.json
# bundle.externalBin）。這些是平台專屬的執行檔（各 ~20MB），刻意不進版控
# （見 .gitignore）——公開 repo 不該塞這麼大的二進位檔，需要打包前先跑這支。
#
# 用法：bash frontend/src-tauri/binaries/fetch-uv.sh [uv 版本，預設 latest]

set -euo pipefail
cd "$(dirname "$0")"

UV_VERSION="${1:-latest}"
if [ "$UV_VERSION" = "latest" ]; then
    UV_VERSION=$(curl -sL https://api.github.com/repos/astral-sh/uv/releases/latest | python3 -c "import json,sys; print(json.load(sys.stdin)['tag_name'])")
fi
echo "抓取 uv $UV_VERSION ..."

declare -A TARGETS=(
    ["uv-aarch64-apple-darwin.tar.gz"]="uv-aarch64-apple-darwin"
    ["uv-x86_64-apple-darwin.tar.gz"]="uv-x86_64-apple-darwin"
    ["uv-x86_64-pc-windows-msvc.zip"]="uv-x86_64-pc-windows-msvc.exe"
)

TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

for asset in "${!TARGETS[@]}"; do
    dest="${TARGETS[$asset]}"
    url="https://github.com/astral-sh/uv/releases/download/${UV_VERSION}/${asset}"
    echo "  $asset -> $dest"
    curl -sL -o "$TMP/$asset" "$url"
    curl -sL -o "$TMP/$asset.sha256" "$url.sha256"
    (cd "$TMP" && shasum -a 256 -c "$asset.sha256")

    case "$asset" in
        *.tar.gz)
            tar xzf "$TMP/$asset" -C "$TMP"
            src_bin="$TMP/${asset%.tar.gz}/uv"
            ;;
        *.zip)
            unzip -o -q "$TMP/$asset" -d "$TMP/${asset%.zip}"
            src_bin="$TMP/${asset%.zip}/uv.exe"
            ;;
    esac
    cp "$src_bin" "./$dest"
    chmod +x "./$dest" 2>/dev/null || true
done

echo "完成。已放置於 frontend/src-tauri/binaries/："
ls -la ./uv-* 2>/dev/null | grep -v "\.sh$"
