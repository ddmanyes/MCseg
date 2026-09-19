# Desktop build and release verification

Desktop version **0.2.2** keeps the 0.2.1 distribution-payload cleanup and adds the H&E color-deconvolution fix (see `releases/0.2.2.md`). It does not change the segmentation model or scientific benchmark values. The source package still uses its own 0.8.0 version.

## Build requirements

Use a clean Git checkout on a native filesystem, Node.js/npm, Rust, and the platform build prerequisites. Build Windows on Windows with the MSVC/C++ build tools; build Apple Silicon macOS on an Apple Silicon Mac with Xcode command-line tools. No paid CI runner is configured or invoked here.

```bash
cd frontend
npm ci
node --test scripts/tests/*.test.mjs
```

### macOS

```bash
bash src-tauri/binaries/fetch-uv.sh 0.12.5
# Optional: point CARGO_TARGET_DIR at a native-filesystem build directory.
npm run tauri:build
```

The macOS build creates the App with Tauri, then a verified compressed DMG with `hdiutil`. It seals the complete App with a local ad-hoc signature and verifies it before creating the image. This is not Developer ID signing or Apple notarization. It does not automate Finder or require icon-layout scripting.

### Windows (PowerShell)

```powershell
npm run tauri:build:windows
```

The Windows script obtains and verifies uv 0.12.5 when the sidecar is missing. When supplying an existing sidecar, independently verify its version and upstream checksum; the build manifest records its hash. Do not reuse an arbitrary binary merely because its filename matches.

## What the build records

`stage-tauri-resources.mjs` packages tracked Python runtime files and YAML configuration, the freshly built frontend, LICENSE, THIRD_PARTY_NOTICES.md, pyproject.toml, and uv.lock. It excludes state.json, tests, caches, hidden files, source maps, and untracked runtime files. It rejects local absolute paths and saved ROIs in the distributable pipeline configuration.

Release staging requires clean tracked sources. `MCSEG_ALLOW_DIRTY_BUILD=1` is for local development only; the release finalizer rejects such manifests. `BUILD_MANIFEST.json` includes source commit, desktop version, platform/architecture, uv sidecar hash, and hashes of packaged files. The release directory contains the installer, `<installer>.sha256`, and `<installer>.build.json`.

The frontend output is generated during `tauri:prepare`. The manifest captures its bytes; it does not establish reproducible binary identity across toolchains or record downloaded model weights.

## Notices

The root LICENSE is included in the package and the installer license configuration. THIRD_PARTY_NOTICES.md includes verbatim available dependency license/notice texts and uv's project licenses. To refresh the dependency inventory after changing locks or build platform:

```bash
cd frontend/src-tauri
cargo metadata --locked --format-version 1 --filter-platform aarch64-apple-darwin > /tmp/mcseg-cargo-metadata.json
cd ../..
python3 scripts/release/collect_notices.py --cargo-metadata /tmp/mcseg-cargo-metadata.json
```

From the repository root, use the appropriate temporary path and `x86_64-pc-windows-msvc` filter on Windows. Commit updated notices before release staging. Some installed packages do not provide standalone license text; the generated inventory names them. This inventory is not a complete audit of every binary's transitive licensing.

## Acceptance before publishing

1. Run packaging tests, the frontend production build, and backend health/config/lifecycle regression tests.
2. Build from clean committed sources and retain the build JSON plus checksums.
3. Mount the DMG read-only or unpack the NSIS installer. Confirm no state.json, tests, cache, secret files, or developer input paths; verify LICENSE/notices and every manifest payload hash.
4. Scan the source and extracted package for credentials.
5. Test a fresh installation on the target platform, starting with no application data directory. Check empty Data Setup/ROI state, setup logs, backend health, shutdown, and a small example analysis. Static payload checks do not replace this step.
6. Publish under a desktop-specific tag pinned to the build source commit; never move the historical v0.8.0 tag. Mark untested platform installers as unavailable/pending rather than implying validation.

Older 0.2.0 installations may already have user data in Application Support / AppData. This cleanup does not delete existing user state automatically. New clean installs receive no preselected dataset or ROI; preserve or reset an existing user's state only with their authorization.
