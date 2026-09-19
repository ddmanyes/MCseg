[CmdletBinding()]
param(
    [string]$UvVersion = "0.12.5"
)

$ErrorActionPreference = "Stop"
Set-StrictMode -Version Latest

# powershell -File does not turn a thrown error into a non-zero exit code, so
# a failed sidecar fetch or bundle step used to surface as "npm run" success
# with no installer produced. Fail loudly instead.
trap {
    Write-Host "[build-windows] FAILED: $_" -ForegroundColor Red
    exit 1
}

if ([Environment]::OSVersion.Platform -ne [PlatformID]::Win32NT) {
    throw "The Windows installer must be built on Windows."
}

$scriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$srcTauriDir = Split-Path -Parent $scriptDir
$frontendDir = Split-Path -Parent $srcTauriDir
$repositoryDir = Split-Path -Parent $frontendDir
$releaseDir = Join-Path $repositoryDir "release"
$uvSidecar = Join-Path $srcTauriDir "binaries\uv-x86_64-pc-windows-msvc.exe"
$tauriCli = Join-Path $frontendDir "node_modules\.bin\tauri.cmd"

# The project Cargo config points macOS/ExFAT builds at an APFS path. Override
# it here so Windows does not resolve /Users/... on the current drive.
$env:CARGO_TARGET_DIR = Join-Path $srcTauriDir "target"

if (-not (Test-Path -LiteralPath $uvSidecar)) {
    Write-Host "[build-windows] Windows uv sidecar is missing; downloading it."
    & (Join-Path $scriptDir "fetch-uv.ps1") -Version $UvVersion
}

if (-not (Test-Path -LiteralPath $uvSidecar)) {
    throw "uv sidecar is still missing after the fetch step: $uvSidecar"
}

if (-not (Test-Path -LiteralPath $tauriCli)) {
    throw "Tauri CLI was not found. Run npm install in frontend first."
}

Write-Host "[build-windows] Sidecar: $(& $uvSidecar --version)"
Write-Host "[build-windows] Cargo target: $env:CARGO_TARGET_DIR"

Push-Location $frontendDir
try {
    & $tauriCli build --bundles nsis
    if ($LASTEXITCODE -ne 0) {
        throw "Tauri Windows NSIS build failed (exit code $LASTEXITCODE)."
    }
}
finally {
    Pop-Location
}

$bundleDir = Join-Path $env:CARGO_TARGET_DIR "release\bundle\nsis"
$installer = Get-ChildItem -LiteralPath $bundleDir -Filter "*-setup.exe" -File |
    Sort-Object LastWriteTime -Descending |
    Select-Object -First 1
if (-not $installer) {
    throw "Build completed, but no NSIS installer was found in $bundleDir."
}

& node (Join-Path $frontendDir "scripts\finalize-release.mjs") $installer.FullName
if ($LASTEXITCODE -ne 0) { throw "Release manifest/checksum generation failed." }
Write-Host "[build-windows] Complete: $releaseDir"
