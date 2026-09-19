[CmdletBinding()]
param(
    [string]$Version = "0.12.5"
)

$ErrorActionPreference = "Stop"
Set-StrictMode -Version Latest

if ([Environment]::OSVersion.Platform -ne [PlatformID]::Win32NT) {
    throw "fetch-uv.ps1 can only run on Windows."
}

$scriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$binariesDir = Split-Path -Parent $scriptDir
$destination = Join-Path $binariesDir "binaries\uv-x86_64-pc-windows-msvc.exe"
$tag = if ($Version.StartsWith("v")) { $Version } else { "v$Version" }
$asset = "uv-x86_64-pc-windows-msvc.zip"
$baseUrl = "https://github.com/astral-sh/uv/releases/download/$tag"
$tempDir = Join-Path ([System.IO.Path]::GetTempPath()) ("mcseg-uv-" + [guid]::NewGuid().ToString("N"))

New-Item -ItemType Directory -Path $tempDir | Out-Null
try {
    $archive = Join-Path $tempDir $asset
    $checksumFile = "$archive.sha256"
    Write-Host "[fetch-uv] Downloading uv $tag ..."
    Invoke-WebRequest -Uri "$baseUrl/$asset" -OutFile $archive
    Invoke-WebRequest -Uri "$baseUrl/$asset.sha256" -OutFile $checksumFile

    $expectedHash = ((Get-Content -Raw $checksumFile).Trim() -split '\s+')[0].ToUpperInvariant()
    $actualHash = (Get-FileHash -Algorithm SHA256 $archive).Hash.ToUpperInvariant()
    if ($actualHash -ne $expectedHash) {
        throw "uv SHA256 mismatch: expected $expectedHash, got $actualHash"
    }

    $extractDir = Join-Path $tempDir "extracted"
    Expand-Archive -LiteralPath $archive -DestinationPath $extractDir
    $uvExecutable = Get-ChildItem -Path $extractDir -Filter "uv.exe" -Recurse -File |
        Select-Object -First 1
    if (-not $uvExecutable) {
        throw "uv.exe was not found in the downloaded archive."
    }

    Copy-Item -LiteralPath $uvExecutable.FullName -Destination $destination -Force
    Write-Host "[fetch-uv] SHA256 verified. Installed: $destination"
}
finally {
    if (Test-Path -LiteralPath $tempDir) {
        Remove-Item -LiteralPath $tempDir -Recurse -Force
    }
}
