"""
MSseg 環境自動設定腳本
功能：
  1. 偵測平台（win32 / darwin / linux）
  2. 修復 Windows 上 macOS XSym .venv symlink 問題
  3. 清除 VIRTUAL_ENV 環境變數衝突
  4. 設定 PYTHONIOENCODING=utf-8（避免 Windows cp950 崩潰）
  5. 偵測 GPU（CUDA / MPS / CPU）並輸出建議
  6. 驗證核心套件（cellpose, scanpy, torch）
"""
from __future__ import annotations

import os
import platform
import subprocess
import sys
from pathlib import Path


def _print(msg: str) -> None:
    print(msg, flush=True)


# ── 1. 基礎環境資訊 ──────────────────────────────────────────────────────────

def check_platform() -> str:
    plat = sys.platform
    py_ver = platform.python_version()
    _print(f"[INFO] Platform : {plat} / Python {py_ver}")
    return plat


# ── 2. Windows .venv XSym 修復 ───────────────────────────────────────────────

def fix_venv_symlink(project_root: Path) -> None:
    """macOS 建立的 .venv XSym 在 Windows 上是普通檔案，非目錄，導致 uv run 失敗。"""
    if sys.platform != "win32":
        return

    venv_path = project_root / ".venv"
    if not venv_path.exists():
        return

    if venv_path.is_dir():
        _print("[OK] .venv is a valid directory.")
        return

    # XSym 檔案：不是目錄也不是正常 symlink
    _print("[WARN] .venv appears to be a macOS XSym file — recreating for Windows...")
    backup = project_root / ".venv_mac_symlink"
    venv_path.rename(backup)
    _print(f"[INFO] Backed up to {backup}")

    result = subprocess.run(
        ["uv", "sync"],
        cwd=str(project_root),
        capture_output=True,
        text=True,
        env={**os.environ, "VIRTUAL_ENV": ""},
    )
    if result.returncode == 0:
        _print("[OK] uv sync completed — new .venv created.")
    else:
        _print(f"[ERROR] uv sync failed:\n{result.stderr.strip()}")
        _print("  Manual fix: cd <project_root> && uv sync")


# ── 3. VIRTUAL_ENV 衝突偵測 ──────────────────────────────────────────────────

def check_virtual_env_conflict() -> None:
    venv = os.environ.get("VIRTUAL_ENV", "")
    if not venv:
        return

    venv_path = Path(venv)
    project_venv = Path.cwd() / ".venv"

    if venv_path.resolve() != project_venv.resolve():
        _print(f"[WARN] VIRTUAL_ENV={venv} conflicts with project .venv")
        _print("  Fix: prefix commands with  VIRTUAL_ENV=\"\"  uv run python ...")
        _print("  Or add to your shell profile:  unset VIRTUAL_ENV")
    else:
        _print("[OK] VIRTUAL_ENV points to project .venv.")


# ── 4. PYTHONIOENCODING ──────────────────────────────────────────────────────

def ensure_utf8_encoding() -> None:
    enc = os.environ.get("PYTHONIOENCODING", "")
    if enc.lower() in ("utf-8", "utf8"):
        _print("[OK] PYTHONIOENCODING=utf-8")
        return

    if sys.platform == "win32":
        _print("[WARN] PYTHONIOENCODING not set — may cause UnicodeEncodeError on Windows")
        _print("  Fix: set PYTHONIOENCODING=utf-8  (or add to .env / shell profile)")
    else:
        _print("[INFO] PYTHONIOENCODING not set (non-Windows, usually fine)")


# ── 5. GPU 偵測 ──────────────────────────────────────────────────────────────

def check_gpu() -> str:
    """Returns: 'cuda' | 'mps' | 'cpu'"""
    try:
        import torch  # noqa: PLC0415
    except ImportError:
        _print("[WARN] torch not installed — cannot check GPU.")
        return "unknown"

    _print(f"[INFO] torch version : {torch.__version__}")

    if torch.cuda.is_available():
        name = torch.cuda.get_device_name(0)
        mem_gb = torch.cuda.get_device_properties(0).total_memory / 1e9
        _print(f"[OK] CUDA GPU : {name}  ({mem_gb:.1f} GB)")
        _print(f"[INFO] CUDA version   : {torch.version.cuda}")
        return "cuda"

    if hasattr(torch.backends, "mps") and torch.backends.mps.is_available():
        _print("[OK] Apple MPS (Metal) GPU available.")
        return "mps"

    # CPU fallback — warn on Windows with NVIDIA card
    _print("[INFO] No GPU acceleration found — running on CPU.")
    if sys.platform == "win32":
        _print("[WARN] On Windows with NVIDIA GPU, ensure torch+cuXXX is installed:")
        _print("  uv add torch torchvision --index-url https://download.pytorch.org/whl/cu128")
    return "cpu"


# ── 6. 核心套件驗證 ──────────────────────────────────────────────────────────

_REQUIRED_PACKAGES = [
    ("cellpose",    "cellpose"),
    ("scanpy",      "scanpy"),
    ("anndata",     "anndata"),
    ("skimage",     "scikit-image"),
    ("zarr",        "zarr"),
    ("tifffile",    "tifffile"),
    ("pandas",      "pandas"),
    ("numpy",       "numpy"),
    ("scipy",       "scipy"),
    ("pyarrow",     "pyarrow"),
    ("yaml",        "pyyaml"),
]


def check_packages() -> list[str]:
    missing = []
    for mod, pkg in _REQUIRED_PACKAGES:
        try:
            __import__(mod)
        except ImportError:
            missing.append(pkg)

    if missing:
        _print(f"[WARN] Missing packages: {', '.join(missing)}")
        _print(f"  Fix: uv add {' '.join(missing)}")
    else:
        _print(f"[OK] All {len(_REQUIRED_PACKAGES)} required packages found.")

    return missing


# ── 主程式 ───────────────────────────────────────────────────────────────────

def main() -> None:
    _print("=" * 56)
    _print("  MSseg Environment Setup Check")
    _print("=" * 56)

    # 找到包含 pyproject.toml 的根目錄（相容 scripts/ 和 skills/scripts/ 兩種位置）
    script_path = Path(__file__).resolve()
    project_root = next(
        (p for p in [script_path.parent, script_path.parent.parent, script_path.parent.parent.parent]
         if (p / "pyproject.toml").exists()),
        script_path.parent,
    )
    _print(f"[INFO] Project root : {project_root}")

    check_platform()
    _print("")

    fix_venv_symlink(project_root)
    _print("")

    check_virtual_env_conflict()
    _print("")

    ensure_utf8_encoding()
    _print("")

    gpu = check_gpu()
    _print("")

    missing = check_packages()
    _print("")

    # Summary
    _print("=" * 56)
    if not missing:
        _print(f"[OK] Environment ready.  GPU={gpu.upper()}")
    else:
        _print(f"[!!] Fix missing packages before running MSseg skills.")
    _print("=" * 56)


if __name__ == "__main__":
    main()
