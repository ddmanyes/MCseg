"""
資料自動發現模組

掃描使用者指定的資料根目錄，自動辨識 10X Genomics 標準輸出：
- Visium HD SpaceRanger 輸出（square_002um, square_008um, H&E image）
- Xenium outs（transcripts.parquet, cells.parquet 等）
"""
from __future__ import annotations

import logging
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

logger = logging.getLogger("pipeline.discovery")

# ── 已知檔案 patterns ────────────────────────────────────────

# H&E 影像：BTF 或大型 TIFF（> 100 MB）
_HE_EXTENSIONS = {".btf", ".tif", ".tiff"}
_HE_MIN_SIZE = 100 * 1024 * 1024  # 100 MB

# Visium HD SpaceRanger binned 目錄的必要檔案
_BINNED_REQUIRED = "filtered_feature_bc_matrix.h5"

# Xenium outs 的必要檔案
_XENIUM_REQUIRED = {"transcripts.parquet", "cells.parquet"}

# 掃描深度限制
_MAX_DEPTH = 6


# ── 資料結構 ─────────────────────────────────────────────────

@dataclass
class DiscoveredFile:
    """單一發現的檔案/目錄"""
    path: str
    label: str
    size_bytes: int = 0
    size_human: str = ""
    exists: bool = True


@dataclass
class DiscoveryResult:
    """掃描結果"""
    data_root: str
    he_image: Optional[DiscoveredFile] = None
    binned_002: Optional[DiscoveredFile] = None
    binned_008: Optional[DiscoveredFile] = None
    xenium_outs: Optional[DiscoveredFile] = None
    extra_files: list[DiscoveredFile] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        result = {"data_root": self.data_root, "warnings": self.warnings}
        for key in ("he_image", "binned_002", "binned_008", "xenium_outs"):
            item = getattr(self, key)
            if item:
                result[key] = {
                    "path": item.path,
                    "label": item.label,
                    "size_bytes": item.size_bytes,
                    "size_human": item.size_human,
                }
            else:
                result[key] = None
        result["extra_files"] = [
            {"path": f.path, "label": f.label, "size_human": f.size_human}
            for f in self.extra_files
        ]
        return result

    def to_paths_dict(self) -> dict[str, str]:
        """轉為 pipeline.yaml paths 更新用的字典（只含找到的項目）"""
        paths = {}
        if self.he_image:
            paths["he_image"] = self.he_image.path
        if self.binned_002:
            paths["binned_002"] = self.binned_002.path
        if self.binned_008:
            paths["binned_008"] = self.binned_008.path
        if self.xenium_outs:
            paths["xenium_outs"] = self.xenium_outs.path
        return paths


# ── 工具函式 ─────────────────────────────────────────────────

def _human_size(nbytes: int) -> str:
    """將位元組轉為人類可讀格式"""
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if abs(nbytes) < 1024.0:
            return f"{nbytes:.1f} {unit}"
        nbytes /= 1024.0
    return f"{nbytes:.1f} PB"


def _dir_size(dirpath: Path) -> int:
    """計算目錄的總檔案大小（頂層）"""
    total = 0
    try:
        for f in dirpath.iterdir():
            if f.is_file() and not f.name.startswith("._"):
                total += f.stat().st_size
    except PermissionError:
        pass
    return total


# ── 核心掃描邏輯 ─────────────────────────────────────────────

def scan_data_root(data_root: str | Path) -> DiscoveryResult:
    """
    掃描資料根目錄，自動辨識 Pipeline 所需的各項輸入檔案。

    Parameters
    ----------
    data_root : 資料根目錄路徑

    Returns
    -------
    DiscoveryResult 結構化結果
    """
    root = Path(os.path.expanduser(str(data_root))).resolve()
    result = DiscoveryResult(data_root=str(root))

    if not root.exists():
        result.warnings.append(f"目錄不存在：{root}")
        return result

    if not root.is_dir():
        result.warnings.append(f"路徑不是目錄：{root}")
        return result

    logger.info(f"開始掃描資料目錄：{root}")

    # 收集候選
    he_candidates: list[tuple[Path, int]] = []
    binned_002_candidates: list[Path] = []
    binned_008_candidates: list[Path] = []
    xenium_candidates: list[Path] = []

    for dirpath, dirnames, filenames in os.walk(str(root)):
        # 深度限制
        depth = str(dirpath).count(os.sep) - str(root).count(os.sep)
        if depth > _MAX_DEPTH:
            dirnames.clear()
            continue

        # 跳過隱藏目錄和結果目錄
        dirnames[:] = [
            d for d in dirnames
            if not d.startswith(".") and d not in ("results", "node_modules", ".git", "__pycache__")
        ]

        current = Path(dirpath)
        fnames_set = set(filenames)

        # 1. 偵測 binned 目錄
        if current.name == "square_002um" and _BINNED_REQUIRED in fnames_set:
            binned_002_candidates.append(current)
        elif current.name == "square_008um" and _BINNED_REQUIRED in fnames_set:
            binned_008_candidates.append(current)

        # 2. 偵測 Xenium outs
        if _XENIUM_REQUIRED.issubset(fnames_set):
            xenium_candidates.append(current)

        # 3. 偵測 H&E 影像
        for fname in filenames:
            if fname.startswith("._"):
                continue
            fpath = current / fname
            suffix = fpath.suffix.lower()
            if suffix in _HE_EXTENSIONS:
                try:
                    fsize = fpath.stat().st_size
                    # BTF 直接加入，大型 TIFF 需 > 100MB
                    if suffix == ".btf" or fsize > _HE_MIN_SIZE:
                        # 排除 morphology.ome.tif（Xenium 的形態影像）
                        if "morphology" not in fname.lower():
                            he_candidates.append((fpath, fsize))
                except OSError:
                    pass

    # ── 選擇最佳候選 ──────────────────────────────────────────

    # H&E：優先 BTF，其次最大的 TIFF
    if he_candidates:
        # BTF 優先
        btf = [c for c in he_candidates if c[0].suffix.lower() == ".btf"]
        if btf:
            best = max(btf, key=lambda c: c[1])
        else:
            best = max(he_candidates, key=lambda c: c[1])
        result.he_image = DiscoveredFile(
            path=str(best[0]),
            label=best[0].name,
            size_bytes=best[1],
            size_human=_human_size(best[1]),
        )
        # 其他候選加入 extra_files
        for c in he_candidates:
            if c[0] != best[0]:
                result.extra_files.append(DiscoveredFile(
                    path=str(c[0]),
                    label=f"(備選 H&E) {c[0].name}",
                    size_bytes=c[1],
                    size_human=_human_size(c[1]),
                ))
    else:
        result.warnings.append("未找到 H&E BTF/TIFF 影像")

    # binned_002
    if binned_002_candidates:
        best = binned_002_candidates[0]
        dsize = _dir_size(best)
        result.binned_002 = DiscoveredFile(
            path=str(best),
            label=f"{best.parent.parent.name}/.../{best.name}",
            size_bytes=dsize,
            size_human=_human_size(dsize),
        )
    else:
        result.warnings.append("未找到 square_002um 目錄（Visium HD 2µm binned）")

    # binned_008
    if binned_008_candidates:
        best = binned_008_candidates[0]
        dsize = _dir_size(best)
        result.binned_008 = DiscoveredFile(
            path=str(best),
            label=f"{best.parent.parent.name}/.../{best.name}",
            size_bytes=dsize,
            size_human=_human_size(dsize),
        )
    else:
        result.warnings.append("未找到 square_008um 目錄（Visium HD 8µm binned）")

    # Xenium outs
    if xenium_candidates:
        best = xenium_candidates[0]
        dsize = _dir_size(best)
        result.xenium_outs = DiscoveredFile(
            path=str(best),
            label=best.name if best.name != "outs" else f"{best.parent.name}/{best.name}",
            size_bytes=dsize,
            size_human=_human_size(dsize),
        )
    else:
        result.warnings.append("未找到 Xenium outs 目錄（可選）")

    # 摘要
    found = sum(1 for x in [result.he_image, result.binned_002, result.binned_008, result.xenium_outs] if x)
    logger.info(f"掃描完成：找到 {found}/4 項目")

    return result
