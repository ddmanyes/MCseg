"""
Stage 5 — Loupe Browser 匯出模組

將 Proseg 分割結果（h5ad + GeoJSON 多邊形）轉換為：
1. .cloupe 檔案（透過 loupepy）
2. .geojson 檔案（匹配條碼後的細胞邊界，供 Loupe Browser 匯入）

設計要點
--------
- 條碼：優先使用 10X 白名單（backend/src/export/10x_whitelist.txt），
  不存在時自動生成符合格式的假條碼。
- 分類欄位限制：唯一值 > 32,000 自動轉為整數 cat codes（Loupe 限制）。
- GeoJSON：使用 ijson 串流讀取，禁止整體 json.load 大型檔案。

移植自：Proseg-Zarr-Integration/scripts/export_to_loupe_merged.py
"""

from __future__ import annotations

import gzip
import json
import logging
import platform
import random
import shutil
from pathlib import Path
from typing import Optional

import anndata
import scanpy as sc

logger = logging.getLogger("pipeline.export.loupe")

# 白名單檔案位置（與本模組同目錄）
_WHITELIST_PATH = Path(__file__).parent / "10x_whitelist.txt"


class LoupeExporter:
    """
    將 Proseg 分割結果匯出為 Loupe Browser 相容格式。

    Parameters
    ----------
    poly_json_path:
        combined_proseg_results_qc.json 路徑（GeoJSON，可以是 .gz 壓縮）。
        若未提供則只產生 .cloupe，不輸出 GeoJSON。
    whitelist_path:
        10X 白名單條碼檔案路徑。若為 None，使用模組預設路徑，
        不存在時自動生成假條碼。
    """

    def __init__(
        self,
        poly_json_path: Optional[str | Path] = None,
        whitelist_path: Optional[str | Path] = None,
    ) -> None:
        self.poly_json_path = Path(poly_json_path) if poly_json_path else None
        self.whitelist_path = (
            Path(whitelist_path) if whitelist_path else _WHITELIST_PATH
        )

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def export(self, h5ad_path: str | Path, output_dir: str | Path) -> Path:
        """
        主要匯出入口。

        Parameters
        ----------
        h5ad_path:
            已過 QC 的 AnnData h5ad 檔案路徑。
        output_dir:
            輸出目錄（會自動建立）。

        Returns
        -------
        Path
            .cloupe 檔案路徑。
        """
        h5ad_path = Path(h5ad_path)
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        output_cloupe = output_dir / (h5ad_path.stem + ".cloupe")
        output_geojson = output_dir / (h5ad_path.stem + ".geojson")

        logger.info("=== Stage 5 Loupe Browser 匯出開始 ===")
        logger.info(f"h5ad: {h5ad_path}")
        logger.info(f"輸出目錄: {output_dir}")

        # 確認 loupepy / ijson 可用
        self._ensure_dependencies()

        # 1. 載入 AnnData
        adata = self._load_adata(h5ad_path)

        # 2. 指派 10X 條碼
        barcodes, original_id_to_barcode = self._assign_barcodes(adata)

        # 3. 前處理 AnnData（分類欄位、obs_names 替換）
        adata = self._prepare_adata(adata, barcodes)

        # 4. 建立 .cloupe
        self._create_cloupe(adata, output_cloupe)

        # 5. 產生 GeoJSON（串流讀取，不撐爆記憶體）
        if self.poly_json_path and self.poly_json_path.exists():
            self._create_geojson(
                original_id_to_barcode=original_id_to_barcode,
                output_geojson=output_geojson,
            )
        else:
            logger.warning(
                f"多邊形 JSON 不存在或未提供：{self.poly_json_path}，跳過 GeoJSON 輸出。"
            )

        logger.info("=== Loupe Browser 匯出完成 ===")
        logger.info(f"1. 開啟 Loupe Browser，載入：{output_cloupe}")
        if output_geojson.exists():
            logger.info(
                f"2. Tools → Import Custom Region → {output_geojson}"
            )
        return output_cloupe

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _ensure_dependencies() -> None:
        """確認 loupepy 和 ijson 已安裝（否則提示使用 uv add）。"""
        missing = []
        try:
            import loupepy  # noqa: F401
        except ImportError:
            missing.append("loupepy")
        try:
            import ijson  # noqa: F401
        except ImportError:
            missing.append("ijson")

        if missing:
            pkg_list = " ".join(missing)
            raise ImportError(
                f"缺少必要套件：{pkg_list}。請執行：uv add {pkg_list}"
            )

    @staticmethod
    def _load_adata(h5ad_path: Path) -> anndata.AnnData:
        """載入 AnnData。"""
        logger.info(f"載入 AnnData：{h5ad_path}")
        adata = sc.read_h5ad(h5ad_path)
        logger.info(f"AnnData 載入完成：{len(adata)} 個細胞，{adata.n_vars} 個基因。")
        return adata

    def _assign_barcodes(
        self, adata: anndata.AnnData
    ) -> tuple[list[str], dict[str, str]]:
        """
        為每個細胞指派 10X 格式條碼。

        優先順序：
        1. 從白名單檔案讀取
        2. 白名單不足時隨機取樣補齊
        3. 白名單不存在時自動生成假條碼

        Returns
        -------
        barcodes:
            長度為 len(adata) 的條碼列表（含 -1 後綴）。
        original_id_to_barcode:
            原始 obs_names → 條碼 的對應字典。
        """
        num_cells = len(adata)
        logger.info(f"為 {num_cells} 個細胞生成 10X 條碼...")

        if self.whitelist_path.exists():
            with open(self.whitelist_path, "r") as f:
                valid_barcodes = [line.strip() for line in f if line.strip()]

            if num_cells <= len(valid_barcodes):
                barcodes = [f"{valid_barcodes[i]}-1" for i in range(num_cells)]
                logger.info(f"從白名單取用 {num_cells} 個條碼。")
            else:
                # 白名單不足，隨機取樣（含重複）
                logger.warning(
                    f"細胞數（{num_cells}）超過白名單長度（{len(valid_barcodes)}），隨機取樣補齊。"
                )
                random.seed(42)
                barcodes = [
                    f"{random.choice(valid_barcodes)}-1" for _ in range(num_cells)
                ]
        else:
            logger.warning(
                f"白名單不存在：{self.whitelist_path}，自動生成假條碼。"
            )
            barcodes = _generate_fake_barcodes(num_cells)

        original_id_to_barcode = {
            orig: bc for orig, bc in zip(adata.obs_names, barcodes)
        }
        # Merge 模式：obs_names 為 "{roi}__cell_{N}"，GeoJSON full_id 為 "{roi}__{N}"
        # 補上別名映射確保 GeoJSON 能正確對應條碼
        for orig, bc in list(original_id_to_barcode.items()):
            if "__cell_" in orig:
                parts = orig.split("__cell_", 1)
                if len(parts) == 2:
                    alias = f"{parts[0]}__{parts[1]}"
                    original_id_to_barcode.setdefault(alias, bc)
        return barcodes, original_id_to_barcode

    @staticmethod
    def _prepare_adata(
        adata: anndata.AnnData, barcodes: list[str]
    ) -> anndata.AnnData:
        """
        前處理 AnnData 使其符合 Loupe 要求：
        1. 確保 Proseg_Cell_ID 欄位存在。
        2. 分類欄位唯一值 > 32,000 轉為整數 cat codes（Loupe 限制）。
        3. 以條碼替換 obs_names。
        """
        num_cells = len(adata)

        if "Proseg_Cell_ID" not in adata.obs:
            adata.obs["Proseg_Cell_ID"] = range(num_cells)

        for col in adata.obs.columns:
            dtype_str = str(adata.obs[col].dtype)
            if adata.obs[col].dtype == "object" or dtype_str == "category":
                n_unique = adata.obs[col].nunique()
                if n_unique > 32_000:
                    logger.warning(
                        f"欄位 '{col}' 有 {n_unique} 個唯一值，超過 Loupe 32k 限制，"
                        "轉換為整數 cat codes。"
                    )
                    adata.obs[col] = adata.obs[col].astype("category").cat.codes
                else:
                    adata.obs[col] = adata.obs[col].astype("category")

        adata.obs_names = barcodes
        return adata

    @staticmethod
    def _create_cloupe(adata: anndata.AnnData, output_cloupe: Path) -> None:
        """呼叫 loupepy 建立 .cloupe 檔案。"""
        import loupepy

        logger.info("初始化 Loupe Converter 環境...")
        loupepy.setup()

        loupe_path: Optional[str] = None

        # Windows 環境：找 loupe_converter.exe
        if platform.system() == "Windows":
            import os

            default_path = os.path.join(
                os.environ.get("APPDATA", ""), "loupe_converter"
            )
            roaming_path = os.path.join(
                os.environ.get("APPDATA", ""), "loupe_converter.exe"
            )
            if os.path.exists(default_path) and not os.path.exists(roaming_path):
                try:
                    shutil.copy(default_path, roaming_path)
                except Exception:
                    pass
            if os.path.exists(roaming_path):
                loupe_path = roaming_path

        logger.info(f"建立 .cloupe 檔案（大型資料集耗時較長）：{output_cloupe}")
        try:
            if loupe_path:
                loupepy.create_loupe_from_anndata(
                    adata,
                    output_cloupe=str(output_cloupe),
                    loupe_converter_path=loupe_path,
                    force=True,
                )
            else:
                loupepy.create_loupe_from_anndata(
                    adata,
                    output_cloupe=str(output_cloupe),
                    force=True,
                )
            logger.info(f".cloupe 建立成功：{output_cloupe}")
        except Exception as exc:
            logger.error(f".cloupe 建立失敗：{exc}")
            raise

    def _create_geojson(
        self,
        original_id_to_barcode: dict[str, str],
        output_geojson: Path,
    ) -> None:
        """
        串流讀取大型 GeoJSON（ijson），將 full_id 映射為條碼，
        輸出 Loupe 相容的 GeoJSON 檔案。

        使用 ijson 串流讀取，避免整體 json.load 大型檔案導致 OOM。
        """
        import ijson

        logger.info(
            f"串流讀取多邊形 JSON（ijson）：{self.poly_json_path}"
        )

        out_features: list[dict] = []

        try:
            # 同時支援 .gz 壓縮和純文字 JSON
            if self.poly_json_path.name.endswith(".gz"):
                f = gzip.open(self.poly_json_path, "rt", encoding="utf-8")
            else:
                f = open(self.poly_json_path, "rt", encoding="utf-8")

            with f:
                features = ijson.items(f, "features.item")
                for feature in features:
                    full_id = feature.get("properties", {}).get("full_id")
                    if full_id not in original_id_to_barcode:
                        continue

                    barcode = original_id_to_barcode[full_id]
                    out_features.append(
                        {
                            "type": "Feature",
                            "properties": {
                                "barcode": barcode,
                                "cell_id": full_id,
                            },
                            "geometry": feature["geometry"],
                        }
                    )

        except (gzip.BadGzipFile, OSError) as exc:
            logger.error(f"讀取 JSON 檔案失敗：{exc}")
            raise
        except Exception as exc:
            logger.error(f"解析 JSON 資料錯誤：{exc}")
            raise

        out_geo = {"type": "FeatureCollection", "features": out_features}
        logger.info(f"寫出 GeoJSON：{output_geojson}（{len(out_features)} 個特徵）")
        with open(output_geojson, "w", encoding="utf-8") as f:
            json.dump(out_geo, f)

        logger.info(f"GeoJSON 寫出完成：{output_geojson}")


# ------------------------------------------------------------------
# 條碼生成工具函式
# ------------------------------------------------------------------

def _generate_fake_barcodes(n: int) -> list[str]:
    """
    生成符合 10X 格式的假條碼（16bp ATCG + -1 後綴）。

    使用固定 seed=42 確保可重現性。
    生成的條碼去重後取前 n 個，確保每個細胞擁有唯一條碼。
    若 n 極大（> 4^16 理論上限）則允許重複。
    """
    bases = "ACGT"
    random.seed(42)
    barcodes: set[str] = set()

    while len(barcodes) < n:
        bc = "".join(random.choices(bases, k=16))
        barcodes.add(bc)

    return [f"{bc}-1" for bc in sorted(barcodes)[:n]]
