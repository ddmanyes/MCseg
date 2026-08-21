"""CLI 全片流程（`backend/src/cli/segment.py`）

`IMPLEMENTATION_PLAN.md` 風險表的 P0-3b：CLI 全片流程原本**完全沒有測試**，
P0-4/P0-6 把 CLI 內部函式搬進 `fullslide/pipeline.py` 時沒有回歸網。

本檔補兩件事：

1. **端到端煙霧測試**：合成 BTF ＋ 少量 bins 走完 `--btf` 與 `--he-crop` 兩條路徑，
   斷言 `mcseg_mask.npy` / `bin_attribution.parquet` / `cells.h5ad` 皆產出。
2. **座標變換回歸**：CLI 必須與 API 走同一條「對位 JSON homography 優先」路徑，
   且以**來源影像**尺寸推導變換（不是裁切後的遮罩尺寸）。
"""
import json
from pathlib import Path

import numpy as np
import pytest


def _write_align(path, *, scale_transform=0.2, scale_images=0.1, serial="H1-X", area="D1"):
    """寫出 Loupe/CytAssist 對位 JSON（schema 與真實檔一致）。"""
    path.write_text(json.dumps({
        "serialNumber": serial,
        "area": area,
        "transform": [scale_transform, 0, 0, 0, scale_transform, 0, 0, 0, 1],
        "cytAssistInfo": {
            "transformImages": [scale_images, 0, 0, 0, scale_images, 0, 0, 0, 1],
        },
    }), encoding="utf-8")
    return path


def _write_spatial(binned_dir, hires_size, scalef, mpp):
    from PIL import Image

    sp = binned_dir / "spatial"
    sp.mkdir(parents=True, exist_ok=True)
    (sp / "scalefactors_json.json").write_text(json.dumps({
        "microns_per_pixel": mpp,
        "tissue_hires_scalef": scalef,
    }))
    Image.new("RGB", hires_size).save(sp / "tissue_hires_image.png")
    return sp


def _tissue_image(size=512, seed=0):
    """整片都算「組織」（灰階 < 220）的合成影像。"""
    rng = np.random.default_rng(seed)
    return rng.integers(30, 180, (size, size, 3), dtype=np.uint8)


@pytest.fixture
def cli_sample(tmp_path):
    """一份最小但完整的 CLI 輸入：BTF ＋ tissue_positions ＋ h5（替身）。

    影像 512×512；bins 均勻散佈在影像內，確保有 bin 落在細胞上。
    """
    import tifffile

    btf = tmp_path / "slide.btf"
    tifffile.imwrite(str(btf), _tissue_image(512), bigtiff=True,
                     tile=(256, 256), photometric="rgb")

    binned = tmp_path / "b002"
    _write_spatial(binned, hires_size=(128, 128), scalef=0.25, mpp=0.2737)

    # bins 均勻鋪滿整張影像，且座標 ≡ 16 (mod 32) —— 正落在 FakeCellposeModel
    # 給 label 的 32px 網格中央，因此任何裁切窗格內都會有 bin 命中細胞。
    rows_cols = [(y, x) for y in range(16, 512, 96) for x in range(16, 512, 96)]
    import pandas as pd
    pd.DataFrame({
        "barcode": [f"BC{i}" for i in range(len(rows_cols))],
        "in_tissue": [1] * len(rows_cols),
        "pxl_row_in_fullres": [rc[0] for rc in rows_cols],
        "pxl_col_in_fullres": [rc[1] for rc in rows_cols],
    }).to_parquet(str(binned / "spatial" / "tissue_positions.parquet"), index=False)

    # aggregate_cells 內部走 sc.read_10x_h5；此處寫 h5ad 並於測試中改導向
    import anndata as ad
    import scipy.sparse as sp_sparse
    h5 = binned / "filtered_feature_bc_matrix.h5"
    adata = ad.AnnData(X=sp_sparse.csr_matrix(
        np.arange(len(rows_cols) * 3, dtype=np.float32).reshape(len(rows_cols), 3)
    ))
    adata.obs_names = [f"BC{i}" for i in range(len(rows_cols))]
    adata.var_names = ["GeneA", "GeneB", "GeneC"]
    adata.write_h5ad(str(h5))

    return {"btf": btf, "binned": binned, "tp": binned / "spatial" / "tissue_positions.parquet",
            "h5": h5, "out": tmp_path / "out"}


@pytest.fixture
def stub_h5_reader(monkeypatch):
    """把 `sc.read_10x_h5` 導向 `read_h5ad`（fixture 寫的是 h5ad 替身）。"""
    import scanpy as sc

    from backend.src.fullslide import pipeline

    monkeypatch.setattr(pipeline.sc if hasattr(pipeline, "sc") else sc,
                        "read_10x_h5", lambda p: sc.read_h5ad(str(p)))
    monkeypatch.setattr(sc, "read_10x_h5", lambda p: sc.read_h5ad(str(p)))


class TestCliSmoke:
    """P0-3b：CLI 端到端煙霧測試（合成資料，不需 GPU/真實模型）"""

    def test_btf_path_produces_mask_attribution_and_h5ad(
        self, cli_sample, fake_cellpose, stub_h5_reader
    ):
        """`--btf` 走完裁切 → 分割 → attribution → 聚合，四項產出齊全。"""
        from backend.src.cli.segment import main

        rc = main([
            "--btf", str(cli_sample["btf"]),
            "--tp", str(cli_sample["tp"]),
            "--h5", str(cli_sample["h5"]),
            "--out", str(cli_sample["out"]),
            "--no-gpu", "--skip-celltypist",
            "--tile-size", "256", "--overlap", "64",
        ])

        out = cli_sample["out"]
        assert rc == 0
        assert (out / "he_crop.tif").exists()
        assert (out / "mcseg_mask.npy").exists()
        assert (out / "bin_attribution.parquet").exists()
        assert (out / "cells.h5ad").exists()

        mask = np.load(str(out / "mcseg_mask.npy"))
        assert mask.shape == (512, 512)
        assert mask.max() > 0, "分割未產生任何細胞"

        import anndata as ad
        cells = ad.read_h5ad(str(out / "cells.h5ad"))
        assert cells.n_obs > 0
        assert cells.n_vars == 3
        # add_centroids 應補上兩套座標系（P0-7）
        for col in ("centroid_x_px", "centroid_x_fullres"):
            assert col in cells.obs.columns

    def test_he_crop_path_skips_crop_step(self, cli_sample, fake_cellpose, stub_h5_reader):
        """`--he-crop` 直接吃已裁切影像，不再產生 he_crop.tif。"""
        import tifffile

        from backend.src.cli.segment import main

        he = cli_sample["btf"].parent / "he_crop_input.tif"
        tifffile.imwrite(str(he), _tissue_image(512))

        rc = main([
            "--he-crop", str(he),
            "--tp", str(cli_sample["tp"]),
            "--h5", str(cli_sample["h5"]),
            "--out", str(cli_sample["out"]),
            "--no-gpu", "--skip-celltypist",
            "--tile-size", "256", "--overlap", "64",
        ])

        assert rc == 0
        assert (cli_sample["out"] / "mcseg_mask.npy").exists()
        assert (cli_sample["out"] / "cells.h5ad").exists()
        assert not (cli_sample["out"] / "he_crop.tif").exists()

    def test_missing_tp_skips_counting_without_failing(self, cli_sample, fake_cellpose):
        """只給影像時仍應分割成功，計數步驟安靜跳過（不得整體失敗）。"""
        from backend.src.cli.segment import main

        rc = main([
            "--btf", str(cli_sample["btf"]),
            "--out", str(cli_sample["out"]),
            "--no-gpu", "--skip-celltypist",
            "--tile-size", "256", "--overlap", "64",
        ])

        assert rc == 0
        assert (cli_sample["out"] / "mcseg_mask.npy").exists()
        assert not (cli_sample["out"] / "cells.h5ad").exists()

    def test_tile_size_flag_reaches_segmenter(self, cli_sample, fake_cellpose, monkeypatch):
        """`--tile-size` / `--overlap` 必須真的傳到分割器。

        兩個旗標原本被 argparse 收下卻從未使用（`step_segment` 硬編碼 1024/128），
        使用者調整後**沒有任何效果也不報錯**。
        """
        from backend.src.cli import segment as cli
        from backend.src.segmentation import cellpose_runner

        seen = {}
        real = cellpose_runner.run_tiled_mcseg_v2

        def spy(*args, **kwargs):
            seen.update(tile_size=kwargs.get("tile_size"), overlap=kwargs.get("overlap"))
            return real(*args, **kwargs)

        monkeypatch.setattr(cellpose_runner, "run_tiled_mcseg_v2", spy)

        cli.main([
            "--btf", str(cli_sample["btf"]),
            "--out", str(cli_sample["out"]),
            "--no-gpu", "--skip-celltypist",
            "--tile-size", "256", "--overlap", "64",
        ])

        assert seen == {"tile_size": 256, "overlap": 64}


class TestCliFullSegSafety:
    """架構深化 P8：CLI 全片路徑接上與 API 共用的 OOM 防護與 MPS batch_size 鉗制"""

    def test_max_load_gb_rejects_oversized_output(self, cli_sample, fake_cellpose):
        """裁切窗格過大（超過 --max-load-gb）須明確報錯，而非默默 OOM。

        `main()` 目前沒有頂層 try/except（與 P8 無關的既有行為，`MemoryError`
        跟其他未捕捉例外一樣會直接冒出，不是回傳非 0 值），故這裡斷言例外本身。
        """
        import pytest

        from backend.src.cli.segment import main

        with pytest.raises(MemoryError, match="請縮小"):
            main([
                "--btf", str(cli_sample["btf"]),
                "--out", str(cli_sample["out"]),
                "--no-gpu", "--skip-celltypist",
                "--tile-size", "256", "--overlap", "64",
                "--max-load-gb", "1e-9",
            ])

        assert not (cli_sample["out"] / "mcseg_mask.npy").exists()

    def test_batch_size_clamped_to_two(self, cli_sample, fake_cellpose, monkeypatch):
        """`--batch-size` 高於 2 時須被鉗制到 2（docs/adr/0005：無條件套用）。"""
        from backend.src.cli import segment as cli
        from backend.src.segmentation import cellpose_runner

        seen = {}
        real = cellpose_runner.run_tiled_mcseg_v2

        def spy(img=None, cfg=None, **kwargs):
            seen["batch_size"] = (cfg or {}).get("batch_size")
            return real(img, cfg, **kwargs)

        monkeypatch.setattr(cellpose_runner, "run_tiled_mcseg_v2", spy)

        cli.main([
            "--btf", str(cli_sample["btf"]),
            "--out", str(cli_sample["out"]),
            "--no-gpu", "--skip-celltypist",
            "--tile-size", "256", "--overlap", "64",
            "--batch-size", "8",
        ])

        assert seen["batch_size"] == 2

    def test_cpsam_flag_still_reaches_segmenter(self, cli_sample, fake_cellpose, monkeypatch):
        """`--cpsam` 仍必須生效——CLI 傳 force_disable_cpsam=False，不受 Web UI

        全片路徑「預設強制關閉 cpsam」的行為影響（兩者的裁切窗格風險不同，
        見 `apply_full_seg_safety_clamp` docstring）。
        """
        from backend.src.cli import segment as cli
        from backend.src.segmentation import cellpose_runner

        seen = {}
        real = cellpose_runner.run_tiled_mcseg_v2

        def spy(img=None, cfg=None, **kwargs):
            seen["use_cpsam"] = (cfg or {}).get("use_cpsam")
            return real(img, cfg, **kwargs)

        monkeypatch.setattr(cellpose_runner, "run_tiled_mcseg_v2", spy)

        cli.main([
            "--btf", str(cli_sample["btf"]),
            "--out", str(cli_sample["out"]),
            "--no-gpu", "--skip-celltypist",
            "--tile-size", "256", "--overlap", "64",
            "--cpsam",
        ])

        assert seen["use_cpsam"] is True


class TestCliExportXenium:
    """`step_export_xenium` 產生的多邊形必須與 Web UI（`export/geometry.py`）同一套過濾規則"""

    def test_filters_small_noise_polygons(self, tmp_path, monkeypatch):
        """單像素雜訊細胞（面積 < 20px）不得出現在 Xenium 匯出的多邊形裡。"""
        from backend.src.cli.segment import step_export_xenium
        from backend.src.export import xenium_exporter

        class _FakeExporter:
            def __init__(self, **kwargs):
                pass

            def export(self, cells_h5ad_path, xen_dir):
                Path(xen_dir).mkdir(parents=True, exist_ok=True)

        monkeypatch.setattr(xenium_exporter, "XeniumExporter", _FakeExporter)

        mask = np.zeros((40, 40), dtype=np.int32)
        mask[5:15, 5:15] = 1       # 大細胞：100 px
        mask[30:32, 30:32] = 2     # 雜訊：4 px < min_area_px=20

        cells_h5ad_path = tmp_path / "cells.h5ad"
        cells_h5ad_path.write_bytes(b"")  # 內容不重要，export() 已被假掉

        step_export_xenium(mask, cells_h5ad_path, tmp_path, pixel_size_um=0.5, he_image_path=None)

        geo = json.loads((tmp_path / "cells_polygons.geojson").read_text(encoding="utf-8"))
        assert len(geo["features"]) == 1
        assert geo["features"][0]["properties"]["cell_id"] == 1


class TestCliLayerIsThin:
    """回歸：GeoJSON 生成邏輯只能有一份（`export/geometry.py`）"""

    def test_no_geojson_generation_logic_duplicated(self):
        src = (Path(__file__).resolve().parents[1] / "src" / "cli" / "segment.py").read_text(
            encoding="utf-8"
        )
        for gone in ("find_contours", "measure.regionprops"):
            assert gone not in src, f"{gone} 不應再出現在 cli/segment.py（改呼叫 export.geometry）"


class TestCliTissuePresets:
    """CLI 的組織參數必須與 Web UI 同源（`config/profiles/`）

    CLI 原本硬編碼一份 `TISSUE_PRESETS`，與 profile 各走一套而**靜默漂移**：
    CRC `voronoi_distance` 8 vs 9、LUAD 有 4 個參數不一致（且 voronoi 在兩個組織
    是反方向差異）。後果是同一個 `crc` preset 在 CLI 與 UI 產出不同遮罩，
    而 README 卻聲稱兩者「use exactly the same engine with consistent
    parameter semantics」。
    """

    PARAM_KEYS = ("dia_small", "dia_mid", "dia_large", "voronoi_distance",
                  "clahe_clip_limit", "flow_threshold", "cellprob_threshold",
                  "min_size", "max_size", "use_hematoxylin")

    def test_no_hardcoded_preset_dict_remains(self):
        """CLI 不得再有硬編碼的參數表（CLAUDE.md §15）。"""
        from backend.src.cli import segment as cli

        assert not hasattr(cli, "TISSUE_PRESETS"), (
            "硬編碼參數表已被 config/profiles/ 取代，重新引入會再次漂移"
        )

    @pytest.mark.parametrize("tissue", ["crc", "luad", "default"])
    def test_preset_matches_profile_yaml(self, tissue):
        """CLI 取到的每個參數都必須等於 profile YAML 的值。"""
        import yaml

        from backend.src.cli.segment import load_tissue_preset

        raw = yaml.safe_load(
            Path(f"config/profiles/{tissue}.yaml").read_text(encoding="utf-8")
        )
        expected = raw["segmentation"]["mcseg_v2"]
        actual = load_tissue_preset(tissue)

        for key in self.PARAM_KEYS:
            assert actual.get(key) == expected.get(key), (
                f"{tissue}.{key}：CLI 取到 {actual.get(key)!r}，"
                f"profile 是 {expected.get(key)!r}"
            )

    def test_available_tissues_ignores_appledouble_files(self, tmp_path, monkeypatch):
        """`._*` AppleDouble 檔不得變成一個幽靈組織選項（CLAUDE.md §4）。

        實測過：在 ExFAT 上寫入 `default.yaml` 會連帶產生 `._default.yaml`，
        使 `--tissue` 的 choices 多出 `._default`。
        """
        from backend.src.cli import segment as cli
        from backend.src.utils import config as cfg_mod

        (tmp_path / "crc.yaml").write_text("x: 1", encoding="utf-8")
        (tmp_path / "._crc.yaml").write_bytes(b"\x00\x05\x16\x07")
        monkeypatch.setattr(cfg_mod, "_PROFILES_DIR", tmp_path)

        assert cli.available_tissues() == ["crc"]

    def test_missing_mcseg_block_fails_loudly(self, tmp_path, monkeypatch):
        """profile 缺 `segmentation.mcseg_v2` 時必須中止，不得靜默跑空參數。

        空參數等同 cellpose 預設值 —— 使用者會拿到一份看似正常、參數卻完全
        不對的遮罩，比直接失敗糟得多。
        """
        from backend.src.cli.segment import load_tissue_preset
        from backend.src.utils import config as cfg_mod

        (tmp_path / "broken.yaml").write_text("profile_name: X", encoding="utf-8")
        monkeypatch.setattr(cfg_mod, "_PROFILES_DIR", tmp_path)

        with pytest.raises(SystemExit, match="缺少 segmentation.mcseg_v2"):
            load_tissue_preset("broken")

    def test_cli_flags_still_override_profile(self, cli_sample, fake_cellpose, monkeypatch):
        """`--voronoi-d` 等旗標仍須蓋過 profile —— 論文基準用 d=8 就靠這條路。"""
        from backend.src.cli import segment as cli
        from backend.src.segmentation import cellpose_runner

        seen = {}
        real = cellpose_runner.run_tiled_mcseg_v2

        def spy(img=None, cfg=None, **kwargs):
            seen.update(cfg or {})
            return real(img, cfg, **kwargs)

        monkeypatch.setattr(cellpose_runner, "run_tiled_mcseg_v2", spy)

        cli.main([
            "--btf", str(cli_sample["btf"]),
            "--out", str(cli_sample["out"]),
            "--tissue", "crc", "--voronoi-d", "8", "--dia-mid", "19",
            "--no-gpu", "--skip-celltypist",
            "--tile-size", "256", "--overlap", "64",
        ])

        assert seen["voronoi_distance"] == 8, "旗標未蓋過 profile 的 9"
        assert seen["dia_mid"] == 19.0
        # 未指定的參數仍應來自 profile
        assert seen["clahe_clip_limit"] == 3.0


class TestCliBinToImageTransform:
    """CLI 必須與 API 共用「對位 JSON 優先」的座標變換

    CLI 原本停在 `resolve_bin_to_mask_scale`（近似分軸縮放）—— P0.5 已判定它在
    SR 畫布相對影像有 padding 時**幾何上是錯的**（把 bin 橫向壓縮進影像寬度）。
    dpcp01 實測該誤差在右緣達約 875 px，足以讓整個 800px 窗格的 bin 全部落空。
    """

    def test_cli_uses_alignment_json_homography(self, cli_sample, monkeypatch):
        """spatial/ 內有可組合的 old/new JSON 時，CLI 須取 homography 而非近似縮放。"""
        from backend.src.cli.segment import _load_cli_config
        from backend.src.fullslide.pipeline import resolve_bin_to_image_transform

        # _load_cli_config() 刻意會讀真實 state.json 的 alignment 設定（CLI 設計上
        # 就是要沿用使用者透過 UI「指定對位 JSON」存的值）——這裡要測的是「spatial/
        # 內自動偵測」這條路徑，必須把即時環境的 extra_alignment_json 隔離掉，
        # 否則使用者若真的透過 UI 設定過（例如指到另一個樣本的對位檔），這個測試
        # 會撿到不屬於本測試的真實檔案而報錯（2026-08-05 實際發生過一次：
        # state.json 指到有勝樣本的 H1-KRDKFYF-...json，serialNumber 跟本測試
        # 自己造的 old.json 對不上）。
        monkeypatch.setattr("backend.src.utils.config.load_config", lambda *a, **k: {"alignment": {}})

        binned = cli_sample["binned"]
        sp = _write_spatial(binned, hires_size=(2870, 6000), scalef=0.25475544, mpp=0.5464)
        _write_align(sp / "old.json", scale_transform=0.2, scale_images=0.2 * 0.5464)
        _write_align(sp / "new.json", scale_transform=0.2, scale_images=0.2 * 0.2732)

        cfg = _load_cli_config()
        cfg.setdefault("paths", {})["binned_002"] = str(binned)

        transform, _, source = resolve_bin_to_image_transform(cfg, (47104, 21504))

        assert transform is not None, "有對位 JSON 卻回退到近似縮放"
        np.testing.assert_allclose(transform, np.diag([2.0, 2.0, 1.0]), rtol=1e-6)
        assert "old.json" in source and "new.json" in source

    def test_cli_config_does_not_inherit_config_file_paths(self, tmp_path):
        """`_load_cli_config` 不得帶入設定檔的 `paths:`。

        CLI 的樣本由引數決定；`config/pipeline.yaml` 與 `state.json` 指向的
        很可能是另一個樣本，混進來會靜默讀錯樣本的 spatial/ 與對位 JSON。
        """
        from backend.src.cli.segment import _load_cli_config

        cfg = _load_cli_config()

        assert "paths" not in cfg
        assert set(cfg) <= {"alignment"}

    def test_crop_derives_transform_from_source_image_not_mask(
        self, cli_sample, fake_cellpose, stub_h5_reader, monkeypatch
    ):
        """裁切時變換須以**來源影像**尺寸推導，與不裁切時一致。

        這是 commit `a338f3e` 在 API 端修掉的同一個缺陷：拿裁切後的遮罩當畫布，
        800×800 窗格會算出 scale (0.071, 0.034)。
        """
        from backend.src.cli import segment as cli

        from backend.src.fullslide import pipeline

        seen = []
        real_fn = pipeline.resolve_bin_to_image_transform

        def spy(config, image_shape):
            seen.append(tuple(image_shape))
            return real_fn(config, image_shape)

        monkeypatch.setattr(pipeline, "resolve_bin_to_image_transform", spy)

        cli.main([
            "--btf", str(cli_sample["btf"]),
            "--crop-y0", "64", "--crop-y1", "320",
            "--btf-col0", "64", "--btf-col1", "320",
            "--tp", str(cli_sample["tp"]),
            "--h5", str(cli_sample["h5"]),
            "--out", str(cli_sample["out"]),
            "--no-gpu", "--skip-celltypist",
            "--tile-size", "256", "--overlap", "64",
        ])

        mask = np.load(str(cli_sample["out"] / "mcseg_mask.npy"))
        assert mask.shape == (256, 256), "遮罩應為裁切窗格尺寸"
        # 但推導變換用的畫布必須是整張影像
        assert seen == [(512, 512)]

    def test_zero_matched_bins_raises_readable_error(self, cli_sample, stub_h5_reader):
        """沒有 bin 對到細胞時須給可讀訊息，不得是空陣列取 max 的 numpy 錯誤。

        座標系錯配（真實案例：近似縮放讓右緣整個窗格的 bin 全部落空）會走到這裡。
        """
        import pandas as pd

        from backend.src.fullslide.pipeline import aggregate_cells

        empty = pd.DataFrame({"barcode": pd.Series([], dtype=str),
                              "cell_id": pd.Series([], dtype="int32")})

        with pytest.raises(ValueError, match="沒有任何 bin 對應到細胞"):
            aggregate_cells(empty, cli_sample["h5"])
