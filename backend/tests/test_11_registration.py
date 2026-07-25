"""Test 11: 對位配準模組（backend/src/registration/）

以合成資料為主；需要真實 dpcp01 對位 JSON 的測試一律 `skipif` 檔案不存在，
確保 Windows / 換樣本時仍可全綠。
"""
import json
from pathlib import Path

import numpy as np
import pytest

# 真實 dpcp01 對位檔（低解析 old + 高解析 new），供組合 homography 的端到端驗證
DPCP01_SPATIAL = Path(
    "/Volumes/SSD/plan_a/tissue sample/raw/binned_outputs/square_002um/spatial"
)
DPCP01_OLD = DPCP01_SPATIAL / "H1-WGR3TC4-D1-fiducials-image-registration.json"
DPCP01_NEW = DPCP01_SPATIAL / "H1-WGR3TC4-D1-fiducials-image-registration_0105.json"


# ── 合成對位 JSON ────────────────────────────────────────────────────────────

def _write_alignment_json(
    path,
    *,
    scale_transform: float = 0.2,
    scale_images: float = 0.1,
    serial_number: str = "H1-TEST01",
    area: str = "D1",
    checksum: str = "deadbeef",
):
    """寫出結構與 Loupe/CytAssist 對位檔相同的最小合成 JSON。"""
    doc = {
        "serialNumber": serial_number,
        "area": area,
        "checksum": checksum,
        "transform": [scale_transform, 0.0, 0.0,
                      0.0, scale_transform, 0.0,
                      0.0, 0.0, 1.0],
        "cytAssistInfo": {
            "transformImages": [scale_images, 0.0, 0.0,
                                0.0, scale_images, 0.0,
                                0.0, 0.0, 1.0],
        },
    }
    Path(path).write_text(json.dumps(doc), encoding="utf-8")
    return path


# ── load_alignment ──────────────────────────────────────────────────────────

class TestLoadAlignment:
    """Loupe / CytAssist 對位 JSON 解析"""

    def test_load_alignment_derives_mpp(self, tmp_path):
        """mpp = scale(transformImages) / scale(transform)。"""
        from backend.src.registration.alignment import load_alignment

        p = _write_alignment_json(
            tmp_path / "a.json", scale_transform=0.2, scale_images=0.1
        )
        al = load_alignment(p)

        assert al.mpp == pytest.approx(0.5, rel=1e-9)
        assert al.serial_number == "H1-TEST01"
        assert al.area == "D1"
        assert al.transform.shape == (3, 3)
        assert al.transform_images.shape == (3, 3)

    def test_load_alignment_rejects_missing_key(self, tmp_path):
        """缺 `cytAssistInfo.transformImages` 時報錯，訊息含檔名但不含完整路徑。"""
        from backend.src.registration.alignment import load_alignment

        p = tmp_path / "broken.json"
        p.write_text(json.dumps({"transform": [1, 0, 0, 0, 1, 0, 0, 0, 1]}), encoding="utf-8")

        with pytest.raises(ValueError) as exc:
            load_alignment(p)
        assert "broken.json" in str(exc.value)
        assert str(tmp_path) not in str(exc.value)


# ── pick_source_alignment ───────────────────────────────────────────────────

class TestPickSourceAlignment:
    """依 scalefactors 的 microns_per_pixel 自動選出 H_old"""

    def test_pick_source_alignment_matches_scalefactors_mpp(self, tmp_path):
        """同 serial 的兩個版本 scale 差 2 倍時，選中與 scalefactors 相符者。"""
        from backend.src.registration.alignment import pick_source_alignment

        low = _write_alignment_json(
            tmp_path / "low.json", scale_transform=0.2, scale_images=0.2 * 0.5465
        )
        high = _write_alignment_json(
            tmp_path / "high.json", scale_transform=0.2, scale_images=0.2 * 0.2732
        )

        source, others = pick_source_alignment([low, high], target_mpp=0.5464)

        assert source.path.name == "low.json"
        assert [o.path.name for o in others] == ["high.json"]

    def test_pick_source_alignment_skips_unparsable(self, tmp_path):
        """`spatial/` 內混有非對位 JSON（如 scalefactors_json.json）時應略過。"""
        from backend.src.registration.alignment import pick_source_alignment

        (tmp_path / "scalefactors_json.json").write_text(
            json.dumps({"microns_per_pixel": 0.5464}), encoding="utf-8"
        )
        good = _write_alignment_json(
            tmp_path / "good.json", scale_transform=0.2, scale_images=0.2 * 0.5465
        )

        source, others = pick_source_alignment(
            sorted(tmp_path.glob("*.json")), target_mpp=0.5464
        )
        assert source.path == good
        assert others == []

    def test_pick_source_alignment_raises_with_candidate_mpps(self, tmp_path):
        """無命中時報錯並列出各候選 mpp，讓使用者判斷該補哪份檔。"""
        from backend.src.registration.alignment import pick_source_alignment

        p = _write_alignment_json(
            tmp_path / "only.json", scale_transform=0.2, scale_images=0.2 * 0.2732
        )
        with pytest.raises(ValueError) as exc:
            pick_source_alignment([p], target_mpp=0.5464)
        assert "0.2732" in str(exc.value)


# ── provenance 驗證 ─────────────────────────────────────────────────────────

class TestValidatePair:
    """拒絕把別片玻片的對位檔組合在一起"""

    def test_validate_pair_accepts_same_slide(self, tmp_path):
        from backend.src.registration.alignment import load_alignment, validate_pair

        a = load_alignment(_write_alignment_json(tmp_path / "a.json"))
        b = load_alignment(
            _write_alignment_json(tmp_path / "b.json", scale_images=0.05)
        )
        validate_pair(a, b)   # 不應 raise

    def test_validate_pair_rejects_serial_mismatch(self, tmp_path):
        from backend.src.registration.alignment import load_alignment, validate_pair

        a = load_alignment(
            _write_alignment_json(tmp_path / "a.json", serial_number="H1-AAA")
        )
        b = load_alignment(
            _write_alignment_json(tmp_path / "b.json", serial_number="H1-BBB")
        )
        with pytest.raises(ValueError) as exc:
            validate_pair(a, b)
        assert "H1-AAA" in str(exc.value) and "H1-BBB" in str(exc.value)

    def test_validate_pair_rejects_area_mismatch(self, tmp_path):
        from backend.src.registration.alignment import load_alignment, validate_pair

        a = load_alignment(_write_alignment_json(tmp_path / "a.json", area="D1"))
        b = load_alignment(_write_alignment_json(tmp_path / "b.json", area="A1"))
        with pytest.raises(ValueError):
            validate_pair(a, b)


# ── compose_bin_to_image ────────────────────────────────────────────────────

class TestComposeBinToImage:
    """組合 homography：SR fullres px → 新影像 px"""

    def test_compose_synthetic_halved_mpp_gives_2x(self, tmp_path):
        """新影像解析度為舊的一半 µm/px → 座標放大 2 倍。"""
        from backend.src.registration.alignment import (
            compose_bin_to_image,
            load_alignment,
        )

        old = load_alignment(
            _write_alignment_json(tmp_path / "old.json", scale_images=0.1)
        )
        new = load_alignment(
            _write_alignment_json(tmp_path / "new.json", scale_images=0.05)
        )
        m = compose_bin_to_image(old, new)

        np.testing.assert_allclose(m, np.diag([2.0, 2.0, 1.0]), atol=1e-9)

    @pytest.mark.skipif(
        not (DPCP01_OLD.exists() and DPCP01_NEW.exists()),
        reason="需要 dpcp01 真實對位 JSON",
    )
    def test_compose_recovers_isotropic_2x(self):
        """真實 dpcp01 兩份對位檔 → 等向 2×、幾乎無旋轉。"""
        from backend.src.registration.alignment import (
            compose_bin_to_image,
            load_alignment,
            validate_pair,
        )

        old = load_alignment(DPCP01_OLD)
        new = load_alignment(DPCP01_NEW)
        validate_pair(old, new)
        m = compose_bin_to_image(old, new)

        sx = np.hypot(m[0, 0], m[1, 0])
        sy = np.hypot(m[0, 1], m[1, 1])
        rot_deg = np.degrees(np.arctan2(m[1, 0], m[0, 0]))

        assert sx == pytest.approx(2.0, abs=1e-3)
        assert sy == pytest.approx(2.0, abs=1e-3)
        assert abs(rot_deg) < 0.01
        # 平移應為個位數 px（實測 ≈ (−0.90, −0.85)）
        assert abs(m[0, 2]) < 5 and abs(m[1, 2]) < 5
