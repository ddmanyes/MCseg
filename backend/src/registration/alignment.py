"""
Loupe / CytAssist 對位 JSON 解析與組合
=====================================

**解決什麼問題**：Space Ranger 的 `pxl_*_in_fullres` 校準到**餵給它的那張影像**。
當使用者因「outs 內的高解析圖不夠清晰」而另外輸出一張更高解析的圖、再用
Loupe Browser 重新對位時，RNA 座標與新影像之間就差一個變換。本模組讀取
兩份對位 JSON 並組出該變換，**不需重跑 Space Ranger**。

JSON 結構（`*_alignment_file.json` 與 `*-fiducials-image-registration.json`
schema 相同），各帶兩個 3×3 homography：

| 欄位 | 映射方向 |
|------|---------|
| `transform` | slide µm → CytAssist px |
| `cytAssistInfo.transformImages` | **影像 px → CytAssist px** |

由此推導該影像的實際解析度：

```text
mpp_image = scale(transformImages) / scale(transform)
```

（dpcp01 實測：無後綴 → 0.5465，與 `scalefactors_json.json` 的
`microns_per_pixel` 0.5464 相符；`_0105` → 0.2732 ＝ 47104×21504 高解析 TIFF。）

再組合兩份即得 bin → 新影像的變換：

```text
新影像 px = inv(H_new.transform_images) @ H_old.transform_images @ (SR fullres px)
```

dpcp01 實測組合結果為等向 scale 2.00001、rot −0.0001°、平移 (−0.90, −0.85) px；
相較舊有的 `resolve_bin_to_mask_scale`（hires 尺寸 ÷ scalef，在 SR 畫布相對
影像有 padding 時會把 bin 橫向壓縮）多命中 13.4% bins / 17.5% 細胞。
"""
from __future__ import annotations

import json
import logging
from dataclasses import dataclass
from pathlib import Path

import numpy as np

logger = logging.getLogger("pipeline.registration.alignment")


@dataclass(frozen=True)
class Alignment:
    """一份對位 JSON 的解析結果。"""

    transform: np.ndarray          # 3×3，slide µm → CytAssist px
    transform_images: np.ndarray   # 3×3，影像 px → CytAssist px
    mpp: float                     # 該影像的 µm/px（由兩個 scale 相除推導）
    serial_number: str
    area: str
    checksum: str
    path: Path

    @property
    def name(self) -> str:
        return self.path.name


def _as_3x3(values, field: str, filename: str) -> np.ndarray:
    """把 JSON 的 9 元素（或 3×3 巢狀）陣列轉為 3×3，並正規化 M[2,2]=1。"""
    arr = np.asarray(values, dtype=float)
    if arr.size != 9:
        raise ValueError(f"{filename}：`{field}` 需為 9 個元素的 3×3 矩陣，實得 {arr.size}")
    m = arr.reshape(3, 3)
    if m[2, 2] == 0 or not np.isfinite(m).all():
        raise ValueError(f"{filename}：`{field}` 不是有效的 homography")
    return m / m[2, 2]


def matrix_scale(m: np.ndarray) -> float:
    """homography 線性部分的等效等向縮放（行列式開根號）。"""
    return float(np.sqrt(abs(np.linalg.det(np.asarray(m, dtype=float)[:2, :2]))))


def load_alignment(path: str | Path) -> Alignment:
    """
    解析一份 Loupe / CytAssist 對位 JSON。

    Raises
    ------
    ValueError
        檔案損壞或缺少必要欄位。訊息含**檔名**但不含完整路徑（可安全回傳前端）。
    """
    path = Path(path)
    try:
        doc = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as e:
        raise ValueError(f"{path.name}：無法讀取或解析對位 JSON（{type(e).__name__}）") from e

    if not isinstance(doc, dict):
        raise ValueError(f"{path.name}：對位 JSON 頂層需為物件")

    if "transform" not in doc:
        raise ValueError(f"{path.name}：對位 JSON 缺少 `transform`")
    cyt = doc.get("cytAssistInfo")
    if not isinstance(cyt, dict) or "transformImages" not in cyt:
        raise ValueError(f"{path.name}：對位 JSON 缺少 `cytAssistInfo.transformImages`")

    transform = _as_3x3(doc["transform"], "transform", path.name)
    transform_images = _as_3x3(cyt["transformImages"], "cytAssistInfo.transformImages", path.name)

    s_slide = matrix_scale(transform)
    if s_slide <= 0:
        raise ValueError(f"{path.name}：`transform` 的縮放為 0，無法推導 µm/px")

    return Alignment(
        transform=transform,
        transform_images=transform_images,
        mpp=matrix_scale(transform_images) / s_slide,
        serial_number=str(doc.get("serialNumber", "")),
        area=str(doc.get("area", "")),
        checksum=str(doc.get("checksum", "")),
        path=path,
    )


def pick_source_alignment(
    paths, target_mpp: float, tol: float = 0.01
) -> tuple[Alignment, list[Alignment]]:
    """
    從候選對位檔中選出 `H_old` —— 即**產生 `pxl_*_in_fullres` 的那張影像**的對位檔。

    判準：推導 mpp 與 `scalefactors_json.json` 的 `microns_per_pixel` 相對誤差
    < `tol`。這同時解掉兩個實際踩過的陷阱：

    1. 把**別片玻片**的註冊檔套進來（mpp 通常不會恰好吻合）；
    2. 同一 serial 存在兩個版本、scale 正好差 2 倍（dpcp01 的 0.5465 vs 0.2732）
       —— 靠檔名或修改時間都選不對。

    無法解析的 JSON（例如 `spatial/` 內的 `scalefactors_json.json`）一律略過，
    不中斷流程。

    Returns
    -------
    tuple[Alignment, list[Alignment]]
        `(source, others)`；`others` 為其餘可解析的候選（保留輸入順序），
        新影像的 `H_new` 通常就在其中。

    Raises
    ------
    ValueError
        無任何候選命中；訊息會列出各候選的推導 mpp 供人工判斷。
    """
    candidates: list[Alignment] = []
    for p in paths:
        try:
            candidates.append(load_alignment(p))
        except ValueError:
            continue   # 非對位 JSON，略過（CLAUDE.md §11 entry 層級容錯）

    if not candidates:
        raise ValueError("找不到任何可解析的對位 JSON（Loupe/CytAssist）")

    best, best_err = None, float("inf")
    for al in candidates:
        err = abs(al.mpp - target_mpp) / target_mpp if target_mpp else float("inf")
        if err < best_err:
            best, best_err = al, err

    if best is None or best_err >= tol:
        detail = "、".join(f"{al.name} → {al.mpp:.4f}" for al in candidates)
        raise ValueError(
            f"沒有對位檔的推導 µm/px 與 scalefactors 的 {target_mpp:.4f} 相符"
            f"（容差 {tol:.1%}）。候選：{detail}"
        )

    others = [al for al in candidates if al.path != best.path]
    logger.info(
        f"選定來源對位檔 {best.name}（mpp {best.mpp:.4f}，"
        f"與 scalefactors {target_mpp:.4f} 差 {best_err:.2%}）"
    )
    return best, others


def validate_pair(h_old: Alignment, h_new: Alignment) -> None:
    """
    確認兩份對位檔出自**同一片玻片的同一區**，否則組合出的變換毫無意義。

    `serialNumber` / `area` 不一致即 raise —— 把別片玻片的對位檔套進來會讓
    RNA 整體錯位卻**不報錯**，是最難察覺的一類故障。

    Raises
    ------
    ValueError
        provenance 不符；訊息含兩邊的值。
    """
    mismatches = []
    if h_old.serial_number != h_new.serial_number:
        mismatches.append(
            f"serialNumber（{h_old.name}={h_old.serial_number!r} vs "
            f"{h_new.name}={h_new.serial_number!r}）"
        )
    if h_old.area != h_new.area:
        mismatches.append(
            f"area（{h_old.name}={h_old.area!r} vs {h_new.name}={h_new.area!r}）"
        )
    if mismatches:
        raise ValueError(
            "兩份對位檔不屬於同一片玻片，拒絕組合：" + "；".join(mismatches)
        )


def compose_bin_to_image(h_old: Alignment, h_new: Alignment) -> np.ndarray:
    """
    組合出「Space Ranger fullres px → 新影像 px」的 3×3 homography。

    ```text
    新影像 px = inv(H_new.transform_images) @ H_old.transform_images @ (SR fullres px)
    ```

    兩者都以 CytAssist 影像為中介座標系，因此無需重跑 Space Ranger。

    Notes
    -----
    `transformImages` 內含的 rot +90° 與鏡射，Space Ranger 在產生
    `pxl_*_in_fullres` 時**已經套入**，兩份相除會互相抵銷 —— 這正是既有
    scale-only 修正能重現 EP 數字的原因。但若有人拿方向不同的影像來分割，
    scale-only 會靜默失敗，而組合式會如實反映該旋轉。
    """
    m = np.linalg.inv(h_new.transform_images) @ h_old.transform_images
    if m[2, 2] == 0 or not np.isfinite(m).all():
        raise ValueError("組合後的 homography 無效（transformImages 可能奇異）")
    return m / m[2, 2]
