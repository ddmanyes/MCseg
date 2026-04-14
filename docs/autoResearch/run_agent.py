"""
run_agent.py — AutoResearch Agent Runner

Iteratively calls Claude Opus to propose, implement, and evaluate cell segmentation
architectures. Adapted from the AutoResearch paradigm (Karpathy 2026).

Usage:
    export ANTHROPIC_API_KEY=your_key_here
    uv run python run_agent.py

Requirements:
    uv add anthropic cellpose scikit-image opencv-python scipy numpy pandas
"""

import json
import time
import shutil
import subprocess
from pathlib import Path
from datetime import datetime

import anthropic

# ── Configuration (edit these for your project) ───────────────────────────────
ROOT         = Path(__file__).parent
SEGMENT_PATH = ROOT / "segment.py"           # the sandbox script the agent modifies
BEST_SEGMENT = ROOT / "segment_best.py"      # best-scoring version (auto-updated)
LOG_PATH     = ROOT / "results" / "experiment_log.jsonl"
MEMORY_PATH  = ROOT / "memory.md"
HISTORY_CSV  = ROOT / "results" / "history_score.csv"
PYTHON_BIN   = "python"                      # replace with path to your venv python
TARGET_SCORE = 0.70                          # stop when this AP@0.5 is reached
MAX_RUNTIME  = 360                           # seconds per experiment (timeout)

# ── System prompt (describes the task to Claude) ─────────────────────────────
# Edit this to match your tissue type, image specs, and available tools.
# See program.md for a human-readable version of the task specification.
SYSTEM_PROMPT = """You are an expert computational biologist and computer vision researcher.
Your task: write improved Python code for `segment.py` that maximises AP@0.5 for
cell instance segmentation on H&E images, evaluated against ground truth masks.

## Scientific context
- Image: H&E RGB numpy array, dtype uint8
- Pixel size: 0.2737 µm/px (adapt to your dataset)
- Ground truth: cytoplasm-inclusive cell boundaries (~5 µm radius beyond nucleus)
- Time budget: < 300 seconds total runtime per experiment

## Available packages
cellpose, scikit-image, opencv-python, scipy, pandas, numpy

## Cellpose v4 API notes
- `channels` parameter is deprecated — pass the image directly
- Default model is `cpsam` (Cellpose-SAM), not cyto3
- After eval with resample=True: always clip mask to image size:
  mask = mask[:img.shape[0], :img.shape[1]]

## Known effective strategies
- CLAHE contrast enhancement (clip=3.0, tile=8) consistently helps
- Voronoi-constrained expansion outperforms expand_labels
- Multi-diameter cyto3 ensemble (13/17/22 px) captures size heterogeneity
- cpsam adds complementary detections at small cell sizes

## Forbidden (causes timeout or errors)
- cv2.watershed, skimage.watershed: always exceeds time budget
- Per-cell morphological loops: O(N_cells) → timeout

## Required function signature (preserve exactly)
```python
def build_and_predict(img, vhd_csv, gt_mask=None) -> np.ndarray:
    \"\"\"Returns 2D int32 array: 0=background, >0=cell instance ID\"\"\"
    ...
    return final_mask
```

## Required __main__ block (preserve exactly)
```python
if __name__ == "__main__":
    import prepare
    img, gt_mask, vhd_csv = prepare.load_data()
    pred_mask = build_and_predict(img, vhd_csv, gt_mask=gt_mask)
    score = prepare.evaluate_iou(pred_mask, gt_mask)
```

## Response format (strict)
Line 1: APPROACH: [one-line strategy description]
Line 2: RATIONALE: [2-3 sentences of scientific reasoning]
Line 3+: Complete Python code starting with import statements. No markdown fences.
Your response MUST end with exactly:
    score = prepare.evaluate_iou(pred_mask, gt_mask)
Output nothing after this line."""


# ── Agent loop utilities ──────────────────────────────────────────────────────

def load_history() -> list[dict]:
    if not LOG_PATH.exists():
        return []
    records = []
    with open(LOG_PATH) as f:
        for line in f:
            line = line.strip()
            if line:
                try:
                    records.append(json.loads(line))
                except json.JSONDecodeError:
                    pass
    return records


def build_user_message(history: list[dict], memory_text: str, reference_code: str) -> str:
    valid = [r for r in history if r.get("score", -1) >= 0]
    best_score = max((r["score"] for r in valid), default=0.0)

    if not valid:
        history_section = "No experiments yet — start fresh."
    else:
        top = sorted(valid, key=lambda r: r["score"], reverse=True)[:8]
        history_section = f"Total experiments: {len(valid)} | Best so far: {best_score:.4f}\n\n"
        history_section += "Top results:\n"
        for r in top:
            history_section += f"  {r['score']:.4f}  {r.get('approach', 'N/A')}\n"
        failed = [r for r in valid if r["score"] < 0.35]
        if failed:
            history_section += "\nFailed approaches (avoid):\n"
            for r in sorted(failed, key=lambda r: r["score"])[:5]:
                history_section += f"  {r['score']:.4f}  {r.get('approach', 'N/A')}\n"

    return f"""Current best AP@0.5: {best_score:.4f}
Target: > {TARGET_SCORE}

## Experiment History
{history_section}

## Research Memory (key insights)
{memory_text}

## Current segment.py (your starting point)
```python
{reference_code}
```

Design the best segment.py you can. Bold architectural changes may be needed to beat
the current best — don't just tweak parameters.

Output APPROACH:, RATIONALE:, then the full Python code."""


def generate_segment(history: list[dict], memory_text: str, reference_code: str) -> str:
    client = anthropic.Anthropic()
    resp = client.messages.create(
        model="claude-opus-4-6",
        max_tokens=8000,
        system=SYSTEM_PROMPT,
        messages=[{"role": "user", "content": build_user_message(history, memory_text, reference_code)}],
    )
    return resp.content[0].text


def parse_response(text: str) -> tuple[str, str, str]:
    lines = text.strip().splitlines()
    approach, rationale, code_lines, in_code = "Unknown", "", [], False
    for line in lines:
        if not in_code:
            if line.startswith("APPROACH:"):
                approach = line[9:].strip()
            elif line.startswith("RATIONALE:"):
                rationale = line[10:].strip()
            elif line.startswith(("import ", "from ", "#!/")):
                in_code = True
                code_lines.append(line)
        else:
            code_lines.append(line)
    ENDING = "    score = prepare.evaluate_iou(pred_mask, gt_mask)"
    code_str = "\n".join(code_lines)
    idx = code_str.rfind(ENDING)
    if idx != -1:
        code_str = code_str[:idx + len(ENDING)]
    return approach, rationale, code_str


def run_experiment(code: str) -> tuple[float | None, str, str]:
    try:
        compile(code, "<segment.py>", "exec")
    except SyntaxError as e:
        return None, "", f"SyntaxError: {e}"
    with open(SEGMENT_PATH, "w") as f:
        f.write(code)
    try:
        result = subprocess.run(
            [PYTHON_BIN, str(SEGMENT_PATH)],
            capture_output=True, text=True, timeout=MAX_RUNTIME,
        )
    except subprocess.TimeoutExpired:
        return None, "", f"TIMEOUT: exceeded {MAX_RUNTIME}s"
    if result.returncode != 0:
        return None, result.stdout, result.stderr[-800:]
    try:
        with open(HISTORY_CSV) as f:
            last_line = f.readlines()[-1].strip()
        score = float(last_line.split(",")[-1])
        return score, result.stdout, ""
    except Exception as e:
        return None, result.stdout, f"Score read failed: {e}"


def save_record(score: float | None, approach: str, rationale: str) -> None:
    LOG_PATH.parent.mkdir(parents=True, exist_ok=True)
    record = {"timestamp": datetime.now().isoformat(),
               "score": score if score is not None else -1.0,
               "approach": approach, "rationale": rationale}
    with open(LOG_PATH, "a") as f:
        f.write(json.dumps(record, ensure_ascii=False) + "\n")


def update_memory(history: list[dict]) -> None:
    client = anthropic.Anthropic()
    valid = [r for r in history if r.get("score", -1) >= 0]
    best = max(valid, key=lambda r: r["score"]) if valid else None
    top_summary = "\n".join(
        f"- score={r['score']:.4f}: {r.get('approach', 'N/A')}"
        for r in sorted(valid, key=lambda r: r["score"], reverse=True)[:12]
    )
    with open(MEMORY_PATH) as f:
        current = f.read()
    resp = client.messages.create(
        model="claude-haiku-4-5-20251001",
        max_tokens=1800,
        messages=[{"role": "user", "content": f"""Update this research memory file.

Current memory:
{current}

Top experiments:
{top_summary}

Best result: score={(best['score'] if best else 0):.4f}, approach={best.get('approach','N/A') if best else 'N/A'}

Instructions: update "Current best" and "What works" sections, add new insights,
remove outdated conclusions. Keep concise and actionable. Return only markdown."""}],
    )
    with open(MEMORY_PATH, "w") as f:
        f.write(resp.content[0].text)
    print("  memory.md updated")


# ── Main loop ─────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    # Initialise memory.md if it doesn't exist
    if not MEMORY_PATH.exists():
        MEMORY_PATH.write_text("# Research Memory\n\n## Current best\nNone yet.\n\n## What works\n- TBD\n\n## What fails\n- TBD\n")

    history = load_history()
    best_score = max((r.get("score", 0) for r in history if r.get("score", -1) >= 0), default=0.0)
    print(f"AutoResearch Agent | History: {len(history)} runs | Best: {best_score:.4f}")

    if not BEST_SEGMENT.exists() and SEGMENT_PATH.exists():
        shutil.copy(SEGMENT_PATH, BEST_SEGMENT)

    i = 0
    while True:
        i += 1
        print(f"\n--- Iteration {i} | Best: {best_score:.4f} ---")

        ref_path = BEST_SEGMENT if BEST_SEGMENT.exists() else SEGMENT_PATH
        reference_code = ref_path.read_text()
        memory_text = MEMORY_PATH.read_text()

        print("  Generating new architecture...")
        try:
            ai_text = generate_segment(history, memory_text, reference_code)
        except Exception as e:
            print(f"  API error: {e}")
            time.sleep(5)
            continue

        approach, rationale, new_code = parse_response(ai_text)
        print(f"  Approach: {approach}")

        if not new_code.strip():
            print("  Could not parse code, skipping")
            save_record(None, approach, "PARSE_FAILED")
            continue

        print(f"  Running experiment (max {MAX_RUNTIME}s)...")
        t0 = time.time()
        score, stdout, stderr = run_experiment(new_code)
        elapsed = time.time() - t0

        if score is None:
            print(f"  FAILED ({elapsed:.0f}s): {stderr[-200:]}")
            save_record(None, approach, f"FAILED | {stderr[:120]}")
            if BEST_SEGMENT.exists():
                shutil.copy(BEST_SEGMENT, SEGMENT_PATH)
            continue

        print(f"  Score: {score:.4f}  ({elapsed:.0f}s)")
        save_record(score, approach, rationale)
        history.append({"score": score, "approach": approach, "rationale": rationale})

        if score > best_score:
            best_score = score
            shutil.copy(SEGMENT_PATH, BEST_SEGMENT)
            print(f"  NEW BEST: {best_score:.4f}")
            with open(ROOT / "results" / "BEST_PARAMS.txt", "a") as f:
                f.write(f"\n[{datetime.now():%Y-%m-%d %H:%M:%S}] {best_score:.4f}\n{approach}\n{rationale}\n")
        else:
            if BEST_SEGMENT.exists():
                shutil.copy(BEST_SEGMENT, SEGMENT_PATH)

        if i % 5 == 0:
            update_memory(history)
            history = load_history()

        if best_score >= TARGET_SCORE:
            print(f"\nTarget reached: AP@0.5 = {best_score:.4f}")
            break

        time.sleep(2)
