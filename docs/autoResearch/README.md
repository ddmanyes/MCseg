# AI-Autonomous Segmentation Discovery (AutoResearch Template)

This directory contains the prompt template, agent runner, and starter script used to develop **MCseg v2** via AI-autonomous architecture search, as described in:

> *MCseg: High-Fidelity Visium HD Cell Segmentation via AI-Autonomous Pipeline Discovery*

## Overview

MCseg v2 was discovered over ~80 overnight cycles using a Claude Opus agent that iteratively proposed, implemented, and evaluated cell segmentation architectures against Xenium single-molecule ground truth (AP@0.5). The agent operated within a strictly sandboxed environment—it could only modify a single Python script—and received the current best implementation, full experiment history, and a distilled research memory at each iteration.

The key insight: by constraining the search to a single evaluable file with a fixed scoring function, the agent could explore bold architectural changes (multi-model ensembles, custom preprocessing, novel boundary strategies) without human guidance, converging on a configuration unlikely to emerge from manual intuition.

## Files

| File | Description |
|------|-------------|
| `program.md` | Human-readable task specification given to the agent at setup |
| `run_agent.py` | Agent runner — the main loop that calls Claude API and executes experiments |
| `segment_template.py` | Starter segmentation script (the sandbox file the agent modifies) |

## How to Adapt for Your Own Segmentation Problem

### 1. Define your evaluation function
The scoring function is the foundation. It must be:
- **Automated**: no human judgment required
- **Fast**: < 5 minutes per cycle
- **Objective**: a single numeric score (higher = better)

Examples: AP@0.5 vs ground-truth masks, F1 score, Dice coefficient.

### 2. Prepare your data
```
your_project/
├── run_agent.py          # copy and adapt from this template
├── prepare.py            # YOUR script: loads image, GT, and calls evaluate()
├── segment.py            # the sandbox file the agent will modify
├── memory.md             # agent's running research notes (auto-updated)
└── results/
    ├── experiment_log.jsonl
    ├── history_score.csv
    └── BEST_PARAMS.txt
```

### 3. Edit `run_agent.py`
Key parameters to set in `run_agent.py`:
```python
PYTHON_BIN = "/path/to/your/venv/bin/python"  # Python with your packages installed
# Edit SYSTEM_PROMPT to describe YOUR tissue type, image specs, and available tools
```

### 4. Edit `program.md`
Describe your task to the agent: what the image looks like, what counts as a cell, what the scoring metric is, and what tools are available.

### 5. Run
```bash
uv run python run_agent.py
```
Leave it running overnight. The agent will iterate automatically, keeping `segment_best.py` updated with the highest-scoring implementation found so far.

## Requirements

```bash
uv add anthropic cellpose scikit-image opencv-python scipy numpy pandas
```

Set your Anthropic API key:
```bash
export ANTHROPIC_API_KEY=your_key_here
```

## Citation

If you use this framework in your work, please cite:

```
[manuscript citation TBD upon acceptance]
```

And the original AutoResearch concept:
```
Karpathy A. AutoResearch: AI agents running research on single-GPU nanochat training automatically.
GitHub 2026. https://github.com/karpathy/autoresearch
```
