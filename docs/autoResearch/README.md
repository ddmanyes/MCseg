# Agent-guided workflow search templates

These files illustrate the development loop behind MCseg. Researchers define candidate operations, reference data, scoring, prompts, and execution constraints; an agent proposes and evaluates workflow variants, followed by researcher review. Routine MCseg use runs the retained cpsam workflow locally and does not invoke this agent or require an external language-model API.

The templates here are an adaptation starting point, not a complete archive of every historical search run. The configured API model in the example runner should not be treated as evidence of the model used in every manuscript experiment.

| File | Purpose |
| --- | --- |
| `program.md` | Example task specification |
| `run_agent.py` | Development loop that calls an external API and executes candidates |
| `segment_template.py` | Starter workflow to adapt to your reference data |

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

Chan, C.-R., Chang, N.-W., Wang, C.-Y., Tan, H.-Y., and Lin, S.-J. (2026). **MCseg: AI agent-guided workflow search for no-code cell segmentation and transcript attribution in spatial transcriptomics.** Manuscript.

The planned bioRxiv reference will be added after posting and DOI assignment. See also the original [AutoResearch project](https://github.com/karpathy/autoresearch).
