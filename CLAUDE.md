# Sorkin / Bombelli Causal Set Lab — Claude Code Guide

## What this project is

Revival of the Bombelli 1987 PhD annealing program. The v1.0.0 revival is complete as a historical reconstruction. The active continuation is **SORKIN-2: algorithmic recoverability in the Bombelli annealer**.

SORKIN-2 is a separate diagnostic continuation. The working question is whether the historical annealer, with its historical energy, move set, schedule, and acceptance rule, can access known-truth causal realizations. It is not a physical embeddability program.

Core rules:

- Distinguish causal realization / embedding existence from annealer accessibility / algorithmic recoverability.
- Do not interpret annealer failure as non-embeddability.
- Do not interpret low final energy as manifoldlikeness.
- Avoid language like "discovers", "proves", "solves", "physical law", "Hauptvermutung test", "embeddability detector", or "manifoldlikeness detector".
- Prefer conservative terms: diagnostic, known-truth case, recoverability, accessibility, schedule sensitivity, basin accessibility, historical energy/move-set limitation.

The primary reference is `Pascal.pdf` (thesis code listing) and `Bombelli_1987_PhD.pdf`.

## Active science (read this before touching anything)

We are working through a numbered phase program. Each phase builds on the previous. Do not skip phases or retroactively change completed benchmark data.

| Phase | What it studies |
|---|---|
| 1 (a–e) | Atlas: success rate vs n, dim, structure |
| 2 (a–g) | Annealing internals: schedule, oracle, init basins, warmup modes |
| 3 (a–f) | PySR symbolic regression on annealing features; exploratory, not a claim |
| 4 (a–d) | Epsilon sweep, survival probe, seed robustness, robustness audit; diagnostic/exploratory |
| 5 | Seed curve morphology; exploratory |

Current frontier: **SORKIN-2 framing note exists** at `docs/SORKIN2_algorithmic_recoverability_note.md`. The next scientific step is a minimal matrix of known-truth cases. Phase 6 is not yet defined.

Phase2F diagnostic: guarded warmup preserves 18/18 small-noise rows on the tested grid.

## Module map

| File | Role |
|---|---|
| `cones.py` | Core: Pascal port of ConesSimulator, input parser, sprinkler |
| `causet_invariants.py` | Order-theoretic invariants (height, width, density, interval profile) |
| `cuda_backend.py` | CUDA shared lib loader (`build/libcones_cuda.so`) |
| `ensemble_scan.py` | Grid sweep over seeds, temps, cooling factors |
| `phase_diagram.py` | Phase map: success rate heatmaps over n × dim |
| `phase_refine.py` | Targeted rescan of specific cells |
| `dimension_sweep.py` | Sweep over embedding dimension |
| `schedule_sweep.py` | Sweep over annealing schedule params |
| `analyze_sweep.py` | Post-processing and aggregation |
| `validation_suite.py` | Cross-phase consistency checks |
| `tools/build_phase*.py` | One script per benchmark phase; generates CSV + MD |
| `legacy/` | Archived exploratory branches and local outputs; provenance only |

## Benchmark conventions

- All benchmark data lives in `benchmarks/foundation/`
- Every phase produces `phaseXY_name.csv` + `phaseXY_name.md`
- The `.md` is the human-readable verdict; the `.csv` is the data
- Never edit CSVs by hand; regenerate via `make regen-phaseXY`
- Seeds 1959, 1962, 1987 are canonical (historically meaningful years)
- Canonical test inputs: `benchmarks/tesis_like_6.in` (fast), `benchmarks/tesis_like_12.in` (medium)

## How to run things

```bash
make test                    # full test suite
make smoke                   # quick CPU sanity check (tesis_like_6.in, dim=2)
make smoke-cuda              # same with CUDA backend
make ensemble                # small ensemble scan to /tmp/
make phase                   # small phase diagram scan to /tmp/
make regen-phase2f           # regenerate a specific phase
```

Backend selection: `--backend auto` prefers CUDA if available, falls back to CPU.

Do not run experiments or regenerate data unless the user explicitly asks for that. For documentation-only work, do not run benchmark targets.

## Tools used daily

- **Python 3.12** — everything except Julia-specific work
- **Julia** — performance-critical runs when explicitly requested
- **PySR** (`pysr`) — existing Phase 3 exploratory artifacts only
- **KAN/GAM** — deferred until a clean multi-family known-truth dataset exists
- **CUDA** — GPU annealing via `libcones_cuda.so`
- **NumPy, SciPy, Matplotlib, Pandas** — standard stack
- **context7** MCP — use for live docs on any of these libraries

## What NOT to touch

- `Pascal.pdf`, `PostScript/`, `biblioteca/` — historical archive, do not edit
- `obj14.bin` — Sorkin's original binary, reference only
- Completed benchmark CSVs — regenerate, never manually edit
- `cones.py` core internals — changes invalidate the historical baseline; document why before touching
- `legacy/` — do not use as a source of claims without a new audit
- Phase3/4A/4B/5 — do not reactivate without an explicit decision

## Multi-agent patterns

When dispatching parallel agents for this project, brief them with:
- Which phase they're working on and what the previous phase found
- The target input file (`tesis_like_6.in` for fast, `tesis_like_12.in` for medium)
- Whether CUDA is available
- The output naming convention (`phaseXY_name.{csv,md}`)
- The SORKIN-2 framing: recoverability/accessibility diagnostics, not physical embeddability claims

Independent tasks suitable for parallel agents:
- Reading existing docs/results before any new script is proposed
- Auditing a single known-truth diagnostic question
- Checking whether an existing generator/test already answers the question

## Scientific interpretation rules

- A single failed run is not evidence of non-embeddability
- Success rate across repeated seeds at fixed (n, dim, schedule) is a recoverability diagnostic, not a physics verdict
- Final energy = 0 means the historical energy found correct causal relations under the current representation
- Final energy >0 means the algorithm did not access zero energy in that run; do not infer non-existence
- Schedule sensitivity is not a physical observable; it is an algorithmic diagnostic to map
- Phase2E/2F are warmup diagnostics
- Phase4C/4D are optimizer-seed and robustness diagnostics
- One question, one target, one file
- No KAN/PySR/GAM until a clean multi-family known-truth dataset exists
