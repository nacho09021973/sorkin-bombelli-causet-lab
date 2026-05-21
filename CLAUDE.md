# Sorkin / Bombelli Causal Set Lab — Claude Code Guide

## What this project is

Revival of the Bombelli 1987 PhD annealing program. The scientific question: which finite causal sets are recoverable by the current optimizer pipeline into low-dimensional Minkowski spacetime, and what controls the transition between easy, hard, and non-embeddable cases. The port is a faithful historical baseline; extensions probe the landscape statistically.

The primary reference is `Pascal.pdf` (thesis code listing) and `Bombelli_1987_PhD.pdf`.

## Active science (read this before touching anything)

We are working through a numbered phase program. Each phase builds on the previous. Do not skip phases or retroactively change completed benchmark data.

| Phase | What it studies |
|---|---|
| 1 (a–e) | Atlas: success rate vs n, dim, structure |
| 2 (a–g) | Annealing internals: schedule, oracle, init basins, warmup modes |
| 3 (a–f) | PySR symbolic regression on annealing features |
| 4 (a–d) | Epsilon sweep, survival probe, seed robustness, robustness audit |
| 5 | Seed curve morphology |

Current frontier: **Phase 5** complete. Phase 6 not yet defined.

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

## Tools used daily

- **Python 3.12** — everything except Julia-specific work
- **Julia** — performance-critical runs, KAN experiments
- **PySR** (`pysr`) — symbolic regression on annealing features (Phase 3)
- **KAN** (Kolmogorov-Arnold Networks) — alternative function approximation
- **CUDA** — GPU annealing via `libcones_cuda.so`
- **NumPy, SciPy, Matplotlib, Pandas** — standard stack
- **context7** MCP — use for live docs on any of these libraries

## What NOT to touch

- `Pascal.pdf`, `PostScript/`, `biblioteca/` — historical archive, do not edit
- `obj14.bin` — Sorkin's original binary, reference only
- Completed benchmark CSVs — regenerate, never manually edit
- `cones.py` core internals — changes invalidate the historical baseline; document why before touching

## Multi-agent patterns

When dispatching parallel agents for this project, brief them with:
- Which phase they're working on and what the previous phase found
- The target input file (`tesis_like_6.in` for fast, `tesis_like_12.in` for medium)
- Whether CUDA is available
- The output naming convention (`phaseXY_name.{csv,md}`)

Independent tasks suitable for parallel agents:
- Regenerating different benchmark phases simultaneously
- Running PySR fits on different feature subsets
- Scanning different regions of the (n, dim) phase diagram

## Scientific interpretation rules

- A single failed run is not evidence of non-embeddability
- Success rate across ≥8 seeds at fixed (n, dim, schedule) is the meaningful observable
- Final energy = 0 means exact embedding found; >0 means optimizer stopped at a local minimum
- Schedule sensitivity is not a physical observable — it is a numerical artifact to be mapped, not minimized blindly
