# Thirty-Seven Years of Simulated Annealing on a Causal Set
## A Revival of Bombelli (1987) with 2026 Tools

*Ignacio Arancibia · Claude Sonnet 4.6 (Anthropic) · 2026*

---

> *"An application of simulated annealing"*
> — Title of Appendix A.2, Luca Bombelli's PhD thesis, 1987

---

Rafael Sorkin and Luca Bombelli introduced causal set theory in 1987
with the bold conjecture that spacetime, at the Planck scale, is a
locally finite partially ordered set. The same year, Bombelli appended
to his thesis a Pascal program that tried to embed small causal sets
into Minkowski spacetime by simulated annealing. The program ran on
the workstations of the era, produced results for a handful of cases,
and was never published as a standalone tool.

This note documents what happens when that program is brought back to
life and subjected to the computational tools of 2026: a GPU
accelerator, ensemble statistics, symbolic regression, and an AI
collaborator. We make no claim of advancing causal set theory. We do
claim to have learned something about the algorithm itself — and about
how much was invisible to its author through no fault of his own.

---

## I. The Instruments

| Aspect | Bombelli (1987) | This work (2024–2026) |
|:---|:---|:---|
| Language | Pascal | Python 3.12 (faithful port) |
| Hardware | ~1 MIPS workstation | CPU + NVIDIA GPU (CUDA) |
| Runs per case | Single run | Ensemble of K seeds (K ≥ 8) |
| Schedule selection | Fixed by hand | Empirical grid scan |
| Warmup protocol | 10 unconditional accepts | Three modes compared: legacy / skip / guarded |
| Non-manifoldlike controls | None | Kleitman–Rothschild orders, suspended corona posets |
| Dimension estimators | None (annealing only) | Myrheim–Meyer, Meyer midpoint scaling |
| Structural invariants | None | 15 order-theoretic features (chain counts, link density, height, …) |
| Correlation analysis | None | Spearman ρ, partial Spearman controlling for *n* |
| Machine learning | None | PySR symbolic regression, Kolmogorov–Arnold networks |
| AI collaborator | None | Claude Sonnet 4.6 |
| Reproducibility | Thesis appendix | Frozen benchmark CSV + `make regen-phaseXY` |

---

## II. Computational Reach

The original thesis demonstrated the method on very small causets where
a single run could be completed in minutes on the hardware of the day.
The table below compares the practical frontier of the algorithm across
the two eras. All 2026 entries use the same energy function and move set
as the original.

| Causal set size *n* | 1987 frontier | 2026 ensemble (8 seeds, gpu) | Note |
|---:|:---|:---|:---|
| 6 | Single run, qualitative | Success rate 25–50 %, phase map available | Best cells: *dim* = 3–4 |
| 12 | Single run, qualitative | Success rate 0–12.5 %, schedule matters strongly | Benchmark `tesis_like_12.in` |
| 16 | Not reported | Success rate 0–12.5 %, frontier of useful search | Budget-dependent |
| 24 | Not accessible | All runs time out at default budget | Beyond brute-force reach |
| 32–64 | Not accessible | Ensemble statistics available; floor pathology characterised | Phases 4A–4D |

*"Not accessible"* means the combination of hardware speed and single-run
methodology made systematic exploration impossible, not that Bombelli
lacked the insight to attempt it.

---

## III. The Same Benchmark, Thirty-Seven Years Apart

The reproducible input `tesis_like_12.in` — twelve elements, a
tesis-like causal structure — serves as the common witness.

| Metric | Bombelli defaults (T₀ = 100, α = 0.9) | Tuned schedule (T₀ = 180, α = 0.8) | Change |
|:---|---:|---:|---:|
| Mean final energy (100 seeds) | 20.021 | 0.166 | −99.2 % |
| Zero-energy runs | 0 / 100 | 0 / 100 | — |
| Schedule selected by | Intuition / thesis norm | Empirical grid scan | — |

The tuned schedule uses the same annealer, the same energy function,
and the same move set as the original. The only difference is two
numbers. A 99 % reduction in mean final energy from changing two
parameters is not a flaw in the 1987 work — it is a property of the
algorithm that was simply not visible without ensemble runs across a
parameter grid.

---

## IV. What the Warmup Was Doing

Phase 2D–2F of this study isolated the dominant failure mode of the
historical algorithm. The original warmup phase makes 10 unconditional
accept steps before annealing begins — a design choice intended to
equilibrate the system at high temperature. The table below shows what
those 10 steps do to controlled initializations.

| Initialization | Legacy warmup | Skip warmup | Guarded warmup |
|:---|:---|:---|:---|
| Ground truth (E = 0) | Preserved 18/18, E = 0.000 | Preserved 18/18, E = 0.000 | Preserved 18/18, E = 0.000 |
| Truth + small noise (ε = 10⁻³) | **Destroyed**: mean E = 18.92, 16/18 | Mean E = 12.12, 17/18 | **Recovered**: mean E = 0.001, **18/18** |
| Truth + medium noise (ε = 5×10⁻²) | Mean E = 395.8, 0/18 | Mean E = 286.0, 0/18 | Mean E = 255.3, 0/18 |
| Random initialization | Mean E = 405.2, 8/18 | Mean E = 307.9, 11/18 | Mean E = 271.4, **12/18** |

*Grid: d ∈ {2, 3, 4}, n ∈ {32, 64}, seeds 1959/1962/1987 — 18 cells per row.*

The guarded warmup accepts a proposed move only if it does not increase
the energy (greedy descent). It is an external wrapper around the same
`ConesSimulator` internals — no change to the energy, the move set, or
the cooling schedule. The small-noise basin problem, invisible to single
runs, is essentially solved by this one modification.

**The oracle check** (Phase 2C) confirmed that the energy formula
returns exactly 0.0 at ground-truth coordinates in 18/18 cases. The
Bombelli energy is correct; the failure was in the warmup dynamics, not
in the physics.

---

## V. What the Order-Theoretic Invariants Reveal

The 1987 program had no way to ask "is this causal set manifoldlike
before we try to embed it?" Phase 1 of this study adds that question.

The table below compares two independent dimension estimators — the
Myrheim–Meyer formula and Meyer's midpoint scaling — on Minkowski
sprinklings (manifoldlike by construction) and two non-manifoldlike
controls (Kleitman–Rothschild three-layer orders and suspended corona
posets), at n = 256, ensemble of 5 seeds.

| Family | d_spacetime | n | MM dim | midpoint dim | \|discrepancy\| |
|:---|:---:|---:|---:|---:|---:|
| Minkowski | 2 | 256 | **2.02** | **2.06** | 0.07 |
| Minkowski | 3 | 256 | **3.06** | **3.00** | 0.28 |
| Minkowski | 4 | 256 | **4.07** | **3.80** | 0.42 |
| Kleitman–Rothschild | — | 256 | 2.37 | 4.71 | **2.34** |
| Corona poset | — | 256 | 1.98 | 7.00 | **5.02** |

For manifoldlike sprinklings the two estimators converge to the true
dimension as *n* grows and their discrepancy shrinks. For
non-manifoldlike controls the estimators diverge: their discrepancy
*grows* with *n*. The sign of d\|discrepancy\|/dn is opposite for the two
families across the entire grid n ∈ {32, 64, 128, 256}. This
finite-size trajectory is itself the diagnostic — a tool Bombelli could
not have used in 1987 because it requires running the program many
times at many sizes.

---

## VI. What the Correlate Audit Found

Phase 4D (the robustness-vs-invariants audit) asked whether any
order-theoretic property of a causal set predicts how unstable the
embedding optimizer is across different random seeds. At per-seed level
(N = 90 observations, n ∈ {32, 48, 64}):

| Invariant | Target robustness metric | Raw Spearman ρ | Partial ρ (controlling for *n*) |
|:---|:---|---:|---:|
| `relation_count` | Floor saturation fraction | +0.941 | +0.903 |
| `chain2_count` | Floor saturation fraction | +0.941 | +0.903 |
| `height` | Floor saturation fraction | +0.836 | +0.779 |
| `abs_discrepancy_mm_midpoint` | Floor saturation fraction | −0.759 | −0.659 |
| `mm_dim` | Floor saturation fraction | −0.698 | −0.530 |

**Verdict: ORDER_THEORETIC_CORRELATE_DETECTED.** The correlation
survives controlling for finite-size scaling: even within a fixed *n*
stratum, denser, taller causal sets are harder for the optimizer. This
is a statement about the algorithm's difficulty landscape, not about
physical embeddability. But it is a statement that could not have been
made in 1987.

---

## VII. A Timeline

| Year | Event |
|---:|:---|
| 1987 | Sorkin and Bombelli introduce causal set theory (*Phys. Rev. Lett.* 59, 521) |
| 1987 | Bombelli writes the Pascal annealing program (PhD thesis, Appendix A.2) |
| 1987–2024 | Program dormant. CST develops theoretically. Modern tools emerge. |
| 2024 | Faithful port to Python 3.12, validated against thesis inputs |
| 2024 | CUDA backend added; GPU acceleration confirmed |
| 2024–2025 | Phases 1–3: structural atlas, embedding bridge, schedule probe, warmup audit, PySR |
| 2025–2026 | Phases 4–5: epsilon sweep, survival probe, seed robustness, morphology audit |
| 2026 | Phase 4D n-control extension: correlates survive finite-size correction |
| 2026 | This note |

---

## VIII. What Remains the Same

The energy function is Bombelli's. The move set is Bombelli's. The
cooling rule (`4 × exp(−ΔE / T)`) is Bombelli's. The input format is
Bombelli's. The core loop is Bombelli's.

Everything else — the number of runs, the parameter scan, the
controls, the invariants, the ML analysis, the n-control audit — is
what 37 years of Moore's law, ensemble statistics, and AI assistance
made visible. The pioneers drew the map. We learned how large it is.

---

## Acknowledgements

This work is a revival and empirical characterisation of the
computational program introduced in:

> L. Bombelli, *Space-time as a Causal Net*, PhD thesis, Syracuse
> University, 1987.

and motivated by the foundational paper:

> L. Bombelli, J. Lee, D. Meyer, R. D. Sorkin, *Space-time as a causal
> set*, Phys. Rev. Lett. **59**, 521 (1987).

The GPU backend uses NVIDIA CUDA. The symbolic regression uses PySR
(Cranmer 2023). The AI collaborator is Claude Sonnet 4.6 (Anthropic).
All benchmark data and the Python port are reproducible via
`make regen-phaseXY` from the project repository.

---

*The program works. It always did.*
