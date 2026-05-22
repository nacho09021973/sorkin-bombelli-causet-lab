# Sorkin-Bombelli Causal Set Lab

This repository contains a reproducible revival of Bombelli's 1987 causal-set
annealer. The v1.0.0 revival is complete as a historical reconstruction.

The current continuation is **SORKIN-2 — algorithmic recoverability in the
Bombelli annealer**. SORKIN-2 studies when the historical annealer recovers
known-truth causal realizations and when it fails to access them.

SORKIN-2 is not a proof of non-embeddability, not a manifoldlikeness
classifier, not a Hauptvermutung test, and not currently an ML project.

The framing note is
[`docs/SORKIN2_algorithmic_recoverability_note.md`](docs/SORKIN2_algorithmic_recoverability_note.md).

## Motivation

The original Bombelli program asks a constructive question: given a causal set
and a target dimension, can a particular simulated annealing procedure find
coordinates whose induced causal relations match the input order?

The guiding distinction is:

1. **Causal realization / embedding existence**: does a zero-energy realization
   exist?
2. **Annealer accessibility / algorithmic recoverability**: can the historical
   energy, move set, schedule, and acceptance rule find it?

Failure to reach zero energy is therefore a diagnostic of this algorithmic
setup. It is not direct evidence that no realization exists.

## Current status

- **Historical revival:** complete and archived in v1.0.0.
- **Current continuation:** SORKIN-2 diagnostic framing is fixed.
- **Next scientific step:** design a minimal matrix of known-truth cases.
- **Provenance freeze:** warmup and robustness diagnostics are committed locally.
- **Exploratory artifacts:** Phase3/4A/4B/5 remain non-claim exploratory
  material unless later formalized. Non-active exploratory branches are archived
  under `legacy/`.
- **Foundation README reconciliation:** may still be pending if
  `benchmarks/foundation/README.md` has unresolved diff.

Implemented so far:

- canonical Minkowski sprinklings;
- Myrheim–Meyer dimension estimator;
- midpoint scaling estimator;
- three-chain abundance observable;
- chain count diagnostics;
- Kleitman–Rothschild controls;
- suspended crown / corona controls;
- finite-size scaling atlases;
- reproducible CSV / Markdown benchmark outputs;
- regression and integrity tests.

The current focus is **not** on optimizing embeddings and **not** on ML. The
current focus is defining known-truth recoverability diagnostics for the
historical annealer.

## Scientific conclusion so far

The current diagnostics separate families of causal sets better than individual
causal sets, but this remains background context for SORKIN-2 rather than a
claim that the annealer classifies physics.

For Minkowski sprinklings, the dimension estimators show convergence behavior with increasing `n`.

For non-manifoldlike controls such as Kleitman–Rothschild and suspended crown / corona posets, the diagnostics show qualitatively different scaling behavior. In particular, Myrheim–Meyer estimates may remain relatively flat while midpoint scaling can drift strongly with scale.

The main current conclusion is conservative:

> The present diagnostics provide ensemble-level separation between manifoldlike sprinklings and non-manifoldlike controls, but they do not yet provide a robust per-causal-set classifier.

The next SORKIN-2 step is narrower: build a minimal known-truth matrix and ask
whether the historical annealer reaches zero energy or exact relation recovery.

## Project philosophy

- Reproducibility over impressive claims.
- Known-truth diagnostics before modeling.
- Negative controls before interpretation.
- Conservative conclusions.
- No temporary scripts without provenance.
- Tests and benchmark outputs should be versioned.
- A good negative result is still a result.
- No KAN/PySR/GAM until a clean multi-family known-truth dataset exists.
- Do not use `legacy/` as a source of claims without a new audit.
- Do not reactivate Phase3/4A/4B/5 without an explicit decision.

## Repository structure

```text
.
├── causet_invariants.py
├── validation_suite.py
├── cones.py
├── tools/
│   ├── build_foundation_benchmarks.py
│   ├── build_phase1c_scaling_atlas.py
│   └── build_phase1d_structural_atlas.py
├── tests/
├── benchmarks/
│   └── foundation/
├── legacy/
├── research_agenda_2026.md
├── results_note_2026.md
├── Makefile
└── README.md
