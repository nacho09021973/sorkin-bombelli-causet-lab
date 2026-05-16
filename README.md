# Sorkin–Bombelli Causal Set Lab

Order-theoretic diagnostics and reproducible benchmarks inspired by the Bombelli–Sorkin causal set embedding program.

This repository is an independent exploratory project on causal sets. It started as a reconstruction of the historical idea of embedding finite causal sets into Minkowski-like coordinates using energy minimization / annealing, and has evolved into a more conservative diagnostic framework:

> Before asking whether an algorithm can embed a causal set, ask whether the partial order itself carries structural evidence of manifoldlike behavior.

## Motivation

Causal set theory suggests that spacetime geometry may be encoded in a locally finite partial order. A central question is therefore not merely whether a numerical optimizer can find coordinates for a given causet, but whether the causet itself contains order-theoretic signatures compatible with a manifoldlike origin.

This project explores that question through small, reproducible benchmarks.

The guiding distinction is:

1. **Order structure**: does the causet have internal invariants compatible with a manifoldlike sprinkling?
2. **Embeddability**: does a low-energy embedding exist?
3. **Algorithmic recovery**: can the historical optimizer actually find it?

Earlier embedding experiments can mix these questions together. This repository tries to separate them.

## Current status

The project currently contains a reproducible diagnostic foundation and a growing structural atlas.

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

The current focus is **not** on optimizing embeddings. The current focus is on identifying which order-theoretic observables carry geometric signal before any embedding algorithm is run.

## Scientific conclusion so far

The current diagnostics separate families of causal sets better than individual causal sets.

For Minkowski sprinklings, the dimension estimators show convergence behavior with increasing `n`.

For non-manifoldlike controls such as Kleitman–Rothschild and suspended crown / corona posets, the diagnostics show qualitatively different scaling behavior. In particular, Myrheim–Meyer estimates may remain relatively flat while midpoint scaling can drift strongly with scale.

The main current conclusion is conservative:

> The present diagnostics provide ensemble-level separation between manifoldlike sprinklings and non-manifoldlike controls, but they do not yet provide a robust per-causal-set classifier.

This repository does **not** claim to solve the causal set Hauptvermutung. It provides a reproducible laboratory for testing which order-theoretic invariants carry useful geometric information.

## Project philosophy

- Reproducibility over impressive claims.
- Order-theoretic diagnostics before embedding.
- Negative controls before interpretation.
- Conservative conclusions.
- No temporary scripts without provenance.
- Tests and benchmark outputs should be versioned.
- A good negative result is still a result.

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
├── research_agenda_2026.md
├── results_note_2026.md
├── Makefile
└── README.md
