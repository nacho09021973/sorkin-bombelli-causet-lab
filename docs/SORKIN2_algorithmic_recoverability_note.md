# SORKIN-2 — Algorithmic recoverability in the Bombelli annealer

## Purpose

SORKIN-2 studies algorithmic recoverability in the historical Bombelli
annealer. It does not treat the annealer as a definitive test of whether a
causal set can be realized or embedded.

The central distinction is:

- causal realization or embedding existence;
- annealer accessibility or algorithmic recoverability.

The first is an existence question. The second is a diagnostic question about a
particular constructive procedure.

## Main Claim

The Bombelli annealer should not be interpreted as a definitive embeddability
test. It is a constructive algorithm with a historical energy, move set,
schedule, and acceptance rule.

Failure to converge to zero energy is evidence about this algorithmic setup. It
is not direct evidence that no causal realization exists.

## Thesis-Based Motivation

The 1987 thesis motivates this narrower reading.

- Bombelli defines energy zero as the state in which the induced causal
  relations are correct.
- The energy measures causal badness of the current configuration.
- The annealer may fail because of local minima, temperature schedule,
  initialization, or move set.
- The thesis notes that the Euclidean displacement and move representation is
  conceptually less satisfying than a formulation that uses the Lorentzian
  structure more directly.

This suggests that the interesting modern question is the gap between
zero-energy existence and annealer accessibility.

## Allowed Interpretation

The following interpretations are in scope:

- optimizer recoverability;
- basin accessibility;
- schedule sensitivity;
- historical energy and move-set limitations;
- diagnostics on known-truth cases.

## Forbidden Interpretation

The following interpretations are out of scope:

- proof of non-embeddability;
- manifoldlikeness classification;
- tests of continuum uniqueness conjectures;
- physical causal-set classification;
- KAN, PySR, or GAM target construction at this stage.

## Why No Modeling Yet

KAN, PySR, and GAM are premature because the current tabular data are derived
probes, mostly one-family, and do not yet form a clean multi-family
known-truth dataset.

Before modeling, SORKIN-2 needs a documented matrix of known-truth cases with
stable provenance and a single diagnostic target.

## Next Step

The next step is to design a minimal matrix of known-truth cases:

- crown or another simple known success case;
- simple Minkowski-like cases with known realization;
- deliberately harder cases;
- one target only: reaches zero energy yes/no, or verified exact relation
  recovery yes/no;
- no machine learning until the matrix exists and has provenance.

## Relation to the 1987 Revival

The v1.0.0 revival remains complete as a historical reconstruction.

SORKIN-2 is a separate interpretive and diagnostic continuation. It does not
replace the revival. It uses the revival to ask a narrower algorithmic question:
when a zero-energy realization is known or expected, when does the historical
annealer recover it, and when does it fail to access it?
