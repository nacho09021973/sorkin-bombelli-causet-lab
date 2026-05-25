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

## Help from bibliography

### Why this section exists

Before programming Kerr diagnostics, we need to separate what the causal-set
black-hole literature already provides from what SORKIN would still have to
invent. SORKIN-2 showed endpoint and selection failures in the historical
Bombelli annealer at N=36, while SORKIN-3 is exploring exact recovery in
Minkowski 1+1D through order-embedding structure. A future SORKIN-4/Kerr line
should therefore begin as a modest ordinal-diagnostic program: generate
known-truth causal sets in black-hole backgrounds, audit order-only
observables, and avoid claiming Kerr reconstruction.

Local bibliography artifacts are stored under
`docs/bibliography/sorkin4_black_holes/`. Dhital (2023) is included from its
public thesis metadata page because the advertised eGrove PDF URL returned HTTP
403 to command-line download.

### Papers and theses read

#### He & Rideout (2009) -- A Causal Set Black Hole

- Citation: Song He and David Rideout, "A Causal Set Black Hole", Classical
  and Quantum Gravity 26, 125015 (2009).
- Link: arXiv:0811.4235, https://arxiv.org/abs/0811.4235; DOI
  https://doi.org/10.1088/0264-9381/26/12/125015.
- Local copy: `docs/bibliography/sorkin4_black_holes/he_rideout_2009_causal_set_black_hole.pdf`.
- Causal-set-specific: yes. The paper is explicitly about computing causal
  relations needed for causal-set sprinklings in a black-hole spacetime.
- Spacetime: Schwarzschild black hole, using Eddington-Finkelstein
  coordinates; not Kerr.
- Numerical sprinkling: yes. The paper describes sprinkling events into a
  bounded Schwarzschild region and then deciding pairwise causal relations.
- Causal-relation algorithm: yes. It gives a practical algorithm for deciding
  whether two Schwarzschild events are causally related by first checking
  sufficient spacelike/timelike bounds and then numerically integrating null
  geodesics for generic pairs.
- What it does: supplies the core Schwarzschild causal-relation machinery that
  flat Minkowski sprinklers lack. It also documents uniform sprinkling in
  Schwarzschild volume measure and implementation details for sorting events
  and testing pairwise relations.
- What helps SORKIN: it is the natural first benchmark before Kerr. A
  SORKIN-4 program should first reproduce a small Schwarzschild causal matrix
  from known coordinates and validate ordinal observables against that matrix.
- What it does not solve: it does not give an order-only estimator of black
  hole mass, horizon position, or spin. It also does not solve Kerr causal
  relations.
- Kerr relevance: limited but important. The authors explicitly say the
  null-geodesic approach should in principle generalize, and mention Kerr as a
  plausible but more complicated target because spherical symmetry is broken
  and closed-timelike-curve regions must be avoided.

#### Dou & Sorkin (2003) -- Black Hole Entropy as Causal Links

- Citation: Djamel Dou and Rafael D. Sorkin, "Black Hole Entropy as Causal
  Links", Foundations of Physics 33, 279--296 (2003).
- Link: arXiv:gr-qc/0302009, https://arxiv.org/abs/gr-qc/0302009; DOI
  https://doi.org/10.1023/A:1023781022519.
- Local copy: `docs/bibliography/sorkin4_black_holes/dou_sorkin_2003_black_hole_entropy_causal_links.pdf`.
- Causal-set-specific: yes. The paper studies black-hole entropy by counting
  irreducible causal relations, i.e. links, in a causal set.
- Spacetime: black-hole horizons with spherical examples, including
  equilibrium Schwarzschild and collapse/far-from-equilibrium cases. The exact
  computations are mainly dimensionally reduced 1+1D/null-surface cases, with
  arguments for four-dimensional area scaling.
- Numerical sprinkling: not as a modern simulation framework. The paper is
  primarily analytic/kinematic and computes expected link counts through
  integrals.
- Causal-relation algorithm: no general pairwise Schwarzschild/Kerr algorithm.
  Its central object is the link condition: x precedes y and the Alexandrov
  interval between them is empty except for x and y.
- What it does: defines horizon-crossing links as candidate "horizon
  molecules" and argues that, with max/min conditions to prevent double
  counting near a hypersurface Sigma, the expected number of such links scales
  with horizon area.
- What helps SORKIN: gives order-theoretic observables that do not require
  metric reconstruction once a causet is available: links, horizon-crossing
  links in oracle mode, maximal/minimal conditions, and area-scaling checks.
- What it does not solve: it does not provide a deployable way to find a
  horizon from the order alone, does not compute Kerr, and does not derive a
  spin estimator.
- Kerr relevance: indirect. The paper names rotating/deformed black holes as
  an obvious direction for future work, but does not carry it out.

#### Homšak & Veroni (2024) -- Boltzmannian state counting for black hole entropy in Causal Set Theory

- Citation: Vid Homšak and Stefano Veroni, "Boltzmannian state counting for
  black hole entropy in Causal Set Theory", Physical Review D 110, 026015
  (2024).
- Link: arXiv:2404.11670, https://arxiv.org/abs/2404.11670; DOI
  https://doi.org/10.1103/PhysRevD.110.026015; code link cited by the paper:
  https://github.com/vidh2000/MSci_Schwarzschild_Causets.
- Local copy: `docs/bibliography/sorkin4_black_holes/homsak_veroni_2024_boltzmannian_state_counting.pdf`.
- Causal-set-specific: yes. It is a numerical causal-set study of
  black-hole thermodynamics.
- Spacetime: Schwarzschild black hole in Eddington-Finkelstein original
  coordinates; not Kerr.
- Numerical sprinkling: yes. The paper reports a highly parallelized C++
  framework for Schwarzschild causal sets, reaching over one million points in
  non-conformally flat spacetime.
- Causal-relation algorithm: yes for Schwarzschild. It builds on He & Rideout
  and discusses a corrected treatment for causal relations inside the horizon,
  then uses null-geodesic arrival conditions and cheap sufficient tests before
  expensive generic pair handling.
- What it does: modernizes Schwarzschild causal-set simulation, counts horizon
  molecules, validates Poisson sprinkling in equal-volume regions, visually and
  analytically checks causal-relation correctness in lower-dimensional limits,
  and tests area scaling of molecule counts.
- What helps SORKIN: provides the best practical benchmark path before Kerr:
  reproduce small Schwarzschild sprinklings, causal matrices, link/molecule
  counts, and area-scaling diagnostics before asking whether any ordinal
  observable can distinguish rotation.
- What it does not solve: it does not infer spin a/M, does not give a Kerr
  sprinkling pipeline, and does not provide an order-only horizon finder or
  metric reconstruction algorithm.
- Kerr relevance: mostly negative evidence for scope. The modern numerical
  literature found here still focuses on Schwarzschild and horizon molecules,
  so Kerr should be treated as a later exploratory extension rather than a
  solved baseline.

#### Eichhorn, Gamito & Stokes (2026) -- Towards black-hole horizons and geodesic focusing in causal sets

- Citation: Astrid Eichhorn, Pedro Gamito, and Nawder Stokes, "Towards
  black-hole horizons and geodesic focusing in causal sets", arXiv:2605.06813
  (2026).
- Link: arXiv:2605.06813, https://arxiv.org/abs/2605.06813.
- Local copy: `docs/bibliography/sorkin4_black_holes/eichhorn_gamito_stokes_2026_horizons_geodesic_focusing.pdf`.
- Causal-set-specific: yes. The paper asks how a black-hole horizon can be
  identified in a causal set.
- Spacetime: toy-model 1+1D Schwarzschild, with discussion of regular
  black-hole examples such as Hayward; not Kerr.
- Numerical sprinkling: yes. It uses finite sprinklings into 1+1D black-hole
  patches, including examples with large point counts, and studies discrete
  horizon diagnostics.
- Causal-relation algorithm: yes for the 1+1D toy setting. It uses analytic
  null geodesics to infer causal relations and then transitive reduction to
  obtain links.
- What it does: gives modern order-level horizon diagnostics beyond molecule
  counting: longest-chain behavior from minimal antichains, ladder tracers of
  null geodesics, a discrete expansion proxy whose sign changes across the
  horizon, and fuzzy ladders for longer null-geodesic tracking.
- What helps SORKIN: this is the closest match to "ordinal diagnostics". It
  suggests concrete order observables to test after Schwarzschild causal-matrix
  generation: minimal-antichain partitions, longest-chain bimodality,
  ladder/fuzzy-ladder structures, and expansion-sign checks.
- What it does not solve: it is explicitly a proof-of-principle, mostly in
  1+1D. It does not provide a Kerr generator, a spin estimator, or a validated
  order-only reconstruction method for 3+1D rotating black holes.
- Kerr relevance: conceptual, not direct. It strengthens the case that
  SORKIN-4 should start with horizon diagnostics in Schwarzschild-like
  settings before trying any Kerr spin observable.

#### Marr (2007) -- Black hole entropy from causal sets

- Citation: Sarah Kathryn Marr, "Black hole entropy from causal sets", PhD
  thesis, Imperial College London (2007).
- Link: https://hdl.handle.net/10044/1/11818.
- Local copy: `docs/bibliography/sorkin4_black_holes/marr_2007_black_hole_entropy_from_causal_sets.pdf`.
- Causal-set-specific: yes. It is a thesis on black-hole entropy from causal
  set structures.
- Spacetime: black-hole entropy and horizon-molecule setting, historically
  downstream of Dou & Sorkin. Secondary reviews cite it for higher-dimensional
  entropy formulas using sub-causal sets rather than only links.
- Numerical sprinkling: historical/deep-reference role; use it for conceptual
  and implementation background rather than as the first executable benchmark.
- Causal-relation algorithm: not the main pairwise Schwarzschild/Kerr
  algorithm for SORKIN. Its relevance is horizon-molecule/statistical
  structure, not Kerr causal-order generation.
- What it does: extends the black-hole entropy from causal sets line and is a
  historical reference for attempts to move beyond simple 1+1D link counting.
- What helps SORKIN: provides background on what "molecule" observables might
  need to become in higher dimensions, which is useful before designing any
  SORKIN-4 invariant list.
- What it does not solve: it does not make Kerr spin recoverable from order
  data and should not be treated as a ready benchmark pipeline.
- Kerr relevance: indirect. It supports a broader horizon-observable design
  space, not a Kerr-specific algorithm.

#### Dhital (2023) -- Black Hole Entropy in the Causal Set Approach

- Citation: Ayush Dhital, "Black Hole Entropy in the Causal Set Approach",
  M.S. thesis, University of Mississippi (2023).
- Link: https://egrove.olemiss.edu/etd/2674/.
- Local copy: no PDF stored. The eGrove metadata page advertises a PDF, but
  command-line download returned HTTP 403; this summary is based on the public
  metadata page and abstract.
- Causal-set-specific: yes. The thesis studies black-hole entropy in the
  causal-set approach.
- Spacetime: Schwarzschild black hole and generalized d-dimensional flat
  black-hole setups according to the public abstract; not Kerr.
- Numerical sprinkling: yes, according to the abstract. It reports computer
  simulations of horizon-crossing links for Schwarzschild and MATLAB numerical
  integration for generalized flat black-hole cases.
- Causal-relation algorithm: not confirmed from the PDF because the PDF was not
  downloaded. Treat as implementation-detail background only until the PDF is
  inspected.
- What it does: appears to provide practical numerical details for
  horizon-crossing links, max/min conditions, and area-scaling checks.
- What helps SORKIN: useful as a secondary implementation reference for
  reproducing Schwarzschild link-count experiments before Kerr.
- What it does not solve: based on the public abstract, it does not address
  Kerr spin recovery or order-only metric reconstruction.
- Kerr relevance: indirect and currently low-confidence until the full thesis
  PDF is obtained.

### Current literature status

- Causal-set black-hole work exists, especially for Schwarzschild causal
  relations, horizon-crossing links, and horizon molecules.
- The strongest concrete algorithmic material found here is Schwarzschild:
  He & Rideout give the pairwise causal-relation algorithm, and Homšak &
  Veroni scale it into a modern numerical framework.
- The newest horizon-diagnostic material found here is Eichhorn, Gamito &
  Stokes (2026), which moves from entropy/link counting toward direct
  order-level horizon diagnostics through chains, ladders, fuzzy ladders, and
  geodesic-focusing proxies.
- I found no mature standard algorithm for "Kerr sprinkling -> estimate spin
  a/M from the order alone" in this minimal search.
- Therefore SORKIN-4/Kerr should be formulated as an exploratory
  ordinal-diagnostic program, not as Kerr reconstruction.

### Practical help for SORKIN-4

- Reproduce Schwarzschild first, using He & Rideout's causal-relation logic as
  the known-truth generator benchmark.
- Use horizon-crossing links and horizon molecules from Dou & Sorkin and
  Homšak & Veroni as order-level observables in oracle mode, where coordinates
  identify inside/outside and the hypersurface only during validation.
- Add area-scaling checks before any Kerr attempt: molecule/link counts should
  scale with horizon area in the Schwarzschild reproduction step.
- Add horizon-localization checks before any Kerr attempt: longest-chain
  bimodality, ladder/fuzzy-ladder availability, and expansion-sign diagnostics
  from Eichhorn, Gamito & Stokes should be tested on small known-truth
  Schwarzschild cases.
- Treat Homšak & Veroni's large-scale C++ framework and validation strategy as
  the performance/verification reference, not as an immediate requirement for
  the Python SORKIN repo.
- Use Marr and Dhital as secondary references for horizon-molecule variants and
  numerical link-counting details, not as evidence that a Kerr spin estimator
  already exists.
- Keep a strict separation between known-truth generation coordinates and the
  final ordinal estimator input.
- Only after Schwarzschild reproduction should Kerr be attempted, with spin
  a/M framed as an exploratory order-only diagnostic target.

### Guardrails

- Do not claim Kerr causal sets are unexplored in an absolute sense unless a
  proper literature review confirms it.
- Do not claim metric reconstruction.
- Do not use coordinates as input for the final ordinal estimator, except in
  oracle/diagnostic mode.
- Separate known-truth generation coordinates from algorithm input.
- Do not present Schwarzschild horizon-molecule success as evidence that Kerr
  spin is recoverable from the order alone.

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
