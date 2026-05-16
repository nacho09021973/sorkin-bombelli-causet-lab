# Phase 2C Oracle Embedding Audit

Ground-truth consistency check for the Minkowski cases in
Phase 2/2B. No annealing is performed. The ground-truth
coordinates from the canonical sprinkler are evaluated
directly against the Bombelli energy, the causal-matrix
convention, and the Lorentz-invariant interval metric.

## Verdict

**ORACLE PASSES**

The ground-truth coordinates are a zero-energy configuration. The causal matrices reconstructed from those coordinates match the stored matrices exactly. The Lorentz-invariant interval residual of the truth against itself is zero. These three checks confirm that the Bombelli energy, the causal-matrix convention, and the interval metric are mutually consistent at the ground-truth embedding. The failure in Phase 2/2B is therefore localized to the optimizer: move set, initialization, or annealing landscape. There is no convention error to fix before running more optimization.

## Protocol

- family: minkowski only.
- dimensions: d_spacetime ∈ {2, 3, 4}.
- sizes: n ∈ {32, 64}.
- case seeds: 1959, 1962, 1987 (first three Phase 1B atlas seeds, same as Phase 2B).
- energy tolerance: |oracle_energy| ≤ 1e-09.
- RMSE tolerance: oracle_interval_rmse ≤ 1e-09.
- causal-matrix pass: total_discordant_pairs = 0.

Three oracle checks per row (no optimisation):

1. **oracle_pass_energy** — Bombelli energy at ground-truth
   coordinates is numerically zero. True iff the energy
   formula and the causal-matrix convention are mutually
   consistent.
2. **oracle_pass_causal_matrix** — causal matrix reconstructed
   from the stored coordinates matches the stored matrix
   bit-for-bit. True iff no floating-point drift or sign
   convention mismatch occurred between construction and
   storage.
3. **oracle_pass_interval_rmse** — Lorentz-invariant RMSE of
   the truth embedding against itself is numerically zero.
   True iff the interval-matrix formula is self-consistent.

## Results

| d | n | seed | oracle E | |oracle E| | RMSE | discordant | E pass | M pass | R pass |
| :---: | ---: | ---: | ---: | ---: | ---: | ---: | :---: | :---: | :---: |
| 2 | 32 | 1959 | 0.0000 | 0.0000 | 0.0000 | 0 | ✓ | ✓ | ✓ |
| 2 | 32 | 1962 | 0.0000 | 0.0000 | 0.0000 | 0 | ✓ | ✓ | ✓ |
| 2 | 32 | 1987 | 0.0000 | 0.0000 | 0.0000 | 0 | ✓ | ✓ | ✓ |
| 2 | 64 | 1959 | 0.0000 | 0.0000 | 0.0000 | 0 | ✓ | ✓ | ✓ |
| 2 | 64 | 1962 | 0.0000 | 0.0000 | 0.0000 | 0 | ✓ | ✓ | ✓ |
| 2 | 64 | 1987 | 0.0000 | 0.0000 | 0.0000 | 0 | ✓ | ✓ | ✓ |
| 3 | 32 | 1959 | 0.0000 | 0.0000 | 0.0000 | 0 | ✓ | ✓ | ✓ |
| 3 | 32 | 1962 | 0.0000 | 0.0000 | 0.0000 | 0 | ✓ | ✓ | ✓ |
| 3 | 32 | 1987 | 0.0000 | 0.0000 | 0.0000 | 0 | ✓ | ✓ | ✓ |
| 3 | 64 | 1959 | 0.0000 | 0.0000 | 0.0000 | 0 | ✓ | ✓ | ✓ |
| 3 | 64 | 1962 | 0.0000 | 0.0000 | 0.0000 | 0 | ✓ | ✓ | ✓ |
| 3 | 64 | 1987 | 0.0000 | 0.0000 | 0.0000 | 0 | ✓ | ✓ | ✓ |
| 4 | 32 | 1959 | 0.0000 | 0.0000 | 0.0000 | 0 | ✓ | ✓ | ✓ |
| 4 | 32 | 1962 | 0.0000 | 0.0000 | 0.0000 | 0 | ✓ | ✓ | ✓ |
| 4 | 32 | 1987 | 0.0000 | 0.0000 | 0.0000 | 0 | ✓ | ✓ | ✓ |
| 4 | 64 | 1959 | 0.0000 | 0.0000 | 0.0000 | 0 | ✓ | ✓ | ✓ |
| 4 | 64 | 1962 | 0.0000 | 0.0000 | 0.0000 | 0 | ✓ | ✓ | ✓ |
| 4 | 64 | 1987 | 0.0000 | 0.0000 | 0.0000 | 0 | ✓ | ✓ | ✓ |

All-three-pass count: 18 / 18.

## Five fixed questions

1. **Does the ground-truth embedding give energy zero?**
   Maximum |oracle_energy| across all 18 cases: 0. Every case satisfies oracle_pass_energy. The Bombelli energy formula returns exactly 0.0 at the ground-truth coordinates, confirming internal consistency between the causal matrix and the energy objective.

2. **Does the reconstructed causal matrix match the stored one?**
   Total discordant pairs across all rows: 0. The causal criterion used by the sprinkler and the criterion used in reconstruction are identical in floating-point. No drift or convention mismatch.

3. **Are the interval residuals zero at the ground truth?**
   Maximum oracle_interval_rmse: 0. The Lorentz-invariant RMSE of a point set against itself is zero, confirming that the interval-matrix formula is self-consistent.

4. **If something fails — what is the failure mode?**
   No failures in this run. All oracle checks pass. There is no evidence of a convention error in the energy formula, the causal criterion, or the interval metric.

5. **If oracle passes — what does that tell us about Phase 2B?**
   The energy objective recognises the ground-truth solution as a zero-energy configuration. The annealing failure in Phase 2/2B is therefore not caused by a broken energy formula. The optimizer does not find the geometry because of move-set or landscape issues, not because the target is wrong. The next diagnostic step is a move-set or initialization audit — not more budget, and not an energy redesign.

## Normalization audit

| d | n | seed | rave_truth | pair_count | original_pairs | reconstructed_pairs |
| :---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 2 | 32 | 1959 | 0.5156 | 496 | 268 | 268 |
| 2 | 32 | 1962 | 0.4983 | 496 | 272 | 272 |
| 2 | 32 | 1987 | 0.4467 | 496 | 225 | 225 |
| 2 | 64 | 1959 | 0.5231 | 2016 | 1013 | 1013 |
| 2 | 64 | 1962 | 0.4619 | 2016 | 956 | 956 |
| 2 | 64 | 1987 | 0.4510 | 2016 | 889 | 889 |
| 3 | 32 | 1959 | 0.4586 | 496 | 49 | 49 |
| 3 | 32 | 1962 | 0.5458 | 496 | 93 | 93 |
| 3 | 32 | 1987 | 0.5100 | 496 | 92 | 92 |
| 3 | 64 | 1959 | 0.4865 | 2016 | 373 | 373 |
| 3 | 64 | 1962 | 0.5116 | 2016 | 413 | 413 |
| 3 | 64 | 1987 | 0.5181 | 2016 | 323 | 323 |
| 4 | 32 | 1959 | 0.5045 | 496 | 47 | 47 |
| 4 | 32 | 1962 | 0.4926 | 496 | 81 | 81 |
| 4 | 32 | 1987 | 0.5129 | 496 | 48 | 48 |
| 4 | 64 | 1959 | 0.4983 | 2016 | 246 | 246 |
| 4 | 64 | 1962 | 0.5023 | 2016 | 249 | 249 |
| 4 | 64 | 1987 | 0.5303 | 2016 | 132 | 132 |

Regenerate via `make regen-phase2c`. Source tool:
`tools/build_phase2c_oracle_embedding_audit.py`.
