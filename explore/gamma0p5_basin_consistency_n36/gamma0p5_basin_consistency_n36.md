# gamma_0p5 basin consistency  N=36

Post-run SORKIN-2 probe.  Runs gamma_0p5 / medium_25_25_8 for N=36
and computes pairwise Jaccard similarity of induced causal orders at
the H2a peak checkpoints (blocks 6–7) and block-1 controls.

## Configuration

- Command: `python3 explore/gamma0p5_basin_consistency_n36/run_gamma0p5_basin_consistency_n36.py`
- Generated at UTC: `2026-05-25T10:48:55+00:00`
- N: `36`, case_seed: `1959`, d_spacetime: `2`
- Schedule: `gamma_0p5` (cooling_factor=0.5)
- Budget: `medium_25_25_8` (warmup=25, anneal=25, blocks=8)
- Optimizer seeds: `1959, 1962, 1987, 2001`
- Target causal relations: `343` of 630 possible pairs.
- Classification uses `causal_f1` against known truth → oracular, not deployable.
- Jaccard threshold basin compartido: `> 0.75`.
- Jaccard threshold basins distintos: `< 0.5`.

## Reproducibility check

- Blocks compared against `explore/schedule_seed_stability_n36/schedule_seed_stability_n36.csv`: `32`.
- max |Δ causal_f1|:   `4.898e-11`  (tolerance 1e-06).
- max |Δ energy_eave|: `4.988e-08`.
- **Reproducible: `True`.**

Run is reproducible within floating-point tolerance.  Jaccard analysis proceeds.

## Jaccard matrix (full 8×8)

Rows/columns ordered: H2a peaks, inconclusive reference, blk-1 controls.
Upper triangle shown (matrix is symmetric; diagonal = 1.000).

| | s1959_b6 | s1987_b6 | s2001_b7 | s1962_b8 | s1959_b1 | s1987_b1 | s2001_b1 | s1962_b1 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| s1959_b6 | 1.000 | 0.251 | 0.256 | 0.253 | 0.186 | 0.173 | 0.198 | 0.219 |
| s1987_b6 | 0.251 | 1.000 | 0.264 | 0.333 | 0.257 | 0.137 | 0.161 | 0.202 |
| s2001_b7 | 0.256 | 0.264 | 1.000 | 0.331 | 0.151 | 0.178 | 0.201 | 0.088 |
| s1962_b8 | 0.253 | 0.333 | 0.331 | 1.000 | 0.230 | 0.162 | 0.241 | 0.174 |
| s1959_b1 | 0.186 | 0.257 | 0.151 | 0.230 | 1.000 | 0.141 | 0.064 | 0.114 |
| s1987_b1 | 0.173 | 0.137 | 0.178 | 0.162 | 0.141 | 1.000 | 0.124 | 0.163 |
| s2001_b1 | 0.198 | 0.161 | 0.201 | 0.241 | 0.064 | 0.124 | 1.000 | 0.119 |
| s1962_b1 | 0.219 | 0.202 | 0.088 | 0.174 | 0.114 | 0.163 | 0.119 | 1.000 |

## Checkpoint details

| checkpoint | type | block | T | causal_f1 | induced_size |
| --- | --- | ---: | ---: | ---: | ---: |
| s1959_b6 | H2a_peak | 6 | 3.12 | 0.5488 | 251 |
| s1987_b6 | H2a_peak | 6 | 3.12 | 0.5487 | 273 |
| s2001_b7 | H2a_peak | 7 | 1.56 | 0.5842 | 215 |
| s1962_b8 | inconclusive | 8 | 0.78 | 0.6096 | 264 |
| s1959_b1 | control_blk1 | 1 | 100.00 | 0.3713 | 158 |
| s1987_b1 | control_blk1 | 1 | 100.00 | 0.2920 | 109 |
| s2001_b1 | control_blk1 | 1 | 100.00 | 0.3760 | 173 |
| s1962_b1 | control_blk1 | 1 | 100.00 | 0.3534 | 155 |

## Diagnostic questions

### Q1 — Do the H2a checkpoints share a causal basin?

**Basins distintos.**  Jaccard medio entre H2a peaks = 0.257 < 0.5.  El F1 parecido en bloques 6–7 no corresponde a configuraciones similares.  Un stopping criterion por temperatura seleccionaría basins distintos según el seed.

### Q2 — Jaccard level between H2a peaks: high, medium, or low?

H2a pairwise Jaccard values: 0.251, 0.256, 0.264.
Mean = **0.257**, min = 0.251, max = 0.264.

### Q3 — Was block-1 Jaccard already high?

**No.**  Mean block-1 cross-seed Jaccard = 0.121 < 0.5.  The hot-start configurations are structurally dissimilar.  The Jaccard elevation at blocks 6–7 (avg 0.257) reflects genuine convergence, not pre-existing similarity.

### Q4 — Does seed 1962 block 8 resemble the H2a peaks?

**Low similarity.**  seed 1962 blk-8 vs H2a peaks: mean Jaccard = 0.305.  seed 1962 final endpoint is structurally different from H2a peak configurations, consistent with reaching a different local minimum.

### Q5 — Does T ≈ 1.5–3.5 have structural support or only scalar F1 support?

**Unclear.**  The Jaccard result (0.257) cannot be cleanly separated from background similarity.  The T ≈ 1.5–3.5 window has only scalar (F1) support in this dataset.

### Q6 — Conservative conclusion

Verdict: **`basins_distintos`**.

Details:
- H2a peak Jaccard (3 seeds, N=36, 1 case): avg=0.257, range=[0.251, 0.264].
- Control blk-1 Jaccard (baseline): avg=0.121.
- seed 1962 blk-8 vs H2a peaks: avg=0.305.

The F1 similarity at blocks 6–7 is not accompanied by structural (causal order) similarity.
A stopping rule targeting T ≈ 3 would recover configurations that differ between seeds.
This closes the temperature-based stopping criterion hypothesis for gamma_0p5 / N=36.

## Conservative interpretation

This is an oracular diagnostic: the Jaccard computation uses the known-truth
induced order, not an observable available without ground truth.

Sample: 3 H2a seeds, 1 case seed (1959), 1 N (36), 1 schedule (gamma_0p5).
No generalization claim is warranted from this probe alone.

## Guardrails

This is a post-run diagnostic only, over one benchmark case with known truth.
It is not an embeddability claim, not a physical gamma claim, not an N-transition claim,
and not proof of general annealer failure.
It is not a deployable checkpoint-selection criterion.
