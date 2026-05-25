# Energy–F1 decoupling diagnostic  N=36

Post-run SORKIN-2 diagnostic.  Reads the oracle checkpoint ceiling CSV
and classifies each (schedule, seed) group by the sign of
`delta_energy_best_minus_final = best_checkpoint_energy_eave − final_energy_eave`
relative to `delta_best_minus_final` (causal F1 gain of the oracle checkpoint).

## Configuration

- Command: `python3 explore/energy_f1_decoupling_n36/run_energy_f1_decoupling_n36.py`
- Generated at UTC: `2026-05-25T09:27:28+00:00`
- Source CSV: `explore/oracle_checkpoint_ceiling_n36/oracle_checkpoint_ceiling_n36.csv`
- Main output CSV: `explore/energy_f1_decoupling_n36/energy_f1_decoupling_n36.csv`
- By-schedule CSV: `explore/energy_f1_decoupling_n36/energy_f1_decoupling_n36_by_schedule.csv`
- SVG: `explore/energy_f1_decoupling_n36/energy_f1_decoupling_n36.svg`
- This script does not run the annealer.
- Classification uses `causal_f1` against known-truth target → oracular, not deployable.

## Failure-mode classification (per group)

| seed | schedule | γ | budget | final F1 | best F1 | ΔF1 | final E | best E | ΔE | blk best | blk minE | best=minE | mode |
| ---: | --- | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- | --- |
| 1959 | gamma_0p5 | 0.5 | medium_25_25_8 | 0.4215 | 0.5488 | +0.1274 | 225.9 | 318.7 | +92.8 | 6 | 8 | no | **H2a** |
| 1962 | gamma_0p5 | 0.5 | medium_25_25_8 | 0.6096 | 0.6096 | +0.0000 | 220.3 | 220.3 | +0.0 | 8 | 8 | yes | inconclusive |
| 1987 | gamma_0p5 | 0.5 | medium_25_25_8 | 0.3732 | 0.5487 | +0.1755 | 215.6 | 314.6 | +98.9 | 6 | 8 | no | **H2a** |
| 2001 | gamma_0p5 | 0.5 | medium_25_25_8 | 0.3750 | 0.5842 | +0.2092 | 234.1 | 267.2 | +33.2 | 7 | 8 | no | **H2a** |
| 1959 | gamma_0p8 | 0.8 | medium_25_25_8 | 0.3742 | 0.4618 | +0.0876 | 438.4 | 493.9 | +55.5 | 7 | 1 | no | **H2a** |
| 1962 | gamma_0p8 | 0.8 | medium_25_25_8 | 0.4644 | 0.4644 | +0.0000 | 482.4 | 482.4 | +0.0 | 8 | 1 | no | inconclusive |
| 1987 | gamma_0p8 | 0.8 | medium_25_25_8 | 0.3872 | 0.4297 | +0.0425 | 447.1 | 548.5 | +101.4 | 5 | 1 | no | **H2a** |
| 2001 | gamma_0p8 | 0.8 | medium_25_25_8 | 0.3503 | 0.4225 | +0.0722 | 488.8 | 509.6 | +20.8 | 7 | 1 | no | **H2a** |
| 1959 | gamma_0p9 | 0.9 | medium_25_25_8 | 0.2976 | 0.3713 | +0.0737 | 509.6 | 318.8 | -190.9 | 1 | 1 | yes | **H2b** |
| 1962 | gamma_0p9 | 0.9 | medium_25_25_8 | 0.2687 | 0.3901 | +0.1214 | 580.4 | 583.9 | +3.4 | 5 | 1 | no | **H2a** |
| 1987 | gamma_0p9 | 0.9 | medium_25_25_8 | 0.3762 | 0.4473 | +0.0711 | 504.9 | 473.3 | -31.7 | 2 | 1 | no | **H2b** |
| 2001 | gamma_0p9 | 0.9 | medium_25_25_8 | 0.2575 | 0.3851 | +0.1276 | 521.4 | 538.6 | +17.2 | 7 | 1 | no | **H2a** |
| 1959 | gamma_0p95 | 0.95 | medium_25_25_8 | 0.2353 | 0.3713 | +0.1360 | 547.2 | 318.8 | -228.4 | 1 | 1 | yes | **H2b** |
| 1962 | gamma_0p95 | 0.95 | medium_25_25_8 | 0.3229 | 0.4038 | +0.0810 | 593.0 | 632.2 | +39.2 | 3 | 1 | no | **H2a** |
| 1987 | gamma_0p95 | 0.95 | medium_25_25_8 | 0.3543 | 0.4473 | +0.0929 | 598.1 | 473.3 | -124.8 | 2 | 1 | no | **H2b** |
| 2001 | gamma_0p95 | 0.95 | medium_25_25_8 | 0.2609 | 0.3760 | +0.1151 | 563.1 | 281.7 | -281.4 | 1 | 1 | yes | **H2b** |

## By-schedule aggregates

| schedule | γ | budget | N | avg ΔE | med ΔE | avg ΔF1 | avg best F1 | avg final F1 | H2a | H2b | inconc | best=minE | avg blk best | avg blk minE |
| --- | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| gamma_0p5 | 0.5 | medium_25_25_8 | 4 | +56.2 | +63.0 | 0.1280 | 0.5728 | 0.4448 | 3 | 0 | 1 | 1 | 6.75 | 8.00 |
| gamma_0p8 | 0.8 | medium_25_25_8 | 4 | +44.4 | +38.1 | 0.0506 | 0.4446 | 0.3940 | 3 | 0 | 1 | 0 | 6.75 | 1.00 |
| gamma_0p9 | 0.9 | medium_25_25_8 | 4 | -50.5 | -14.1 | 0.0984 | 0.3984 | 0.3000 | 2 | 2 | 0 | 1 | 3.75 | 1.00 |
| gamma_0p95 | 0.95 | medium_25_25_8 | 4 | -148.9 | -176.6 | 0.1062 | 0.3996 | 0.2933 | 1 | 3 | 0 | 2 | 1.75 | 1.00 |

## Diagnostic questions

### Q1 — Is there a clear separation between fast and slow schedules?

**Yes, with conservative support.**  Fast schedules (gamma_0p5, gamma_0p8): 6/8 groups → H2a, avg ΔE = +50.3.  Slow schedules (gamma_0p9, gamma_0p95): 5/8 groups → H2b, avg ΔE = -99.7.  The sign of ΔE is predominantly positive for fast schedules and predominantly negative for slow schedules.  With only 4 groups per schedule the separation is suggestive, not conclusive.

### Q2 — How many groups fall into H2a?

**9 of 16 groups** are classified as `H2a_over_annealing_candidate`.

### Q3 — How many groups fall into H2b?

**5 of 16 groups** are classified as `H2b_escape_or_nonconvergence_candidate`.
**2 of 16 groups** are `inconclusive` (best == final or ΔE = 0).

### Q4 — Is the classification by sign of ΔE robust or marginal?

3 group(s) have |ΔE| < 5 energy units (modes: ['inconclusive', 'inconclusive', 'H2a_over_annealing_candidate']).  These are near the boundary and their classification is fragile.

### Q5 — Is it still reasonable to treat the 16 groups as a single population?

**No.**  The failure-mode split correlates strongly with cooling speed.  Pooling all 16 groups flattens two qualitatively distinct dynamics: one in which the annealer over-cools through the good causal region (H2a) and one in which it visits then escapes that region (H2b).  Treating the 16 groups as a homogeneous sample would produce averages that describe neither regime accurately.

### Q6 — What does this imply for the next probe?

If the H2a/H2b separation is confirmed:

- **H2a probes** should focus on the temperature window during which causal F1 peaks and whether it corresponds to a stable or transient basin.  A trajectory-dump probe for fast schedules (gamma_0p5 / gamma_0p8) recording causal F1 at each block across multiple seeds would establish whether the good-causal temperature range is reproducible.

- **H2b probes** should focus on why the annealer escapes from the initially-good state.  This is a different failure: the energy landscape is not monotone in causal quality even at the start.  A probe recording which relations flip between the best block and later blocks would characterise the escape mechanism.

- **Do not design a single stopping criterion** that targets both regimes simultaneously.  The two regimes require different interventions.

## Conservative interpretation

This diagnostic uses `causal_f1` against the known-truth partial order to identify the oracle-best checkpoint.  It is therefore not a deployable selection criterion for truth-free cases.

The H2a / H2b classification is a hypothesis with `N_groups = 16` and only 4 seeds per schedule.  It may reflect true regime structure or schedule-specific noise.  The classification should be treated as a diagnostic candidate, not a confirmed physical finding.

The diagnostic does not distinguish between:
  - a fundamentally misaligned energy function (same function for all regimes),
  - schedule-specific convergence behaviour (different dynamics per cooling rate).
Both could produce the observed pattern.  Next probes should hold the energy function fixed (it is historical) and vary the observable: trajectory timing, relation-flip rate, or seed-to-seed configuration similarity.

## Guardrails

This is a post-run diagnostic only, over benchmark cases with known truth.
It is not an embeddability claim, not a physical gamma claim, not an N-transition claim,
and not proof of general annealer failure.
It is not a deployable criterion for truth-free cases.
