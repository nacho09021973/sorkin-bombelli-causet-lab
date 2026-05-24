# SORKIN-2 exploratory probe index

Status: exploratory inventory; not a result paper; not conservative documentation.

## Current probe inventory

### gamma_n24_probe

- status: completed exploratory N=24 gamma probe plus post-analysis.
- N: 24.
- gamma coverage: 11-point logarithmic gamma grid from 0.5 to 0.9.
- T0: 100.
- seeds: case seed 1959; optimizer seed 1987.
- outputs present: README.md, runner script, CSV output, generated Markdown summary, SVG plot, ranking diagnostic CSV/Markdown, and energy-vs-RMSE SVG.
- main exploratory lesson: final energy selected gamma=0.562373, interval RMSE selected gamma=0.848624, rankings disagree, all success flags are false, and this does not establish a robust gamma optimum.

### gamma_grid_stress_n24

- status: completed exploratory N=24 gamma-grid stress test.
- N: 24.
- gamma coverage: 193 gamma grid points across 8 grid families.
- T0: 100.
- seeds: case seed 1959; optimizer seed 1987.
- outputs present: README.md, runner script, CSV output, generated Markdown summary, and SVG summary figure.
- main exploratory lesson: no robust gamma optimum is established; winners are grid- and metric-dependent, and all success flags are false.

### causal_order_metric_probe_n24

- status: completed exploratory N=24 causal-order metric probe.
- N: 24.
- gamma coverage: 0.55, 0.73, 0.82.
- T0: 100.
- seeds: case seed 1959; optimizer seed 1987.
- outputs present: runner script, CSV output, generated Markdown summary, and SVG figure.
- main exploratory lesson: final energy and causal F1 selected gamma=0.55, interval RMSE selected gamma=0.73, and no exact/success outcome was reported.

### causal_order_seed_stability_n24

- status: completed exploratory N=24 seed-stability probe.
- N: 24.
- gamma coverage: 0.50, 0.53, 0.55, 0.57, 0.60.
- T0: 100.
- seeds: case seed 1959; optimizer seeds 1959, 1962, 1987, 2001, 2026.
- outputs present: runner script, CSV output, generated Markdown summary, and SVG figure.
- main exploratory lesson: gamma=0.50 was best by mean, median, and worst-case causal F1 across the tested optimizer seeds, with no exact/success outcome.

### causal_order_budget_scaling_n24

- status: completed exploratory N=24 budget/accessibility probe in the existing artifact.
- N: 24.
- gamma coverage: fixed gamma=0.5.
- T0: 100.
- seeds: case seed 1959; optimizer seeds 1959, 1962, 1987, 2001, 2026.
- outputs present: runner script, runner backup files, CSV output, generated Markdown summary, SVG figure, and log files.
- main exploratory lesson: at gamma=0.5, the long_50_50_16 budget rows show near-exact causal-order recovery, with mean F1 0.993610, mean recall 0.992405, mean precision 0.994935, 6 missing relations total, 4 extra relations total, exact 2/5, and success 5/5.

### gamma_n36_probe

- status: script-only/deferred exploratory N=36 attempt.
- N: 36.
- gamma coverage: runner script states the same 11-point grid from 0.5 to 0.9.
- T0: 100.
- seeds: runner script states case seed 1959 and optimizer seed 1987.
- outputs present: README.md and run_gamma_n36_probe.py only; no CSV, no generated Markdown summary, no plot, and no completed numerical result.
- main exploratory lesson: no completed N=36 numerical result is present; this folder should not be cited as evidence for an N=36 outcome.

## Current lessons

- N=24 gamma-only scans do not currently support a robust privileged gamma near 0.8.
- Different metrics select different gamma values.
- The strongest current exploratory signal is budget/accessibility, not gamma alone.
- At N=24 and gamma=0.5, increasing the annealer budget to long_50_50_16 produced near-exact causal-order recovery in the existing artifact.
- This suggests that recoverability may depend on an effective combination of gamma, T0, budget, topology, and seed.
- The earlier gamma resonance hypothesis should remain alive as a hypothesis generator, but should not be treated as supported by the N=24 probes.

## Open questions

- Does N=36 show the same budget/accessibility behavior?
- Is there any gamma window once budget is controlled?
- Is gamma=0.8 useful only under specific budget constraints?
- Can a thermal mobility quantity combine gamma and budget better than gamma alone?
- Are failures dominated by missing relations, extra relations, or metric mismatch?

## Guardrails

- This file is exploratory.
- It is not evidence for a physical constant.
- It is not evidence for embeddability.
- It is not a causal-set theorem.
- Do not cite gamma=0.8 as privileged from these probes.
- Any promoted claim must cite named runs, CSVs, figures, and scripts.
- Any future N=36 claim must be based on completed N=36 outputs, not on the deferred script-only folder.
