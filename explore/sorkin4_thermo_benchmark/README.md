# SORKIN-4 Thermo Benchmark

This directory contains conservative **Level-A** closed-form guardrails for
Schwarzschild/Kerr horizon geometry and standard thermodynamic scalars.

Scope:

- closed-form identity and trend checks only;
- regression protection for formulas already used in the benchmark sequence;
- no claim of new Hawking/Bekenstein physics.

Level-A vs Level-B:

- Level-A here means algebraic/closed-form consistency checks.
- Level-B would require discrete causal-set pipeline observables (for example:
  links, horizon crossings, horizon molecules, or horizon counts) to
  reconstruct continuum horizon laws without hard-coding those laws as outputs.

Fixed-M sweep caveat:

At fixed `M`, sweeping Kerr spin compares different stationary solutions.
Decreases in area/entropy along that sweep are not physical dynamical
evolution and do not violate area-theorem/second-law statements.
