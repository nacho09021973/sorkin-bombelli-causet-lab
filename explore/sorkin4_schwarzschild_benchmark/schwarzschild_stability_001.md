# S4-SCHW-STABILITY-001: Exterior Schwarzschild Stability Sweep

Generated: 2026-05-27T12:20:48.235630+00:00

## Purpose

Exploratory diagnostic. Sweeps `N`, `seed`, and `margin` to check whether
the exterior Schwarzschild benchmark remains well-behaved outside the single
frozen reference case (N=12, seed=1959, margin=0.35).

## Sweep Parameters

- N: [12, 16]
- seed: [1959, 1960, 1961]
- margin: [0.25, 0.35, 0.5]
- Total cells: 18

## Aggregate Results

| Check | Result |
|-------|--------|
| Q1: all order checks pass | True |
| Q2: turning branch no invasion | True |
| Q3: min turning gap (all cells) | 0.14016744672036555 |
| Q3: max turning gap (all cells) | 0.1666184302832776 |
| Q3: gap consistent positive | True |
| Q4: near-inner True assertions (total) | 0 |
| Q4: near-outer True assertions (total) | 1 |
| **all_checks_pass** | **True** |

## Per-Cell Results

| N | seed | margin | true | undecided | antisym | transitive | gaps | min_gap | local_pass |
|---|------|--------|------|-----------|---------|------------|------|---------|------------|
| 12 | 1959 | 0.25 | 1 | 0 | True | True | 16 | 0.1666 | True |
| 12 | 1959 | 0.35 | 1 | 0 | True | True | 16 | 0.1658 | True |
| 12 | 1959 | 0.50 | 1 | 0 | True | True | 16 | 0.1645 | True |
| 12 | 1960 | 0.25 | 0 | 1 | True | True | 28 | 0.1480 | True |
| 12 | 1960 | 0.35 | 0 | 1 | True | True | 28 | 0.1469 | True |
| 12 | 1960 | 0.50 | 0 | 1 | True | True | 32 | 0.1453 | True |
| 12 | 1961 | 0.25 | 3 | 0 | True | True | 32 | 0.1444 | True |
| 12 | 1961 | 0.35 | 3 | 0 | True | True | 32 | 0.1428 | True |
| 12 | 1961 | 0.50 | 3 | 0 | True | True | 33 | 0.1402 | True |
| 16 | 1959 | 0.25 | 1 | 0 | True | True | 37 | 0.1666 | True |
| 16 | 1959 | 0.35 | 1 | 0 | True | True | 37 | 0.1658 | True |
| 16 | 1959 | 0.50 | 1 | 0 | True | True | 37 | 0.1645 | True |
| 16 | 1960 | 0.25 | 1 | 2 | True | True | 56 | 0.1480 | True |
| 16 | 1960 | 0.35 | 1 | 2 | True | True | 56 | 0.1469 | True |
| 16 | 1960 | 0.50 | 0 | 3 | True | True | 62 | 0.1453 | True |
| 16 | 1961 | 0.25 | 4 | 1 | True | True | 57 | 0.1444 | True |
| 16 | 1961 | 0.35 | 3 | 2 | True | True | 57 | 0.1428 | True |
| 16 | 1961 | 0.50 | 3 | 2 | True | True | 58 | 0.1402 | True |

## Notes

- `gaps` = number of outgoing pairs for which a valid turning-branch gap was computed.
- `min_gap` = minimum of (phi_turning_min − phi_direct_max) over those pairs; a positive
  value means the turning branch cannot reach any target that the direct branch can reach.
- Q4 counts near-margin True assertions using a 10 % domain-width threshold; these are
  informational only and do not affect `all_checks_pass`.
- This diagnostic does not claim embeddability, manifoldlikeness, or physical correctness.
  It only reports algorithmic stability of the partial He-Rideout exterior model across
  the swept parameters.
