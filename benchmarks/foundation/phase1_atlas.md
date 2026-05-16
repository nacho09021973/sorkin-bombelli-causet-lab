# Phase 1 Order-Theoretic Atlas

Comparison of two independent order-theoretic dimension
estimators on canonical Minkowski-diamond sprinklings
(manifoldlike) and Kleitman-Rothschild three-level posets
(non-manifoldlike controls). No embedding, energy, or
annealing is invoked in this comparison.

Columns:

- **family**: `minkowski` or `kleitman_rothschild`.
- **d_spacetime**: the sprinkling dimension for Minkowski
  causets; blank for KR (no embedding dimension is defined).
- **n**: number of events.
- **mean_mm**: ensemble mean of the Myrheim-Meyer dimension
  over all seeds in the row.
- **mean_midpoint**: ensemble mean of Meyer's midpoint
  scaling dimension.
- **mean_discrepancy**: ensemble mean of ``|mm - midpoint|``.

Manifoldlike rows should show `mean_mm ~ mean_midpoint ~ d_spacetime`.
Non-manifoldlike rows show the diagnostic separation.

| family | d | n | seeds | mean MM | mean midpoint | mean discrepancy |
| --- | :---: | ---: | ---: | ---: | ---: | ---: |
| kleitman_rothschild | - | 16 | 5 | 2.3962 | 1.4797 | 0.9164 |
| kleitman_rothschild | - | 32 | 5 | 2.3431 | 2.2582 | 0.1684 |
| kleitman_rothschild | - | 64 | 5 | 2.3611 | 3.0846 | 0.7235 |
| minkowski | 2 | 16 | 5 | 1.9248 | 1.4721 | 0.4905 |
| minkowski | 2 | 32 | 5 | 1.9301 | 1.6546 | 0.2754 |
| minkowski | 2 | 64 | 5 | 2.0144 | 1.9074 | 0.1574 |
| minkowski | 3 | 16 | 5 | 3.0662 | 1.6654 | 1.3355 |
| minkowski | 3 | 32 | 5 | 3.3164 | 1.7428 | 1.5735 |
| minkowski | 3 | 64 | 5 | 3.2151 | 2.3550 | 0.8855 |
| minkowski | 4 | 16 | 5 | 3.7443 | 1.0000 | 2.6563 |
| minkowski | 4 | 32 | 5 | 3.7065 | 1.5546 | 2.1519 |
| minkowski | 4 | 64 | 5 | 3.9174 | 2.4853 | 1.4321 |

Regenerate with `make regen-phase1`. The underlying tool is
`tools/build_phase1_atlas.py`.
