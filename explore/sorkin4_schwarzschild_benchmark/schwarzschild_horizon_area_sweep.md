# Schwarzschild Horizon Area Sweep

This is an aligned radial-strand SORKIN-4 diagnostic.
It is not a 4D area-law measurement.

- Mass values: [0.75, 1.0, 1.25, 1.5, 1.75, 2.0]
- Seeds: 1..40
- N exterior/interior: 16/8
- Mean horizon-link counts monotone non-decreasing: True
- Failed order checks: 0

| M | A=16piM^2 | mean horizon links | std | min | max | nonzero seeds |
|---:|---:|---:|---:|---:|---:|---:|
| 0.75 | 28.2743 | 2.775 | 2.37158 | 0 | 10 | 33 |
| 1.00 | 50.2655 | 3.375 | 2.6428 | 0 | 11 | 36 |
| 1.25 | 78.5398 | 3.8 | 2.63818 | 0 | 12 | 37 |
| 1.50 | 113.097 | 4.575 | 2.78287 | 0 | 12 | 39 |
| 1.75 | 153.938 | 4.9 | 2.498 | 1 | 12 | 40 |
| 2.00 | 201.062 | 5.05 | 2.74727 | 1 | 12 | 40 |

Interpretation:

- `horizon_crossing_links` counts exterior-to-interior links in the transitive reduction.
- The aligned radial setup gives a controlled 1+1-like diagnostic with no non-radial undecided pairs.
- The result checks monotonic response to increasing formal horizon area, not the 4D area coefficient.
