# Known-truth N=20

This directory freezes a small manufactured known-truth case for SORKIN-2
recoverability diagnostics. It is not a community golden table.

Configuration:

- `N = 20`
- `family = minkowski`
- `d_spacetime = 2`
- `case_seed = 1959`
- geometry: Minkowski 1+1 unit causal diamond
- generator: `validation_suite.sprinkle_minkowski_diamond`

The generator writes:

- `case_metadata.json`
- `ground_truth_coords.csv`
- `target_pairs.csv`
- `target_matrix.csv`
- `target_invariants.json`
- `checksums_sha256.txt`

Regenerate the frozen artifacts with:

```bash
python3 explore/known_truth_n20/generate_known_truth_n20.py
```

The target causal order and coordinates are deterministic for the constants
above. These files are intended as an internal known-truth benchmark for
annealer diagnostics, not as evidence of embeddability beyond construction.
