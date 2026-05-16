"""Regression tests for the foundation benchmark fixtures.

The foundation benchmark is a frozen grid of canonical Minkowski-diamond
sprinklings together with their precomputed order-theoretic invariants
(``benchmarks/foundation/invariants.json``). These tests verify that
the current implementation of the sprinkler and the invariants
exactly reproduces the frozen values, catching any unintended drift
in the scientific protocol.
"""

from __future__ import annotations

import json
import unittest
from pathlib import Path

import causet_invariants
import validation_suite as vs


ROOT = Path(__file__).resolve().parents[1]
FOUNDATION = ROOT / "benchmarks" / "foundation"


class FoundationBenchmarkIntegrityTests(unittest.TestCase):
    def setUp(self) -> None:
        self.payload = json.loads(
            (FOUNDATION / "invariants.json").read_text(encoding="utf-8")
        )

    def test_cells_cover_expected_grid(self) -> None:
        cells = self.payload["cells"]
        self.assertEqual(len(cells), 45)
        dims = {c["d_spacetime"] for c in cells}
        sizes = {c["n"] for c in cells}
        self.assertEqual(dims, {2, 3, 4})
        self.assertEqual(sizes, {16, 32, 64})

    def test_frozen_invariants_match_recomputation(self) -> None:
        for cell in self.payload["cells"]:
            d = cell["d_spacetime"]
            n = cell["n"]
            seed = cell["seed"]
            matrix, _ = vs.sprinkle_minkowski_diamond(
                n=n, seed=seed, d_spacetime=d
            )
            fp = causet_invariants.invariants_fingerprint(matrix)
            # JSON serialization turns int keys in chain_counts into
            # strings; normalize the recomputed fingerprint for
            # comparison.
            fp_normalized = dict(fp)
            fp_normalized["chain_counts"] = {
                str(k): v for k, v in fp["chain_counts"].items()
            }
            self.assertEqual(
                fp_normalized,
                cell["fingerprint"],
                msg=(
                    f"Foundation cell (d={d}, n={n}, seed={seed}) "
                    f"invariants drifted from frozen values"
                ),
            )

    def test_mean_myrheim_meyer_per_dimension_is_consistent(self) -> None:
        # Aggregate MM dim estimates per d_spacetime across all cells
        # and verify that the mean is within Monte Carlo tolerance of
        # the true dimension. This is a sanity check on the entire
        # pipeline: sprinkler, ordering fraction, and MM inversion all
        # have to be correct for this to pass.
        by_dim: dict[int, list[float]] = {}
        for cell in self.payload["cells"]:
            by_dim.setdefault(cell["d_spacetime"], []).append(
                cell["fingerprint"]["myrheim_meyer_dim"]
            )
        for d, ests in by_dim.items():
            mean = sum(ests) / len(ests)
            self.assertAlmostEqual(
                mean, float(d), delta=0.6,
                msg=(
                    f"d_spacetime={d}: mean MM dim {mean:.3f} too far "
                    f"from {d}"
                ),
            )


if __name__ == "__main__":
    unittest.main()
