"""Integrity tests for the Phase 1B finite-size scaling atlas.

The artifact at ``benchmarks/foundation/phase1b_scaling_atlas.csv``
must exist with the expected schema, must cover every cell of the
``(family, target_dim if applicable, n, seed)`` grid declared in
``tools/build_phase1b_scaling_atlas.py``, and must contain only
finite numeric values in the ``mm_dim`` and ``midpoint_dim``
columns. NaN or infinity in those columns at this size range
would indicate the estimators are silently failing on a cell.
"""

from __future__ import annotations

import csv
import math
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ATLAS_CSV = ROOT / "benchmarks" / "foundation" / "phase1b_scaling_atlas.csv"
ATLAS_MD = ROOT / "benchmarks" / "foundation" / "phase1b_scaling_atlas.md"

EXPECTED_HEADERS = (
    "family",
    "d_spacetime",
    "n",
    "seed",
    "mm_dim",
    "midpoint_dim",
    "abs_discrepancy",
)

EXPECTED_FAMILIES = {"minkowski", "kleitman_rothschild"}
EXPECTED_SIZES = {32, 64, 128, 256}
EXPECTED_MINKOWSKI_DIMS = {"2", "3", "4"}
EXPECTED_SEEDS = {1959, 1962, 1987, 2009, 2026}


class Phase1BScalingAtlasIntegrityTests(unittest.TestCase):
    def setUp(self) -> None:
        self.assertTrue(
            ATLAS_CSV.exists(),
            msg=(
                f"missing atlas CSV at {ATLAS_CSV}; "
                "run `make regen-phase1b`"
            ),
        )
        with ATLAS_CSV.open(newline="", encoding="utf-8") as fh:
            self.rows = list(csv.DictReader(fh))

    def test_markdown_exists(self) -> None:
        self.assertTrue(
            ATLAS_MD.exists(),
            msg=(
                f"missing atlas markdown at {ATLAS_MD}; "
                "run `make regen-phase1b`"
            ),
        )

    def test_columns_match_expected(self) -> None:
        with ATLAS_CSV.open(newline="", encoding="utf-8") as fh:
            reader = csv.reader(fh)
            header = tuple(next(reader))
        self.assertEqual(header, EXPECTED_HEADERS)

    def test_row_count_and_families(self) -> None:
        # 4 sizes x 5 seeds x (3 Minkowski d's + 1 KR) = 80 rows.
        self.assertEqual(len(self.rows), 80)
        self.assertEqual(
            {r["family"] for r in self.rows}, EXPECTED_FAMILIES
        )

    def test_covers_full_grid(self) -> None:
        # Every Minkowski (d, n, seed) cell must be present.
        minkowski_keys = {
            (r["d_spacetime"], int(r["n"]), int(r["seed"]))
            for r in self.rows if r["family"] == "minkowski"
        }
        expected_minkowski = {
            (d, n, s)
            for d in EXPECTED_MINKOWSKI_DIMS
            for n in EXPECTED_SIZES
            for s in EXPECTED_SEEDS
        }
        self.assertEqual(minkowski_keys, expected_minkowski)
        # Every KR (n, seed) cell must be present.
        kr_keys = {
            (int(r["n"]), int(r["seed"]))
            for r in self.rows if r["family"] == "kleitman_rothschild"
        }
        expected_kr = {(n, s) for n in EXPECTED_SIZES for s in EXPECTED_SEEDS}
        self.assertEqual(kr_keys, expected_kr)

    def test_no_nan_or_inf_in_dimension_columns(self) -> None:
        # At n >= 32 both estimators should always produce finite
        # values on both families. A NaN here means a generator or
        # estimator silently failed on a cell.
        for row in self.rows:
            for col in ("mm_dim", "midpoint_dim"):
                raw = row[col]
                self.assertNotIn(
                    raw.lower(), {"nan", "inf", "-inf"},
                    msg=(
                        f"{col}={raw} for cell "
                        f"family={row['family']}, d={row['d_spacetime']}, "
                        f"n={row['n']}, seed={row['seed']}"
                    ),
                )
                value = float(raw)
                self.assertTrue(
                    math.isfinite(value),
                    msg=(
                        f"non-finite {col}={value} for cell "
                        f"family={row['family']}, n={row['n']}, "
                        f"seed={row['seed']}"
                    ),
                )


if __name__ == "__main__":
    unittest.main()
