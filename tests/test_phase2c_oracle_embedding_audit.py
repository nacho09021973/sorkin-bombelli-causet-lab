"""Integrity tests for the Phase 2C oracle embedding audit artifact."""

from __future__ import annotations

import csv
import math
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ORACLE_CSV = ROOT / "benchmarks" / "foundation" / "phase2c_oracle_embedding_audit.csv"
ORACLE_MD = ROOT / "benchmarks" / "foundation" / "phase2c_oracle_embedding_audit.md"

EXPECTED_HEADERS = (
    "family",
    "target_dim",
    "n",
    "seed",
    "d_spatial",
    "pair_count",
    "oracle_energy",
    "oracle_energy_abs",
    "rave_truth",
    "oracle_interval_rmse",
    "original_pair_count",
    "reconstructed_pair_count",
    "false_positive_pairs",
    "false_negative_pairs",
    "total_discordant_pairs",
    "oracle_pass_energy",
    "oracle_pass_causal_matrix",
    "oracle_pass_interval_rmse",
    "notes",
)


class Phase2COracleAuditTests(unittest.TestCase):
    def setUp(self) -> None:
        self.assertTrue(
            ORACLE_CSV.exists(),
            msg=(
                f"missing Phase 2C CSV at {ORACLE_CSV}; "
                "run `make regen-phase2c`"
            ),
        )
        with ORACLE_CSV.open(newline="", encoding="utf-8") as fh:
            self.rows = list(csv.DictReader(fh))

    def test_markdown_exists(self) -> None:
        self.assertTrue(
            ORACLE_MD.exists(),
            msg=(
                f"missing Phase 2C markdown at {ORACLE_MD}; "
                "run `make regen-phase2c`"
            ),
        )

    def test_columns_match_expected(self) -> None:
        with ORACLE_CSV.open(newline="", encoding="utf-8") as fh:
            reader = csv.reader(fh)
            header = tuple(next(reader))
        self.assertEqual(header, EXPECTED_HEADERS)

    def test_has_at_least_one_row(self) -> None:
        self.assertGreater(len(self.rows), 0)

    def test_only_minkowski_family(self) -> None:
        families = {row["family"] for row in self.rows}
        self.assertEqual(
            families,
            {"minkowski"},
            msg=(
                "Phase 2C must not include KR or corona controls; "
                f"got families={families}"
            ),
        )

    def test_no_kr_or_corona(self) -> None:
        for row in self.rows:
            self.assertNotIn(
                row["family"],
                {"kleitman_rothschild", "corona_poset"},
                msg=f"unexpected control family in oracle: {row['family']}",
            )

    def test_oracle_energy_is_finite(self) -> None:
        for row in self.rows:
            raw = row["oracle_energy"]
            self.assertNotIn(raw.lower(), {"nan", "inf", "-inf", "na"})
            value = float(raw)
            self.assertTrue(math.isfinite(value), msg=f"oracle_energy={raw}")

    def test_oracle_interval_rmse_is_finite(self) -> None:
        for row in self.rows:
            raw = row["oracle_interval_rmse"]
            self.assertNotIn(raw.lower(), {"nan", "inf", "-inf", "na"})
            value = float(raw)
            self.assertTrue(math.isfinite(value), msg=f"oracle_interval_rmse={raw}")

    def test_total_discordant_pairs_is_nonneg_int(self) -> None:
        for row in self.rows:
            value = int(row["total_discordant_pairs"])
            self.assertGreaterEqual(value, 0)

    def test_zero_discordant_implies_pass_causal_matrix(self) -> None:
        for row in self.rows:
            disc = int(row["total_discordant_pairs"])
            flag = row["oracle_pass_causal_matrix"]
            if disc == 0:
                self.assertEqual(
                    flag,
                    "true",
                    msg=(
                        f"total_discordant_pairs=0 but oracle_pass_causal_matrix={flag!r} "
                        f"for {row}"
                    ),
                )

    def test_pass_flags_are_boolean(self) -> None:
        bool_cols = (
            "oracle_pass_energy",
            "oracle_pass_causal_matrix",
            "oracle_pass_interval_rmse",
        )
        for row in self.rows:
            for col in bool_cols:
                self.assertIn(
                    row[col],
                    {"true", "false"},
                    msg=f"{col}={row[col]!r} is not boolean",
                )

    def test_d_spatial_is_d_minus_one(self) -> None:
        for row in self.rows:
            d = int(row["target_dim"])
            ds = int(row["d_spatial"])
            self.assertEqual(ds, d - 1)

    def test_pair_count_matches_n(self) -> None:
        for row in self.rows:
            n = int(row["n"])
            expected = n * (n - 1) // 2
            self.assertEqual(int(row["pair_count"]), expected)

    def test_discordant_counts_consistent(self) -> None:
        for row in self.rows:
            fp = int(row["false_positive_pairs"])
            fn = int(row["false_negative_pairs"])
            total = int(row["total_discordant_pairs"])
            self.assertEqual(fp + fn, total)


if __name__ == "__main__":
    unittest.main()
