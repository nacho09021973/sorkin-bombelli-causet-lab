"""Regression tests for S4-KERR-K12-EQUATORIAL-WINDING-SECTOR-AUDIT-001."""

from __future__ import annotations

import csv
import json
import math
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
KERR_DIR = ROOT / "explore" / "sorkin4_kerr_benchmark"
PREFIX = "kerr_k12_equatorial_winding_sector_audit_001_n12_seed1959"

CSV_PATH = KERR_DIR / f"{PREFIX}.csv"
JSON_PATH = KERR_DIR / f"{PREFIX}.json"
MD_PATH = KERR_DIR / f"{PREFIX}.md"
PNG_PATH = KERR_DIR / f"{PREFIX}.png"
README_PATH = KERR_DIR / "README.md"


def _load_csv() -> list[dict[str, str]]:
    with CSV_PATH.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def _load_json() -> dict:
    return json.loads(JSON_PATH.read_text(encoding="utf-8"))


def _to_bool(value: str | bool) -> bool:
    if isinstance(value, bool):
        return value
    return value.strip().lower() == "true"


class KerrK12ArtifactTests(unittest.TestCase):
    def test_artifacts_exist(self):
        assert CSV_PATH.exists(), f"Missing {CSV_PATH}"
        assert JSON_PATH.exists(), f"Missing {JSON_PATH}"
        assert MD_PATH.exists(), f"Missing {MD_PATH}"
        assert PNG_PATH.exists(), f"Missing {PNG_PATH}"

    def test_csv_json_parse(self):
        assert len(_load_csv()) > 0
        assert "cases" in _load_json()

    def test_png_nonzero_size(self):
        assert PNG_PATH.stat().st_size > 0

    def test_summary_all_checks_pass(self):
        assert _load_json()["global_summary"]["all_checks_pass"] is True

    def test_synthetic_generated_and_recovered_positive(self):
        s = _load_json()["global_summary"]
        assert s["synthetic_winding_targets_generated"] > 0
        assert s["synthetic_winding_targets_recovered"] > 0

    def test_non_advisory_non_unresolved_pass(self):
        for c in _load_json()["cases"]:
            if not c["advisory_only"] and not c["unresolved"]:
                assert c["all_checks_pass"] is True, c["case_id"]

    def test_no_sprinkling_and_no_global_relations(self):
        for c in _load_json()["cases"]:
            assert c["no_sprinkling_pair_used"] is True, c["case_id"]
            assert c["no_global_causal_relations_decided"] is True, c["case_id"]

    def test_accepted_cases_invariants(self):
        for c in _load_json()["cases"]:
            if c["all_checks_pass"]:
                assert c["all_points_exterior"] is True, c["case_id"]
                assert c["null_condition_pass"] is True, c["case_id"]
                assert c["constants_consistency_pass"] is True, c["case_id"]

    def test_schwarzschild_b0_radial_sector_zero(self):
        for c in _load_json()["cases"]:
            if c["case_id"] == "k12_schw_b0_outgoing":
                assert abs(c["delta_phi_raw"]) <= 1.0e-12
                assert c["best_sector_m"] == 0

    def test_recovered_targets_sector_match(self):
        for c in _load_json()["cases"]:
            if (not c["advisory_only"]) and c["synthetic_winding_sector_recovered"]:
                assert c["correct_sector_recovered"] is True, c["case_id"]
                assert c["best_sector_m"] == c["correct_sector_m"], c["case_id"]
                assert abs(c["best_sector_residual"]) <= 1.0e-6, c["case_id"]

    def test_caveats_in_md_and_readme(self):
        text = (
            MD_PATH.read_text(encoding="utf-8").lower()
            + "\n"
            + README_PATH.read_text(encoding="utf-8").lower()
        )
        required = [
            "synthetic targets",
            "does not use sprinkling event pairs",
            "does not decide causal reachability",
            "does not implement a production kerr causal classifier",
            "correct_sector_recovered is not physical reachability",
        ]
        for needle in required:
            assert needle in text, needle

    def test_positive_spin_global_causal_accounting(self):
        s = _load_json()["global_summary"]
        assert s["global_true_relations"] == 0
        assert s["global_false_relations"] == 0
        assert s["global_undecided_pairs"] == 66

    def test_csv_boolean_fields(self):
        for r in _load_csv():
            assert _to_bool(r["no_sprinkling_pair_used"]), r["case_id"]
            assert _to_bool(r["no_global_causal_relations_decided"]), r["case_id"]
