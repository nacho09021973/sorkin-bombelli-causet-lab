"""Regression tests for S4-KERR-K11-EQUATORIAL-SHOOTING-SANDBOX-001."""

from __future__ import annotations

import csv
import json
import math
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
KERR_DIR = ROOT / "explore" / "sorkin4_kerr_benchmark"
PREFIX = "kerr_k11_equatorial_shooting_sandbox_001_n12_seed1959"

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


class KerrK11ArtifactTests(unittest.TestCase):
    def test_artifacts_exist(self):
        assert CSV_PATH.exists(), f"Missing {CSV_PATH}"
        assert JSON_PATH.exists(), f"Missing {JSON_PATH}"
        assert MD_PATH.exists(), f"Missing {MD_PATH}"
        assert PNG_PATH.exists(), f"Missing {PNG_PATH}"

    def test_csv_json_parseable(self):
        rows = _load_csv()
        payload = _load_json()
        assert len(rows) > 0
        assert "cases" in payload

    def test_png_nonzero_size(self):
        assert PNG_PATH.stat().st_size > 0

    def test_summary_all_checks_pass(self):
        assert _load_json()["global_summary"]["all_checks_pass"] is True

    def test_summary_target_counts(self):
        summary = _load_json()["global_summary"]
        assert summary["synthetic_targets_generated"] > 0
        assert summary["synthetic_targets_hit"] > 0

    def test_non_advisory_non_unresolved_cases_pass(self):
        for case in _load_json()["cases"]:
            if not case["advisory_only"] and not case["unresolved"]:
                assert case["all_checks_pass"] is True, case["case_id"]

    def test_target_generated_and_no_sprinkling(self):
        for case in _load_json()["cases"]:
            assert case["target_was_forward_generated"] is True, case["case_id"]
            assert case["no_sprinkling_pair_used"] is True, case["case_id"]

    def test_all_points_exterior(self):
        for case in _load_json()["cases"]:
            assert case["all_points_exterior"] is True, case["case_id"]

    def test_non_advisory_accepted_null_and_constants(self):
        for case in _load_json()["cases"]:
            if case["all_checks_pass"] and not case["advisory_only"]:
                assert case["null_condition_pass"] is True, case["case_id"]
                assert case["constants_consistency_pass"] is True, case["case_id"]

    def test_schwarzschild_radial_known_answer_gate(self):
        for case in _load_json()["cases"]:
            if case["case_type"] == "schwarzschild_radial_known_answer_shooting":
                assert case["schwarzschild_radial_limit_pass"] is True, case["case_id"]

    def test_bshoot_converged_accuracy(self):
        for case in _load_json()["cases"]:
            if "b_shooting" in case["case_type"] and case["solver_converged"]:
                assert case["recovered_b_error"] is not None, case["case_id"]
                assert case["recovered_b_error"] <= 1.0e-5, case["case_id"]
                assert abs(case["endpoint_phi_residual"]) <= 1.0e-6, case["case_id"]

    def test_synthetic_hit_weighted_residual(self):
        for case in _load_json()["cases"]:
            if case["synthetic_target_hit"]:
                assert case["endpoint_weighted_residual"] <= 1.0, case["case_id"]

    def test_caveats_in_md_and_readme(self):
        text = (
            MD_PATH.read_text(encoding="utf-8").lower()
            + "\n"
            + README_PATH.read_text(encoding="utf-8").lower()
        )
        required = [
            "synthetic known-answer",
            "does not use sprinkling event pairs",
            "does not decide causal reachability",
            "does not implement a production kerr causal classifier",
        ]
        for needle in required:
            assert needle in text, needle

    def test_positive_spin_global_causal_accounting(self):
        summary = _load_json()["global_summary"]
        assert summary["global_true_relations"] == 0
        assert summary["global_false_relations"] == 0
        assert summary["global_undecided_pairs"] == 66

    def test_csv_flags(self):
        for row in _load_csv():
            assert _to_bool(row["target_was_forward_generated"]), row["case_id"]
            assert _to_bool(row["no_sprinkling_pair_used"]), row["case_id"]
            assert _to_bool(row["all_points_exterior"]), row["case_id"]
