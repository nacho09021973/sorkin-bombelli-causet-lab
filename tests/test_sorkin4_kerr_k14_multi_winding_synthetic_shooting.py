"""Regression tests for S4-KERR-K14-MULTI-WINDING-SYNTHETIC-SHOOTING-001."""

from __future__ import annotations

import csv
import json
import math
import subprocess
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
KERR_DIR = ROOT / "explore" / "sorkin4_kerr_benchmark"
SCRIPT_PATH = KERR_DIR / "audit_kerr_k14_multi_winding_synthetic_shooting_001.py"
PREFIX = "kerr_k14_multi_winding_synthetic_shooting_001_n12_seed1959"
K13B_JSON = KERR_DIR / "kerr_k13b_near_photon_whirling_probe_001_n12_seed1959.json"

CSV_PATH = KERR_DIR / f"{PREFIX}.csv"
JSON_PATH = KERR_DIR / f"{PREFIX}.json"
MD_PATH = KERR_DIR / f"{PREFIX}.md"
PNG_PATH = KERR_DIR / f"{PREFIX}.png"

PHI_TOL = 1.0e-5
WEIGHTED_TOL = 1.0
LAMBDA_TOL = 1.0e-5


def _load_csv() -> list[dict[str, str]]:
    with CSV_PATH.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def _load_json() -> dict:
    return json.loads(JSON_PATH.read_text(encoding="utf-8"))


def _to_bool(v: str | bool) -> bool:
    if isinstance(v, bool):
        return v
    return v.strip().lower() == "true"


class KerrK14MultiWindingSyntheticShootingTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run(["python3", str(SCRIPT_PATH)], check=True, cwd=str(ROOT))

    def test_artifacts_exist(self):
        assert CSV_PATH.exists()
        assert JSON_PATH.exists()
        assert MD_PATH.exists()
        assert PNG_PATH.exists()

    def test_csv_json_parse(self):
        assert len(_load_csv()) > 0
        assert "cases" in _load_json()

    def test_png_nonzero(self):
        assert PNG_PATH.stat().st_size > 0

    def test_summary_all_checks_pass(self):
        assert _load_json()["global_summary"]["all_checks_pass"] is True

    def test_whirling_selected_and_recovered(self):
        s = _load_json()["global_summary"]
        assert s["whirling_targets_selected"] > 0
        assert s["whirling_targets_recovered"] > 0

    def test_max_abs_delta_phi_target_gt_pi(self):
        assert _load_json()["global_summary"]["max_abs_delta_phi_target"] > math.pi

    def test_conditional_pi_and_2pi_recovered(self):
        k13b = json.loads(K13B_JSON.read_text(encoding="utf-8"))
        accepted = [
            c for c in k13b["cases"]
            if (not c["advisory_only"]) and (not c["unresolved"]) and c["all_checks_pass"]
        ]
        has_pi = any(c["abs_delta_phi"] > math.pi for c in accepted)
        has_2pi = any(c["abs_delta_phi"] > 2.0 * math.pi for c in accepted)
        s = _load_json()["global_summary"]
        if has_pi:
            assert s["any_pi_whirling_target_recovered"] is True
        if has_2pi and any(c["abs_delta_phi_target"] > 2.0 * math.pi for c in _load_json()["cases"] if c["source_was_k13b_whirling"]):
            assert s["any_2pi_whirling_target_recovered"] is True

    def test_non_advisory_non_unresolved_pass(self):
        for c in _load_json()["cases"]:
            if not c["advisory_only"] and not c["unresolved"]:
                assert c["all_checks_pass"] is True, c["case_id"]

    def test_target_forward_generated_and_guardrails(self):
        for c in _load_json()["cases"]:
            assert c["target_was_forward_generated"] is True, c["case_id"]
            assert c["no_sprinkling_pair_used"] is True, c["case_id"]
            assert c["no_global_causal_relations_decided"] is True, c["case_id"]
            assert c["no_causal_classifier_introduced"] is True, c["case_id"]

    def test_accepted_invariants(self):
        for c in _load_json()["cases"]:
            if c["all_checks_pass"]:
                assert c["all_points_exterior"] is True, c["case_id"]
                assert c["null_condition_pass"] is True, c["case_id"]
                assert c["constants_consistency_pass"] is True, c["case_id"]

    def test_recovered_whirling_tolerances_and_sector(self):
        for c in _load_json()["cases"]:
            if c["physical_whirling_synthetic_target_recovered"]:
                assert abs(c["endpoint_phi_residual_sector_adjusted"]) <= PHI_TOL, c["case_id"]
                assert c["endpoint_weighted_residual"] <= WEIGHTED_TOL, c["case_id"]
                assert c["correct_sector_recovered"] is True, c["case_id"]

    def test_lambda_mode_error_tolerance(self):
        for c in _load_json()["cases"]:
            if c["recovery_mode"] == "lambda_fixed_b" and c["synthetic_known_answer_recovered"]:
                assert c["recovered_lambda_error"] is not None, c["case_id"]
                assert c["recovered_lambda_error"] <= LAMBDA_TOL, c["case_id"]

    def test_md_caveats(self):
        text = MD_PATH.read_text(encoding="utf-8").lower()
        required = [
            "synthetic known-answer",
            "not causal reachability",
            "does not use sprinkling event pairs",
            "does not implement a production kerr causal classifier",
        ]
        for needle in required:
            assert needle in text, needle

    def test_global_causal_accounting(self):
        s = _load_json()["global_summary"]
        assert s["global_true_relations"] == 0
        assert s["global_false_relations"] == 0
        assert s["global_undecided_pairs"] == 66

    def test_csv_guardrail_flags(self):
        for r in _load_csv():
            assert _to_bool(r["no_sprinkling_pair_used"]), r["case_id"]
            assert _to_bool(r["no_global_causal_relations_decided"]), r["case_id"]
            assert _to_bool(r["no_causal_classifier_introduced"]), r["case_id"]


if __name__ == "__main__":
    unittest.main()
