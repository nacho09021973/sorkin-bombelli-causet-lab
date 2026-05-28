"""Regression tests for S4-KERR-K17-CONTROLLED-CANDIDATE-PAIR-SANDBOX-001."""

from __future__ import annotations

import csv
import json
import subprocess
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
KERR_DIR = ROOT / "explore" / "sorkin4_kerr_benchmark"
SCRIPT_PATH = KERR_DIR / "audit_kerr_k17_controlled_candidate_pair_sandbox_001.py"
PREFIX = "kerr_k17_controlled_candidate_pair_sandbox_001_n12_seed1959"

CSV_PATH = KERR_DIR / f"{PREFIX}.csv"
JSON_PATH = KERR_DIR / f"{PREFIX}.json"
MD_PATH = KERR_DIR / f"{PREFIX}.md"
PNG_PATH = KERR_DIR / f"{PREFIX}.png"

W_TOL = 1.0e-3


def _load_csv() -> list[dict[str, str]]:
    with CSV_PATH.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def _load_json() -> dict:
    return json.loads(JSON_PATH.read_text(encoding="utf-8"))


def _to_bool(v: str | bool) -> bool:
    if isinstance(v, bool):
        return v
    return v.strip().lower() == "true"


class KerrK17ControlledCandidatePairSandboxTests(unittest.TestCase):
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

    def test_basic_summary(self):
        s = _load_json()["global_summary"]
        assert s["total_candidate_pairs"] > 0
        assert s["cloud_events_available"] > 0
        assert s["negative_controls_passed"] > 0
        assert s["candidate_hits"] + s["candidate_misses"] + s["candidate_undecided"] == s["total_candidate_pairs"]

    def test_non_negative_cloud_pair_flags(self):
        for c in _load_json()["cases"]:
            if c["pair_type"] != "negative_control":
                assert c["event_A_from_cloud"] is True, c["case_id"]
                assert c["event_B_from_cloud"] is True, c["case_id"]
                assert c["controlled_candidate_pair"] is True, c["case_id"]
                assert c["no_arbitrary_full_cloud_classification"] is True, c["case_id"]

    def test_global_guardrail_flags(self):
        for c in _load_json()["cases"]:
            assert c["no_sprinkling_reachability_claimed"] is True, c["case_id"]
            assert c["no_global_causal_relations_decided"] is True, c["case_id"]
            assert c["no_production_classifier_introduced"] is True, c["case_id"]

    def test_candidate_hit_invariants(self):
        for c in _load_json()["cases"]:
            if c["candidate_hit"]:
                assert c["endpoint_weighted_residual"] <= W_TOL, c["case_id"]
                assert c["all_points_exterior"] is True, c["case_id"]
                assert c["null_condition_pass"] is True, c["case_id"]
                assert c["constants_consistency_pass"] is True, c["case_id"]
                assert c["all_checks_pass"] is True, c["case_id"]

    def test_candidate_miss_wording(self):
        # enforce label-only contract via artifact md
        text = MD_PATH.read_text(encoding="utf-8").lower()
        assert "candidate_miss is not proof of spacelike separation" in text
        assert "causal_false" not in text

    def test_md_caveats(self):
        text = MD_PATH.read_text(encoding="utf-8").lower()
        required = [
            "controlled candidate pairs",
            "does not classify the full cloud",
            "candidate_hit is a sandbox numerical recovery, not physical/global causal reachability",
            "candidate_undecided is the default conservative result",
        ]
        for needle in required:
            assert needle in text, needle

    def test_global_accounting(self):
        s = _load_json()["global_summary"]
        assert s["global_true_relations"] == 0
        assert s["global_false_relations"] == 0
        assert s["global_undecided_pairs"] == 66

    def test_csv_guardrails(self):
        for r in _load_csv():
            assert _to_bool(r["no_sprinkling_reachability_claimed"]), r["case_id"]
            assert _to_bool(r["no_global_causal_relations_decided"]), r["case_id"]
            assert _to_bool(r["no_production_classifier_introduced"]), r["case_id"]


if __name__ == "__main__":
    unittest.main()
