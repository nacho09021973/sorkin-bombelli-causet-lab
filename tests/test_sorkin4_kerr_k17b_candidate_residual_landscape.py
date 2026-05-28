"""Regression tests for S4-KERR-K17B-CANDIDATE-RESIDUAL-LANDSCAPE-001."""

from __future__ import annotations

import csv
import json
import subprocess
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
KERR_DIR = ROOT / "explore" / "sorkin4_kerr_benchmark"
SCRIPT_PATH = KERR_DIR / "audit_kerr_k17b_candidate_residual_landscape_001.py"
PREFIX = "kerr_k17b_candidate_residual_landscape_001_n12_seed1959"
K17_JSON = KERR_DIR / "kerr_k17_controlled_candidate_pair_sandbox_001_n12_seed1959.json"

CSV_PATH = KERR_DIR / f"{PREFIX}.csv"
JSON_PATH = KERR_DIR / f"{PREFIX}.json"
MD_PATH = KERR_DIR / f"{PREFIX}.md"
PNG_PATH = KERR_DIR / f"{PREFIX}.png"


def _load_csv() -> list[dict[str, str]]:
    with CSV_PATH.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def _load_json() -> dict:
    return json.loads(JSON_PATH.read_text(encoding="utf-8"))


class KerrK17bCandidateResidualLandscapeTests(unittest.TestCase):
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

    def test_summary_pass(self):
        assert _load_json()["global_summary"]["all_checks_pass"] is True

    def test_pair_count_matches_k17(self):
        k17_total = json.loads(K17_JSON.read_text(encoding="utf-8"))["global_summary"]["total_candidate_pairs"]
        assert _load_json()["global_summary"]["total_pairs_analyzed"] == k17_total

    def test_k17_original_zero_hits(self):
        s = _load_json()["global_summary"]
        assert s["original_candidate_hits"] == 0
        assert s["original_candidate_undecided"] > 0

    def test_best_residual_not_worse(self):
        s = _load_json()["global_summary"]
        assert s["best_residual_expanded"] <= s["best_residual_original"]

    def test_no_causal_claim_flags(self):
        for c in _load_json()["cases"]:
            assert c["no_causal_claim_introduced"] is True
            assert c["no_production_classifier_introduced"] is True

    def test_md_near_hit_caveat(self):
        text = MD_PATH.read_text(encoding="utf-8").lower()
        assert "diagnostic_near_hit is not reachability" in text


if __name__ == "__main__":
    unittest.main()
