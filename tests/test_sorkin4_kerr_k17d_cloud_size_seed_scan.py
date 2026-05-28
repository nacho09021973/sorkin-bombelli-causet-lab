"""Regression tests for S4-KERR-K17D-CLOUD-SIZE-SEED-SCAN-001."""

from __future__ import annotations

import csv
import json
import subprocess
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
KERR_DIR = ROOT / "explore" / "sorkin4_kerr_benchmark"
SCRIPT_PATH = KERR_DIR / "audit_kerr_k17d_cloud_size_seed_scan_001.py"
PREFIX = "kerr_k17d_cloud_size_seed_scan_001"

CSV_PATH = KERR_DIR / f"{PREFIX}.csv"
JSON_PATH = KERR_DIR / f"{PREFIX}.json"
MD_PATH = KERR_DIR / f"{PREFIX}.md"
PNG_PATH = KERR_DIR / f"{PREFIX}.png"


def _load_csv() -> list[dict[str, str]]:
    with CSV_PATH.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def _load_json() -> dict:
    return json.loads(JSON_PATH.read_text(encoding="utf-8"))


def _to_bool(v: str | bool) -> bool:
    if isinstance(v, bool):
        return v
    return v.strip().lower() == "true"


class KerrK17dCloudSizeSeedScanTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        subprocess.run(
            ["python3", str(SCRIPT_PATH), "--full"],
            check=True,
            cwd=str(ROOT),
        )

    def test_artifacts_exist(self):
        assert CSV_PATH.exists()
        assert JSON_PATH.exists()
        assert MD_PATH.exists()
        assert PNG_PATH.exists()

    def test_csv_json_parse(self):
        assert len(_load_csv()) > 0
        payload = _load_json()
        assert "selected_pairs" in payload
        assert "cell_summaries" in payload
        assert "global_summary" in payload
        assert "caveats" in payload

    def test_png_nonzero(self):
        assert PNG_PATH.stat().st_size > 0

    def test_cell_summaries_27(self):
        assert len(_load_json()["cell_summaries"]) == 27

    def test_global_summary_keys(self):
        gs = _load_json()["global_summary"]
        for key in (
            "median_best_residual_per_N",
            "monotone_N_decrease_flag",
            "N_at_which_W_TOL_first_crossed",
            "N_at_which_10x_W_TOL_first_crossed",
            "recommendation",
        ):
            assert key in gs, key

    def test_cell_summary_rejection_counts(self):
        for cs in _load_json()["cell_summaries"]:
            for k in (
                "n_rejected_time_short",
                "n_rejected_radial_proxy",
                "n_rejected_angular_large",
                "n_selection_unresolved",
            ):
                assert k in cs, (cs["N"], cs["seed"], cs["spin_a"], k)

    def test_no_causal_claim_flags(self):
        for r in _load_csv():
            assert _to_bool(r["no_causal_claim_introduced"]), r["case_id"]
            assert _to_bool(r["no_production_classifier_introduced"]), r["case_id"]
            assert _to_bool(r["no_global_causal_relations_decided"]), r["case_id"]

    def test_md_caveats(self):
        text = MD_PATH.read_text(encoding="utf-8").lower()
        assert "selection_candidate is not reachability" in text
        assert "rejected_by_selection is not proof of spacelike separation" in text
        assert "residual_probe_pass is not causal reachability" in text

    def test_summary_all_checks_pass(self):
        assert _load_json()["global_summary"]["all_checks_pass"] is True


if __name__ == "__main__":
    unittest.main()
