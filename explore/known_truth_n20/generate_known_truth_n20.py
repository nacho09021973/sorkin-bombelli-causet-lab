#!/usr/bin/env python3
"""Freeze a deterministic N=20 known-truth Minkowski sprinkling case.

This is a manufactured known-truth reference for SORKIN-2 diagnostics.
It does not use any external community golden table.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import sys
from pathlib import Path
from typing import Iterable


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import causet_invariants  # noqa: E402
import validation_suite as vs  # noqa: E402


OUT_DIR = Path(__file__).resolve().parent
CASE_METADATA_PATH = OUT_DIR / "case_metadata.json"
COORDS_PATH = OUT_DIR / "ground_truth_coords.csv"
TARGET_PAIRS_PATH = OUT_DIR / "target_pairs.csv"
TARGET_MATRIX_PATH = OUT_DIR / "target_matrix.csv"
INVARIANTS_PATH = OUT_DIR / "target_invariants.json"
CHECKSUMS_PATH = OUT_DIR / "checksums_sha256.txt"

N = 20
D_SPACETIME = 2
CASE_SEED = 1959
FAMILY = "minkowski"
CASE_ID = f"minkowski_n{N}_d{D_SPACETIME}_seed{CASE_SEED}"
GENERATOR = "validation_suite.sprinkle_minkowski_diamond"


def _pairs(matrix: list[list[bool]]) -> list[tuple[int, int]]:
    return [
        (i, j)
        for i in range(len(matrix) - 1)
        for j in range(i + 1, len(matrix))
        if matrix[i][j]
    ]


def _write_coords(points: list[tuple[float, ...]]) -> None:
    spatial_headers = [f"x{k}" for k in range(1, D_SPACETIME)]
    with COORDS_PATH.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["event_index", "t", *spatial_headers])
        for idx, point in enumerate(points):
            writer.writerow([idx, *[f"{value:.17g}" for value in point]])


def _write_pairs(pairs: Iterable[tuple[int, int]]) -> None:
    with TARGET_PAIRS_PATH.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["i", "j"])
        writer.writerows(pairs)


def _write_matrix(matrix: list[list[bool]]) -> None:
    headers = ["i", *[f"j{j}" for j in range(len(matrix))]]
    with TARGET_MATRIX_PATH.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(headers)
        for i, row in enumerate(matrix):
            writer.writerow([i, *[1 if value else 0 for value in row]])


def _json_default(value: object) -> object:
    if isinstance(value, float):
        if math.isfinite(value):
            return value
        return str(value)
    return value


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> int:
    matrix, points = vs.sprinkle_minkowski_diamond(
        n=N,
        seed=CASE_SEED,
        d_spacetime=D_SPACETIME,
    )
    pairs = _pairs(matrix)

    metadata = {
        "case_id": CASE_ID,
        "family": FAMILY,
        "n": N,
        "d_spacetime": D_SPACETIME,
        "case_seed": CASE_SEED,
        "geometry": "Minkowski 1+1 unit causal diamond",
        "generator": GENERATOR,
        "truth_status": "manufactured_known_truth_by_fixed_sprinkling_seed",
        "notes": [
            "This is not an external community golden table.",
            "The target order and coordinates are frozen from validation_suite.sprinkle_minkowski_diamond.",
        ],
    }
    CASE_METADATA_PATH.write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    _write_coords(points)
    _write_pairs(pairs)
    _write_matrix(matrix)

    fingerprint = causet_invariants.invariants_fingerprint(matrix)
    invariant_payload = {
        "case_id": CASE_ID,
        "total_relations_target": len(pairs),
        "ordering_fraction": causet_invariants.ordering_fraction(matrix),
        "longest_chain": causet_invariants.height(matrix),
        "layer_rank_summary": causet_invariants.antichain_profile(matrix),
        "myrheim_meyer_dimension": causet_invariants.myrheim_meyer_dimension(matrix),
        "chain_counts": causet_invariants.chain_counts(matrix, k_max=4),
        "link_count": causet_invariants.link_count(matrix),
        "invariants_fingerprint": fingerprint,
        "omitted_existing_invariants": [
            "interval abundance: no public reusable interval-abundance helper found",
        ],
    }
    INVARIANTS_PATH.write_text(
        json.dumps(invariant_payload, indent=2, sort_keys=True, default=_json_default) + "\n",
        encoding="utf-8",
    )

    checksum_paths = [
        CASE_METADATA_PATH,
        COORDS_PATH,
        TARGET_PAIRS_PATH,
        TARGET_MATRIX_PATH,
        INVARIANTS_PATH,
    ]
    CHECKSUMS_PATH.write_text(
        "".join(f"{_sha256(path)}  {path.name}\n" for path in checksum_paths),
        encoding="utf-8",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
