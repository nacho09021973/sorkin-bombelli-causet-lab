#!/usr/bin/env python3
"""Build the Phase 1D structural atlas.

Phase 1D keeps the Phase 1C families and finite-size grid, but adds a
third order-theoretic diagnostic: 3-chain abundance. The observable is
reported as both raw chain counts and the normalized fraction of
triples that form a 3-chain. It is intentionally not converted into a
dimension estimator; the aim is to test whether one more independent
order statistic improves ensemble-level separation without invoking
coordinates, energy, embedding, fitting, or optimization.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import causet_invariants as ci  # noqa: E402
from tools.build_phase1_atlas import (  # noqa: E402
    _discrepancy,
    _estimate_dimensions,
    _format_field,
)
from tools.build_phase1b_scaling_atlas import (  # noqa: E402
    FOUNDATION,
    SEEDS,
    SIZES,
    SPACETIME_DIMS,
)
from tools.build_phase1c_scaling_atlas import CONTROL_GENERATORS  # noqa: E402
import validation_suite as vs  # noqa: E402


CSV_HEADERS = (
    "family",
    "target_dim",
    "n",
    "seed",
    "mm_dim",
    "midpoint_dim",
    "abs_discrepancy_mm_midpoint",
    "chain2_count",
    "chain3_count",
    "chain3_abundance",
)


def _row_for_matrix(
    *,
    family: str,
    target_dim: int | str,
    n: int,
    seed: int,
    matrix,
) -> dict:
    mm, mid = _estimate_dimensions(matrix)
    counts = ci.chain_counts(matrix, k_max=3)
    return {
        "family": family,
        "target_dim": target_dim,
        "n": n,
        "seed": seed,
        "mm_dim": mm,
        "midpoint_dim": mid,
        "abs_discrepancy_mm_midpoint": _discrepancy(mm, mid),
        "chain2_count": counts[2],
        "chain3_count": counts[3],
        "chain3_abundance": ci.three_chain_abundance(matrix),
    }


def build_rows() -> list[dict]:
    """Return one row per Phase 1D grid cell."""

    rows: list[dict] = []
    for n in SIZES:
        for seed in SEEDS:
            for d in SPACETIME_DIMS:
                matrix, _ = vs.sprinkle_minkowski_diamond(
                    n=n, seed=seed, d_spacetime=d
                )
                rows.append(
                    _row_for_matrix(
                        family="minkowski",
                        target_dim=d,
                        n=n,
                        seed=seed,
                        matrix=matrix,
                    )
                )
            for family, generator in CONTROL_GENERATORS:
                rows.append(
                    _row_for_matrix(
                        family=family,
                        target_dim="",
                        n=n,
                        seed=seed,
                        matrix=generator(n=n, seed=seed),
                    )
                )
    return rows


def write_csv(rows: list[dict], path: Path) -> None:
    lines = [",".join(CSV_HEADERS)]
    for row in rows:
        lines.append(",".join(_format_field(row[h]) for h in CSV_HEADERS))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _ensemble_mean(values) -> float:
    finite = [v for v in values if isinstance(v, (int, float)) and math.isfinite(v)]
    if not finite:
        return float("nan")
    return sum(finite) / len(finite)


def _ensemble_std(values) -> float:
    finite = [v for v in values if isinstance(v, (int, float)) and math.isfinite(v)]
    if len(finite) < 2:
        return float("nan")
    mu = sum(finite) / len(finite)
    variance = sum((v - mu) ** 2 for v in finite) / len(finite)
    return math.sqrt(variance)


def summarize(rows: list[dict]) -> list[dict]:
    """Aggregate by ``(family, target_dim, n)`` over seeds."""

    buckets: dict = {}
    for row in rows:
        key = (row["family"], row["target_dim"], row["n"])
        buckets.setdefault(key, []).append(row)

    out: list[dict] = []
    for key in sorted(buckets, key=lambda kv: (kv[0], str(kv[1]), kv[2])):
        rs = buckets[key]
        family, target_dim, n = key
        out.append({
            "family": family,
            "target_dim": target_dim,
            "n": n,
            "seeds": len(rs),
            "mean_mm": _ensemble_mean([r["mm_dim"] for r in rs]),
            "std_mm": _ensemble_std([r["mm_dim"] for r in rs]),
            "mean_midpoint": _ensemble_mean([r["midpoint_dim"] for r in rs]),
            "std_midpoint": _ensemble_std([r["midpoint_dim"] for r in rs]),
            "mean_abs_discrepancy": _ensemble_mean(
                [r["abs_discrepancy_mm_midpoint"] for r in rs]
            ),
            "mean_chain3_abundance": _ensemble_mean(
                [r["chain3_abundance"] for r in rs]
            ),
            "std_chain3_abundance": _ensemble_std(
                [r["chain3_abundance"] for r in rs]
            ),
        })
    return out


def write_markdown(summary: list[dict], path: Path) -> None:
    lines = [
        "# Phase 1D Structural Atlas",
        "",
        "Ensemble statistics for the Phase 1C finite-size grid with",
        "one additional order-theoretic observable: 3-chain abundance.",
        "The raw CSV records `chain2_count`, `chain3_count`, and",
        "`chain3_abundance = chain3_count / binom(n, 3)`.",
        "",
        "This is not a calibrated dimension estimator. It is a minimal",
        "structural statistic asking whether higher-order chain",
        "abundance helps separate Minkowski sprinklings from",
        "Kleitman-Rothschild and suspended corona controls without",
        "coordinates, energy, embedding, fitting, or optimization.",
        "",
        "| family | target d | n | seeds | mean MM | std MM | mean midpoint | std midpoint | mean \\|disc\\| | mean C3 abundance | std C3 abundance |",
        "| --- | :---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in summary:
        d_display = row["target_dim"] if row["target_dim"] != "" else "-"
        lines.append(
            "| {family} | {d} | {n} | {seeds} | "
            "{mm} | {smm} | {mid} | {smid} | {disc} | {c3} | {sc3} |".format(
                family=row["family"],
                d=d_display,
                n=row["n"],
                seeds=row["seeds"],
                mm=_format_field(row["mean_mm"]),
                smm=_format_field(row["std_mm"]),
                mid=_format_field(row["mean_midpoint"]),
                smid=_format_field(row["std_midpoint"]),
                disc=_format_field(row["mean_abs_discrepancy"]),
                c3=_format_field(row["mean_chain3_abundance"]),
                sc3=_format_field(row["std_chain3_abundance"]),
            )
        )
    lines += [
        "",
        "Reading by family:",
        "",
        "- Minkowski d=2: high relation density and high C3 abundance;",
        "  the observable is expected to overlap corona controls at",
        "  larger n because both sit near effective dimension two by",
        "  ordering fraction.",
        "- Minkowski d=3,d=4: C3 abundance falls rapidly with dimension",
        "  and gives an independent structural scale beside MM and",
        "  midpoint.",
        "- Kleitman-Rothschild: intermediate C3 abundance with growing",
        "  MM/midpoint discrepancy.",
        "- Corona: high C3 abundance but a sharply growing midpoint",
        "  discrepancy, so C3 alone is not a manifoldness classifier.",
        "",
        "Regenerate via `make regen-phase1d`. Source tool:",
        "`tools/build_phase1d_structural_atlas.py`.",
        "",
        "For interpretation see the *Phase 1D* section of",
        "`results_note_2026.md`.",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    rows = build_rows()
    FOUNDATION.mkdir(parents=True, exist_ok=True)
    write_csv(rows, FOUNDATION / "phase1d_structural_atlas.csv")
    summary = summarize(rows)
    write_markdown(summary, FOUNDATION / "phase1d_structural_atlas.md")
    print(
        f"Wrote {len(rows)} structural rows ({len(summary)} summary rows) "
        f"to {FOUNDATION}"
    )


if __name__ == "__main__":
    main()
