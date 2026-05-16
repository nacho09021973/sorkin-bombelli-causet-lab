#!/usr/bin/env python3
"""Build the Phase 1B finite-size scaling atlas.

Extends the Phase 1 order-theoretic atlas across a range of sizes
``n in {32, 64, 128, 256}`` for the same families
(canonical Minkowski sprinklings at ``d_spacetime in {2, 3, 4}`` and
Kleitman-Rothschild non-manifoldlike controls), with the same fixed
seeds. The physical question is *finite-size scaling*: does the
separation between manifoldlike sprinklings and non-manifoldlike
controls improve systematically as ``n`` grows?

This phase does **not** invoke any embedding algorithm, energy
function, or optimizer. The single artifact it produces is a
reproducible table at ``benchmarks/foundation/phase1b_scaling_atlas.{csv,md}``.

The Phase 1 helpers ``_estimate_dimensions``, ``_discrepancy``, and
``_format_field`` are reused so that the two atlases share their
numerical conventions verbatim.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from tools.build_phase1_atlas import (  # noqa: E402
    _discrepancy,
    _estimate_dimensions,
    _format_field,
)
import validation_suite as vs  # noqa: E402


FOUNDATION = ROOT / "benchmarks" / "foundation"


# Phase 1B grid. Same five seeds as Phase 1 so each Phase 1 cell at
# ``n in {32, 64}`` appears as a subset of the Phase 1B grid; no
# reduction is applied at ``n = 256`` since per-cell cost stays
# under a few seconds on stdlib Python.
SIZES = (32, 64, 128, 256)
SEEDS = (1959, 1962, 1987, 2009, 2026)
SPACETIME_DIMS = (2, 3, 4)


def build_rows() -> list[dict]:
    """Return one row per ``(family, target_dim, n, seed)`` cell.

    The Minkowski rows record ``d_spacetime``; the KR rows leave
    it blank because the construction is intrinsically
    non-embeddable.
    """

    rows: list[dict] = []
    for n in SIZES:
        for seed in SEEDS:
            for d in SPACETIME_DIMS:
                matrix, _ = vs.sprinkle_minkowski_diamond(
                    n=n, seed=seed, d_spacetime=d
                )
                mm, mid = _estimate_dimensions(matrix)
                rows.append({
                    "family": "minkowski",
                    "d_spacetime": d,
                    "n": n,
                    "seed": seed,
                    "mm_dim": mm,
                    "midpoint_dim": mid,
                    "abs_discrepancy": _discrepancy(mm, mid),
                })
            kr_matrix = vs.generate_kleitman_rothschild(n=n, seed=seed)
            mm, mid = _estimate_dimensions(kr_matrix)
            rows.append({
                "family": "kleitman_rothschild",
                "d_spacetime": "",
                "n": n,
                "seed": seed,
                "mm_dim": mm,
                "midpoint_dim": mid,
                "abs_discrepancy": _discrepancy(mm, mid),
            })
    return rows


_CSV_HEADERS = (
    "family",
    "d_spacetime",
    "n",
    "seed",
    "mm_dim",
    "midpoint_dim",
    "abs_discrepancy",
)


def write_csv(rows: list[dict], path: Path) -> None:
    lines = [",".join(_CSV_HEADERS)]
    for row in rows:
        lines.append(",".join(_format_field(row[h]) for h in _CSV_HEADERS))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _ensemble_mean(values) -> float:
    finite = [v for v in values if isinstance(v, float) and math.isfinite(v)]
    if not finite:
        return float("nan")
    return sum(finite) / len(finite)


def _ensemble_std(values) -> float:
    """Population standard deviation over the finite entries."""

    finite = [v for v in values if isinstance(v, float) and math.isfinite(v)]
    if len(finite) < 2:
        return float("nan")
    mu = sum(finite) / len(finite)
    variance = sum((v - mu) ** 2 for v in finite) / len(finite)
    return math.sqrt(variance)


def summarize(rows: list[dict]) -> list[dict]:
    """Aggregate by ``(family, d_spacetime, n)`` over the seed axis."""

    buckets: dict = {}
    for row in rows:
        key = (row["family"], row["d_spacetime"], row["n"])
        buckets.setdefault(key, []).append(row)

    out: list[dict] = []
    for key in sorted(
        buckets,
        key=lambda kv: (kv[0], str(kv[1]), kv[2]),
    ):
        rs = buckets[key]
        family, d, n = key
        mm_vals = [r["mm_dim"] for r in rs]
        mid_vals = [r["midpoint_dim"] for r in rs]
        disc_vals = [r["abs_discrepancy"] for r in rs]
        out.append({
            "family": family,
            "d_spacetime": d,
            "n": n,
            "seeds": len(rs),
            "mean_mm": _ensemble_mean(mm_vals),
            "std_mm": _ensemble_std(mm_vals),
            "mean_midpoint": _ensemble_mean(mid_vals),
            "std_midpoint": _ensemble_std(mid_vals),
            "mean_abs_discrepancy": _ensemble_mean(disc_vals),
        })
    return out


def write_markdown(summary: list[dict], path: Path) -> None:
    lines = [
        "# Phase 1B Finite-Size Scaling Atlas",
        "",
        "Ensemble statistics of two order-theoretic dimension",
        "estimators across `n` in {32, 64, 128, 256} for canonical",
        "Minkowski-diamond sprinklings (d_spacetime in {2, 3, 4}) and",
        "Kleitman-Rothschild three-level non-manifoldlike controls.",
        "Each cell aggregates five seeds: 1959, 1962, 1987, 2009, 2026.",
        "",
        "Columns:",
        "",
        "- `family`, `d`, `n`, `seeds`: the cell descriptor.",
        "- `mean_mm`, `std_mm`: Myrheim-Meyer dim, ensemble mean and",
        "  population standard deviation (zero in the small-ensemble",
        "  case the row has fewer than two finite values).",
        "- `mean_midpoint`, `std_midpoint`: Meyer's midpoint scaling",
        "  dim, same convention.",
        "- `mean |disc|`: ensemble mean of `|mm - midpoint|`.",
        "",
        "What to look for:",
        "",
        "- Manifoldlike (Minkowski): `mean_mm` should converge toward",
        "  `d_spacetime` as `n` grows; `std_mm` should shrink.",
        "- Non-manifoldlike (KR): `mean_mm` should remain approximately",
        "  flat in `n`.",
        "- If `mean_midpoint` also converges toward the true `d` for",
        "  Minkowski while staying disjoint from KR, the two",
        "  estimators jointly support per-causet classification. If",
        "  not, the diagnostic is ensemble-level only.",
        "",
        "| family | d | n | seeds | mean MM | std MM | mean midpoint | std midpoint | mean \\|disc\\| |",
        "| --- | :---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in summary:
        d_display = row["d_spacetime"] if row["d_spacetime"] != "" else "-"
        lines.append(
            "| {family} | {d} | {n} | {seeds} | "
            "{mm} | {sm} | {mid} | {smid} | {disc} |".format(
                family=row["family"],
                d=d_display,
                n=row["n"],
                seeds=row["seeds"],
                mm=_format_field(row["mean_mm"]),
                sm=_format_field(row["std_mm"]),
                mid=_format_field(row["mean_midpoint"]),
                smid=_format_field(row["std_midpoint"]),
                disc=_format_field(row["mean_abs_discrepancy"]),
            )
        )
    lines += [
        "",
        "Regenerate via `make regen-phase1b`. Source tool:",
        "`tools/build_phase1b_scaling_atlas.py`.",
        "",
        "For the conservative interpretation across these numbers see",
        "the *Phase 1B* section of `results_note_2026.md`.",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    rows = build_rows()
    FOUNDATION.mkdir(parents=True, exist_ok=True)
    write_csv(rows, FOUNDATION / "phase1b_scaling_atlas.csv")
    summary = summarize(rows)
    write_markdown(summary, FOUNDATION / "phase1b_scaling_atlas.md")
    print(
        f"Wrote {len(rows)} atlas rows ({len(summary)} summary rows) "
        f"to {FOUNDATION}"
    )


if __name__ == "__main__":
    main()
