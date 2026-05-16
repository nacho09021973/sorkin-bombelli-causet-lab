#!/usr/bin/env python3
"""Build the Phase 1C two-control finite-size scaling atlas.

Phase 1B showed ensemble-level separation between canonical
Minkowski-diamond sprinklings and Kleitman-Rothschild controls.
Phase 1C repeats the same finite-size grid with a second
non-manifoldlike control: suspended corona/crown posets. The point is
to check whether the signature "Myrheim-Meyer roughly flat, midpoint
growing with log n" is specific to the KR generator or persists on a
different cheap non-manifoldlike family.
"""

from __future__ import annotations

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
from tools.build_phase1b_scaling_atlas import (  # noqa: E402
    FOUNDATION,
    SEEDS,
    SIZES,
    SPACETIME_DIMS,
    _CSV_HEADERS,
    summarize,
    write_csv,
)
import validation_suite as vs  # noqa: E402


CONTROL_GENERATORS = (
    ("kleitman_rothschild", vs.generate_kleitman_rothschild),
    ("corona_poset", vs.generate_corona_poset),
)


def build_rows() -> list[dict]:
    """Return one row per Phase 1C grid cell."""

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
            for family, generator in CONTROL_GENERATORS:
                matrix = generator(n=n, seed=seed)
                mm, mid = _estimate_dimensions(matrix)
                rows.append({
                    "family": family,
                    "d_spacetime": "",
                    "n": n,
                    "seed": seed,
                    "mm_dim": mm,
                    "midpoint_dim": mid,
                    "abs_discrepancy": _discrepancy(mm, mid),
                })
    return rows


def write_markdown(summary: list[dict], path: Path) -> None:
    lines = [
        "# Phase 1C Two-Control Scaling Atlas",
        "",
        "Ensemble statistics of two order-theoretic dimension",
        "estimators across `n` in {32, 64, 128, 256} for canonical",
        "Minkowski-diamond sprinklings (d_spacetime in {2, 3, 4}) and",
        "two non-manifoldlike controls: Kleitman-Rothschild three-level",
        "posets and suspended corona/crown posets. Each cell aggregates",
        "five seeds: 1959, 1962, 1987, 2009, 2026.",
        "",
        "Columns:",
        "",
        "- `family`, `d`, `n`, `seeds`: the cell descriptor.",
        "- `mean_mm`, `std_mm`: Myrheim-Meyer dim, ensemble mean and",
        "  population standard deviation.",
        "- `mean_midpoint`, `std_midpoint`: Meyer's midpoint scaling",
        "  dim, same convention.",
        "- `mean |disc|`: ensemble mean of `|mm - midpoint|`.",
        "",
        "What to look for:",
        "",
        "- Manifoldlike (Minkowski): `mean_mm` should converge toward",
        "  `d_spacetime` as `n` grows, and the two estimators should",
        "  move toward agreement at larger `n`.",
        "- Non-manifoldlike controls: `mean_mm` should be comparatively",
        "  flat within a family, while `mean_midpoint` grows with",
        "  finite-size scale rather than converging to the same target.",
        "- Agreement of KR and corona trends would support the claim that",
        "  the Phase 1B signature is not KR-specific.",
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
        "Regenerate via `make regen-phase1c`. Source tool:",
        "`tools/build_phase1c_scaling_atlas.py`.",
        "",
        "For interpretation see the *Phase 1C* section of",
        "`results_note_2026.md`.",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    rows = build_rows()
    FOUNDATION.mkdir(parents=True, exist_ok=True)
    write_csv(rows, FOUNDATION / "phase1c_scaling_atlas.csv")
    summary = summarize(rows)
    write_markdown(summary, FOUNDATION / "phase1c_scaling_atlas.md")
    print(
        f"Wrote {len(rows)} atlas rows ({len(summary)} summary rows) "
        f"to {FOUNDATION}"
    )


if __name__ == "__main__":
    main()
