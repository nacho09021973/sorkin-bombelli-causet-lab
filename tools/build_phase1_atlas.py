#!/usr/bin/env python3
"""Build the Phase 1 order-theoretic atlas.

For each ``(n, seed)`` cell in the foundation grid this tool computes
two independent order-theoretic dimension estimators on both a
canonical Minkowski-diamond sprinkling and a Kleitman-Rothschild
non-manifoldlike control of matching size and seed:

- :func:`causet_invariants.myrheim_meyer_dimension` (inversion of the
  closed-form ordering-fraction formula).
- :func:`causet_invariants.midpoint_scaling_dimension` (Meyer's
  log-ratio of interval cardinalities).

The output is a CSV and a markdown table at
``benchmarks/foundation/phase1_atlas.{csv,md}``. The CSV is the
authoritative artifact; the markdown is a human-readable summary.

The point of the atlas is **diagnostic**, not optimization. No
embedding, energy, or annealing is invoked. The question the atlas
answers is whether order-theoretic invariants alone can separate
manifoldlike causets (sprinklings) from non-manifoldlike controls
(KR posets) of the same size.
"""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import causet_invariants  # noqa: E402
import validation_suite as vs  # noqa: E402


FOUNDATION = ROOT / "benchmarks" / "foundation"


# We deliberately reuse the foundation grid sizes and seeds. The KR
# control does not have a meaningful ``d_spacetime``, so the grid is
# (n, seed) for each family. We match every Minkowski sprinkling at
# every ``d_spacetime`` against KR controls at the same ``(n, seed)``;
# the KR cells are computed once per ``(n, seed)``.
SIZES = (16, 32, 64)
SEEDS = (1959, 1962, 1987, 2009, 2026)
SPACETIME_DIMS = (2, 3, 4)


def _estimate_dimensions(matrix) -> tuple[float, float]:
    mm = causet_invariants.myrheim_meyer_dimension(matrix)
    mid = causet_invariants.midpoint_scaling_dimension(matrix)
    return mm, mid


def _discrepancy(mm: float, mid: float) -> float:
    """Absolute difference between the two estimators.

    Returns ``float('nan')`` if either input is ``nan`` or infinite.
    """

    import math
    if not math.isfinite(mm) or not math.isfinite(mid):
        return float("nan")
    return abs(mm - mid)


def build_rows() -> list[dict]:
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
                    "discrepancy": _discrepancy(mm, mid),
                })
            kr_matrix = vs.generate_kleitman_rothschild(n=n, seed=seed)
            mm, mid = _estimate_dimensions(kr_matrix)
            rows.append({
                "family": "kleitman_rothschild",
                "d_spacetime": "",  # not defined for non-manifoldlike
                "n": n,
                "seed": seed,
                "mm_dim": mm,
                "midpoint_dim": mid,
                "discrepancy": _discrepancy(mm, mid),
            })
    return rows


def _format_field(value) -> str:
    import math
    if isinstance(value, float):
        if math.isnan(value):
            return "nan"
        if math.isinf(value):
            return "inf"
        return f"{value:.4f}"
    return str(value)


def write_csv(rows: list[dict], path: Path) -> None:
    headers = (
        "family",
        "d_spacetime",
        "n",
        "seed",
        "mm_dim",
        "midpoint_dim",
        "discrepancy",
    )
    lines = [",".join(headers)]
    for row in rows:
        lines.append(",".join(_format_field(row[h]) for h in headers))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def summarize(rows: list[dict]) -> list[dict]:
    """Aggregate by ``(family, d_spacetime if applicable, n)`` over seeds."""

    import math
    buckets: dict = {}
    for row in rows:
        key = (row["family"], row["d_spacetime"], row["n"])
        buckets.setdefault(key, []).append(row)
    summary: list[dict] = []
    for (family, d, n), rs in sorted(
        buckets.items(),
        key=lambda kv: (kv[0][0], str(kv[0][1]), kv[0][2]),
    ):
        mm_values = [r["mm_dim"] for r in rs if math.isfinite(r["mm_dim"])]
        mid_values = [r["midpoint_dim"] for r in rs if math.isfinite(r["midpoint_dim"])]
        disc_values = [r["discrepancy"] for r in rs if math.isfinite(r["discrepancy"])]
        summary.append({
            "family": family,
            "d_spacetime": d,
            "n": n,
            "seeds": len(rs),
            "mean_mm": sum(mm_values) / len(mm_values) if mm_values else float("nan"),
            "mean_midpoint": sum(mid_values) / len(mid_values) if mid_values else float("nan"),
            "mean_discrepancy": sum(disc_values) / len(disc_values) if disc_values else float("nan"),
        })
    return summary


def write_markdown(summary: list[dict], path: Path) -> None:
    lines = [
        "# Phase 1 Order-Theoretic Atlas",
        "",
        "Comparison of two independent order-theoretic dimension",
        "estimators on canonical Minkowski-diamond sprinklings",
        "(manifoldlike) and Kleitman-Rothschild three-level posets",
        "(non-manifoldlike controls). No embedding, energy, or",
        "annealing is invoked in this comparison.",
        "",
        "Columns:",
        "",
        "- **family**: `minkowski` or `kleitman_rothschild`.",
        "- **d_spacetime**: the sprinkling dimension for Minkowski",
        "  causets; blank for KR (no embedding dimension is defined).",
        "- **n**: number of events.",
        "- **mean_mm**: ensemble mean of the Myrheim-Meyer dimension",
        "  over all seeds in the row.",
        "- **mean_midpoint**: ensemble mean of Meyer's midpoint",
        "  scaling dimension.",
        "- **mean_discrepancy**: ensemble mean of ``|mm - midpoint|``.",
        "",
        "Manifoldlike rows should show `mean_mm ~ mean_midpoint ~ d_spacetime`.",
        "Non-manifoldlike rows show the diagnostic separation.",
        "",
        "| family | d | n | seeds | mean MM | mean midpoint | mean discrepancy |",
        "| --- | :---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in summary:
        lines.append(
            "| {family} | {d} | {n} | {seeds} | "
            "{mm} | {mid} | {disc} |".format(
                family=row["family"],
                d=row["d_spacetime"] if row["d_spacetime"] != "" else "-",
                n=row["n"],
                seeds=row["seeds"],
                mm=_format_field(row["mean_mm"]),
                mid=_format_field(row["mean_midpoint"]),
                disc=_format_field(row["mean_discrepancy"]),
            )
        )
    lines.append("")
    lines.append(
        "Regenerate with `make regen-phase1`. The underlying tool is"
    )
    lines.append("`tools/build_phase1_atlas.py`.")
    lines.append("")
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    rows = build_rows()
    FOUNDATION.mkdir(parents=True, exist_ok=True)
    write_csv(rows, FOUNDATION / "phase1_atlas.csv")
    summary = summarize(rows)
    write_markdown(summary, FOUNDATION / "phase1_atlas.md")
    print(
        f"Wrote {len(rows)} atlas rows ({len(summary)} summary rows) "
        f"to {FOUNDATION}"
    )


if __name__ == "__main__":
    main()
