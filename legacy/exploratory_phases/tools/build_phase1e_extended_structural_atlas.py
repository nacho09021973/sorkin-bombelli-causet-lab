#!/usr/bin/env python3
"""Phase 1E — extended structural atlas for the Phase 3F dataset.

Phase 1D recorded mm_dim, midpoint_dim, chain2/3 counts and chain3
abundance over 3 seeds × 2 sizes × 3 dims (Minkowski).  Phase 3 ablations
in 3A–3E exposed two limitations:

1. The invariants needed by Phase 3 (link_count, ordering_fraction,
   relation_count, height) came from a separate JSON file that did not
   cover n=128.

2. Three seeds is too few for intra-stratum residual analysis: the
   target variance within (noise, n, d, warmup_mode) was either degenerate
   or singleton-dominated.

Phase 1E fixes both at once: a single self-contained CSV that includes
the full invariant battery used by Phase 3F, over the expanded grid

    SEEDS = 15 deterministic integers (first 5 match Phase 1B for backward
            compatibility; remaining 10 are historical year integers
            chosen for memorability).
    SIZES = (32, 64, 128)
    DIMS  = (2, 3, 4)

Only the Minkowski family is recorded — KR/corona controls were used in
Phase 1D for ensemble separation diagnostics and are not needed for the
order-theoretic residual analysis of Phase 3F.

Output columns
--------------
family, target_dim, n, seed,
mm_dim, midpoint_dim, abs_discrepancy_mm_midpoint,
chain2_count, chain3_count, chain3_abundance, chain4_count,
link_count, link_density, relation_count, ordering_fraction, height

Output
------
benchmarks/foundation/phase1e_extended_structural_atlas.csv
benchmarks/foundation/phase1e_extended_structural_atlas.md
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
import validation_suite as vs  # noqa: E402


FOUNDATION = ROOT / "benchmarks" / "foundation"

# 15 deterministic seeds; the first 5 are the Phase 1B/1D seeds so any
# downstream join across older atlases remains consistent.  The
# remaining 10 are historical year integers (chosen for human memorability).
EXTENDED_SEEDS: tuple[int, ...] = (
    1959, 1962, 1987, 2009, 2026,           # legacy block (Phase 1B/1D)
    1812, 1848, 1871, 1905, 1929,
    1945, 1968, 1989, 2001, 2017,
)
SIZES: tuple[int, ...]    = (32, 64, 128)
SPACETIME_DIMS: tuple[int, ...] = (2, 3, 4)


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
    "chain4_count",
    "link_count",
    "link_density",
    "relation_count",
    "ordering_fraction",
    "height",
)


def _row_for_cell(d: int, n: int, seed: int) -> dict:
    matrix, _ = vs.sprinkle_minkowski_diamond(n=n, seed=seed, d_spacetime=d)
    mm, mid = _estimate_dimensions(matrix)
    chains = ci.chain_counts(matrix, k_max=4)
    rc  = ci.relation_count(matrix)
    of  = ci.ordering_fraction(matrix)
    lc  = ci.link_count(matrix)
    ht  = ci.height(matrix)
    return {
        "family": "minkowski",
        "target_dim": d,
        "n": n,
        "seed": seed,
        "mm_dim": mm,
        "midpoint_dim": mid,
        "abs_discrepancy_mm_midpoint": _discrepancy(mm, mid),
        "chain2_count": chains[2],
        "chain3_count": chains[3],
        "chain3_abundance": ci.three_chain_abundance(matrix),
        "chain4_count": chains[4],
        "link_count": lc,
        "link_density": lc / n if n > 0 else 0.0,
        "relation_count": rc,
        "ordering_fraction": of,
        "height": ht,
    }


def build_rows() -> list[dict]:
    rows: list[dict] = []
    total = len(EXTENDED_SEEDS) * len(SIZES) * len(SPACETIME_DIMS)
    done = 0
    for n in SIZES:
        for seed in EXTENDED_SEEDS:
            for d in SPACETIME_DIMS:
                rows.append(_row_for_cell(d=d, n=n, seed=seed))
                done += 1
                if done % 30 == 0:
                    print(f"  {done}/{total} cells", flush=True)
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
    var = sum((v - mu) ** 2 for v in finite) / len(finite)
    return math.sqrt(var)


def summarize(rows: list[dict]) -> list[dict]:
    buckets: dict = {}
    for row in rows:
        key = (row["target_dim"], row["n"])
        buckets.setdefault(key, []).append(row)
    out = []
    for key in sorted(buckets):
        d, n = key
        rs = buckets[key]
        out.append({
            "target_dim": d,
            "n": n,
            "seeds": len(rs),
            "mean_mm":             _ensemble_mean([r["mm_dim"] for r in rs]),
            "std_mm":              _ensemble_std ([r["mm_dim"] for r in rs]),
            "mean_midpoint":       _ensemble_mean([r["midpoint_dim"] for r in rs]),
            "std_midpoint":        _ensemble_std ([r["midpoint_dim"] for r in rs]),
            "mean_abs_disc":       _ensemble_mean([r["abs_discrepancy_mm_midpoint"] for r in rs]),
            "std_abs_disc":        _ensemble_std ([r["abs_discrepancy_mm_midpoint"] for r in rs]),
            "mean_link_density":   _ensemble_mean([r["link_density"] for r in rs]),
            "std_link_density":    _ensemble_std ([r["link_density"] for r in rs]),
            "mean_ordering_frac":  _ensemble_mean([r["ordering_fraction"] for r in rs]),
            "std_ordering_frac":   _ensemble_std ([r["ordering_fraction"] for r in rs]),
        })
    return out


def write_markdown(summary: list[dict], n_rows: int, path: Path) -> None:
    lines = [
        "# Phase 1E — Extended Structural Atlas",
        "",
        "Self-contained invariant table for the Phase 3F PySR ablation.",
        "Replaces the join across `phase1d_structural_atlas.csv` +",
        "`invariants.json` with a single CSV that covers the expanded",
        "grid (15 seeds, n ∈ {32, 64, 128}, d ∈ {2, 3, 4}) and includes",
        "every invariant the Phase 3 panels need.",
        "",
        "Minkowski family only — controls (KR, corona) belong to Phase 1D.",
        "",
        f"- Total rows: {n_rows}",
        f"- Seeds ({len(EXTENDED_SEEDS)}): {', '.join(str(s) for s in EXTENDED_SEEDS)}",
        f"- Sizes: {', '.join(str(s) for s in SIZES)}",
        f"- Spacetime dims: {', '.join(str(s) for s in SPACETIME_DIMS)}",
        "",
        "## Ensemble statistics by (d, n)",
        "",
        "| d | n | seeds | mean MM | std MM | mean midpoint | std midpoint "
        "| mean \\|disc\\| | std \\|disc\\| | mean link_density | mean ordering_frac |",
        "| :---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in summary:
        lines.append(
            "| {d} | {n} | {seeds} | {mm} | {smm} | {mid} | {smid} | "
            "{disc} | {sdisc} | {ld} | {of} |".format(
                d=row["target_dim"], n=row["n"], seeds=row["seeds"],
                mm=_format_field(row["mean_mm"]),
                smm=_format_field(row["std_mm"]),
                mid=_format_field(row["mean_midpoint"]),
                smid=_format_field(row["std_midpoint"]),
                disc=_format_field(row["mean_abs_disc"]),
                sdisc=_format_field(row["std_abs_disc"]),
                ld=_format_field(row["mean_link_density"]),
                of=_format_field(row["mean_ordering_frac"]),
            )
        )
    lines += [
        "",
        "Regenerate via `make regen-phase1e`.",
        "Source: `tools/build_phase1e_extended_structural_atlas.py`.",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    print(f"Phase 1E: {len(EXTENDED_SEEDS) * len(SIZES) * len(SPACETIME_DIMS)} cells.")
    rows = build_rows()
    FOUNDATION.mkdir(parents=True, exist_ok=True)
    write_csv(rows, FOUNDATION / "phase1e_extended_structural_atlas.csv")
    write_markdown(summarize(rows), len(rows),
                   FOUNDATION / "phase1e_extended_structural_atlas.md")
    print(f"Wrote {len(rows)} rows to {FOUNDATION / 'phase1e_extended_structural_atlas.csv'}.")


if __name__ == "__main__":
    main()
