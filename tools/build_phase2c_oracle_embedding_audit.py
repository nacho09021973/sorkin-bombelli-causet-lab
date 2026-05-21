#!/usr/bin/env python3
"""Build the Phase 2C oracle embedding audit.

Phase 2C answers the question left open after Phase 2B: given that
increasing the annealer budget does not close the Minkowski energy gap,
is the failure in the energy definition / causal-matrix convention /
interval metric, or is it purely in the optimizer?

The oracle protocol is:

1. For each Minkowski case (no annealing), take the ground-truth
   coordinates directly from the sprinkler.
2. Evaluate the Bombelli energy at those coordinates via
   :func:`validation_suite.bombelli_energy_at`. If the energy
   function and causal-matrix convention are mutually consistent the
   result must be zero.
3. Reconstruct the causal matrix from the ground-truth coordinates
   using the same criterion as the sprinkler. Count discordant pairs
   against the stored causal matrix to check for floating-point drift
   or convention mismatches.
4. Evaluate :func:`validation_suite.interval_rmse` comparing the
   ground-truth embedding to itself. This must be zero by definition;
   a nonzero result would indicate a broken residual formula.
5. Record pass/fail flags and the energy-function normalization
   constant for auditing.

Nothing is optimized. No random state is touched after the sprinkler.
The result is a pure algebraic consistency check:

- oracle_pass_energy: Bombelli energy at ground truth is
  numerically compatible with zero.
- oracle_pass_causal_matrix: no discordant pairs between the stored
  causal matrix and the one reconstructed from the stored coordinates.
- oracle_pass_interval_rmse: the Lorentz-invariant RMSE of a point
  set against itself is numerically zero.

If all three pass, the failure in Phase 2/2B is localized to the
optimizer (move set, initialization, or annealing landscape). If any
fails, there is a convention or formula inconsistency that must be
fixed before running more annealing.

Phase 2C covers the same (d, n, seed) grid as Phase 2B: d in {2,3,4},
n in {32, 64}, seeds = the first three Phase 1B atlas seeds.
Non-manifoldlike controls are excluded; they have no ground-truth
coordinates.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from tools.build_phase1_atlas import _format_field  # noqa: E402
from tools.build_phase1b_scaling_atlas import FOUNDATION, SEEDS  # noqa: E402
import validation_suite as vs  # noqa: E402


SPACETIME_DIMS = (2, 3, 4)
SIZES = (32, 64)
PROBE_SEEDS = SEEDS[:3]  # (1959, 1962, 1987) — same slice as Phase 2B

# Thresholds for the three oracle pass criteria.
# Energy and interval RMSE should be identically zero; the small
# tolerances absorb any floating-point edge cases without weakening
# the conclusion when the result is, as expected, exactly zero.
ENERGY_TOL = 1e-9
RMSE_TOL = 1e-9


CSV_HEADERS = (
    "family",
    "target_dim",
    "n",
    "seed",
    "d_spatial",
    "pair_count",
    "oracle_energy",
    "oracle_energy_abs",
    "rave_truth",
    "oracle_interval_rmse",
    "original_pair_count",
    "reconstructed_pair_count",
    "false_positive_pairs",
    "false_negative_pairs",
    "total_discordant_pairs",
    "oracle_pass_energy",
    "oracle_pass_causal_matrix",
    "oracle_pass_interval_rmse",
    "notes",
)


def _reconstruct_causal_matrix(
    points: list[vs.Coord],
    n: int,
) -> list[list[bool]]:
    """Rebuild the causal matrix from ground-truth coordinates.

    Uses the identical floating-point criterion as
    :func:`validation_suite.sprinkle_minkowski_diamond`:
    pair (i, j) is causal iff ``dt * dt >= dx_sq``.
    """

    matrix: list[list[bool]] = [[False] * n for _ in range(n)]
    for i in range(n - 1):
        ti = points[i][0]
        xi = points[i][1:]
        for j in range(i + 1, n):
            tj = points[j][0]
            xj = points[j][1:]
            dt = tj - ti
            dx_sq = sum((b - a) * (b - a) for a, b in zip(xi, xj))
            if dt * dt >= dx_sq:
                matrix[i][j] = True
    return matrix


def _count_relations(matrix: list[list[bool]], n: int) -> int:
    return sum(1 for i in range(n - 1) for j in range(i + 1, n) if matrix[i][j])


def _discordant_pairs(
    z_orig: list[list[bool]],
    z_recon: list[list[bool]],
    n: int,
) -> tuple[int, int]:
    """Return (false_positives, false_negatives).

    false_positive: orig=True, recon=False (stored as causal but not
    recoverable from the stored coordinates).
    false_negative: orig=False, recon=True (stored as non-causal but
    coordinates imply causal).
    """

    fp = fn = 0
    for i in range(n - 1):
        for j in range(i + 1, n):
            if z_orig[i][j] and not z_recon[i][j]:
                fp += 1
            elif not z_orig[i][j] and z_recon[i][j]:
                fn += 1
    return fp, fn


def _oracle_row(d_spacetime: int, n: int, seed: int) -> dict:
    d_spatial = d_spacetime - 1
    pair_count = n * (n - 1) // 2

    matrix, points = vs.sprinkle_minkowski_diamond(
        n=n, seed=seed, d_spacetime=d_spacetime
    )

    # 1. Bombelli energy at ground-truth coordinates (no annealing).
    rave_truth = sum(p[0] for p in points) / n
    oracle_energy = vs.bombelli_energy_at(matrix, points, d_spatial=d_spatial)
    oracle_energy_abs = abs(oracle_energy)

    # 2. Lorentz-invariant interval RMSE of the truth against itself.
    oracle_interval_rmse = vs.interval_rmse(points, points)

    # 3. Causal matrix reconstruction from stored coordinates.
    z_recon = _reconstruct_causal_matrix(points, n)
    original_pair_count = _count_relations(matrix, n)
    reconstructed_pair_count = _count_relations(z_recon, n)
    fp, fn = _discordant_pairs(matrix, z_recon, n)
    total_discordant = fp + fn

    # Pass criteria.
    pass_energy = oracle_energy_abs <= ENERGY_TOL
    pass_matrix = total_discordant == 0
    pass_rmse = oracle_interval_rmse <= RMSE_TOL

    # Notes: report per-case anomalies if any.
    anomalies: list[str] = []
    if not pass_energy:
        anomalies.append(
            f"energy={oracle_energy:.6g} exceeds tolerance {ENERGY_TOL}"
        )
    if not pass_matrix:
        anomalies.append(
            f"{fp} false-positive and {fn} false-negative pairs"
        )
    if not pass_rmse:
        anomalies.append(
            f"interval_rmse={oracle_interval_rmse:.6g} exceeds {RMSE_TOL}"
        )
    notes = "; ".join(anomalies) if anomalies else "all oracle checks pass"

    return {
        "family": "minkowski",
        "target_dim": d_spacetime,
        "n": n,
        "seed": seed,
        "d_spatial": d_spatial,
        "pair_count": pair_count,
        "oracle_energy": oracle_energy,
        "oracle_energy_abs": oracle_energy_abs,
        "rave_truth": rave_truth,
        "oracle_interval_rmse": oracle_interval_rmse,
        "original_pair_count": original_pair_count,
        "reconstructed_pair_count": reconstructed_pair_count,
        "false_positive_pairs": fp,
        "false_negative_pairs": fn,
        "total_discordant_pairs": total_discordant,
        "oracle_pass_energy": pass_energy,
        "oracle_pass_causal_matrix": pass_matrix,
        "oracle_pass_interval_rmse": pass_rmse,
        "notes": notes,
    }


def build_rows() -> list[dict]:
    rows: list[dict] = []
    for d in SPACETIME_DIMS:
        for n in SIZES:
            for seed in PROBE_SEEDS:
                rows.append(_oracle_row(d, n, seed))
    return rows


def _fmt(value) -> str:
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, str):
        return value
    if isinstance(value, int):
        return str(value)
    return _format_field(value)


def write_csv(rows: list[dict], path: Path) -> None:
    lines = [",".join(CSV_HEADERS)]
    for row in rows:
        lines.append(",".join(_fmt(row[h]) for h in CSV_HEADERS))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


# ------------------------------------------------------------------
# Markdown report
# ------------------------------------------------------------------


def _verdict(rows: list[dict]) -> tuple[str, str]:
    """Return (label, explanation) summarising the oracle verdict."""

    all_energy = all(r["oracle_pass_energy"] for r in rows)
    all_matrix = all(r["oracle_pass_causal_matrix"] for r in rows)
    all_rmse = all(r["oracle_pass_interval_rmse"] for r in rows)

    if all_energy and all_matrix and all_rmse:
        return (
            "ORACLE PASSES",
            "The ground-truth coordinates are a zero-energy "
            "configuration. The causal matrices reconstructed from "
            "those coordinates match the stored matrices exactly. "
            "The Lorentz-invariant interval residual of the truth "
            "against itself is zero. These three checks confirm that "
            "the Bombelli energy, the causal-matrix convention, and "
            "the interval metric are mutually consistent at the "
            "ground-truth embedding. The failure in Phase 2/2B is "
            "therefore localized to the optimizer: move set, "
            "initialization, or annealing landscape. There is no "
            "convention error to fix before running more optimization.",
        )

    parts: list[str] = []
    if not all_energy:
        parts.append("energy_oracle_fails")
    if not all_matrix:
        parts.append("causal_matrix_mismatch")
    if not all_rmse:
        parts.append("rmse_oracle_fails")
    label = "ORACLE FAILS (" + ", ".join(parts) + ")"
    explanation = (
        "At least one oracle check failed. This indicates a "
        "convention or formula inconsistency between the ground-truth "
        "coordinates, the stored causal matrix, the Bombelli energy "
        "objective, or the Lorentz-invariant interval metric. "
        "Running more annealing will not resolve an inconsistency at "
        "this level. The next step is an audit of the relevant "
        "convention (energy formula, causal criterion, or interval "
        "sign convention) before any further optimization."
    )
    return label, explanation


def write_markdown(rows: list[dict], path: Path) -> None:
    verdict_label, verdict_text = _verdict(rows)

    lines = [
        "# Phase 2C Oracle Embedding Audit",
        "",
        "Ground-truth consistency check for the Minkowski cases in",
        "Phase 2/2B. No annealing is performed. The ground-truth",
        "coordinates from the canonical sprinkler are evaluated",
        "directly against the Bombelli energy, the causal-matrix",
        "convention, and the Lorentz-invariant interval metric.",
        "",
        "## Verdict",
        "",
        f"**{verdict_label}**",
        "",
        verdict_text,
        "",
        "## Protocol",
        "",
        f"- family: minkowski only.",
        f"- dimensions: d_spacetime ∈ {{2, 3, 4}}.",
        f"- sizes: n ∈ {{32, 64}}.",
        f"- case seeds: {', '.join(str(s) for s in PROBE_SEEDS)} "
        "(first three Phase 1B atlas seeds, same as Phase 2B).",
        f"- energy tolerance: |oracle_energy| ≤ {ENERGY_TOL}.",
        f"- RMSE tolerance: oracle_interval_rmse ≤ {RMSE_TOL}.",
        "- causal-matrix pass: total_discordant_pairs = 0.",
        "",
        "Three oracle checks per row (no optimisation):",
        "",
        "1. **oracle_pass_energy** — Bombelli energy at ground-truth",
        "   coordinates is numerically zero. True iff the energy",
        "   formula and the causal-matrix convention are mutually",
        "   consistent.",
        "2. **oracle_pass_causal_matrix** — causal matrix reconstructed",
        "   from the stored coordinates matches the stored matrix",
        "   bit-for-bit. True iff no floating-point drift or sign",
        "   convention mismatch occurred between construction and",
        "   storage.",
        "3. **oracle_pass_interval_rmse** — Lorentz-invariant RMSE of",
        "   the truth embedding against itself is numerically zero.",
        "   True iff the interval-matrix formula is self-consistent.",
        "",
        "## Results",
        "",
        "| d | n | seed | oracle E | |oracle E| | RMSE | discordant | E pass | M pass | R pass |",
        "| :---: | ---: | ---: | ---: | ---: | ---: | ---: | :---: | :---: | :---: |",
    ]
    for r in rows:
        lines.append(
            "| {d} | {n} | {seed} | {oe} | {oea} | {rmse} | {disc} "
            "| {pe} | {pm} | {pr} |".format(
                d=r["target_dim"],
                n=r["n"],
                seed=r["seed"],
                oe=_format_field(r["oracle_energy"]),
                oea=_format_field(r["oracle_energy_abs"]),
                rmse=_format_field(r["oracle_interval_rmse"]),
                disc=r["total_discordant_pairs"],
                pe="✓" if r["oracle_pass_energy"] else "✗",
                pm="✓" if r["oracle_pass_causal_matrix"] else "✗",
                pr="✓" if r["oracle_pass_interval_rmse"] else "✗",
            )
        )

    n_pass = sum(
        1 for r in rows
        if r["oracle_pass_energy"]
        and r["oracle_pass_causal_matrix"]
        and r["oracle_pass_interval_rmse"]
    )

    lines += [
        "",
        f"All-three-pass count: {n_pass} / {len(rows)}.",
        "",
        "## Five fixed questions",
        "",
        "1. **Does the ground-truth embedding give energy zero?**",
    ]
    energy_vals = [r["oracle_energy"] for r in rows]
    max_e = max(abs(e) for e in energy_vals)
    lines.append(
        f"   Maximum |oracle_energy| across all {len(rows)} cases: "
        f"{max_e:.6g}. "
        + (
            "Every case satisfies oracle_pass_energy. The Bombelli "
            "energy formula returns exactly 0.0 at the ground-truth "
            "coordinates, confirming internal consistency between the "
            "causal matrix and the energy objective."
            if all(r["oracle_pass_energy"] for r in rows) else
            "At least one case fails. The energy is not zero at the "
            "ground truth; there is a mismatch between the causal "
            "matrix and the energy formula."
        )
    )
    lines += [
        "",
        "2. **Does the reconstructed causal matrix match the stored one?**",
    ]
    total_disc = sum(r["total_discordant_pairs"] for r in rows)
    lines.append(
        f"   Total discordant pairs across all rows: {total_disc}. "
        + (
            "The causal criterion used by the sprinkler and the "
            "criterion used in reconstruction are identical in "
            "floating-point. No drift or convention mismatch."
            if total_disc == 0 else
            f"Discordant pairs detected. This means the stored causal "
            f"matrix cannot be reproduced from the stored coordinates "
            f"using the project's causal criterion. The source of the "
            f"mismatch must be audited before any further embedding "
            f"work."
        )
    )
    lines += [
        "",
        "3. **Are the interval residuals zero at the ground truth?**",
    ]
    max_rmse = max(r["oracle_interval_rmse"] for r in rows)
    lines.append(
        f"   Maximum oracle_interval_rmse: {max_rmse:.6g}. "
        + (
            "The Lorentz-invariant RMSE of a point set against itself "
            "is zero, confirming that the interval-matrix formula is "
            "self-consistent."
            if all(r["oracle_pass_interval_rmse"] for r in rows) else
            "Nonzero self-RMSE detected. The interval formula has a "
            "numerical inconsistency at the ground truth."
        )
    )
    lines += [
        "",
        "4. **If something fails — what is the failure mode?**",
        "   " + (
            "No failures in this run. All oracle checks pass. There is "
            "no evidence of a convention error in the energy formula, "
            "the causal criterion, or the interval metric."
            if n_pass == len(rows) else
            "See the discordant-pairs columns and the 'notes' field in "
            "the CSV for per-case failure attribution."
        ),
        "",
        "5. **If oracle passes — what does that tell us about Phase 2B?**",
        "   " + (
            "The energy objective recognises the ground-truth solution "
            "as a zero-energy configuration. The annealing failure in "
            "Phase 2/2B is therefore not caused by a broken energy "
            "formula. The optimizer does not recover the oracle "
            "configuration under the current move set / landscape because "
            "of move-set or landscape issues, not because the target "
            "is wrong. The next diagnostic step is a move-set or "
            "initialization audit — not more budget, and not an energy "
            "redesign."
            if n_pass == len(rows) else
            "Oracle does not fully pass. Phase 2B conclusions are "
            "confounded by the formula inconsistency. Fix the "
            "convention before interpreting annealing outcomes."
        ),
        "",
        "## Normalization audit",
        "",
        "| d | n | seed | rave_truth | pair_count | original_pairs | reconstructed_pairs |",
        "| :---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for r in rows:
        lines.append(
            "| {d} | {n} | {seed} | {rave} | {pc} | {op} | {rp} |".format(
                d=r["target_dim"],
                n=r["n"],
                seed=r["seed"],
                rave=_format_field(r["rave_truth"]),
                pc=r["pair_count"],
                op=r["original_pair_count"],
                rp=r["reconstructed_pair_count"],
            )
        )

    lines += [
        "",
        "Regenerate via `make regen-phase2c`. Source tool:",
        "`tools/build_phase2c_oracle_embedding_audit.py`.",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    rows = build_rows()
    FOUNDATION.mkdir(parents=True, exist_ok=True)
    write_csv(rows, FOUNDATION / "phase2c_oracle_embedding_audit.csv")
    write_markdown(rows, FOUNDATION / "phase2c_oracle_embedding_audit.md")

    n_all_pass = sum(
        1 for r in rows
        if r["oracle_pass_energy"]
        and r["oracle_pass_causal_matrix"]
        and r["oracle_pass_interval_rmse"]
    )
    verdict_label, _ = _verdict(rows)
    print(
        f"Wrote {len(rows)} Phase 2C oracle rows to {FOUNDATION}. "
        f"All-three-pass: {n_all_pass}/{len(rows)}. "
        f"Verdict: {verdict_label}."
    )


if __name__ == "__main__":
    main()
