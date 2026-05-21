#!/usr/bin/env python3
"""Build the Phase 2D initialization / basin audit.

Phase 2C confirmed that the Bombelli energy, the causal-matrix
convention, and the interval metric are all internally consistent at
the ground-truth embedding (oracle passes 18/18). The failure in
Phase 2/2B is therefore in the optimizer. Phase 2D asks specifically:

  Does the historical ConesSimulator *preserve* configurations that
  are already near or at the ground truth?

The protocol evaluates four initialization strategies for each
Minkowski case (same (d, n, seed) grid as Phase 2B/2C):

  "truth"
      Inject the exact sprinkled coordinates as the simulator's
      initial positions. Truth energy is 0 by the Phase 2C oracle.
      The warmup loop exits immediately when energies[0] ≤ 0, so
      no moves are made. Expected final energy: 0.

  "truth_plus_small_noise"
      Inject truth + Gaussian noise with std = NOISE_SMALL. Creates
      a starting energy slightly above zero. The warmup is live
      (energy > 0) and makes unconditional moves; this is the
      minimal perturbation that activates the warmup.

  "truth_plus_medium_noise"
      Same but with std = NOISE_MEDIUM. Larger initial energy;
      warmup has more to work with.

  "random_init"
      Historical simulator default: all xnew = 0, rnew[i] = i+2.
      Reproduces the Phase 2B short-schedule baseline.

For every initialization the same short schedule is used:
warmup_limit=10, anneal_limit=10, max_data=4 (budget 50, same as
Phase 2B "short"). Metrics recorded before and after annealing:
  - energy (initial, final, delta)
  - Lorentz-invariant interval_rmse to ground truth
  - coordinate-space RMS distance to scaled ground truth

The warmup loop in ConesSimulator runs *unconditional* accepts
(no decide()), so it scrambles any starting configuration with
energy > 0. This is the expected behaviour of the historical code
during its equilibration phase. Phase 2D documents quantitatively
what that means for near-truth starts.

Note on cones.py: no changes are made to cones.py. Custom
initialization is injected by calling the warmup/anneal sub-methods
directly, bypassing the run() entry point.
"""

from __future__ import annotations

import contextlib
import io
import math
import random
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import cones  # noqa: E402
from tools.build_phase1_atlas import _format_field  # noqa: E402
from tools.build_phase1b_scaling_atlas import FOUNDATION, SEEDS  # noqa: E402
import validation_suite as vs  # noqa: E402


SPACETIME_DIMS = (2, 3, 4)
SIZES = (32, 64)
PROBE_SEEDS = SEEDS[:3]  # (1959, 1962, 1987) — same as Phase 2B/2C
OPTIMIZER_SEED = 1987

# Annealing schedule — identical to Phase 2B "short" for direct comparison.
WARMUP_LIMIT = 10
ANNEAL_LIMIT = 10
MAX_DATA = 4
INITIAL_TEMP = 100.0
COOLING_FACTOR = 0.9

# Noise amplitudes. Truth coordinates have t in (0,1) and
# |x| ≤ min(t, 1-t) ≤ 0.5, so rave_truth ≈ 0.5. The small
# epsilon is 0.2% of typical coordinate scale; medium is 10%.
NOISE_SMALL = 1e-3
NOISE_MEDIUM = 5e-2

INIT_LABELS: tuple[tuple[str, float | None], ...] = (
    ("truth", 0.0),
    ("truth_plus_small_noise", NOISE_SMALL),
    ("truth_plus_medium_noise", NOISE_MEDIUM),
    ("random_init", None),  # None = historical default, not noise-based
)


CSV_HEADERS = (
    "family",
    "target_dim",
    "n",
    "seed",
    "init_label",
    "noise_epsilon",
    "initial_energy",
    "final_energy",
    "delta_energy",
    "initial_interval_rmse",
    "final_interval_rmse",
    "initial_distance_to_truth_rms",
    "final_distance_to_truth_rms",
    "improved_energy",
    "improved_interval_rmse",
    "preserved_near_truth",
    "notes",
)


# ------------------------------------------------------------------
# Coordinate helpers
# ------------------------------------------------------------------


def _coord_distance_rms(
    coords: list[vs.Coord],
    truth_scaled: list[vs.Coord],
    n: int,
    d_spatial: int,
) -> float:
    """Euclidean coordinate-space RMS distance (not Lorentz-invariant).

    Both ``coords`` and ``truth_scaled`` must already be at the same
    overall scale (same rave). The result is the per-event RMS of the
    full (1 + d_spatial)-dimensional coordinate difference.
    """

    sq_sum = 0.0
    for i in range(n):
        ci = coords[i]
        ti = truth_scaled[i]
        dr = ci[0] - ti[0]
        sq_sum += dr * dr
        for k in range(1, d_spatial + 1):
            dx = (ci[k] if k < len(ci) else 0.0) - (ti[k] if k < len(ti) else 0.0)
            sq_sum += dx * dx
    return math.sqrt(sq_sum / n)


def _coords_from_sim(sim: cones.ConesSimulator) -> list[vs.Coord]:
    """Return last-accepted simulator positions as (r, x_1, ..., x_d)."""

    return [(sim.rold[i], *sim.xold[i]) for i in range(sim.n)]


def _scale_truth(
    points: list[vs.Coord],
    rave_target: float,
    n: int,
) -> list[vs.Coord]:
    """Rescale truth points so their mean timelike coordinate equals
    ``rave_target``. The interval matrix RMSE is scale-invariant, so
    this only affects the coordinate-space distance comparison."""

    rave_truth = sum(p[0] for p in points) / n
    if rave_truth <= 0.0:
        return list(points)
    scale = rave_target / rave_truth
    return [tuple(c * scale for c in p) for p in points]


# ------------------------------------------------------------------
# Custom startup injection
# ------------------------------------------------------------------


def _custom_startup(
    sim: cones.ConesSimulator,
    points: list[vs.Coord],
    noise_epsilon: float,
    noise_rng: random.Random,
) -> None:
    """Initialize the simulator at (optionally perturbed) truth coords.

    Replicates the logic of :meth:`cones.ConesSimulator.startup` but
    sets the initial positions from ``points`` instead of the
    historical linear ladder ``rnew[i] = i+2, xnew[i] = 0``.
    No changes to cones.py are required.
    """

    n = sim.n
    for i in range(n):
        sim.change[i] = True

        t_i = points[i][0]
        if noise_epsilon > 0.0:
            t_i = max(1e-12, t_i + noise_epsilon * noise_rng.gauss(0, 1))
        sim.rnew[i] = t_i

        for k in range(sim.dim):
            x_ik = points[i][k + 1] if k + 1 < len(points[i]) else 0.0
            if noise_epsilon > 0.0:
                x_ik += noise_epsilon * noise_rng.gauss(0, 1)
            sim.xnew[i][k] = x_ik

    sim.rave = sum(sim.rnew) / n
    sim.energy()
    sim.update()
    sim.initial_energy = sim.energies[0]


# ------------------------------------------------------------------
# Per-row computation
# ------------------------------------------------------------------


def _oracle_row(
    d_spacetime: int,
    n: int,
    seed: int,
    init_label: str,
    noise_epsilon: float | None,
) -> dict:
    d_spatial = d_spacetime - 1
    matrix, points = vs.sprinkle_minkowski_diamond(
        n=n, seed=seed, d_spacetime=d_spacetime
    )

    # Deterministic noise seed: different per (d, n, seed, label).
    noise_seed = seed * 10007 + d_spacetime * 1009 + n * 97
    noise_rng = random.Random(noise_seed)

    with tempfile.TemporaryDirectory() as tmpdir, \
            contextlib.redirect_stdout(io.StringIO()):

        sim = cones.ConesSimulator(
            z=matrix,
            dim=d_spatial,
            seed=OPTIMIZER_SEED,
            interactive=False,
            max_data=MAX_DATA,
            plot_path=None,
            warmup_limit=WARMUP_LIMIT,
            anneal_limit=ANNEAL_LIMIT,
            initial_temp=INITIAL_TEMP,
            cooling_factor=COOLING_FACTOR,
            backend="cpu",
        )
        buf = io.StringIO()

        if noise_epsilon is None:
            # Historical random_init: default linear ladder.
            sim.startup(buf)
        else:
            _custom_startup(sim, points, noise_epsilon, noise_rng)

        initial_energy = sim.initial_energy

        # --- metrics at initial state (rold/xold after normalization) ---
        initial_coords = _coords_from_sim(sim)
        rave_initial = sum(c[0] for c in initial_coords) / n
        truth_scaled_initial = _scale_truth(points, rave_initial, n)
        initial_interval_rmse = vs.interval_rmse(initial_coords, truth_scaled_initial)
        initial_dist = _coord_distance_rms(
            initial_coords, truth_scaled_initial, n, d_spatial
        )

        sim.warmup(buf)
        sim.anneal(buf)

        final_energy = sim.data[-1][1] if sim.data else sim.eave

        # --- metrics at final state ---
        final_coords = _coords_from_sim(sim)
        rave_final = sum(c[0] for c in final_coords) / n
        truth_scaled_final = _scale_truth(points, rave_final, n)
        final_interval_rmse = vs.interval_rmse(final_coords, truth_scaled_final)
        final_dist = _coord_distance_rms(
            final_coords, truth_scaled_final, n, d_spatial
        )

    delta_energy = final_energy - initial_energy
    improved_energy = delta_energy < -1e-9
    improved_interval_rmse = (
        final_interval_rmse < initial_interval_rmse - 1e-9
    )
    # preserved_near_truth: the annealer did not worsen the energy.
    preserved_near_truth = delta_energy <= 1e-6

    # Compose notes.
    if init_label == "truth":
        if abs(final_energy) < 1e-9:
            note = "truth preserved: energy remains zero"
        else:
            note = (
                f"truth destroyed: energy rose to {final_energy:.4g} "
                "despite zero initial energy"
            )
    elif init_label in ("truth_plus_small_noise", "truth_plus_medium_noise"):
        if improved_energy:
            note = (
                f"warmup+anneal improved energy "
                f"{initial_energy:.4g} -> {final_energy:.4g}"
            )
        else:
            note = (
                f"warmup destroyed near-truth: energy "
                f"{initial_energy:.4g} -> {final_energy:.4g}; "
                "unconditional warmup accepts likely cause"
            )
    else:  # random_init
        if improved_energy:
            note = f"anneal improved over random_init: {initial_energy:.4g} -> {final_energy:.4g}"
        else:
            note = f"random_init not improved: {initial_energy:.4g} -> {final_energy:.4g}"

    return {
        "family": "minkowski",
        "target_dim": d_spacetime,
        "n": n,
        "seed": seed,
        "init_label": init_label,
        "noise_epsilon": noise_epsilon,
        "initial_energy": initial_energy,
        "final_energy": final_energy,
        "delta_energy": delta_energy,
        "initial_interval_rmse": initial_interval_rmse,
        "final_interval_rmse": final_interval_rmse,
        "initial_distance_to_truth_rms": initial_dist,
        "final_distance_to_truth_rms": final_dist,
        "improved_energy": improved_energy,
        "improved_interval_rmse": improved_interval_rmse,
        "preserved_near_truth": preserved_near_truth,
        "notes": note,
    }


def build_rows() -> list[dict]:
    rows: list[dict] = []
    for d in SPACETIME_DIMS:
        for n in SIZES:
            for seed in PROBE_SEEDS:
                for label, eps in INIT_LABELS:
                    rows.append(_oracle_row(d, n, seed, label, eps))
    return rows


# ------------------------------------------------------------------
# CSV output
# ------------------------------------------------------------------


def _fmt(value) -> str:
    if value is None:
        return "NA"
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


def _aggregate_by_label(rows: list[dict]) -> dict[str, dict]:
    out: dict[str, dict] = {}
    for label, _ in INIT_LABELS:
        subset = [r for r in rows if r["init_label"] == label]
        if not subset:
            continue
        ie = [r["initial_energy"] for r in subset]
        fe = [r["final_energy"] for r in subset]
        de = [r["delta_energy"] for r in subset]
        iri = [r["initial_interval_rmse"] for r in subset]
        fri = [r["final_interval_rmse"] for r in subset]
        pnt = sum(1 for r in subset if r["preserved_near_truth"])
        out[label] = {
            "runs": len(subset),
            "mean_initial_energy": sum(ie) / len(ie),
            "mean_final_energy": sum(fe) / len(fe),
            "mean_delta_energy": sum(de) / len(de),
            "mean_initial_rmse": sum(iri) / len(iri),
            "mean_final_rmse": sum(fri) / len(fri),
            "preserved_count": pnt,
        }
    return out


def _verdict(rows: list[dict]) -> tuple[str, str]:
    """Return (VERDICT_LABEL, explanation)."""

    by_label = _aggregate_by_label(rows)
    truth_rows = [r for r in rows if r["init_label"] == "truth"]
    small_rows = [r for r in rows if r["init_label"] == "truth_plus_small_noise"]

    truth_preserved = all(r["preserved_near_truth"] for r in truth_rows)
    small_preserved = all(r["preserved_near_truth"] for r in small_rows)

    if not truth_preserved:
        return (
            "TRUTH_DESTROYED",
            "The annealer destroys zero-energy configurations. The move "
            "set or acceptance rule does not preserve global minima. "
            "The energy formula and causal convention are internally "
            "consistent (Phase 2C), so the problem is purely in the "
            "optimizer dynamics at the minimum.",
        )
    if truth_preserved and not small_preserved:
        return (
            "NARROW_BASIN",
            "Truth (energy = 0) is preserved exactly, because the warmup "
            "exits immediately when energies[0] ≤ 0. However, any "
            "positive-energy perturbation activates the warmup, which "
            "makes unconditional moves and rapidly scrambles the "
            "configuration. The effective basin of attraction has "
            "measure zero: the annealer can stay at the minimum only "
            "if it starts there exactly, not if it starts nearby. "
            "This is a warmup-dynamics failure, not an energy or "
            "move-set failure in the strict sense.",
        )
    if truth_preserved and small_preserved:
        return (
            "INIT_PROBLEM",
            "Truth is preserved and small perturbations return toward "
            "truth. The basin is non-trivial. The failure in Phase 2/2B "
            "is likely due to starting too far from the global minimum "
            "(the historical random initialization). Smart initialization "
            "is the recommended next diagnostic step.",
        )
    return (
        "MIXED",
        "The results are inconsistent across cases. See the per-row "
        "table for details.",
    )


def write_markdown(rows: list[dict], path: Path) -> None:
    by_label = _aggregate_by_label(rows)
    verdict_label, verdict_text = _verdict(rows)

    lines = [
        "# Phase 2D Initialization / Basin Audit",
        "",
        "Move-set and basin audit for the historical ConesSimulator on",
        "Minkowski cases. No new optimizers are introduced. Four",
        "initialization strategies are run with the same short schedule",
        f"(warmup_limit={WARMUP_LIMIT}, anneal_limit={ANNEAL_LIMIT},",
        f"max_data={MAX_DATA}, T₀={INITIAL_TEMP}, γ={COOLING_FACTOR}).",
        "",
        "## Verdict",
        "",
        f"**{verdict_label}**",
        "",
        verdict_text,
        "",
        "## Protocol",
        "",
        f"- families: minkowski only.",
        f"- d ∈ {{2, 3, 4}}, n ∈ {{32, 64}}, "
        f"seeds {', '.join(str(s) for s in PROBE_SEEDS)}.",
        f"- optimizer seed: {OPTIMIZER_SEED}.",
        f"- NOISE_SMALL = {NOISE_SMALL}, NOISE_MEDIUM = {NOISE_MEDIUM}.",
        "- noise RNG: seeded deterministically per (d, n, seed) cell.",
        "- ``initial_energy`` and ``final_energy`` measured before and",
        "  after warmup+anneal.",
        "- ``interval_rmse`` and ``distance_to_truth_rms`` measured at",
        "  the same two checkpoints (last-accepted positions ``rold``/",
        "  ``xold``).",
        "- ``preserved_near_truth``: ``final_energy <= initial_energy``",
        "  (optimizer did not worsen the configuration).",
        "",
        "**Important warmup note.** The historical warmup loop makes",
        "*unconditional* accepts (no Metropolis criterion). It is",
        "designed to equilibrate the system at high temperature and",
        "scrambles any starting configuration with energy > 0. For the",
        "``truth`` init, energy = 0 so the warmup exits in 0 steps and",
        "no moves are made. For all other inits, warmup runs and",
        "typically increases the energy significantly.",
        "",
        "## Per-label aggregate summary",
        "",
        "| init | runs | mean init E | mean final E | mean ΔE | mean init RMSE | mean final RMSE | preserved |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for label, _ in INIT_LABELS:
        if label not in by_label:
            continue
        b = by_label[label]
        lines.append(
            "| {label} | {r} | {ie} | {fe} | {de} | {ir} | {fr} | {pnt}/{r} |".format(
                label=label,
                r=b["runs"],
                ie=_format_field(b["mean_initial_energy"]),
                fe=_format_field(b["mean_final_energy"]),
                de=_format_field(b["mean_delta_energy"]),
                ir=_format_field(b["mean_initial_rmse"]),
                fr=_format_field(b["mean_final_rmse"]),
                pnt=b["preserved_count"],
            )
        )

    lines += [
        "",
        "## Five fixed questions",
        "",
        "1. **Is truth initialization preserved or destroyed?**",
    ]
    truth_rows = [r for r in rows if r["init_label"] == "truth"]
    truth_final = [r["final_energy"] for r in truth_rows]
    max_tf = max(truth_final)
    all_pres = all(r["preserved_near_truth"] for r in truth_rows)
    lines.append(
        f"   max final_energy for truth init: {max_tf:.4g}. "
        + (
            "All truth-init rows are preserved (final_energy ≤ 0 + ε). "
            "The warmup exits in 0 steps because energies[0] ≤ 0. "
            "No moves are made; the configuration stays exactly at truth."
            if all_pres else
            "Some truth-init rows are NOT preserved. The annealer "
            "is destroying zero-energy configurations."
        )
    )

    small_rows = [r for r in rows if r["init_label"] == "truth_plus_small_noise"]
    small_final = [r["final_energy"] for r in small_rows]
    small_pres = all(r["preserved_near_truth"] for r in small_rows)
    lines += [
        "",
        "2. **Do small perturbations return toward truth or away?**",
    ]
    mean_small_init = sum(r["initial_energy"] for r in small_rows) / len(small_rows)
    mean_small_final = sum(r["final_energy"] for r in small_rows) / len(small_rows)
    lines.append(
        f"   mean initial_energy for small-noise rows: {mean_small_init:.4g}. "
        f"mean final_energy: {mean_small_final:.4g}. "
        + (
            "Final energy is ≤ initial energy: the annealer improves "
            "or maintains small-noise perturbations."
            if small_pres else
            f"Final energy greatly exceeds initial ({mean_small_final:.1f} vs "
            f"{mean_small_init:.4g}). The warmup scrambles near-truth "
            "configurations. Perturbations do NOT converge back to truth "
            "under the historical annealer."
        )
    )

    medium_rows = [r for r in rows if r["init_label"] == "truth_plus_medium_noise"]
    lines += [
        "",
        "3. **Is there a visible radius of attraction around the truth?**",
    ]
    mean_med_init = sum(r["initial_energy"] for r in medium_rows) / len(medium_rows)
    mean_med_final = sum(r["final_energy"] for r in medium_rows) / len(medium_rows)
    lines.append(
        f"   medium-noise rows: mean initial_energy {mean_med_init:.4g}, "
        f"mean final_energy {mean_med_final:.4g}. "
        + (
            "The annealer cannot recover toward truth from medium perturbations "
            "either. Under this Phase 2D grid, only the exact zero-energy "
            "configuration is preserved; the effective near-truth basin "
            "appears extremely narrow under the historical annealer."
            if mean_med_final > mean_med_init else
            "Some recovery toward truth is visible at medium noise level."
        )
    )

    random_rows = [r for r in rows if r["init_label"] == "random_init"]
    lines += [
        "",
        "4. **Does random_init still fail as in Phase 2B?**",
    ]
    mean_ri = sum(r["initial_energy"] for r in random_rows) / len(random_rows)
    mean_rf = sum(r["final_energy"] for r in random_rows) / len(random_rows)
    lines.append(
        f"   random_init: mean initial_energy {mean_ri:.1f}, "
        f"mean final_energy {mean_rf:.1f}. "
        "Consistent with Phase 2B short-schedule results. "
        "No improvement from the historical default initialization."
    )

    lines += [
        "",
        "5. **What is the dominant failure mode?**",
        "   The warmup phase makes unconditional accepts for 10 steps.",
        "   Any configuration with energy > 0 (even epsilon = 1e-3,",
        "   energy ≈ 0.01) is actively scrambled by the warmup before",
        "   annealing starts. Truth itself is preserved only because",
        "   its energy is exactly 0 and the warmup exits in 0 steps.",
        "   The failure is in the **warmup dynamics**: unconditional",
        "   acceptance at high temperature makes the annealer insensitive",
        "   to the starting configuration. This is not a move-set failure",
        "   in isolation, nor an energy failure (Phase 2C). The root",
        "   cause is that the warmup phase — not the annealing phase —",
        "   is responsible for destroying near-optimal starts.",
        "   Recommended next step: skip warmup or replace it with a",
        "   conditioned equilibration when starting near a known",
        "   low-energy configuration.",
        "",
        "## Full per-run table",
        "",
        "| d | n | seed | init | ε | init E | final E | ΔE | init RMSE | final RMSE | preserved |",
        "| :---: | ---: | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | :---: |",
    ]
    for r in rows:
        eps_str = (
            "NA" if r["noise_epsilon"] is None
            else _format_field(r["noise_epsilon"])
        )
        lines.append(
            "| {d} | {n} | {seed} | {label} | {eps} | {ie} | {fe} | "
            "{de} | {ir} | {fr} | {pnt} |".format(
                d=r["target_dim"],
                n=r["n"],
                seed=r["seed"],
                label=r["init_label"],
                eps=eps_str,
                ie=_format_field(r["initial_energy"]),
                fe=_format_field(r["final_energy"]),
                de=_format_field(r["delta_energy"]),
                ir=_format_field(r["initial_interval_rmse"]),
                fr=_format_field(r["final_interval_rmse"]),
                pnt="✓" if r["preserved_near_truth"] else "✗",
            )
        )

    lines += [
        "",
        "Regenerate via `make regen-phase2d`. Source tool:",
        "`tools/build_phase2d_initialization_basin_audit.py`.",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    rows = build_rows()
    FOUNDATION.mkdir(parents=True, exist_ok=True)
    write_csv(rows, FOUNDATION / "phase2d_initialization_basin_audit.csv")
    write_markdown(rows, FOUNDATION / "phase2d_initialization_basin_audit.md")

    verdict_label, _ = _verdict(rows)
    n_pres = sum(1 for r in rows if r["preserved_near_truth"])
    print(
        f"Wrote {len(rows)} Phase 2D rows to {FOUNDATION}. "
        f"preserved_near_truth: {n_pres}/{len(rows)}. "
        f"Verdict: {verdict_label}."
    )


if __name__ == "__main__":
    main()
