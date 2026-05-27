#!/usr/bin/env python3
"""S4-KERR-K3-LOCAL-CONE-001: Kerr equatorial local null-cone diagnostic.

This is the first genuine Kerr-geometry diagnostic in the SORKIN-4 program.
It computes Boyer-Lindquist equatorial metric coefficients at sampled exterior
points and classifies small-displacement local intervals by the sign of ds².

BOUNDARY (read before extending this module):
  local_timelike_candidate  ≠  causal_relation = True
  local_spacelike_candidate ≠  causal_relation = False

The labels are local metric-sign diagnostics ONLY.  They do NOT:
  - imply global causal reachability between the two events,
  - integrate Kerr null geodesics of any kind,
  - decide prograde or retrograde causal relations,
  - constitute a Kerr causal solver of any kind.

For a=0, the Schwarzschild BL reduction is verified and the a=0 global causal
control counts from K1/K2 are preserved.  For a>0, all global causal pairs
remain undecided by the same invariant as K1/K2.
"""

from __future__ import annotations

import csv
import json
import math
import sys
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Optional

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

sys.path.insert(0, str(Path(__file__).resolve().parent))

from run_kerr_minimal_benchmark import (          # noqa: E402
    Event,
    R_MAX,
    T_MIN,
    T_MAX,
    build_relation_states,
    count_true_relations,
    generate_exterior_events,
    kerr_horizon_radius,
    signed_delta_phi,
)
from explore.sorkin4_schwarzschild_benchmark import (  # noqa: E402
    run_schwarzschild_minimal_benchmark as schwarz,
)


# ---------------------------------------------------------------------------
# Audit identity
# ---------------------------------------------------------------------------

AUDIT_ID       = "S4-KERR-K3-LOCAL-CONE-001"
OUT_DIR        = Path(__file__).resolve().parent
DEFAULT_OUT_PREFIX = "kerr_k3_local_cone_001_n12_seed1959"
COMMAND        = (
    "python3 explore/sorkin4_kerr_benchmark/audit_kerr_k3_local_cone_001.py"
)

# ---------------------------------------------------------------------------
# Sweep parameters
# ---------------------------------------------------------------------------

DEFAULT_N      = 12
DEFAULT_SEED   = 1959
DEFAULT_MASS   = 1.0
DEFAULT_MARGIN = schwarz.EXTERIOR_MARGIN           # 0.35, same as K1/K2
DEFAULT_SPINS  = (0.0, 0.25, 0.5, 0.75)

# ---------------------------------------------------------------------------
# Local-pair filter thresholds (Adjustment 1 from design review)
#
# Only pairs satisfying BOTH conditions are classified as local samples:
#   |dr|   <= DR_LOCAL_THRESHOLD   (Boyer-Lindquist coordinate distance)
#   |dphi| <= DPHI_LOCAL_THRESHOLD (shortest azimuthal separation, radians)
#
# Pairs outside these thresholds are counted as local_skipped_pair_count.
# ---------------------------------------------------------------------------

DR_LOCAL_THRESHOLD   = 1.0    # BL coordinate units
DPHI_LOCAL_THRESHOLD = 0.5    # radians (~28.6°)

# ---------------------------------------------------------------------------
# Numerical tolerances
# ---------------------------------------------------------------------------

METRIC_TOL    = 1.0e-12    # Schwarzschild-reduction check tolerance
LOCAL_DS2_EPS = 1.0e-9     # scale factor for ds² classifier

# ---------------------------------------------------------------------------
# CSV schema (one row per spin value)
# ---------------------------------------------------------------------------

CSV_FIELDS = (
    "spin_a",
    "M",
    "N",
    "seed",
    "margin",
    "r_plus",
    "r_erg_eq",
    "r_min_observed",
    "all_points_exterior",
    "min_Delta",
    "min_g_rr",
    "min_g_phiphi",
    "max_abs_g_tphi",
    "local_evaluated_pair_count",
    "local_skipped_pair_count",
    "local_max_abs_dr",
    "local_max_abs_dphi",
    "local_max_abs_dt",
    "local_timelike_count",
    "local_nullish_count",
    "local_spacelike_count",
    "global_true_relations",
    "global_false_relations",
    "global_undecided_pairs",
    "schwarzschild_reduction_pass",
    "frame_dragging_sign_pass",
    "all_checks_pass",
)


# ---------------------------------------------------------------------------
# Metric helpers
# ---------------------------------------------------------------------------

def equatorial_metric_at_r(r: float, mass: float, spin: float) -> dict[str, float]:
    """Boyer-Lindquist equatorial metric coefficients at theta=pi/2, signature -+++.

    Equatorial simplifications (theta=pi/2 → sin=1, cos=0):
      Sigma_eq = r²
      Delta    = r² - 2Mr + a²
      g_tt     = -(1 - 2M/r)
      g_tphi   = -2Ma/r
      g_rr     = r²/Delta
      g_phiphi = r² + a² + 2Ma²/r

    These reduce to standard Schwarzschild BL coefficients when a=0:
      g_tt     = -(1 - 2M/r)
      g_tphi   = 0
      g_rr     = 1/(1 - 2M/r)
      g_phiphi = r²

    Raises ValueError if Delta <= 0 (point inside or on the Kerr horizon).
    Delta is computed explicitly here (Adjustment 2 from design review) so
    that the pre-check is visible in the call stack.
    """
    delta = r * r - 2.0 * mass * r + spin * spin
    if delta <= 0.0:
        raise ValueError(
            f"equatorial_metric_at_r: Delta={delta:.6g} <= 0 "
            f"at r={r:.6g}, M={mass}, a={spin}; "
            "point is inside or on the Kerr horizon."
        )
    return {
        "g_tt":     -(1.0 - 2.0 * mass / r),
        "g_tphi":   -2.0 * mass * spin / r,
        "g_rr":     r * r / delta,
        "g_phiphi": r * r + spin * spin + 2.0 * mass * spin * spin / r,
        "delta":    delta,
    }


def is_local_pair(p: Event, q: Event) -> bool:
    """Return True iff (p, q) satisfies both local-displacement thresholds."""
    return (
        abs(q.r - p.r) <= DR_LOCAL_THRESHOLD
        and abs(signed_delta_phi(p.phi, q.phi)) <= DPHI_LOCAL_THRESHOLD
    )


def local_interval_ds2_equatorial(
    p: Event,
    q: Event,
    mass: float,
    spin: float,
) -> tuple[float, dict[str, float]]:
    """ds² at the midpoint radius for an equatorial pair in BL coordinates.

    ds² = g_tt dt² + 2 g_tphi dt dphi + g_rr dr² + g_phiphi dphi²

    dtheta = 0 by construction (all K3 events are at theta = pi/2).
    The metric is evaluated at r_mid = (r_i + r_j) / 2.

    Returns (ds2, metric_at_midpoint).
    """
    r_mid  = 0.5 * (p.r + q.r)
    metric = equatorial_metric_at_r(r_mid, mass, spin)
    dt     = q.t - p.t
    dr     = q.r - p.r
    dphi   = signed_delta_phi(p.phi, q.phi)
    ds2 = (
        metric["g_tt"]     * dt   * dt
        + 2.0 * metric["g_tphi"] * dt   * dphi
        + metric["g_rr"]   * dr   * dr
        + metric["g_phiphi"] * dphi * dphi
    )
    return ds2, metric


def classify_local_interval(
    ds2:        float,
    scale:      float,
    eps_factor: float = LOCAL_DS2_EPS,
) -> str:
    """Classify a local ds² value as one of three local-diagnostic labels.

    The tolerance is scale-adaptive: tol = eps_factor * max(1, |scale|).

    Returns one of:
      'timelike_local_candidate'   — ds² < -tol
      'nullish_local_candidate'    — |ds²| <= tol
      'spacelike_local_candidate'  — ds² > tol

    These labels are NOT global causal relation decisions.
    """
    tol = eps_factor * max(1.0, abs(scale))
    if ds2 < -tol:
        return "timelike_local_candidate"
    if ds2 > tol:
        return "spacelike_local_candidate"
    return "nullish_local_candidate"


# ---------------------------------------------------------------------------
# Per-spin case runner
# ---------------------------------------------------------------------------

def run_spin_case(
    n:      int,
    seed:   int,
    mass:   float,
    spin:   float,
    margin: float,
) -> dict[str, Any]:
    """Run one (spin, N, seed, margin) K3 diagnostic cell.

    Returns a dict with all CSV_FIELDS plus internal '_pair_details' and
    '_events' keys (stripped before writing to CSV).
    """
    r_plus = kerr_horizon_radius(mass, spin)
    r_min  = r_plus + margin

    # Generate equatorial events (theta = pi/2 forced)
    events = generate_exterior_events(n, seed, r_min, equatorial=True)

    # --- Adjustment 2: compute Delta explicitly per event before metric call ---
    event_deltas:   list[float] = []
    event_g_rr:     list[float] = []
    event_g_phiphi: list[float] = []
    event_g_tphi:   list[float] = []

    for ev in events:
        # Explicit Delta first — traceable in the artifact even if sampling changes.
        delta = ev.r * ev.r - 2.0 * mass * ev.r + spin * spin
        event_deltas.append(delta)
        m = equatorial_metric_at_r(ev.r, mass, spin)
        event_g_rr.append(m["g_rr"])
        event_g_phiphi.append(m["g_phiphi"])
        event_g_tphi.append(m["g_tphi"])

    # --- Adjustment 3: causal accounting per K1/K2 rules ---
    # equatorial_scaffold mode: a=0 → Schwarzschild control counts;
    #                           a>0 → all pairs undecided (by construction).
    matrix, states = build_relation_states(events, mass, spin, "equatorial_scaffold")
    possible_pairs  = n * (n - 1) // 2
    true_relations  = count_true_relations(matrix)
    false_relations = sum(
        1 for i in range(n - 1) for j in range(i + 1, n)
        if states[i][j] is False
    )
    undecided_pairs = sum(
        1 for i in range(n - 1) for j in range(i + 1, n)
        if states[i][j] is None
    )

    # --- Local cone classification (Adjustment 1: small-displacement only) ---
    local_evaluated:   int   = 0
    local_skipped:     int   = 0
    local_timelike:    int   = 0
    local_nullish:     int   = 0
    local_spacelike:   int   = 0
    max_abs_dr:        float = 0.0
    max_abs_dphi:      float = 0.0
    max_abs_dt:        float = 0.0
    pair_details: list[dict[str, Any]] = []

    for i in range(n - 1):
        for j in range(i + 1, n):
            p, q = events[i], events[j]
            if not is_local_pair(p, q):
                local_skipped += 1
                continue

            dr_val   = q.r - p.r
            dphi_val = signed_delta_phi(p.phi, q.phi)
            dt_val   = q.t - p.t

            max_abs_dr   = max(max_abs_dr,   abs(dr_val))
            max_abs_dphi = max(max_abs_dphi, abs(dphi_val))
            max_abs_dt   = max(max_abs_dt,   abs(dt_val))

            ds2, metric_mid = local_interval_ds2_equatorial(p, q, mass, spin)

            # Scale-adaptive tolerance: key metric terms at midpoint
            scale = max(
                1.0,
                abs(metric_mid["g_tt"]     * dt_val  * dt_val),
                abs(metric_mid["g_phiphi"] * dphi_val * dphi_val),
                abs(2.0 * metric_mid["g_tphi"] * dt_val * dphi_val),
            )
            label = classify_local_interval(ds2, scale)

            if label == "timelike_local_candidate":
                local_timelike  += 1
            elif label == "nullish_local_candidate":
                local_nullish   += 1
            else:
                local_spacelike += 1
            local_evaluated += 1

            pair_details.append({
                "i":     i,
                "j":     j,
                "dr":    dr_val,
                "dphi":  dphi_val,
                "dt":    dt_val,
                "r_mid": 0.5 * (p.r + q.r),
                "ds2":   ds2,
                "label": label,
            })

    # --- Metric sanity checks ---
    r_min_observed      = min(ev.r for ev in events)
    all_points_exterior = r_min_observed > r_min

    # Schwarzschild reduction check (a=0 only)
    schwarzschild_reduction_pass = False
    if abs(spin) <= 0.0:
        g_tphi_zero = all(abs(g) <= METRIC_TOL for g in event_g_tphi)
        g_rr_correct = all(
            abs(event_g_rr[k] - 1.0 / (1.0 - 2.0 * mass / events[k].r)) <= METRIC_TOL
            for k in range(n)
        )
        g_phiphi_correct = all(
            abs(event_g_phiphi[k] - events[k].r ** 2) <= METRIC_TOL
            for k in range(n)
        )
        schwarzschild_reduction_pass = g_tphi_zero and g_rr_correct and g_phiphi_correct

    # Frame-dragging sign check (a>0 only): g_tphi = -2Ma/r < 0 for a>0
    frame_dragging_sign_pass = False
    if abs(spin) > 0.0:
        frame_dragging_sign_pass = all(g < 0.0 for g in event_g_tphi)

    checks_list = [
        all_points_exterior,
        min(event_deltas) > 0.0,
        min(event_g_rr)   > 0.0,
        min(event_g_phiphi) > 0.0,
    ]
    if abs(spin) <= 0.0:
        checks_list.append(schwarzschild_reduction_pass)
    else:
        checks_list.append(frame_dragging_sign_pass)

    return {
        # ---- CSV / JSON public fields ----
        "spin_a":                    spin,
        "M":                         mass,
        "N":                         n,
        "seed":                      seed,
        "margin":                    margin,
        "r_plus":                    r_plus,
        "r_erg_eq":                  2.0 * mass,
        "r_min_observed":            r_min_observed,
        "all_points_exterior":       all_points_exterior,
        "min_Delta":                 min(event_deltas),
        "min_g_rr":                  min(event_g_rr),
        "min_g_phiphi":              min(event_g_phiphi),
        "max_abs_g_tphi":            max(abs(g) for g in event_g_tphi),
        "local_evaluated_pair_count": local_evaluated,
        "local_skipped_pair_count":   local_skipped,
        "local_max_abs_dr":          max_abs_dr,
        "local_max_abs_dphi":        max_abs_dphi,
        "local_max_abs_dt":          max_abs_dt,
        "local_timelike_count":      local_timelike,
        "local_nullish_count":       local_nullish,
        "local_spacelike_count":     local_spacelike,
        "global_true_relations":     true_relations,
        "global_false_relations":    false_relations,
        "global_undecided_pairs":    undecided_pairs,
        "schwarzschild_reduction_pass": schwarzschild_reduction_pass,
        "frame_dragging_sign_pass":  frame_dragging_sign_pass,
        "all_checks_pass":           all(checks_list),
        # ---- Internal: not written to CSV ----
        "_pair_details":             pair_details,
        "_events":                   [asdict(ev) for ev in events],
    }


# ---------------------------------------------------------------------------
# Full audit runner
# ---------------------------------------------------------------------------

def run_audit(
    n:      int              = DEFAULT_N,
    seed:   int              = DEFAULT_SEED,
    mass:   float            = DEFAULT_MASS,
    spins:  tuple[float,...] = DEFAULT_SPINS,
    margin: float            = DEFAULT_MARGIN,
) -> dict[str, Any]:
    """Run K3 for all spins and return the full payload dict."""
    if mass != 1.0:
        raise ValueError("K3 diagnostic is fixed to M=1")
    if any(abs(a) >= mass for a in spins):
        raise ValueError("K3 requires |a| < M for all spins")

    rows = [run_spin_case(n, seed, mass, spin, margin) for spin in spins]
    all_pass = all(row["all_checks_pass"] for row in rows)

    aggregate: dict[str, Any] = {
        "audit":            AUDIT_ID,
        "benchmark":        "S4-K3 Kerr equatorial local null-cone diagnostic",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "command":          COMMAND,
        "N":                n,
        "seed":             seed,
        "M":                mass,
        "spins":            list(spins),
        "margin":           margin,
        "theta":            math.pi / 2.0,
        "possible_pairs":   n * (n - 1) // 2,
        "all_checks_pass":  all_pass,
        "positive_spin_cases_all_undecided": all(
            row["global_true_relations"]  == 0
            and row["global_false_relations"] == 0
            and row["global_undecided_pairs"] == n * (n - 1) // 2
            for row in rows
            if row["spin_a"] > 0.0
        ),
        "scope_note": (
            "K3 computes equatorial BL metric coefficients and classifies "
            "small-displacement local intervals by the sign of ds². "
            "Local labels are metric-sign diagnostics ONLY — they do not "
            "imply global causal reachability, do not integrate Kerr geodesics, "
            "and do not constitute a Kerr causal solver of any kind."
        ),
    }

    return {
        "aggregate": aggregate,
        "rows":      [_public_row(r) for r in rows],
        "cases":     [_case_payload(r) for r in rows],
    }


def _public_row(row: dict[str, Any]) -> dict[str, Any]:
    """Return a copy of a row without internal underscore-prefixed keys."""
    return {k: v for k, v in row.items() if not k.startswith("_")}


def _case_payload(row: dict[str, Any]) -> dict[str, Any]:
    """Return the per-case JSON payload (events + pair details)."""
    return {
        "spin_a":          row["spin_a"],
        "all_checks_pass": row["all_checks_pass"],
        "events":          row["_events"],
        "pair_details":    row["_pair_details"],
    }


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------

def _fmt(value: object) -> str:
    if isinstance(value, bool):
        return "true" if value else "false"
    if value is None:
        return ""
    if isinstance(value, float):
        return f"{value:.12g}"
    return str(value)


def write_outputs(
    payload:    dict[str, Any],
    out_prefix: str = DEFAULT_OUT_PREFIX,
) -> tuple[Path, Path, Path]:
    csv_path  = OUT_DIR / f"{out_prefix}.csv"
    json_path = OUT_DIR / f"{out_prefix}.json"
    md_path   = OUT_DIR / f"{out_prefix}.md"
    rows      = payload["rows"]
    aggregate = payload["aggregate"]

    # CSV — one row per spin value
    with csv_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=CSV_FIELDS)
        writer.writeheader()
        for row in rows:
            writer.writerow({f: _fmt(row.get(f, "")) for f in CSV_FIELDS})

    # JSON — full payload
    json_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    # Markdown
    _write_md(md_path, rows, aggregate)

    return csv_path, json_path, md_path


def _write_md(
    md_path:   Path,
    rows:      list[dict[str, Any]],
    aggregate: dict[str, Any],
) -> None:
    lines = [
        f"# {AUDIT_ID}: Kerr Equatorial Local Null-Cone Diagnostic",
        "",
        f"Generated: {aggregate['generated_at_utc']}",
        "",
        "## What this is",
        "",
        "This is a **local metric-sign diagnostic**, not a global Kerr causal solver.",
        "",
        "It computes Boyer-Lindquist equatorial metric coefficients and evaluates",
        "the quadratic interval `ds²` at the midpoint radius for small-displacement",
        "equatorial pairs. Each evaluated pair is labelled as one of:",
        "",
        "- `timelike_local_candidate`  — `ds² < -tol`",
        "- `nullish_local_candidate`   — `|ds²| ≤ tol`",
        "- `spacelike_local_candidate` — `ds² > tol`",
        "",
        "**It does NOT:**",
        "",
        "- Establish null geodesic connectivity between the two events.",
        "- Integrate Kerr geodesics of any kind.",
        "- Decide prograde or retrograde causal relations.",
        "- Constitute a Kerr causal solver of any kind.",
        "",
        "This is only the first local consistency check before any Kerr causal inference.",
        "",
        "## Parameters",
        "",
        f"- M = {aggregate['M']}, theta = pi/2 (equatorial), spins = {aggregate['spins']}",
        f"- N = {aggregate['N']}, seed = {aggregate['seed']}, margin = {aggregate['margin']}",
        f"- Local-pair filter: |dr| ≤ {DR_LOCAL_THRESHOLD} (BL units), "
        f"|dphi| ≤ {DPHI_LOCAL_THRESHOLD} rad (~{math.degrees(DPHI_LOCAL_THRESHOLD):.1f}°)",
        "",
        "## Summary",
        "",
        "| Check | Result |",
        "|-------|--------|",
        f"| **all_checks_pass** | **{aggregate['all_checks_pass']}** |",
        f"| positive_spin_cases_all_undecided | {aggregate['positive_spin_cases_all_undecided']} |",
        "",
        "## Per-Spin Results",
        "",
        "| a | r_+ | min_Δ | min_g_rr | min_g_φφ | evaluated | timelike | nullish"
        " | spacelike | global_true | undecided | checks |",
        "|---|-----|-------|----------|----------|-----------|----------|---------|"
        "-----------|-------------|-----------|--------|",
    ]
    for r in rows:
        lines.append(
            "| {spin_a:.2f} | {r_plus:.4f} | {min_Delta:.4f} | {min_g_rr:.4f} | "
            "{min_g_phiphi:.4f} | {local_evaluated_pair_count} | "
            "{local_timelike_count} | {local_nullish_count} | {local_spacelike_count} | "
            "{global_true_relations} | {global_undecided_pairs} | {all_checks_pass} |".format(**r)
        )
    lines += [
        "",
        "## Interpretation",
        "",
        "- `evaluated` = pairs passing both displacement filters "
        f"(|dr|≤{DR_LOCAL_THRESHOLD}, |dphi|≤{DPHI_LOCAL_THRESHOLD} rad).",
        "- `global_true` = global causal assertions. For a=0 these are the Schwarzschild",
        "  control counts. For a>0 this is always 0 (all pairs remain undecided globally).",
        "- `undecided` = global causal pairs not decided by this diagnostic.",
        "- **Local labels are not global causal decisions.** A `timelike_local_candidate`",
        "  label means `ds² < 0` at the midpoint radius under the local BL metric.",
        "  It does not imply global causal reachability.",
        "",
    ]
    md_path.write_text("\n".join(lines), encoding="utf-8")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main() -> None:
    print(f"Running {AUDIT_ID}")
    payload  = run_audit()
    agg      = payload["aggregate"]
    csv_path, json_path, md_path = write_outputs(payload)
    print(f"all_checks_pass={agg['all_checks_pass']}")
    for row in payload["rows"]:
        print(
            f"  a={row['spin_a']:.2f}  r_plus={row['r_plus']:.6g}"
            f"  evaluated={row['local_evaluated_pair_count']}"
            f"  timelike={row['local_timelike_count']}"
            f"  nullish={row['local_nullish_count']}"
            f"  spacelike={row['local_spacelike_count']}"
            f"  global_true={row['global_true_relations']}"
            f"  undecided={row['global_undecided_pairs']}"
            f"  checks={row['all_checks_pass']}"
        )
    print(f"wrote {csv_path}")
    print(f"wrote {json_path}")
    print(f"wrote {md_path}")


if __name__ == "__main__":
    main()
