# S4-KERR-K3-LOCAL-CONE-001 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add `audit_kerr_k3_local_cone_001.py` — the first genuinely Kerr-geometry diagnostic — that computes equatorial Boyer-Lindquist metric coefficients and classifies small-displacement local intervals as timelike/nullish/spacelike, while preserving the K1/K2 invariant that for a>0 no global causal relations are decided.

**Architecture:** One new audit module in `explore/sorkin4_kerr_benchmark/`; imports helpers from the existing `run_kerr_minimal_benchmark` without modifying it; one test file in `tests/`; three permanent artifacts (CSV/JSON/MD). README gets a single K3 section appended.

**Tech Stack:** Python 3.12, stdlib only (math, csv, json, random, dataclasses). Imports `run_kerr_minimal_benchmark as kerr` and `run_schwarzschild_minimal_benchmark as schwarz`.

---

## Task 1: Write the failing tests

**Files:**
- Create: `tests/test_sorkin4_kerr_k3_local_cone.py`

- [ ] **Step 1.1 — Write the test file**

```python
# tests/test_sorkin4_kerr_k3_local_cone.py
"""Tests for S4-KERR-K3-LOCAL-CONE-001.

Covers:
1.  Artifact files exist and are parseable.
2.  all_checks_pass is True in the JSON aggregate.
3.  a=0 schwarzschild_reduction_pass is True.
4.  For a=0, max_abs_g_tphi == 0 (within tolerance).
5.  For a>0, frame_dragging_sign_pass is True.
6.  For all spins, all_points_exterior is True.
7.  For all spins, min_Delta > 0, min_g_rr > 0, min_g_phiphi > 0.
8.  For a>0, global_true_relations=0, global_false_relations=0,
    global_undecided_pairs = N*(N-1)/2.
9.  local_timelike_count + local_nullish_count + local_spacelike_count
    == local_evaluated_pair_count (NOT necessarily N*(N-1)/2).
10. Existing K1/K2 tests: the test runner covers them in the validation
    command; no explicit re-test here to avoid duplication.
"""
from __future__ import annotations

import json
import math
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "explore" / "sorkin4_kerr_benchmark"))

import audit_kerr_k3_local_cone_001 as audit

ARTIFACT_DIR = ROOT / "explore" / "sorkin4_kerr_benchmark"
ARTIFACT_JSON = ARTIFACT_DIR / "kerr_k3_local_cone_001_n12_seed1959.json"
ARTIFACT_CSV  = ARTIFACT_DIR / "kerr_k3_local_cone_001_n12_seed1959.csv"
ARTIFACT_MD   = ARTIFACT_DIR / "kerr_k3_local_cone_001_n12_seed1959.md"

METRIC_TOL = 1.0e-12


class KerrK3ArtifactTests(unittest.TestCase):
    """Guard the committed frozen artifacts."""

    def setUp(self) -> None:
        for path in (ARTIFACT_JSON, ARTIFACT_CSV, ARTIFACT_MD):
            self.assertTrue(
                path.exists(),
                msg=(
                    f"Missing K3 artifact {path}. Run: "
                    "python explore/sorkin4_kerr_benchmark/"
                    "audit_kerr_k3_local_cone_001.py"
                ),
            )
        with ARTIFACT_JSON.open(encoding="utf-8") as fh:
            self.payload = json.load(fh)

    # Test 1 — files exist & parseable: covered by setUp

    # Test 2 — aggregate all_checks_pass
    def test_all_checks_pass(self) -> None:
        self.assertTrue(self.payload["aggregate"]["all_checks_pass"])

    # Tests 3-9 on frozen rows
    def test_per_row_checks(self) -> None:
        rows = self.payload["rows"]
        n = self.payload["aggregate"]["N"]
        possible_pairs = n * (n - 1) // 2

        for row in rows:
            spin = row["spin_a"]

            # 6: all_points_exterior
            self.assertTrue(
                row["all_points_exterior"],
                msg=f"all_points_exterior failed for a={spin}",
            )
            # 7: metric positivity
            self.assertGreater(row["min_Delta"], 0.0,   msg=f"min_Delta ≤ 0 for a={spin}")
            self.assertGreater(row["min_g_rr"], 0.0,    msg=f"min_g_rr ≤ 0 for a={spin}")
            self.assertGreater(row["min_g_phiphi"], 0.0,msg=f"min_g_phiphi ≤ 0 for a={spin}")

            # 9: local count integrity
            local_sum = (
                row["local_timelike_count"]
                + row["local_nullish_count"]
                + row["local_spacelike_count"]
            )
            self.assertEqual(
                local_sum,
                row["local_evaluated_pair_count"],
                msg=f"local count mismatch for a={spin}: {local_sum} != {row['local_evaluated_pair_count']}",
            )

            if abs(spin) <= 0.0:
                # Test 3: Schwarzschild reduction
                self.assertTrue(
                    row["schwarzschild_reduction_pass"],
                    msg="schwarzschild_reduction_pass failed for a=0",
                )
                # Test 4: g_tphi = 0 for a=0
                self.assertLessEqual(
                    row["max_abs_g_tphi"],
                    METRIC_TOL,
                    msg=f"max_abs_g_tphi={row['max_abs_g_tphi']} > tol for a=0",
                )
            else:
                # Test 5: frame dragging sign for a>0
                self.assertTrue(
                    row["frame_dragging_sign_pass"],
                    msg=f"frame_dragging_sign_pass failed for a={spin}",
                )
                # Test 8: global causal accounting for a>0
                self.assertEqual(row["global_true_relations"], 0,   msg=f"a={spin}")
                self.assertEqual(row["global_false_relations"], 0,  msg=f"a={spin}")
                self.assertEqual(
                    row["global_undecided_pairs"],
                    possible_pairs,
                    msg=f"a={spin}",
                )


class KerrK3ComputationTests(unittest.TestCase):
    """In-process tests: run computation without relying on artifact files."""

    def setUp(self) -> None:
        self.payload = audit.run_audit()
        self.rows = self.payload["rows"]
        self.aggregate = self.payload["aggregate"]

    def test_aggregate_all_checks_pass(self) -> None:
        self.assertTrue(self.aggregate["all_checks_pass"])

    def test_spins_and_n(self) -> None:
        self.assertEqual(self.aggregate["spins"], [0.0, 0.25, 0.5, 0.75])
        self.assertEqual(self.aggregate["M"], 1.0)
        self.assertEqual(self.aggregate["N"], 12)
        self.assertEqual(len(self.rows), 4)

    def test_per_row_local_count_integrity(self) -> None:
        for row in self.rows:
            local_sum = (
                row["local_timelike_count"]
                + row["local_nullish_count"]
                + row["local_spacelike_count"]
            )
            self.assertEqual(local_sum, row["local_evaluated_pair_count"])

    def test_a0_schwarzschild_reduction(self) -> None:
        a0_row = next(r for r in self.rows if r["spin_a"] == 0.0)
        self.assertTrue(a0_row["schwarzschild_reduction_pass"])
        self.assertLessEqual(a0_row["max_abs_g_tphi"], METRIC_TOL)

    def test_positive_spin_frame_dragging(self) -> None:
        for row in self.rows:
            if row["spin_a"] > 0.0:
                self.assertTrue(row["frame_dragging_sign_pass"], msg=f"a={row['spin_a']}")

    def test_positive_spin_global_undecided(self) -> None:
        n = self.aggregate["N"]
        possible = n * (n - 1) // 2
        for row in self.rows:
            if row["spin_a"] > 0.0:
                self.assertEqual(row["global_true_relations"], 0)
                self.assertEqual(row["global_false_relations"], 0)
                self.assertEqual(row["global_undecided_pairs"], possible)

    def test_equatorial_metric_at_r_a0(self) -> None:
        """a=0 must give exact Schwarzschild equatorial BL coefficients."""
        r, M = 4.0, 1.0
        m = audit.equatorial_metric_at_r(r, M, 0.0)
        self.assertAlmostEqual(m["g_tt"],     -(1.0 - 2*M/r), places=14)
        self.assertAlmostEqual(m["g_tphi"],   0.0,            places=14)
        self.assertAlmostEqual(m["g_rr"],     1.0/(1.0-2*M/r),places=14)
        self.assertAlmostEqual(m["g_phiphi"], r*r,            places=14)

    def test_equatorial_metric_at_r_frame_dragging_sign(self) -> None:
        """g_tphi < 0 for a>0, M>0 at any exterior r."""
        for a in (0.25, 0.5, 0.75):
            m = audit.equatorial_metric_at_r(4.0, 1.0, a)
            self.assertLess(m["g_tphi"], 0.0, msg=f"g_tphi >= 0 for a={a}")

    def test_local_pair_filter(self) -> None:
        """Only pairs with |dr| <= DR and |dphi| <= DPHI are evaluated."""
        from explore.sorkin4_kerr_benchmark.run_kerr_minimal_benchmark import Event
        import math
        # Two events with dr=0.2, dphi=0.1 — both within threshold
        p = Event(index=0, t=0.0, r=3.0, theta=math.pi/2, phi=0.0)
        q = Event(index=1, t=1.0, r=3.2, theta=math.pi/2, phi=0.1)
        self.assertTrue(audit.is_local_pair(p, q))
        # dr=2.0 exceeds DR_LOCAL_THRESHOLD=1.0
        q_far = Event(index=2, t=1.0, r=5.5, theta=math.pi/2, phi=0.1)
        self.assertFalse(audit.is_local_pair(p, q_far))


if __name__ == "__main__":
    unittest.main()
```

- [ ] **Step 1.2 — Run the test, verify it fails with ImportError**

```bash
cd /home/ignac/sorkin
python3 -m unittest tests.test_sorkin4_kerr_k3_local_cone 2>&1 | head -10
```
Expected: `ModuleNotFoundError: No module named 'audit_kerr_k3_local_cone_001'`

---

## Task 2: Implement the audit module

**Files:**
- Create: `explore/sorkin4_kerr_benchmark/audit_kerr_k3_local_cone_001.py`

- [ ] **Step 2.1 — Write the module**

```python
#!/usr/bin/env python3
"""S4-KERR-K3-LOCAL-CONE-001: Kerr equatorial local null-cone diagnostic.

This is the first genuine Kerr-geometry diagnostic in the SORKIN-4 program.
It computes Boyer-Lindquist metric coefficients at sampled equatorial exterior
points and classifies small-displacement local intervals by the sign of ds².

BOUNDARY (read before extending this module):
  local_timelike_candidate  ≠  causal_relation = True
  local_spacelike_candidate ≠  causal_relation = False

The labels are local metric-sign diagnostics ONLY.  They do NOT:
  - imply global causal reachability between the two events,
  - integrate Kerr null geodesics,
  - decide prograde or retrograde causal relations,
  - constitute a Kerr causal solver of any kind.

For a=0, the Schwarzschild reduction is verified and the a=0 control counts
from K1/K2 are preserved.  For a>0, all global causal pairs remain undecided.
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

AUDIT_ID = "S4-KERR-K3-LOCAL-CONE-001"
OUT_DIR = Path(__file__).resolve().parent
DEFAULT_OUT_PREFIX = "kerr_k3_local_cone_001_n12_seed1959"
COMMAND = (
    "python3 explore/sorkin4_kerr_benchmark/audit_kerr_k3_local_cone_001.py"
)

# ---------------------------------------------------------------------------
# Sweep parameters
# ---------------------------------------------------------------------------

DEFAULT_N      = 12
DEFAULT_SEED   = 1959
DEFAULT_MASS   = 1.0
DEFAULT_MARGIN = schwarz.EXTERIOR_MARGIN          # 0.35, same as K1/K2
DEFAULT_SPINS  = (0.0, 0.25, 0.5, 0.75)

# ---------------------------------------------------------------------------
# Local-pair filter thresholds (Adjustment 1)
# ---------------------------------------------------------------------------

DR_LOCAL_THRESHOLD   = 1.0   # |dr| in BL coords
DPHI_LOCAL_THRESHOLD = 0.5   # |dphi| in radians (~28.6°)

# ---------------------------------------------------------------------------
# Numerical tolerances
# ---------------------------------------------------------------------------

METRIC_TOL    = 1.0e-12   # for Schwarzschild-reduction checks
LOCAL_DS2_EPS = 1.0e-9    # scale factor for ds² classifier

# ---------------------------------------------------------------------------
# CSV schema
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

    Equatorial simplifications (theta=pi/2, sin=1, cos=0):
      Sigma = r²
      Delta = r² - 2Mr + a²
      g_tt     = -(1 - 2M/r)
      g_tphi   = -2Ma/r
      g_rr     = r²/Delta
      g_phiphi = r² + a² + 2Ma²/r

    Raises ValueError if Delta <= 0 (point inside or on horizon).
    """
    delta = r * r - 2.0 * mass * r + spin * spin
    if delta <= 0.0:
        raise ValueError(
            f"equatorial_metric_at_r: Delta={delta:.6g} <= 0 at r={r:.6g}, "
            f"M={mass}, a={spin}; point is inside or on the Kerr horizon."
        )
    return {
        "g_tt":     -(1.0 - 2.0 * mass / r),
        "g_tphi":   -2.0 * mass * spin / r,
        "g_rr":     r * r / delta,
        "g_phiphi": r * r + spin * spin + 2.0 * mass * spin * spin / r,
        "delta":    delta,
    }


def is_local_pair(p: Event, q: Event) -> bool:
    """Return True if (p, q) is within the local-displacement thresholds."""
    dr   = abs(q.r - p.r)
    dphi = abs(signed_delta_phi(p.phi, q.phi))
    return dr <= DR_LOCAL_THRESHOLD and dphi <= DPHI_LOCAL_THRESHOLD


def local_interval_ds2_equatorial(
    p: Event,
    q: Event,
    mass: float,
    spin: float,
) -> tuple[float, dict[str, float]]:
    """Compute ds² at the midpoint radius for an equatorial pair.

    ds² = g_tt dt² + 2 g_tphi dt dphi + g_rr dr² + g_phiphi dphi²

    dtheta = 0 by construction (all K3 events are at theta=pi/2).

    Returns (ds2, metric_at_midpoint).
    """
    r_mid = 0.5 * (p.r + q.r)
    metric = equatorial_metric_at_r(r_mid, mass, spin)
    dt    = q.t   - p.t
    dr    = q.r   - p.r
    dphi  = signed_delta_phi(p.phi, q.phi)
    ds2 = (
        metric["g_tt"]     * dt   * dt
        + 2.0 * metric["g_tphi"] * dt   * dphi
        + metric["g_rr"]   * dr   * dr
        + metric["g_phiphi"] * dphi * dphi
    )
    return ds2, metric


def classify_local_interval(
    ds2: float,
    scale: float,
    eps_factor: float = LOCAL_DS2_EPS,
) -> str:
    """Classify a local ds² value as timelike / nullish / spacelike.

    The tolerance is scale-adaptive: tol = eps_factor * max(1, |scale|).

    Returns one of:
      'timelike_local_candidate'   — ds² < -tol
      'nullish_local_candidate'    — |ds²| <= tol
      'spacelike_local_candidate'  — ds² >  tol
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
    n: int,
    seed: int,
    mass: float,
    spin: float,
    margin: float,
) -> dict[str, Any]:
    """Run one (spin, N, seed, margin) K3 case and return a result dict."""
    r_plus = kerr_horizon_radius(mass, spin)
    r_min  = r_plus + margin

    # Generate equatorial events (theta=pi/2 for all)
    events = generate_exterior_events(n, seed, r_min, equatorial=True)

    # --- Adjustment 2: compute Delta explicitly for each event ---
    event_deltas:   list[float] = []
    event_g_rr:     list[float] = []
    event_g_phiphi: list[float] = []
    event_g_tphi:   list[float] = []

    for ev in events:
        delta = ev.r * ev.r - 2.0 * mass * ev.r + spin * spin
        event_deltas.append(delta)
        m = equatorial_metric_at_r(ev.r, mass, spin)
        event_g_rr.append(m["g_rr"])
        event_g_phiphi.append(m["g_phiphi"])
        event_g_tphi.append(m["g_tphi"])

    # --- Adjustment 3: causal accounting per K1/K2 rules ---
    # For a=0: use Schwarzschild control (equatorial_scaffold reproduces it).
    # For a>0: build_relation_states with "equatorial_scaffold" leaves all
    #          pairs undecided by construction.
    matrix, states = build_relation_states(events, mass, spin, "equatorial_scaffold")
    possible_pairs  = n * (n - 1) // 2
    true_relations  = count_true_relations(matrix)
    false_relations = sum(
        1 for i in range(n - 1) for j in range(i + 1, n) if states[i][j] is False
    )
    undecided_pairs = sum(
        1 for i in range(n - 1) for j in range(i + 1, n) if states[i][j] is None
    )

    # --- Local cone classification (Adjustment 1: small-displacement only) ---
    local_evaluated   = 0
    local_skipped     = 0
    local_timelike    = 0
    local_nullish     = 0
    local_spacelike   = 0
    max_abs_dr        = 0.0
    max_abs_dphi      = 0.0
    max_abs_dt        = 0.0
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

            # Scale for adaptive tolerance
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
                "i": i,
                "j": j,
                "dr": dr_val,
                "dphi": dphi_val,
                "dt": dt_val,
                "r_mid": 0.5 * (p.r + q.r),
                "ds2": ds2,
                "label": label,
            })

    # --- Metric checks ---
    r_min_observed     = min(ev.r for ev in events)
    all_points_exterior = r_min_observed > r_min

    # Schwarzschild reduction (a=0 only)
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

    # Frame-dragging sign (a>0 only): g_tphi = -2Ma/r < 0
    frame_dragging_sign_pass = False
    if abs(spin) > 0.0:
        frame_dragging_sign_pass = all(g < 0.0 for g in event_g_tphi)

    # Assemble per-spin checks
    checks_list = [
        all_points_exterior,
        min(event_deltas) > 0.0,
        min(event_g_rr) > 0.0,
        min(event_g_phiphi) > 0.0,
    ]
    if abs(spin) <= 0.0:
        checks_list.append(schwarzschild_reduction_pass)
    else:
        checks_list.append(frame_dragging_sign_pass)

    row: dict[str, Any] = {
        "spin_a":                   spin,
        "M":                        mass,
        "N":                        n,
        "seed":                     seed,
        "margin":                   margin,
        "r_plus":                   r_plus,
        "r_erg_eq":                 2.0 * mass,
        "r_min_observed":           r_min_observed,
        "all_points_exterior":      all_points_exterior,
        "min_Delta":                min(event_deltas),
        "min_g_rr":                 min(event_g_rr),
        "min_g_phiphi":             min(event_g_phiphi),
        "max_abs_g_tphi":           max(abs(g) for g in event_g_tphi),
        "local_evaluated_pair_count": local_evaluated,
        "local_skipped_pair_count":   local_skipped,
        "local_max_abs_dr":         max_abs_dr,
        "local_max_abs_dphi":       max_abs_dphi,
        "local_max_abs_dt":         max_abs_dt,
        "local_timelike_count":     local_timelike,
        "local_nullish_count":      local_nullish,
        "local_spacelike_count":    local_spacelike,
        "global_true_relations":    true_relations,
        "global_false_relations":   false_relations,
        "global_undecided_pairs":   undecided_pairs,
        "schwarzschild_reduction_pass": schwarzschild_reduction_pass,
        "frame_dragging_sign_pass": frame_dragging_sign_pass,
        "all_checks_pass":          all(checks_list),
        # Internal: not written to CSV
        "_pair_details":            pair_details,
        "_events":                  [asdict(ev) for ev in events],
    }
    return row


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
    """Run K3 for all spins and return the full payload."""
    if mass != 1.0:
        raise ValueError("K3 is fixed to M=1")
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
            "and do not constitute a Kerr causal solver."
        ),
    }

    return {
        "aggregate": aggregate,
        "rows":      [_public_row(r) for r in rows],
        "cases":     [_case_payload(r) for r in rows],
    }


def _public_row(row: dict[str, Any]) -> dict[str, Any]:
    """Strip internal keys from a row for CSV/JSON rows list."""
    return {k: v for k, v in row.items() if not k.startswith("_")}


def _case_payload(row: dict[str, Any]) -> dict[str, Any]:
    return {
        "spin_a":       row["spin_a"],
        "all_checks_pass": row["all_checks_pass"],
        "events":       row["_events"],
        "pair_details": row["_pair_details"],
    }


# ---------------------------------------------------------------------------
# I/O
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

    # CSV — one row per spin
    with csv_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=CSV_FIELDS)
        writer.writeheader()
        for row in rows:
            writer.writerow({f: _fmt(row.get(f, "")) for f in CSV_FIELDS})

    # JSON
    json_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    # MD
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
        "the quadratic interval `ds²` for small-displacement equatorial pairs.",
        "Each pair is labelled as `timelike_local_candidate`, `nullish_local_candidate`,",
        "or `spacelike_local_candidate`.",
        "",
        "**It does NOT:**",
        "",
        "- Establish null geodesic connectivity between the two events.",
        "- Integrate Kerr geodesics of any kind.",
        "- Decide prograde or retrograde causal relations.",
        "- Constitute a Kerr causal solver of any kind.",
        "",
        "This is the first local consistency check before any Kerr causal inference.",
        "",
        "## Parameters",
        "",
        f"- M = {aggregate['M']}, theta = pi/2, spins = {aggregate['spins']}",
        f"- N = {aggregate['N']}, seed = {aggregate['seed']}, margin = {aggregate['margin']}",
        f"- Local-pair filter: |dr| ≤ {DR_LOCAL_THRESHOLD}, |dphi| ≤ {DPHI_LOCAL_THRESHOLD} rad",
        "",
        "## Results",
        "",
        f"| Check | Result |",
        f"|-------|--------|",
        f"| **all_checks_pass** | **{aggregate['all_checks_pass']}** |",
        f"| positive_spin_cases_all_undecided | {aggregate['positive_spin_cases_all_undecided']} |",
        "",
        "## Per-Spin Table",
        "",
        "| a | r_+ | min_Δ | min_g_rr | min_g_φφ | evaluated | timelike | nullish | spacelike | true | undecided | checks |",
        "|---|-----|-------|----------|----------|-----------|----------|---------|-----------|------|-----------|--------|",
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
        "- `evaluated` = pairs passing the local-displacement filter "
        f"(|dr|≤{DR_LOCAL_THRESHOLD}, |dphi|≤{DPHI_LOCAL_THRESHOLD} rad).",
        "- `true` = global causal assertions (a=0 Schwarzschild control only; 0 for a>0).",
        "- `undecided` = global causal pairs not decided by this diagnostic.",
        "- Local labels are NOT global causal decisions.",
        "",
    ]
    md_path.write_text("\n".join(lines), encoding="utf-8")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main() -> None:
    print(f"Running {AUDIT_ID}")
    payload   = run_audit()
    agg       = payload["aggregate"]
    csv_path, json_path, md_path = write_outputs(payload)
    print(f"all_checks_pass={agg['all_checks_pass']}")
    for row in payload["rows"]:
        print(
            f"  a={row['spin_a']:.2f} r_plus={row['r_plus']:.6g} "
            f"evaluated={row['local_evaluated_pair_count']} "
            f"timelike={row['local_timelike_count']} "
            f"nullish={row['local_nullish_count']} "
            f"spacelike={row['local_spacelike_count']} "
            f"global_true={row['global_true_relations']} "
            f"undecided={row['global_undecided_pairs']} "
            f"checks={row['all_checks_pass']}"
        )
    print(f"wrote {csv_path}")
    print(f"wrote {json_path}")
    print(f"wrote {md_path}")


if __name__ == "__main__":
    main()
```

- [ ] **Step 2.2 — Run the computation tests (artifact tests will still fail)**

```bash
cd /home/ignac/sorkin
python3 -m unittest tests.test_sorkin4_kerr_k3_local_cone.KerrK3ComputationTests -v 2>&1
```
Expected: all KerrK3ComputationTests pass; KerrK3ArtifactTests skipped/not run.

---

## Task 3: Generate the permanent artifacts

**Files:**
- Generate: `explore/sorkin4_kerr_benchmark/kerr_k3_local_cone_001_n12_seed1959.{csv,json,md}`

- [ ] **Step 3.1 — Run the audit script**

```bash
cd /home/ignac/sorkin
python3 explore/sorkin4_kerr_benchmark/audit_kerr_k3_local_cone_001.py
```
Expected output ends with:
```
all_checks_pass=True
  a=0.00 ...
  a=0.25 ...
  a=0.50 ...
  a=0.75 ...
wrote .../kerr_k3_local_cone_001_n12_seed1959.csv
wrote .../kerr_k3_local_cone_001_n12_seed1959.json
wrote .../kerr_k3_local_cone_001_n12_seed1959.md
```

- [ ] **Step 3.2 — Run full test class including artifact tests**

```bash
cd /home/ignac/sorkin
python3 -m unittest tests.test_sorkin4_kerr_k3_local_cone -v 2>&1
```
Expected: all tests pass (0 failures).

---

## Task 4: Update the Kerr benchmark README

**Files:**
- Modify: `explore/sorkin4_kerr_benchmark/README.md`

- [ ] **Step 4.1 — Append the K3 section**

Open `explore/sorkin4_kerr_benchmark/README.md` and append after the K2 block:

```markdown
The K3 local cone diagnostic is `audit_kerr_k3_local_cone_001.py`.
It is still conservative: no Kerr causal relations are decided.
It fixes `theta = pi/2` and introduces the first genuine Kerr-geometry
computation: the Boyer-Lindquist equatorial metric coefficients.

K3 computes for each sub-extremal spin in `(0.0, 0.25, 0.5, 0.75)`:
- `Delta`, `g_tt`, `g_tphi`, `g_rr`, `g_phiphi` at each event point.
- For small-displacement equatorial pairs (`|dr| ≤ 1.0`, `|dphi| ≤ 0.5 rad`),
  the local quadratic interval `ds²` at the midpoint radius.
- Classification of each evaluated pair as `timelike_local_candidate`,
  `nullish_local_candidate`, or `spacelike_local_candidate`.

K3 freezes these controls:

- `a=0.0` Schwarzschild reduction: `g_tphi = 0`, `g_rr = 1/(1-2M/r)`,
  `g_phiphi = r²`, exactly.
- `a>0` frame-dragging sign: `g_tphi = -2Ma/r < 0`.
- `a>0` global causal accounting: all unordered pairs remain undecided
  with true `0`, false `0`, undecided `N*(N-1)/2`.
- Local labels are **metric-sign diagnostics only**.  They do not imply
  null geodesic connectivity or any global causal reachability.

The K3 artifact is
`kerr_k3_local_cone_001_n12_seed1959.{csv,json,md}`.
```

---

## Task 5: Full validation run and commit

- [ ] **Step 5.1 — Run the full validation command from the spec**

```bash
cd /home/ignac/sorkin
python3 -m unittest \
  tests.test_sorkin4_kerr_k3_local_cone \
  tests.test_sorkin4_kerr_k2_equatorial_diagnostic \
  tests.test_sorkin4_kerr_k1_control \
  tests.test_sorkin4_kerr_l0_scaffold \
  tests.test_sorkin4_schw_stability_001 \
  tests.test_sorkin4_exterior_turning_asymptotic_audit \
  -v 2>&1
```
Expected: all tests pass, 0 failures.

- [ ] **Step 5.2 — Verify git diff shows only the intended new files**

```bash
cd /home/ignac/sorkin
git status
git diff HEAD --stat
```
Expected: only these paths appear:
```
explore/sorkin4_kerr_benchmark/audit_kerr_k3_local_cone_001.py  (new)
explore/sorkin4_kerr_benchmark/kerr_k3_local_cone_001_n12_seed1959.csv  (new)
explore/sorkin4_kerr_benchmark/kerr_k3_local_cone_001_n12_seed1959.json  (new)
explore/sorkin4_kerr_benchmark/kerr_k3_local_cone_001_n12_seed1959.md  (new)
explore/sorkin4_kerr_benchmark/README.md  (modified)
tests/test_sorkin4_kerr_k3_local_cone.py  (new)
docs/superpowers/specs/2026-05-27-kerr-k3-local-cone-design.md  (new)
docs/superpowers/plans/2026-05-27-kerr-k3-local-cone.md  (new)
```
No changes to `run_schwarzschild_minimal_benchmark.py`, `run_kerr_minimal_benchmark.py`,
`cones.py`, or any existing benchmark CSV/JSON.

- [ ] **Step 5.3 — Commit**

```bash
cd /home/ignac/sorkin
git add \
  explore/sorkin4_kerr_benchmark/audit_kerr_k3_local_cone_001.py \
  explore/sorkin4_kerr_benchmark/kerr_k3_local_cone_001_n12_seed1959.csv \
  explore/sorkin4_kerr_benchmark/kerr_k3_local_cone_001_n12_seed1959.json \
  explore/sorkin4_kerr_benchmark/kerr_k3_local_cone_001_n12_seed1959.md \
  explore/sorkin4_kerr_benchmark/README.md \
  tests/test_sorkin4_kerr_k3_local_cone.py \
  docs/superpowers/specs/2026-05-27-kerr-k3-local-cone-design.md \
  docs/superpowers/plans/2026-05-27-kerr-k3-local-cone.md
git commit -m "sorkin4: add Kerr K3 local cone diagnostic

Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>"
```

---

## Self-Review

**Spec coverage:**
- ✅ BL equatorial metric implemented via `equatorial_metric_at_r`
- ✅ Delta computed explicitly before `kerr_metric_components` call (Adj. 2)
- ✅ Local-pair filter with |dr|/|dphi| thresholds + count fields (Adj. 1)
- ✅ a=0 Schwarzschild control preserved in causal accounting (Adj. 3)
- ✅ a>0 all pairs undecided globally
- ✅ `schwarzschild_reduction_pass` checks g_tphi=0, g_rr=1/(1-2M/r), g_phiphi=r²
- ✅ `frame_dragging_sign_pass` checks g_tphi < 0 for a>0
- ✅ All 27 CSV fields from spec present
- ✅ MD states the four "does NOT" boundaries explicitly
- ✅ All 10 test requirements covered across the two test classes
- ✅ README updated conservatively
- ✅ Commit message matches spec

**Type/name consistency:**
- `equatorial_metric_at_r` called consistently in `run_spin_case` and tests
- `is_local_pair` called in `run_spin_case` and tests
- `DR_LOCAL_THRESHOLD`, `DPHI_LOCAL_THRESHOLD` exported at module level (used in tests and MD)
- `_public_row` strips `_pair_details` and `_events` — CSV fields don't reference them
- `CSV_FIELDS` tuple matches all `row` keys in `_public_row`

**Placeholder scan:** No TBD/TODO/fill-in-later in any step. ✅
