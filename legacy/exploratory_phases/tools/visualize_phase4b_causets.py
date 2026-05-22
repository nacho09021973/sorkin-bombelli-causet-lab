#!/usr/bin/env python3
"""Generate Phase 4B causal-set visual audit SVGs.

This is a provenance/visual-audit tool, not a new simulation phase.  It
selects near-mean, best, and worst seeds from the persisted Phase 4B
per-seed CSV, reproduces the corresponding causets deterministically,
and writes Hasse/projection SVGs via ``cones.write_causet_svg``.

The causet depends on (n, target_dim, seed).  Epsilon is used only to
select which seed/loss rows to inspect.
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import cones  # noqa: E402
import validation_suite as vs  # noqa: E402


FOUNDATION = ROOT / "benchmarks" / "foundation"
AGGREGATE_CSV = FOUNDATION / "phase4b_survival_probe.csv"
PER_EPSILON_CSV = FOUNDATION / "phase4b_survival_probe_per_epsilon.csv"
PER_SEED_CSV = FOUNDATION / "phase4b_survival_probe_per_seed.csv"
VISUAL_DIR = FOUNDATION / "phase4b_causet_svgs"
VISUAL_INDEX_CSV = FOUNDATION / "phase4b_causet_visual_index.csv"
VISUAL_LINK_AUDIT_CSV = FOUNDATION / "phase4b_visual_link_audit.csv"
INTRINSIC_POSET_AUDIT_CSV = FOUNDATION / "phase4b_intrinsic_poset_audit.csv"
PHASE4B_MD = FOUNDATION / "phase4b_survival_probe.md"

PRIORITY_CELLS: tuple[tuple[int, int], ...] = (
    (32, 4),
    (48, 3),
    (48, 4),
    (64, 4),
    (32, 2),
)

VISUAL_INDEX_HEADERS = (
    "n",
    "target_dim",
    "epsilon",
    "role",
    "seed",
    "loss",
    "mean_loss_cell_epsilon",
    "abs_delta_from_mean",
    "svg_path",
    "curve_shape_cell",
    "survival_label_cell",
    "borderline_v_like_cell",
)

VISUAL_LINK_AUDIT_HEADERS = (
    "n",
    "target_dim",
    "epsilon",
    "role",
    "seed",
    "loss",
    "n_links_Hasse",
    "curve_shape_cell",
    "survival_label_cell",
    "borderline_v_like_cell",
)

INTRINSIC_POSET_AUDIT_HEADERS = (
    "n",
    "target_dim",
    "epsilon",
    "role",
    "seed",
    "loss",
    "ordering_fraction",
    "chain3_abundance",
    "n_links_Hasse",
    "dim_discrepancy_rel_midpoint",
    "curve_shape_cell",
    "survival_label_cell",
    "borderline_v_like_cell",
)


def _read_csv(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def _fmt_float_for_filename(value: float) -> str:
    return f"{value:g}".replace(".", "p").replace("-", "m")


def _selected_epsilon_by_cell(aggregate_rows: list[dict]) -> dict[tuple[int, int], float]:
    out: dict[tuple[int, int], float] = {}
    for row in aggregate_rows:
        cell = (int(row["n"]), int(row["target_dim"]))
        if cell in PRIORITY_CELLS:
            out[cell] = float(row["epsilon_at_min"])
    return out


def select_visual_rows(
    per_seed_rows: list[dict],
    aggregate_rows: list[dict],
    priority_cells: tuple[tuple[int, int], ...] = PRIORITY_CELLS,
) -> list[dict]:
    """Select near-mean, best, and worst valid seeds for each priority cell."""
    eps_by_cell = _selected_epsilon_by_cell(aggregate_rows)
    selected: list[dict] = []

    for cell in priority_cells:
        if cell not in eps_by_cell:
            continue
        n, d = cell
        eps = eps_by_cell[cell]
        rows = [
            r for r in per_seed_rows
            if int(r["n"]) == n
            and int(r["target_dim"]) == d
            and math.isclose(float(r["epsilon"]), eps, rel_tol=0.0, abs_tol=1e-12)
            and r["valid"].lower() == "true"
        ]
        if not rows:
            continue

        mean_loss = float(rows[0]["mean_loss_cell_epsilon"])
        role_rows = {
            "near_mean": min(rows, key=lambda r: abs(float(r["loss"]) - mean_loss)),
            "best": min(rows, key=lambda r: float(r["loss"])),
            "worst": max(rows, key=lambda r: float(r["loss"])),
        }
        for role in ("near_mean", "best", "worst"):
            row = dict(role_rows[role])
            loss = float(row["loss"])
            row["role"] = role
            row["abs_delta_from_mean"] = abs(loss - mean_loss)
            selected.append(row)
    return selected


def _svg_filename(row: dict) -> str:
    eps = _fmt_float_for_filename(float(row["epsilon"]))
    return (
        f"phase4b_n{int(row['n'])}_d{int(row['target_dim'])}"
        f"_eps{eps}_seed{int(row['seed'])}_{row['role']}.svg"
    )


def generate_svgs(selected_rows: list[dict], out_dir: Path = VISUAL_DIR) -> list[dict]:
    out_dir.mkdir(parents=True, exist_ok=True)
    index_rows: list[dict] = []
    for row in selected_rows:
        n = int(row["n"])
        d = int(row["target_dim"])
        seed = int(row["seed"])
        matrix, points = vs.sprinkle_minkowski_diamond(
            n=n,
            seed=seed,
            d_spacetime=d,
        )
        svg_path = out_dir / _svg_filename(row)
        cones.write_causet_svg(points, svg_path, matrix)
        try:
            stored_svg_path = str(svg_path.relative_to(ROOT))
        except ValueError:
            stored_svg_path = str(svg_path)
        index_rows.append({
            "n": n,
            "target_dim": d,
            "epsilon": float(row["epsilon"]),
            "role": row["role"],
            "seed": seed,
            "loss": float(row["loss"]),
            "mean_loss_cell_epsilon": float(row["mean_loss_cell_epsilon"]),
            "abs_delta_from_mean": float(row["abs_delta_from_mean"]),
            "svg_path": stored_svg_path,
            "curve_shape_cell": row["curve_shape_cell"],
            "survival_label_cell": row["survival_label_cell"],
            "borderline_v_like_cell": row["borderline_v_like_cell"],
        })
    return index_rows


def build_link_audit_rows(selected_rows: list[dict]) -> list[dict]:
    """Count Hasse links for the selected Phase 4B visual audit causets."""
    out: list[dict] = []
    for row in selected_rows:
        n = int(row["n"])
        d = int(row["target_dim"])
        seed = int(row["seed"])
        matrix, _ = vs.sprinkle_minkowski_diamond(
            n=n,
            seed=seed,
            d_spacetime=d,
        )
        out.append({
            "n": n,
            "target_dim": d,
            "epsilon": float(row["epsilon"]),
            "role": row["role"],
            "seed": seed,
            "loss": float(row["loss"]),
            "n_links_Hasse": len(cones.transitive_reduction(matrix)),
            "curve_shape_cell": row["curve_shape_cell"],
            "survival_label_cell": row["survival_label_cell"],
            "borderline_v_like_cell": row["borderline_v_like_cell"],
        })
    return out


def build_intrinsic_poset_audit_rows(
    selected_rows: list[dict],
    link_rows: list[dict],
) -> list[dict]:
    """Join selected seed provenance to Hasse link counts."""
    links_by_key = {
        (
            int(r["n"]),
            int(r["target_dim"]),
            float(r["epsilon"]),
            str(r["role"]),
            int(r["seed"]),
        ): r
        for r in link_rows
    }
    out: list[dict] = []
    for row in selected_rows:
        key = (
            int(row["n"]),
            int(row["target_dim"]),
            float(row["epsilon"]),
            str(row["role"]),
            int(row["seed"]),
        )
        link_row = links_by_key[key]
        out.append({
            "n": int(row["n"]),
            "target_dim": int(row["target_dim"]),
            "epsilon": float(row["epsilon"]),
            "role": row["role"],
            "seed": int(row["seed"]),
            "loss": float(row["loss"]),
            "ordering_fraction": float(row["ordering_fraction"]),
            "chain3_abundance": float(row["chain3_abundance"]),
            "n_links_Hasse": int(link_row["n_links_Hasse"]),
            "dim_discrepancy_rel_midpoint": float(row["dim_discrepancy_rel_midpoint"]),
            "curve_shape_cell": row["curve_shape_cell"],
            "survival_label_cell": row["survival_label_cell"],
            "borderline_v_like_cell": row["borderline_v_like_cell"],
        })
    return out


def write_visual_index(index_rows: list[dict], path: Path = VISUAL_INDEX_CSV) -> None:
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=VISUAL_INDEX_HEADERS)
        writer.writeheader()
        for row in index_rows:
            writer.writerow(row)


def write_link_audit_csv(
    link_rows: list[dict],
    path: Path = VISUAL_LINK_AUDIT_CSV,
) -> None:
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=VISUAL_LINK_AUDIT_HEADERS)
        writer.writeheader()
        for row in link_rows:
            writer.writerow(row)


def write_intrinsic_poset_audit_csv(
    intrinsic_rows: list[dict],
    path: Path = INTRINSIC_POSET_AUDIT_CSV,
) -> None:
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=INTRINSIC_POSET_AUDIT_HEADERS)
        writer.writeheader()
        for row in intrinsic_rows:
            writer.writerow(row)


def _link_audit_markdown_table(link_rows: list[dict]) -> list[str]:
    lines = [
        "| n | target_dim | epsilon | role | seed | loss | n_links_Hasse | curve_shape | survival_label | borderline |",
        "| ---: | :---: | ---: | --- | ---: | ---: | ---: | --- | --- | :---: |",
    ]
    for row in sorted(
        link_rows,
        key=lambda r: (int(r["n"]), int(r["target_dim"]), str(r["role"])),
    ):
        lines.append(
            f"| {row['n']} | {row['target_dim']} | {float(row['epsilon']):.4g} "
            f"| {row['role']} | {row['seed']} | {float(row['loss']):.4g} "
            f"| {row['n_links_Hasse']} | {row['curve_shape_cell']} "
            f"| {row['survival_label_cell']} | {row['borderline_v_like_cell']} |"
        )
    return lines


def _intrinsic_poset_markdown_table(intrinsic_rows: list[dict]) -> list[str]:
    lines = [
        "| n | target_dim | role | seed | loss | ordering_fraction | chain3_abundance | n_links_Hasse | dim_disc_rel |",
        "| ---: | :---: | --- | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in sorted(
        intrinsic_rows,
        key=lambda r: (int(r["n"]), int(r["target_dim"]), str(r["role"])),
    ):
        lines.append(
            f"| {row['n']} | {row['target_dim']} | {row['role']} "
            f"| {row['seed']} | {float(row['loss']):.4g} "
            f"| {float(row['ordering_fraction']):.4g} "
            f"| {float(row['chain3_abundance']):.4g} "
            f"| {row['n_links_Hasse']} "
            f"| {float(row['dim_discrepancy_rel_midpoint']):.4g} |"
        )
    return lines


def update_markdown_visual_section(
    md_path: Path = PHASE4B_MD,
    link_rows: list[dict] | None = None,
    intrinsic_rows: list[dict] | None = None,
) -> None:
    if not md_path.exists():
        return
    text = md_path.read_text(encoding="utf-8")
    marker = "## Visual audit provenance"
    section = "\n".join([
        marker,
        "",
        "The SVG files in `phase4b_causet_svgs/` are Hasse/projection diagrams: edges are links from the transitive reduction, not all transitive relations.",
        "",
        "The SVG projection uses the `(t, x1)` plane, so for `target_dim=3` or `target_dim=4` it does not represent the full spatial geometry.",
        "",
        "Seeds are selected from `phase4b_survival_probe_per_seed.csv` as `near_mean`, `best`, and `worst` at each priority cell's `epsilon_at_min`; they are not arbitrary examples.",
        "",
        "This is an exploratory visual audit only and is not confirmatory evidence.",
        "",
    ])
    link_marker = "## Visual link audit"
    link_section = ""
    if link_rows is not None:
        link_section = "\n".join([
            link_marker,
            "",
            "`phase4b_visual_link_audit.csv` tabulates `loss` against `n_links_Hasse` for the selected visual audit seeds. This is exploratory and does not infer causality.",
            "",
            *_link_audit_markdown_table(link_rows),
            "",
        ])
    intrinsic_marker = "## Intrinsic poset audit"
    intrinsic_section = ""
    if intrinsic_rows is not None:
        intrinsic_section = "\n".join([
            intrinsic_marker,
            "",
            "`phase4b_intrinsic_poset_audit.csv` compares `loss` against intrinsic partial-order observables for the selected visual audit seeds, not against the SVG projection.",
            "",
            "No robust scalar predictor is established: there is no stable monotonic relationship between `loss` and any single observable among `ordering_fraction`, `n_links_Hasse`, or `chain3_abundance`.",
            "",
            "Some cells show local alignments, especially `(48,3)`, but other visualized cells contradict a simple one-variable explanation.",
            "",
            "The global exploratory outcome remains **MIXED**, and no physical causality is inferred from this audit.",
            "",
            *_intrinsic_poset_markdown_table(intrinsic_rows),
            "",
        ])
    if marker in text:
        before = text.split(marker, 1)[0].rstrip()
        after = text.split(marker, 1)[1]
        next_marker = after.find("\n## ")
        if next_marker >= 0:
            rest = after[next_marker:].lstrip("\n")
            new_text = before + "\n\n" + section + rest
        else:
            new_text = before + "\n\n" + section
    else:
        insert_before = "\n## Global exploratory outcome"
        if insert_before in text:
            new_text = text.replace(insert_before, "\n" + section + insert_before, 1)
        else:
            new_text = text.rstrip() + "\n\n" + section
    if link_section:
        if link_marker in new_text:
            before = new_text.split(link_marker, 1)[0].rstrip()
            after = new_text.split(link_marker, 1)[1]
            next_marker = after.find("\n## ")
            if next_marker >= 0:
                rest = after[next_marker:].lstrip("\n")
                new_text = before + "\n\n" + link_section + rest
            else:
                new_text = before + "\n\n" + link_section
        else:
            insert_before = "\n## Global exploratory outcome"
            if insert_before in new_text:
                new_text = new_text.replace(
                    insert_before,
                    "\n" + link_section + insert_before,
                    1,
                )
            else:
                new_text = new_text.rstrip() + "\n\n" + link_section
    if intrinsic_section:
        if intrinsic_marker in new_text:
            before = new_text.split(intrinsic_marker, 1)[0].rstrip()
            after = new_text.split(intrinsic_marker, 1)[1]
            next_marker = after.find("\n## ")
            if next_marker >= 0:
                rest = after[next_marker:].lstrip("\n")
                new_text = before + "\n\n" + intrinsic_section + rest
            else:
                new_text = before + "\n\n" + intrinsic_section
        else:
            insert_before = "\n## Global exploratory outcome"
            if insert_before in new_text:
                new_text = new_text.replace(
                    insert_before,
                    "\n" + intrinsic_section + insert_before,
                    1,
                )
            else:
                new_text = new_text.rstrip() + "\n\n" + intrinsic_section
    md_path.write_text(new_text, encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=VISUAL_DIR)
    parser.add_argument("--index-csv", type=Path, default=VISUAL_INDEX_CSV)
    parser.add_argument("--link-audit-csv", type=Path, default=VISUAL_LINK_AUDIT_CSV)
    parser.add_argument(
        "--intrinsic-audit-csv",
        type=Path,
        default=INTRINSIC_POSET_AUDIT_CSV,
    )
    args = parser.parse_args()

    aggregate_rows = _read_csv(AGGREGATE_CSV)
    # Read per-epsilon to fail early if provenance is incomplete; selection
    # itself is driven by per-seed rows joined to epsilon_at_min.
    _read_csv(PER_EPSILON_CSV)
    per_seed_rows = _read_csv(PER_SEED_CSV)

    selected = select_visual_rows(per_seed_rows, aggregate_rows)
    index_rows = generate_svgs(selected, args.output_dir)
    link_rows = build_link_audit_rows(selected)
    intrinsic_rows = build_intrinsic_poset_audit_rows(selected, link_rows)
    write_visual_index(index_rows, args.index_csv)
    write_link_audit_csv(link_rows, args.link_audit_csv)
    write_intrinsic_poset_audit_csv(intrinsic_rows, args.intrinsic_audit_csv)
    update_markdown_visual_section(PHASE4B_MD, link_rows, intrinsic_rows)

    print(f"Selected rows: {len(selected)}")
    print(f"SVGs written:  {len(index_rows)}")
    print(f"Index CSV:     {args.index_csv}")
    print(f"Link audit CSV: {args.link_audit_csv}")
    print(f"Intrinsic audit CSV: {args.intrinsic_audit_csv}")


if __name__ == "__main__":
    main()
