#!/usr/bin/env python3
"""Compare independently generated PPP-AR integer constraints in Ginan traces.

The comparison currently canonicalises the common wide-lane-like four-term
constraints produced by the receiver/signal single-difference experiment:

    (N_L1 - N_L2)_satellite_a - (N_L1 - N_L2)_satellite_b = integer.

Reference satellites may differ between runs.  Each epoch/receiver graph is
therefore reduced to all implied satellite-pair differences before comparison.
Two-term single-signal satellite differences are also canonicalised.  Any
remaining constraint form is counted but deliberately not interpreted.
Agreement is restart-consistency evidence, not proof that the common integers
are correct.
"""

from __future__ import annotations

import argparse
import json
import re
from collections import defaultdict, deque
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path


EPOCH_RE = re.compile(r"^fixAndHoldAmbiguities:\s+(.+?)\s*$")
TERM_RE = re.compile(r"([+-]\d+)\s+A\(([^,]+),(G\d{2}),([^)]+)\)")
RHS_RE = re.compile(r"=\s*([+-]?\d+(?:\.\d+)?)\s*$")


@dataclass(frozen=True)
class Constraint:
    epoch: str
    coefficients: dict[tuple[str, str, str], int]
    rhs: Fraction


@dataclass(frozen=True)
class WideLaneEdge:
    epoch: str
    receiver: str
    satellite_a: str
    satellite_b: str
    difference_a_minus_b: Fraction


@dataclass(frozen=True)
class DifferenceEdge:
    epoch: str
    receiver: str
    coordinate: str
    satellite_a: str
    satellite_b: str
    difference_a_minus_b: Fraction


def parse_trace(path: Path) -> list[Constraint]:
    constraints: list[Constraint] = []
    epoch: str | None = None
    with path.open("r", encoding="utf-8", errors="replace") as stream:
        for raw_line in stream:
            epoch_match = EPOCH_RE.match(raw_line.strip())
            if epoch_match:
                epoch = epoch_match.group(1)
                continue
            if "Applying:" not in raw_line:
                continue
            if epoch is None:
                raise ValueError(f"{path}: Applying row precedes AR epoch marker")
            rhs_match = RHS_RE.search(raw_line)
            if rhs_match is None:
                raise ValueError(f"{path}: malformed Applying row: {raw_line.rstrip()}")
            coefficients: dict[tuple[str, str, str], int] = defaultdict(int)
            for coefficient, receiver, satellite, signal in TERM_RE.findall(raw_line):
                coefficients[(receiver, satellite, signal)] += int(coefficient)
            coefficients = {
                key: value for key, value in coefficients.items() if value != 0
            }
            if not coefficients:
                raise ValueError(f"{path}: Applying row has no ambiguity terms")
            constraints.append(
                Constraint(
                    epoch=epoch,
                    coefficients=coefficients,
                    rhs=Fraction(rhs_match.group(1)),
                )
            )
    return constraints


def as_wide_lane_edge(constraint: Constraint) -> WideLaneEdge | None:
    receivers = {key[0] for key in constraint.coefficients}
    satellites = sorted({key[1] for key in constraint.coefficients})
    signals = {key[2] for key in constraint.coefficients}
    if len(receivers) != 1 or len(satellites) != 2 or signals != {"L1C", "L2W"}:
        return None
    if len(constraint.coefficients) != 4:
        return None

    wide_lane_coefficients: dict[str, int] = {}
    receiver = next(iter(receivers))
    for satellite in satellites:
        l1 = constraint.coefficients.get((receiver, satellite, "L1C"), 0)
        l2 = constraint.coefficients.get((receiver, satellite, "L2W"), 0)
        if l1 != -l2 or abs(l1) != 1:
            return None
        wide_lane_coefficients[satellite] = l1

    satellite_a, satellite_b = satellites
    coefficient_a = wide_lane_coefficients[satellite_a]
    coefficient_b = wide_lane_coefficients[satellite_b]
    if coefficient_a != -coefficient_b:
        return None

    return WideLaneEdge(
        epoch=constraint.epoch,
        receiver=receiver,
        satellite_a=satellite_a,
        satellite_b=satellite_b,
        difference_a_minus_b=constraint.rhs / coefficient_a,
    )


def as_difference_edge(constraint: Constraint) -> DifferenceEdge | None:
    wide_lane = as_wide_lane_edge(constraint)
    if wide_lane is not None:
        return DifferenceEdge(
            epoch=wide_lane.epoch,
            receiver=wide_lane.receiver,
            coordinate="WL:L1C-L2W",
            satellite_a=wide_lane.satellite_a,
            satellite_b=wide_lane.satellite_b,
            difference_a_minus_b=wide_lane.difference_a_minus_b,
        )

    if len(constraint.coefficients) != 2:
        return None
    receivers = {key[0] for key in constraint.coefficients}
    satellites = sorted({key[1] for key in constraint.coefficients})
    signals = {key[2] for key in constraint.coefficients}
    if len(receivers) != 1 or len(satellites) != 2 or len(signals) != 1:
        return None

    receiver = next(iter(receivers))
    signal = next(iter(signals))
    satellite_a, satellite_b = satellites
    coefficient_a = constraint.coefficients.get((receiver, satellite_a, signal), 0)
    coefficient_b = constraint.coefficients.get((receiver, satellite_b, signal), 0)
    if abs(coefficient_a) != 1 or coefficient_a != -coefficient_b:
        return None
    return DifferenceEdge(
        epoch=constraint.epoch,
        receiver=receiver,
        coordinate=signal,
        satellite_a=satellite_a,
        satellite_b=satellite_b,
        difference_a_minus_b=constraint.rhs / coefficient_a,
    )


def implied_differences(
    edges: list[WideLaneEdge] | list[DifferenceEdge],
) -> tuple[dict[tuple[str, str], Fraction], list[str]]:
    adjacency: dict[str, list[tuple[str, Fraction]]] = defaultdict(list)
    for edge in edges:
        adjacency[edge.satellite_a].append(
            (edge.satellite_b, edge.difference_a_minus_b)
        )
        adjacency[edge.satellite_b].append(
            (edge.satellite_a, -edge.difference_a_minus_b)
        )

    potentials: dict[str, Fraction] = {}
    components: dict[str, int] = {}
    conflicts: list[str] = []
    component = 0
    for root in sorted(adjacency):
        if root in potentials:
            continue
        component += 1
        potentials[root] = Fraction(0)
        components[root] = component
        queue = deque([root])
        while queue:
            current = queue.popleft()
            for neighbour, current_minus_neighbour in adjacency[current]:
                expected = potentials[current] - current_minus_neighbour
                if neighbour in potentials:
                    if potentials[neighbour] != expected:
                        conflicts.append(
                            f"{current}-{neighbour}: expected {expected}, "
                            f"found {potentials[neighbour]}"
                        )
                    continue
                potentials[neighbour] = expected
                components[neighbour] = component
                queue.append(neighbour)

    differences: dict[tuple[str, str], Fraction] = {}
    satellites = sorted(potentials)
    for i, satellite_a in enumerate(satellites):
        for satellite_b in satellites[i + 1 :]:
            if components[satellite_a] != components[satellite_b]:
                continue
            differences[(satellite_a, satellite_b)] = (
                potentials[satellite_a] - potentials[satellite_b]
            )
    return differences, conflicts


def index_wide_lane_edges(
    constraints: list[Constraint],
) -> tuple[dict[tuple[str, str], list[WideLaneEdge]], int]:
    indexed: dict[tuple[str, str], list[WideLaneEdge]] = defaultdict(list)
    wide_lane_count = 0
    for constraint in constraints:
        edge = as_wide_lane_edge(constraint)
        if edge is None:
            continue
        indexed[(edge.epoch, edge.receiver)].append(edge)
        wide_lane_count += 1
    return indexed, wide_lane_count


def index_difference_edges(
    constraints: list[Constraint],
) -> tuple[dict[tuple[str, str, str], list[DifferenceEdge]], int]:
    indexed: dict[tuple[str, str, str], list[DifferenceEdge]] = defaultdict(list)
    recognised_count = 0
    for constraint in constraints:
        edge = as_difference_edge(constraint)
        if edge is None:
            continue
        indexed[(edge.epoch, edge.receiver, edge.coordinate)].append(edge)
        recognised_count += 1
    return indexed, recognised_count


def compare_constraints(path_a: Path, path_b: Path) -> dict[str, object]:
    constraints_a = parse_trace(path_a)
    constraints_b = parse_trace(path_b)
    wide_lane_index_a, wide_lane_count_a = index_wide_lane_edges(constraints_a)
    wide_lane_index_b, wide_lane_count_b = index_wide_lane_edges(constraints_b)
    index_a, recognised_count_a = index_difference_edges(constraints_a)
    index_b, recognised_count_b = index_difference_edges(constraints_b)

    common_keys = sorted(set(index_a) & set(index_b))
    shared_pair_comparisons = 0
    shared_wide_lane_pair_comparisons = 0
    shared_single_signal_pair_comparisons = 0
    mismatches: list[dict[str, object]] = []
    internal_conflicts: list[dict[str, object]] = []

    for key in common_keys:
        differences_a, conflicts_a = implied_differences(index_a[key])
        differences_b, conflicts_b = implied_differences(index_b[key])
        for source, conflicts in (("a", conflicts_a), ("b", conflicts_b)):
            for conflict in conflicts:
                internal_conflicts.append(
                    {"epoch": key[0], "receiver": key[1], "source": source, "detail": conflict}
                )
        for pair in sorted(set(differences_a) & set(differences_b)):
            shared_pair_comparisons += 1
            if key[2] == "WL:L1C-L2W":
                shared_wide_lane_pair_comparisons += 1
            else:
                shared_single_signal_pair_comparisons += 1
            if differences_a[pair] == differences_b[pair]:
                continue
            mismatches.append(
                {
                    "epoch": key[0],
                    "receiver": key[1],
                    "coordinate": key[2],
                    "satellite_pair": list(pair),
                    "run_a_difference": str(differences_a[pair]),
                    "run_b_difference": str(differences_b[pair]),
                }
            )

    accepted_epochs_a = {constraint.epoch for constraint in constraints_a}
    accepted_epochs_b = {constraint.epoch for constraint in constraints_b}
    return {
        "run_a": str(path_a.resolve()),
        "run_b": str(path_b.resolve()),
        "run_a_constraint_rows": len(constraints_a),
        "run_b_constraint_rows": len(constraints_b),
        "run_a_wide_lane_rows": wide_lane_count_a,
        "run_b_wide_lane_rows": wide_lane_count_b,
        "run_a_single_signal_rows": recognised_count_a - wide_lane_count_a,
        "run_b_single_signal_rows": recognised_count_b - wide_lane_count_b,
        "run_a_unrecognised_rows": len(constraints_a) - recognised_count_a,
        "run_b_unrecognised_rows": len(constraints_b) - recognised_count_b,
        "run_a_accepted_epochs": len(accepted_epochs_a),
        "run_b_accepted_epochs": len(accepted_epochs_b),
        "overlapping_accepted_epochs": len(accepted_epochs_a & accepted_epochs_b),
        "overlapping_epoch_receiver_coordinate_groups": len(common_keys),
        "shared_implied_pair_comparisons": shared_pair_comparisons,
        "shared_implied_wide_lane_pair_comparisons": shared_wide_lane_pair_comparisons,
        "shared_implied_single_signal_pair_comparisons": shared_single_signal_pair_comparisons,
        "integer_difference_mismatch_count": len(mismatches),
        "internal_graph_conflict_count": len(internal_conflicts),
        "integer_differences_consistent_on_shared_pairs": (
            shared_pair_comparisons > 0
            and not mismatches
            and not internal_conflicts
        ),
        "mismatch_examples": mismatches[:20],
        "internal_conflict_examples": internal_conflicts[:20],
        "claim_limit": (
            "Agreement covers only shared implied L1C/L2W wide-lane and single-signal "
            "satellite-pair constraints. It is restart-consistency evidence, not an "
            "independent correct-integer certificate."
        ),
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("trace_a", type=Path)
    parser.add_argument("trace_b", type=Path)
    parser.add_argument("--json-output", type=Path)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    report = compare_constraints(args.trace_a, args.trace_b)
    payload = json.dumps(report, indent=2, ensure_ascii=False) + "\n"
    print(payload, end="")
    if args.json_output:
        args.json_output.write_text(payload, encoding="utf-8")
    return 2 if report["integer_difference_mismatch_count"] or report["internal_graph_conflict_count"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
