#!/usr/bin/env python3
"""Audit the exact rank and signal structure of submitted PPP-AR rows.

The trace reports constraints in the original ambiguity-state coordinates.
This tool uses exact rational Gaussian elimination, so rank conclusions do not
depend on a floating-point tolerance.  It deliberately audits each epoch on
its own; rows are not accumulated across epochs because ambiguity arcs may
change after slips or state resets.
"""

from __future__ import annotations

import argparse
import json
import re
from collections import defaultdict
from dataclasses import dataclass, field
from fractions import Fraction
from pathlib import Path
from statistics import median

from compare_pppar_integer_constraints import (
    Constraint,
    EPOCH_RE,
    RHS_RE,
    TERM_RE,
    as_difference_edge,
    as_wide_lane_edge,
)


INTEGER_COORDINATES_RE = re.compile(
    r"PPP_AR\s+INTEGER_COORDINATES\s+original=(\d+)\s+integer=(\d+)"
)
RECEIVER_SD_RE = re.compile(
    r"PPP_AR\s+RECEIVER_SD\s+receiver=(\S+)\s+system=(\S+)\s+"
    r"signal=(\S+)\s+reference=(\S+)\s+members=(\d+)\s+status=(\S+)"
)
SUBMISSION_RE = re.compile(
    r"PPP_AR\s+PSEUDOOBS_SUBMISSION\s+rows=(\d+)\s+status=(\S+)"
)


@dataclass
class EpochAudit:
    epoch: str
    integer_coordinate_count: int | None = None
    original_ambiguity_count: int | None = None
    target_by_receiver: dict[str, int] = field(default_factory=lambda: defaultdict(int))
    constraints: list[Constraint] = field(default_factory=list)
    submitted_rows: int | None = None
    submission_status: str | None = None


def exact_rank(rows: list[dict[tuple[str, str, str], int]]) -> int:
    """Return exact row rank for sparse integer rows."""
    columns = sorted({column for row in rows for column in row})
    if not rows or not columns:
        return 0
    matrix = [
        [Fraction(row.get(column, 0)) for column in columns]
        for row in rows
    ]
    rank = 0
    for column in range(len(columns)):
        pivot = next(
            (row for row in range(rank, len(matrix)) if matrix[row][column] != 0),
            None,
        )
        if pivot is None:
            continue
        matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
        pivot_value = matrix[rank][column]
        matrix[rank] = [value / pivot_value for value in matrix[rank]]
        for row in range(len(matrix)):
            if row == rank or matrix[row][column] == 0:
                continue
            factor = matrix[row][column]
            matrix[row] = [
                value - factor * pivot_entry
                for value, pivot_entry in zip(matrix[row], matrix[rank])
            ]
        rank += 1
        if rank == len(matrix):
            break
    return rank


def parse_trace(path: Path) -> list[EpochAudit]:
    epochs: list[EpochAudit] = []
    current: EpochAudit | None = None
    with path.open("r", encoding="utf-8", errors="replace") as stream:
        for raw_line in stream:
            stripped = raw_line.strip()
            epoch_match = EPOCH_RE.match(stripped)
            if epoch_match:
                current = EpochAudit(epoch=epoch_match.group(1))
                epochs.append(current)
                continue
            if current is None:
                continue

            coordinates_match = INTEGER_COORDINATES_RE.search(raw_line)
            if coordinates_match:
                current.original_ambiguity_count = int(coordinates_match.group(1))
                current.integer_coordinate_count = int(coordinates_match.group(2))
                continue

            receiver_match = RECEIVER_SD_RE.search(raw_line)
            if receiver_match:
                receiver, _system, _signal, _reference, members, status = (
                    receiver_match.groups()
                )
                if status == "INTEGER_COORDINATES_CREATED":
                    current.target_by_receiver[receiver] += max(0, int(members) - 1)
                continue

            if "Applying:" in raw_line:
                rhs_match = RHS_RE.search(raw_line)
                if rhs_match is None:
                    raise ValueError(f"{path}: malformed Applying row: {stripped}")
                coefficients: dict[tuple[str, str, str], int] = defaultdict(int)
                for coefficient, receiver, satellite, signal in TERM_RE.findall(raw_line):
                    coefficients[(receiver, satellite, signal)] += int(coefficient)
                coefficients = {
                    key: value for key, value in coefficients.items() if value != 0
                }
                if not coefficients:
                    raise ValueError(f"{path}: Applying row has no ambiguity terms")
                current.constraints.append(
                    Constraint(
                        epoch=current.epoch,
                        coefficients=coefficients,
                        rhs=Fraction(rhs_match.group(1)),
                    )
                )
                continue

            submission_match = SUBMISSION_RE.search(raw_line)
            if submission_match:
                current.submitted_rows = int(submission_match.group(1))
                current.submission_status = submission_match.group(2)

    return epochs


def audit_trace(path: Path) -> dict[str, object]:
    epochs = parse_trace(path)
    submitted = [epoch for epoch in epochs if epoch.submitted_rows is not None]
    details: list[dict[str, object]] = []
    total_wide_lane = 0
    total_single_signal = 0
    total_other = 0
    total_rank = 0
    duplicate_rows = 0
    complete_rank_epochs = 0
    deficits: list[int] = []

    for epoch in submitted:
        if epoch.submitted_rows != len(epoch.constraints):
            raise ValueError(
                f"{path}: {epoch.epoch}: submission reports {epoch.submitted_rows} rows "
                f"but {len(epoch.constraints)} Applying rows were parsed"
            )
        rows = [constraint.coefficients for constraint in epoch.constraints]
        rank = exact_rank(rows)
        wide_lane_rows = [
            constraint.coefficients
            for constraint in epoch.constraints
            if as_wide_lane_edge(constraint) is not None
        ]
        single_signal_rows = [
            constraint.coefficients
            for constraint in epoch.constraints
            if as_wide_lane_edge(constraint) is None
            and as_difference_edge(constraint) is not None
        ]
        other_count = len(rows) - len(wide_lane_rows) - len(single_signal_rows)
        target = epoch.integer_coordinate_count
        deficit = None if target is None else target - rank
        if deficit is not None:
            deficits.append(deficit)
            complete_rank_epochs += int(deficit == 0)

        receiver_rows: dict[str, list[dict[tuple[str, str, str], int]]] = defaultdict(list)
        cross_receiver_rows = 0
        for row in rows:
            receivers = {key[0] for key in row}
            if len(receivers) == 1:
                receiver_rows[next(iter(receivers))].append(row)
            else:
                cross_receiver_rows += 1
        receiver_details = {
            receiver: {
                "target_rank": epoch.target_by_receiver.get(receiver),
                "submitted_rows": len(receiver_rows.get(receiver, [])),
                "submitted_rank": exact_rank(receiver_rows.get(receiver, [])),
            }
            for receiver in sorted(set(epoch.target_by_receiver) | set(receiver_rows))
        }

        total_wide_lane += len(wide_lane_rows)
        total_single_signal += len(single_signal_rows)
        total_other += other_count
        total_rank += rank
        duplicate_rows += len(rows) - rank
        details.append(
            {
                "epoch": epoch.epoch,
                "original_ambiguity_count": epoch.original_ambiguity_count,
                "target_integer_rank": target,
                "submitted_rows": len(rows),
                "submitted_rank": rank,
                "rank_deficit": deficit,
                "wide_lane_rows": len(wide_lane_rows),
                "wide_lane_rank": exact_rank(wide_lane_rows),
                "single_signal_rows": len(single_signal_rows),
                "single_signal_rank": exact_rank(single_signal_rows),
                "other_rows": other_count,
                "cross_receiver_rows": cross_receiver_rows,
                "receiver_components": receiver_details,
                "submission_status": epoch.submission_status,
            }
        )

    return {
        "trace": str(path.resolve()),
        "scope": "per_epoch_exact_integer_row_rank_no_cross_epoch_accumulation",
        "ar_attempt_epoch_count": len(epochs),
        "submitted_epoch_count": len(submitted),
        "submitted_row_count": sum(len(epoch.constraints) for epoch in submitted),
        "submitted_rank_sum": total_rank,
        "duplicate_or_dependent_row_count": duplicate_rows,
        "wide_lane_row_count": total_wide_lane,
        "single_signal_row_count": total_single_signal,
        "other_row_count": total_other,
        "complete_rank_epoch_count": complete_rank_epochs,
        "rank_deficit_min": min(deficits) if deficits else None,
        "rank_deficit_median": median(deficits) if deficits else None,
        "rank_deficit_max": max(deficits) if deficits else None,
        "epochs": details,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("trace", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = audit_trace(args.trace)
    payload = json.dumps(report, indent=2, sort_keys=True)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(payload + "\n", encoding="utf-8")
    print(payload)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
