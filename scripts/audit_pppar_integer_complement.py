#!/usr/bin/env python3
"""Audit diagnostic-only conditional PPP-AR integer-complement probes.

The audit verifies both the trace protocol and the exact integer rank after
mapping accepted stage-two rows back to the original ambiguity states.  It
does not treat a diagnostic probe as a submitted or scientifically validated
integer constraint.
"""

from __future__ import annotations

import argparse
import json
import re
from dataclasses import dataclass, field
from fractions import Fraction
from pathlib import Path

from audit_pppar_integer_rank import exact_rank
from compare_pppar_integer_constraints import (
    Constraint,
    EPOCH_RE,
    RHS_RE,
    TERM_RE,
    as_difference_edge,
    as_wide_lane_edge,
)


DIAGNOSTIC_RE = re.compile(
    r"PPP_AR\s+INTEGER_COMPLEMENT_DIAGNOSTIC\s+"
    r"stage1_rows=(\d+)\s+remaining_coordinates=(\d+)\s+"
    r"stage2_probe_rows=(\d+)\s+combined_independent_rows=(\d+)\s+"
    r"target_integer_rank=(\d+)\s+stage2_status=(\S+)\s+"
    r"stage2_success_rate=(\S+)\s+stage2_minimum_decorrelated=(\d+)\s+"
    r"stage2_selected_decorrelated=(\d+)\s+"
    r"stage2_integer_candidates=(\d+)\s+stage2_best_squared_norm=(\S+)\s+"
    r"stage2_second_squared_norm=(\S+)\s+stage2_ratio=(\S+)\s+"
    r"action=(\S+)"
)
ROW_RE = re.compile(
    r"PPP_AR\s+INTEGER_COMPLEMENT_ROW\s+stage=2\s+row=(\d+)\s+"
    r"rhs=(\S+)\s+support=(\d+)\s+status=(\S+)\s+"
    r"terms=(.*?)action=(\S+)"
)


@dataclass
class ComplementDiagnostic:
    epoch: str
    stage1_rows: int
    remaining_coordinates: int
    stage2_probe_rows: int
    combined_independent_rows: int
    target_integer_rank: int
    stage2_status: str
    stage2_success_rate: float
    stage2_minimum_decorrelated: int
    stage2_selected_decorrelated: int
    stage2_integer_candidates: int
    stage2_best_squared_norm: float
    stage2_second_squared_norm: float
    stage2_ratio: float
    action: str
    stage1_constraints: list[Constraint] = field(default_factory=list)
    stage2_constraints: list[Constraint] = field(default_factory=list)


def _parse_coefficients(text: str) -> dict[tuple[str, str, str], int]:
    coefficients: dict[tuple[str, str, str], int] = {}
    for coefficient, receiver, satellite, signal in TERM_RE.findall(text):
        key = (receiver, satellite, signal)
        coefficients[key] = coefficients.get(key, 0) + int(coefficient)
    return {key: value for key, value in coefficients.items() if value != 0}


def parse_trace(path: Path) -> list[ComplementDiagnostic]:
    epoch: str | None = None
    diagnostics: list[ComplementDiagnostic] = []
    by_epoch: dict[str, ComplementDiagnostic] = {}
    stage1_by_epoch: dict[str, list[Constraint]] = {}

    with path.open("r", encoding="utf-8", errors="replace") as stream:
        for raw_line in stream:
            epoch_match = EPOCH_RE.match(raw_line.strip())
            if epoch_match:
                epoch = epoch_match.group(1)
                stage1_by_epoch.setdefault(epoch, [])
                continue
            if epoch is None:
                continue

            if "Applying:" in raw_line:
                rhs_match = RHS_RE.search(raw_line)
                coefficients = _parse_coefficients(raw_line)
                if rhs_match is None or not coefficients:
                    raise ValueError(f"{path}: malformed stage-one row: {raw_line.strip()}")
                stage1_by_epoch[epoch].append(
                    Constraint(epoch, coefficients, Fraction(rhs_match.group(1)))
                )
                continue

            diagnostic_match = DIAGNOSTIC_RE.search(raw_line)
            if diagnostic_match:
                if epoch in by_epoch:
                    raise ValueError(f"{path}: duplicate complement diagnostic at {epoch}")
                values = diagnostic_match.groups()
                diagnostic = ComplementDiagnostic(
                    epoch=epoch,
                    stage1_rows=int(values[0]),
                    remaining_coordinates=int(values[1]),
                    stage2_probe_rows=int(values[2]),
                    combined_independent_rows=int(values[3]),
                    target_integer_rank=int(values[4]),
                    stage2_status=values[5],
                    stage2_success_rate=float(values[6]),
                    stage2_minimum_decorrelated=int(values[7]),
                    stage2_selected_decorrelated=int(values[8]),
                    stage2_integer_candidates=int(values[9]),
                    stage2_best_squared_norm=float(values[10]),
                    stage2_second_squared_norm=float(values[11]),
                    stage2_ratio=float(values[12]),
                    action=values[13],
                )
                diagnostics.append(diagnostic)
                by_epoch[epoch] = diagnostic
                continue

            row_match = ROW_RE.search(raw_line)
            if row_match:
                if epoch not in by_epoch:
                    raise ValueError(f"{path}: complement row precedes diagnostic at {epoch}")
                row_index, rhs, support, status, terms, action = row_match.groups()
                diagnostic = by_epoch[epoch]
                if int(row_index) != len(diagnostic.stage2_constraints):
                    raise ValueError(f"{path}: non-contiguous complement row index at {epoch}")
                if status != "INTEGER_MAPPED" or action != "PROBE_ONLY_NOT_SUBMITTED":
                    raise ValueError(f"{path}: invalid complement row status/action at {epoch}")
                coefficients = _parse_coefficients(terms)
                if len(coefficients) != int(support) or not coefficients:
                    raise ValueError(f"{path}: complement support mismatch at {epoch}")
                diagnostic.stage2_constraints.append(
                    Constraint(epoch, coefficients, Fraction(rhs))
                )

    for diagnostic in diagnostics:
        diagnostic.stage1_constraints = stage1_by_epoch.get(diagnostic.epoch, [])
    return diagnostics


def audit_trace(
    path: Path,
    success_threshold: float = 0.9999,
    ratio_threshold: float = 3.0,
) -> dict[str, object]:
    diagnostics = parse_trace(path)
    accepted = 0
    insufficient = 0
    stage2_wide_lane = 0
    stage2_single_signal = 0
    stage2_other = 0
    details: list[dict[str, object]] = []

    for item in diagnostics:
        if item.action != "PROBE_ONLY_NOT_SUBMITTED":
            raise ValueError(f"{path}: diagnostic action is not probe-only at {item.epoch}")
        if item.remaining_coordinates != item.target_integer_rank - item.stage1_rows:
            raise ValueError(f"{path}: remaining-coordinate mismatch at {item.epoch}")
        if item.combined_independent_rows != item.stage1_rows + item.stage2_probe_rows:
            raise ValueError(f"{path}: combined-row arithmetic mismatch at {item.epoch}")
        if item.combined_independent_rows > item.target_integer_rank:
            raise ValueError(f"{path}: combined rows exceed target rank at {item.epoch}")
        if len(item.stage1_constraints) != item.stage1_rows:
            raise ValueError(f"{path}: stage-one row count mismatch at {item.epoch}")
        if len(item.stage2_constraints) != item.stage2_probe_rows:
            raise ValueError(f"{path}: stage-two row count mismatch at {item.epoch}")

        stage1_rank = exact_rank([row.coefficients for row in item.stage1_constraints])
        stage2_rank = exact_rank([row.coefficients for row in item.stage2_constraints])
        combined_rank = exact_rank(
            [row.coefficients for row in item.stage1_constraints + item.stage2_constraints]
        )
        if stage1_rank != item.stage1_rows or stage2_rank != item.stage2_probe_rows:
            raise ValueError(f"{path}: dependent rows within a stage at {item.epoch}")
        if combined_rank != item.combined_independent_rows:
            raise ValueError(f"{path}: stage-two rows are not independent of stage one at {item.epoch}")

        is_accepted = item.stage2_status == "RESOLVED_RATIO_ACCEPTED"
        if is_accepted:
            accepted += 1
            if item.stage2_probe_rows <= 0:
                raise ValueError(f"{path}: accepted stage two has no rows at {item.epoch}")
            if item.stage2_success_rate < success_threshold:
                raise ValueError(f"{path}: accepted stage two is below success threshold at {item.epoch}")
            if item.stage2_ratio < ratio_threshold or item.stage2_integer_candidates < 2:
                raise ValueError(f"{path}: accepted stage two is below ratio requirements at {item.epoch}")
        elif item.stage2_status == "INSUFFICIENT_DECORRELATED_AMBIGUITIES":
            insufficient += 1
            if (
                item.stage2_probe_rows != 0
                or item.stage2_selected_decorrelated
                >= item.stage2_minimum_decorrelated
            ):
                raise ValueError(f"{path}: inconsistent insufficient-dimension result at {item.epoch}")
        elif item.stage2_probe_rows != 0:
            raise ValueError(f"{path}: rejected stage two reports fixed rows at {item.epoch}")

        wide_lane = sum(as_wide_lane_edge(row) is not None for row in item.stage2_constraints)
        single_signal = sum(
            as_wide_lane_edge(row) is None and as_difference_edge(row) is not None
            for row in item.stage2_constraints
        )
        other = item.stage2_probe_rows - wide_lane - single_signal
        stage2_wide_lane += wide_lane
        stage2_single_signal += single_signal
        stage2_other += other
        details.append(
            {
                "epoch": item.epoch,
                "stage1_rows": item.stage1_rows,
                "remaining_coordinates": item.remaining_coordinates,
                "stage2_probe_rows": item.stage2_probe_rows,
                "combined_independent_rows": item.combined_independent_rows,
                "target_integer_rank": item.target_integer_rank,
                "stage2_status": item.stage2_status,
                "stage2_success_rate": item.stage2_success_rate,
                "stage2_minimum_decorrelated": item.stage2_minimum_decorrelated,
                "stage2_selected_decorrelated": item.stage2_selected_decorrelated,
                "stage2_integer_candidates": item.stage2_integer_candidates,
                "stage2_best_squared_norm": item.stage2_best_squared_norm,
                "stage2_second_squared_norm": item.stage2_second_squared_norm,
                "stage2_ratio": item.stage2_ratio,
                "stage1_exact_rank": stage1_rank,
                "stage2_exact_rank": stage2_rank,
                "combined_exact_rank": combined_rank,
                "stage2_wide_lane_rows": wide_lane,
                "stage2_single_signal_rows": single_signal,
                "stage2_other_rows": other,
                "action": item.action,
            }
        )

    return {
        "trace": str(path.resolve()),
        "scope": "diagnostic_only_conditional_integer_complement_exact_rank",
        "success_threshold": success_threshold,
        "ratio_threshold": ratio_threshold,
        "diagnostic_epoch_count": len(diagnostics),
        "accepted_probe_epoch_count": accepted,
        "insufficient_dimension_epoch_count": insufficient,
        "stage2_probe_row_count": sum(item.stage2_probe_rows for item in diagnostics),
        "stage2_wide_lane_row_count": stage2_wide_lane,
        "stage2_single_signal_row_count": stage2_single_signal,
        "stage2_other_row_count": stage2_other,
        "all_actions_probe_only_not_submitted": all(
            item.action == "PROBE_ONLY_NOT_SUBMITTED" for item in diagnostics
        ),
        "epochs": details,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("trace", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--success-threshold", type=float, default=0.9999)
    parser.add_argument("--ratio-threshold", type=float, default=3.0)
    args = parser.parse_args()

    report = audit_trace(args.trace, args.success_threshold, args.ratio_threshold)
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(rendered, encoding="utf-8")
    print(rendered, end="")


if __name__ == "__main__":
    main()
