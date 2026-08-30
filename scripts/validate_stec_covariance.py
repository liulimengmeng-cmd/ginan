#!/usr/bin/env python3
"""Validate Ginan GINAN_STEC_COVARIANCE_V2/V3 epoch blocks.

The validator treats every META/STATE/COV block as an independently auditable
posterior covariance epoch.  It does not infer successful ambiguity fixing from
the posterior_stage field; that field records only where in the processing
pipeline the covariance was captured.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Iterator, TextIO

import numpy as np


@dataclass
class StateRow:
    local_index: int
    site: str
    satellite: str
    state_number: int
    filter_index: int
    estimate_tecu: float
    variance_tecu2: float


@dataclass
class EpochBlock:
    schema: str
    gps_week: int
    gps_tow_text: str
    status: str
    state_count: int
    upper_triangle_count: int
    writer_max_abs_asymmetry_tecu2: float
    posterior_stage: str
    ar_routine_invoked: bool
    ar_eligible_ambiguity_count: int
    ar_integer_ambiguity_coordinate_count: int
    ar_receiver_single_difference_applied: bool
    ar_receiver_datum_group_count: int
    ar_dropped_singleton_group_count: int
    ar_resolved_combination_count: int
    ar_pseudoobservations_submitted: bool
    ar_mode: str
    ar_configured_success_rate_threshold: float
    ar_configured_solution_ratio_threshold: float
    ar_diagnostic_status: str
    ar_selected_decorrelated_ambiguity_count: int
    ar_integer_candidate_count: int
    ar_bootstrapped_success_rate: float
    ar_best_squared_norm: float
    ar_second_squared_norm: float
    ar_solution_ratio: float
    states: dict[int, StateRow] = field(default_factory=dict)
    covariance: dict[tuple[int, int], float] = field(default_factory=dict)


def _require_length(row: list[str], expected: int, line_number: int) -> None:
    if len(row) != expected:
        raise ValueError(
            f"line {line_number}: {row[0] if row else 'EMPTY'} expects "
            f"{expected} fields, got {len(row)}"
        )


def _parse_binary_flag(value: str, field_name: str, line_number: int) -> bool:
    if value not in {"0", "1"}:
        raise ValueError(f"line {line_number}: {field_name} expects 0 or 1, got {value!r}")
    return value == "1"


def iter_epoch_blocks(stream: TextIO) -> Iterator[EpochBlock]:
    schema: str | None = None
    current: EpochBlock | None = None

    for line_number, raw_line in enumerate(stream, start=1):
        stripped = raw_line.strip()
        if stripped in {
            "# GINAN_STEC_COVARIANCE_V2",
            "# GINAN_STEC_COVARIANCE_V3",
        }:
            detected = stripped.removeprefix("# ")
            if schema is not None and schema != detected:
                raise ValueError(
                    f"line {line_number}: mixed covariance schemas {schema} and {detected}"
                )
            schema = detected
            continue
        if not raw_line.strip() or raw_line.startswith("#"):
            continue

        row = next(csv.reader([raw_line]))
        record_type = row[0]

        if record_type == "META":
            if schema is None:
                raise ValueError(f"line {line_number}: META precedes covariance schema marker")
            if schema == "GINAN_STEC_COVARIANCE_V2":
                _require_length(row, 22, line_number)
                integer_coordinate_count = int(row[9])
                receiver_single_difference_applied = False
                receiver_datum_group_count = 0
                dropped_singleton_group_count = 0
                resolved_index = 10
            else:
                _require_length(row, 26, line_number)
                integer_coordinate_count = int(row[10])
                receiver_single_difference_applied = _parse_binary_flag(
                    row[11], "ar_receiver_single_difference_applied", line_number
                )
                receiver_datum_group_count = int(row[12])
                dropped_singleton_group_count = int(row[13])
                resolved_index = 14
            if current is not None:
                yield current
            current = EpochBlock(
                schema=schema,
                gps_week=int(row[1]),
                gps_tow_text=row[2],
                status=row[3],
                state_count=int(row[4]),
                upper_triangle_count=int(row[5]),
                writer_max_abs_asymmetry_tecu2=float(row[6]),
                posterior_stage=row[7],
                ar_routine_invoked=_parse_binary_flag(
                    row[8], "ar_routine_invoked", line_number
                ),
                ar_eligible_ambiguity_count=int(row[9]),
                ar_integer_ambiguity_coordinate_count=integer_coordinate_count,
                ar_receiver_single_difference_applied=receiver_single_difference_applied,
                ar_receiver_datum_group_count=receiver_datum_group_count,
                ar_dropped_singleton_group_count=dropped_singleton_group_count,
                ar_resolved_combination_count=int(row[resolved_index]),
                ar_pseudoobservations_submitted=_parse_binary_flag(
                    row[resolved_index + 1], "ar_pseudoobservations_submitted", line_number
                ),
                ar_mode=row[resolved_index + 2],
                ar_configured_success_rate_threshold=float(row[resolved_index + 3]),
                ar_configured_solution_ratio_threshold=float(row[resolved_index + 4]),
                ar_diagnostic_status=row[resolved_index + 5],
                ar_selected_decorrelated_ambiguity_count=int(row[resolved_index + 6]),
                ar_integer_candidate_count=int(row[resolved_index + 7]),
                ar_bootstrapped_success_rate=float(row[resolved_index + 8]),
                ar_best_squared_norm=float(row[resolved_index + 9]),
                ar_second_squared_norm=float(row[resolved_index + 10]),
                ar_solution_ratio=float(row[resolved_index + 11]),
            )
            continue

        if record_type == "STATE":
            _require_length(row, 10, line_number)
        elif record_type == "COV":
            _require_length(row, 6, line_number)
        else:
            raise ValueError(f"line {line_number}: unknown record type {record_type!r}")

        if current is None:
            raise ValueError(f"line {line_number}: {record_type} precedes META")

        if int(row[1]) != current.gps_week or row[2] != current.gps_tow_text:
            raise ValueError(f"line {line_number}: epoch key differs from current META")

        if record_type == "STATE":
            state = StateRow(
                local_index=int(row[3]),
                site=row[4],
                satellite=row[5],
                state_number=int(row[6]),
                filter_index=int(row[7]),
                estimate_tecu=float(row[8]),
                variance_tecu2=float(row[9]),
            )
            if state.local_index in current.states:
                raise ValueError(f"line {line_number}: duplicate STATE index {state.local_index}")
            current.states[state.local_index] = state
        elif record_type == "COV":
            pair = (int(row[3]), int(row[4]))
            if pair in current.covariance:
                raise ValueError(f"line {line_number}: duplicate COV pair {pair}")
            current.covariance[pair] = float(row[5])

    if schema is None:
        raise ValueError("missing GINAN_STEC_COVARIANCE_V2/V3 schema marker")
    if current is not None:
        yield current


def validate_epoch(
    epoch: EpochBlock,
    psd_absolute_tolerance: float = 1e-12,
    psd_relative_tolerance: float = 1e-10,
) -> dict[str, object]:
    errors: list[str] = []
    n = epoch.state_count
    expected_entries = n * (n + 1) // 2

    if epoch.status != "OK":
        errors.append(f"writer_status={epoch.status}")
    if epoch.ar_eligible_ambiguity_count < 0:
        errors.append("ar_eligible_ambiguity_count_negative")
    if epoch.ar_integer_ambiguity_coordinate_count < 0:
        errors.append("ar_integer_ambiguity_coordinate_count_negative")
    if (
        epoch.ar_integer_ambiguity_coordinate_count
        > epoch.ar_eligible_ambiguity_count
    ):
        errors.append("ar_integer_ambiguity_coordinate_count_exceeds_eligible_count")
    if epoch.ar_receiver_datum_group_count < 0:
        errors.append("ar_receiver_datum_group_count_negative")
    if epoch.ar_dropped_singleton_group_count < 0:
        errors.append("ar_dropped_singleton_group_count_negative")
    if epoch.schema == "GINAN_STEC_COVARIANCE_V3" and (
        epoch.ar_integer_ambiguity_coordinate_count
        + epoch.ar_receiver_datum_group_count
        + epoch.ar_dropped_singleton_group_count
        != epoch.ar_eligible_ambiguity_count
    ):
        errors.append("ar_integer_coordinate_accounting_inconsistent")
    if (
        epoch.ar_receiver_single_difference_applied
        and epoch.ar_receiver_datum_group_count <= 0
    ):
        errors.append("ar_receiver_single_difference_without_datum_group")
    if (
        not epoch.ar_receiver_single_difference_applied
        and epoch.ar_receiver_datum_group_count != 0
    ):
        errors.append("ar_receiver_datum_group_without_single_difference")
    if epoch.ar_resolved_combination_count < 0:
        errors.append("ar_resolved_combination_count_negative")
    if (
        epoch.ar_resolved_combination_count
        > epoch.ar_integer_ambiguity_coordinate_count
    ):
        errors.append("ar_resolved_combination_count_exceeds_integer_coordinate_count")
    if epoch.ar_selected_decorrelated_ambiguity_count < 0:
        errors.append("ar_selected_decorrelated_ambiguity_count_negative")
    if (
        epoch.ar_selected_decorrelated_ambiguity_count
        > epoch.ar_integer_ambiguity_coordinate_count
    ):
        errors.append("ar_selected_decorrelated_count_exceeds_integer_coordinate_count")
    if epoch.ar_integer_candidate_count < 0:
        errors.append("ar_integer_candidate_count_negative")
    if not epoch.ar_routine_invoked and (
        epoch.ar_eligible_ambiguity_count != 0
        or epoch.ar_integer_ambiguity_coordinate_count != 0
        or epoch.ar_receiver_single_difference_applied
        or epoch.ar_receiver_datum_group_count != 0
        or epoch.ar_dropped_singleton_group_count != 0
        or epoch.ar_resolved_combination_count != 0
        or epoch.ar_pseudoobservations_submitted
    ):
        errors.append("ar_evidence_present_without_routine_invocation")
    if not epoch.ar_routine_invoked and (
        epoch.ar_diagnostic_status != "NOT_RUN"
        or epoch.ar_selected_decorrelated_ambiguity_count != 0
        or epoch.ar_integer_candidate_count != 0
        or epoch.ar_bootstrapped_success_rate != -1
        or epoch.ar_best_squared_norm != -1
        or epoch.ar_second_squared_norm != -1
        or epoch.ar_solution_ratio != -1
    ):
        errors.append("ar_diagnostics_present_without_routine_invocation")
    if epoch.ar_pseudoobservations_submitted and epoch.ar_resolved_combination_count <= 0:
        errors.append("ar_pseudoobservations_submitted_without_resolved_combinations")
    if not math.isfinite(epoch.ar_configured_success_rate_threshold):
        errors.append("ar_configured_success_rate_threshold_nonfinite")
    elif not 0 <= epoch.ar_configured_success_rate_threshold <= 1:
        errors.append("ar_configured_success_rate_threshold_out_of_range")
    if not math.isfinite(epoch.ar_configured_solution_ratio_threshold):
        errors.append("ar_configured_solution_ratio_threshold_nonfinite")
    elif epoch.ar_configured_solution_ratio_threshold <= 0:
        errors.append("ar_configured_solution_ratio_threshold_not_positive")
    for field_name, value in (
        ("ar_bootstrapped_success_rate", epoch.ar_bootstrapped_success_rate),
        ("ar_best_squared_norm", epoch.ar_best_squared_norm),
        ("ar_second_squared_norm", epoch.ar_second_squared_norm),
        ("ar_solution_ratio", epoch.ar_solution_ratio),
    ):
        if not math.isfinite(value):
            errors.append(f"{field_name}_nonfinite")
    if epoch.ar_bootstrapped_success_rate != -1 and not (
        0 <= epoch.ar_bootstrapped_success_rate <= 1
    ):
        errors.append("ar_bootstrapped_success_rate_out_of_range")
    for field_name, value in (
        ("ar_best_squared_norm", epoch.ar_best_squared_norm),
        ("ar_second_squared_norm", epoch.ar_second_squared_norm),
        ("ar_solution_ratio", epoch.ar_solution_ratio),
    ):
        if value < 0 and value != -1:
            errors.append(f"{field_name}_invalid_negative_sentinel")
    if epoch.ar_routine_invoked and epoch.ar_diagnostic_status == "NOT_RUN":
        errors.append("ar_routine_invoked_without_diagnostic_status")
    if epoch.ar_diagnostic_status == "SUCCESS_RATE_BELOW_THRESHOLD" and not (
        0 <= epoch.ar_bootstrapped_success_rate
        < epoch.ar_configured_success_rate_threshold
    ):
        errors.append("ar_success_rate_failure_status_inconsistent")
    if epoch.ar_diagnostic_status == "RATIO_BELOW_THRESHOLD" and not (
        0 <= epoch.ar_solution_ratio < epoch.ar_configured_solution_ratio_threshold
    ):
        errors.append("ar_ratio_failure_status_inconsistent")
    if n <= 0:
        errors.append("state_count_not_positive")
    if epoch.upper_triangle_count != expected_entries:
        errors.append(
            f"upper_triangle_count={epoch.upper_triangle_count},expected={expected_entries}"
        )
    if sorted(epoch.states) != list(range(n)):
        errors.append("state_catalog_incomplete_or_noncontiguous")
    if len(epoch.covariance) != expected_entries:
        errors.append(f"covariance_entries={len(epoch.covariance)},expected={expected_entries}")

    matrix = np.full((n, n), np.nan, dtype=float)
    for (row, column), value in epoch.covariance.items():
        if row < 0 or column < row or column >= n:
            errors.append(f"invalid_upper_triangle_index=({row},{column})")
            continue
        matrix[row, column] = value
        matrix[column, row] = value

    estimates = np.array(
        [epoch.states[index].estimate_tecu for index in range(n) if index in epoch.states],
        dtype=float,
    )
    finite = estimates.size == n and np.isfinite(estimates).all() and np.isfinite(matrix).all()
    if not finite:
        errors.append("nonfinite_estimate_or_covariance")

    diagonal_error = math.inf
    minimum_eigenvalue = math.nan
    maximum_eigenvalue = math.nan
    condition_number = math.inf
    psd_tolerance = math.nan

    if finite:
        state_variances = np.array(
            [epoch.states[index].variance_tecu2 for index in range(n)], dtype=float
        )
        diagonal_error = float(np.max(np.abs(np.diag(matrix) - state_variances)))
        scale = max(1.0, float(np.max(np.abs(matrix))))
        diagonal_tolerance = psd_absolute_tolerance + psd_relative_tolerance * scale
        if diagonal_error > diagonal_tolerance:
            errors.append(
                f"state_variance_diagonal_mismatch={diagonal_error:.17g},"
                f"tolerance={diagonal_tolerance:.17g}"
            )

        eigenvalues = np.linalg.eigvalsh(matrix)
        minimum_eigenvalue = float(eigenvalues[0])
        maximum_eigenvalue = float(eigenvalues[-1])
        psd_tolerance = psd_absolute_tolerance + psd_relative_tolerance * scale
        if minimum_eigenvalue < -psd_tolerance:
            errors.append(
                f"not_positive_semidefinite_min_eigenvalue={minimum_eigenvalue:.17g},"
                f"tolerance={psd_tolerance:.17g}"
            )

        positive = eigenvalues[eigenvalues > psd_tolerance]
        if positive.size:
            condition_number = float(positive[-1] / positive[0])

    return {
        "gps_week": epoch.gps_week,
        "gps_tow": float(epoch.gps_tow_text),
        "posterior_stage": epoch.posterior_stage,
        "ar_routine_invoked": epoch.ar_routine_invoked,
        "ar_eligible_ambiguity_count": epoch.ar_eligible_ambiguity_count,
        "ar_integer_ambiguity_coordinate_count": epoch.ar_integer_ambiguity_coordinate_count,
        "ar_receiver_single_difference_applied": epoch.ar_receiver_single_difference_applied,
        "ar_receiver_datum_group_count": epoch.ar_receiver_datum_group_count,
        "ar_dropped_singleton_group_count": epoch.ar_dropped_singleton_group_count,
        "ar_resolved_combination_count": epoch.ar_resolved_combination_count,
        "ar_pseudoobservations_submitted": epoch.ar_pseudoobservations_submitted,
        "ar_mode": epoch.ar_mode,
        "ar_configured_success_rate_threshold": epoch.ar_configured_success_rate_threshold,
        "ar_configured_solution_ratio_threshold": epoch.ar_configured_solution_ratio_threshold,
        "ar_diagnostic_status": epoch.ar_diagnostic_status,
        "ar_selected_decorrelated_ambiguity_count": epoch.ar_selected_decorrelated_ambiguity_count,
        "ar_integer_candidate_count": epoch.ar_integer_candidate_count,
        "ar_bootstrapped_success_rate": epoch.ar_bootstrapped_success_rate,
        "ar_best_squared_norm": epoch.ar_best_squared_norm,
        "ar_second_squared_norm": epoch.ar_second_squared_norm,
        "ar_solution_ratio": epoch.ar_solution_ratio,
        "writer_status": epoch.status,
        "state_count": n,
        "upper_triangle_count": epoch.upper_triangle_count,
        "writer_max_abs_asymmetry_tecu2": epoch.writer_max_abs_asymmetry_tecu2,
        "diagonal_error_tecu2": diagonal_error,
        "minimum_eigenvalue_tecu2": minimum_eigenvalue,
        "maximum_eigenvalue_tecu2": maximum_eigenvalue,
        "psd_tolerance_tecu2": psd_tolerance,
        "positive_spectrum_condition_number": condition_number,
        "valid": not errors,
        "errors": errors,
    }


def validate_file(
    path: Path,
    psd_absolute_tolerance: float = 1e-12,
    psd_relative_tolerance: float = 1e-10,
) -> dict[str, object]:
    epoch_results: list[dict[str, object]] = []
    schema: str | None = None
    with path.open("r", encoding="utf-8", newline="") as stream:
        for epoch in iter_epoch_blocks(stream):
            schema = epoch.schema
            epoch_results.append(
                validate_epoch(epoch, psd_absolute_tolerance, psd_relative_tolerance)
            )

    valid_epochs = sum(bool(epoch["valid"]) for epoch in epoch_results)
    return {
        "schema": schema,
        "input": str(path.resolve()),
        "epoch_count": len(epoch_results),
        "valid_epoch_count": valid_epochs,
        "invalid_epoch_count": len(epoch_results) - valid_epochs,
        "all_valid": bool(epoch_results) and valid_epochs == len(epoch_results),
        "epochs": epoch_results,
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="GINAN_STEC_COVARIANCE_V2/V3 file")
    parser.add_argument("--json-output", type=Path, help="write the full validation report")
    parser.add_argument("--psd-absolute-tolerance", type=float, default=1e-12)
    parser.add_argument("--psd-relative-tolerance", type=float, default=1e-10)
    return parser


def main(argv: Iterable[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        report = validate_file(
            args.input,
            psd_absolute_tolerance=args.psd_absolute_tolerance,
            psd_relative_tolerance=args.psd_relative_tolerance,
        )
    except (OSError, ValueError, csv.Error) as error:
        print(f"STEC covariance validation error: {error}", file=sys.stderr)
        return 2

    summary = {
        key: report[key]
        for key in ("schema", "input", "epoch_count", "valid_epoch_count", "invalid_epoch_count", "all_valid")
    }
    print(json.dumps(summary, ensure_ascii=False, indent=2))

    if args.json_output:
        args.json_output.parent.mkdir(parents=True, exist_ok=True)
        args.json_output.write_text(
            json.dumps(report, ensure_ascii=False, indent=2, allow_nan=True) + "\n",
            encoding="utf-8",
        )

    return 0 if report["all_valid"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
