#!/usr/bin/env python3
"""Validate a raw STEC V3 sidecar against its satellite-difference V1 sidecar.

The validator reconstructs the sparse satellite-difference transform ``D`` at
every epoch and independently checks ``x_sd = D x`` and
``C_sd = D C D.T``.  It also audits deterministic datum construction, complete
upper-triangle serialization, diagonal consistency, symmetry, and positive
semidefiniteness.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import sys
import tempfile
from dataclasses import dataclass, field
from decimal import Decimal, InvalidOperation
from pathlib import Path
from typing import Iterator, Sequence, TextIO

import numpy as np


RAW_SCHEMA = "GINAN_STEC_COVARIANCE_V3"
SD_SCHEMA = "GINAN_STEC_SATELLITE_DIFFERENCE_COVARIANCE_V1"
REPORT_SCHEMA = "GINAN_STEC_SATELLITE_DIFFERENCE_COVARIANCE_AUDIT_V1"
MAX_AUDIT_STATES = 10000


class ValidationInputError(ValueError):
    """Raised when a sidecar is structurally unparseable."""


@dataclass(frozen=True)
class RawState:
    local_index: int
    site: str
    satellite: str
    state_number: int
    filter_index: int
    estimate_tecu: float
    variance_tecu2: float


@dataclass
class RawEpoch:
    gps_week: int
    gps_tow_text: str
    status: str
    state_count: int
    upper_triangle_count: int
    max_abs_asymmetry_tecu2: float
    posterior_stage: str
    states: dict[int, RawState] = field(default_factory=dict)
    covariance: dict[tuple[int, int], float] = field(default_factory=dict)


@dataclass(frozen=True)
class DatumRow:
    datum_index: int
    site: str
    constellation: str
    state_number: int
    reference_satellite: str
    member_count: int
    status: str


@dataclass(frozen=True)
class TransformRow:
    difference_local_index: int
    source_local_index: int
    coefficient: float


@dataclass(frozen=True)
class DifferenceState:
    local_index: int
    datum_index: int
    site: str
    constellation: str
    target_satellite: str
    reference_satellite: str
    state_number: int
    target_source_local_index: int
    reference_source_local_index: int
    estimate_tecu: float
    variance_tecu2: float


@dataclass
class DifferenceEpoch:
    gps_week: int
    gps_tow_text: str
    status: str
    source_state_count: int
    difference_state_count: int
    datum_count: int
    singleton_datum_count: int
    upper_triangle_count: int
    max_abs_asymmetry_tecu2: float
    posterior_stage: str
    datums: dict[int, DatumRow] = field(default_factory=dict)
    transforms: list[TransformRow] = field(default_factory=list)
    states: dict[int, DifferenceState] = field(default_factory=dict)
    covariance: dict[tuple[int, int], float] = field(default_factory=dict)


def _require_length(row: list[str], expected: int, line_number: int) -> None:
    if len(row) != expected:
        record = row[0] if row else "EMPTY"
        raise ValidationInputError(
            f"line {line_number}: {record} expects {expected} fields, got {len(row)}"
        )


def _check_record_epoch(
    row: list[str],
    gps_week: int,
    gps_tow_text: str,
    line_number: int,
) -> None:
    if int(row[1]) != gps_week or row[2] != gps_tow_text:
        raise ValidationInputError(
            f"line {line_number}: record epoch differs from current META"
        )


def iter_raw_epochs(stream: TextIO) -> Iterator[RawEpoch]:
    schema_seen = False
    current: RawEpoch | None = None
    for line_number, raw_line in enumerate(stream, start=1):
        stripped = raw_line.strip()
        if stripped == f"# {RAW_SCHEMA}":
            schema_seen = True
            continue
        if not stripped or raw_line.startswith("#"):
            continue
        row = next(csv.reader([raw_line]))
        record = row[0]
        if record == "META":
            _require_length(row, 26, line_number)
            if not schema_seen:
                raise ValidationInputError(
                    f"line {line_number}: META precedes {RAW_SCHEMA} marker"
                )
            if current is not None:
                yield current
            current = RawEpoch(
                gps_week=int(row[1]),
                gps_tow_text=row[2],
                status=row[3],
                state_count=int(row[4]),
                upper_triangle_count=int(row[5]),
                max_abs_asymmetry_tecu2=float(row[6]),
                posterior_stage=row[7],
            )
            continue
        if current is None:
            raise ValidationInputError(f"line {line_number}: {record} precedes META")
        _check_record_epoch(row, current.gps_week, current.gps_tow_text, line_number)
        if record == "STATE":
            _require_length(row, 10, line_number)
            state = RawState(
                local_index=int(row[3]),
                site=row[4],
                satellite=row[5],
                state_number=int(row[6]),
                filter_index=int(row[7]),
                estimate_tecu=float(row[8]),
                variance_tecu2=float(row[9]),
            )
            if state.local_index in current.states:
                raise ValidationInputError(
                    f"line {line_number}: duplicate raw STATE index {state.local_index}"
                )
            current.states[state.local_index] = state
        elif record == "COV":
            _require_length(row, 6, line_number)
            pair = (int(row[3]), int(row[4]))
            if pair in current.covariance:
                raise ValidationInputError(
                    f"line {line_number}: duplicate raw COV pair {pair}"
                )
            current.covariance[pair] = float(row[5])
        else:
            raise ValidationInputError(
                f"line {line_number}: unknown raw record type {record!r}"
            )
    if not schema_seen:
        raise ValidationInputError(f"missing {RAW_SCHEMA} schema marker")
    if current is not None:
        yield current


def iter_difference_epochs(stream: TextIO) -> Iterator[DifferenceEpoch]:
    schema_seen = False
    current: DifferenceEpoch | None = None
    transform_pairs: set[tuple[int, int]] = set()
    for line_number, raw_line in enumerate(stream, start=1):
        stripped = raw_line.strip()
        if stripped == f"# {SD_SCHEMA}":
            schema_seen = True
            continue
        if not stripped or raw_line.startswith("#"):
            continue
        row = next(csv.reader([raw_line]))
        record = row[0]
        if record == "META":
            _require_length(row, 11, line_number)
            if not schema_seen:
                raise ValidationInputError(
                    f"line {line_number}: META precedes {SD_SCHEMA} marker"
                )
            if current is not None:
                yield current
            current = DifferenceEpoch(
                gps_week=int(row[1]),
                gps_tow_text=row[2],
                status=row[3],
                source_state_count=int(row[4]),
                difference_state_count=int(row[5]),
                datum_count=int(row[6]),
                singleton_datum_count=int(row[7]),
                upper_triangle_count=int(row[8]),
                max_abs_asymmetry_tecu2=float(row[9]),
                posterior_stage=row[10],
            )
            transform_pairs = set()
            continue
        if current is None:
            raise ValidationInputError(f"line {line_number}: {record} precedes META")
        _check_record_epoch(row, current.gps_week, current.gps_tow_text, line_number)
        if record == "DATUM":
            _require_length(row, 10, line_number)
            datum = DatumRow(
                datum_index=int(row[3]),
                site=row[4],
                constellation=row[5],
                state_number=int(row[6]),
                reference_satellite=row[7],
                member_count=int(row[8]),
                status=row[9],
            )
            if datum.datum_index in current.datums:
                raise ValidationInputError(
                    f"line {line_number}: duplicate DATUM index {datum.datum_index}"
                )
            current.datums[datum.datum_index] = datum
        elif record == "TRANSFORM":
            _require_length(row, 6, line_number)
            transform = TransformRow(
                difference_local_index=int(row[3]),
                source_local_index=int(row[4]),
                coefficient=float(row[5]),
            )
            pair = (transform.difference_local_index, transform.source_local_index)
            if pair in transform_pairs:
                raise ValidationInputError(
                    f"line {line_number}: duplicate TRANSFORM pair {pair}"
                )
            transform_pairs.add(pair)
            current.transforms.append(transform)
        elif record == "SD_STATE":
            _require_length(row, 14, line_number)
            state = DifferenceState(
                local_index=int(row[3]),
                datum_index=int(row[4]),
                site=row[5],
                constellation=row[6],
                target_satellite=row[7],
                reference_satellite=row[8],
                state_number=int(row[9]),
                target_source_local_index=int(row[10]),
                reference_source_local_index=int(row[11]),
                estimate_tecu=float(row[12]),
                variance_tecu2=float(row[13]),
            )
            if state.local_index in current.states:
                raise ValidationInputError(
                    f"line {line_number}: duplicate SD_STATE index {state.local_index}"
                )
            current.states[state.local_index] = state
        elif record == "SD_COV":
            _require_length(row, 6, line_number)
            pair = (int(row[3]), int(row[4]))
            if pair in current.covariance:
                raise ValidationInputError(
                    f"line {line_number}: duplicate SD_COV pair {pair}"
                )
            current.covariance[pair] = float(row[5])
        else:
            raise ValidationInputError(
                f"line {line_number}: unknown satellite-difference record type {record!r}"
            )
    if not schema_seen:
        raise ValidationInputError(f"missing {SD_SCHEMA} schema marker")
    if current is not None:
        yield current


def _epoch_key(gps_week: int, gps_tow_text: str) -> tuple[int, Decimal]:
    if gps_week < 0:
        raise ValidationInputError(f"negative gps_week {gps_week}")
    try:
        tow = Decimal(gps_tow_text)
    except InvalidOperation as error:
        raise ValidationInputError(f"invalid gps_tow {gps_tow_text!r}") from error
    if not tow.is_finite():
        raise ValidationInputError(f"nonfinite gps_tow {gps_tow_text!r}")
    if tow < 0 or tow >= Decimal(604800):
        raise ValidationInputError(f"gps_tow outside [0,604800): {gps_tow_text!r}")
    return gps_week, tow


def _key_text(key: tuple[int, Decimal]) -> str:
    return f"{key[0]}:{key[1]}"


def _matrix_from_upper_triangle(
    count: int,
    declared_count: int,
    entries: dict[tuple[int, int], float],
    prefix: str,
    errors: list[str],
) -> np.ndarray | None:
    if count <= 0:
        errors.append(f"{prefix}_state_count_not_positive={count}")
        return None
    if count > MAX_AUDIT_STATES:
        errors.append(f"{prefix}_state_count_exceeds_audit_limit={count}")
        return None
    expected_count = count * (count + 1) // 2
    if declared_count != expected_count:
        errors.append(
            f"{prefix}_upper_triangle_count={declared_count},expected={expected_count}"
        )
    expected_pairs = {(row, column) for row in range(count) for column in range(row, count)}
    actual_pairs = set(entries)
    missing = sorted(expected_pairs - actual_pairs)
    extra = sorted(actual_pairs - expected_pairs)
    if missing:
        errors.append(f"{prefix}_upper_triangle_missing_pairs={len(missing)}:{missing[:5]}")
    if extra:
        errors.append(f"{prefix}_invalid_upper_triangle_pairs={len(extra)}:{extra[:5]}")
    if len(entries) != expected_count:
        errors.append(
            f"{prefix}_covariance_entry_count={len(entries)},expected={expected_count}"
        )
    matrix = np.full((count, count), np.nan, dtype=float)
    for (row, column), value in entries.items():
        if (row, column) not in expected_pairs:
            continue
        matrix[row, column] = value
        matrix[column, row] = value
    return matrix


def _array_comparison(
    actual: np.ndarray,
    expected: np.ndarray,
    absolute_tolerance: float,
    relative_tolerance: float,
) -> tuple[bool, float | None, float | None]:
    if actual.shape != expected.shape or not np.isfinite(actual).all() or not np.isfinite(expected).all():
        return False, None, None
    difference = np.abs(actual - expected)
    if difference.size == 0:
        return True, 0.0, 0.0
    limits = absolute_tolerance + relative_tolerance * np.abs(expected)
    normalized = difference / np.maximum(limits, np.finfo(float).tiny)
    return (
        bool(np.all(difference <= limits)),
        float(np.max(difference)),
        float(np.max(normalized)),
    )


def _psd_metrics(
    matrix: np.ndarray,
    absolute_tolerance: float,
    relative_tolerance: float,
) -> tuple[float | None, float | None, bool]:
    if not np.isfinite(matrix).all():
        return None, None, False
    scale = max(1.0, float(np.max(np.abs(matrix))))
    tolerance = absolute_tolerance + relative_tolerance * scale
    minimum_eigenvalue = float(np.linalg.eigvalsh((matrix + matrix.T) / 2)[0])
    return minimum_eigenvalue, tolerance, minimum_eigenvalue >= -tolerance


def _satellite_number(satellite: str) -> int:
    if len(satellite) < 2:
        return sys.maxsize
    suffix = satellite[1:]
    if not suffix.isdigit():
        return sys.maxsize
    return int(suffix)


def _satellite_order(state: RawState) -> tuple[int, str]:
    return _satellite_number(state.satellite), state.satellite


def _expected_groups(raw: RawEpoch) -> list[tuple[tuple[str, str, int], list[int]]]:
    groups: dict[tuple[str, str, int], list[int]] = {}
    for source_index in range(raw.state_count):
        state = raw.states[source_index]
        constellation = state.satellite[:1]
        groups.setdefault((state.site, constellation, state.state_number), []).append(source_index)
    result: list[tuple[tuple[str, str, int], list[int]]] = []
    for key in sorted(groups):
        members = sorted(groups[key], key=lambda index: _satellite_order(raw.states[index]))
        result.append((key, members))
    return result


def validate_epoch_pair(
    raw: RawEpoch,
    difference: DifferenceEpoch,
    absolute_tolerance: float,
    relative_tolerance: float,
    psd_absolute_tolerance: float,
    psd_relative_tolerance: float,
) -> dict[str, object]:
    errors: list[str] = []
    metrics: dict[str, float | None] = {
        "raw_diagonal_max_abs_error_tecu2": None,
        "sd_state_transform_max_abs_error_tecu": None,
        "sd_state_transform_max_normalized_error": None,
        "sd_covariance_transform_max_abs_error_tecu2": None,
        "sd_covariance_transform_max_normalized_error": None,
        "sd_diagonal_max_abs_error_tecu2": None,
        "sd_predicted_diagonal_max_abs_error_tecu2": None,
        "raw_symmetry_max_abs_error_tecu2": None,
        "sd_symmetry_max_abs_error_tecu2": None,
        "raw_minimum_eigenvalue_tecu2": None,
        "raw_psd_tolerance_tecu2": None,
        "sd_minimum_eigenvalue_tecu2": None,
        "sd_psd_tolerance_tecu2": None,
    }

    if raw.status != "OK":
        errors.append(f"raw_writer_status={raw.status}")
    if difference.status != "OK":
        errors.append(f"sd_writer_status={difference.status}")
    if raw.posterior_stage != difference.posterior_stage:
        errors.append(
            f"posterior_stage_mismatch=raw:{raw.posterior_stage},sd:{difference.posterior_stage}"
        )
    if difference.source_state_count != raw.state_count:
        errors.append(
            f"source_state_count_mismatch=raw:{raw.state_count},sd:{difference.source_state_count}"
        )

    n = raw.state_count
    m = difference.difference_state_count
    if sorted(raw.states) != list(range(n)):
        errors.append("raw_state_catalog_incomplete_or_noncontiguous")
    if sorted(difference.states) != list(range(m)):
        errors.append("sd_state_catalog_incomplete_or_noncontiguous")
    if difference.datum_count != len(difference.datums):
        errors.append(
            f"datum_count={difference.datum_count},records={len(difference.datums)}"
        )
    if sorted(difference.datums) != list(range(difference.datum_count)):
        errors.append("datum_catalog_incomplete_or_noncontiguous")

    raw_matrix = _matrix_from_upper_triangle(
        n,
        raw.upper_triangle_count,
        raw.covariance,
        "raw",
        errors,
    )
    sd_matrix = _matrix_from_upper_triangle(
        m,
        difference.upper_triangle_count,
        difference.covariance,
        "sd",
        errors,
    )

    raw_catalog_complete = sorted(raw.states) == list(range(n))
    sd_catalog_complete = sorted(difference.states) == list(range(m))
    raw_estimate: np.ndarray | None = None
    if raw_catalog_complete:
        raw_estimate = np.array(
            [raw.states[index].estimate_tecu for index in range(n)], dtype=float
        )
        raw_variances = np.array(
            [raw.states[index].variance_tecu2 for index in range(n)], dtype=float
        )
        if not np.isfinite(raw_estimate).all() or not np.isfinite(raw_variances).all():
            errors.append("raw_state_nonfinite")
        if raw_matrix is not None and np.isfinite(raw_matrix).all():
            valid, maximum, _ = _array_comparison(
                raw_variances,
                np.diag(raw_matrix),
                absolute_tolerance,
                relative_tolerance,
            )
            metrics["raw_diagonal_max_abs_error_tecu2"] = maximum
            if not valid:
                errors.append(f"raw_state_variance_diagonal_mismatch={maximum}")

    expected_groups: list[tuple[tuple[str, str, int], list[int]]] = []
    expected_difference_rows: list[dict[str, object]] = []
    expected_transform: np.ndarray | None = None
    if raw_catalog_complete:
        expected_groups = _expected_groups(raw)
        expected_singletons = sum(len(members) == 1 for _, members in expected_groups)
        expected_difference_count = sum(len(members) - 1 for _, members in expected_groups)
        if difference.datum_count != len(expected_groups):
            errors.append(
                f"datum_group_count={difference.datum_count},expected={len(expected_groups)}"
            )
        if difference.singleton_datum_count != expected_singletons:
            errors.append(
                "singleton_datum_count="
                f"{difference.singleton_datum_count},expected={expected_singletons}"
            )
        if m != expected_difference_count:
            errors.append(
                f"difference_state_count={m},expected={expected_difference_count}"
            )
        if m != n - len(expected_groups):
            errors.append(
                f"difference_group_accounting={m},expected_source_minus_datums={n - len(expected_groups)}"
            )

        for datum_index, (group_key, members) in enumerate(expected_groups):
            site, constellation, state_number = group_key
            reference_source_index = members[0]
            reference_satellite = raw.states[reference_source_index].satellite
            satellites = [raw.states[index].satellite for index in members]
            if len(set(satellites)) != len(satellites):
                errors.append(f"datum_{datum_index}_duplicate_satellites={satellites}")
            expected_status = (
                "DIFFERENCES_CREATED" if len(members) >= 2 else "SINGLETON_DROPPED"
            )
            actual_datum = difference.datums.get(datum_index)
            expected_datum = (
                site,
                constellation,
                state_number,
                reference_satellite,
                len(members),
                expected_status,
            )
            if actual_datum is None:
                errors.append(f"missing_datum_index={datum_index}")
            else:
                actual = (
                    actual_datum.site,
                    actual_datum.constellation,
                    actual_datum.state_number,
                    actual_datum.reference_satellite,
                    actual_datum.member_count,
                    actual_datum.status,
                )
                if actual != expected_datum:
                    errors.append(
                        f"datum_{datum_index}_deterministic_definition_mismatch="
                        f"actual:{actual},expected:{expected_datum}"
                    )
            for target_source_index in members[1:]:
                expected_difference_rows.append(
                    {
                        "datum_index": datum_index,
                        "site": site,
                        "constellation": constellation,
                        "target_satellite": raw.states[target_source_index].satellite,
                        "reference_satellite": reference_satellite,
                        "state_number": state_number,
                        "target_source_local_index": target_source_index,
                        "reference_source_local_index": reference_source_index,
                    }
                )

        if (
            len(expected_difference_rows) == m
            and 0 < m <= MAX_AUDIT_STATES
            and 0 < n <= MAX_AUDIT_STATES
        ):
            expected_transform = np.zeros((m, n), dtype=float)
            for local_index, expected in enumerate(expected_difference_rows):
                expected_transform[
                    local_index, int(expected["target_source_local_index"])
                ] = 1
                expected_transform[
                    local_index, int(expected["reference_source_local_index"])
                ] = -1

    datum_member_sum = sum(datum.member_count for datum in difference.datums.values())
    if datum_member_sum != difference.source_state_count:
        errors.append(
            f"datum_member_count_sum={datum_member_sum},source_state_count={difference.source_state_count}"
        )

    transform: np.ndarray | None = None
    transform_structure_valid = 0 < m <= MAX_AUDIT_STATES and 0 < n <= MAX_AUDIT_STATES
    if len(difference.transforms) != 2 * m:
        errors.append(
            f"transform_term_count={len(difference.transforms)},expected={2 * m}"
        )
        transform_structure_valid = False
    terms_by_row: dict[int, list[TransformRow]] = {}
    for term in difference.transforms:
        terms_by_row.setdefault(term.difference_local_index, []).append(term)
        if not 0 <= term.difference_local_index < m:
            errors.append(f"transform_invalid_difference_index={term.difference_local_index}")
            transform_structure_valid = False
        if not 0 <= term.source_local_index < n:
            errors.append(f"transform_invalid_source_index={term.source_local_index}")
            transform_structure_valid = False
        if not math.isfinite(term.coefficient):
            errors.append("transform_nonfinite_coefficient")
            transform_structure_valid = False

    for local_index in range(max(0, m)):
        terms = terms_by_row.get(local_index, [])
        coefficients = sorted(term.coefficient for term in terms)
        source_indices = [term.source_local_index for term in terms]
        if len(terms) != 2 or coefficients != [-1.0, 1.0] or len(set(source_indices)) != 2:
            errors.append(
                f"transform_row_{local_index}_must_have_two_distinct_plus_minus_one_terms="
                f"terms:{[(term.source_local_index, term.coefficient) for term in terms]}"
            )
            transform_structure_valid = False

    if transform_structure_valid:
        transform = np.zeros((m, n), dtype=float)
        for term in difference.transforms:
            transform[term.difference_local_index, term.source_local_index] = term.coefficient
        if expected_transform is not None and not np.array_equal(transform, expected_transform):
            errors.append("transform_differs_from_deterministic_datum_transform")

    if sd_catalog_complete and expected_difference_rows:
        for local_index, expected in enumerate(expected_difference_rows):
            if local_index not in difference.states:
                continue
            state = difference.states[local_index]
            actual = {
                "datum_index": state.datum_index,
                "site": state.site,
                "constellation": state.constellation,
                "target_satellite": state.target_satellite,
                "reference_satellite": state.reference_satellite,
                "state_number": state.state_number,
                "target_source_local_index": state.target_source_local_index,
                "reference_source_local_index": state.reference_source_local_index,
            }
            if actual != expected:
                errors.append(
                    f"sd_state_{local_index}_deterministic_definition_mismatch="
                    f"actual:{actual},expected:{expected}"
                )

    if transform is not None and raw_estimate is not None and sd_catalog_complete:
        sd_estimate = np.array(
            [difference.states[index].estimate_tecu for index in range(m)], dtype=float
        )
        predicted_estimate = transform @ raw_estimate
        valid, maximum, normalized = _array_comparison(
            sd_estimate,
            predicted_estimate,
            absolute_tolerance,
            relative_tolerance,
        )
        metrics["sd_state_transform_max_abs_error_tecu"] = maximum
        metrics["sd_state_transform_max_normalized_error"] = normalized
        if not valid:
            errors.append(f"sd_state_not_equal_D_x_max_abs_error={maximum}")

    predicted_covariance: np.ndarray | None = None
    if (
        transform is not None
        and raw_matrix is not None
        and sd_matrix is not None
        and np.isfinite(raw_matrix).all()
        and np.isfinite(sd_matrix).all()
    ):
        predicted_covariance = transform @ raw_matrix @ transform.T
        valid, maximum, normalized = _array_comparison(
            sd_matrix,
            predicted_covariance,
            absolute_tolerance,
            relative_tolerance,
        )
        metrics["sd_covariance_transform_max_abs_error_tecu2"] = maximum
        metrics["sd_covariance_transform_max_normalized_error"] = normalized
        if not valid:
            errors.append(f"sd_covariance_not_equal_D_C_Dt_max_abs_error={maximum}")

    if sd_catalog_complete and sd_matrix is not None and np.isfinite(sd_matrix).all():
        sd_variances = np.array(
            [difference.states[index].variance_tecu2 for index in range(m)], dtype=float
        )
        valid, maximum, _ = _array_comparison(
            sd_variances,
            np.diag(sd_matrix),
            absolute_tolerance,
            relative_tolerance,
        )
        metrics["sd_diagonal_max_abs_error_tecu2"] = maximum
        if not valid:
            errors.append(f"sd_state_variance_diagonal_mismatch={maximum}")
        if predicted_covariance is not None:
            valid, maximum, _ = _array_comparison(
                sd_variances,
                np.diag(predicted_covariance),
                absolute_tolerance,
                relative_tolerance,
            )
            metrics["sd_predicted_diagonal_max_abs_error_tecu2"] = maximum
            if not valid:
                errors.append(f"sd_state_variance_predicted_diagonal_mismatch={maximum}")

    for prefix, matrix, writer_asymmetry in (
        ("raw", raw_matrix, raw.max_abs_asymmetry_tecu2),
        ("sd", sd_matrix, difference.max_abs_asymmetry_tecu2),
    ):
        if not math.isfinite(writer_asymmetry) or writer_asymmetry < 0:
            errors.append(f"{prefix}_writer_asymmetry_invalid={writer_asymmetry}")
        if matrix is None or not np.isfinite(matrix).all():
            errors.append(f"{prefix}_covariance_nonfinite_or_incomplete")
            continue
        symmetry_valid, symmetry_error, _ = _array_comparison(
            matrix,
            matrix.T,
            absolute_tolerance,
            relative_tolerance,
        )
        metrics[f"{prefix}_symmetry_max_abs_error_tecu2"] = symmetry_error
        if not symmetry_valid:
            errors.append(f"{prefix}_covariance_not_symmetric={symmetry_error}")
        scale = max(1.0, float(np.max(np.abs(matrix))))
        writer_limit = absolute_tolerance + relative_tolerance * scale
        if math.isfinite(writer_asymmetry) and writer_asymmetry > writer_limit:
            errors.append(
                f"{prefix}_writer_asymmetry={writer_asymmetry},tolerance={writer_limit}"
            )
        minimum_eigenvalue, psd_tolerance, psd_valid = _psd_metrics(
            matrix,
            psd_absolute_tolerance,
            psd_relative_tolerance,
        )
        metrics[f"{prefix}_minimum_eigenvalue_tecu2"] = minimum_eigenvalue
        metrics[f"{prefix}_psd_tolerance_tecu2"] = psd_tolerance
        if not psd_valid:
            errors.append(
                f"{prefix}_not_positive_semidefinite_min_eigenvalue="
                f"{minimum_eigenvalue},tolerance={psd_tolerance}"
            )

    return {
        "gps_week": raw.gps_week,
        "gps_tow": float(Decimal(raw.gps_tow_text)),
        "posterior_stage": raw.posterior_stage,
        "source_state_count": n,
        "difference_state_count": m,
        "datum_count": difference.datum_count,
        "singleton_datum_count": difference.singleton_datum_count,
        "transform_term_count": len(difference.transforms),
        **metrics,
        "valid": not errors,
        "errors": errors,
    }


def validate_files(
    raw_path: Path,
    difference_path: Path,
    absolute_tolerance: float = 1e-10,
    relative_tolerance: float = 1e-10,
    psd_absolute_tolerance: float = 1e-12,
    psd_relative_tolerance: float = 1e-10,
) -> dict[str, object]:
    for name, value in (
        ("absolute_tolerance", absolute_tolerance),
        ("relative_tolerance", relative_tolerance),
        ("psd_absolute_tolerance", psd_absolute_tolerance),
        ("psd_relative_tolerance", psd_relative_tolerance),
    ):
        if not math.isfinite(value) or value < 0:
            raise ValidationInputError(f"{name} must be finite and nonnegative")

    with raw_path.open("r", encoding="utf-8", newline="") as stream:
        raw_epochs = list(iter_raw_epochs(stream))
    with difference_path.open("r", encoding="utf-8", newline="") as stream:
        difference_epochs = list(iter_difference_epochs(stream))

    raw_by_key: dict[tuple[int, Decimal], RawEpoch] = {}
    difference_by_key: dict[tuple[int, Decimal], DifferenceEpoch] = {}
    for epoch in raw_epochs:
        key = _epoch_key(epoch.gps_week, epoch.gps_tow_text)
        if key in raw_by_key:
            raise ValidationInputError(f"duplicate raw epoch {_key_text(key)}")
        raw_by_key[key] = epoch
    for epoch in difference_epochs:
        key = _epoch_key(epoch.gps_week, epoch.gps_tow_text)
        if key in difference_by_key:
            raise ValidationInputError(f"duplicate SD epoch {_key_text(key)}")
        difference_by_key[key] = epoch

    raw_keys = set(raw_by_key)
    difference_keys = set(difference_by_key)
    missing_sd = sorted(raw_keys - difference_keys)
    missing_raw = sorted(difference_keys - raw_keys)
    file_errors: list[str] = []
    if missing_sd:
        file_errors.append(
            f"epochs_missing_from_sd={[ _key_text(key) for key in missing_sd ]}"
        )
    if missing_raw:
        file_errors.append(
            f"epochs_missing_from_raw={[ _key_text(key) for key in missing_raw ]}"
        )
    if not raw_epochs:
        file_errors.append("raw_sidecar_has_no_epochs")
    if not difference_epochs:
        file_errors.append("sd_sidecar_has_no_epochs")

    epoch_reports = [
        validate_epoch_pair(
            raw_by_key[key],
            difference_by_key[key],
            absolute_tolerance,
            relative_tolerance,
            psd_absolute_tolerance,
            psd_relative_tolerance,
        )
        for key in sorted(raw_keys & difference_keys)
    ]
    valid_epoch_count = sum(bool(epoch["valid"]) for epoch in epoch_reports)
    all_valid = (
        not file_errors
        and bool(epoch_reports)
        and valid_epoch_count == len(epoch_reports)
        and len(raw_epochs) == len(difference_epochs) == len(epoch_reports)
    )
    return {
        "schema": REPORT_SCHEMA,
        "raw_schema": RAW_SCHEMA,
        "satellite_difference_schema": SD_SCHEMA,
        "raw_input": str(raw_path.resolve()),
        "satellite_difference_input": str(difference_path.resolve()),
        "tolerances": {
            "absolute": absolute_tolerance,
            "relative": relative_tolerance,
            "psd_absolute": psd_absolute_tolerance,
            "psd_relative": psd_relative_tolerance,
        },
        "raw_epoch_count": len(raw_epochs),
        "satellite_difference_epoch_count": len(difference_epochs),
        "matched_epoch_count": len(epoch_reports),
        "valid_epoch_count": valid_epoch_count,
        "invalid_epoch_count": len(epoch_reports) - valid_epoch_count,
        "all_valid": all_valid,
        "file_errors": file_errors,
        "epochs": epoch_reports,
    }


def write_json_atomic(path: Path, payload: dict[str, object]) -> None:
    path = path.resolve()
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary_name: str | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            dir=path.parent,
            prefix=path.name + ".",
            suffix=".tmp",
            delete=False,
        ) as stream:
            temporary_name = stream.name
            json.dump(
                payload,
                stream,
                ensure_ascii=False,
                indent=2,
                sort_keys=True,
                allow_nan=False,
            )
            stream.write("\n")
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary_name, path)
        temporary_name = None
    finally:
        if temporary_name is not None:
            Path(temporary_name).unlink(missing_ok=True)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("raw_v3", type=Path)
    parser.add_argument("satellite_difference_v1", type=Path)
    parser.add_argument("output_json", type=Path)
    parser.add_argument("--absolute-tolerance", type=float, default=1e-10)
    parser.add_argument("--relative-tolerance", type=float, default=1e-10)
    parser.add_argument("--psd-absolute-tolerance", type=float, default=1e-12)
    parser.add_argument("--psd-relative-tolerance", type=float, default=1e-10)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        report = validate_files(
            args.raw_v3,
            args.satellite_difference_v1,
            args.absolute_tolerance,
            args.relative_tolerance,
            args.psd_absolute_tolerance,
            args.psd_relative_tolerance,
        )
        status = 0 if report["all_valid"] else 1
    except (OSError, UnicodeError, csv.Error, ValidationInputError, ValueError) as error:
        report = {
            "schema": REPORT_SCHEMA,
            "raw_schema": RAW_SCHEMA,
            "satellite_difference_schema": SD_SCHEMA,
            "raw_input": str(args.raw_v3.resolve()),
            "satellite_difference_input": str(args.satellite_difference_v1.resolve()),
            "all_valid": False,
            "fatal_error": str(error),
            "epochs": [],
        }
        status = 2
    write_json_atomic(args.output_json, report)
    summary = {
        key: report.get(key)
        for key in (
            "raw_epoch_count",
            "satellite_difference_epoch_count",
            "matched_epoch_count",
            "valid_epoch_count",
            "invalid_epoch_count",
            "all_valid",
            "fatal_error",
        )
        if key in report
    }
    print(json.dumps(summary, ensure_ascii=False, indent=2, allow_nan=False))
    return status


if __name__ == "__main__":
    raise SystemExit(main())
