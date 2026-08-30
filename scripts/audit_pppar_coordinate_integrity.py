#!/usr/bin/env python3
"""Report coordinate-integrity evidence for a primary PPP-AR and FLOAT trace.

The audit retains the last *complete* REC_POS X/Y/Z vector for each
epoch/receiver, compares primary and FLOAT solutions, and compares both with
station coordinates in an IGS CRD SINEX file.  Optional thresholds are
screening limits only: neither a non-exceedance nor this report alone
establishes scientifically valid PPP-AR.
"""

from __future__ import annotations

import argparse
import bisect
import json
import math
import os
import re
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Iterable


Vector = tuple[float, float, float]
PositionKey = tuple[datetime, str]
PendingKey = tuple[str, datetime, str]
AXES = ("X", "Y", "Z")
OUTSIDE_STATES_BLOCK = "OUTSIDE_STATES_BLOCK"
AXIS_ALIASES = {"X": "X", "Y": "Y", "Z": "Z", "0": "X", "1": "Y", "2": "Z"}
CONSTRAINT_MARKER_RE = re.compile(r"^\s*fixAndHoldAmbiguities:\s+(.+?)\s*$")
APPLYING_RE = re.compile(r"^\s*Applying:")
AMBIGUITY_RECEIVER_RE = re.compile(r"A\(([^,\s]+),")
PSEUDOOBS_DESIGN_RE = re.compile(
    r"PPP_AR\s+PSEUDOOBS_DESIGN\s+rows=(\d+).*?status=([A-Z0-9_]+)"
)
PSEUDOOBS_NOISE_RE = re.compile(
    r"PPP_AR\s+PSEUDOOBS_NOISE\s+rows=(\d+).*?status=([A-Z0-9_]+)"
)
PSEUDOOBS_SUBMISSION_RE = re.compile(
    r"PPP_AR\s+PSEUDOOBS_SUBMISSION\s+rows=(\d+).*?"
    r"status=FILTER_CALL_RETURNED_SUBMITTED_UNVERIFIED"
)


@dataclass(frozen=True)
class PositionSample:
    epoch: datetime
    receiver: str
    ecef_m: Vector
    completed_at_line: int
    source_state_block: str | None


@dataclass(frozen=True)
class SinexReference:
    receiver: str
    point: str
    solution: str
    reference_epoch: datetime | None
    valid_from: datetime | None
    valid_to: datetime | None
    ecef_m: Vector


@dataclass(frozen=True)
class SinexSelection:
    reference: SinexReference | None
    status: str
    candidate_count: int
    reference_age_days: float | None


@dataclass(frozen=True)
class Thresholds:
    ar_float_m: float | None = None
    sinex_m: float | None = None
    jump_m: float | None = None


@dataclass
class TraceParseResult:
    samples: dict[PositionKey, PositionSample]
    pseudoobs_attempt_rows_by_epoch: Counter[datetime]
    pseudoobs_attempt_receiver_rows_by_epoch: dict[datetime, Counter[str]]
    submitted_unverified_rows_by_epoch: Counter[datetime]
    metadata: dict[str, object]


@dataclass
class SinexParseResult:
    references_by_receiver: dict[str, list[SinexReference]]
    metadata: dict[str, object]


def parse_trace_epoch(text: str) -> datetime:
    """Parse the UTC-like epoch format emitted by Ginan trace output."""

    value = text.strip()
    try:
        epoch = datetime.fromisoformat(value.replace("Z", "+00:00"))
    except ValueError as iso_error:
        # Python 3.10's fromisoformat is stricter than newer Python releases
        # for the two-digit fractional seconds emitted by this Ginan build.
        # strptime's %f accepts one through six digits and keeps the WSL
        # experiment runtime consistent with the Windows unit-test runtime.
        normalized = value.removesuffix("Z")
        epoch = None
        for pattern in ("%Y-%m-%d %H:%M:%S.%f", "%Y-%m-%d %H:%M:%S"):
            try:
                epoch = datetime.strptime(normalized, pattern)
                break
            except ValueError:
                continue
        if epoch is None:
            raise ValueError(f"unsupported Ginan trace epoch {value!r}") from iso_error
    if epoch.tzinfo is None:
        epoch = epoch.replace(tzinfo=timezone.utc)
    return epoch.astimezone(timezone.utc)


def epoch_text(epoch: datetime) -> str:
    return epoch.astimezone(timezone.utc).isoformat(timespec="microseconds").replace(
        "+00:00", "Z"
    )


def finite_float(text: str) -> float:
    value = float(text.replace("D", "E").replace("d", "e"))
    if not math.isfinite(value):
        raise ValueError("non-finite numeric value")
    return value


def parse_rec_pos_row(line: str) -> tuple[datetime, str, str, float] | None:
    """Parse one tab-delimited KF state row, or return None for other rows."""

    fields = line.rstrip("\r\n").split("\t")
    if not fields or fields[0].strip() != "*":
        return None
    try:
        type_index = next(
            index for index, field in enumerate(fields) if field.strip() == "REC_POS"
        )
    except StopIteration:
        return None
    # KFKey renders as type, satellite, receiver, coordinate.  REC_POS has an
    # empty satellite field, which is preserved by splitting on tabs.
    if len(fields) <= type_index + 4 or len(fields) <= 2:
        raise ValueError("malformed REC_POS row")
    receiver = fields[type_index + 2].strip().upper()
    axis = AXIS_ALIASES.get(fields[type_index + 3].strip().upper())
    if not receiver or axis is None:
        raise ValueError("malformed REC_POS receiver or coordinate")
    return (
        parse_trace_epoch(fields[2]),
        receiver,
        axis,
        finite_float(fields[type_index + 4].strip()),
    )


def _fragment(
    key: PendingKey,
    components: dict[str, float],
    reason: str,
) -> dict[str, object]:
    return {
        "state_block": key[0],
        "receiver": key[2],
        "epoch": epoch_text(key[1]),
        "present_components": sorted(components),
        "missing_components": sorted(set(AXES) - set(components)),
        "reason": reason,
    }


def normalize_state_block(name: str) -> str:
    value = name.strip().lstrip("+").lstrip("-")
    if value == OUTSIDE_STATES_BLOCK:
        return value
    return value if value.startswith("STATES/") else f"STATES/{value}"


def select_state_block(
    available_blocks: set[str], role: str, requested_state_block: str | None
) -> tuple[str | None, str]:
    if requested_state_block is not None:
        requested = normalize_state_block(requested_state_block)
        if requested in available_blocks:
            return requested, "selected_explicitly"
        return None, "requested_state_block_missing"
    expected = {
        "primary": {"STATES/AR", "STATES/AR_RTS"},
        "float": {"STATES/PPP", "STATES/PPP_RTS"},
    }[role]
    candidates = sorted(available_blocks & expected)
    if len(candidates) == 1:
        return candidates[0], "selected_unique_role_appropriate_block"
    if not candidates:
        return None, "role_appropriate_state_block_missing"
    return None, "ambiguous_multiple_role_appropriate_state_blocks"


def parse_trace(
    path: Path,
    role: str,
    requested_state_block: str | None = None,
) -> TraceParseResult:
    """Parse REC_POS vectors and epochs where integer constraints were applied.

    A second complete X/Y/Z group for the same epoch/receiver replaces the
    earlier group.  A repeated component before a group is complete starts a
    new group, so coordinates from two partial groups are never mixed.
    """

    if role not in {"primary", "float"}:
        raise ValueError("trace role must be 'primary' or 'float'")
    samples_by_block: dict[str, dict[PositionKey, PositionSample]] = defaultdict(dict)
    pending: dict[PendingKey, dict[str, float]] = {}
    complete_occurrences: Counter[PendingKey] = Counter()
    pseudoobs_attempt_rows_by_epoch: Counter[datetime] = Counter()
    pseudoobs_attempt_receiver_rows_by_epoch: dict[datetime, Counter[str]] = defaultdict(
        Counter
    )
    design_validated_rows_by_epoch: Counter[datetime] = Counter()
    noise_validated_rows_by_epoch: Counter[datetime] = Counter()
    submitted_unverified_rows_by_epoch: Counter[datetime] = Counter()
    fragments: list[dict[str, object]] = []
    component_rows = 0
    malformed_rows = 0
    current_constraint_epoch: datetime | None = None
    constraint_attempt_markers = 0
    current_state_block: str | None = None

    with path.open("r", encoding="utf-8", errors="replace") as stream:
        for line_number, line in enumerate(stream, start=1):
            stripped = line.strip()
            if stripped.startswith("+STATES/"):
                current_state_block = stripped[1:]
                continue
            if stripped.startswith("-STATES/"):
                closing_block = stripped[1:]
                for pending_key, components in pending.items():
                    if pending_key[0] == closing_block and components:
                        fragments.append(
                            _fragment(
                                pending_key,
                                components,
                                "state_block_closed_before_complete",
                            )
                        )
                        pending[pending_key] = {}
                current_state_block = None
                continue
            marker = CONSTRAINT_MARKER_RE.match(line)
            if marker:
                current_constraint_epoch = parse_trace_epoch(marker.group(1))
                constraint_attempt_markers += 1
                continue
            if APPLYING_RE.match(line) and current_constraint_epoch is not None:
                receivers = set(AMBIGUITY_RECEIVER_RE.findall(line))
                if not receivers:
                    continue
                pseudoobs_attempt_rows_by_epoch[current_constraint_epoch] += 1
                for receiver in receivers:
                    pseudoobs_attempt_receiver_rows_by_epoch[current_constraint_epoch][
                        receiver.upper()
                    ] += 1

            design = PSEUDOOBS_DESIGN_RE.search(line)
            if design and current_constraint_epoch is not None:
                if design.group(2) == "MATCHES_Z_TIMES_D":
                    design_validated_rows_by_epoch[current_constraint_epoch] += int(
                        design.group(1)
                    )
                else:
                    current_constraint_epoch = None
            noise = PSEUDOOBS_NOISE_RE.search(line)
            if noise and current_constraint_epoch is not None:
                if noise.group(2) == "INDEPENDENT_DIAGONAL":
                    noise_validated_rows_by_epoch[current_constraint_epoch] += int(
                        noise.group(1)
                    )
                else:
                    current_constraint_epoch = None
            submission = PSEUDOOBS_SUBMISSION_RE.search(line)
            if submission and current_constraint_epoch is not None:
                submitted_unverified_rows_by_epoch[current_constraint_epoch] += int(
                    submission.group(1)
                )
                current_constraint_epoch = None

            if "REC_POS" not in line:
                continue
            try:
                parsed = parse_rec_pos_row(line)
            except ValueError:
                malformed_rows += 1
                continue
            if parsed is None:
                continue
            epoch, receiver, axis, value = parsed
            component_rows += 1
            state_block = current_state_block or OUTSIDE_STATES_BLOCK
            pending_key = (state_block, epoch, receiver)
            components = pending.setdefault(pending_key, {})
            if axis in components:
                fragments.append(
                    _fragment(
                        pending_key,
                        components,
                        "component_repeated_before_complete",
                    )
                )
                components = {}
                pending[pending_key] = components
            components[axis] = value
            if set(components) != set(AXES):
                continue
            complete_occurrences[pending_key] += 1
            samples_by_block[state_block][(epoch, receiver)] = PositionSample(
                epoch=epoch,
                receiver=receiver,
                ecef_m=(components["X"], components["Y"], components["Z"]),
                completed_at_line=line_number,
                source_state_block=current_state_block,
            )
            pending[pending_key] = {}

    for key, components in sorted(pending.items()):
        if components:
            fragments.append(_fragment(key, components, "end_of_file_before_complete"))

    available_blocks = set(samples_by_block)
    selected_block, selection_status = select_state_block(
        available_blocks, role, requested_state_block
    )
    samples = samples_by_block.get(selected_block, {}) if selected_block else {}
    selected_fragments = [
        fragment for fragment in fragments if fragment["state_block"] == selected_block
    ]
    duplicate_groups = sum(
        max(0, count - 1)
        for key, count in complete_occurrences.items()
        if key[0] == selected_block
    )
    selected_complete_occurrences = sum(
        count
        for key, count in complete_occurrences.items()
        if key[0] == selected_block
    )
    return TraceParseResult(
        samples=samples,
        pseudoobs_attempt_rows_by_epoch=pseudoobs_attempt_rows_by_epoch,
        pseudoobs_attempt_receiver_rows_by_epoch=dict(
            pseudoobs_attempt_receiver_rows_by_epoch
        ),
        submitted_unverified_rows_by_epoch=submitted_unverified_rows_by_epoch,
        metadata={
            "path": path.as_posix(),
            "role": role,
            "state_block_selection_status": selection_status,
            "selected_state_block": selected_block,
            "requested_state_block": (
                normalize_state_block(requested_state_block)
                if requested_state_block is not None
                else None
            ),
            "available_state_blocks": {
                block: len(block_samples)
                for block, block_samples in sorted(samples_by_block.items())
            },
            "rec_pos_component_row_count": component_rows,
            "malformed_rec_pos_row_count": malformed_rows,
            "complete_vector_occurrence_count_all_blocks": sum(
                complete_occurrences.values()
            ),
            "complete_vector_occurrence_count": selected_complete_occurrences,
            "retained_epoch_receiver_vector_count": len(samples),
            "duplicate_complete_vector_count": duplicate_groups,
            "discarded_incomplete_fragment_count": len(selected_fragments),
            "discarded_incomplete_fragments": selected_fragments,
            "discarded_incomplete_fragments_all_blocks_count": len(fragments),
            "constraint_attempt_marker_count": constraint_attempt_markers,
            "visible_integer_pseudoobservation_attempt_row_count": sum(
                pseudoobs_attempt_rows_by_epoch.values()
            ),
            "design_validated_pseudoobservation_row_count": sum(
                design_validated_rows_by_epoch.values()
            ),
            "noise_validated_pseudoobservation_row_count": sum(
                noise_validated_rows_by_epoch.values()
            ),
            "submitted_unverified_pseudoobservation_row_count": sum(
                submitted_unverified_rows_by_epoch.values()
            ),
            "duplicate_selection_rule": (
                "within the selected state block, retain the last complete X/Y/Z "
                "group in file order for each epoch/receiver; never merge across "
                "a repeated component or across state blocks"
            ),
        },
    )


def parse_sinex_epoch(text: str) -> datetime | None:
    fields = text.split(":")
    if len(fields) != 3:
        raise ValueError(f"unsupported SINEX epoch {text!r}")
    year, day_of_year, seconds = (int(field) for field in fields)
    if day_of_year == 0:
        return None
    # Match the SINEX 2.02 convention used by Ginan's own parser.
    full_year = 1900 + year if year > 50 else 2000 + year
    return datetime(full_year, 1, 1, tzinfo=timezone.utc) + timedelta(
        days=day_of_year - 1, seconds=seconds
    )


def parse_sinex(path: Path) -> SinexParseResult:
    """Parse complete STAX/STAY/STAZ records from SOLUTION/ESTIMATE."""

    components: dict[tuple[str, str, str, datetime | None], dict[str, float]] = {}
    duplicate_component_rows = 0
    estimate_rows = 0
    malformed_rows = 0
    section: str | None = None
    validity_by_solution: dict[
        tuple[str, str, str], tuple[datetime | None, datetime | None]
    ] = {}
    malformed_epoch_rows = 0
    with path.open("r", encoding="ascii", errors="replace") as stream:
        for line in stream:
            stripped = line.strip()
            if stripped == "+SOLUTION/ESTIMATE":
                section = "ESTIMATE"
                continue
            if stripped == "+SOLUTION/EPOCHS":
                section = "EPOCHS"
                continue
            if stripped in {"-SOLUTION/ESTIMATE", "-SOLUTION/EPOCHS"}:
                section = None
                continue
            if section is None or not stripped or stripped.startswith("*"):
                continue
            fields = stripped.split()
            if section == "EPOCHS":
                if len(fields) < 7:
                    malformed_epoch_rows += 1
                    continue
                try:
                    key = (fields[0].upper(), fields[1], fields[2])
                    validity_by_solution[key] = (
                        parse_sinex_epoch(fields[4]),
                        parse_sinex_epoch(fields[5]),
                    )
                except (ValueError, IndexError):
                    malformed_epoch_rows += 1
                continue
            if len(fields) < 10 or fields[1] not in {"STAX", "STAY", "STAZ"}:
                continue
            estimate_rows += 1
            try:
                receiver = fields[2].upper()
                point = fields[3]
                solution = fields[4]
                reference_epoch = parse_sinex_epoch(fields[5])
                value = finite_float(fields[8])
            except (ValueError, IndexError):
                malformed_rows += 1
                continue
            key = (receiver, point, solution, reference_epoch)
            axis = fields[1][-1]
            record = components.setdefault(key, {})
            if axis in record:
                duplicate_component_rows += 1
            record[axis] = value

    references_by_receiver: dict[str, list[SinexReference]] = defaultdict(list)
    incomplete: list[dict[str, object]] = []
    for key, record in sorted(
        components.items(),
        key=lambda item: (
            item[0][0],
            item[0][1],
            item[0][2],
            item[0][3] or datetime.min.replace(tzinfo=timezone.utc),
        ),
    ):
        receiver, point, solution, reference_epoch = key
        if set(record) != set(AXES):
            incomplete.append(
                {
                    "receiver": receiver,
                    "point": point,
                    "solution": solution,
                    "reference_epoch": (
                        epoch_text(reference_epoch) if reference_epoch else None
                    ),
                    "missing_components": sorted(set(AXES) - set(record)),
                }
            )
            continue
        valid_from, valid_to = validity_by_solution.get(
            (receiver, point, solution), (None, None)
        )
        references_by_receiver[receiver].append(
            SinexReference(
                receiver=receiver,
                point=point,
                solution=solution,
                reference_epoch=reference_epoch,
                valid_from=valid_from,
                valid_to=valid_to,
                ecef_m=(record["X"], record["Y"], record["Z"]),
            )
        )

    return SinexParseResult(
        references_by_receiver=dict(references_by_receiver),
        metadata={
            "path": path.as_posix(),
            "station_coordinate_estimate_row_count": estimate_rows,
            "malformed_station_coordinate_row_count": malformed_rows,
            "solution_epoch_record_count": len(validity_by_solution),
            "malformed_solution_epoch_row_count": malformed_epoch_rows,
            "duplicate_station_component_row_count": duplicate_component_rows,
            "complete_reference_count": sum(map(len, references_by_receiver.values())),
            "incomplete_reference_count": len(incomplete),
            "incomplete_references": incomplete,
            "reference_selection_rule": (
                "require exactly one SOLUTION/EPOCHS validity interval containing the "
                "trace epoch, then use that solution's coordinate reference epoch"
            ),
        },
    )


def select_sinex_reference(
    references: list[SinexReference],
    epoch: datetime,
    max_reference_age_days: float | None = None,
) -> SinexSelection:
    active = [
        reference
        for reference in references
        if reference.valid_from is not None
        and reference.valid_to is not None
        and reference.valid_from <= epoch <= reference.valid_to
    ]
    if not active:
        return SinexSelection(None, "no_solution_valid_at_trace_epoch", len(references), None)
    if len(active) > 1:
        return SinexSelection(None, "multiple_solutions_valid_at_trace_epoch", len(active), None)
    reference = active[0]
    if reference.reference_epoch is None:
        return SinexSelection(None, "valid_solution_reference_epoch_unknown", 1, None)
    age_days = abs((reference.reference_epoch - epoch).total_seconds()) / 86400.0
    if max_reference_age_days is not None and age_days > max_reference_age_days:
        return SinexSelection(None, "reference_age_exceeds_configured_limit", 1, age_days)
    return SinexSelection(reference, "unique_time_valid_solution", 1, age_days)


def subtract(left: Vector, right: Vector) -> Vector:
    return tuple(a - b for a, b in zip(left, right))  # type: ignore[return-value]


def norm(vector: Vector) -> float:
    return math.sqrt(sum(component * component for component in vector))


def ecef_delta_to_enu(delta: Vector, reference: Vector) -> Vector | None:
    """Rotate an ECEF delta to local ENU at a WGS-84 reference position."""

    x, y, z = reference
    if norm(reference) == 0:
        return None
    semi_major = 6378137.0
    flattening = 1 / 298.257223563
    semi_minor = semi_major * (1 - flattening)
    eccentricity_sq = flattening * (2 - flattening)
    second_eccentricity_sq = (semi_major**2 - semi_minor**2) / semi_minor**2
    horizontal = math.hypot(x, y)
    longitude = math.atan2(y, x)
    if horizontal == 0:
        latitude = math.copysign(math.pi / 2, z)
    else:
        theta = math.atan2(z * semi_major, horizontal * semi_minor)
        latitude = math.atan2(
            z + second_eccentricity_sq * semi_minor * math.sin(theta) ** 3,
            horizontal - eccentricity_sq * semi_major * math.cos(theta) ** 3,
        )
    sin_lat, cos_lat = math.sin(latitude), math.cos(latitude)
    sin_lon, cos_lon = math.sin(longitude), math.cos(longitude)
    dx, dy, dz = delta
    return (
        -sin_lon * dx + cos_lon * dy,
        -sin_lat * cos_lon * dx - sin_lat * sin_lon * dy + cos_lat * dz,
        cos_lat * cos_lon * dx + cos_lat * sin_lon * dy + sin_lat * dz,
    )


def vector_json(vector: Vector) -> dict[str, float]:
    return {axis.lower(): value for axis, value in zip(AXES, vector)}


def enu_json(vector: Vector) -> dict[str, float]:
    return {axis: value for axis, value in zip(("east", "north", "up"), vector)}


def comparison_json(
    estimate: Vector,
    reference: Vector,
    threshold: float | None,
    include_enu: bool,
) -> dict[str, object]:
    delta = subtract(estimate, reference)
    distance = norm(delta)
    result: dict[str, object] = {
        "ecef_delta_m": vector_json(delta),
        "distance_3d_m": distance,
        "threshold_m": threshold,
        "exceeds_threshold": distance > threshold if threshold is not None else None,
    }
    if include_enu:
        enu = ecef_delta_to_enu(delta, reference)
        result["enu_delta_m"] = enu_json(enu) if enu is not None else None
    return result


def statistics(values: Iterable[float]) -> dict[str, float | int | None]:
    materialized = list(values)
    if not materialized:
        return {"count": 0, "maximum_m": None, "rms_m": None}
    return {
        "count": len(materialized),
        "maximum_m": max(materialized),
        "rms_m": math.sqrt(sum(value * value for value in materialized) / len(materialized)),
    }


def _key_json(key: PositionKey) -> dict[str, str]:
    return {"receiver": key[1], "epoch": epoch_text(key[0])}


def _metric_value(record: dict[str, object] | None) -> float | None:
    if record is None:
        return None
    return float(record["distance_3d_m"])


def audit_coordinate_integrity(
    primary_trace: Path,
    float_trace: Path,
    sinex_path: Path,
    thresholds: Thresholds | None = None,
    primary_state_block: str | None = None,
    float_state_block: str | None = None,
    max_sinex_reference_age_days: float | None = None,
) -> dict[str, object]:
    thresholds = thresholds or Thresholds()
    primary = parse_trace(primary_trace, "primary", primary_state_block)
    floating = parse_trace(float_trace, "float", float_state_block)
    sinex = parse_sinex(sinex_path)

    all_keys = sorted(set(primary.samples) | set(floating.samples))
    comparisons: list[dict[str, object]] = []
    metric_values: dict[str, list[float]] = defaultdict(list)
    receiver_metric_values: dict[str, dict[str, list[float]]] = defaultdict(
        lambda: defaultdict(list)
    )
    attempt_epochs = sorted(primary.pseudoobs_attempt_rows_by_epoch)
    submitted_epochs = sorted(primary.submitted_unverified_rows_by_epoch)
    receiver_attempt_epochs: dict[str, list[datetime]] = defaultdict(list)
    for attempt_epoch, receiver_counts in primary.pseudoobs_attempt_receiver_rows_by_epoch.items():
        for receiver in receiver_counts:
            receiver_attempt_epochs[receiver].append(attempt_epoch)
    for epochs in receiver_attempt_epochs.values():
        epochs.sort()
    receiver_submitted_epochs = {
        receiver: [
            epoch
            for epoch in epochs
            if epoch in primary.submitted_unverified_rows_by_epoch
        ]
        for receiver, epochs in receiver_attempt_epochs.items()
    }
    exceedance_counts: Counter[str] = Counter()

    for key in all_keys:
        epoch, receiver = key
        primary_sample = primary.samples.get(key)
        float_sample = floating.samples.get(key)
        references = sinex.references_by_receiver.get(receiver, [])
        selection = select_sinex_reference(
            references,
            epoch,
            max_reference_age_days=max_sinex_reference_age_days,
        )
        reference = selection.reference
        missing: list[str] = []
        if primary_sample is None:
            missing.append("primary_trace")
        if float_sample is None:
            missing.append("float_trace")
        if reference is None:
            missing.append("sinex_reference")

        ar_vs_float = None
        if primary_sample is not None and float_sample is not None:
            ar_vs_float = comparison_json(
                primary_sample.ecef_m,
                float_sample.ecef_m,
                thresholds.ar_float_m,
                include_enu=False,
            )

        primary_vs_sinex = None
        float_vs_sinex = None
        if reference is not None and primary_sample is not None:
            primary_vs_sinex = comparison_json(
                primary_sample.ecef_m,
                reference.ecef_m,
                thresholds.sinex_m,
                include_enu=True,
            )
        if reference is not None and float_sample is not None:
            float_vs_sinex = comparison_json(
                float_sample.ecef_m,
                reference.ecef_m,
                thresholds.sinex_m,
                include_enu=True,
            )

        for metric, result in (
            ("ar_vs_float_3d", ar_vs_float),
            ("primary_vs_sinex_3d", primary_vs_sinex),
            ("float_vs_sinex_3d", float_vs_sinex),
        ):
            value = _metric_value(result)
            if value is not None:
                metric_values[metric].append(value)
                receiver_metric_values[receiver][metric].append(value)
            if result is not None and result["exceeds_threshold"] is True:
                exceedance_counts[metric] += 1

        attempt_index = bisect.bisect_right(attempt_epochs, epoch)
        latest_any_attempt = (
            attempt_epochs[attempt_index - 1] if attempt_index else None
        )
        station_attempt_epochs = receiver_attempt_epochs.get(receiver, [])
        station_attempt_index = bisect.bisect_right(station_attempt_epochs, epoch)
        latest_station_attempt = (
            station_attempt_epochs[station_attempt_index - 1]
            if station_attempt_index
            else None
        )
        station_submitted_epochs = receiver_submitted_epochs.get(receiver, [])
        station_submitted_index = bisect.bisect_right(station_submitted_epochs, epoch)
        latest_station_submission = (
            station_submitted_epochs[station_submitted_index - 1]
            if station_submitted_index
            else None
        )
        reference_json = None
        if reference is not None:
            reference_json = {
                "receiver": reference.receiver,
                "point": reference.point,
                "solution": reference.solution,
                "reference_epoch": (
                    epoch_text(reference.reference_epoch)
                    if reference.reference_epoch is not None
                    else None
                ),
                "valid_from": (
                    epoch_text(reference.valid_from) if reference.valid_from else None
                ),
                "valid_to": epoch_text(reference.valid_to) if reference.valid_to else None,
                "ecef_m": vector_json(reference.ecef_m),
            }
        comparisons.append(
            {
                "receiver": receiver,
                "epoch": epoch_text(epoch),
                "primary_ecef_m": (
                    vector_json(primary_sample.ecef_m) if primary_sample else None
                ),
                "primary_source": (
                    {
                        "state_block": primary_sample.source_state_block,
                        "completed_at_line": primary_sample.completed_at_line,
                    }
                    if primary_sample
                    else None
                ),
                "float_ecef_m": (
                    vector_json(float_sample.ecef_m) if float_sample else None
                ),
                "float_source": (
                    {
                        "state_block": float_sample.source_state_block,
                        "completed_at_line": float_sample.completed_at_line,
                    }
                    if float_sample
                    else None
                ),
                "sinex_reference": reference_json,
                "sinex_reference_selection": {
                    "status": selection.status,
                    "candidate_count": selection.candidate_count,
                    "reference_age_days": selection.reference_age_days,
                    "maximum_reference_age_days": max_sinex_reference_age_days,
                },
                "ar_vs_float": ar_vs_float,
                "primary_vs_sinex": primary_vs_sinex,
                "float_vs_sinex": float_vs_sinex,
                "any_integer_pseudoobservation_attempt_at_epoch": (
                    primary.pseudoobs_attempt_rows_by_epoch.get(epoch, 0) > 0
                ),
                "integer_pseudoobservation_attempt_for_receiver_at_epoch": (
                    primary.pseudoobs_attempt_receiver_rows_by_epoch.get(epoch, {}).get(
                        receiver, 0
                    )
                    > 0
                ),
                "integer_pseudoobservation_attempt_rows_for_receiver_at_epoch": (
                    primary.pseudoobs_attempt_receiver_rows_by_epoch.get(epoch, {}).get(
                        receiver, 0
                    )
                ),
                "submitted_unverified_at_epoch": (
                    primary.submitted_unverified_rows_by_epoch.get(epoch, 0) > 0
                ),
                "submitted_unverified_for_receiver_at_epoch": (
                    primary.submitted_unverified_rows_by_epoch.get(epoch, 0) > 0
                    and primary.pseudoobs_attempt_receiver_rows_by_epoch.get(
                        epoch, {}
                    ).get(receiver, 0)
                    > 0
                ),
                "latest_any_integer_pseudoobservation_attempt_epoch": (
                    epoch_text(latest_any_attempt) if latest_any_attempt else None
                ),
                "latest_integer_pseudoobservation_attempt_epoch_for_receiver": (
                    epoch_text(latest_station_attempt) if latest_station_attempt else None
                ),
                "latest_submitted_unverified_epoch_for_receiver": (
                    epoch_text(latest_station_submission)
                    if latest_station_submission
                    else None
                ),
                "missing": missing,
            }
        )

    jumps: list[dict[str, object]] = []
    primary_by_receiver: dict[str, list[PositionSample]] = defaultdict(list)
    for sample in primary.samples.values():
        primary_by_receiver[sample.receiver].append(sample)
    for receiver, samples in sorted(primary_by_receiver.items()):
        samples.sort(key=lambda sample: sample.epoch)
        for previous, current in zip(samples, samples[1:]):
            delta = subtract(current.ecef_m, previous.ecef_m)
            distance = norm(delta)
            lower = bisect.bisect_right(attempt_epochs, previous.epoch)
            upper = bisect.bisect_right(attempt_epochs, current.epoch)
            interval_attempts = attempt_epochs[lower:upper]
            station_attempt_epochs = receiver_attempt_epochs.get(receiver, [])
            station_lower = bisect.bisect_right(
                station_attempt_epochs, previous.epoch
            )
            station_upper = bisect.bisect_right(station_attempt_epochs, current.epoch)
            station_interval_attempts = station_attempt_epochs[station_lower:station_upper]
            station_submitted_epochs = receiver_submitted_epochs.get(receiver, [])
            submitted_lower = bisect.bisect_right(
                station_submitted_epochs, previous.epoch
            )
            submitted_upper = bisect.bisect_right(station_submitted_epochs, current.epoch)
            station_interval_submissions = station_submitted_epochs[
                submitted_lower:submitted_upper
            ]
            exceeds = distance > thresholds.jump_m if thresholds.jump_m is not None else None
            jump = {
                "receiver": receiver,
                "previous_epoch": epoch_text(previous.epoch),
                "current_epoch": epoch_text(current.epoch),
                "interval_seconds": (current.epoch - previous.epoch).total_seconds(),
                "ecef_delta_m": vector_json(delta),
                "distance_3d_m": distance,
                "integer_pseudoobservation_attempts_in_interval": [
                    {
                        "epoch": epoch_text(event),
                        "attempt_row_count": primary.pseudoobs_attempt_rows_by_epoch[event],
                    }
                    for event in interval_attempts
                ],
                "receiver_attempts_in_interval": [
                    {
                        "epoch": epoch_text(event),
                        "attempt_row_count": (
                            primary.pseudoobs_attempt_receiver_rows_by_epoch[event][receiver]
                        ),
                    }
                    for event in station_interval_attempts
                ],
                "receiver_submitted_unverified_events_in_interval": [
                    {
                        "epoch": epoch_text(event),
                        "submitted_unverified_row_count": (
                            primary.submitted_unverified_rows_by_epoch[event]
                        ),
                    }
                    for event in station_interval_submissions
                ],
                "jump_threshold_m": thresholds.jump_m,
                "exceeds_jump_threshold": exceeds,
            }
            jumps.append(jump)
            metric_values["primary_epoch_jump_3d"].append(distance)
            receiver_metric_values[receiver]["primary_epoch_jump_3d"].append(distance)
            if exceeds is True:
                exceedance_counts["primary_epoch_jump_3d"] += 1
                if station_interval_attempts:
                    exceedance_counts[
                        "receiver_attempt_interval_primary_epoch_jump_3d"
                    ] += 1
                if station_interval_submissions:
                    exceedance_counts[
                        "submitted_unverified_interval_primary_epoch_jump_3d"
                    ] += 1

    metrics = (
        "ar_vs_float_3d",
        "primary_vs_sinex_3d",
        "float_vs_sinex_3d",
        "primary_epoch_jump_3d",
    )
    receivers = sorted({key[1] for key in all_keys})
    configured = any(
        threshold is not None
        for threshold in (thresholds.ar_float_m, thresholds.sinex_m, thresholds.jump_m)
    )
    evaluated_counts = {
        metric: len(metric_values[metric])
        for metric in metrics
    }
    primary_keys = set(primary.samples)
    float_keys = set(floating.samples)
    trace_receivers = {key[1] for key in all_keys}
    common_trace_key_count = len(primary_keys & float_keys)
    trace_key_union_count = len(primary_keys | float_keys)
    comparison_coverage = (
        common_trace_key_count / trace_key_union_count if trace_key_union_count else None
    )
    sinex_selection_failures = Counter(
        row["sinex_reference_selection"]["status"]
        for row in comparisons
        if row["sinex_reference_selection"]["status"] != "unique_time_valid_solution"
    )
    data_quality_issues: list[dict[str, object]] = []

    def add_quality_issue(code: str, count: int = 1) -> None:
        if count:
            data_quality_issues.append({"code": code, "count": count})

    if not str(primary.metadata["state_block_selection_status"]).startswith("selected_"):
        add_quality_issue("primary_state_block_not_selected")
    if not str(floating.metadata["state_block_selection_status"]).startswith("selected_"):
        add_quality_issue("float_state_block_not_selected")
    add_quality_issue(
        "primary_malformed_rec_pos_rows",
        int(primary.metadata["malformed_rec_pos_row_count"]),
    )
    add_quality_issue(
        "float_malformed_rec_pos_rows",
        int(floating.metadata["malformed_rec_pos_row_count"]),
    )
    add_quality_issue(
        "primary_incomplete_rec_pos_fragments",
        int(primary.metadata["discarded_incomplete_fragment_count"]),
    )
    add_quality_issue(
        "float_incomplete_rec_pos_fragments",
        int(floating.metadata["discarded_incomplete_fragment_count"]),
    )
    add_quality_issue("trace_keys_missing_primary", len(float_keys - primary_keys))
    add_quality_issue("trace_keys_missing_float", len(primary_keys - float_keys))
    add_quality_issue(
        "sinex_malformed_station_coordinate_rows",
        int(sinex.metadata["malformed_station_coordinate_row_count"]),
    )
    add_quality_issue(
        "sinex_malformed_solution_epoch_rows",
        int(sinex.metadata["malformed_solution_epoch_row_count"]),
    )
    add_quality_issue(
        "sinex_incomplete_references",
        int(sinex.metadata["incomplete_reference_count"]),
    )
    for status, count in sorted(sinex_selection_failures.items()):
        add_quality_issue(f"sinex_selection_{status}", count)
    if attempt_row_count := int(
        primary.metadata["visible_integer_pseudoobservation_attempt_row_count"]
    ):
        submitted_count = int(
            primary.metadata["submitted_unverified_pseudoobservation_row_count"]
        )
        if submitted_count != attempt_row_count:
            add_quality_issue(
                "pseudoobservation_attempt_submission_row_count_mismatch",
                abs(attempt_row_count - submitted_count),
            )
    configured_without_values: list[str] = []
    if thresholds.ar_float_m is not None and not evaluated_counts["ar_vs_float_3d"]:
        configured_without_values.append("ar_vs_float_3d")
    if thresholds.sinex_m is not None:
        if not evaluated_counts["primary_vs_sinex_3d"]:
            configured_without_values.append("primary_vs_sinex_3d")
        if not evaluated_counts["float_vs_sinex_3d"]:
            configured_without_values.append("float_vs_sinex_3d")
    if thresholds.jump_m is not None and not evaluated_counts["primary_epoch_jump_3d"]:
        configured_without_values.append("primary_epoch_jump_3d")

    configured_metric_count = (
        int(thresholds.ar_float_m is not None)
        + 2 * int(thresholds.sinex_m is not None)
        + int(thresholds.jump_m is not None)
    )
    if not configured:
        screening_status = "not_evaluated_no_thresholds_configured"
    elif len(configured_without_values) == configured_metric_count:
        screening_status = "not_evaluated_no_comparable_values"
    elif exceedance_counts:
        screening_status = "configured_threshold_exceedance_observed"
    elif configured_without_values or data_quality_issues:
        screening_status = "no_exceedance_in_evaluated_values_data_quality_issues_present"
    else:
        screening_status = "no_configured_threshold_exceedance_observed"

    marker_count = int(primary.metadata["constraint_attempt_marker_count"])
    attempt_row_count = int(
        primary.metadata["visible_integer_pseudoobservation_attempt_row_count"]
    )
    submitted_row_count = int(
        primary.metadata["submitted_unverified_pseudoobservation_row_count"]
    )
    if submitted_row_count:
        pseudoobs_detection_status = "submitted_unverified_telemetry_observed"
    elif attempt_row_count:
        pseudoobs_detection_status = "attempt_rows_observed_without_submission_telemetry"
    elif marker_count:
        pseudoobs_detection_status = "attempt_markers_observed_no_visible_attempt_rows"
    else:
        pseudoobs_detection_status = "no_attempt_markers_observed"

    return {
        "schema": "GINAN_PPPAR_COORDINATE_INTEGRITY_AUDIT_V1",
        "inputs": {
            "primary_trace": primary_trace.as_posix(),
            "float_trace": float_trace.as_posix(),
            "igs_crd_sinex": sinex_path.as_posix(),
            "primary_state_block": primary.metadata["selected_state_block"],
            "float_state_block": floating.metadata["selected_state_block"],
            "maximum_sinex_reference_age_days": max_sinex_reference_age_days,
        },
        "definitions": {
            "ar_vs_float_delta": "primary ECEF minus FLOAT ECEF",
            "sinex_delta": "trace estimate ECEF minus selected SINEX ECEF",
            "trace_epoch_matching": "exact receiver and normalized UTC epoch",
            "local_frame": "ENU rotation at the selected SINEX ECEF using WGS-84",
            "jump_interval": "consecutive retained primary epochs for one receiver",
            "sinex_epoch_handling": (
                "published SINEX coordinates are used directly; no station-velocity "
                "propagation is inferred"
            ),
        },
        "parsing": {
            "primary_trace": primary.metadata,
            "float_trace": floating.metadata,
            "sinex": sinex.metadata,
        },
        "integer_pseudoobservation_events": [
            {
                "epoch": epoch_text(epoch),
                "attempt_row_count": row_count,
                "submitted_unverified_row_count": (
                    primary.submitted_unverified_rows_by_epoch.get(epoch, 0)
                ),
                "submission_status": (
                    "filter_call_returned_submitted_unverified"
                    if epoch in primary.submitted_unverified_rows_by_epoch
                    else "attempt_only_not_confirmed_submitted"
                ),
                "receivers": sorted(
                    primary.pseudoobs_attempt_receiver_rows_by_epoch.get(epoch, {})
                ),
                "receiver_attempt_row_counts": dict(
                    sorted(
                        primary.pseudoobs_attempt_receiver_rows_by_epoch.get(
                            epoch, {}
                        ).items()
                    )
                ),
            }
            for epoch, row_count in sorted(primary.pseudoobs_attempt_rows_by_epoch.items())
        ],
        "integer_pseudoobservation_detection": {
            "status": pseudoobs_detection_status,
            "rule": (
                "Applying rows are attempts only; submitted_unverified requires "
                "PSEUDOOBS_SUBMISSION status=FILTER_CALL_RETURNED_SUBMITTED_UNVERIFIED"
            ),
            "visibility_limit": (
                "telemetry depends on trace verbosity and code version; submitted_unverified "
                "does not assert that an integer is correct or scientifically accepted"
            ),
        },
        "thresholds_m": {
            "ar_vs_float_3d": thresholds.ar_float_m,
            "trace_vs_sinex_3d": thresholds.sinex_m,
            "primary_epoch_jump_3d": thresholds.jump_m,
        },
        "screening": {
            "status": screening_status,
            "evaluated_value_counts": evaluated_counts,
            "configured_metrics_without_evaluable_values": configured_without_values,
            "exceedance_counts": dict(sorted(exceedance_counts.items())),
            "scientific_interpretation": (
                "Threshold exceedances are screening flags, not diagnoses. Absence of an "
                "exceedance is not PPP-AR acceptance. Independent restarts/reinitialization, "
                "FLOAT controls, and integer-datum validation remain necessary. Temporal "
                "association with a submitted-unverified event does not establish causation."
            ),
        },
        "data_quality": {
            "status": (
                "issues_observed_inconclusive"
                if data_quality_issues
                else "no_reported_data_quality_issues"
            ),
            "issues": data_quality_issues,
            "primary_float_exact_epoch_coverage_fraction": comparison_coverage,
            "common_trace_epoch_receiver_count": common_trace_key_count,
            "trace_epoch_receiver_union_count": trace_key_union_count,
            "interpretation": (
                "Data-quality status is not a scientific pass. Any listed issue limits "
                "or prevents interpretation of threshold screening."
            ),
        },
        "station_epoch_comparisons": comparisons,
        "primary_coordinate_jumps": jumps,
        "statistics": {
            "overall": {metric: statistics(metric_values[metric]) for metric in metrics},
            "by_receiver": {
                receiver: {
                    metric: statistics(receiver_metric_values[receiver][metric])
                    for metric in metrics
                }
                for receiver in receivers
            },
        },
        "missing_data": {
            "missing_from_primary_trace": [
                _key_json(key) for key in sorted(float_keys - primary_keys)
            ],
            "missing_from_float_trace": [
                _key_json(key) for key in sorted(primary_keys - float_keys)
            ],
            "receivers_without_complete_sinex_reference": sorted(
                receiver
                for receiver in trace_receivers
                if receiver not in sinex.references_by_receiver
            ),
            "sinex_reference_unavailable_for_station_epoch": [
                {
                    "receiver": row["receiver"],
                    "epoch": row["epoch"],
                    "selection_status": row["sinex_reference_selection"]["status"],
                }
                for row in comparisons
                if row["sinex_reference"] is None
            ],
        },
    }


def nonnegative_threshold(text: str) -> float:
    try:
        value = finite_float(text)
    except ValueError as error:
        raise argparse.ArgumentTypeError(str(error)) from error
    if value < 0:
        raise argparse.ArgumentTypeError("threshold must be non-negative")
    return value


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("primary_trace", type=Path, help="candidate PPP-AR trace")
    parser.add_argument("float_trace", type=Path, help="independent FLOAT-control trace")
    parser.add_argument("igs_crd_sinex", type=Path, help="IGS CRD SINEX reference")
    parser.add_argument("--output", type=Path, help="write JSON here instead of stdout")
    parser.add_argument("--ar-float-threshold-m", type=nonnegative_threshold)
    parser.add_argument("--sinex-threshold-m", type=nonnegative_threshold)
    parser.add_argument("--jump-threshold-m", type=nonnegative_threshold)
    parser.add_argument(
        "--primary-state-block",
        help="explicit state block, for example AR or AR_RTS; default requires a unique AR block",
    )
    parser.add_argument(
        "--float-state-block",
        help=(
            "explicit state block, for example PPP or PPP_RTS; default requires "
            "a unique PPP block"
        ),
    )
    parser.add_argument(
        "--max-sinex-reference-age-days",
        type=nonnegative_threshold,
        help="optional maximum absolute age from SINEX reference epoch",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    payload = audit_coordinate_integrity(
        args.primary_trace,
        args.float_trace,
        args.igs_crd_sinex,
        Thresholds(
            ar_float_m=args.ar_float_threshold_m,
            sinex_m=args.sinex_threshold_m,
            jump_m=args.jump_threshold_m,
        ),
        primary_state_block=args.primary_state_block,
        float_state_block=args.float_state_block,
        max_sinex_reference_age_days=args.max_sinex_reference_age_days,
    )
    serialized = json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n"
    if args.output is None:
        print(serialized, end="")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        partial = args.output.with_name(args.output.name + ".part")
        partial.write_text(serialized, encoding="utf-8")
        os.replace(partial, args.output)
    # This audit deliberately reports evidence and screening flags.  It does
    # not turn a lack of threshold exceedances into a scientific pass verdict.
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
