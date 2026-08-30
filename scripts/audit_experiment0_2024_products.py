#!/usr/bin/env python3
"""Audit the frozen central-day precise products for Experiment 0.

The script is deliberately offline and read-only with respect to the input
tree.  It scans the WUM final SP3, CLK and Bias-SINEX files already fetched by
``fetch_experiment0_2024_data.py`` and writes one JSON audit using an atomic
same-directory replacement.

A GPS satellite belongs to the product-availability intersection (``eligible``)
when all of the following hold for the central UTC day:

* at least one valid SP3 position is present;
* at least one valid satellite clock is present;
* C1C, C2W, L1C and L2W OSB intervals each cover the full day.

The audit separately computes the stricter intersection requiring every
5-minute SP3 and every 30-second CLK epoch.  This distinction is intentional:
an analysis-centre clock file can contain a short product-wide outage without
changing its satellite roster.  Such outages remain validation failures and
are reported per satellite rather than silently dropping the affected
satellites from the product-availability intersection.

The frozen expected eligible intersections are G02-G09,G11-G32 for the quiet
day and G02-G32 for the storm day.  An eligible-set mismatch, metadata mismatch
or full-day coverage defect is recorded in the JSON and results in exit status
2.
"""

from __future__ import annotations

import argparse
import gzip
import json
import math
import os
import re
import sys
import tempfile
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from datetime import date, datetime, timedelta, timezone
from pathlib import Path
from typing import Iterable, TextIO


SCHEMA = "GINAN_EXPERIMENT0_2024_PRODUCT_AUDIT_V1"
GPS_SATELLITES = tuple(f"G{prn:02d}" for prn in range(1, 33))
REQUIRED_OBSERVABLES = ("C1C", "C2W", "L1C", "L2W")
EXPECTED_APC_MODEL = "IGS20_2303.ATX"
EXPECTED_GPS_CLOCK_REFERENCE = ("C1W", "C2W")
SP3_INTERVAL_SECONDS = 300
CLK_INTERVAL_SECONDS = 30

EXPERIMENTS = {
    "quiet": {
        "day": date(2024, 5, 8),
        "doy": 129,
        "eligible": tuple(f"G{prn:02d}" for prn in range(2, 33) if prn != 10),
    },
    "storm": {
        "day": date(2024, 5, 11),
        "doy": 132,
        "eligible": tuple(f"G{prn:02d}" for prn in range(2, 33)),
    },
}

GPS_PRN_RE = re.compile(r"^G\d{2}$")


class ProductAuditError(RuntimeError):
    """Raised when an input product cannot be located or parsed."""


@dataclass
class EpochRecords:
    valid: dict[str, Counter[datetime]] = field(
        default_factory=lambda: defaultdict(Counter)
    )
    raw_count: Counter[str] = field(default_factory=Counter)
    invalid_count: Counter[str] = field(default_factory=Counter)
    off_grid_count: Counter[str] = field(default_factory=Counter)
    outside_day_count: Counter[str] = field(default_factory=Counter)


@dataclass(frozen=True)
class BiasInterval:
    start: datetime
    end: datetime


@dataclass
class BiasRecords:
    metadata: dict[str, object] = field(default_factory=dict)
    intervals: dict[str, dict[str, list[BiasInterval]]] = field(
        default_factory=lambda: defaultdict(lambda: defaultdict(list))
    )
    raw_count: dict[str, Counter[str]] = field(
        default_factory=lambda: defaultdict(Counter)
    )
    invalid_count: dict[str, Counter[str]] = field(
        default_factory=lambda: defaultdict(Counter)
    )


def open_text(path: Path) -> TextIO:
    """Open a plain or gzip-compressed ASCII GNSS product as text."""

    if path.suffix.lower() == ".gz":
        return gzip.open(path, "rt", encoding="ascii", errors="replace")
    return path.open("rt", encoding="ascii", errors="replace")


def utc_day_bounds(day: date) -> tuple[datetime, datetime]:
    start = datetime(day.year, day.month, day.day, tzinfo=timezone.utc)
    return start, start + timedelta(days=1)


def expected_epochs(day: date, interval_seconds: int) -> tuple[datetime, ...]:
    start, end = utc_day_bounds(day)
    count = int((end - start).total_seconds()) // interval_seconds
    return tuple(start + timedelta(seconds=index * interval_seconds) for index in range(count))


def datetime_from_calendar_fields(fields: Iterable[str]) -> datetime:
    values = list(fields)
    if len(values) != 6:
        raise ValueError(f"expected six calendar fields, received {values}")
    year, month, day, hour, minute = (int(value) for value in values[:5])
    second = float(values[5])
    base = datetime(year, month, day, hour, minute, tzinfo=timezone.utc)
    return base + timedelta(seconds=second)


def parse_bias_time_tag(tag: str) -> datetime | None:
    match = re.fullmatch(r"(\d{4}):(\d{3}):(\d{5})", tag)
    if not match:
        raise ValueError(f"invalid Bias-SINEX time tag: {tag}")
    year, doy, seconds = (int(value) for value in match.groups())
    if year == 0 or doy == 0:
        return None
    return datetime(year, 1, 1, tzinfo=timezone.utc) + timedelta(
        days=doy - 1, seconds=seconds
    )


def is_valid_sp3_position(tokens: list[str]) -> bool:
    if len(tokens) < 4:
        return False
    try:
        coordinates = tuple(float(value) for value in tokens[1:4])
    except ValueError:
        return False
    if not all(math.isfinite(value) and abs(value) < 999999 for value in coordinates):
        return False
    return any(abs(value) > 1e-9 for value in coordinates)


def is_valid_clock(tokens: list[str]) -> bool:
    if len(tokens) < 10:
        return False
    try:
        value_count = int(tokens[8])
        clock_value = float(tokens[9])
    except ValueError:
        return False
    return value_count >= 1 and math.isfinite(clock_value) and abs(clock_value) < 999999


def parse_sp3(path: Path, day: date) -> EpochRecords:
    records = EpochRecords()
    start, end = utc_day_bounds(day)
    expected = set(expected_epochs(day, SP3_INTERVAL_SECONDS))
    current_epoch: datetime | None = None

    with open_text(path) as stream:
        for line_number, line in enumerate(stream, start=1):
            if line.startswith("*"):
                tokens = line.split()
                try:
                    current_epoch = datetime_from_calendar_fields(tokens[1:7])
                except (ValueError, IndexError) as error:
                    raise ProductAuditError(
                        f"invalid SP3 epoch at {path}:{line_number}: {line.rstrip()}"
                    ) from error
                continue
            if not line.startswith("PG"):
                continue
            if current_epoch is None:
                raise ProductAuditError(
                    f"SP3 position appears before an epoch at {path}:{line_number}"
                )
            tokens = line.split()
            satellite = tokens[0][1:] if tokens else ""
            if satellite not in GPS_SATELLITES:
                continue
            if not (start <= current_epoch < end):
                records.outside_day_count[satellite] += 1
                continue
            records.raw_count[satellite] += 1
            if not is_valid_sp3_position(tokens):
                records.invalid_count[satellite] += 1
            elif current_epoch not in expected:
                records.off_grid_count[satellite] += 1
            else:
                records.valid[satellite][current_epoch] += 1
    return records


def parse_clk(path: Path, day: date) -> EpochRecords:
    records = EpochRecords()
    start, end = utc_day_bounds(day)
    expected = set(expected_epochs(day, CLK_INTERVAL_SECONDS))

    with open_text(path) as stream:
        for line_number, line in enumerate(stream, start=1):
            if not line.startswith("AS G"):
                continue
            tokens = line.split()
            if len(tokens) < 9:
                raise ProductAuditError(
                    f"truncated satellite CLK record at {path}:{line_number}: {line.rstrip()}"
                )
            satellite = tokens[1]
            if satellite not in GPS_SATELLITES:
                continue
            try:
                epoch = datetime_from_calendar_fields(tokens[2:8])
            except ValueError as error:
                raise ProductAuditError(
                    f"invalid CLK epoch at {path}:{line_number}: {line.rstrip()}"
                ) from error
            if not (start <= epoch < end):
                records.outside_day_count[satellite] += 1
                continue
            records.raw_count[satellite] += 1
            if not is_valid_clock(tokens):
                records.invalid_count[satellite] += 1
            elif epoch not in expected:
                records.off_grid_count[satellite] += 1
            else:
                records.valid[satellite][epoch] += 1
    return records


def parse_bia(path: Path) -> BiasRecords:
    records = BiasRecords()
    in_description = False

    with open_text(path) as stream:
        for line_number, line in enumerate(stream, start=1):
            stripped = line.strip()
            if stripped == "+BIAS/DESCRIPTION":
                in_description = True
                continue
            if stripped == "-BIAS/DESCRIPTION":
                in_description = False
                continue
            if in_description and stripped and not stripped.startswith("*"):
                tokens = stripped.split()
                key = tokens[0]
                values = tokens[1:]
                if key == "SATELLITE_CLOCK_REFERENCE_OBSERVABLES" and len(values) >= 3:
                    references = records.metadata.setdefault(
                        "satellite_clock_reference_observables", {}
                    )
                    assert isinstance(references, dict)
                    references[values[0]] = values[1:]
                elif key in {
                    "APC_MODEL",
                    "BIAS_MODE",
                    "DETERMINATION_METHOD",
                    "TIME_SYSTEM",
                }:
                    records.metadata[key.lower()] = " ".join(values)
                elif key in {"OBSERVATION_SAMPLING", "PARAMETER_SPACING"}:
                    try:
                        records.metadata[key.lower()] = int(values[0])
                    except (ValueError, IndexError) as error:
                        raise ProductAuditError(
                            f"invalid {key} metadata at {path}:{line_number}"
                        ) from error
                continue

            if not line.startswith(" OSB"):
                continue
            tokens = line.split()
            if len(tokens) < 9 or not GPS_PRN_RE.fullmatch(tokens[2]):
                continue
            satellite = tokens[2]
            observable = tokens[3]
            if observable not in REQUIRED_OBSERVABLES:
                continue
            records.raw_count[satellite][observable] += 1
            try:
                interval_start = parse_bias_time_tag(tokens[4])
                interval_end = parse_bias_time_tag(tokens[5])
                estimate = float(tokens[7])
            except (ValueError, IndexError) as error:
                records.invalid_count[satellite][observable] += 1
                continue
            if (
                interval_start is None
                or interval_end is None
                or interval_end <= interval_start
                or not math.isfinite(estimate)
            ):
                records.invalid_count[satellite][observable] += 1
                continue
            records.intervals[satellite][observable].append(
                BiasInterval(interval_start, interval_end)
            )
    return records


def group_missing_epochs(
    missing: Iterable[datetime], interval_seconds: int
) -> list[dict[str, object]]:
    epochs = sorted(missing)
    if not epochs:
        return []
    groups: list[list[datetime]] = [[epochs[0]]]
    step = timedelta(seconds=interval_seconds)
    for epoch in epochs[1:]:
        if epoch - groups[-1][-1] == step:
            groups[-1].append(epoch)
        else:
            groups.append([epoch])
    return [
        {
            "first_epoch": group[0].isoformat(),
            "last_epoch": group[-1].isoformat(),
            "count": len(group),
        }
        for group in groups
    ]


def summarize_epoch_records(
    records: EpochRecords, day: date, interval_seconds: int
) -> tuple[dict[str, dict[str, object]], set[str]]:
    expected_tuple = expected_epochs(day, interval_seconds)
    expected = set(expected_tuple)
    summaries: dict[str, dict[str, object]] = {}
    complete: set[str] = set()

    for satellite in GPS_SATELLITES:
        counter = records.valid[satellite]
        valid_epochs = set(counter)
        missing = expected - valid_epochs
        duplicates = sum(max(0, count - 1) for count in counter.values())
        valid_record_count = sum(counter.values())
        summary = {
            "raw_record_count": records.raw_count[satellite],
            "valid_record_count": valid_record_count,
            "unique_valid_epoch_count": len(valid_epochs),
            "invalid_record_count": records.invalid_count[satellite],
            "duplicate_epoch_count": duplicates,
            "off_grid_record_count": records.off_grid_count[satellite],
            "outside_day_record_count": records.outside_day_count[satellite],
            "missing_epoch_count": len(missing),
            "gaps": group_missing_epochs(missing, interval_seconds),
            "first_valid_epoch": min(valid_epochs).isoformat() if valid_epochs else None,
            "last_valid_epoch": max(valid_epochs).isoformat() if valid_epochs else None,
        }
        is_complete = (
            len(missing) == 0
            and duplicates == 0
            and records.invalid_count[satellite] == 0
            and records.off_grid_count[satellite] == 0
            and records.raw_count[satellite] == len(expected_tuple)
            and valid_record_count == len(expected_tuple)
        )
        summary["complete"] = is_complete
        if is_complete:
            complete.add(satellite)
        summaries[satellite] = summary
    return summaries, complete


def merge_intervals(intervals: Iterable[BiasInterval]) -> list[BiasInterval]:
    ordered = sorted(intervals, key=lambda interval: (interval.start, interval.end))
    merged: list[BiasInterval] = []
    for interval in ordered:
        if not merged or interval.start > merged[-1].end:
            merged.append(interval)
        elif interval.end > merged[-1].end:
            merged[-1] = BiasInterval(merged[-1].start, interval.end)
    return merged


def interval_gaps(
    intervals: Iterable[BiasInterval], start: datetime, end: datetime
) -> list[dict[str, object]]:
    clipped = [
        BiasInterval(max(interval.start, start), min(interval.end, end))
        for interval in intervals
        if interval.end > start and interval.start < end
    ]
    merged = merge_intervals(clipped)
    cursor = start
    gaps: list[dict[str, object]] = []
    for interval in merged:
        if interval.start > cursor:
            gaps.append(
                {
                    "start": cursor.isoformat(),
                    "end_exclusive": interval.start.isoformat(),
                    "duration_seconds": int((interval.start - cursor).total_seconds()),
                }
            )
        cursor = max(cursor, interval.end)
    if cursor < end:
        gaps.append(
            {
                "start": cursor.isoformat(),
                "end_exclusive": end.isoformat(),
                "duration_seconds": int((end - cursor).total_seconds()),
            }
        )
    return gaps


def summarize_bias_records(
    records: BiasRecords, day: date
) -> tuple[dict[str, dict[str, object]], set[str]]:
    start, end = utc_day_bounds(day)
    summaries: dict[str, dict[str, object]] = {}
    complete: set[str] = set()

    for satellite in GPS_SATELLITES:
        observable_summaries: dict[str, object] = {}
        satellite_complete = True
        for observable in REQUIRED_OBSERVABLES:
            intervals = records.intervals[satellite][observable]
            gaps = interval_gaps(intervals, start, end)
            invalid_count = records.invalid_count[satellite][observable]
            covers_full_day = not gaps and invalid_count == 0 and bool(intervals)
            observable_summaries[observable] = {
                "raw_record_count": records.raw_count[satellite][observable],
                "valid_interval_count": len(intervals),
                "invalid_record_count": invalid_count,
                "covers_full_day": covers_full_day,
                "coverage_gaps": gaps,
            }
            satellite_complete = satellite_complete and covers_full_day
        if satellite_complete:
            complete.add(satellite)
        summaries[satellite] = {
            "observables": observable_summaries,
            "complete": satellite_complete,
        }
    return summaries, complete


def find_product(root: Path, label: str, filename: str) -> Path:
    directory = root / label / "products"
    candidates = [directory / filename, directory / f"{filename}.gz"]
    found = [path for path in candidates if path.is_file()]
    if len(found) != 1:
        rendered = ", ".join(str(path) for path in candidates)
        raise ProductAuditError(
            f"expected exactly one {label} product for {filename}; checked {rendered}"
        )
    return found[0]


def relative_or_absolute(path: Path, root: Path) -> str:
    try:
        return path.relative_to(root).as_posix()
    except ValueError:
        return str(path)


def audit_experiment(root: Path, label: str) -> dict[str, object]:
    specification = EXPERIMENTS[label]
    day = specification["day"]
    doy = specification["doy"]
    assert isinstance(day, date)
    assert isinstance(doy, int)

    prefix = f"WUM0MGXFIN_2024{doy:03d}0000_"
    sp3_path = find_product(root, label, prefix + "01D_05M_ORB.SP3")
    clk_path = find_product(root, label, prefix + "01D_30S_CLK.CLK")
    bia_path = find_product(root, label, prefix + "01D_01D_OSB.BIA")

    sp3_records = parse_sp3(sp3_path, day)
    clk_records = parse_clk(clk_path, day)
    bia_records = parse_bia(bia_path)
    sp3_summary, sp3_complete = summarize_epoch_records(
        sp3_records, day, SP3_INTERVAL_SECONDS
    )
    clk_summary, clk_complete = summarize_epoch_records(
        clk_records, day, CLK_INTERVAL_SECONDS
    )
    bia_summary, bia_complete = summarize_bias_records(bia_records, day)

    sp3_available = {
        satellite
        for satellite, summary in sp3_summary.items()
        if int(summary["unique_valid_epoch_count"]) > 0
    }
    clk_available = {
        satellite
        for satellite, summary in clk_summary.items()
        if int(summary["unique_valid_epoch_count"]) > 0
    }
    eligible = sp3_available & clk_available & bia_complete
    full_day_complete = sp3_complete & clk_complete & bia_complete
    expected = set(specification["eligible"])
    missing = sorted(expected - eligible)
    unexpected = sorted(eligible - expected)
    eligible_with_coverage_defects = sorted(eligible - full_day_complete)

    clock_references = bia_records.metadata.get(
        "satellite_clock_reference_observables", {}
    )
    gps_clock_reference: list[str] = []
    if isinstance(clock_references, dict):
        value = clock_references.get("G", [])
        if isinstance(value, list):
            gps_clock_reference = value
    apc_model = str(bia_records.metadata.get("apc_model", ""))
    bias_mode = str(bia_records.metadata.get("bias_mode", ""))
    metadata_validation = {
        "apc_model_matches_expected": apc_model.upper() == EXPECTED_APC_MODEL,
        "gps_clock_reference_matches_expected": tuple(gps_clock_reference)
        == EXPECTED_GPS_CLOCK_REFERENCE,
        "bias_mode_is_absolute": bias_mode.upper() == "ABSOLUTE",
    }
    metadata_validation["passed"] = all(metadata_validation.values())
    intersection_passed = not missing and not unexpected
    full_day_completeness_passed = not eligible_with_coverage_defects

    return {
        "date": day.isoformat(),
        "files": {
            "sp3": relative_or_absolute(sp3_path, root),
            "clk": relative_or_absolute(clk_path, root),
            "bia": relative_or_absolute(bia_path, root),
        },
        "sp3": {
            "nominal_interval_seconds": SP3_INTERVAL_SECONDS,
            "expected_epoch_count": len(expected_epochs(day, SP3_INTERVAL_SECONDS)),
            "available_satellites": sorted(sp3_available),
            "complete_satellites": sorted(sp3_complete),
            "satellites": sp3_summary,
        },
        "clk": {
            "nominal_interval_seconds": CLK_INTERVAL_SECONDS,
            "expected_epoch_count": len(expected_epochs(day, CLK_INTERVAL_SECONDS)),
            "available_satellites": sorted(clk_available),
            "complete_satellites": sorted(clk_complete),
            "satellites": clk_summary,
        },
        "bia": {
            "required_observables": list(REQUIRED_OBSERVABLES),
            "metadata": bia_records.metadata,
            "complete_satellites": sorted(bia_complete),
            "satellites": bia_summary,
        },
        "eligible_gps_satellites": sorted(eligible),
        "full_day_complete_gps_satellites": sorted(full_day_complete),
        "expected_eligible_gps_satellites": sorted(expected),
        "validation": {
            "eligible_intersection_passed": intersection_passed,
            "missing_expected_satellites": missing,
            "unexpected_satellites": unexpected,
            "full_day_completeness_passed": full_day_completeness_passed,
            "eligible_satellites_with_coverage_defects": (
                eligible_with_coverage_defects
            ),
            "metadata": metadata_validation,
            "passed": (
                intersection_passed
                and full_day_completeness_passed
                and bool(metadata_validation["passed"])
            ),
        },
    }


def audit_data_root(root: Path) -> dict[str, object]:
    resolved_root = root.resolve()
    if not resolved_root.is_dir():
        raise ProductAuditError(f"input data root is not a directory: {resolved_root}")
    experiments = {
        label: audit_experiment(resolved_root, label) for label in EXPERIMENTS
    }
    passed = all(
        bool(experiment["validation"]["passed"])
        for experiment in experiments.values()
    )
    return {
        "schema": SCHEMA,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "input_root": str(resolved_root),
        "eligibility_rule": (
            "central-day SP3 and CLK product presence intersected with full-day "
            "C1C/C2W/L1C/L2W Bias-SINEX coverage"
        ),
        "full_day_completeness_rule": (
            "eligible satellites must additionally contain every 5-minute SP3 "
            "epoch and every 30-second CLK epoch without invalid, duplicate or "
            "off-grid records"
        ),
        "experiments": experiments,
        "validation_passed": passed,
    }


def write_json_atomic(payload: dict[str, object], output_path: Path) -> None:
    output_path = output_path.resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    temporary_name: str | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            dir=output_path.parent,
            prefix=f".{output_path.name}.",
            suffix=".tmp",
            delete=False,
            newline="\n",
        ) as stream:
            temporary_name = stream.name
            json.dump(payload, stream, indent=2, sort_keys=True)
            stream.write("\n")
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary_name, output_path)
        temporary_name = None
    finally:
        if temporary_name is not None:
            Path(temporary_name).unlink(missing_ok=True)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("data_root", type=Path, help="frozen Experiment 0 input root")
    parser.add_argument(
        "--output",
        type=Path,
        help=(
            "audit JSON path (default: DATA_ROOT/experiment0_2024_product_audit.json)"
        ),
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    root = args.data_root.resolve()
    output = (
        args.output.resolve()
        if args.output is not None
        else root / "experiment0_2024_product_audit.json"
    )
    try:
        payload = audit_data_root(root)
        write_json_atomic(payload, output)
    except (OSError, ProductAuditError) as error:
        print(f"product audit failed: {error}", file=sys.stderr)
        return 1

    for label, experiment in payload["experiments"].items():
        validation = experiment["validation"]
        eligible = ",".join(experiment["eligible_gps_satellites"])
        print(
            f"{label}: passed={validation['passed']} eligible={eligible}",
            flush=True,
        )
        if validation["missing_expected_satellites"]:
            print(
                f"{label}: missing={','.join(validation['missing_expected_satellites'])}",
                flush=True,
            )
        if validation["unexpected_satellites"]:
            print(
                f"{label}: unexpected={','.join(validation['unexpected_satellites'])}",
                flush=True,
            )
        if validation["eligible_satellites_with_coverage_defects"]:
            print(
                f"{label}: coverage-defects="
                f"{','.join(validation['eligible_satellites_with_coverage_defects'])}",
                flush=True,
            )
    print(f"audit JSON: {output}", flush=True)
    return 0 if payload["validation_passed"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
