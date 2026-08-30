#!/usr/bin/env python3
"""Audit the frozen Experiment 0 2024 train/held-out station split.

Usage::

    python scripts/validate_experiment0_2024_split.py \
        input_manifest.json quiet_model.yaml storm_model.yaml split_audit.json

Only the explicit ``inputs.gnss_observations.rnx_inputs`` list in each
day-level YAML is inspected.  Included YAML files are deliberately not
resolved: the quiet/storm overlays are the authority for observation inputs.
The audit is written atomically and the process exits with status 0 only when
all split invariants hold.
"""

from __future__ import annotations

import argparse
import ast
import json
import os
import re
import tempfile
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Sequence


EXPECTED_DATES = {
    "quiet": "2024-05-08",
    "storm": "2024-05-11",
}

EXPECTED_TRAIN = (
    "ARMC",
    "BATH",
    "HOB2",
    "MCHL",
    "MOBS",
    "STR1",
    "SYDN",
    "TID1",
    "WGGA",
    "YARR",
)

EXPECTED_HELDOUT = (
    "STR2",
    "BALL",
    "PARK",
)

RINEX_OBSERVATION_NAME = re.compile(
    r"^(?P<station>[A-Z0-9]{4})[A-Z0-9]{5}_[RSU]_"
    r"\d{11}_01D_30S_MO\.(?:RNX|CRX)(?:\.GZ)?$",
    re.IGNORECASE,
)

YAML_MAPPING_KEY = re.compile(r"^(?P<indent> *)(?P<key>[A-Za-z_][A-Za-z0-9_-]*):(?P<value>.*)$")
YAML_LIST_ITEM = re.compile(r"^(?P<indent> *)-\s*(?P<value>.+?)\s*$")


class AuditInputError(ValueError):
    """Raised when an audit input cannot be interpreted unambiguously."""


def _strip_yaml_comment(value: str) -> str:
    quote: str | None = None
    escaped = False
    for index, character in enumerate(value):
        if escaped:
            escaped = False
            continue
        if character == "\\" and quote == '"':
            escaped = True
            continue
        if character in ("'", '"'):
            if quote is None:
                quote = character
            elif quote == character:
                quote = None
            continue
        if character == "#" and quote is None:
            return value[:index].rstrip()
    return value.rstrip()


def _parse_yaml_scalar(value: str, path: Path, line_number: int) -> str:
    value = _strip_yaml_comment(value).strip()
    if not value:
        raise AuditInputError(f"empty rnx_inputs item: {path}:{line_number}")
    if value[0] in ("'", '"'):
        try:
            parsed = ast.literal_eval(value)
        except (SyntaxError, ValueError) as error:
            raise AuditInputError(
                f"invalid quoted rnx_inputs item: {path}:{line_number}"
            ) from error
        if not isinstance(parsed, str):
            raise AuditInputError(f"rnx_inputs item is not a string: {path}:{line_number}")
        return parsed
    return value


def explicit_rnx_inputs(path: Path) -> list[str]:
    """Extract the one explicit day-level ``rnx_inputs`` scalar list."""

    lines = path.read_text(encoding="utf-8").splitlines()
    mapping_stack: list[tuple[int, str]] = []
    target_indent: int | None = None
    target_count = 0
    values: list[str] = []

    for line_number, raw_line in enumerate(lines, start=1):
        if "\t" in raw_line[: len(raw_line) - len(raw_line.lstrip())]:
            raise AuditInputError(f"tab indentation is not supported: {path}:{line_number}")
        stripped = raw_line.strip()
        if not stripped or stripped.startswith("#"):
            continue

        indent = len(raw_line) - len(raw_line.lstrip(" "))
        if target_indent is not None:
            if indent > target_indent:
                item = YAML_LIST_ITEM.match(raw_line)
                if item:
                    values.append(_parse_yaml_scalar(item.group("value"), path, line_number))
                    continue
                raise AuditInputError(
                    f"rnx_inputs must contain only scalar list items: {path}:{line_number}"
                )
            target_indent = None

        mapping = YAML_MAPPING_KEY.match(raw_line)
        if not mapping:
            continue
        indent = len(mapping.group("indent"))
        while mapping_stack and mapping_stack[-1][0] >= indent:
            mapping_stack.pop()
        key = mapping.group("key")
        path_keys = [entry[1] for entry in mapping_stack] + [key]
        inline_value = _strip_yaml_comment(mapping.group("value")).strip()
        if path_keys == ["inputs", "gnss_observations", "rnx_inputs"]:
            target_count += 1
            if inline_value:
                raise AuditInputError(
                    f"rnx_inputs must use a block list, not an inline value: {path}:{line_number}"
                )
            target_indent = indent
        mapping_stack.append((indent, key))

    if target_count != 1:
        raise AuditInputError(
            f"expected one explicit inputs.gnss_observations.rnx_inputs block in {path}; "
            f"found {target_count}"
        )
    if not values:
        raise AuditInputError(f"rnx_inputs is empty: {path}")
    return values


def station_from_rinex_path(value: str) -> str:
    name = value.replace("\\", "/").rsplit("/", 1)[-1]
    match = RINEX_OBSERVATION_NAME.fullmatch(name)
    if not match:
        raise AuditInputError(f"not a frozen 01D/30S RINEX3 observation path: {value}")
    return match.group("station").upper()


def _station_list(manifest: dict[str, object], key: str, errors: list[str]) -> list[str]:
    value = manifest.get(key)
    if not isinstance(value, list) or not all(isinstance(item, str) for item in value):
        errors.append(f"manifest.{key} must be a list of station strings")
        return []
    stations = [item.upper() for item in value]
    duplicates = sorted(station for station, count in Counter(stations).items() if count > 1)
    if duplicates:
        errors.append(f"manifest.{key} contains duplicate stations: {duplicates}")
    return stations


def _manifest_observations(
    manifest: dict[str, object],
    label: str,
    errors: list[str],
) -> list[str]:
    files = manifest.get("files")
    if not isinstance(files, list):
        errors.append("manifest.files must be a list")
        return []
    stations: list[str] = []
    for index, record in enumerate(files):
        if not isinstance(record, dict):
            errors.append(f"manifest.files[{index}] must be an object")
            continue
        if record.get("experiment") != label:
            continue
        if "observation" not in str(record.get("role", "")).lower():
            continue
        relative_path = record.get("relative_path")
        if not isinstance(relative_path, str):
            errors.append(f"manifest.files[{index}].relative_path must be a string")
            continue
        try:
            stations.append(station_from_rinex_path(relative_path))
        except AuditInputError as error:
            errors.append(f"manifest.files[{index}]: {error}")
    return stations


def build_audit(
    manifest_path: Path,
    quiet_yaml: Path,
    storm_yaml: Path,
) -> dict[str, object]:
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if not isinstance(manifest, dict):
        raise AuditInputError("input manifest root must be a JSON object")

    errors: list[str] = []
    if manifest.get("schema") != "GINAN_EXPERIMENT0_2024_INPUT_MANIFEST_V1":
        errors.append(
            "manifest.schema must be GINAN_EXPERIMENT0_2024_INPUT_MANIFEST_V1; "
            f"found {manifest.get('schema')}"
        )
    train = _station_list(manifest, "model_stations", errors)
    heldout = _station_list(manifest, "heldout_stations", errors)
    train_set = set(train)
    heldout_set = set(heldout)
    intersection = sorted(train_set & heldout_set)

    if len(train) != 10 or len(train_set) != 10:
        errors.append(f"manifest model split must contain 10 unique stations; found {len(train_set)}")
    if len(heldout) != 3 or len(heldout_set) != 3:
        errors.append(f"manifest held-out split must contain 3 unique stations; found {len(heldout_set)}")
    if intersection:
        errors.append(f"manifest model and held-out stations overlap: {intersection}")
    if train_set != set(EXPECTED_TRAIN):
        errors.append(
            "manifest model stations differ from the frozen Experiment 0 set: "
            f"missing={sorted(set(EXPECTED_TRAIN) - train_set)}, "
            f"unexpected={sorted(train_set - set(EXPECTED_TRAIN))}"
        )
    if heldout_set != set(EXPECTED_HELDOUT):
        errors.append(
            "manifest held-out stations differ from the frozen Experiment 0 set: "
            f"missing={sorted(set(EXPECTED_HELDOUT) - heldout_set)}, "
            f"unexpected={sorted(heldout_set - set(EXPECTED_HELDOUT))}"
        )

    all_manifest_stations = _station_list(manifest, "stations", errors)
    if set(all_manifest_stations) != train_set | heldout_set:
        errors.append("manifest.stations does not equal model_stations union heldout_stations")

    dates = manifest.get("dates")
    if dates != EXPECTED_DATES:
        errors.append(f"manifest.dates must equal {EXPECTED_DATES}; found {dates}")

    yaml_paths = {
        "quiet": quiet_yaml,
        "storm": storm_yaml,
    }
    per_day: dict[str, dict[str, object]] = {}
    for label, yaml_path in yaml_paths.items():
        try:
            rnx_inputs = explicit_rnx_inputs(yaml_path)
        except (AuditInputError, OSError, UnicodeError) as error:
            errors.append(f"{label} YAML: {error}")
            rnx_inputs = []

        input_stations: list[str] = []
        invalid_inputs: list[str] = []
        for value in rnx_inputs:
            try:
                input_stations.append(station_from_rinex_path(value))
            except AuditInputError:
                invalid_inputs.append(value)

        counts = Counter(input_stations)
        duplicate_stations = sorted(station for station, count in counts.items() if count > 1)
        station_set = set(input_stations)
        heldout_inputs = [
            value
            for value in rnx_inputs
            if station_from_rinex_path(value) in heldout_set
        ] if not invalid_inputs else [
            value
            for value in rnx_inputs
            if value not in invalid_inputs and station_from_rinex_path(value) in heldout_set
        ]
        missing_train = sorted(train_set - station_set)
        unexpected_stations = sorted(station_set - train_set)

        manifest_stations = _manifest_observations(manifest, label, errors)
        manifest_counts = Counter(manifest_stations)
        manifest_duplicates = sorted(
            station for station, count in manifest_counts.items() if count > 1
        )
        manifest_missing = sorted((train_set | heldout_set) - set(manifest_stations))
        manifest_unexpected = sorted(set(manifest_stations) - (train_set | heldout_set))

        if len(rnx_inputs) != 10:
            errors.append(f"{label} YAML must contain exactly 10 rnx_inputs; found {len(rnx_inputs)}")
        if invalid_inputs:
            errors.append(f"{label} YAML contains invalid RINEX paths: {invalid_inputs}")
        if duplicate_stations:
            errors.append(f"{label} YAML contains duplicate stations: {duplicate_stations}")
        if heldout_inputs:
            errors.append(f"{label} YAML contains held-out observations: {heldout_inputs}")
        if missing_train or unexpected_stations:
            errors.append(
                f"{label} YAML station set differs from manifest model stations: "
                f"missing={missing_train}, unexpected={unexpected_stations}"
            )
        if manifest_duplicates or manifest_missing or manifest_unexpected:
            errors.append(
                f"manifest {label} observations must contain each of the 13 frozen stations once: "
                f"duplicates={manifest_duplicates}, missing={manifest_missing}, "
                f"unexpected={manifest_unexpected}"
            )

        per_day[label] = {
            "yaml": str(yaml_path.resolve()),
            "rnx_input_count": len(rnx_inputs),
            "rnx_inputs": rnx_inputs,
            "train": sorted(station_set),
            "heldout_inputs": heldout_inputs,
            "heldout_input_count": len(heldout_inputs),
            "duplicate_stations": duplicate_stations,
            "invalid_inputs": invalid_inputs,
            "missing_train": missing_train,
            "unexpected_stations": unexpected_stations,
            "manifest_observation_stations": sorted(set(manifest_stations)),
        }

    quiet_train = set(per_day["quiet"]["train"])
    storm_train = set(per_day["storm"]["train"])
    roles_identical = quiet_train == storm_train == train_set == set(EXPECTED_TRAIN)
    if not roles_identical:
        errors.append(
            "quiet and storm model roles are not identical to the frozen train set: "
            f"quiet={sorted(quiet_train)}, storm={sorted(storm_train)}"
        )

    heldout_input_count = sum(
        int(day["heldout_input_count"])
        for day in per_day.values()
    )

    return {
        "schema": "GINAN_EXPERIMENT0_2024_SPLIT_AUDIT_V1",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "valid": not errors,
        "manifest": str(manifest_path.resolve()),
        "dates": EXPECTED_DATES,
        "train": sorted(train_set),
        "heldout": sorted(heldout_set),
        "intersection": intersection,
        "heldout_input_count": heldout_input_count,
        "roles_identical": roles_identical,
        "per_day": per_day,
        "errors": errors,
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
            json.dump(payload, stream, indent=2, sort_keys=True)
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
    parser.add_argument("manifest", type=Path)
    parser.add_argument("quiet_yaml", type=Path)
    parser.add_argument("storm_yaml", type=Path)
    parser.add_argument("output_json", type=Path)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        audit = build_audit(args.manifest, args.quiet_yaml, args.storm_yaml)
    except (AuditInputError, OSError, UnicodeError, json.JSONDecodeError) as error:
        audit = {
            "schema": "GINAN_EXPERIMENT0_2024_SPLIT_AUDIT_V1",
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "valid": False,
            "manifest": str(args.manifest.resolve()),
            "dates": EXPECTED_DATES,
            "train": [],
            "heldout": [],
            "intersection": [],
            "heldout_input_count": 0,
            "roles_identical": False,
            "per_day": {},
            "errors": [str(error)],
        }
    write_json_atomic(args.output_json, audit)
    print(json.dumps(audit, indent=2, sort_keys=True))
    return 0 if audit["valid"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
