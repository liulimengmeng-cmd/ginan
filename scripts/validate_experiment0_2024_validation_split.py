#!/usr/bin/env python3
"""Audit the frozen Experiment 0 model/validation station separation.

The model audit remains authoritative for the ten-station PEA runs.  This
second guard additionally requires the validation-only YAMLs to contain
exactly STR2, BALL and PARK, with exact day-specific paths from the frozen
input manifest and no model-station observation.
"""

from __future__ import annotations

import argparse
import json
import os
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Sequence

from validate_experiment0_2024_split import (
    AuditInputError,
    EXPECTED_HELDOUT,
    EXPECTED_TRAIN,
    build_audit,
    explicit_rnx_inputs,
    station_from_rinex_path,
)


def _manifest_heldout_paths(
    manifest: dict[str, object], label: str
) -> set[str]:
    heldout = set(EXPECTED_HELDOUT)
    paths: set[str] = set()
    files = manifest.get("files")
    if not isinstance(files, list):
        raise AuditInputError("manifest.files must be a list")
    for record in files:
        if not isinstance(record, dict) or record.get("experiment") != label:
            continue
        if "observation" not in str(record.get("role", "")).lower():
            continue
        relative_path = record.get("relative_path")
        if not isinstance(relative_path, str):
            raise AuditInputError("manifest observation relative_path must be a string")
        if station_from_rinex_path(relative_path) in heldout:
            paths.add(relative_path.replace("\\", "/"))
    return paths


def build_validation_audit(
    manifest_path: Path,
    quiet_model_yaml: Path,
    storm_model_yaml: Path,
    quiet_heldout_yaml: Path,
    storm_heldout_yaml: Path,
) -> dict[str, object]:
    model = build_audit(manifest_path, quiet_model_yaml, storm_model_yaml)
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if not isinstance(manifest, dict):
        raise AuditInputError("input manifest root must be a JSON object")

    errors = [] if model["valid"] else ["model_split_audit_failed"]
    train = set(EXPECTED_TRAIN)
    heldout = set(EXPECTED_HELDOUT)
    per_day: dict[str, dict[str, object]] = {}
    for label, yaml_path in (
        ("quiet", quiet_heldout_yaml),
        ("storm", storm_heldout_yaml),
    ):
        values = explicit_rnx_inputs(yaml_path)
        normalized = [value.replace("\\", "/") for value in values]
        stations = [station_from_rinex_path(value) for value in values]
        station_set = set(stations)
        expected_paths = _manifest_heldout_paths(manifest, label)
        actual_paths = set(normalized)
        duplicates = sorted({station for station in stations if stations.count(station) > 1})
        model_inputs = sorted(station_set & train)
        missing_heldout = sorted(heldout - station_set)
        unexpected = sorted(station_set - heldout)
        missing_paths = sorted(expected_paths - actual_paths)
        unexpected_paths = sorted(actual_paths - expected_paths)

        if len(values) != 3:
            errors.append(f"{label} heldout YAML must contain exactly 3 rnx_inputs")
        if duplicates:
            errors.append(f"{label} heldout YAML has duplicate stations: {duplicates}")
        if model_inputs:
            errors.append(f"{label} heldout YAML contains model stations: {model_inputs}")
        if missing_heldout or unexpected:
            errors.append(
                f"{label} heldout station set mismatch: missing={missing_heldout}, "
                f"unexpected={unexpected}"
            )
        if missing_paths or unexpected_paths:
            errors.append(
                f"{label} heldout paths differ from manifest: missing={missing_paths}, "
                f"unexpected={unexpected_paths}"
            )

        per_day[label] = {
            "yaml": str(yaml_path.resolve()),
            "rnx_input_count": len(values),
            "rnx_inputs": normalized,
            "stations": sorted(station_set),
            "model_input_count": len(model_inputs),
            "model_inputs": model_inputs,
            "missing_heldout": missing_heldout,
            "unexpected_stations": unexpected,
            "paths_match_manifest": not missing_paths and not unexpected_paths,
        }

    roles_identical = (
        set(per_day["quiet"]["stations"])
        == set(per_day["storm"]["stations"])
        == heldout
    )
    if not roles_identical:
        errors.append("quiet/storm heldout roles are not identical")

    return {
        "schema": "GINAN_EXPERIMENT0_2024_VALIDATION_SPLIT_AUDIT_V1",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "valid": not errors,
        "manifest": str(manifest_path.resolve()),
        "model_split_valid": bool(model["valid"]),
        "model_stations": sorted(train),
        "heldout_stations": sorted(heldout),
        "intersection": sorted(train & heldout),
        "roles_identical": roles_identical,
        "per_day": per_day,
        "errors": errors,
    }


def _atomic_write_json(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="\n") as stream:
            json.dump(payload, stream, indent=2)
            stream.write("\n")
        os.replace(temporary_name, path)
    except BaseException:
        try:
            os.unlink(temporary_name)
        except FileNotFoundError:
            pass
        raise


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("quiet_model_yaml", type=Path)
    parser.add_argument("storm_model_yaml", type=Path)
    parser.add_argument("quiet_heldout_yaml", type=Path)
    parser.add_argument("storm_heldout_yaml", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args(argv)

    try:
        audit = build_validation_audit(
            args.manifest,
            args.quiet_model_yaml,
            args.storm_model_yaml,
            args.quiet_heldout_yaml,
            args.storm_heldout_yaml,
        )
    except (AuditInputError, OSError, UnicodeError, json.JSONDecodeError) as error:
        print(f"validation split audit input error: {error}")
        return 2
    _atomic_write_json(args.output, audit)
    print(json.dumps(audit, indent=2))
    return 0 if audit["valid"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
