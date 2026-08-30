#!/usr/bin/env python3
"""Blind freeze/evaluate driver for the Experiment 0 2024 float pivot.

The driver is intentionally separate from the scientific evaluator.  It adds
an externally hashable freeze, immutable input receipts, exact held-out run
plans, and no-replace publication around
``experiment0_2024_spatial_validation.build_freeze`` / ``build_report``.

Freeze plan JSON (paths may be absolute or repository-relative)::

  {
    "schema": "GINAN_EXPERIMENT0_2024_BLIND_PLAN_V2",
    "repository": "/path/to/ginan",
    "static_inputs": {
      "code.driver": "scripts/experiment0_2024_blind_driver.py",
      "code.spatial": "scripts/experiment0_2024_spatial_validation.py",
      "code.parser": "scripts/validate_stec_satellite_difference_covariance.py",
      "code.provenance_helper": "scripts/experiment0_2024_blind_provenance.py",
      "manifest.input": "/data/input_manifest.json",
      "audit.rinex_qc": "/data/rinex_qc.json",
      "audit.split": "/data/split_audit.json",
      "audit.product": "/data/experiment0_2024_product_audit.json",
      "audit.validation_split": "/data/validation_split_audit.json",
      "yaml.quiet_model": "...yaml",
      "yaml.storm_model": "...yaml",
      "yaml.quiet_heldout": "...yaml",
      "yaml.storm_heldout": "...yaml",
      "yaml.include.common": "...yaml"
    },
    "model": {
      "quiet": {"stec": "...", "raw_cov": "...", "sd": "...",
                  "raw_audit": "...", "sd_audit": "...", "console": "..."},
      "storm": {"...": "same fields"}
    },
    "pea": {
      "binary": "bin/pea", "commit": "1296966",
      "environment": {"LD_LIBRARY_PATH": "/path/to/lib"},
      "quiet_argv": ["./bin/pea", "-y", "...quiet...yaml", "-a", "EXP0_OUTPUT_ROOT:/new/quiet"],
      "storm_argv": ["./bin/pea", "-y", "...storm...yaml", "-a", "EXP0_OUTPUT_ROOT:/new/storm"]
    },
    "heldout": {
      "quiet": {
        "output_root": "/new/quiet", "stec": "...STEC",
        "raw_cov": "...STEC.COV", "sd": "...STEC.SD.COV",
        "raw_audit": "...raw.audit.json", "sd_audit": "...sd.audit.json",
        "console": "...console.log", "run_receipt": "...run.receipt.json"
      },
      "storm": {"...": "same fields"}
    },
    "audit_policy": {
      "quiet": {"expected_invalid_epochs": []},
      "storm": {"expected_invalid_epochs": [{"gps_week": 2313, "gps_tow": 531420}]}
    },
    "freeze_output": "/data/freeze.json",
    "freeze_receipt_output": "/data/freeze.receipt.json",
    "report_output": "/data/report.json",
    "evaluation_receipt_output": "/data/report.receipt.json",
    "evaluation_claim_output": "/data/report.claim.json"
  }

All ``yaml.*`` paths must equal the recursively discovered include closure of
the four named entry points.  Relative include paths use the repository root,
matching these registered Ginan configurations.  The PEA environment is the
exact declared environment map; no implicit inherited variables are claimed.

Security boundary: formal execution assumes ext4 and trusted parent
directories that are not concurrently renamed or replaced.  The ``freeze``
subcommand does not run PEA and does not read held-out outputs.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import subprocess
import sys
import tempfile
import time
from datetime import datetime, timezone
from decimal import Decimal, InvalidOperation
from pathlib import Path
from typing import Mapping, Sequence

import yaml

import audit_experiment0_2024_products as product_audit
import audit_experiment0_2024_rinex as rinex_audit
import experiment0_2024_blind_provenance as provenance
import experiment0_2024_spatial_validation as spatial
import validate_experiment0_2024_split as split_audit
import validate_experiment0_2024_validation_split as validation_split_audit
import validate_stec_satellite_difference_covariance as sd_parser


PLAN_SCHEMA = "GINAN_EXPERIMENT0_2024_BLIND_PLAN_V2"
FREEZE_RECEIPT_SCHEMA = "GINAN_EXPERIMENT0_2024_BLIND_FREEZE_RECEIPT_V2"
RUN_RECEIPT_SCHEMA = "GINAN_EXPERIMENT0_2024_BLIND_RUN_RECEIPT_V1"
EVALUATION_RECEIPT_SCHEMA = "GINAN_EXPERIMENT0_2024_BLIND_EVALUATION_RECEIPT_V2"
EVALUATION_CLAIM_SCHEMA = "GINAN_EXPERIMENT0_2024_BLIND_EVALUATION_CLAIM_V1"
EXPECTED_POSTERIOR_STAGE = "FILTER_POSTERIOR_NO_EPOCH_AR"
DAY_LABELS = ("quiet", "storm")
MODEL_FIELDS = (
    "stec",
    "raw_cov",
    "sd",
    "raw_audit",
    "sd_audit",
    "console",
)
HELDOUT_ARTIFACT_FIELDS = MODEL_FIELDS
HELDOUT_FIELDS = HELDOUT_ARTIFACT_FIELDS + ("run_receipt",)
EXPECTED_EPOCH_COUNT = 2880
MAX_STORM_NO_STATES = 10
EXPECTED_MANIFEST_FILE_COUNT = 53
EXPECTED_MANIFEST_SCHEMA = "GINAN_EXPERIMENT0_2024_INPUT_MANIFEST_V1"
EXPECTED_VALIDATION_SPLIT_SCHEMA = "GINAN_EXPERIMENT0_2024_VALIDATION_SPLIT_AUDIT_V1"
EXPECTED_RINEX_QC_SCHEMA = "GINAN_EXPERIMENT0_2024_RINEX_QC_V1"
EXPECTED_SPLIT_AUDIT_SCHEMA = "GINAN_EXPERIMENT0_2024_SPLIT_AUDIT_V1"
EXPECTED_PRODUCT_AUDIT_SCHEMA = "GINAN_EXPERIMENT0_2024_PRODUCT_AUDIT_V1"
EXPECTED_DATES = {"quiet": "2024-05-08", "storm": "2024-05-11"}
EXPECTED_DAY_START_KEYS = {
    "quiet": (2313, Decimal(259200)),
    "storm": (2313, Decimal(518400)),
}
EXPECTED_MODEL_STATIONS = {
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
}
EXPECTED_HELDOUT_STATIONS = {"BALL", "PARK", "STR2"}
EXPECTED_MANIFEST_ENTRY_FIELDS = {
    "experiment",
    "relative_path",
    "role",
    "sha256",
    "size",
    "source",
    "source_metadata",
}
TIME_BINARY = Path("/usr/bin/time")
REQUIRED_STATIC_ROLES = {
    "code.driver",
    "code.spatial",
    "code.parser",
    "code.provenance_helper",
    "code.raw_validator",
    "code.rinex_auditor",
    "code.split_validator",
    "code.product_auditor",
    "code.validation_split_validator",
    "manifest.input",
    "audit.rinex_qc",
    "audit.split",
    "audit.product",
    "audit.validation_split",
    "yaml.quiet_model",
    "yaml.storm_model",
    "yaml.quiet_heldout",
    "yaml.storm_heldout",
    "yaml.quiet_model_inputs",
    "yaml.storm_model_inputs",
}
YAML_ENTRY_ROLES = {
    "quiet_model": "yaml.quiet_model",
    "storm_model": "yaml.storm_model",
    "quiet_heldout": "yaml.quiet_heldout",
    "storm_heldout": "yaml.storm_heldout",
}
ACTUAL_CODE_PATHS = {
    "code.driver": Path(__file__).resolve(),
    "code.spatial": Path(spatial.__file__).resolve(),
    "code.parser": Path(sd_parser.__file__).resolve(),
    "code.provenance_helper": Path(provenance.__file__).resolve(),
    "code.raw_validator": (
        Path(__file__).resolve().parent / "validate_stec_covariance.py"
    ),
    "code.rinex_auditor": Path(rinex_audit.__file__).resolve(),
    "code.split_validator": Path(split_audit.__file__).resolve(),
    "code.product_auditor": Path(product_audit.__file__).resolve(),
    "code.validation_split_validator": Path(validation_split_audit.__file__).resolve(),
}
_SHA256_PATTERN = re.compile(r"^[0-9a-fA-F]{64}$")
_COMMIT_PATTERN = re.compile(r"^[0-9a-fA-F]{7,40}$")
_EXIT_STATUS_PATTERN = re.compile(r"\bExit\s+status\s*[:=]\s*(-?\d+)\b", re.I)


class BlindDriverError(RuntimeError):
    """Raised when the blind execution contract is violated."""


def _absolute(path: str | os.PathLike[str], base: Path | None = None) -> Path:
    candidate = Path(path)
    if not candidate.is_absolute():
        if base is None:
            base = Path.cwd()
        candidate = base / candidate
    return Path(os.path.abspath(candidate))


def _mapping(value: object, label: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise BlindDriverError(f"{label} must be a JSON object")
    return value


def _string(value: object, label: str) -> str:
    if not isinstance(value, str) or not value:
        raise BlindDriverError(f"{label} must be a non-empty string")
    return value


def _reject_json_constant(value: str) -> object:
    raise ValueError(f"non-standard JSON constant {value}")


def _load_json_object(path: Path, label: str) -> dict[str, object]:
    before = provenance.record_regular_file(path)
    try:
        encoded = path.read_bytes()
        payload = json.loads(
            encoded.decode("utf-8"), parse_constant=_reject_json_constant
        )
    except (OSError, UnicodeError, json.JSONDecodeError, ValueError) as exc:
        raise BlindDriverError(f"cannot read {label} {path}: {exc}") from exc
    after = provenance.record_regular_file(path)
    if (
        before != after
        or len(encoded) != before["size"]
        or hashlib.sha256(encoded).hexdigest() != before["sha256"]
    ):
        raise BlindDriverError(f"{label} changed while it was read: {path}")
    if not isinstance(payload, dict):
        raise BlindDriverError(f"{label} must contain a JSON object")
    return payload


def _normalise_string_map(value: object, label: str) -> dict[str, str]:
    source = _mapping(value, label)
    output: dict[str, str] = {}
    for key, item in source.items():
        if not isinstance(key, str) or not key:
            raise BlindDriverError(f"{label} contains an invalid key")
        output[key] = _string(item, f"{label}.{key}")
    return output


def _normalise_argv(value: object, label: str) -> list[str]:
    if not isinstance(value, list) or not value:
        raise BlindDriverError(f"{label} must be a non-empty JSON array")
    if any(not isinstance(item, str) or not item for item in value):
        raise BlindDriverError(f"{label} must contain only non-empty strings")
    return list(value)


def _normalise_registered_epoch_list(
    value: object, label: str
) -> list[dict[str, object]]:
    if not isinstance(value, list):
        raise BlindDriverError(f"{label} must be a JSON array")
    output: list[dict[str, object]] = []
    keys: set[tuple[int, Decimal]] = set()
    for index, item in enumerate(value):
        entry = _mapping(item, f"{label}[{index}]")
        if set(entry) != {"gps_week", "gps_tow"}:
            raise BlindDriverError(f"{label}[{index}] must contain gps_week/gps_tow")
        week = entry.get("gps_week")
        if not isinstance(week, int) or isinstance(week, bool) or week < 0:
            raise BlindDriverError(f"{label}[{index}].gps_week is invalid")
        try:
            tow = Decimal(str(entry.get("gps_tow")))
        except InvalidOperation as exc:
            raise BlindDriverError(f"{label}[{index}].gps_tow is invalid") from exc
        if not tow.is_finite() or tow < 0 or tow >= Decimal(604800):
            raise BlindDriverError(f"{label}[{index}].gps_tow is invalid")
        key = (week, tow)
        if key in keys:
            raise BlindDriverError(f"{label} repeats epoch {week}:{tow}")
        keys.add(key)
        output.append({"gps_week": week, "gps_tow": str(tow)})
    return output


def _exact_aliases(argv: Sequence[str], label: str) -> dict[str, str]:
    aliases: dict[str, str] = {}
    index = 0
    while index < len(argv):
        if argv[index] != "-a":
            index += 1
            continue
        if index + 1 >= len(argv):
            raise BlindDriverError(f"{label} has trailing -a without KEY:VALUE")
        declaration = argv[index + 1]
        if ":" not in declaration:
            raise BlindDriverError(f"{label} alias is not KEY:VALUE: {declaration!r}")
        key, value = declaration.split(":", 1)
        if not key or not value:
            raise BlindDriverError(f"{label} alias is not KEY:VALUE: {declaration!r}")
        if key in aliases:
            raise BlindDriverError(f"{label} repeats alias {key}")
        aliases[key] = value
        index += 2
    return aliases


def _repo_path(value: object, repository: Path, label: str) -> Path:
    return _absolute(_string(value, label), repository)


def _manifest_expanded_inputs(manifest_path: Path) -> dict[str, str]:
    manifest = _load_json_object(manifest_path, "input manifest")
    if manifest.get("schema") != EXPECTED_MANIFEST_SCHEMA:
        raise BlindDriverError(
            f"input manifest schema must be {EXPECTED_MANIFEST_SCHEMA}"
        )
    if manifest.get("dates") != EXPECTED_DATES:
        raise BlindDriverError(
            f"input manifest dates must be {EXPECTED_DATES}, got {manifest.get('dates')!r}"
        )
    if set(manifest.get("model_stations", [])) != EXPECTED_MODEL_STATIONS:
        raise BlindDriverError("input manifest model station set mismatch")
    if set(manifest.get("heldout_stations", [])) != EXPECTED_HELDOUT_STATIONS:
        raise BlindDriverError("input manifest held-out station set mismatch")
    if set(manifest.get("stations", [])) != (
        EXPECTED_MODEL_STATIONS | EXPECTED_HELDOUT_STATIONS
    ):
        raise BlindDriverError("input manifest complete station set mismatch")
    if manifest.get("observation_policy") != {
        "constellation": "GPS",
        "nominal_interval_seconds": 30,
        "signals": ["C1C", "L1C", "C2W", "L2W"],
    }:
        raise BlindDriverError("input manifest observation policy mismatch")
    if manifest.get("precise_product_policy") != {
        "analysis_center": "Wuhan University",
        "antenna_calibration": "shared/products/igs20_2303.atx",
        "bias_apc_model": "IGS20_2303.ATX",
        "series": "WUM0MGXFIN",
    }:
        raise BlindDriverError("input manifest precise-product policy mismatch")
    files = manifest.get("files")
    if not isinstance(files, list) or len(files) != EXPECTED_MANIFEST_FILE_COUNT:
        raise BlindDriverError(
            f"input manifest must contain exactly {EXPECTED_MANIFEST_FILE_COUNT} files"
        )
    root = manifest_path.resolve().parent
    inputs: dict[str, str] = {}
    seen_paths: set[Path] = set()
    for index, value in enumerate(files):
        entry = _mapping(value, f"input_manifest.files[{index}]")
        if set(entry) != EXPECTED_MANIFEST_ENTRY_FIELDS:
            raise BlindDriverError(
                f"input_manifest.files[{index}] fields differ from the registered schema"
            )
        if (
            entry.get("experiment") not in {"quiet", "storm", "shared", "selection"}
            or not isinstance(entry.get("role"), str)
            or not entry.get("role")
            or not isinstance(entry.get("source"), str)
            or not entry.get("source")
            or not isinstance(entry.get("source_metadata"), Mapping)
        ):
            raise BlindDriverError(
                f"input_manifest.files[{index}] has invalid provenance fields"
            )
        relative_text = _string(
            entry.get("relative_path"), f"input_manifest.files[{index}].relative_path"
        )
        relative = Path(relative_text)
        if relative.is_absolute() or ".." in relative.parts:
            raise BlindDriverError(
                f"manifest relative_path escapes manifest root: {relative_text!r}"
            )
        path = root / relative
        record = provenance.record_regular_file(path)
        resolved_path = Path(
            _string(record["resolved_path"], "manifest file resolved_path")
        )
        try:
            resolved_path.relative_to(root)
        except ValueError as exc:
            raise BlindDriverError(
                f"manifest relative_path escapes manifest root: {relative_text!r}"
            ) from exc
        if resolved_path in seen_paths:
            raise BlindDriverError(f"duplicate manifest file: {relative_text}")
        seen_paths.add(resolved_path)
        expected_size = entry.get("size")
        expected_sha = entry.get("sha256")
        if (
            not isinstance(expected_size, int)
            or expected_size < 0
            or not isinstance(expected_sha, str)
            or not _SHA256_PATTERN.fullmatch(expected_sha)
        ):
            raise BlindDriverError(
                f"manifest file has invalid size/SHA-256: {relative_text}"
            )
        if record["size"] != expected_size or record["sha256"] != expected_sha.lower():
            raise BlindDriverError(
                f"manifest file identity mismatch: {relative_text}; "
                f"expected size={expected_size}, sha256={expected_sha.lower()}, "
                f"got size={record['size']}, sha256={record['sha256']}"
            )
        inputs[f"manifest.file.{index:02d}:{relative.as_posix()}"] = str(path)
    return inputs


def _validate_validation_split_contract(
    plan: Mapping[str, object],
) -> dict[str, object]:
    static = _mapping(plan["static_inputs"], "static_inputs")
    path = Path(_string(static["audit.validation_split"], "audit.validation_split"))
    audit = _load_json_object(path, "validation split audit")
    if audit.get("schema") != EXPECTED_VALIDATION_SPLIT_SCHEMA:
        raise BlindDriverError("validation split audit schema mismatch")
    manifest_path = Path(_string(static["manifest.input"], "manifest.input")).resolve()
    if (
        audit.get("valid") is not True
        or audit.get("model_split_valid") is not True
        or audit.get("roles_identical") is not True
        or audit.get("errors") != []
        or audit.get("intersection") != []
        or Path(_string(audit.get("manifest"), "validation_split.manifest")).resolve()
        != manifest_path
        or set(audit.get("model_stations", [])) != EXPECTED_MODEL_STATIONS
        or set(audit.get("heldout_stations", [])) != EXPECTED_HELDOUT_STATIONS
    ):
        raise BlindDriverError("validation split audit failed the registered split")
    per_day = _mapping(audit.get("per_day"), "validation_split.per_day")
    for day in DAY_LABELS:
        day_audit = _mapping(per_day.get(day), f"validation_split.per_day.{day}")
        expected_yaml = Path(
            _string(static[f"yaml.{day}_heldout"], f"yaml.{day}_heldout")
        ).resolve()
        if (
            Path(
                _string(day_audit.get("yaml"), f"validation_split.{day}.yaml")
            ).resolve()
            != expected_yaml
            or day_audit.get("rnx_input_count") != len(EXPECTED_HELDOUT_STATIONS)
            or set(day_audit.get("stations", [])) != EXPECTED_HELDOUT_STATIONS
            or day_audit.get("model_input_count") != 0
            or day_audit.get("model_inputs") != []
            or day_audit.get("missing_heldout") != []
            or day_audit.get("unexpected_stations") != []
            or day_audit.get("paths_match_manifest") is not True
        ):
            raise BlindDriverError(
                f"validation split audit failed registered {day} held-out contract"
            )
    regenerated = validation_split_audit.build_validation_audit(
        manifest_path,
        Path(_string(static["yaml.quiet_model_inputs"], "yaml.quiet_model_inputs")),
        Path(_string(static["yaml.storm_model_inputs"], "yaml.storm_model_inputs")),
        Path(_string(static["yaml.quiet_heldout"], "yaml.quiet_heldout")),
        Path(_string(static["yaml.storm_heldout"], "yaml.storm_heldout")),
    )
    stored_without_time = dict(audit)
    regenerated_without_time = dict(regenerated)
    stored_without_time.pop("generated_utc", None)
    regenerated_without_time.pop("generated_utc", None)
    if regenerated_without_time != stored_without_time:
        raise BlindDriverError(
            "validation split audit does not reproduce from frozen manifest/YAML"
        )
    return {
        "schema": EXPECTED_VALIDATION_SPLIT_SCHEMA,
        "model_stations": sorted(EXPECTED_MODEL_STATIONS),
        "heldout_stations": sorted(EXPECTED_HELDOUT_STATIONS),
        "valid": True,
    }


def _validate_registered_audit_contracts(
    plan: Mapping[str, object],
) -> dict[str, object]:
    static = _mapping(plan["static_inputs"], "static_inputs")
    manifest_path = Path(_string(static["manifest.input"], "manifest.input")).resolve()
    manifest_record = provenance.record_regular_file(manifest_path)
    data_root = manifest_path.parent

    rinex = _load_json_object(
        Path(_string(static["audit.rinex_qc"], "audit.rinex_qc")), "RINEX QC audit"
    )
    rinex_files = rinex.get("files")
    if (
        rinex.get("schema") != EXPECTED_RINEX_QC_SCHEMA
        or rinex.get("pass") is not True
        or rinex.get("expected_file_count") != 26
        or rinex.get("audited_file_count") != 26
        or not isinstance(rinex_files, list)
        or len(rinex_files) != 26
        or any(
            not isinstance(entry, Mapping) or entry.get("pass") is not True
            for entry in rinex_files
        )
        or Path(
            _string(rinex.get("input_manifest"), "rinex_qc.input_manifest")
        ).resolve()
        != manifest_path
        or rinex.get("input_manifest_sha256") != manifest_record["sha256"]
    ):
        raise BlindDriverError("RINEX QC audit failed the registered 26-file policy")
    regenerated_rinex = rinex_audit.build_audit(data_root)
    if regenerated_rinex != rinex:
        raise BlindDriverError(
            "RINEX QC audit does not reproduce from the frozen manifest/files"
        )

    split = _load_json_object(
        Path(_string(static["audit.split"], "audit.split")), "model split audit"
    )
    split_days = _mapping(split.get("per_day"), "split.per_day")
    if (
        split.get("schema") != EXPECTED_SPLIT_AUDIT_SCHEMA
        or split.get("valid") is not True
        or split.get("dates") != EXPECTED_DATES
        or set(split.get("train", [])) != EXPECTED_MODEL_STATIONS
        or set(split.get("heldout", [])) != EXPECTED_HELDOUT_STATIONS
        or split.get("intersection") != []
        or split.get("roles_identical") is not True
        or split.get("errors") != []
        or Path(_string(split.get("manifest"), "split.manifest")).resolve()
        != manifest_path
    ):
        raise BlindDriverError("model split audit failed the registered 10+3 split")
    for day in DAY_LABELS:
        day_split = _mapping(split_days.get(day), f"split.per_day.{day}")
        expected_yaml = Path(
            _string(static[f"yaml.{day}_model_inputs"], f"yaml.{day}_model_inputs")
        ).resolve()
        if (
            Path(_string(day_split.get("yaml"), f"split.{day}.yaml")).resolve()
            != expected_yaml
            or day_split.get("rnx_input_count") != len(EXPECTED_MODEL_STATIONS)
            or set(day_split.get("train", [])) != EXPECTED_MODEL_STATIONS
            or day_split.get("heldout_input_count") != 0
            or day_split.get("heldout_inputs") != []
            or day_split.get("duplicate_stations") != []
            or day_split.get("invalid_inputs") != []
            or day_split.get("missing_train") != []
            or day_split.get("unexpected_stations") != []
        ):
            raise BlindDriverError(
                f"model split audit failed registered {day} model contract"
            )
    regenerated_split = split_audit.build_audit(
        manifest_path,
        Path(_string(static["yaml.quiet_model_inputs"], "yaml.quiet_model_inputs")),
        Path(_string(static["yaml.storm_model_inputs"], "yaml.storm_model_inputs")),
    )
    split_without_time = dict(split)
    regenerated_split_without_time = dict(regenerated_split)
    split_without_time.pop("generated_utc", None)
    regenerated_split_without_time.pop("generated_utc", None)
    if regenerated_split_without_time != split_without_time:
        raise BlindDriverError(
            "model split audit does not reproduce from frozen manifest/YAML"
        )

    product = _load_json_object(
        Path(_string(static["audit.product"], "audit.product")), "product audit"
    )
    experiments = _mapping(product.get("experiments"), "product.experiments")
    if (
        product.get("schema") != EXPECTED_PRODUCT_AUDIT_SCHEMA
        or product.get("validation_passed") is not False
        or set(experiments) != set(DAY_LABELS)
    ):
        raise BlindDriverError("product audit overall registered status mismatch")
    expected_defect_counts = {"quiet": 1, "storm": 31}
    for day in DAY_LABELS:
        experiment = _mapping(experiments.get(day), f"product.experiments.{day}")
        validation = _mapping(
            experiment.get("validation"), f"product.experiments.{day}.validation"
        )
        metadata = _mapping(
            validation.get("metadata"), f"product.experiments.{day}.metadata"
        )
        defects = validation.get("eligible_satellites_with_coverage_defects")
        if (
            validation.get("eligible_intersection_passed") is not True
            or validation.get("full_day_completeness_passed") is not False
            or validation.get("passed") is not False
            or validation.get("missing_expected_satellites") != []
            or validation.get("unexpected_satellites") != []
            or metadata.get("passed") is not True
            or not isinstance(defects, list)
            or len(defects) != expected_defect_counts[day]
        ):
            raise BlindDriverError(
                f"product audit {day} status differs from registered CLK-gap policy"
            )
    regenerated_product = product_audit.audit_data_root(data_root)
    product_without_time = dict(product)
    regenerated_product_without_time = dict(regenerated_product)
    product_without_time.pop("generated_utc", None)
    regenerated_product_without_time.pop("generated_utc", None)
    if regenerated_product_without_time != product_without_time:
        raise BlindDriverError(
            "product audit does not reproduce from the frozen product files"
        )
    return {
        "rinex_qc": {"pass": True, "file_count": 26},
        "model_split": {"valid": True, "model": 10, "heldout": 3},
        "product": {
            "metadata_and_eligibility_passed": True,
            "known_clk_gap_defect_counts": expected_defect_counts,
            "full_day_completeness_passed": False,
        },
    }


def _discover_yaml_closure(
    repository: Path, entrypoints: Mapping[str, Path]
) -> set[Path]:
    visited: set[Path] = set()
    active: set[Path] = set()

    def visit(path: Path) -> None:
        path = path.resolve(strict=True)
        if path in active:
            raise BlindDriverError(f"YAML include cycle detected at {path}")
        if path in visited:
            return
        provenance.record_regular_file(path)
        active.add(path)
        try:
            document = yaml.safe_load(path.read_text(encoding="utf-8"))
        except (OSError, UnicodeError, yaml.YAMLError) as exc:
            raise BlindDriverError(f"cannot parse YAML {path}: {exc}") from exc
        if document is None:
            document = {}
        if not isinstance(document, Mapping):
            raise BlindDriverError(f"YAML root is not a mapping: {path}")
        inputs = document.get("inputs", {})
        if inputs is None:
            inputs = {}
        if not isinstance(inputs, Mapping):
            raise BlindDriverError(f"inputs is not a mapping in {path}")
        includes = inputs.get("include_yamls", [])
        if includes is None:
            includes = []
        if not isinstance(includes, list) or any(
            not isinstance(item, str) or not item for item in includes
        ):
            raise BlindDriverError(f"inputs.include_yamls is invalid in {path}")
        for include in includes:
            included_path = _absolute(include, repository)
            try:
                included_path.resolve(strict=True).relative_to(repository)
            except (FileNotFoundError, ValueError) as exc:
                raise BlindDriverError(
                    f"YAML include must exist inside repository: {include!r} from {path}"
                ) from exc
            visit(included_path)
        active.remove(path)
        visited.add(path)

    for entrypoint in entrypoints.values():
        visit(entrypoint)
    return visited


def _normalise_freeze_plan(plan_path: Path) -> dict[str, object]:
    raw = _load_json_object(plan_path, "freeze plan")
    if raw.get("schema") != PLAN_SCHEMA:
        raise BlindDriverError(f"freeze plan schema must be {PLAN_SCHEMA}")
    repository = _repo_path(raw.get("repository"), Path.cwd(), "repository")
    try:
        repository = repository.resolve(strict=True)
    except FileNotFoundError as exc:
        raise BlindDriverError(f"repository does not exist: {repository}") from exc

    raw_static = _normalise_string_map(raw.get("static_inputs"), "static_inputs")
    missing_roles = sorted(REQUIRED_STATIC_ROLES - set(raw_static))
    if missing_roles:
        raise BlindDriverError(f"static_inputs missing roles: {missing_roles}")
    static_inputs = {
        role: str(_repo_path(path, repository, f"static_inputs.{role}"))
        for role, path in raw_static.items()
    }
    for role, actual in ACTUAL_CODE_PATHS.items():
        declared = Path(static_inputs[role]).resolve(strict=True)
        if declared != actual:
            raise BlindDriverError(
                f"{role} must name the executing code {actual}, got {declared}"
            )

    yaml_entries = {
        label: Path(static_inputs[role]) for label, role in YAML_ENTRY_ROLES.items()
    }
    discovered_yaml = _discover_yaml_closure(repository, yaml_entries)
    declared_yaml = {
        Path(path).resolve(strict=True)
        for role, path in static_inputs.items()
        if role.startswith("yaml.")
    }
    if discovered_yaml != declared_yaml:
        missing = sorted(str(path) for path in discovered_yaml - declared_yaml)
        extra = sorted(str(path) for path in declared_yaml - discovered_yaml)
        raise BlindDriverError(
            f"declared YAML inputs are not the complete include closure: "
            f"missing={missing}, extra={extra}"
        )

    raw_model = _mapping(raw.get("model"), "model")
    model: dict[str, dict[str, str]] = {}
    for day in DAY_LABELS:
        day_model = _mapping(raw_model.get(day), f"model.{day}")
        model[day] = {
            field: str(
                _repo_path(day_model.get(field), repository, f"model.{day}.{field}")
            )
            for field in MODEL_FIELDS
        }

    raw_pea = _mapping(raw.get("pea"), "pea")
    pea_binary = _repo_path(raw_pea.get("binary"), repository, "pea.binary")
    pea_commit = _string(raw_pea.get("commit"), "pea.commit").lower()
    if not _COMMIT_PATTERN.fullmatch(pea_commit):
        raise BlindDriverError("pea.commit must be a 7-40 digit hexadecimal commit")
    environment = _normalise_string_map(raw_pea.get("environment"), "pea.environment")
    pea: dict[str, object] = {
        "binary": str(pea_binary),
        "commit": pea_commit,
        "environment": environment,
        "environment_semantics": "exact declared map; no implicit inherited variables claimed",
        "working_directory": str(repository),
    }

    raw_heldout = _mapping(raw.get("heldout"), "heldout")
    data_root = Path(static_inputs["manifest.input"]).resolve(strict=True).parent
    heldout: dict[str, dict[str, object]] = {}
    for day in DAY_LABELS:
        day_plan = _mapping(raw_heldout.get(day), f"heldout.{day}")
        output_root = _repo_path(
            day_plan.get("output_root"), repository, f"heldout.{day}.output_root"
        )
        artifacts = {
            field: _repo_path(day_plan.get(field), repository, f"heldout.{day}.{field}")
            for field in HELDOUT_FIELDS
        }
        for field, artifact in artifacts.items():
            try:
                artifact.relative_to(output_root)
            except ValueError as exc:
                raise BlindDriverError(
                    f"heldout.{day}.{field} must be inside output_root {output_root}"
                ) from exc
        argv = _normalise_argv(raw_pea.get(f"{day}_argv"), f"pea.{day}_argv")
        executable = _absolute(argv[0], repository)
        if executable.resolve(strict=False) != pea_binary.resolve(strict=False):
            raise BlindDriverError(
                f"pea.{day}_argv[0] does not name pea.binary: {argv[0]}"
            )
        yaml_entry = Path(static_inputs[f"yaml.{day}_heldout"])
        yaml_arguments = [
            _absolute(argv[index + 1], repository).resolve(strict=False)
            for index, argument in enumerate(argv[:-1])
            if argument == "-y"
        ]
        if yaml_arguments != [yaml_entry.resolve(strict=True)]:
            raise BlindDriverError(
                f"pea.{day}_argv must contain exactly one '-y {yaml_entry}'"
            )
        aliases = _exact_aliases(argv, f"pea.{day}_argv")
        required_aliases = {
            "EXP0_DATA_ROOT": str(data_root),
            "EXP0_OUTPUT_ROOT": str(output_root),
        }
        for alias, expected in required_aliases.items():
            if aliases.get(alias) != expected:
                raise BlindDriverError(
                    f"pea.{day}_argv must declare exactly one {alias}:{expected}; "
                    f"got {aliases.get(alias)!r}"
                )
        if set(aliases) != set(required_aliases):
            raise BlindDriverError(
                f"pea.{day}_argv aliases must be exactly {sorted(required_aliases)}"
            )
        if (
            sum(argument == "-a" for argument in argv) != len(aliases)
            or len(argv) != 1 + 2 + 2 * len(aliases)
            or any(
                argument not in {"-y", "-a"}
                for argument in argv[1:]
                if not argument.startswith("EXP0_")
                and _absolute(argument, repository).resolve(strict=False)
                != yaml_entry.resolve(strict=True)
            )
        ):
            raise BlindDriverError(
                f"pea.{day}_argv contains an extra or malformed argument"
            )
        pea[f"{day}_argv"] = argv
        heldout[day] = {
            "output_root": str(output_root),
            **{field: str(path) for field, path in artifacts.items()},
        }

    quiet_root = Path(heldout["quiet"]["output_root"])
    storm_root = Path(heldout["storm"]["output_root"])
    if quiet_root == storm_root:
        raise BlindDriverError("quiet and storm held-out output roots must differ")

    raw_audit_policy = _mapping(raw.get("audit_policy"), "audit_policy")
    audit_policy = {
        day: {
            "expected_invalid_epochs": _normalise_registered_epoch_list(
                _mapping(raw_audit_policy.get(day), f"audit_policy.{day}").get(
                    "expected_invalid_epochs"
                ),
                f"audit_policy.{day}.expected_invalid_epochs",
            )
        }
        for day in DAY_LABELS
    }

    freeze_output = _repo_path(raw.get("freeze_output"), repository, "freeze_output")
    freeze_receipt_output = _repo_path(
        raw.get("freeze_receipt_output"), repository, "freeze_receipt_output"
    )
    report_output = _repo_path(raw.get("report_output"), repository, "report_output")
    evaluation_receipt_output = _repo_path(
        raw.get("evaluation_receipt_output"),
        repository,
        "evaluation_receipt_output",
    )
    evaluation_claim_output = _repo_path(
        raw.get("evaluation_claim_output"), repository, "evaluation_claim_output"
    )
    final_paths = {
        freeze_output,
        freeze_receipt_output,
        report_output,
        evaluation_receipt_output,
        evaluation_claim_output,
    }
    if len(final_paths) != 5:
        raise BlindDriverError("freeze/report and receipt paths must all differ")

    return {
        "schema": PLAN_SCHEMA,
        "plan_path": str(plan_path.resolve(strict=True)),
        "repository": str(repository),
        "static_inputs": static_inputs,
        "model": model,
        "pea": pea,
        "heldout": heldout,
        "audit_policy": audit_policy,
        "freeze_output": str(freeze_output),
        "freeze_receipt_output": str(freeze_receipt_output),
        "report_output": str(report_output),
        "evaluation_receipt_output": str(evaluation_receipt_output),
        "evaluation_claim_output": str(evaluation_claim_output),
    }


def _freeze_inputs(plan: Mapping[str, object]) -> dict[str, str]:
    inputs = dict(_mapping(plan["static_inputs"], "static_inputs"))
    manifest_path = Path(_string(inputs["manifest.input"], "manifest.input"))
    inputs.update(_manifest_expanded_inputs(manifest_path))
    inputs["plan.freeze"] = _string(plan["plan_path"], "plan_path")
    pea = _mapping(plan["pea"], "pea")
    inputs["executable.pea"] = _string(pea["binary"], "pea.binary")
    inputs["runtime.time"] = str(TIME_BINARY)
    model = _mapping(plan["model"], "model")
    for day in DAY_LABELS:
        day_model = _mapping(model[day], f"model.{day}")
        for field in MODEL_FIELDS:
            inputs[f"model.{day}.{field}"] = _string(
                day_model[field], f"model.{day}.{field}"
            )
    return inputs


def _freeze_planned_outputs(plan: Mapping[str, object]) -> dict[str, str]:
    outputs = {
        "freeze": _string(plan["freeze_output"], "freeze_output"),
        "freeze_receipt": _string(
            plan["freeze_receipt_output"], "freeze_receipt_output"
        ),
        "report": _string(plan["report_output"], "report_output"),
        "evaluation_receipt": _string(
            plan["evaluation_receipt_output"], "evaluation_receipt_output"
        ),
        "evaluation_claim": _string(
            plan["evaluation_claim_output"], "evaluation_claim_output"
        ),
    }
    heldout = _mapping(plan["heldout"], "heldout")
    for day in DAY_LABELS:
        day_plan = _mapping(heldout[day], f"heldout.{day}")
        outputs[f"heldout.{day}.output_root"] = _string(
            day_plan["output_root"], f"heldout.{day}.output_root"
        )
        for field in HELDOUT_FIELDS:
            outputs[f"heldout.{day}.{field}"] = _string(
                day_plan[field], f"heldout.{day}.{field}"
            )
    return outputs


def freeze_from_plan(plan_path: Path) -> dict[str, object]:
    """Build and publish one blind freeze plus its receipt."""

    plan = _normalise_freeze_plan(plan_path)
    planned_outputs = _freeze_planned_outputs(plan)
    provenance.assert_planned_outputs_absent(planned_outputs)
    pre = provenance.begin_provenance(
        plan["repository"], _freeze_inputs(plan), planned_outputs
    )

    model = _mapping(plan["model"], "model")
    quiet_model = _mapping(model["quiet"], "model.quiet")
    static_inputs = _mapping(plan["static_inputs"], "static_inputs")
    model_validation = {
        day: _validate_day_audits(
            day,
            _mapping(model[day], f"model.{day}"),
            f"model.{day}",
            _registered_invalid_epoch_keys(plan, day),
        )
        for day in DAY_LABELS
    }
    split_validation = _validate_validation_split_contract(plan)
    registered_audit_validation = _validate_registered_audit_contracts(plan)
    payload = spatial.build_freeze(
        Path(_string(quiet_model["stec"], "model.quiet.stec")),
        Path(_string(quiet_model["sd"], "model.quiet.sd")),
        Path(_string(static_inputs["audit.product"], "audit.product")),
        Path(
            _string(static_inputs["audit.validation_split"], "audit.validation_split")
        ),
    )
    if payload.get("heldout_output_read") is not False:
        raise BlindDriverError(
            "scientific freeze did not declare heldout_output_read=false"
        )
    completed = provenance.complete_provenance(pre, _freeze_inputs(plan))
    payload = dict(payload)
    payload["blind_run_plan"] = plan
    payload["blind_provenance"] = completed
    payload["blind_model_validation"] = model_validation
    payload["blind_split_validation"] = split_validation
    payload["blind_registered_audit_validation"] = registered_audit_validation

    freeze_path = Path(_string(plan["freeze_output"], "freeze_output"))
    freeze_record = provenance.atomic_write_json_no_replace(freeze_path, payload)
    receipt_payload = {
        "schema": FREEZE_RECEIPT_SCHEMA,
        "freeze": freeze_record,
        "run_plan": plan,
        "provenance": completed,
    }
    receipt_path = Path(_string(plan["freeze_receipt_output"], "freeze_receipt_output"))
    receipt_record = provenance.atomic_write_json_no_replace(
        receipt_path, receipt_payload
    )
    return {
        "freeze": freeze_record,
        "freeze_receipt": receipt_record,
        "heldout_output_read": False,
    }


def _anchored_freeze(
    freeze_path: Path, receipt_path: Path, expected_sha256: str
) -> tuple[dict[str, object], dict[str, object], dict[str, object]]:
    if not _SHA256_PATTERN.fullmatch(expected_sha256):
        raise BlindDriverError("expected freeze SHA-256 must contain 64 hex digits")
    expected_sha256 = expected_sha256.lower()
    freeze_record = provenance.record_regular_file(freeze_path)
    if freeze_record["sha256"] != expected_sha256:
        raise BlindDriverError(
            f"freeze SHA-256 mismatch: expected {expected_sha256}, "
            f"got {freeze_record['sha256']}"
        )
    freeze = _load_json_object(freeze_path, "freeze")
    receipt = _load_json_object(receipt_path, "freeze receipt")
    if receipt.get("schema") != FREEZE_RECEIPT_SCHEMA:
        raise BlindDriverError("invalid freeze receipt schema")
    receipt_freeze = _mapping(receipt.get("freeze"), "receipt.freeze")
    if receipt_freeze.get("sha256") != expected_sha256:
        raise BlindDriverError("external freeze SHA-256 does not match freeze receipt")
    if (
        Path(_string(receipt_freeze.get("path"), "receipt.freeze.path"))
        != freeze_path.resolve()
    ):
        raise BlindDriverError("freeze path does not match freeze receipt")
    plan = _mapping(receipt.get("run_plan"), "receipt.run_plan")
    if (
        Path(_string(plan.get("freeze_output"), "run_plan.freeze_output"))
        != freeze_path.resolve()
    ):
        raise BlindDriverError("freeze path differs from frozen run plan")
    if (
        Path(
            _string(plan.get("freeze_receipt_output"), "run_plan.freeze_receipt_output")
        )
        != receipt_path.resolve()
    ):
        raise BlindDriverError("freeze receipt path differs from frozen run plan")
    if freeze.get("blind_run_plan") != plan:
        raise BlindDriverError("freeze and receipt run plans differ")
    frozen_provenance = _mapping(receipt.get("provenance"), "receipt.provenance")
    if freeze.get("blind_provenance") != frozen_provenance:
        raise BlindDriverError("freeze and receipt provenance differ")
    if frozen_provenance.get("inputs_unchanged") is not True:
        raise BlindDriverError("freeze receipt does not prove unchanged inputs")
    return freeze, dict(receipt), dict(plan)


def _verify_frozen_static_state(
    plan: Mapping[str, object], receipt: Mapping[str, object]
) -> None:
    frozen = _mapping(receipt.get("provenance"), "receipt.provenance")
    post = _mapping(frozen.get("inputs_post"), "provenance.inputs_post")
    provenance.assert_input_snapshot_unchanged(post)
    current_git = provenance.git_identity(_string(plan["repository"], "repository"))
    if current_git != frozen.get("git"):
        raise BlindDriverError("git identity differs from frozen receipt")
    current_runtime = provenance.runtime_identity()
    if current_runtime != frozen.get("runtime"):
        raise BlindDriverError("Python/NumPy runtime differs from frozen receipt")


def _frozen_file_record(
    freeze_receipt: Mapping[str, object], role: str
) -> Mapping[str, object]:
    frozen = _mapping(freeze_receipt.get("provenance"), "receipt.provenance")
    inputs = _mapping(frozen.get("inputs_post"), "provenance.inputs_post")
    return _mapping(inputs.get(role), f"provenance.inputs_post.{role}")


def _execute_logged(
    command: Sequence[str],
    working_directory: Path,
    environment: Mapping[str, str],
    console_path: Path,
    *,
    exclusive: bool,
) -> dict[str, object]:
    console_path.parent.mkdir(parents=True, exist_ok=True)
    mode = "xb" if exclusive else "ab"
    started = datetime.now(timezone.utc)
    started_monotonic = time.monotonic()
    with console_path.open(mode) as stream:
        completed = subprocess.run(
            list(command),
            cwd=working_directory,
            env=dict(environment),
            stdout=stream,
            stderr=subprocess.STDOUT,
            check=False,
        )
    ended = datetime.now(timezone.utc)
    return {
        "command": list(command),
        "working_directory": str(working_directory.resolve()),
        "started_utc": started.isoformat(),
        "ended_utc": ended.isoformat(),
        "duration_seconds": time.monotonic() - started_monotonic,
        "returncode": completed.returncode,
    }


def _publish_existing_file_no_replace(source: Path, destination: Path) -> None:
    if os.path.lexists(destination):
        raise FileExistsError(f"planned output already exists: {destination}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    try:
        os.link(source, destination, follow_symlinks=False)
    except BaseException:
        source.unlink(missing_ok=True)
        raise
    source.unlink(missing_ok=True)


def _capture_available_artifacts(
    day_plan: Mapping[str, object], day: str
) -> dict[str, object]:
    records: dict[str, object] = {}
    for field in HELDOUT_ARTIFACT_FIELDS:
        path = Path(_string(day_plan[field], f"heldout.{day}.{field}"))
        if not os.path.lexists(path):
            records[field] = {"path": str(path), "status": "ABSENT"}
            continue
        try:
            records[field] = provenance.record_regular_file(path)
        except Exception as exc:  # failure receipt must survive a malformed artifact
            records[field] = {
                "path": str(path),
                "status": "UNREADABLE",
                "error_type": type(exc).__name__,
                "error": str(exc),
            }
    return records


def _write_post_pea_failure_receipt(
    *,
    path: Path,
    day: str,
    freeze_path: Path,
    freeze_sha256: str,
    pea: Mapping[str, object],
    argv: Sequence[str],
    environment: Mapping[str, str],
    binary_record: Mapping[str, object],
    time_record: Mapping[str, object],
    pea_execution: Mapping[str, object],
    validators: Mapping[str, object],
    day_plan: Mapping[str, object],
    stage: str,
    error: BaseException | str,
    status: str = "VALIDATION_FAILED",
) -> dict[str, object]:
    error_text = str(error)
    error_type = type(error).__name__ if isinstance(error, BaseException) else "Error"
    payload = {
        "schema": RUN_RECEIPT_SCHEMA,
        "status": status,
        "day": day,
        "freeze": {
            "path": str(freeze_path.resolve()),
            "sha256": freeze_sha256.lower(),
        },
        "pea": {
            "commit": pea["commit"],
            "argv": list(argv),
            "environment": dict(environment),
            "binary": dict(binary_record),
            "time_binary": dict(time_record),
        },
        "execution": dict(pea_execution),
        "validators": dict(validators),
        "failure": {"stage": stage, "error_type": error_type, "error": error_text},
        "available_artifacts": _capture_available_artifacts(day_plan, day),
    }
    return provenance.atomic_write_json_no_replace(path, payload)


def run_frozen_day(
    freeze_path: Path,
    freeze_receipt_path: Path,
    expected_freeze_sha256: str,
    day: str,
) -> dict[str, object]:
    """Execute one exact frozen held-out PEA command and publish its receipt."""

    if day not in DAY_LABELS:
        raise BlindDriverError(f"day must be one of {DAY_LABELS}")
    _, freeze_receipt, plan = _anchored_freeze(
        freeze_path.resolve(), freeze_receipt_path.resolve(), expected_freeze_sha256
    )
    _verify_frozen_static_state(plan, freeze_receipt)
    heldout = _mapping(plan["heldout"], "heldout")
    day_plan = _mapping(heldout[day], f"heldout.{day}")
    root = Path(_string(day_plan["output_root"], f"heldout.{day}.output_root"))
    planned = {"output_root": str(root)}
    planned.update(
        {
            field: _string(day_plan[field], f"heldout.{day}.{field}")
            for field in HELDOUT_FIELDS
        }
    )
    provenance.assert_planned_outputs_absent(planned)
    if not root.parent.is_dir():
        raise BlindDriverError(f"held-out output parent does not exist: {root.parent}")

    pea = _mapping(plan["pea"], "pea")
    argv = _normalise_argv(pea[f"{day}_argv"], f"pea.{day}_argv")
    environment = _normalise_string_map(pea["environment"], "pea.environment")
    working_directory = Path(_string(pea["working_directory"], "pea.working_directory"))
    time_binary = TIME_BINARY
    time_record = provenance.record_regular_file(time_binary)
    if time_record != _frozen_file_record(freeze_receipt, "runtime.time"):
        raise BlindDriverError("/usr/bin/time identity differs from the freeze")
    console = Path(_string(day_plan["console"], f"heldout.{day}.console"))
    run_receipt_path = Path(
        _string(day_plan["run_receipt"], f"heldout.{day}.run_receipt")
    )

    try:
        root.mkdir()
    except FileExistsError as exc:
        raise BlindDriverError(
            f"held-out output root was claimed concurrently: {root}"
        ) from exc

    temporary_descriptor, temporary_console_name = tempfile.mkstemp(
        prefix=f".{day}.pea-console.", suffix=".tmp", dir=root.parent
    )
    os.close(temporary_descriptor)
    temporary_console = Path(temporary_console_name)
    pea_command = [str(time_binary), "-v", "--", *argv]
    pea_execution: dict[str, object] = {
        "command": pea_command,
        "working_directory": str(working_directory.resolve()),
        "returncode": None,
    }
    try:
        pea_execution = _execute_logged(
            pea_command,
            working_directory,
            environment,
            temporary_console,
            exclusive=False,
        )
        if not root.is_dir():
            raise BlindDriverError(f"PEA output root is not a directory: {root}")
        _publish_existing_file_no_replace(temporary_console, console)
    except Exception as exc:
        temporary_console.unlink(missing_ok=True)
        receipt_record = _write_post_pea_failure_receipt(
            path=run_receipt_path,
            day=day,
            freeze_path=freeze_path,
            freeze_sha256=expected_freeze_sha256,
            pea=pea,
            argv=argv,
            environment=environment,
            binary_record=_frozen_file_record(freeze_receipt, "executable.pea"),
            time_record=time_record,
            pea_execution=pea_execution,
            validators={},
            day_plan=day_plan,
            stage="PEA_EXECUTION_OR_CONSOLE_PUBLICATION",
            error=exc,
            status="EXECUTION_FAILED",
        )
        raise BlindDriverError(
            f"{day} execution/console publication failed; failure receipt written to "
            f"{receipt_record['path']}: {exc}"
        ) from exc
    except BaseException:
        temporary_console.unlink(missing_ok=True)
        raise
    if pea_execution["returncode"] != 0:
        failed = {
            "schema": RUN_RECEIPT_SCHEMA,
            "status": "PEA_FAILED",
            "day": day,
            "freeze": {
                "path": str(freeze_path.resolve()),
                "sha256": expected_freeze_sha256.lower(),
            },
            "pea": {
                "commit": pea["commit"],
                "argv": argv,
                "environment": environment,
                "binary": dict(_frozen_file_record(freeze_receipt, "executable.pea")),
                "time_binary": time_record,
            },
            "execution": pea_execution,
        }
        receipt_record = provenance.atomic_write_json_no_replace(
            run_receipt_path, failed
        )
        raise BlindDriverError(
            f"{day} PEA returned {pea_execution['returncode']}; failure receipt "
            f"written to {receipt_record['path']}"
        )

    raw_cov = Path(_string(day_plan["raw_cov"], f"heldout.{day}.raw_cov"))
    sd_path = Path(_string(day_plan["sd"], f"heldout.{day}.sd"))
    raw_audit = Path(_string(day_plan["raw_audit"], f"heldout.{day}.raw_audit"))
    sd_audit = Path(_string(day_plan["sd_audit"], f"heldout.{day}.sd_audit"))
    provenance.assert_planned_outputs_absent(
        {"raw_audit": raw_audit, "sd_audit": sd_audit, "run_receipt": run_receipt_path}
    )
    raw_audit.parent.mkdir(parents=True, exist_ok=True)
    sd_audit.parent.mkdir(parents=True, exist_ok=True)
    raw_validator = Path(
        _string(
            _mapping(plan["static_inputs"], "static_inputs")["code.raw_validator"],
            "code.raw_validator",
        )
    )
    sd_validator = Path(
        _string(
            _mapping(plan["static_inputs"], "static_inputs")["code.parser"],
            "code.parser",
        )
    )
    validator_executions: dict[str, object] = {}
    try:
        raw_execution = _execute_logged(
            [
                sys.executable,
                str(raw_validator),
                str(raw_cov),
                "--json-output",
                str(raw_audit),
            ],
            working_directory,
            environment,
            console,
            exclusive=False,
        )
        validator_executions["raw"] = raw_execution
        sd_execution = _execute_logged(
            [
                sys.executable,
                str(sd_validator),
                str(raw_cov),
                str(sd_path),
                str(sd_audit),
            ],
            working_directory,
            environment,
            console,
            exclusive=False,
        )
        validator_executions["sd"] = sd_execution
    except Exception as exc:
        receipt_record = _write_post_pea_failure_receipt(
            path=run_receipt_path,
            day=day,
            freeze_path=freeze_path,
            freeze_sha256=expected_freeze_sha256,
            pea=pea,
            argv=argv,
            environment=environment,
            binary_record=_frozen_file_record(freeze_receipt, "executable.pea"),
            time_record=time_record,
            pea_execution=pea_execution,
            validators=validator_executions,
            day_plan=day_plan,
            stage="VALIDATOR_LAUNCH_OR_LOGGING",
            error=exc,
        )
        raise BlindDriverError(
            f"{day} validator launch/logging failed; failure receipt written to "
            f"{receipt_record['path']}: {exc}"
        ) from exc
    if raw_execution["returncode"] not in (0, 1) or sd_execution["returncode"] not in (
        0,
        1,
    ):
        error = (
            f"{day} validator fatal return codes: raw={raw_execution['returncode']}, "
            f"sd={sd_execution['returncode']}"
        )
        receipt_record = _write_post_pea_failure_receipt(
            path=run_receipt_path,
            day=day,
            freeze_path=freeze_path,
            freeze_sha256=expected_freeze_sha256,
            pea=pea,
            argv=argv,
            environment=environment,
            binary_record=_frozen_file_record(freeze_receipt, "executable.pea"),
            time_record=time_record,
            pea_execution=pea_execution,
            validators=validator_executions,
            day_plan=day_plan,
            stage="VALIDATOR_EXECUTION",
            error=error,
        )
        raise BlindDriverError(
            f"{error}; failure receipt written to {receipt_record['path']}"
        )
    try:
        validation = {
            "audits": _validate_day_audits(
                day,
                day_plan,
                f"heldout.{day}",
                _registered_invalid_epoch_keys(plan, day),
            ),
            "console": _successful_console(console, f"{day} console"),
            "sd_stage": _posterior_stage_summary(sd_path, f"{day} SD sidecar"),
        }
        frozen = _mapping(freeze_receipt["provenance"], "freeze_receipt.provenance")
        provenance.assert_input_snapshot_unchanged(
            _mapping(frozen["inputs_post"], "freeze_receipt.inputs_post")
        )
        artifact_paths = {
            field: _string(day_plan[field], f"heldout.{day}.{field}")
            for field in HELDOUT_ARTIFACT_FIELDS
        }
        artifacts = provenance.capture_input_snapshot(artifact_paths)
    except Exception as exc:
        receipt_record = _write_post_pea_failure_receipt(
            path=run_receipt_path,
            day=day,
            freeze_path=freeze_path,
            freeze_sha256=expected_freeze_sha256,
            pea=pea,
            argv=argv,
            environment=environment,
            binary_record=_frozen_file_record(freeze_receipt, "executable.pea"),
            time_record=time_record,
            pea_execution=pea_execution,
            validators=validator_executions,
            day_plan=day_plan,
            stage="POST_PEA_VALIDATION",
            error=exc,
        )
        raise BlindDriverError(
            f"{day} post-PEA validation failed; failure receipt written to "
            f"{receipt_record['path']}: {exc}"
        ) from exc
    receipt = {
        "schema": RUN_RECEIPT_SCHEMA,
        "status": "SUCCESS",
        "day": day,
        "freeze": {
            "path": str(freeze_path.resolve()),
            "sha256": expected_freeze_sha256.lower(),
        },
        "pea": {
            "commit": pea["commit"],
            "argv": argv,
            "environment": environment,
            "binary": dict(_frozen_file_record(freeze_receipt, "executable.pea")),
            "time_binary": time_record,
        },
        "execution": pea_execution,
        "validators": validator_executions,
        "validation": validation,
        "artifacts": artifacts,
    }
    receipt_record = provenance.atomic_write_json_no_replace(run_receipt_path, receipt)
    return {"day": day, "status": "SUCCESS", "run_receipt": receipt_record}


def _verify_day_run_receipt(
    day: str,
    plan: Mapping[str, object],
    freeze_receipt: Mapping[str, object],
    expected_freeze_sha256: str,
    expected_run_receipt_sha256: str,
) -> dict[str, object]:
    heldout = _mapping(plan["heldout"], "heldout")
    day_plan = _mapping(heldout[day], f"heldout.{day}")
    receipt_path = Path(_string(day_plan["run_receipt"], f"heldout.{day}.run_receipt"))
    receipt_record = provenance.record_regular_file(receipt_path)
    if (
        not _SHA256_PATTERN.fullmatch(expected_run_receipt_sha256)
        or receipt_record["sha256"] != expected_run_receipt_sha256.lower()
    ):
        raise BlindDriverError(f"{day} run receipt SHA-256 mismatch")
    receipt = _load_json_object(receipt_path, f"{day} run receipt")
    if (
        receipt.get("schema") != RUN_RECEIPT_SCHEMA
        or receipt.get("status") != "SUCCESS"
        or receipt.get("day") != day
    ):
        raise BlindDriverError(f"invalid {day} run receipt status/schema/day")
    anchor = _mapping(receipt.get("freeze"), f"{day}.run_receipt.freeze")
    if anchor.get("sha256") != expected_freeze_sha256.lower() or Path(
        _string(anchor.get("path"), f"{day}.run_receipt.freeze.path")
    ) != Path(_string(plan["freeze_output"], "freeze_output")):
        raise BlindDriverError(f"{day} run receipt freeze SHA-256 mismatch")
    pea_receipt = _mapping(receipt.get("pea"), f"{day}.run_receipt.pea")
    pea_plan = _mapping(plan["pea"], "pea")
    if (
        pea_receipt.get("commit") != pea_plan.get("commit")
        or pea_receipt.get("argv") != pea_plan.get(f"{day}_argv")
        or pea_receipt.get("environment") != pea_plan.get("environment")
        or pea_receipt.get("binary")
        != _frozen_file_record(freeze_receipt, "executable.pea")
    ):
        raise BlindDriverError(f"{day} run receipt PEA identity/argv/env mismatch")
    expected_time = _frozen_file_record(freeze_receipt, "runtime.time")
    if pea_receipt.get("time_binary") != expected_time:
        raise BlindDriverError(f"{day} run receipt time binary mismatch")
    execution = _mapping(receipt.get("execution"), f"{day}.run_receipt.execution")
    working_directory = Path(
        _string(pea_plan["working_directory"], "pea.working_directory")
    ).resolve()
    expected_pea_command = [str(TIME_BINARY), "-v", "--", *pea_plan[f"{day}_argv"]]
    if (
        execution.get("returncode") != 0
        or execution.get("command") != expected_pea_command
        or execution.get("working_directory") != str(working_directory)
    ):
        raise BlindDriverError(
            f"{day} run receipt has wrong PEA command/cwd/return code"
        )
    validators = _mapping(receipt.get("validators"), f"{day}.run_receipt.validators")
    static = _mapping(plan["static_inputs"], "static_inputs")
    expected_validator_commands = {
        "raw": [
            sys.executable,
            _string(static["code.raw_validator"], "code.raw_validator"),
            _string(day_plan["raw_cov"], f"heldout.{day}.raw_cov"),
            "--json-output",
            _string(day_plan["raw_audit"], f"heldout.{day}.raw_audit"),
        ],
        "sd": [
            sys.executable,
            _string(static["code.parser"], "code.parser"),
            _string(day_plan["raw_cov"], f"heldout.{day}.raw_cov"),
            _string(day_plan["sd"], f"heldout.{day}.sd"),
            _string(day_plan["sd_audit"], f"heldout.{day}.sd_audit"),
        ],
    }
    for validator in ("raw", "sd"):
        validation_run = _mapping(
            validators.get(validator), f"{day}.run_receipt.validators.{validator}"
        )
        if (
            validation_run.get("returncode") not in (0, 1)
            or validation_run.get("command") != expected_validator_commands[validator]
            or validation_run.get("working_directory") != str(working_directory)
        ):
            raise BlindDriverError(
                f"{day} run receipt has wrong {validator} validator command/cwd/status"
            )
    artifacts = _mapping(receipt.get("artifacts"), f"{day}.run_receipt.artifacts")
    if set(artifacts) != set(HELDOUT_ARTIFACT_FIELDS):
        raise BlindDriverError(f"{day} run receipt artifact roles are incomplete")
    for field in HELDOUT_ARTIFACT_FIELDS:
        frozen_artifact = _mapping(artifacts[field], f"{day}.artifacts.{field}")
        expected_path = Path(
            _string(day_plan[field], f"heldout.{day}.{field}")
        ).resolve()
        if (
            Path(_string(frozen_artifact.get("path"), f"{day}.{field}.path"))
            != expected_path
        ):
            raise BlindDriverError(f"{day} run receipt {field} path mismatch")
        current = provenance.record_regular_file(expected_path)
        if current != frozen_artifact:
            raise BlindDriverError(f"{day} artifact changed after run: {field}")
    validation = {
        "audits": _validate_day_audits(
            day,
            day_plan,
            f"heldout.{day}",
            _registered_invalid_epoch_keys(plan, day),
        ),
        "console": _successful_console(
            Path(_string(day_plan["console"], f"heldout.{day}.console")),
            f"{day} console",
        ),
        "sd_stage": _posterior_stage_summary(
            Path(_string(day_plan["sd"], f"heldout.{day}.sd")),
            f"{day} SD sidecar",
        ),
    }
    return {
        "receipt": receipt_record,
        "freeze_sha256": expected_freeze_sha256.lower(),
        "pea_sha256": _mapping(pea_receipt["binary"], f"{day}.run_receipt.pea.binary")[
            "sha256"
        ],
        "validation": validation,
    }


def _heldout_inputs(plan: Mapping[str, object]) -> dict[str, str]:
    inputs: dict[str, str] = {}
    heldout = _mapping(plan["heldout"], "heldout")
    for day in DAY_LABELS:
        day_plan = _mapping(heldout[day], f"heldout.{day}")
        root = Path(_string(day_plan["output_root"], f"heldout.{day}.output_root"))
        if not root.is_dir():
            raise BlindDriverError(f"held-out output root is not a directory: {root}")
        for field in HELDOUT_FIELDS:
            path = Path(_string(day_plan[field], f"heldout.{day}.{field}"))
            try:
                path.relative_to(root)
            except ValueError as exc:
                raise BlindDriverError(
                    f"heldout.{day}.{field} escaped frozen output root"
                ) from exc
            inputs[f"heldout.{day}.{field}"] = str(path)
    return inputs


def _audit_integer(audit: Mapping[str, object], field: str, label: str) -> int:
    value = audit.get(field)
    if not isinstance(value, int) or isinstance(value, bool) or value < 0:
        raise BlindDriverError(f"{label}.{field} must be a nonnegative integer")
    return value


def _audit_epoch_key(epoch: Mapping[str, object], label: str) -> tuple[int, Decimal]:
    week = epoch.get("gps_week")
    if not isinstance(week, int) or isinstance(week, bool) or week < 0:
        raise BlindDriverError(f"{label} has invalid gps_week")
    try:
        tow = Decimal(str(epoch.get("gps_tow")))
    except InvalidOperation as exc:
        raise BlindDriverError(f"{label} has invalid gps_tow") from exc
    if not tow.is_finite() or tow < 0 or tow >= Decimal(604800):
        raise BlindDriverError(f"{label} has invalid gps_tow")
    return week, tow


def _audit_epochs(
    audit: Mapping[str, object], label: str
) -> tuple[list[Mapping[str, object]], set[tuple[int, Decimal]]]:
    raw_epochs = audit.get("epochs")
    if not isinstance(raw_epochs, list) or len(raw_epochs) != EXPECTED_EPOCH_COUNT:
        raise BlindDriverError(
            f"{label}.epochs must contain exactly {EXPECTED_EPOCH_COUNT} entries"
        )
    epochs: list[Mapping[str, object]] = []
    keys: set[tuple[int, Decimal]] = set()
    invalid: set[tuple[int, Decimal]] = set()
    for index, value in enumerate(raw_epochs):
        epoch = _mapping(value, f"{label}.epochs[{index}]")
        key = _audit_epoch_key(epoch, f"{label}.epochs[{index}]")
        if key in keys:
            raise BlindDriverError(f"{label} has duplicate epoch {key}")
        keys.add(key)
        if epoch.get("posterior_stage") != EXPECTED_POSTERIOR_STAGE:
            raise BlindDriverError(
                f"{label} epoch {key} has wrong posterior stage "
                f"{epoch.get('posterior_stage')!r}"
            )
        if epoch.get("valid") is not True:
            if epoch.get("valid") is not False:
                raise BlindDriverError(f"{label} epoch {key} has non-boolean valid")
            invalid.add(key)
        epochs.append(epoch)
    return epochs, invalid


def _format_epoch_keys(keys: set[tuple[int, Decimal]]) -> list[str]:
    return [f"{week}:{tow}" for week, tow in sorted(keys)]


def _registered_invalid_epoch_keys(
    plan: Mapping[str, object], day: str
) -> set[tuple[int, Decimal]]:
    policy = _mapping(plan["audit_policy"], "audit_policy")
    day_policy = _mapping(policy[day], f"audit_policy.{day}")
    entries = day_policy.get("expected_invalid_epochs")
    normalised = _normalise_registered_epoch_list(
        entries, f"audit_policy.{day}.expected_invalid_epochs"
    )
    return {
        (int(entry["gps_week"]), Decimal(str(entry["gps_tow"]))) for entry in normalised
    }


def _validate_day_audits(
    day: str,
    files: Mapping[str, object],
    label_prefix: str,
    expected_invalid: set[tuple[int, Decimal]],
) -> dict[str, object]:
    if day not in DAY_LABELS:
        raise BlindDriverError(f"unknown day policy: {day}")
    raw_cov = Path(_string(files["raw_cov"], f"{label_prefix}.raw_cov")).resolve()
    sd_path = Path(_string(files["sd"], f"{label_prefix}.sd")).resolve()
    raw_path = Path(_string(files["raw_audit"], f"{label_prefix}.raw_audit"))
    sd_audit_path = Path(_string(files["sd_audit"], f"{label_prefix}.sd_audit"))
    raw = _load_json_object(raw_path, f"{label_prefix} raw audit")
    sd = _load_json_object(sd_audit_path, f"{label_prefix} SD audit")
    if raw.get("schema") != sd_parser.RAW_SCHEMA:
        raise BlindDriverError(f"{label_prefix} raw audit schema mismatch")
    if (
        sd.get("schema") != sd_parser.REPORT_SCHEMA
        or sd.get("raw_schema") != sd_parser.RAW_SCHEMA
        or sd.get("satellite_difference_schema") != sd_parser.SD_SCHEMA
    ):
        raise BlindDriverError(f"{label_prefix} SD audit schema mismatch")
    if (
        Path(_string(raw.get("input"), f"{label_prefix}.raw.input")).resolve()
        != raw_cov
    ):
        raise BlindDriverError(f"{label_prefix} raw audit names the wrong raw sidecar")
    if (
        Path(_string(sd.get("raw_input"), f"{label_prefix}.sd.raw_input")).resolve()
        != raw_cov
        or Path(
            _string(
                sd.get("satellite_difference_input"),
                f"{label_prefix}.sd.satellite_difference_input",
            )
        ).resolve()
        != sd_path
    ):
        raise BlindDriverError(f"{label_prefix} SD audit names the wrong sidecars")

    raw_epochs, raw_invalid = _audit_epochs(raw, f"{label_prefix}.raw")
    sd_epochs, sd_invalid = _audit_epochs(sd, f"{label_prefix}.sd")
    raw_epoch_keys = [
        _audit_epoch_key(epoch, f"{label_prefix}.raw") for epoch in raw_epochs
    ]
    sd_epoch_keys = [
        _audit_epoch_key(epoch, f"{label_prefix}.sd") for epoch in sd_epochs
    ]
    if raw_epoch_keys != sd_epoch_keys:
        raise BlindDriverError(
            f"{label_prefix} raw/SD complete epoch key order differs"
        )
    start_week, start_tow = EXPECTED_DAY_START_KEYS[day]
    expected_grid = [
        (start_week, start_tow + Decimal(index * 30))
        for index in range(EXPECTED_EPOCH_COUNT)
    ]
    if raw_epoch_keys != expected_grid:
        raise BlindDriverError(
            f"{label_prefix} epochs are not the registered full-day 30-second grid"
        )
    raw_count = _audit_integer(raw, "epoch_count", f"{label_prefix}.raw")
    raw_valid_count = _audit_integer(raw, "valid_epoch_count", f"{label_prefix}.raw")
    raw_invalid_count = _audit_integer(
        raw, "invalid_epoch_count", f"{label_prefix}.raw"
    )
    if (
        raw_count != EXPECTED_EPOCH_COUNT
        or raw_valid_count + raw_invalid_count != raw_count
        or raw_invalid_count != len(raw_invalid)
    ):
        raise BlindDriverError(f"{label_prefix} raw audit counts are inconsistent")
    sd_counts = {
        field: _audit_integer(sd, field, f"{label_prefix}.sd")
        for field in (
            "raw_epoch_count",
            "satellite_difference_epoch_count",
            "matched_epoch_count",
            "valid_epoch_count",
            "invalid_epoch_count",
        )
    }
    if (
        sd_counts["raw_epoch_count"] != EXPECTED_EPOCH_COUNT
        or sd_counts["satellite_difference_epoch_count"] != EXPECTED_EPOCH_COUNT
        or sd_counts["matched_epoch_count"] != EXPECTED_EPOCH_COUNT
        or sd_counts["valid_epoch_count"] + sd_counts["invalid_epoch_count"]
        != EXPECTED_EPOCH_COUNT
        or sd_counts["invalid_epoch_count"] != len(sd_invalid)
        or sd.get("file_errors") != []
    ):
        raise BlindDriverError(f"{label_prefix} SD audit counts/files are inconsistent")
    if raw_invalid != sd_invalid:
        raise BlindDriverError(
            f"{label_prefix} raw/SD invalid epoch keys differ: "
            f"raw={_format_epoch_keys(raw_invalid)}, sd={_format_epoch_keys(sd_invalid)}"
        )

    invalid_count = len(raw_invalid)
    if raw_invalid != expected_invalid:
        raise BlindDriverError(
            f"{label_prefix} invalid epoch keys differ from the registered policy: "
            f"expected={_format_epoch_keys(expected_invalid)}, "
            f"actual={_format_epoch_keys(raw_invalid)}"
        )
    if day == "quiet" and invalid_count != 0:
        raise BlindDriverError(f"{label_prefix} quiet audit must be 2880/2880 valid")
    if day == "storm" and invalid_count > MAX_STORM_NO_STATES:
        raise BlindDriverError(
            f"{label_prefix} storm invalid epoch count {invalid_count} exceeds "
            f"{MAX_STORM_NO_STATES}"
        )
    expected_all_valid = invalid_count == 0
    if (
        raw.get("all_valid") is not expected_all_valid
        or sd.get("all_valid") is not expected_all_valid
    ):
        raise BlindDriverError(f"{label_prefix} all_valid flags disagree with counts")

    raw_by_key = {
        _audit_epoch_key(epoch, f"{label_prefix}.raw"): epoch for epoch in raw_epochs
    }
    sd_by_key = {
        _audit_epoch_key(epoch, f"{label_prefix}.sd"): epoch for epoch in sd_epochs
    }
    for key in raw_invalid:
        raw_epoch = raw_by_key[key]
        sd_epoch = sd_by_key[key]
        raw_errors = raw_epoch.get("errors")
        sd_errors = sd_epoch.get("errors")
        if (
            raw_epoch.get("writer_status") != "NO_STATES"
            or raw_epoch.get("state_count") != 0
            or not isinstance(raw_errors, list)
            or "writer_status=NO_STATES" not in raw_errors
            or sd_epoch.get("source_state_count") != 0
            or sd_epoch.get("difference_state_count") != 0
            or not isinstance(sd_errors, list)
            or "raw_writer_status=NO_STATES" not in sd_errors
            or "sd_writer_status=NO_STATES" not in sd_errors
        ):
            raise BlindDriverError(
                f"{label_prefix} invalid epoch {key} is not registered NO_STATES"
            )
    return {
        "day_policy": day,
        "raw_audit": str(raw_path.resolve()),
        "sd_audit": str(sd_audit_path.resolve()),
        "epoch_count": EXPECTED_EPOCH_COUNT,
        "valid_epoch_count": EXPECTED_EPOCH_COUNT - invalid_count,
        "invalid_epoch_count": invalid_count,
        "invalid_epoch_keys": _format_epoch_keys(raw_invalid),
        "all_valid": expected_all_valid,
    }


def _successful_console(path: Path, label: str) -> dict[str, object]:
    finished = False
    exit_statuses: list[int] = []
    try:
        with path.open("r", encoding="utf-8", errors="strict") as stream:
            for line in stream:
                if "pea finished" in line.lower():
                    finished = True
                exit_statuses.extend(
                    int(match.group(1)) for match in _EXIT_STATUS_PATTERN.finditer(line)
                )
    except (OSError, UnicodeError) as exc:
        raise BlindDriverError(f"cannot scan {label} {path}: {exc}") from exc
    if (
        not finished
        or not exit_statuses
        or any(status != 0 for status in exit_statuses)
    ):
        raise BlindDriverError(
            f"{label} lacks PEA finished plus only Exit status 0: "
            f"finished={finished}, exit_statuses={exit_statuses}"
        )
    return {
        "path": str(path.resolve()),
        "pea_finished": True,
        "exit_statuses": exit_statuses,
    }


def _posterior_stage_summary(path: Path, label: str) -> dict[str, object]:
    stages: set[str] = set()
    epoch_count = 0
    try:
        with path.open("r", encoding="utf-8", errors="strict", newline="") as stream:
            for epoch in sd_parser.iter_difference_epochs(stream):
                epoch_count += 1
                stages.add(epoch.posterior_stage)
    except (OSError, UnicodeError, sd_parser.ValidationInputError, ValueError) as exc:
        raise BlindDriverError(f"cannot scan {label} posterior stage: {exc}") from exc
    if epoch_count == 0 or stages != {EXPECTED_POSTERIOR_STAGE}:
        raise BlindDriverError(
            f"{label} posterior stage mismatch: epochs={epoch_count}, "
            f"stages={sorted(stages)}, expected={EXPECTED_POSTERIOR_STAGE}"
        )
    return {
        "path": str(path.resolve()),
        "epoch_count": epoch_count,
        "posterior_stage": EXPECTED_POSTERIOR_STAGE,
    }


def evaluate_frozen(
    freeze_path: Path,
    freeze_receipt_path: Path,
    expected_freeze_sha256: str,
    expected_quiet_run_receipt_sha256: str,
    expected_storm_run_receipt_sha256: str,
    report_output: Path,
    evaluation_receipt_output: Path,
) -> dict[str, object]:
    """Validate the anchor, evaluate exactly once, and publish report/receipt."""

    freeze, freeze_receipt, plan = _anchored_freeze(
        freeze_path.resolve(), freeze_receipt_path.resolve(), expected_freeze_sha256
    )
    expected_report = Path(_string(plan["report_output"], "report_output"))
    expected_receipt = Path(
        _string(plan["evaluation_receipt_output"], "evaluation_receipt_output")
    )
    expected_claim = Path(
        _string(plan["evaluation_claim_output"], "evaluation_claim_output")
    )
    if report_output.resolve(strict=False) != expected_report.resolve(
        strict=False
    ) or evaluation_receipt_output.resolve(strict=False) != expected_receipt.resolve(
        strict=False
    ):
        raise BlindDriverError(
            "evaluate output paths must exactly equal the frozen report/receipt paths"
        )
    provenance.assert_planned_outputs_absent(
        {
            "report": expected_report,
            "evaluation_receipt": expected_receipt,
            "evaluation_claim": expected_claim,
        }
    )
    report_output = expected_report
    evaluation_receipt_output = expected_receipt
    expected_run_receipt_hashes = {
        "quiet": expected_quiet_run_receipt_sha256,
        "storm": expected_storm_run_receipt_sha256,
    }
    for day, digest in expected_run_receipt_hashes.items():
        if not _SHA256_PATTERN.fullmatch(digest):
            raise BlindDriverError(f"expected {day} run receipt SHA-256 is invalid")
        expected_run_receipt_hashes[day] = digest.lower()
    claim_record = provenance.atomic_write_json_no_replace(
        expected_claim,
        {
            "schema": EVALUATION_CLAIM_SCHEMA,
            "freeze": {
                "path": str(freeze_path.resolve()),
                "sha256": expected_freeze_sha256.lower(),
            },
            "expected_run_receipt_sha256": expected_run_receipt_hashes,
            "claimed_utc": datetime.now(timezone.utc).isoformat(),
        },
    )
    _verify_frozen_static_state(plan, freeze_receipt)
    run_receipts = {
        day: _verify_day_run_receipt(
            day,
            plan,
            freeze_receipt,
            expected_freeze_sha256,
            expected_run_receipt_hashes[day],
        )
        for day in DAY_LABELS
    }
    heldout_inputs = _heldout_inputs(plan)

    frozen_provenance = _mapping(
        freeze_receipt["provenance"], "freeze_receipt.provenance"
    )
    frozen_inputs = _mapping(
        frozen_provenance["inputs_post"], "freeze_receipt.inputs_post"
    )
    evaluation_inputs = {
        f"frozen.{role}": _string(record.get("path"), f"frozen.{role}.path")
        for role, record_value in frozen_inputs.items()
        for record in [_mapping(record_value, f"frozen.{role}")]
    }
    evaluation_inputs.update(heldout_inputs)
    evaluation_inputs["anchor.freeze"] = str(freeze_path.resolve())
    evaluation_inputs["anchor.freeze_receipt"] = str(freeze_receipt_path.resolve())
    evaluation_pre = provenance.begin_provenance(
        _string(plan["repository"], "repository"),
        evaluation_inputs,
        {"report": report_output, "evaluation_receipt": evaluation_receipt_output},
    )

    heldout = _mapping(plan["heldout"], "heldout")
    validation = {
        day: _mapping(run_receipts[day]["validation"], f"run_receipts.{day}")
        for day in DAY_LABELS
    }

    model = _mapping(plan["model"], "model")
    quiet_model = _mapping(model["quiet"], "model.quiet")
    storm_model = _mapping(model["storm"], "model.storm")
    quiet_heldout = _mapping(heldout["quiet"], "heldout.quiet")
    storm_heldout = _mapping(heldout["storm"], "heldout.storm")
    static_inputs = _mapping(plan["static_inputs"], "static_inputs")
    # Scientific evaluation is deliberately a single call after all guards.
    report = spatial.build_report(
        freeze_path,
        Path(_string(static_inputs["audit.product"], "audit.product")),
        Path(_string(quiet_model["stec"], "model.quiet.stec")),
        Path(_string(quiet_model["sd"], "model.quiet.sd")),
        Path(_string(quiet_heldout["stec"], "heldout.quiet.stec")),
        Path(_string(quiet_heldout["sd"], "heldout.quiet.sd")),
        Path(_string(storm_model["stec"], "model.storm.stec")),
        Path(_string(storm_model["sd"], "model.storm.sd")),
        Path(_string(storm_heldout["stec"], "heldout.storm.stec")),
        Path(_string(storm_heldout["sd"], "heldout.storm.sd")),
    )
    completed = provenance.complete_provenance(evaluation_pre, evaluation_inputs)
    report = dict(report)
    report["blind_validation"] = validation
    report["blind_anchor"] = {
        "expected_freeze_sha256": expected_freeze_sha256.lower(),
        "freeze_receipt": str(freeze_receipt_path.resolve()),
        "expected_run_receipt_sha256": expected_run_receipt_hashes,
        "evaluation_claim": claim_record,
    }
    report["blind_provenance"] = completed
    report_record = provenance.atomic_write_json_no_replace(report_output, report)
    receipt = {
        "schema": EVALUATION_RECEIPT_SCHEMA,
        "freeze": {
            "path": str(freeze_path.resolve()),
            "sha256": expected_freeze_sha256.lower(),
        },
        "report": report_record,
        "validation": validation,
        "run_receipts": run_receipts,
        "evaluation_claim": claim_record,
        "provenance": completed,
    }
    receipt_record = provenance.atomic_write_json_no_replace(
        evaluation_receipt_output, receipt
    )
    return {"report": report_record, "evaluation_receipt": receipt_record}


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    subparsers = parser.add_subparsers(dest="command", required=True)
    freeze_parser = subparsers.add_parser(
        "freeze", help="build a blind freeze from one complete plan JSON"
    )
    freeze_parser.add_argument("--plan", type=Path, required=True)

    run_parser = subparsers.add_parser(
        "run", help="execute one quiet/storm held-out command from an anchored freeze"
    )
    run_parser.add_argument("--freeze", type=Path, required=True)
    run_parser.add_argument("--freeze-receipt", type=Path, required=True)
    run_parser.add_argument("--expected-freeze-sha256", required=True)
    run_parser.add_argument("--day", choices=DAY_LABELS, required=True)

    evaluate_parser = subparsers.add_parser(
        "evaluate", help="evaluate exact held-out paths anchored by a freeze SHA-256"
    )
    evaluate_parser.add_argument("--freeze", type=Path, required=True)
    evaluate_parser.add_argument("--freeze-receipt", type=Path, required=True)
    evaluate_parser.add_argument("--expected-freeze-sha256", required=True)
    evaluate_parser.add_argument("--expected-quiet-run-receipt-sha256", required=True)
    evaluate_parser.add_argument("--expected-storm-run-receipt-sha256", required=True)
    evaluate_parser.add_argument("--report-output", type=Path, required=True)
    evaluate_parser.add_argument(
        "--evaluation-receipt-output", type=Path, required=True
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        if args.command == "freeze":
            result = freeze_from_plan(args.plan)
        elif args.command == "run":
            result = run_frozen_day(
                args.freeze,
                args.freeze_receipt,
                args.expected_freeze_sha256,
                args.day,
            )
        else:
            result = evaluate_frozen(
                args.freeze,
                args.freeze_receipt,
                args.expected_freeze_sha256,
                args.expected_quiet_run_receipt_sha256,
                args.expected_storm_run_receipt_sha256,
                args.report_output,
                args.evaluation_receipt_output,
            )
    except (
        BlindDriverError,
        provenance.ProvenanceError,
        spatial.SpatialInputError,
        OSError,
        UnicodeError,
        json.JSONDecodeError,
        yaml.YAMLError,
    ) as exc:
        print(f"blind driver error: {exc}", file=sys.stderr)
        return 2
    print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
