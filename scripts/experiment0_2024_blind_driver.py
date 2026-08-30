#!/usr/bin/env python3
"""Blind freeze/evaluate driver for the Experiment 0 2024 float pivot.

The driver is intentionally separate from the scientific evaluator.  It adds
an externally hashable freeze, immutable input receipts, exact held-out run
plans, and no-replace publication around
``experiment0_2024_spatial_validation.build_freeze`` / ``build_report``.

Freeze plan JSON (paths may be absolute or repository-relative)::

  {
    "schema": "GINAN_EXPERIMENT0_2024_BLIND_PLAN_V1",
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
      "quiet": {"stec": "...", "sd": "...", "audit": "...", "console": "..."},
      "storm": {"stec": "...", "sd": "...", "audit": "...", "console": "..."}
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
        "console": "...console.log"
      },
      "storm": {"...": "same fields"}
    },
    "freeze_output": "/data/freeze.json",
    "freeze_receipt_output": "/data/freeze.receipt.json"
  }

All ``yaml.*`` paths must equal the recursively discovered include closure of
the four named entry points.  Relative include paths use the repository root,
matching these registered Ginan configurations.  The PEA environment is the
exact declared environment map; no implicit inherited variables are claimed.

Security boundary: formal execution assumes ext4 and trusted parent
directories that are not concurrently renamed or replaced.  This module does
not run PEA and does not read held-out outputs during ``freeze``.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
from pathlib import Path
from typing import Mapping, Sequence

import yaml

import experiment0_2024_blind_provenance as provenance
import experiment0_2024_spatial_validation as spatial
import validate_stec_satellite_difference_covariance as sd_parser


PLAN_SCHEMA = "GINAN_EXPERIMENT0_2024_BLIND_PLAN_V1"
FREEZE_RECEIPT_SCHEMA = "GINAN_EXPERIMENT0_2024_BLIND_FREEZE_RECEIPT_V1"
EVALUATION_RECEIPT_SCHEMA = "GINAN_EXPERIMENT0_2024_BLIND_EVALUATION_RECEIPT_V1"
EXPECTED_POSTERIOR_STAGE = "FILTER_POSTERIOR_NO_EPOCH_AR"
DAY_LABELS = ("quiet", "storm")
MODEL_FIELDS = ("stec", "sd", "audit", "console")
HELDOUT_FIELDS = (
    "stec",
    "raw_cov",
    "sd",
    "raw_audit",
    "sd_audit",
    "console",
)
REQUIRED_STATIC_ROLES = {
    "code.driver",
    "code.spatial",
    "code.parser",
    "code.provenance_helper",
    "manifest.input",
    "audit.rinex_qc",
    "audit.split",
    "audit.product",
    "audit.validation_split",
    "yaml.quiet_model",
    "yaml.storm_model",
    "yaml.quiet_heldout",
    "yaml.storm_heldout",
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


def _load_json_object(path: Path, label: str) -> dict[str, object]:
    provenance.record_regular_file(path)
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise BlindDriverError(f"cannot read {label} {path}: {exc}") from exc
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


def _repo_path(value: object, repository: Path, label: str) -> Path:
    return _absolute(_string(value, label), repository)


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
        if yaml_entry.resolve(strict=True) not in yaml_arguments:
            raise BlindDriverError(f"pea.{day}_argv must contain '-y {yaml_entry}'")
        if not any(str(output_root) in argument for argument in argv):
            raise BlindDriverError(
                f"pea.{day}_argv does not contain exact output root {output_root}"
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

    freeze_output = _repo_path(raw.get("freeze_output"), repository, "freeze_output")
    freeze_receipt_output = _repo_path(
        raw.get("freeze_receipt_output"), repository, "freeze_receipt_output"
    )
    if freeze_output == freeze_receipt_output:
        raise BlindDriverError("freeze and freeze receipt paths must differ")

    return {
        "schema": PLAN_SCHEMA,
        "plan_path": str(plan_path.resolve(strict=True)),
        "repository": str(repository),
        "static_inputs": static_inputs,
        "model": model,
        "pea": pea,
        "heldout": heldout,
        "freeze_output": str(freeze_output),
        "freeze_receipt_output": str(freeze_receipt_output),
    }


def _freeze_inputs(plan: Mapping[str, object]) -> dict[str, str]:
    inputs = dict(_mapping(plan["static_inputs"], "static_inputs"))
    inputs["plan.freeze"] = _string(plan["plan_path"], "plan_path")
    pea = _mapping(plan["pea"], "pea")
    inputs["executable.pea"] = _string(pea["binary"], "pea.binary")
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


def _audit_all_valid(path: Path, label: str) -> dict[str, object]:
    audit = _load_json_object(path, label)
    if audit.get("all_valid") is not True:
        raise BlindDriverError(f"{label} does not declare all_valid=true: {path}")
    return {
        "path": str(path.resolve()),
        "schema": audit.get("schema"),
        "all_valid": True,
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
    report_output: Path,
    evaluation_receipt_output: Path,
) -> dict[str, object]:
    """Validate the anchor, evaluate exactly once, and publish report/receipt."""

    provenance.assert_planned_outputs_absent(
        {"report": report_output, "evaluation_receipt": evaluation_receipt_output}
    )
    freeze, freeze_receipt, plan = _anchored_freeze(
        freeze_path.resolve(), freeze_receipt_path.resolve(), expected_freeze_sha256
    )
    _verify_frozen_static_state(plan, freeze_receipt)
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
    validation: dict[str, object] = {}
    for day in DAY_LABELS:
        day_plan = _mapping(heldout[day], f"heldout.{day}")
        validation[day] = {
            "raw_audit": _audit_all_valid(
                Path(_string(day_plan["raw_audit"], f"heldout.{day}.raw_audit")),
                f"{day} raw audit",
            ),
            "sd_audit": _audit_all_valid(
                Path(_string(day_plan["sd_audit"], f"heldout.{day}.sd_audit")),
                f"{day} SD audit",
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

    evaluate_parser = subparsers.add_parser(
        "evaluate", help="evaluate exact held-out paths anchored by a freeze SHA-256"
    )
    evaluate_parser.add_argument("--freeze", type=Path, required=True)
    evaluate_parser.add_argument("--freeze-receipt", type=Path, required=True)
    evaluate_parser.add_argument("--expected-freeze-sha256", required=True)
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
        else:
            result = evaluate_frozen(
                args.freeze,
                args.freeze_receipt,
                args.expected_freeze_sha256,
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
