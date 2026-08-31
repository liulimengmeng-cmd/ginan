#!/usr/bin/env python3
"""Run and audit the corrected 2024 Experiment 1 PPP-AR matrix.

The driver deliberately separates an execution-success receipt from a scientific
claim.  It freezes the exact PEA binary, verifies all 23 registered inputs,
runs five conditions serially, validates covariance sidecars, and publishes
receipts with atomic no-replace semantics.

Examples (from the repository root under WSL)::

    python3 scripts/experiment1_2024_corrected_driver.py run --dry-run
    python3 scripts/experiment1_2024_corrected_driver.py run
    python3 scripts/experiment1_2024_corrected_driver.py audit \
        --batch-root /home/rx/GINAN/inputData/outputs/exp1_2024_arfix_v1_...
    python3 scripts/experiment1_2024_corrected_driver.py audit \
        --batch-root /home/rx/GINAN/inputData/outputs/exp1_2024_arfix_v1_... \
        --write-audit

``--dry-run`` verifies inputs and constructs the complete plan, but creates no
directory and never starts PEA.  By default any live process whose ``comm`` is
``pea`` blocks execution.  The explicit ``--allow-existing-pea-pid`` exception
is intentionally narrow: the complete live PEA set must equal that one PID,
MemAvailable must be at least 3 GiB, SwapFree at least 8 GiB, and new PEA runs
are wrapped in a recorded nice value (19 by default) plus idle-class ionice.
The resource snapshot and the concurrency exception are recorded in every
receipt.

The ``audit`` subcommand is read-only by default and prints its report to
standard output.  It writes a no-replace JSON artifact only when
``--write-audit`` is supplied; ``--output`` may then select a non-default path.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import stat
import subprocess
import sys
import tempfile
import time
from dataclasses import asdict, dataclass
from datetime import datetime, timedelta, timezone
from decimal import Decimal, InvalidOperation
from pathlib import Path
from typing import Iterable, Mapping, Sequence


MANIFEST_SCHEMA = "GINAN_EXPERIMENT1_2024_INPUT_MANIFEST_V1"
PLAN_SCHEMA = "GINAN_EXPERIMENT1_2024_CORRECTED_PLAN_V1"
RUN_RECEIPT_SCHEMA = "GINAN_EXPERIMENT1_2024_CORRECTED_RUN_RECEIPT_V1"
BATCH_RECEIPT_SCHEMA = "GINAN_EXPERIMENT1_2024_CORRECTED_BATCH_RECEIPT_V1"
BATCH_AUDIT_SCHEMA = "GINAN_EXPERIMENT1_2024_CORRECTED_BATCH_AUDIT_V1"
TRACE_AUDIT_SCHEMA = "GINAN_EXPERIMENT1_2024_AR_TRACE_AUDIT_V1"
EXPECTED_INPUT_COUNT = 23
MIN_AVAILABLE_MEMORY_BYTES = 3 * 1024**3
MIN_SWAP_FREE_BYTES = 8 * 1024**3
FORBIDDEN_DYNAMIC_LOADER_ENV = ("LD_PRELOAD", "LD_AUDIT")
EXECUTION_ENV_ALLOWLIST = (
    "HOME",
    "LANG",
    "LC_ALL",
    "LC_CTYPE",
    "LOGNAME",
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "PATH",
    "PYTHONHASHSEED",
    "SHELL",
    "TMPDIR",
    "TZ",
    "USER",
)
GPS_EPOCH = datetime(1980, 1, 6)
HASH_RE = re.compile(r"^[0-9a-f]{64}$")
BATCH_ID_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")
EPOCH_RE = re.compile(r"^fixAndHoldAmbiguities:\s+(.+?)\s*$")
FEEDBACK_RE = re.compile(
    r"PPP_AR INTEGER_FEEDBACK rows=(\d+) integer_coordinates=(\d+) "
    r"original_ambiguities=(\d+) status=([A-Z0-9_]+)"
)
DESIGN_RE = re.compile(
    r"PPP_AR PSEUDOOBS_DESIGN rows=(\d+) original_ambiguities=(\d+) "
    r"(?:z_rows=(\d+) z_columns=(\d+) )?status=([A-Z0-9_]+)"
)
NOISE_RE = re.compile(
    r"PPP_AR PSEUDOOBS_NOISE rows=(\d+)(?: variance=([^\s]+))? "
    r"status=([A-Z0-9_]+)"
)
SUBMISSION_RE = re.compile(
    r"PPP_AR PSEUDOOBS_SUBMISSION rows=(\d+) status=([A-Z0-9_]+)"
)


@dataclass(frozen=True)
class RunSpec:
    order: int
    condition: str
    directory: str
    config: str
    epochs: int
    receiver_sd_expected: bool | None


RUN_SPECS = (
    RunSpec(
        1,
        "float",
        "01_float",
        "Docs/stecCovarianceFeasibility/experiment_1_float_control.yaml",
        480,
        None,
    ),
    RunSpec(
        2,
        "primary",
        "02_primary",
        "Docs/stecCovarianceFeasibility/experiment_1_receiver_sd_ar.yaml",
        480,
        True,
    ),
    RunSpec(
        3,
        "phase_osb_disabled",
        "03_phase_osb_disabled",
        "Docs/stecCovarianceFeasibility/experiment_1_phase_osb_disabled_control.yaml",
        480,
        True,
    ),
    RunSpec(
        4,
        "undifferenced",
        "04_undifferenced",
        "Docs/stecCovarianceFeasibility/experiment_1_undifferenced_ar_control.yaml",
        480,
        False,
    ),
    RunSpec(
        5,
        "restart_0030",
        "05_restart_0030",
        "Docs/stecCovarianceFeasibility/experiment_1_restart_0030.yaml",
        420,
        True,
    ),
)


class DriverError(RuntimeError):
    """A reproducibility or validation invariant failed."""


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _json_bytes(payload: object) -> bytes:
    return (
        json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n"
    ).encode("utf-8")


def _sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def record_file(path: Path | str) -> dict[str, object]:
    """Hash a stable regular file and return its pathname identity."""

    absolute = Path(path).expanduser().absolute()
    try:
        before = absolute.stat()
    except FileNotFoundError as exc:
        raise DriverError(f"required file does not exist: {absolute}") from exc
    if not stat.S_ISREG(before.st_mode):
        raise DriverError(f"required path is not a regular file: {absolute}")
    digest = hashlib.sha256()
    with absolute.open("rb") as stream:
        while True:
            chunk = stream.read(1024 * 1024)
            if not chunk:
                break
            digest.update(chunk)
    after = absolute.stat()
    identity_before = (before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns)
    identity_after = (after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns)
    if identity_before != identity_after:
        raise DriverError(f"file changed while hashing: {absolute}")
    return {
        "path": str(absolute),
        "resolved_path": str(absolute.resolve()),
        "size": int(after.st_size),
        "mtime_ns": int(after.st_mtime_ns),
        "mode": stat.S_IMODE(after.st_mode),
        "dev": int(after.st_dev),
        "inode": int(after.st_ino),
        "sha256": digest.hexdigest(),
    }


def _same_content_record(left: Mapping[str, object], right: Mapping[str, object]) -> bool:
    return left.get("size") == right.get("size") and left.get("sha256") == right.get(
        "sha256"
    )


def build_effective_environment(ld_library_path: str) -> dict[str, str]:
    """Freeze a small, non-secret execution environment and reject loader injection."""

    forbidden = [name for name in FORBIDDEN_DYNAMIC_LOADER_ENV if os.environ.get(name)]
    if forbidden:
        raise DriverError(
            "dynamic-loader injection variables are not permitted: " + ", ".join(forbidden)
        )
    environment = {
        name: os.environ[name] for name in EXECUTION_ENV_ALLOWLIST if name in os.environ
    }
    environment["LD_LIBRARY_PATH"] = ld_library_path
    return dict(sorted(environment.items()))


def verify_effective_environment(environment: Mapping[str, object]) -> None:
    allowed = set(EXECUTION_ENV_ALLOWLIST) | {"LD_LIBRARY_PATH"}
    if not environment or not set(environment).issubset(allowed):
        raise DriverError("execution environment contains an unregistered variable")
    if any(name in environment for name in FORBIDDEN_DYNAMIC_LOADER_ENV):
        raise DriverError("execution environment contains a dynamic-loader injection variable")
    if not isinstance(environment.get("LD_LIBRARY_PATH"), str) or not environment.get(
        "LD_LIBRARY_PATH"
    ):
        raise DriverError("execution environment has no LD_LIBRARY_PATH")
    if not all(isinstance(key, str) and isinstance(value, str) for key, value in environment.items()):
        raise DriverError("execution environment must contain only string pairs")


def atomic_write_bytes_no_replace(path: Path, payload: bytes) -> dict[str, object]:
    """Atomically publish bytes and fail if the destination already exists."""

    destination = path.absolute()
    destination.parent.mkdir(parents=True, exist_ok=True)
    if os.path.lexists(destination):
        raise FileExistsError(f"output already exists: {destination}")
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{destination.name}.", suffix=".tmp", dir=destination.parent
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "wb") as stream:
            stream.write(payload)
            stream.flush()
            os.fsync(stream.fileno())
        try:
            os.link(temporary, destination, follow_symlinks=False)
        except FileExistsError:
            raise
        except OSError as exc:
            raise DriverError(
                f"filesystem cannot atomically publish {destination} without replace: {exc}"
            ) from exc
    finally:
        temporary.unlink(missing_ok=True)
    return record_file(destination)


def atomic_write_json_no_replace(path: Path, payload: object) -> dict[str, object]:
    return atomic_write_bytes_no_replace(path, _json_bytes(payload))


def _load_object(path: Path, label: str) -> dict[str, object]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise DriverError(f"cannot read {label} {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise DriverError(f"{label} must contain one JSON object: {path}")
    return value


def _safe_relative_path(value: object, label: str) -> Path:
    if not isinstance(value, str) or not value:
        raise DriverError(f"{label} must be a non-empty string")
    path = Path(value)
    if path.is_absolute() or ".." in path.parts:
        raise DriverError(f"{label} must be a safe relative path: {value!r}")
    return path


def verify_input_manifest(manifest_path: Path, data_root: Path) -> dict[str, object]:
    manifest = _load_object(manifest_path, "input manifest")
    if manifest.get("schema") != MANIFEST_SCHEMA:
        raise DriverError(f"input manifest schema must be {MANIFEST_SCHEMA}")
    files = manifest.get("files")
    if not isinstance(files, list) or len(files) != EXPECTED_INPUT_COUNT:
        raise DriverError(
            f"input manifest must contain exactly {EXPECTED_INPUT_COUNT} files"
        )
    seen_roles: set[str] = set()
    seen_paths: set[Path] = set()
    records: list[dict[str, object]] = []
    for index, raw in enumerate(files):
        if not isinstance(raw, dict) or set(raw) != {
            "role",
            "relative_path",
            "size",
            "sha256",
        }:
            raise DriverError(f"input manifest files[{index}] has invalid fields")
        role = raw["role"]
        if not isinstance(role, str) or not role or role in seen_roles:
            raise DriverError(f"input manifest files[{index}] has duplicate/invalid role")
        relative = _safe_relative_path(raw["relative_path"], f"files[{index}].relative_path")
        if relative in seen_paths:
            raise DriverError(f"input manifest repeats path {relative}")
        expected_size = raw["size"]
        expected_hash = raw["sha256"]
        if not isinstance(expected_size, int) or expected_size < 0:
            raise DriverError(f"input manifest files[{index}] has invalid size")
        if not isinstance(expected_hash, str) or not HASH_RE.fullmatch(expected_hash):
            raise DriverError(f"input manifest files[{index}] has invalid SHA-256")
        actual = record_file(data_root / relative)
        if actual["size"] != expected_size or actual["sha256"] != expected_hash:
            raise DriverError(
                f"input mismatch for {relative}: expected size={expected_size}, "
                f"sha256={expected_hash}; got size={actual['size']}, "
                f"sha256={actual['sha256']}"
            )
        records.append(
            {
                "role": role,
                "relative_path": relative.as_posix(),
                "expected_size": expected_size,
                "expected_sha256": expected_hash,
                "actual": actual,
            }
        )
        seen_roles.add(role)
        seen_paths.add(relative)
    manifest_record = record_file(manifest_path)
    snapshot_digest = _sha256_bytes(_json_bytes(records))
    return {
        "schema": MANIFEST_SCHEMA,
        "manifest": manifest_record,
        "data_root": str(data_root.absolute()),
        "file_count": len(records),
        "all_match": True,
        "snapshot_sha256": snapshot_digest,
        "files": records,
    }


def active_pea_processes(proc_root: Path = Path("/proc")) -> list[dict[str, object]]:
    processes: list[dict[str, object]] = []
    if not proc_root.is_dir():
        return processes
    for child in proc_root.iterdir():
        if not child.name.isdigit():
            continue
        try:
            comm = (child / "comm").read_text(encoding="utf-8", errors="replace").strip()
        except OSError:
            continue
        if comm != "pea":
            continue
        try:
            stat_text = (child / "stat").read_text(encoding="utf-8", errors="strict")
        except FileNotFoundError:
            continue
        except OSError as exc:
            raise DriverError(f"cannot read process identity for PID {child.name}: {exc}") from exc
        close_parenthesis = stat_text.rfind(")")
        stat_fields = stat_text[close_parenthesis + 1 :].split()
        if close_parenthesis < 0 or len(stat_fields) <= 19:
            raise DriverError(f"invalid /proc stat identity for PID {child.name}")
        try:
            starttime_ticks = int(stat_fields[19])
        except ValueError as exc:
            raise DriverError(f"invalid start time for PID {child.name}") from exc
        command = ""
        try:
            command = (child / "cmdline").read_bytes().replace(b"\0", b" ").decode(
                "utf-8", errors="replace"
            ).strip()
        except OSError:
            pass
        processes.append(
            {
                "pid": int(child.name),
                "starttime_ticks": starttime_ticks,
                "comm": comm,
                "command": command,
            }
        )
    return sorted(processes, key=lambda item: int(item["pid"]))


def resource_snapshot(meminfo_path: Path = Path("/proc/meminfo")) -> dict[str, object]:
    values: dict[str, int] = {}
    try:
        for line in meminfo_path.read_text(encoding="utf-8").splitlines():
            match = re.match(r"^(MemAvailable|SwapFree):\s+(\d+)\s+kB$", line)
            if match:
                values[match.group(1)] = int(match.group(2)) * 1024
    except OSError as exc:
        raise DriverError(f"cannot read resource state {meminfo_path}: {exc}") from exc
    if set(values) != {"MemAvailable", "SwapFree"}:
        raise DriverError(f"resource state is incomplete in {meminfo_path}")
    return {
        "captured_utc": _utc_now(),
        "mem_available_bytes": values["MemAvailable"],
        "swap_free_bytes": values["SwapFree"],
    }


def concurrency_gate(
    allow_existing_pea_pid: int | None,
    concurrent_nice: int = 19,
    *,
    proc_root: Path = Path("/proc"),
    meminfo_path: Path = Path("/proc/meminfo"),
) -> dict[str, object]:
    if concurrent_nice < 0 or concurrent_nice > 19:
        raise DriverError("concurrent nice value must be between 0 and 19")
    active = active_pea_processes(proc_root)
    if not active:
        if allow_existing_pea_pid is not None:
            raise DriverError(
                "--allow-existing-pea-pid names a PID that is not an active pea process"
            )
        return {
            "allowed": True,
            "concurrent_existing_pea": False,
            "active_pea_processes": [],
            "resource_snapshot": resource_snapshot(meminfo_path),
        }
    active_pids = {int(item["pid"]) for item in active}
    if allow_existing_pea_pid is None:
        raise DriverError(f"live PEA process blocks experiment: {sorted(active_pids)}")
    if active_pids != {allow_existing_pea_pid}:
        raise DriverError(
            "--allow-existing-pea-pid requires the complete live PEA set to equal "
            f"{{{allow_existing_pea_pid}}}; found {sorted(active_pids)}"
        )
    resources = resource_snapshot(meminfo_path)
    if int(resources["mem_available_bytes"]) < MIN_AVAILABLE_MEMORY_BYTES:
        raise DriverError("concurrent PEA exception denied: MemAvailable is below 3 GiB")
    if int(resources["swap_free_bytes"]) < MIN_SWAP_FREE_BYTES:
        raise DriverError("concurrent PEA exception denied: SwapFree is below 8 GiB")
    return {
        "allowed": True,
        "concurrent_existing_pea": True,
        "allowed_existing_pea_pid": allow_existing_pea_pid,
        "active_pea_processes": active,
        "resource_snapshot": resources,
        "execution_priority": {
            "nice_adjustment": concurrent_nice,
            "ionice_class_name": "idle",
            "ionice_class_number": 3,
        },
    }


def build_execution_policy(
    concurrency: Mapping[str, object],
    *,
    concurrent_nice: int,
    time_binary: Path,
    nice_binary: Path,
    ionice_binary: Path,
) -> dict[str, object]:
    """Freeze the exact wrapper policy independently of a resource snapshot."""

    wrapper_records: dict[str, dict[str, object]] = {
        "time": record_file(time_binary),
    }
    command_prefix = [str(Path(wrapper_records["time"]["path"])), "-v", "--"]
    policy: dict[str, object] = {
        "mode": "normal",
        "command_prefix": command_prefix,
        "wrapper_records": wrapper_records,
        "priority": None,
        "binary_path_trust_boundary": (
            "unique no-replace frozen path with pre-run and post-run hashing; local single-user "
            "non-adversarial workspace (same-uid malicious replacement is outside this audit)"
        ),
    }
    if not concurrency.get("concurrent_existing_pea"):
        if concurrency.get("execution_priority") is not None:
            raise DriverError("normal preflight unexpectedly contains an execution priority")
        return policy

    expected_priority = {
        "nice_adjustment": concurrent_nice,
        "ionice_class_name": "idle",
        "ionice_class_number": 3,
    }
    if concurrency.get("execution_priority") != expected_priority:
        raise DriverError("preflight execution priority differs from the requested policy")
    wrapper_records["nice"] = record_file(nice_binary)
    wrapper_records["ionice"] = record_file(ionice_binary)
    command_prefix.extend(
        [
            str(Path(wrapper_records["nice"]["path"])),
            "-n",
            str(concurrent_nice),
            str(Path(wrapper_records["ionice"]["path"])),
            "-c",
            "3",
        ]
    )
    policy["mode"] = "concurrent_low_priority"
    policy["priority"] = expected_priority
    return policy


def verify_execution_policy(policy: Mapping[str, object]) -> None:
    wrapper_records = policy.get("wrapper_records")
    command_prefix = policy.get("command_prefix")
    if not isinstance(wrapper_records, dict) or not isinstance(command_prefix, list):
        raise DriverError("execution policy is incomplete")
    if "time" not in wrapper_records:
        raise DriverError("execution policy does not register /usr/bin/time")
    expected_roles = {"time"}
    for role, record in wrapper_records.items():
        if not isinstance(record, dict):
            raise DriverError(f"execution wrapper record is invalid: {role}")
        _verify_record(record, f"execution wrapper {role}")
    expected = [str(Path(str(wrapper_records["time"]["path"]))), "-v", "--"]
    if policy.get("mode") == "concurrent_low_priority":
        priority = policy.get("priority")
        if not isinstance(priority, dict):
            raise DriverError("concurrent execution policy has no priority settings")
        if "nice" not in wrapper_records or "ionice" not in wrapper_records:
            raise DriverError("concurrent execution wrappers are incomplete")
        if "nice_adjustment" not in priority:
            raise DriverError("concurrent execution policy has no nice adjustment")
        nice_adjustment = priority.get("nice_adjustment")
        ionice_class_number = priority.get("ionice_class_number")
        if type(nice_adjustment) is not int or not 0 <= nice_adjustment <= 19:
            raise DriverError("concurrent nice adjustment must be an integer from 0 to 19")
        if priority.get("ionice_class_name") != "idle" or type(
            ionice_class_number
        ) is not int or ionice_class_number != 3:
            raise DriverError("concurrent execution policy has an invalid ionice mapping")
        expected_roles.update({"nice", "ionice"})
        expected.extend(
            [
                str(Path(str(wrapper_records["nice"]["path"]))),
                "-n",
                str(nice_adjustment),
                str(Path(str(wrapper_records["ionice"]["path"]))),
                "-c",
                str(ionice_class_number),
            ]
        )
    elif policy.get("mode") != "normal" or policy.get("priority") is not None:
        raise DriverError("execution policy mode and priority are inconsistent")
    if command_prefix != expected:
        raise DriverError("execution policy command prefix is inconsistent")
    if set(wrapper_records) != expected_roles:
        raise DriverError("execution policy wrapper role set is inconsistent")


def verify_gate_matches_policy(
    gate: Mapping[str, object], policy: Mapping[str, object]
) -> None:
    if gate.get("allowed") is not True:
        raise DriverError("resource gate is not explicitly allowed")
    active = gate.get("active_pea_processes")
    resources = gate.get("resource_snapshot")
    if not isinstance(active, list) or not isinstance(resources, dict):
        raise DriverError("resource gate process or resource snapshot is invalid")
    try:
        active_pids = {
            int(process["pid"])
            for process in active
            if isinstance(process, dict) and "pid" in process
        }
    except (TypeError, ValueError) as exc:
        raise DriverError("resource gate contains an invalid PEA PID") from exc
    if len(active_pids) != len(active):
        raise DriverError("resource gate contains duplicate or malformed PEA processes")
    if any(
        type(process.get("starttime_ticks")) is not int
        or int(process["starttime_ticks"]) <= 0
        for process in active
    ):
        raise DriverError("resource gate contains an invalid PEA process start time")
    concurrent = bool(gate.get("concurrent_existing_pea"))
    if concurrent != (policy.get("mode") == "concurrent_low_priority"):
        raise DriverError("resource gate concurrency mode differs from the frozen plan")
    if concurrent:
        try:
            allowed_pid = int(gate["allowed_existing_pea_pid"])
            mem_available = int(resources["mem_available_bytes"])
            swap_free = int(resources["swap_free_bytes"])
        except (KeyError, TypeError, ValueError) as exc:
            raise DriverError("concurrent resource gate is incomplete") from exc
        if active_pids != {allowed_pid}:
            raise DriverError("concurrent resource gate active PID set is not the allowed singleton")
        if mem_available < MIN_AVAILABLE_MEMORY_BYTES:
            raise DriverError("concurrent resource gate is below the 3 GiB memory threshold")
        if swap_free < MIN_SWAP_FREE_BYTES:
            raise DriverError("concurrent resource gate is below the 8 GiB swap threshold")
        if gate.get("execution_priority") != policy.get("priority"):
            raise DriverError("resource gate priority differs from the frozen plan")
    else:
        if active_pids:
            raise DriverError("normal resource gate contains active PEA processes")
        if gate.get("allowed_existing_pea_pid") is not None:
            raise DriverError("normal resource gate unexpectedly allows an existing PEA PID")
        if gate.get("execution_priority") is not None:
            raise DriverError("normal resource gate unexpectedly contains a priority")


def verify_allowed_process_identity(
    planned_gate: Mapping[str, object], current_gate: Mapping[str, object]
) -> None:
    if not planned_gate.get("concurrent_existing_pea"):
        return
    planned = planned_gate.get("active_pea_processes")
    current = current_gate.get("active_pea_processes")
    if not isinstance(planned, list) or not isinstance(current, list):
        raise DriverError("allowed PEA process identity record is invalid")
    identity_fields = ("pid", "starttime_ticks", "comm", "command")
    planned_identities = [
        tuple(process.get(field) for field in identity_fields)
        for process in planned
        if isinstance(process, dict)
    ]
    current_identities = [
        tuple(process.get(field) for field in identity_fields)
        for process in current
        if isinstance(process, dict)
    ]
    if planned_identities != current_identities:
        raise DriverError("allowed PEA PID was reused or its process identity changed")


def _run_capture(command: Sequence[str], cwd: Path, *, binary: bool = False) -> bytes | str:
    completed = subprocess.run(
        list(command),
        cwd=cwd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    if completed.returncode != 0:
        raise DriverError(
            f"command failed ({completed.returncode}): {' '.join(command)}: "
            f"{completed.stderr.decode(errors='replace').strip()}"
        )
    return completed.stdout if binary else completed.stdout.decode("utf-8").strip()


def git_identity(repository: Path) -> dict[str, object]:
    head = str(_run_capture(["git", "rev-parse", "HEAD"], repository))
    branch = str(_run_capture(["git", "branch", "--show-current"], repository))
    status_bytes = _run_capture(
        ["git", "status", "--porcelain=v1", "-z"], repository, binary=True
    )
    diff_bytes = _run_capture(
        ["git", "diff", "--binary", "HEAD", "--"], repository, binary=True
    )
    assert isinstance(status_bytes, bytes) and isinstance(diff_bytes, bytes)
    return {
        "repository": str(repository.resolve()),
        "head": head,
        "branch": branch,
        "dirty": bool(status_bytes),
        "status_sha256": _sha256_bytes(status_bytes),
        "status_entries": [
            item.decode("utf-8", errors="replace")
            for item in status_bytes.split(b"\0")
            if item
        ],
        "tracked_diff_sha256": _sha256_bytes(diff_bytes),
    }


def _copy_binary_no_replace(source: Path, destination: Path) -> dict[str, object]:
    if os.path.lexists(destination):
        raise FileExistsError(f"frozen binary already exists: {destination}")
    source_record = record_file(source)
    destination.parent.mkdir(parents=True, exist_ok=True)
    with source.open("rb") as input_stream, destination.open("xb") as output_stream:
        while True:
            chunk = input_stream.read(1024 * 1024)
            if not chunk:
                break
            output_stream.write(chunk)
        output_stream.flush()
        os.fsync(output_stream.fileno())
    # Remove write permission after publication.  The plan also records the
    # remaining single-user, non-adversarial pathname trust boundary.
    destination.chmod(0o555)
    frozen_record = record_file(destination)
    if not _same_content_record(source_record, frozen_record):
        raise DriverError("frozen PEA binary differs from the build output")
    return frozen_record


def freeze_planned_binary(
    planned_record: Mapping[str, object], source: Path, destination: Path
) -> dict[str, object]:
    """Freeze PEA only if its bytes still match the already-published plan."""

    immediately_before = record_file(source)
    if not _same_content_record(planned_record, immediately_before):
        raise DriverError("PEA binary changed after planning and before freezing")
    frozen = _copy_binary_no_replace(source, destination)
    if not _same_content_record(planned_record, frozen):
        raise DriverError("frozen PEA binary does not match the planned binary")
    immediately_after = record_file(source)
    if not _same_content_record(planned_record, immediately_after):
        raise DriverError("PEA binary changed while it was being frozen")
    return frozen


def _build_id(path: Path) -> str | None:
    try:
        completed = subprocess.run(
            ["readelf", "-n", str(path)],
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            check=False,
            text=True,
        )
    except OSError:
        return None
    if completed.returncode != 0:
        return None
    match = re.search(r"Build ID:\s*([0-9a-fA-F]+)", completed.stdout)
    return match.group(1).lower() if match else None


def _default_batch_id(git: Mapping[str, object], binary: Mapping[str, object]) -> str:
    timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    return (
        f"exp1_2024_arfix_v1_{timestamp}_head{str(git['head'])[:12]}_"
        f"bin{str(binary['sha256'])[:12]}"
    )


def _config_records(repository: Path) -> dict[str, dict[str, object]]:
    return {
        spec.condition: record_file(repository / spec.config) for spec in RUN_SPECS
    }


def _included_yaml_paths(config: Path, repository: Path) -> list[Path]:
    """Read the narrow ``inputs.include_yamls`` list used by Ginan configs."""

    included: list[Path] = []
    include_indent: int | None = None
    try:
        lines = config.read_text(encoding="utf-8").splitlines()
    except (OSError, UnicodeError) as exc:
        raise DriverError(f"cannot read YAML include closure from {config}: {exc}") from exc
    for raw in lines:
        content = raw.split("#", 1)[0].rstrip()
        if not content.strip():
            continue
        indent = len(content) - len(content.lstrip(" "))
        stripped = content.strip()
        if stripped == "include_yamls:":
            include_indent = indent
            continue
        if include_indent is None:
            continue
        if indent <= include_indent:
            include_indent = None
            continue
        match = re.match(r"^-\s+(.+?)\s*$", stripped)
        if not match:
            continue
        value = match.group(1).strip().strip("'\"")
        path = Path(value)
        if not path.is_absolute():
            path = repository / path
        included.append(path.resolve())
    return included


def config_include_closure(config: Path, repository: Path) -> dict[str, object]:
    """Capture one deterministic recursive include closure, rejecting cycles."""

    ordered: list[dict[str, object]] = []
    visited: set[Path] = set()
    active: set[Path] = set()

    def visit(path: Path) -> None:
        resolved = path.resolve()
        if resolved in active:
            raise DriverError(f"YAML include cycle detected at {resolved}")
        if resolved in visited:
            return
        active.add(resolved)
        record = record_file(resolved)
        ordered.append(record)
        for included in _included_yaml_paths(resolved, repository):
            visit(included)
        active.remove(resolved)
        visited.add(resolved)

    visit(config)
    return {
        "root": str(config.resolve()),
        "file_count": len(ordered),
        "files": ordered,
        "snapshot_sha256": _sha256_bytes(_json_bytes(ordered)),
    }


def verify_config_closure(closure: Mapping[str, object]) -> None:
    files = closure.get("files")
    if not isinstance(files, list) or not files:
        raise DriverError("configuration include closure is empty")
    current: list[dict[str, object]] = []
    for record in files:
        if not isinstance(record, dict):
            raise DriverError("configuration include closure contains an invalid record")
        actual = record_file(Path(str(record.get("path", ""))))
        if not _same_content_record(record, actual):
            raise DriverError(f"configuration include changed: {record.get('path')}")
        current.append(actual)
    # mtime/inode are intentionally excluded from the closure digest comparison.
    planned_content = [
        {"path": item["path"], "size": item["size"], "sha256": item["sha256"]}
        for item in files
    ]
    current_content = [
        {"path": item["path"], "size": item["size"], "sha256": item["sha256"]}
        for item in current
    ]
    if planned_content != current_content:
        raise DriverError("configuration include closure content changed")


def _same_config_closure(
    left: Mapping[str, object], right: Mapping[str, object]
) -> bool:
    """Compare only the immutable path/size/hash identity of two closures."""

    def content(closure: Mapping[str, object]) -> list[tuple[str, int, str]] | None:
        files = closure.get("files")
        if not isinstance(files, list):
            return None
        result: list[tuple[str, int, str]] = []
        for record in files:
            if not isinstance(record, dict):
                return None
            try:
                result.append(
                    (
                        str(record["path"]),
                        int(record["size"]),
                        str(record["sha256"]),
                    )
                )
            except (KeyError, TypeError, ValueError):
                return None
        return result

    return (
        str(left.get("root")) == str(right.get("root"))
        and int(left.get("file_count", -1)) == int(right.get("file_count", -2))
        and content(left) == content(right)
    )


def build_plan(
    *,
    repository: Path,
    data_root: Path,
    output_parent: Path,
    manifest_path: Path,
    pea_binary: Path,
    batch_id: str | None,
    allow_existing_pea_pid: int | None,
    concurrency: Mapping[str, object],
    execution_policy: Mapping[str, object],
    ld_library_path: str,
) -> dict[str, object]:
    repository = repository.resolve()
    data_root = data_root.resolve()
    output_parent = output_parent.resolve()
    if not repository.is_dir() or not output_parent.is_dir():
        raise DriverError("repository and output parent must already exist")
    source = git_identity(repository)
    binary = record_file(pea_binary)
    manifest_snapshot = verify_input_manifest(manifest_path, data_root)
    configs = _config_records(repository)
    chosen_id = batch_id or _default_batch_id(source, binary)
    if not BATCH_ID_RE.fullmatch(chosen_id):
        raise DriverError(f"invalid batch id: {chosen_id!r}")
    batch_root = output_parent / chosen_id
    frozen_binary = batch_root / "provenance" / "pea"
    runs: list[dict[str, object]] = []
    for spec in RUN_SPECS:
        output_root = batch_root / spec.directory
        closure = config_include_closure(repository / spec.config, repository)
        argv = [
            str(frozen_binary),
            "-y",
            str(repository / spec.config),
            "-n",
            str(spec.epochs),
            "-a",
            f"EXP1_DATA_ROOT:{data_root}",
            "-a",
            f"EXP1_OUTPUT_ROOT:{output_root}",
        ]
        runs.append(
            {
                **asdict(spec),
                "config_record": configs[spec.condition],
                "config_include_closure": closure,
                "output_root": str(output_root),
                "argv": argv,
            }
        )
    static_paths = {
        "driver": Path(__file__).resolve(),
        "manifest": manifest_path.resolve(),
        "python_interpreter": Path(sys.executable).resolve(),
        "raw_validator": repository / "scripts/validate_stec_covariance.py",
        "sd_validator": repository
        / "scripts/validate_stec_satellite_difference_covariance.py",
        "restart_comparator": repository
        / "scripts/compare_pppar_integer_constraints.py",
    }
    return {
        "schema": PLAN_SCHEMA,
        "created_utc": _utc_now(),
        "batch_id": chosen_id,
        "batch_root": str(batch_root),
        "repository": str(repository),
        "data_root": str(data_root),
        "output_parent": str(output_parent),
        "source": source,
        "build_binary": binary,
        "frozen_binary": str(frozen_binary),
        "input_snapshot_pre": manifest_snapshot,
        "static_files": {role: record_file(path) for role, path in static_paths.items()},
        "environment": build_effective_environment(ld_library_path),
        "allow_existing_pea_pid": allow_existing_pea_pid,
        "concurrency_preflight": dict(concurrency),
        "execution_policy": dict(execution_policy),
        "runs": runs,
    }


def _execute_logged(
    command: Sequence[str], cwd: Path, environment: Mapping[str, str], log_path: Path, mode: str
) -> dict[str, object]:
    started_utc = _utc_now()
    started = time.monotonic()
    with log_path.open(mode) as stream:
        completed = subprocess.run(
            list(command),
            cwd=cwd,
            env=dict(environment),
            stdout=stream,
            stderr=subprocess.STDOUT,
            check=False,
        )
        stream.flush()
        os.fsync(stream.fileno())
    return {
        "command": list(command),
        "working_directory": str(cwd.resolve()),
        "started_utc": started_utc,
        "ended_utc": _utc_now(),
        "duration_seconds": time.monotonic() - started,
        "returncode": completed.returncode,
    }


def _one(paths: Iterable[Path], label: str) -> Path:
    candidates = sorted({path.resolve() for path in paths})
    if len(candidates) != 1:
        raise DriverError(f"expected one {label}, found {len(candidates)}: {candidates}")
    return candidates[0]


def discover_outputs(output_root: Path) -> dict[str, Path]:
    ionstec = output_root / "ionstec"
    if not ionstec.is_dir():
        raise DriverError(f"missing ionstec directory: {ionstec}")
    traces = list(output_root.glob("Network-*.trace"))
    return {
        "trace_forward": _one(
            (path for path in traces if "_smoothed" not in path.name), "forward trace"
        ),
        "trace_smoothed": _one(
            (path for path in traces if "_smoothed" in path.name), "smoothed trace"
        ),
        "raw_cov_forward": _one(ionstec.glob("*.STEC.COV"), "forward raw covariance"),
        "raw_cov_smoothed": _one(
            ionstec.glob("*.STEC_smoothed.COV"), "smoothed raw covariance"
        ),
        "sd_cov_forward": _one(
            ionstec.glob("*.STEC.SD.COV"), "forward satellite-difference covariance"
        ),
        "sd_cov_smoothed": _one(
            ionstec.glob("*.STEC.SD_smoothed.COV"),
            "smoothed satellite-difference covariance",
        ),
    }


def _int_field(row: Sequence[str], index: int, label: str) -> int:
    try:
        return int(row[index])
    except (IndexError, ValueError) as exc:
        raise DriverError(f"invalid {label} in covariance META row") from exc


def _gps_key_from_calendar(value: str) -> tuple[int, Decimal]:
    text = value.strip()
    try:
        parsed = datetime.strptime(text, "%Y-%m-%d %H:%M:%S.%f")
    except ValueError:
        parsed = datetime.strptime(text, "%Y-%m-%d %H:%M:%S")
    delta = parsed - GPS_EPOCH
    total_seconds = Decimal(delta.days * 86400 + delta.seconds) + Decimal(
        delta.microseconds
    ) / Decimal(1_000_000)
    week = int(total_seconds // Decimal(604800))
    tow = total_seconds - Decimal(week * 604800)
    return week, tow.normalize()


def _covariance_meta(path: Path) -> list[dict[str, object]]:
    meta: list[dict[str, object]] = []
    with path.open("r", encoding="utf-8", newline="") as stream:
        first = stream.readline().strip()
        if first != "# GINAN_STEC_COVARIANCE_V3":
            raise DriverError(f"expected V3 covariance sidecar: {path}")
        reader = csv.reader(stream)
        for row in reader:
            if not row or row[0] != "META":
                continue
            try:
                tow = Decimal(row[2]).normalize()
            except (IndexError, InvalidOperation) as exc:
                raise DriverError(f"invalid GPS time in {path}") from exc
            meta.append(
                {
                    "key": (_int_field(row, 1, "gps_week"), tow),
                    "writer_status": row[3],
                    "ar_invoked": _int_field(row, 8, "ar_routine_invoked"),
                    "eligible": _int_field(row, 9, "eligible ambiguities"),
                    "integer_coordinates": _int_field(row, 10, "integer coordinates"),
                    "receiver_sd": _int_field(row, 11, "receiver SD flag"),
                    "datum_groups": _int_field(row, 12, "datum groups"),
                    "dropped_singletons": _int_field(row, 13, "dropped singletons"),
                    "resolved": _int_field(row, 14, "resolved combinations"),
                    "pseudo_submitted": _int_field(row, 15, "pseudo submitted flag"),
                    "mode": row[16],
                    "diagnostic": row[19],
                    "selected": _int_field(row, 20, "selected ambiguities"),
                    "candidate_count": _int_field(row, 21, "candidate count"),
                }
            )
    return meta


def _trace_feedback_events(
    path: Path,
) -> dict[tuple[int, Decimal], dict[str, list[dict[str, object]]]]:
    current: tuple[int, Decimal] | None = None
    events: dict[tuple[int, Decimal], dict[str, list[dict[str, object]]]] = {}
    with path.open("r", encoding="utf-8", errors="replace") as stream:
        for raw in stream:
            epoch = EPOCH_RE.match(raw.strip())
            if epoch:
                current = _gps_key_from_calendar(epoch.group(1))
            feedback_matches = list(FEEDBACK_RE.finditer(raw))
            design_matches = list(DESIGN_RE.finditer(raw))
            noise_matches = list(NOISE_RE.finditer(raw))
            submission_matches = list(SUBMISSION_RE.finditer(raw))
            if (
                feedback_matches
                or design_matches
                or noise_matches
                or submission_matches
            ) and current is None:
                raise DriverError("integer feedback evidence precedes AR epoch marker")
            assert current is not None or not (
                feedback_matches or design_matches or noise_matches or submission_matches
            )
            if current is None:
                continue
            bucket = events.setdefault(
                current,
                {"feedback": [], "design": [], "noise": [], "submission": []},
            )
            for match in feedback_matches:
                bucket["feedback"].append(
                    {
                        "rows": int(match.group(1)),
                        "integer_coordinates": int(match.group(2)),
                        "original_ambiguities": int(match.group(3)),
                        "status": match.group(4),
                    }
                )
            for match in design_matches:
                bucket["design"].append(
                    {
                        "rows": int(match.group(1)),
                        "original_ambiguities": int(match.group(2)),
                        "z_rows": int(match.group(3)) if match.group(3) else None,
                        "z_columns": int(match.group(4)) if match.group(4) else None,
                        "status": match.group(5),
                    }
                )
            for match in noise_matches:
                bucket["noise"].append(
                    {
                        "rows": int(match.group(1)),
                        "variance": match.group(2),
                        "status": match.group(3),
                    }
                )
            for match in submission_matches:
                bucket["submission"].append(
                    {"rows": int(match.group(1)), "status": match.group(2)}
                )
    return events


def audit_trace(
    *, condition: str, expected_epochs: int, receiver_sd_expected: bool | None,
    raw_covariance: Path, trace: Path,
) -> dict[str, object]:
    meta = _covariance_meta(raw_covariance)
    feedback = _trace_feedback_events(trace)
    accounting_violations: list[str] = []
    mode_violations: list[str] = []
    candidate_violations: list[str] = []
    noise_violations: list[str] = []
    resolved_epochs = 0
    full_rank_epochs = 0
    max_resolved = 0
    max_candidates = 0
    meta_keys = {entry["key"] for entry in meta}
    for entry in meta:
        key = entry["key"]
        invoked = int(entry["ar_invoked"])
        eligible = int(entry["eligible"])
        integer_coordinates = int(entry["integer_coordinates"])
        datum_groups = int(entry["datum_groups"])
        dropped = int(entry["dropped_singletons"])
        resolved = int(entry["resolved"])
        max_resolved = max(max_resolved, resolved)
        pseudo = int(entry["pseudo_submitted"])
        candidates = int(entry["candidate_count"])
        max_candidates = max(max_candidates, candidates)
        if condition == "float":
            if invoked != 0 or resolved != 0 or pseudo != 0 or entry["mode"] != "OFF":
                mode_violations.append(f"{key}: FLOAT invoked/submitted AR")
        else:
            if invoked != 1:
                mode_violations.append(f"{key}: AR routine was not invoked")
            if receiver_sd_expected is True:
                if int(entry["receiver_sd"]) != 1:
                    mode_violations.append(f"{key}: receiver SD was not applied")
                if integer_coordinates + datum_groups + dropped != eligible:
                    accounting_violations.append(
                        f"{key}: {integer_coordinates}+{datum_groups}+{dropped}!={eligible}"
                    )
            elif receiver_sd_expected is False and int(entry["receiver_sd"]) != 0:
                mode_violations.append(f"{key}: undifferenced control applied receiver SD")
        if candidates > 2:
            candidate_violations.append(f"{key}: candidate_count={candidates}")
        events = feedback.get(
            key, {"feedback": [], "design": [], "noise": [], "submission": []}
        )
        if resolved > 0:
            resolved_epochs += 1
            if resolved == integer_coordinates:
                full_rank_epochs += 1
            if pseudo != 1:
                noise_violations.append(f"{key}: resolved={resolved} but pseudo flag={pseudo}")
            expected_statuses = {
                "feedback": "Z_TIMES_D_MAPPED",
                "design": "MATCHES_Z_TIMES_D",
                "noise": "INDEPENDENT_DIAGONAL",
                "submission": "FILTER_CALL_RETURNED_SUBMITTED_UNVERIFIED",
            }
            for evidence_type, expected_status in expected_statuses.items():
                evidence = events[evidence_type]
                if len(evidence) != 1:
                    noise_violations.append(
                        f"{key}: expected one {evidence_type} event, found {len(evidence)}"
                    )
                    continue
                event = evidence[0]
                if event["rows"] != resolved or event["status"] != expected_status:
                    noise_violations.append(
                        f"{key}: resolved={resolved}, {evidence_type}={event}"
                    )
                if evidence_type in {"feedback", "design"} and int(
                    event["original_ambiguities"]
                ) != eligible:
                    noise_violations.append(
                        f"{key}: {evidence_type} original ambiguity count differs from META"
                    )
                if evidence_type == "feedback" and int(
                    event["integer_coordinates"]
                ) != integer_coordinates:
                    noise_violations.append(
                        f"{key}: feedback integer-coordinate count differs from META"
                    )
        elif any(events.values()):
            noise_violations.append(
                f"{key}: integer-feedback evidence exists without resolved rows"
            )
    for key in set(feedback) - meta_keys:
        if any(feedback[key].values()):
            noise_violations.append(
                f"{key}: trace integer-feedback evidence has no covariance META row"
            )
    if condition == "undifferenced" and full_rank_epochs:
        mode_violations.append(
            f"undifferenced negative control produced {full_rank_epochs} complete-rank epochs"
        )
    checks = {
        "expected_epoch_count": len(meta) == expected_epochs,
        "receiver_datum_accounting": not accounting_violations,
        "condition_mode": not mode_violations,
        "lambda_candidate_limit": not candidate_violations,
        "integer_feedback_design_noise_submission": not noise_violations,
        "undifferentiated_no_complete_rank": not (
            condition == "undifferenced" and full_rank_epochs > 0
        ),
    }
    coverage_checks = {
        "multirow_runtime_coverage": condition != "primary" or max_resolved >= 2,
    }
    return {
        "schema": TRACE_AUDIT_SCHEMA,
        "condition": condition,
        "expected_epoch_count": expected_epochs,
        "epoch_count": len(meta),
        "resolved_epoch_count": resolved_epochs,
        "complete_rank_epoch_count": full_rank_epochs,
        "maximum_resolved_combination_count": max_resolved,
        "maximum_integer_candidate_count": max_candidates,
        "feedback_evidence_epoch_count": sum(
            1 for value in feedback.values() if value["feedback"]
        ),
        "design_evidence_epoch_count": sum(
            1 for value in feedback.values() if value["design"]
        ),
        "noise_evidence_epoch_count": sum(
            1 for value in feedback.values() if value["noise"]
        ),
        "submission_evidence_epoch_count": sum(
            1 for value in feedback.values() if value["submission"]
        ),
        "checks": checks,
        "coverage_checks": coverage_checks,
        "multirow_runtime_coverage": coverage_checks["multirow_runtime_coverage"],
        "all_valid": all(checks.values()),
        "violations": {
            "accounting": accounting_violations[:20],
            "mode": mode_violations[:20],
            "candidate": candidate_violations[:20],
            "noise": noise_violations[:20],
        },
        "claim_limit": (
            "This audit proves runtime exposure of the candidate-pool, mapped-design, "
            "independent-noise, and filter-call-returned submission invariants only. It does "
            "not establish filter acceptance or correctness of the submitted integer feedback."
        ),
    }


def _validator_commands(
    repository: Path, outputs: Mapping[str, Path], audit_dir: Path
) -> list[tuple[str, list[str], Path]]:
    raw_validator = repository / "scripts/validate_stec_covariance.py"
    sd_validator = repository / "scripts/validate_stec_satellite_difference_covariance.py"
    return [
        (
            "forward_raw",
            [
                sys.executable,
                str(raw_validator),
                str(outputs["raw_cov_forward"]),
                "--json-output",
                str(audit_dir / "forward.raw.json"),
            ],
            audit_dir / "forward.raw.json",
        ),
        (
            "forward_sd",
            [
                sys.executable,
                str(sd_validator),
                str(outputs["raw_cov_forward"]),
                str(outputs["sd_cov_forward"]),
                str(audit_dir / "forward.sd.json"),
            ],
            audit_dir / "forward.sd.json",
        ),
        (
            "smoothed_raw",
            [
                sys.executable,
                str(raw_validator),
                str(outputs["raw_cov_smoothed"]),
                "--json-output",
                str(audit_dir / "smoothed.raw.json"),
            ],
            audit_dir / "smoothed.raw.json",
        ),
        (
            "smoothed_sd",
            [
                sys.executable,
                str(sd_validator),
                str(outputs["raw_cov_smoothed"]),
                str(outputs["sd_cov_smoothed"]),
                str(audit_dir / "smoothed.sd.json"),
            ],
            audit_dir / "smoothed.sd.json",
        ),
    ]


def _available_files(root: Path) -> dict[str, object]:
    result: dict[str, object] = {}
    for path in sorted(root.rglob("*")):
        if not path.is_file() or path.name == "run.receipt.json":
            continue
        try:
            result[path.relative_to(root).as_posix()] = record_file(path)
        except Exception as exc:  # preserve failure evidence without hiding original error
            result[path.relative_to(root).as_posix()] = {
                "path": str(path),
                "error": str(exc),
            }
    return result


def _publish_run_receipt(output_root: Path, payload: dict[str, object]) -> dict[str, object]:
    return atomic_write_json_no_replace(output_root / "run.receipt.json", payload)


def execute_condition(
    *,
    plan: Mapping[str, object],
    run: Mapping[str, object],
    frozen_binary: Path,
    gate: Mapping[str, object],
) -> dict[str, object]:
    repository = Path(str(plan["repository"]))
    output_root = Path(str(run["output_root"]))
    try:
        output_root.mkdir()
    except FileExistsError as exc:
        raise DriverError(f"run output directory already exists: {output_root}") from exc
    audit_dir = output_root / "audit"
    audit_dir.mkdir()
    console = output_root / "run.console.log"
    exit_code_path = output_root / "run.exit_code"
    verify_effective_environment(dict(plan["environment"]))
    environment = {str(k): str(v) for k, v in dict(plan["environment"]).items()}
    pea_argv = list(run["argv"])
    pea_argv[0] = str(frozen_binary)
    policy = plan.get("execution_policy")
    if not isinstance(policy, dict):
        raise DriverError("frozen execution policy is missing from the plan")
    verify_execution_policy(policy)
    verify_gate_matches_policy(gate, policy)
    command = [*list(policy["command_prefix"]), *pea_argv]
    execution: dict[str, object] = {
        "command": command,
        "working_directory": str(repository),
        "returncode": None,
    }
    validators: dict[str, object] = {}
    validation: dict[str, object] = {}
    status = "DRIVER_FAILED"
    failure: dict[str, object] | None = None
    try:
        execution = _execute_logged(command, repository, environment, console, "xb")
        atomic_write_bytes_no_replace(
            exit_code_path, f"{execution['returncode']}\n".encode("ascii")
        )
        if execution["returncode"] != 0:
            status = "PEA_FAILED"
            failure = {"stage": "PEA", "error": f"return code {execution['returncode']}"}
        else:
            outputs = discover_outputs(output_root)
            validator_log = audit_dir / "validators.console.log"
            planned_python = dict(plan["static_files"])["python_interpreter"]
            if not isinstance(planned_python, dict):
                raise DriverError("planned Python interpreter record is invalid")
            validator_specs = _validator_commands(repository, outputs, audit_dir)
            for index, (name, validator_command, report_path) in enumerate(validator_specs):
                _verify_record(planned_python, "planned Python interpreter")
                validator_command[0] = str(planned_python["path"])
                tool_role = "raw_validator" if name.endswith("_raw") else "sd_validator"
                planned_tool = dict(plan["static_files"])[tool_role]
                if not isinstance(planned_tool, dict):
                    raise DriverError(f"planned validator record is invalid: {tool_role}")
                _verify_record(planned_tool, f"planned validator {tool_role}")
                if Path(validator_command[1]).resolve() != Path(
                    str(planned_tool["path"])
                ).resolve():
                    raise DriverError(f"validator command path differs from plan: {name}")
                if os.path.lexists(report_path):
                    raise DriverError(f"validator output already exists: {report_path}")
                result = _execute_logged(
                    validator_command,
                    repository,
                    environment,
                    validator_log,
                    "xb" if index == 0 else "ab",
                )
                _verify_record(planned_tool, f"post-run validator {tool_role}")
                _verify_record(planned_python, "post-run Python interpreter")
                result["tool_role"] = tool_role
                result["tool"] = record_file(Path(str(planned_tool["path"])))
                result["python_interpreter"] = record_file(
                    Path(str(planned_python["path"]))
                )
                validators[name] = result
            validator_reports = {
                name: _load_object(report_path, f"{name} validator report")
                for name, _, report_path in validator_specs
            }
            trace_audit = audit_trace(
                condition=str(run["condition"]),
                expected_epochs=int(run["epochs"]),
                receiver_sd_expected=run["receiver_sd_expected"],
                raw_covariance=outputs["raw_cov_forward"],
                trace=outputs["trace_forward"],
            )
            trace_audit_path = audit_dir / "ar_trace_audit.json"
            atomic_write_json_no_replace(trace_audit_path, trace_audit)
            validator_valid = all(
                result["returncode"] == 0 for result in validators.values()
            ) and all(bool(report.get("all_valid")) for report in validator_reports.values())
            validation = {
                "validator_reports": validator_reports,
                "trace_audit": trace_audit,
                "all_valid": validator_valid and bool(trace_audit["all_valid"]),
            }
            status = "SUCCESS" if validation["all_valid"] else "VALIDATION_FAILED"
            if status != "SUCCESS":
                failure = {"stage": "POSTPROCESSING", "error": "one or more audits failed"}
    except Exception as exc:
        status = "DRIVER_FAILED"
        failure = {
            "stage": "EXECUTION_OR_POSTPROCESSING",
            "error_type": type(exc).__name__,
            "error": str(exc),
        }
    payload: dict[str, object] = {
        "schema": RUN_RECEIPT_SCHEMA,
        "status": status,
        "batch_id": plan["batch_id"],
        "condition": run["condition"],
        "order": run["order"],
        "expected_epochs": run["epochs"],
        "config": run["config_record"],
        "config_include_closure": run["config_include_closure"],
        "aliases": {
            "EXP1_DATA_ROOT": plan["data_root"],
            "EXP1_OUTPUT_ROOT": str(output_root),
        },
        "source": plan["source"],
        "input_manifest": {
            "manifest": dict(plan["input_snapshot_pre"])["manifest"],
            "snapshot_sha256": dict(plan["input_snapshot_pre"])["snapshot_sha256"],
            "file_count": EXPECTED_INPUT_COUNT,
        },
        "pea": {
            "binary": record_file(frozen_binary),
            "elf_build_id": _build_id(frozen_binary),
            "argv": pea_argv,
            "environment": plan["environment"],
        },
        "planned_execution_policy": policy,
        "concurrency": dict(gate),
        "execution": execution,
        "validators": validators,
        "validation": validation,
        "artifacts": _available_files(output_root),
    }
    if failure is not None:
        payload["failure"] = failure
    receipt_record = _publish_run_receipt(output_root, payload)
    return {"status": status, "receipt": receipt_record, "payload": payload}


def _restart_comparison(
    plan: Mapping[str, object], run_results: Sequence[Mapping[str, object]], batch_root: Path
) -> dict[str, object]:
    by_condition = {
        str(result["payload"]["condition"]): result for result in run_results  # type: ignore[index]
    }
    primary_root = Path(str(next(run for run in plan["runs"] if run["condition"] == "primary")["output_root"]))  # type: ignore[index]
    restart_root = Path(str(next(run for run in plan["runs"] if run["condition"] == "restart_0030")["output_root"]))  # type: ignore[index]
    if by_condition["primary"]["status"] != "SUCCESS" or by_condition["restart_0030"]["status"] != "SUCCESS":
        return {"status": "SKIPPED_FAILED_RUN"}
    primary_trace = discover_outputs(primary_root)["trace_forward"]
    restart_trace = discover_outputs(restart_root)["trace_forward"]
    report_path = batch_root / "audit" / "restart_integer_comparison.json"
    log_path = batch_root / "audit" / "restart_integer_comparison.console.log"
    planned_comparator = dict(plan["static_files"])["restart_comparator"]
    planned_python = dict(plan["static_files"])["python_interpreter"]
    if not isinstance(planned_comparator, dict):
        raise DriverError("planned restart comparator record is invalid")
    if not isinstance(planned_python, dict):
        raise DriverError("planned Python interpreter record is invalid")
    _verify_record(planned_comparator, "planned restart comparator")
    _verify_record(planned_python, "planned Python interpreter")
    comparator = Path(str(planned_comparator["path"]))
    command = [
        str(planned_python["path"]),
        str(comparator),
        str(primary_trace),
        str(restart_trace),
        "--json-output",
        str(report_path),
    ]
    verify_effective_environment(dict(plan["environment"]))
    environment = {str(k): str(v) for k, v in dict(plan["environment"]).items()}
    execution = _execute_logged(command, Path(str(plan["repository"])), environment, log_path, "xb")
    _verify_record(planned_comparator, "post-run restart comparator")
    _verify_record(planned_python, "post-run Python interpreter")
    comparator_record = record_file(comparator)
    python_record = record_file(Path(str(planned_python["path"])))
    report = _load_object(report_path, "restart comparison")
    shared = int(report.get("shared_implied_pair_comparisons", 0))
    consistent = bool(report.get("integer_differences_consistent_on_shared_pairs"))
    return {
        "status": (
            "CONSISTENT_SHARED_FIXES"
            if execution["returncode"] == 0 and shared > 0 and consistent
            else "INCONCLUSIVE_NO_SHARED_FIXES"
            if execution["returncode"] == 0 and shared == 0
            else "INCONSISTENT_OR_FAILED"
        ),
        "execution": execution,
        "tool": comparator_record,
        "python_interpreter": python_record,
        "report": report,
        "artifacts": {
            "report": record_file(report_path),
            "console": record_file(log_path),
        },
    }


def _same_snapshot(left: Mapping[str, object], right: Mapping[str, object]) -> bool:
    return left.get("snapshot_sha256") == right.get("snapshot_sha256") and left.get(
        "manifest", {}
    ).get("sha256") == right.get("manifest", {}).get("sha256")  # type: ignore[union-attr]


def _batch_scientific_assessment(
    run_results: Sequence[Mapping[str, object]], comparison: Mapping[str, object]
) -> dict[str, object]:
    summaries: dict[str, Mapping[str, object]] = {}
    for result in run_results:
        payload = result.get("payload")
        if not isinstance(payload, dict):
            continue
        validation = payload.get("validation")
        trace = validation.get("trace_audit") if isinstance(validation, dict) else None
        if isinstance(trace, dict):
            summaries[str(payload.get("condition"))] = trace
    primary = int(summaries.get("primary", {}).get("resolved_epoch_count", 0))
    primary_max_resolved = int(
        summaries.get("primary", {}).get("maximum_resolved_combination_count", 0)
    )
    phase_off = int(
        summaries.get("phase_osb_disabled", {}).get("resolved_epoch_count", 0)
    )
    undifferenced_full = int(
        summaries.get("undifferenced", {}).get("complete_rank_epoch_count", 0)
    )
    flags: list[str] = []
    if primary == 0:
        flags.append("NO_PRIMARY_INTEGER_FEEDBACK")
    if primary_max_resolved < 2:
        flags.append("PRIMARY_MULTIROW_RUNTIME_COVERAGE_NOT_ESTABLISHED")
    if phase_off > 0:
        flags.append("PHASE_OSB_DISABLED_SUBMITTED_FEEDBACK_ROWS")
    if undifferenced_full > 0:
        flags.append("UNDIFFERENCED_COMPLETE_RANK_STOP_RULE_VIOLATION")
    if comparison.get("status") != "CONSISTENT_SHARED_FIXES":
        flags.append("RESTART_SHARED_FIXES_NOT_ESTABLISHED")
    return {
        "status": "DEFECT_FIX_RUNTIME_EVIDENCE_COMPLETE" if not flags else "INCONCLUSIVE",
        "primary_resolved_epoch_count": primary,
        "primary_maximum_resolved_combination_count": primary_max_resolved,
        "primary_multirow_runtime_coverage": primary_max_resolved >= 2,
        "phase_osb_disabled_resolved_epoch_count": phase_off,
        "phase_osb_disabled_to_primary_ratio": (
            phase_off / primary if primary > 0 else None
        ),
        "undifferenced_complete_rank_epoch_count": undifferenced_full,
        "flags": flags,
        "required_followup": (
            "Run a separately labelled non-scientific multirow pipeline probe; do not treat "
            "single-row formal-data feedback as coverage of the covariance-shape defect."
            if primary_max_resolved < 2
            else None
        ),
        "claim_limit": (
            "Any phase-OSB-disabled submitted feedback keeps OSB/sign sensitivity unresolved; "
            "the driver does not silently promote an execution-success batch."
        ),
    }


def run_batch(args: argparse.Namespace) -> dict[str, object]:
    repository = args.repository.resolve()
    data_root = args.data_root.resolve()
    output_parent = args.output_parent.resolve()
    manifest = args.input_manifest.resolve()
    pea_binary = args.pea_binary.resolve()
    preflight = concurrency_gate(args.allow_existing_pea_pid, args.concurrent_nice)
    execution_policy = build_execution_policy(
        preflight,
        concurrent_nice=args.concurrent_nice,
        time_binary=args.time_binary,
        nice_binary=args.nice_binary,
        ionice_binary=args.ionice_binary,
    )
    verify_execution_policy(execution_policy)
    verify_gate_matches_policy(preflight, execution_policy)
    plan = build_plan(
        repository=repository,
        data_root=data_root,
        output_parent=output_parent,
        manifest_path=manifest,
        pea_binary=pea_binary,
        batch_id=args.batch_id,
        allow_existing_pea_pid=args.allow_existing_pea_pid,
        concurrency=preflight,
        execution_policy=execution_policy,
        ld_library_path=args.ld_library_path,
    )
    verify_effective_environment(dict(plan["environment"]))
    if plan["allow_existing_pea_pid"] != preflight.get("allowed_existing_pea_pid"):
        raise DriverError("requested allowed PEA PID differs from the preflight gate")
    if args.dry_run:
        return {
            "status": "DRY_RUN",
            "pea_started": False,
            "batch_directory_created": False,
            "plan": plan,
        }
    batch_root = Path(str(plan["batch_root"]))
    try:
        batch_root.mkdir()
    except FileExistsError as exc:
        raise DriverError(f"batch directory already exists: {batch_root}") from exc
    (batch_root / "provenance").mkdir()
    (batch_root / "audit").mkdir()
    plan_record = atomic_write_json_no_replace(batch_root / "batch.plan.json", plan)
    frozen_binary = Path(str(plan["frozen_binary"]))
    frozen_record = freeze_planned_binary(plan["build_binary"], pea_binary, frozen_binary)
    frozen_record["elf_build_id"] = _build_id(frozen_binary)
    run_results: list[dict[str, object]] = []
    batch_status = "FAILED"
    comparison: dict[str, object] = {"status": "NOT_RUN"}
    failure: dict[str, object] | None = None
    try:
        for run in plan["runs"]:  # type: ignore[assignment]
            gate = concurrency_gate(args.allow_existing_pea_pid, args.concurrent_nice)
            verify_execution_policy(execution_policy)
            verify_gate_matches_policy(gate, execution_policy)
            verify_allowed_process_identity(preflight, gate)
            verify_config_closure(run["config_include_closure"])
            if not _same_content_record(frozen_record, record_file(frozen_binary)):
                raise DriverError("frozen binary changed during the batch")
            result = execute_condition(
                plan=plan,
                run=run,
                frozen_binary=frozen_binary,
                gate=gate,
            )
            run_results.append(result)
            if result["status"] != "SUCCESS":
                failure = {
                    "stage": str(run["condition"]),
                    "error": f"run status {result['status']}",
                }
                break
        if len(run_results) == len(RUN_SPECS) and all(
            result["status"] == "SUCCESS" for result in run_results
        ):
            comparison = _restart_comparison(plan, run_results, batch_root)
            batch_status = (
                "SUCCESS"
                if comparison["status"] != "INCONSISTENT_OR_FAILED"
                else "AUDIT_FAILED"
            )
        input_post = verify_input_manifest(manifest, data_root)
        if not _same_snapshot(plan["input_snapshot_pre"], input_post):
            raise DriverError("input snapshot changed during the batch")
        if not _same_content_record(frozen_record, record_file(frozen_binary)):
            raise DriverError("frozen binary changed during the batch")
    except Exception as exc:
        batch_status = "FAILED"
        failure = {
            "stage": "BATCH_DRIVER",
            "error_type": type(exc).__name__,
            "error": str(exc),
        }
        input_post = {"status": "NOT_VERIFIED", "error": str(exc)}
    scientific_assessment = _batch_scientific_assessment(run_results, comparison)
    receipt_payload: dict[str, object] = {
        "schema": BATCH_RECEIPT_SCHEMA,
        "status": batch_status,
        "batch_id": plan["batch_id"],
        "batch_root": str(batch_root),
        "plan": plan_record,
        "source": plan["source"],
        "binary": frozen_record,
        "concurrency_preflight": preflight,
        "execution_policy": execution_policy,
        "input_snapshot_pre": plan["input_snapshot_pre"],
        "input_snapshot_post": input_post,
        "runs": [
            {
                "condition": result["payload"]["condition"],
                "status": result["status"],
                "receipt": result["receipt"],
            }
            for result in run_results
        ],
        "restart_comparison": comparison,
        "scientific_assessment": scientific_assessment,
        "claim_limit": (
            "This five-run batch tests the corrected ambiguity machinery. It does not include "
            "inverted/double phase-OSB or deliberate wrong-fix controls and therefore cannot "
            "establish fixed-STEC or PPP-AR scientific acceptance."
        ),
    }
    if failure is not None:
        receipt_payload["failure"] = failure
    batch_receipt = atomic_write_json_no_replace(
        batch_root / "batch.receipt.json", receipt_payload
    )
    result = {
        "status": batch_status,
        "batch_root": str(batch_root),
        "batch_receipt": batch_receipt,
        "run_count": len(run_results),
    }
    if batch_status != "SUCCESS":
        raise DriverError(json.dumps(result, sort_keys=True))
    return result


def _verify_record(record: Mapping[str, object], label: str) -> None:
    path = Path(str(record.get("path", "")))
    current = record_file(path)
    if not _same_content_record(record, current):
        raise DriverError(f"{label} changed: {path}")
    if "mode" in record and record.get("mode") != current.get("mode"):
        raise DriverError(f"{label} mode changed: {path}")


def audit_batch(args: argparse.Namespace) -> dict[str, object]:
    batch_root = args.batch_root.resolve()
    plan = _load_object(batch_root / "batch.plan.json", "batch plan")
    receipt = _load_object(batch_root / "batch.receipt.json", "batch receipt")
    if plan.get("schema") != PLAN_SCHEMA or receipt.get("schema") != BATCH_RECEIPT_SCHEMA:
        raise DriverError("batch plan/receipt schema mismatch")
    if Path(str(plan.get("batch_root"))).resolve() != batch_root:
        raise DriverError("batch root does not match the frozen plan")
    execution_policy = plan.get("execution_policy")
    if not isinstance(execution_policy, dict):
        raise DriverError("batch plan does not contain an execution policy")
    if receipt.get("execution_policy") != execution_policy:
        raise DriverError("batch receipt execution policy differs from the plan")
    if receipt.get("concurrency_preflight") != plan.get("concurrency_preflight"):
        raise DriverError("batch receipt concurrency preflight differs from the plan")
    verify_execution_policy(execution_policy)
    plan_preflight = plan.get("concurrency_preflight")
    if not isinstance(plan_preflight, dict):
        raise DriverError("batch plan concurrency preflight is invalid")
    verify_gate_matches_policy(plan_preflight, execution_policy)
    if plan.get("allow_existing_pea_pid") != plan_preflight.get("allowed_existing_pea_pid"):
        raise DriverError("planned allowed PEA PID differs from the preflight gate")
    planned_environment = plan.get("environment")
    if not isinstance(planned_environment, dict):
        raise DriverError("batch plan execution environment is invalid")
    verify_effective_environment(planned_environment)
    static_files = plan.get("static_files")
    if not isinstance(static_files, dict) or set(static_files) != {
        "driver",
        "manifest",
        "python_interpreter",
        "raw_validator",
        "sd_validator",
        "restart_comparator",
    }:
        raise DriverError("batch plan static file set is invalid")
    for role, record in static_files.items():
        if not isinstance(record, dict):
            raise DriverError(f"batch plan static file record is invalid: {role}")
        _verify_record(record, f"planned static file {role}")
    plan_record = receipt.get("plan")
    if not isinstance(plan_record, dict):
        raise DriverError("batch receipt does not register its plan")
    _verify_record(plan_record, "batch plan")
    frozen_receipt = receipt.get("binary")
    if not isinstance(frozen_receipt, dict):
        raise DriverError("batch receipt does not register the frozen binary")
    _verify_record(frozen_receipt, "frozen binary")
    if int(frozen_receipt.get("mode", -1)) & 0o222:
        raise DriverError("frozen binary is writable")
    if Path(str(frozen_receipt.get("path"))).resolve() != Path(
        str(plan.get("frozen_binary"))
    ).resolve():
        raise DriverError("frozen binary path differs from the batch plan")
    planned_build_binary = plan.get("build_binary")
    if not isinstance(planned_build_binary, dict) or not _same_content_record(
        planned_build_binary, frozen_receipt
    ):
        raise DriverError("frozen binary bytes differ from the planned build binary")
    current_inputs = verify_input_manifest(
        Path(str(plan["static_files"]["manifest"]["path"])),  # type: ignore[index]
        Path(str(plan["data_root"])),
    )
    if not _same_snapshot(plan["input_snapshot_pre"], current_inputs):
        raise DriverError("current input snapshot differs from the batch plan")
    run_summaries: list[dict[str, object]] = []
    science_inputs: list[dict[str, object]] = []
    all_valid = receipt.get("status") == "SUCCESS"
    planned_runs = plan.get("runs")
    if not isinstance(planned_runs, list) or len(planned_runs) != len(RUN_SPECS):
        raise DriverError("batch plan does not contain all five runs")
    expected_conditions = [spec.condition for spec in RUN_SPECS]
    if [str(item.get("condition")) for item in planned_runs if isinstance(item, dict)] != (
        expected_conditions
    ):
        raise DriverError("batch plan does not preserve the required serial run order")
    planned_by_condition: dict[str, Mapping[str, object]] = {}
    for planned_run in planned_runs:
        if not isinstance(planned_run, dict):
            raise DriverError("invalid run entry in batch plan")
        condition = str(planned_run.get("condition"))
        if condition in planned_by_condition:
            raise DriverError(f"duplicate planned condition: {condition}")
        planned_by_condition[condition] = planned_run
        closure = planned_run.get("config_include_closure")
        if not isinstance(closure, dict):
            raise DriverError(f"configuration closure missing for {condition}")
        verify_config_closure(closure)
    if set(planned_by_condition) != set(expected_conditions):
        raise DriverError("batch plan condition set differs from the five-run protocol")
    registered_runs = receipt.get("runs")
    if not isinstance(registered_runs, list) or len(registered_runs) != len(RUN_SPECS):
        raise DriverError("batch receipt does not contain all five runs")
    if [str(item.get("condition")) for item in registered_runs if isinstance(item, dict)] != (
        expected_conditions
    ):
        raise DriverError("batch receipt does not preserve the required serial run order")
    for registered in registered_runs:
        if not isinstance(registered, dict):
            raise DriverError("invalid run registration in batch receipt")
        condition = str(registered.get("condition"))
        planned_run = planned_by_condition.get(condition)
        if planned_run is None:
            raise DriverError(f"unplanned run condition in batch receipt: {condition}")
        receipt_record = registered.get("receipt")
        if not isinstance(receipt_record, dict):
            raise DriverError("run receipt record is missing")
        _verify_record(receipt_record, f"{registered.get('condition')} receipt")
        run_receipt = _load_object(Path(str(receipt_record["path"])), "run receipt")
        if run_receipt.get("schema") != RUN_RECEIPT_SCHEMA:
            raise DriverError("run receipt schema mismatch")
        if run_receipt.get("condition") != condition:
            raise DriverError(f"run receipt condition mismatch: {condition}")
        if run_receipt.get("planned_execution_policy") != execution_policy:
            raise DriverError(f"execution policy mismatch in {condition} receipt")
        run_gate = run_receipt.get("concurrency")
        if not isinstance(run_gate, dict):
            raise DriverError(f"resource gate missing from {condition} receipt")
        verify_gate_matches_policy(run_gate, execution_policy)
        verify_allowed_process_identity(plan_preflight, run_gate)
        if run_gate.get("allowed_existing_pea_pid") != plan_preflight.get(
            "allowed_existing_pea_pid"
        ):
            raise DriverError(f"allowed PEA PID differs from the plan for {condition}")
        planned_closure = planned_run.get("config_include_closure")
        receipt_closure = run_receipt.get("config_include_closure")
        if not isinstance(planned_closure, dict) or not isinstance(receipt_closure, dict):
            raise DriverError(f"configuration closure missing from {condition} receipt")
        if not _same_config_closure(planned_closure, receipt_closure):
            raise DriverError(f"configuration closure mismatch for {condition}")
        pea = run_receipt.get("pea")
        if not isinstance(pea, dict) or not isinstance(pea.get("binary"), dict):
            raise DriverError(f"frozen binary record missing from {condition} receipt")
        if not _same_content_record(frozen_receipt, pea["binary"]):
            raise DriverError(f"frozen binary mismatch in {condition} receipt")
        if Path(str(pea["binary"].get("path"))).resolve() != Path(
            str(plan["frozen_binary"])
        ).resolve():
            raise DriverError(f"frozen binary path mismatch in {condition} receipt")
        planned_argv = list(planned_run.get("argv", []))
        if not planned_argv or Path(str(planned_argv[0])).resolve() != Path(
            str(plan["frozen_binary"])
        ).resolve():
            raise DriverError(f"planned PEA argv does not use the frozen binary for {condition}")
        if pea.get("argv") != planned_argv:
            raise DriverError(f"PEA argv differs from the frozen plan for {condition}")
        if pea.get("environment") != planned_environment:
            raise DriverError(f"PEA environment differs from the frozen plan for {condition}")
        execution = run_receipt.get("execution")
        if not isinstance(execution, dict):
            raise DriverError(f"execution record missing from {condition} receipt")
        expected_command = [*list(execution_policy["command_prefix"]), *planned_argv]
        if execution.get("command") != expected_command:
            raise DriverError(f"execution wrapper/command differs from the plan for {condition}")
        validators = run_receipt.get("validators")
        expected_validator_roles = {
            "forward_raw": "raw_validator",
            "forward_sd": "sd_validator",
            "smoothed_raw": "raw_validator",
            "smoothed_sd": "sd_validator",
        }
        if not isinstance(validators, dict) or set(validators) != set(
            expected_validator_roles
        ):
            raise DriverError(f"validator execution set is invalid for {condition}")
        planned_python = static_files["python_interpreter"]
        for validator_name, tool_role in expected_validator_roles.items():
            validator = validators[validator_name]
            planned_tool = static_files[tool_role]
            if not isinstance(validator, dict) or not isinstance(planned_tool, dict):
                raise DriverError(f"validator record is invalid: {condition}:{validator_name}")
            if validator.get("tool_role") != tool_role:
                raise DriverError(f"validator tool role mismatch: {condition}:{validator_name}")
            tool = validator.get("tool")
            interpreter = validator.get("python_interpreter")
            if not isinstance(tool, dict) or not isinstance(interpreter, dict):
                raise DriverError(f"validator tool identity is missing: {condition}:{validator_name}")
            if not _same_content_record(planned_tool, tool) or Path(
                str(tool.get("path"))
            ).resolve() != Path(str(planned_tool["path"])).resolve():
                raise DriverError(f"validator bytes/path mismatch: {condition}:{validator_name}")
            if not isinstance(planned_python, dict) or not _same_content_record(
                planned_python, interpreter
            ) or Path(str(interpreter.get("path"))).resolve() != Path(
                str(planned_python["path"])
            ).resolve():
                raise DriverError(
                    f"validator Python interpreter mismatch: {condition}:{validator_name}"
                )
            validator_command = validator.get("command")
            if not isinstance(validator_command, list) or len(validator_command) < 2:
                raise DriverError(f"validator command is invalid: {condition}:{validator_name}")
            if Path(str(validator_command[0])).resolve() != Path(
                str(planned_python["path"])
            ).resolve() or Path(str(validator_command[1])).resolve() != Path(
                str(planned_tool["path"])
            ).resolve():
                raise DriverError(f"validator command identity mismatch: {condition}:{validator_name}")
            if validator.get("returncode") != 0:
                raise DriverError(f"validator returned nonzero: {condition}:{validator_name}")
        artifacts = run_receipt.get("artifacts")
        if not isinstance(artifacts, dict):
            raise DriverError("run receipt artifact map is missing")
        for role, record in artifacts.items():
            if not isinstance(record, dict) or "sha256" not in record:
                raise DriverError(f"artifact record is invalid: {role}")
            _verify_record(record, f"{registered.get('condition')}:{role}")
        exit_code_record = artifacts.get("run.exit_code")
        if not isinstance(exit_code_record, dict):
            raise DriverError(f"run.exit_code is missing for {condition}")
        try:
            published_exit_code = int(
                Path(str(exit_code_record["path"])).read_text(encoding="ascii").strip()
            )
            receipt_exit_code = int(run_receipt.get("execution", {}).get("returncode"))  # type: ignore[union-attr]
        except (KeyError, OSError, TypeError, ValueError) as exc:
            raise DriverError(f"invalid exit-code receipt for {condition}") from exc
        if published_exit_code != receipt_exit_code:
            raise DriverError(f"exit-code receipt mismatch for {condition}")
        if receipt_exit_code != 0:
            raise DriverError(f"successful run receipt has nonzero exit code for {condition}")
        trace_audit = run_receipt.get("validation", {}).get("trace_audit", {})  # type: ignore[union-attr]
        valid = run_receipt.get("status") == "SUCCESS" and bool(
            run_receipt.get("validation", {}).get("all_valid")  # type: ignore[union-attr]
        )
        all_valid &= valid
        science_inputs.append(
            {
                "payload": {
                    "condition": condition,
                    "validation": {"trace_audit": trace_audit},
                }
            }
        )
        run_summaries.append(
            {
                "condition": registered.get("condition"),
                "valid": valid,
                "resolved_epoch_count": trace_audit.get("resolved_epoch_count"),
                "complete_rank_epoch_count": trace_audit.get("complete_rank_epoch_count"),
                "maximum_resolved_combination_count": trace_audit.get(
                    "maximum_resolved_combination_count"
                ),
                "multirow_runtime_coverage": trace_audit.get(
                    "multirow_runtime_coverage"
                ),
                "maximum_integer_candidate_count": trace_audit.get(
                    "maximum_integer_candidate_count"
                ),
            }
        )
    comparison = receipt.get("restart_comparison", {})
    comparison_status = comparison.get("status") if isinstance(comparison, dict) else None
    if not isinstance(comparison, dict):
        raise DriverError("restart comparison receipt is invalid")
    comparison_tool = comparison.get("tool")
    comparison_python = comparison.get("python_interpreter")
    planned_comparator = static_files["restart_comparator"]
    planned_python = static_files["python_interpreter"]
    if not all(
        isinstance(record, dict)
        for record in (
            comparison_tool,
            comparison_python,
            planned_comparator,
            planned_python,
        )
    ):
        raise DriverError("restart comparison tool identity is missing")
    if not _same_content_record(planned_comparator, comparison_tool) or Path(
        str(comparison_tool.get("path"))
    ).resolve() != Path(str(planned_comparator["path"])).resolve():
        raise DriverError("restart comparator bytes/path differ from the plan")
    if not _same_content_record(planned_python, comparison_python) or Path(
        str(comparison_python.get("path"))
    ).resolve() != Path(str(planned_python["path"])).resolve():
        raise DriverError("restart comparator Python interpreter differs from the plan")
    comparison_execution = comparison.get("execution")
    if not isinstance(comparison_execution, dict):
        raise DriverError("restart comparison execution record is missing")
    comparison_command = comparison_execution.get("command")
    if not isinstance(comparison_command, list) or len(comparison_command) < 2:
        raise DriverError("restart comparison command is invalid")
    if Path(str(comparison_command[0])).resolve() != Path(
        str(planned_python["path"])
    ).resolve() or Path(str(comparison_command[1])).resolve() != Path(
        str(planned_comparator["path"])
    ).resolve():
        raise DriverError("restart comparison command identity differs from the plan")
    if comparison_execution.get("returncode") != 0:
        raise DriverError("restart comparison returned nonzero")
    comparison_artifacts = comparison.get("artifacts")
    if not isinstance(comparison_artifacts, dict):
        raise DriverError("restart comparison artifacts are missing")
    for role, record in comparison_artifacts.items():
        if not isinstance(record, dict):
            raise DriverError(f"restart comparison artifact is invalid: {role}")
        _verify_record(record, f"restart_comparison:{role}")
    derived_scientific_assessment = _batch_scientific_assessment(
        science_inputs, comparison if isinstance(comparison, dict) else {}
    )
    registered_scientific_assessment = receipt.get("scientific_assessment")
    scientific_assessment_matches_receipt = (
        isinstance(registered_scientific_assessment, dict)
        and registered_scientific_assessment == derived_scientific_assessment
    )
    all_valid &= scientific_assessment_matches_receipt
    scientific_status = (
        str(derived_scientific_assessment["status"]) if all_valid else "AUDIT_FAILED"
    )
    report = {
        "schema": BATCH_AUDIT_SCHEMA,
        "status": "PASS" if all_valid else "FAIL",
        "batch_id": plan["batch_id"],
        "batch_root": str(batch_root),
        "inputs_match": True,
        "runs": run_summaries,
        "restart_comparison_status": comparison_status,
        "scientific_status": scientific_status,
        "scientific_assessment": derived_scientific_assessment,
        "scientific_assessment_matches_receipt": scientific_assessment_matches_receipt,
        "claim_limit": receipt.get("claim_limit"),
    }
    write_audit = bool(getattr(args, "write_audit", False))
    dry_run = bool(getattr(args, "dry_run", False))
    output_arg = getattr(args, "output", None)
    if output_arg is not None and not write_audit:
        raise DriverError("--output requires --write-audit")
    if dry_run and write_audit:
        raise DriverError("--dry-run cannot be combined with --write-audit")
    if not write_audit:
        return {"status": report["status"], "audit_written": False, "report": report}
    output = output_arg or (batch_root / "audit" / "batch.audit.json")
    output_record = atomic_write_json_no_replace(output, report)
    return {
        "status": report["status"],
        "audit_written": True,
        "audit": output_record,
        "report": report,
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    run_parser = subparsers.add_parser("run", help="verify, freeze and run all five conditions")
    run_parser.add_argument("--repository", type=Path, default=Path.cwd())
    run_parser.add_argument(
        "--data-root", type=Path, default=Path("/home/rx/GINAN/inputData")
    )
    run_parser.add_argument(
        "--output-parent",
        type=Path,
        default=Path("/home/rx/GINAN/inputData/outputs"),
    )
    run_parser.add_argument(
        "--input-manifest",
        type=Path,
        default=Path(
            "Docs/stecCovarianceFeasibility/experiment_1_2024_input_manifest.json"
        ),
    )
    run_parser.add_argument("--pea-binary", type=Path, default=Path("bin/pea"))
    run_parser.add_argument("--time-binary", type=Path, default=Path("/usr/bin/time"))
    run_parser.add_argument("--nice-binary", type=Path, default=Path("/usr/bin/nice"))
    run_parser.add_argument("--ionice-binary", type=Path, default=Path("/usr/bin/ionice"))
    run_parser.add_argument(
        "--ld-library-path", default="/home/rx/.local/boost-1.82/lib"
    )
    run_parser.add_argument("--batch-id")
    run_parser.add_argument("--allow-existing-pea-pid", type=int)
    run_parser.add_argument(
        "--concurrent-nice",
        type=int,
        choices=range(0, 20),
        default=19,
        help="nice value for this batch when --allow-existing-pea-pid is used",
    )
    run_parser.add_argument("--dry-run", action="store_true")
    audit_parser = subparsers.add_parser(
        "audit", help="read-only rehash and audit of an existing batch"
    )
    audit_parser.add_argument("--batch-root", type=Path, required=True)
    audit_parser.add_argument(
        "--write-audit",
        action="store_true",
        help="publish a no-replace JSON audit artifact (default: read-only stdout)",
    )
    audit_parser.add_argument(
        "--output", type=Path, help="output path used only with --write-audit"
    )
    audit_parser.add_argument(
        "--dry-run",
        action="store_true",
        help="deprecated read-only spelling retained for compatibility",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        result = run_batch(args) if args.command == "run" else audit_batch(args)
    except (DriverError, OSError, ValueError, TypeError) as exc:
        print(f"experiment1 corrected driver error: {exc}", file=sys.stderr)
        return 2
    print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
