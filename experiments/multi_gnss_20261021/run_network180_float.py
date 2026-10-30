#!/usr/bin/env python3
"""Prepare/execute the original 180-station cold G/E/C static-coordinate FLOAT.

Preparation never starts PEA. Native paths are supplied explicitly so the same
driver works on comput and WSL. A parser dry-run, reviewed coordinate overlay,
resource approval and multi-system checkpoint acceptance are execution gates.
"""
from __future__ import annotations

import argparse
import copy
from datetime import datetime, timezone
import json
import math
import os
from pathlib import Path
import re
import signal
import subprocess
import sys
import time

import run_stages as stages

STATION_CONFIG_SHA256 = "4470e50d4139c6ccdc3abd6e177848519561a070f17c58d348e857606d5c4aba"
CHECKPOINT_EPOCHS = ["2024-07-17 01:59:30", "2024-07-17 02:00:00"]
MODE = "STATIC_ESTIMATED_GEC_FLOAT_180"
GIB = 1024 ** 3
RESOURCE_LIMITS = {"maximum_owned_tree_rss_bytes": 85 * GIB,
    "minimum_host_available_memory_bytes": 12 * GIB, "minimum_output_free_disk_bytes": 20 * GIB,
    "sample_interval_seconds": 10, "sigint_epoch_boundary_grace_seconds": 120}


def read_station_selection(path, audit):
    if stages.file_sha(path) != STATION_CONFIG_SHA256:
        raise RuntimeError("Original 180 estimation-station YAML SHA-256 mismatch")
    if audit["source_station_config"]["sha256"] != STATION_CONFIG_SHA256:
        raise RuntimeError("Input audit does not belong to the original service station list")
    names = re.findall(r"^    - (.+\.rnx)$", path.read_text(), re.MULTILINE)
    ids = [name[:4] for name in names]
    if len(names) != 180 or len(set(names)) != 180 or len(set(ids)) != 180:
        raise RuntimeError("Require all original 180 unique estimation station files/IDs")
    if any(Path(name).name != name or re.search(r"[\s*?\[\]]", name) for name in names):
        raise RuntimeError("Station files must be exact basenames, not paths/globs")
    header_names = [row["filename"] for row in audit["station_headers"]]
    if len(header_names) != 180 or set(header_names) != set(names):
        raise RuntimeError("180-station header audit differs from the frozen estimation list")
    return names


def coordinate_approval(overlay_path, review_path):
    review = json.loads(review_path.read_text())
    required = ("approved_for_static_coordinate_estimation", "all_service_receivers_estimated",
                "no_pseudo_coordinate_constraints")
    if any(review.get(key) is not True for key in required) or review.get("mode") != "STATIC_ESTIMATED":
        raise RuntimeError("Coordinate review is missing a static estimated-coordinate gate")
    if review.get("coordinate_overlay_sha256") != stages.file_sha(overlay_path):
        raise RuntimeError("Coordinate overlay differs from its reviewed SHA-256")
    value = json.loads(overlay_path.read_text())
    if set(value) - {"estimation_parameters", "processing_options", "receiver_options", "outputs"}:
        raise RuntimeError("Coordinate overlay changes unrelated top-level configuration")
    pos = value.get("estimation_parameters", {}).get("receivers", {}).get("global", {}).get("pos", {})
    if pos.get("estimated") != [True] or pos.get("process_noise") != [0]:
        raise RuntimeError("Coordinate overlay must explicitly estimate static positions with original Q=0")
    if pos.get("sigma") != [1]:
        raise RuntimeError("Coordinate overlay must preserve the original 1 m position prior")
    receiver = value["estimation_parameters"]["receivers"]["global"]
    if receiver.get("pos_rate", {}).get("estimated") != [False]:
        raise RuntimeError("Static coordinate trial must explicitly disable receiver position rates")
    if value.get("processing_options", {}).get("minimum_constraints", {}).get("enable") is not False:
        raise RuntimeError("Minimum-constraint pseudo position fixing must remain disabled")
    return review, value


def network_overlay(root, input_root, navigation, audit, names):
    adapted = copy.deepcopy(audit)
    # Existing reviewed D overlay supplies G/E/C signals/satellite membership.
    # Its selection is private here; existing 12-station defaults are unchanged.
    adapted["selected_12_station_candidates"] = [name[:4] for name in names]
    value = stages.overlay("D", root, adapted)
    value["inputs"]["gnss_observations"].update(
        gnss_observations_root=str(input_root / "data") + "/", rnx_inputs=names)
    value["inputs"]["satellite_data"]["nav_files"] = [str(navigation)]
    epoch = value["processing_options"]["epoch_control"]
    epoch.update(start_epoch="2024-07-17 00:00:00", end_epoch="2024-07-17 02:00:00", epoch_interval=30)
    product = value["processing_options"]["gnss_general"]["zhang_pppar"]
    product.update(deterministic_checkpoint=True, checkpoint_runtime_id="E29-GEC180-FREECOORD-20261021",
        checkpoint_output_directory=str(root / "checkpoints"), checkpoint_capture_epochs=CHECKPOINT_EPOCHS,
        checkpoint_restore_path="")
    value["outputs"]["metadata"]["config_description"] = "multi_gnss_180_freecoord_20261021"
    return value


def identity(path, **extra):
    path = Path(path).resolve(strict=True)
    if any(c.isspace() for c in str(path)):
        raise RuntimeError("Production effective-path summary requires paths without whitespace")
    return dict(path=str(path), bytes=path.stat().st_size, sha256=stages.file_sha(path), **extra)


def pin_inputs(input_root, navigation, names, configs):
    return {"selection_policy": "ORIGINAL_180_ESTIMATION_LIST_EXACT_CLI_REPLACEMENT",
        "expected_unique_station_count": 180,
        "observations": [identity(input_root / "data" / name, station=name[:4]) for name in names],
        "navigation": [identity(navigation)], "configurations": [identity(path) for path in configs],
        "cnv123_propagation_claim": False, "scientific_acceptance": False}


def exact_input_cli(pinned):
    paths = [row["path"] for row in pinned["observations"]]
    nav = [row["path"] for row in pinned["navigation"]]
    stations = [row["station"] for row in pinned["observations"]]
    if len(paths) != 180 or len(set(paths)) != 180 or len(set(stations)) != 180 or len(nav) != 1:
        raise RuntimeError("Need exactly 180 unique service receivers and one mixed NAV")
    return ["--rnx_inputs", *paths, "--nav_files", *nav]


def command(binary, configs, pinned, dry_run=False):
    result = [str(binary), "-y", *map(str, configs), "-d", "case", *exact_input_cli(pinned)]
    return result + (["--dry-run"] if dry_run else [])


def verify_pins(plan):
    exact_input_cli(plan["pinned_input_contract"])
    for group in ("observations", "navigation", "configurations"):
        for row in plan["pinned_input_contract"][group]:
            path = Path(row["path"])
            if not path.is_file() or path.stat().st_size != row["bytes"] or stages.file_sha(path) != row["sha256"]:
                raise RuntimeError("Pinned input changed: " + str(path))
    if stages.file_sha(plan["binary"]) != plan["binary_sha256"]:
        raise RuntimeError("Pinned production binary changed")
    for row in [*plan["reviewed_coordinate_files"], plan["input_audit"], *plan["source_identity"]]:
        if stages.file_sha(row["path"]) != row["sha256"]:
            raise RuntimeError("Reviewed coordinate/audit/driver source changed")


def verify_effective_summary(text, pinned):
    found = stages.parse_effective_input_summary(text)
    for key, group, count in (("rnx_inputs", "observations", 180), ("nav_files", "navigation", 1)):
        actual, expected = found[key], {row["path"] for row in pinned[group]}
        if len(actual) != count or len(set(actual)) != count or set(actual) != expected:
            raise RuntimeError("Actual parser effective " + key + " differs from pinned 180+1 inputs")
    if any(value for key, value in found.items() if key not in ("rnx_inputs", "nav_files")):
        raise RuntimeError("Unexpected extra observation/navigation streams")
    return {"passed": True, "actual_unique_observation_files": 180, "actual_navigation_files": 1,
            "scope": "Real parser only; no numerical observations or scientific acceptance"}


def prepare(binary, audit_path, root, frozen_config, input_root, navigation, coordinate_overlay,
            coordinate_review, timeout=172800, float_threads=12, omp_threads=8):
    if not math.isfinite(timeout) or timeout <= 0:
        raise ValueError("Timeout must be finite and positive")
    if not all(isinstance(count, int) and 1 <= count <= 1024 for count in (float_threads, omp_threads)):
        raise ValueError("Thread counts must be positive integers at most 1024")
    audit = json.loads(audit_path.read_text())
    if audit.get("missing_files"):
        raise RuntimeError("Audited service input files are missing")
    binary, input_root, frozen_config = (Path(path).resolve(strict=True) for path in
                                        (binary, input_root, frozen_config))
    names = read_station_selection(frozen_config / stages.CONFIG_NAMES[1], audit)
    review, coordinate_value = coordinate_approval(coordinate_overlay, coordinate_review)
    root = root.resolve()
    root.mkdir(parents=True, exist_ok=False)
    for directory in ("config", "outputs", "checkpoints"):
        (root / directory).mkdir()
    overlay_path = root / "config/network180_overlay.yaml"
    stages.write_json(overlay_path, network_overlay(root, input_root, navigation, audit, names))
    coordinate_copy = root / "config/reviewed_static_coordinate_overlay.yaml"
    stages.write_json(coordinate_copy, coordinate_value)
    configs = [frozen_config / name for name in stages.CONFIG_NAMES] + [overlay_path, coordinate_copy]
    pinned = pin_inputs(input_root, navigation, names, configs)
    env = stages.environment("D")
    env.update(ZHANG_R51_FLOAT_THREADS=str(float_threads), OMP_NUM_THREADS=str(omp_threads),
               ZHANG_MULTI_GNSS_FLOAT_CHECKPOINT="1")
    if any("RESTORE" in key or "REUSE_FLOAT" in key for key in env):
        raise RuntimeError("Cold network trial cannot import an old checkpoint")
    plan = {"schema_version": 1, "stage": "NETWORK180_GEC_FREECOORD_FLOAT", "mode": MODE,
        "state": "PREPARED_NOT_RUN", "run_root": str(root), "binary": str(binary),
        "binary_sha256": stages.file_sha(binary), "working_directory": str(input_root),
        "configs": list(map(str, configs)), "input_audit": identity(audit_path),
        "original_station_config_sha256": STATION_CONFIG_SHA256,
        "station_ids": [name[:4] for name in names], "station_selection_is_not_validation_network": True,
        "start_epoch": "2024-07-17 00:00:00", "end_epoch": "2024-07-17 02:00:00",
        "expected_epoch_count": 241, "checkpoint_capture_epochs": CHECKPOINT_EPOCHS,
        "cold_start": True, "reuse_gps_checkpoint": False, "ratio": 2, "all_nis": False,
        "epoch_only_ar": True, "ar_execution": False, "timeout_seconds": timeout,
        "coordinate_review": review, "reviewed_coordinate_files": [identity(coordinate_overlay), identity(coordinate_review)],
        "scientific_acceptance": False, "pinned_input_contract": pinned,
        "execution_gate_note": "Requires root-reviewed multi-system checkpoint roundtrip and host-bound resource approval; preparation is not acceptance",
        "effective_input_preflight": {"status": "NOT_RUN", "required_before_solver_launch": True},
        "environment": {key: value for key, value in env.items()
                        if key.startswith("ZHANG_") or key.endswith("_NUM_THREADS")},
        "source_identity": [identity(Path(__file__)), identity(Path(stages.__file__))]}
    plan["runtime_supervision"] = dict(RESOURCE_LIMITS,
        policy="Owned process tree only; first SIGINT for epoch-barrier stop, then bounded TERM/KILL escalation",
        fatal_qc_policy="Only final authoritative PPP transaction rejection/terminal fatal from this unique run's logs; ordinary retries and station SPP failures are not fatal")
    plan["command"] = ["/usr/bin/time", "-v", *command(binary, configs, pinned)]
    stages.write_json(root / "stage_plan.json", plan)
    return plan, env


def preflight(plan, env):
    root = Path(plan["run_root"])
    verify_pins(plan)
    cmd = command(plan["binary"], plan["configs"], plan["pinned_input_contract"], True)
    result = subprocess.run(cmd, cwd=plan["working_directory"], env=env, text=True,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=180)
    log = root / "config_preflight.log"
    log.write_text(result.stdout)
    report = {"command": cmd, "exit_code": result.returncode, "log": str(log), "log_sha256": stages.file_sha(log)}
    if result.returncode or "Dry-run successful: configuration and sanity checks passed" not in result.stdout:
        report.update(passed=False, reason="Real PEA parser/configuration rejected trial; see log, do not bypass checkpoint gate")
    else:
        report.update(verify_effective_summary(result.stdout, plan["pinned_input_contract"]))
    plan["effective_input_preflight"] = report
    if not report["passed"]:
        plan["state"] = "CONFIGURATION_PREFLIGHT_FAILED_NOT_RUN"
    stages.write_json(root / "stage_plan.json", plan)
    verify_pins(plan)
    if not report["passed"]:
        raise RuntimeError(report["reason"])
    return report


def execution_approval(plan, path):
    review = json.loads(path.read_text())
    required = ("build_passed", "mandatory_tests_passed", "approved_for_execution",
                "multi_system_checkpoint_roundtrip_passed", "static_coordinate_configuration_review_passed",
                "approved_for_180_float_trial")
    if any(review.get(key) is not True for key in required):
        raise RuntimeError("Execution blocked: missing build/checkpoint/coordinate/resource approval")
    if review.get("binary_sha256") != plan["binary_sha256"]:
        raise RuntimeError("Execution approval binary SHA-256 mismatch")
    host = os.uname().nodename
    if review.get("host") != host:
        raise RuntimeError("Resource approval is for a different native host: " + host)
    if not all(isinstance(review.get(key), int) and review[key] > 0 for key in
               ("minimum_available_memory_bytes", "minimum_free_disk_bytes")):
        raise RuntimeError("Host resource approval must supply positive memory/disk reserve gates")
    memory = Path("/proc/meminfo").read_text()
    available = int(re.search(r"^MemAvailable:\s+(\d+)\s+kB$", memory, re.MULTILINE)[1]) * 1024
    stat = os.statvfs(plan["run_root"])
    free_disk = stat.f_bavail * stat.f_frsize
    if available < review["minimum_available_memory_bytes"] or free_disk < review["minimum_free_disk_bytes"]:
        raise RuntimeError("Native host does not meet the reviewed RAM/disk reserve gates")
    return dict(review, measured_available_memory_bytes=available, measured_free_disk_bytes=free_disk)


def proc_identity(pid):
    """Linux PID birth identity. A recycled PID must never receive our signal."""
    text = Path(f"/proc/{pid}/stat").read_text()
    fields = text.rsplit(")", 1)[1].split()
    return {"pid": pid, "parent_pid": int(fields[1]), "start_ticks": int(fields[19])}


def owned_process_tree(pid):
    records = {}
    for directory in Path("/proc").iterdir():
        if not directory.name.isdecimal():
            continue
        try:
            row = proc_identity(int(directory.name))
            status = (directory / "status").read_text()
            match = re.search(r"^VmRSS:\s+(\d+)\s+kB$", status, re.MULTILINE)
            row["rss_bytes"] = int(match[1]) * 1024 if match else 0
            records[row["pid"]] = row
        except (OSError, ValueError, IndexError):
            continue  # An unrelated /proc process may end during enumeration.
    if pid not in records:
        raise RuntimeError("Owned child identity unavailable during live resource sample")
    selected = {pid}
    while True:
        expanded = selected | {row["pid"] for row in records.values() if row["parent_pid"] in selected}
        if expanded == selected:
            break
        selected = expanded
    return [records[item] for item in sorted(selected)]


def resource_sample(pid, root):
    tree = owned_process_tree(pid)
    memory = Path("/proc/meminfo").read_text()
    match = re.search(r"^MemAvailable:\s+(\d+)\s+kB$", memory, re.MULTILINE)
    if not match:
        raise RuntimeError("Linux MemAvailable is unavailable")
    stat = os.statvfs(root)
    return {"sampled_utc": datetime.now(timezone.utc).isoformat(), "processes": tree,
        "owned_tree_rss_bytes": sum(row["rss_bytes"] for row in tree),
        "host_available_memory_bytes": int(match[1]) * 1024,
        "output_free_disk_bytes": stat.f_bavail * stat.f_frsize}


def resource_violations(sample):
    reasons = []
    if sample["owned_tree_rss_bytes"] >= RESOURCE_LIMITS["maximum_owned_tree_rss_bytes"]:
        reasons.append("OWNED_TREE_RSS_AT_OR_ABOVE_85_GIB")
    if sample["host_available_memory_bytes"] < RESOURCE_LIMITS["minimum_host_available_memory_bytes"]:
        reasons.append("HOST_MEM_AVAILABLE_BELOW_12_GIB")
    if sample["output_free_disk_bytes"] < RESOURCE_LIMITS["minimum_output_free_disk_bytes"]:
        reasons.append("OUTPUT_DISK_AVAILABLE_BELOW_20_GIB")
    return reasons


def authoritative_failure(line):
    # These are emitted only after the authoritative /PPP factor has been
    # rejected, or after mainOncePerEpoch returned false. Neither per-measurement
    # QC warnings nor an exhausted single-station SPP are a network rollback.
    return ("PPP AUTHORITATIVE MEASUREMENT TRANSACTION REJECTED" in line or
            "AUTHORITATIVE_PPP_FACTOR_TRANSACTION_FAILED:" in line)


class CurrentRunLogTail:
    def __init__(self, root):
        self.root = Path(root).resolve()
        self.offsets = {}

    def final_failures(self):
        failures = []
        paths = [self.root / "runner.log", *(self.root / "outputs").glob("Network-*.TRACE")]
        for path in paths:
            if not path.is_file() or path.is_symlink() or self.root not in path.resolve().parents:
                continue
            offset, pending = self.offsets.get(str(path), (0, b""))
            if path.stat().st_size < offset:
                raise RuntimeError("Current-run log was truncated during resource supervision")
            with path.open("rb") as stream:
                stream.seek(offset)
                data = pending + stream.read()
                offset = stream.tell()
            lines = data.splitlines(keepends=True)
            pending = b""
            for line in lines:
                # Fatal TRACE can be the final un-newlined record. It is safe
                # to recognize a complete literal marker before line flush.
                text = line.decode("utf-8", errors="replace")
                if authoritative_failure(text):
                    failures.append({"path": str(path), "line": text.strip()})
                elif not line.endswith((b"\n", b"\r")):
                    pending = line
            self.offsets[str(path)] = (offset, pending)
        return failures


def stop_owned_child(child, birth_identity):
    report = {"policy": "SIGINT epoch barrier, bounded 120 s, then TERM 10 s then KILL owned group",
              "signals": [], "requested_utc": datetime.now(timezone.utc).isoformat()}
    for sig, grace in ((signal.SIGINT, 120), (signal.SIGTERM, 10), (signal.SIGKILL, 10)):
        if child.poll() is not None:
            break
        current = proc_identity(child.pid)
        if current["start_ticks"] != birth_identity["start_ticks"] or os.getpgid(child.pid) != child.pid:
            raise RuntimeError("Refuse to signal a changed/recycled/non-owned process group")
        os.killpg(child.pid, sig)
        report["signals"].append({"signal": int(sig), "sent_utc": datetime.now(timezone.utc).isoformat(),
                                   "maximum_wait_seconds": grace})
        try:
            child.wait(timeout=grace)
        except subprocess.TimeoutExpired:
            continue
    if child.poll() is None:
        raise RuntimeError("Owned child did not terminate after bounded signal escalation")
    report["exit_code"] = child.returncode
    report["finished_utc"] = datetime.now(timezone.utc).isoformat()
    return report


def execute(plan, env, approval_path):
    if sys.platform != "linux":
        raise RuntimeError("Execute only on the reviewed native Linux compute host")
    approved = execution_approval(plan, approval_path)
    preflight(plan, env)
    verify_pins(plan)
    root = Path(plan["run_root"])
    plan.update(state="RUNNING", execution_approval=approved, started_utc=datetime.now(timezone.utc).isoformat())
    stages.write_json(root / "stage_plan.json", plan)
    started = time.monotonic()
    timed_out = False
    supervision_reasons, final_failures, stop_report = [], [], None
    log_tail = CurrentRunLogTail(root)
    with (root / "runner.log").open("wb") as stream:
        child = subprocess.Popen(plan["command"], cwd=plan["working_directory"], env=env,
                                 stdout=stream, stderr=subprocess.STDOUT, start_new_session=True)
        plan["pid"] = child.pid
        stages.write_json(root / "stage_plan.json", plan)
        birth_identity = proc_identity(child.pid)
        with (root / "runtime_resource_samples.jsonl").open("w", buffering=1) as sample_log:
            while child.poll() is None:
                if time.monotonic() - started >= plan["timeout_seconds"]:
                    timed_out = True
                    supervision_reasons.append("EXECUTION_TIMEOUT")
                try:
                    sample = resource_sample(child.pid, root)
                    if sample["processes"][0]["pid"] != child.pid:
                        # Children can have lower numeric PIDs after PID wrap.
                        own = next(row for row in sample["processes"] if row["pid"] == child.pid)
                    else:
                        own = sample["processes"][0]
                    if own["start_ticks"] != birth_identity["start_ticks"]:
                        raise RuntimeError("Owned process identity changed")
                    sample_log.write(json.dumps(sample) + "\n")
                    supervision_reasons.extend(resource_violations(sample))
                    final_failures.extend(log_tail.final_failures())
                    if final_failures:
                        supervision_reasons.append("FINAL_AUTHORITATIVE_FLOAT_TRANSACTION_REJECTED")
                except (OSError, RuntimeError, ValueError, IndexError) as error:
                    if child.poll() is None:
                        supervision_reasons.append("RESOURCE_OR_CURRENT_LOG_SUPERVISION_UNAVAILABLE:" + str(error))
                if supervision_reasons:
                    stop_report = stop_owned_child(child, birth_identity)
                    break
                try:
                    child.wait(timeout=RESOURCE_LIMITS["sample_interval_seconds"])
                except subprocess.TimeoutExpired:
                    pass
        child.wait()
    final_failures.extend(log_tail.final_failures())
    succeeded = child.returncode == 0 and not timed_out and not supervision_reasons and not final_failures
    plan.update(state="COMPLETED_REQUIRES_SCIENTIFIC_AND_CHECKPOINT_AUDIT" if succeeded else "EXECUTION_FAILED",
        exit_code=child.returncode, timed_out=timed_out, elapsed_seconds=time.monotonic() - started,
        finished_utc=datetime.now(timezone.utc).isoformat(), scientific_acceptance=False,
        supervision_stop_reasons=supervision_reasons, authoritative_float_failures=final_failures,
        stop_signal_report=stop_report)
    stages.write_json(root / "stage_plan.json", plan)
    verify_pins(plan)
    return succeeded


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for option in ("binary", "audit", "run-root", "frozen-config-dir", "input-root", "navigation",
                   "coordinate-overlay", "coordinate-review"):
        parser.add_argument("--" + option, type=Path, required=True)
    parser.add_argument("--timeout-seconds", type=float, default=172800)
    parser.add_argument("--float-threads", type=int, default=12)
    parser.add_argument("--omp-threads", type=int, default=8)
    parser.add_argument("--build-approval", type=Path)
    parser.add_argument("--preflight", action="store_true")
    parser.add_argument("--execute", action="store_true")
    args = parser.parse_args()
    if args.execute and args.build_approval is None:
        parser.error("--execute requires --build-approval; prepare is not permission to run")
    plan, env = prepare(args.binary, args.audit, args.run_root, args.frozen_config_dir, args.input_root,
        args.navigation, args.coordinate_overlay, args.coordinate_review, args.timeout_seconds,
        args.float_threads, args.omp_threads)
    if args.execute:
        try:
            passed = execute(plan, env, args.build_approval)
        except Exception as error:
            plan.update(state="EXECUTION_FAILED", failure_reason=str(error),
                        finished_utc=datetime.now(timezone.utc).isoformat(), scientific_acceptance=False)
            stages.write_json(Path(plan["run_root"]) / "stage_plan.json", plan)
            raise
        raise SystemExit(0 if passed else 1)
    if args.preflight:
        preflight(plan, env)
    print(json.dumps(plan, indent=2))


if __name__ == "__main__":
    main()
