#!/usr/bin/env python3
"""One explicitly approved, bounded comput first-epoch FLOAT diagnostic.

Uses the separately frozen v2 identity/cleanup helpers, never the unfinished
local-only sampling changes. No estimator, input, or original-run modification.
"""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import importlib
import json
import os
from pathlib import Path
import platform
import signal
import subprocess
import sys
import time

GIB = 1024 ** 3
EPOCH = "2024-07-17 00:00:00"
SERVER_LIMITS = {"maximum_address_space_bytes": 64 * GIB,
    "maximum_owned_tree_rss_bytes": 40 * GIB,
    "minimum_host_available_memory_bytes": 16 * GIB,
    "minimum_output_free_disk_bytes": 20 * GIB, "sample_interval_seconds": 1,
    "term_grace_seconds": 1, "kill_grace_seconds": 1}
FROZEN_PROFILE_SHA = "5bba49764d500cb54a1dc8c8c944a01119a274566acff3b970597e14e7856df2"


def utc_now():
    return datetime.now(timezone.utc).isoformat()


def load_helpers(directory):
    import hashlib
    directory = directory.resolve(strict=True)
    candidate = directory / "run_local_first_epoch_profile.py"
    if hashlib.sha256(candidate.read_bytes()).hexdigest() != FROZEN_PROFILE_SHA:
        raise RuntimeError("Server diagnostic requires the reviewed frozen v2 profile helper")
    sys.path.insert(0, str(directory))
    modules = [importlib.import_module(name) for name in
               ("run_stages", "run_network180_float", "run_local_first_epoch_profile")]
    if any(Path(module.__file__).resolve().parent != directory for module in modules):
        raise RuntimeError("Mixed live/frozen diagnostic helpers are forbidden")
    stages, network, profile = modules
    profile.LIMITS = dict(SERVER_LIMITS)
    return stages, network, profile


def debugger_prefix(root):
    return ["/lib64/ld-linux-x86-64.so.2", "--library-path", str(root / "lib"),
            str(root / "bin/gdb"), "-nx", "-nh", "--data-directory=" + str(root / "share/gdb")]


def numerical_command(plan, network):
    return ["/usr/bin/time", "-v", *plan["debugger_prefix"], "--quiet", "--batch",
        "--return-child-result", "-ex", "set pagination off", "-ex", "set confirm off",
        "-ex", "set startup-with-shell off", "-ex", "set print thread-events off",
        "-iex", "set auto-load python-scripts off", "-iex", "set auto-load local-gdbinit off",
        "-iex", "set auto-load safe-path /nonexistent-ginan-debugger-autoload",
        "-iex", "set debuginfod enabled off", "-ex", "handle SIGINT stop print nopass",
        "-ex", "run", "-ex", "thread 1", "-ex", "bt 32", "-ex", "thread apply all bt 8",
        "--args", *network.command(plan["binary"], plan["configs"], plan["pinned_input_contract"])]


def sampling_action(text, svd_at, capture_at, now):
    if "ZHANG_FLOAT_PHASE_END phase=FLOAT_RANK_JACOBI_SVD " in text:
        return "SVD_COMPLETED_BEFORE_STACK_SAMPLE"
    if svd_at is not None and capture_at is None and now - svd_at >= 15:
        return "REQUEST_SVD_STACK"
    return None


def read_argv(pid):
    return [part.decode(errors="surrogateescape") for part in
            (Path("/proc") / str(pid) / "cmdline").read_bytes().split(b"\0") if part]


def interrupt_debugger(records, prefix, profile):
    """PIDfd + exact planned argv + saved birth/SID, not a generic loader name."""
    for saved in records:
        descriptor = None
        try:
            argv = read_argv(saved["pid"])
            if argv[:len(prefix)] != prefix:
                continue
            descriptor = os.pidfd_open(saved["pid"])
            actual = profile.proc_record(saved["pid"])
            if actual is None or any(actual[key] != saved[key] for key in ("start_ticks", "session")):
                continue
            if read_argv(saved["pid"])[:len(prefix)] != prefix:
                continue
            signal.pidfd_send_signal(descriptor, signal.SIGINT, None, 0)
            return {"event": "SVD_15_SECOND_DEBUGGER_STACK_REQUEST", "utc": utc_now(),
                "pid": saved["pid"], "start_ticks": saved["start_ticks"], "session": saved["session"],
                "actual_command": saved["command"], "verified_argv_prefix": prefix, "signal": "SIGINT"}
        except ProcessLookupError:
            continue
        finally:
            if descriptor is not None:
                os.close(descriptor)
    raise RuntimeError("Exact owned portable debugger identity not found; refuse broad SIGINT")


def sample_resources(root, started, processes, profile):
    for row in processes:
        try:
            text = (Path("/proc") / str(row["pid"]) / "stat").read_text()
            fields = text[text.rfind(")") + 2:].split()
            if int(fields[19]) == row["start_ticks"]:
                row.update(cpu_user_ticks=int(fields[11]), cpu_system_ticks=int(fields[12]),
                    thread_count=int(fields[17]), clock_ticks_per_second=os.sysconf("SC_CLK_TCK"))
        except (OSError, ValueError, IndexError):
            row["cpu_sample_unavailable"] = True
    mem = profile.linux_memory()
    disk = os.statvfs(root)
    return {"elapsed_seconds": time.monotonic() - started, "utc": utc_now(), "processes": processes,
        "owned_tree_rss_bytes": sum(row["rss_bytes"] for row in processes),
        "host_available_memory_bytes": mem["MemAvailable"],
        "output_free_disk_bytes": disk.f_bavail * disk.f_frsize,
        "swap_used_bytes": mem["SwapTotal"] - mem["SwapFree"]}


def semantic_changes(old, new, path=()):
    if isinstance(old, dict) and isinstance(new, dict) and old.keys() == new.keys():
        return [change for key in old for change in semantic_changes(old[key], new[key], path + (key,))]
    return [] if old == new else [{"path": list(path), "old": old, "new": new}]


def configuration_equivalence(plan, reference, network):
    old_rows = reference["pinned_input_contract"]["configurations"]
    new_rows = plan["pinned_input_contract"]["configurations"]
    if len(old_rows) != 8 or len(new_rows) != 8:
        raise RuntimeError("Original and diagnostic configurations must have exactly eight layers")
    unchanged = []
    for index in [*range(6), 7]:
        old, new = old_rows[index], new_rows[index]
        actual_old = network.identity(Path(old["path"]))
        if ((actual_old["bytes"], actual_old["sha256"]) != (old["bytes"], old["sha256"]) or
                (new["bytes"], new["sha256"]) != (old["bytes"], old["sha256"])):
            raise RuntimeError("Frozen/coordinate configuration SHA differs from original server run")
        unchanged.append({"layer": index + 1, "sha256": new["sha256"]})
    for old, new in zip(reference["reviewed_coordinate_files"], plan["reviewed_coordinate_files"]):
        if network.identity(Path(old["path"]))["sha256"] != old["sha256"] or new["sha256"] != old["sha256"]:
            raise RuntimeError("Coordinate source/review SHA differs from original server run")
    original_overlay = network.identity(Path(old_rows[6]["path"]))
    if original_overlay["sha256"] != old_rows[6]["sha256"]:
        raise RuntimeError("Original server overlay changed since its pinned plan")
    changes = semantic_changes(json.loads(Path(old_rows[6]["path"]).read_text()),
                               json.loads(Path(new_rows[6]["path"]).read_text()))
    product = ("processing_options", "gnss_general", "zhang_pppar")
    allowed = {("processing_options", "epoch_control", "end_epoch"),
        product + ("checkpoint_capture_epochs",), product + ("checkpoint_output_directory",),
        product + ("product_filename",), product + ("product_covariance_filename",),
        ("outputs", "outputs_root"), ("outputs", "metadata", "config_description")}
    if any(tuple(change["path"]) not in allowed for change in changes):
        raise RuntimeError("Diagnostic overlay changes a model/threshold outside the isolation allowlist: " + json.dumps(changes))
    return {"passed": True, "unchanged_layers": unchanged, "allowed_overlay_changes": changes,
            "scope": "Only end/capture epoch and isolated output/identity fields may differ"}


def prepare(args, stages, network, profile):
    manifest = json.loads(args.source_manifest.read_text())
    if (manifest["binary_sha256"] != stages.file_sha(args.binary) or
            manifest["archive_sha256"] != stages.file_sha(args.source_archive)):
        raise RuntimeError("Actual producer/source archive SHA binding failed")
    expected = {row["path"]: row["sha256"] for row in manifest["source_files"]}
    for name in ("run_stages.py", "run_network180_float.py", "run_local_first_epoch_profile.py"):
        if stages.file_sha(args.frozen_drivers / name) != expected["experiments/multi_gnss_20261021/" + name]:
            raise RuntimeError("Frozen helper does not match full source archive: " + name)
    reference = json.loads(args.reference_plan.read_text())
    if not (reference["ratio"] == 2 and reference["all_nis"] is False and
            reference["ar_execution"] is False and reference["epoch_only_ar"] is True):
        raise RuntimeError("Original server reference is not reviewed ratio2/NIS-off FLOAT")
    plan, env = network.prepare(args.binary, args.audit, args.run_root, args.frozen_config,
        args.input_root, args.navigation, args.coordinate_overlay, args.coordinate_review,
        timeout=900, float_threads=28, omp_threads=28)
    root = Path(plan["run_root"])
    overlay = Path(plan["configs"][-2])
    value = json.loads(overlay.read_text())
    value["processing_options"]["epoch_control"]["end_epoch"] = EPOCH
    value["processing_options"]["gnss_general"]["zhang_pppar"]["checkpoint_capture_epochs"] = [EPOCH]
    stages.write_json(overlay, value)
    names = network.read_station_selection(args.frozen_config / stages.CONFIG_NAMES[1],
                                          json.loads(args.audit.read_text()))
    plan["pinned_input_contract"] = network.pin_inputs(args.input_root, args.navigation, names, plan["configs"])
    plan["configuration_equivalence"] = configuration_equivalence(plan, reference, network)
    current = {Path(row["path"]).name: row for group in ("observations", "navigation")
               for row in plan["pinned_input_contract"][group]}
    for group in ("observations", "navigation"):
        for row in reference["pinned_input_contract"][group]:
            actual = current.get(Path(row["path"]).name)
            if actual is None or (actual["bytes"], actual["sha256"]) != (row["bytes"], row["sha256"]):
                raise RuntimeError("Physical input differs from original server run: " + row["path"])
    verified = []
    for row in json.loads(args.input_manifest.read_text())["files"]:
        actual = current.get(Path(row["path"]).name) or network.identity(Path(row["path"]))
        if (actual["bytes"], actual["sha256"]) != (row["bytes"], row["sha256"]):
            raise RuntimeError("Model/input differs from server reference manifest: " + row["path"])
        verified.append(actual)
    for key, value in reference["environment"].items():
        if env.get(key) != value:
            raise RuntimeError("Original server numerical environment changed: " + key)
    env["ZHANG_FLOAT_PROFILE"] = "1"
    env["DEBUGINFOD_URLS"] = ""
    if str(args.gdb_root) in env.get("LD_LIBRARY_PATH", ""):
        raise RuntimeError("Debugger library path must not leak into numerical inferior")
    private_tools = [network.identity(path) for path in sorted(args.gdb_root.rglob("*")) if path.is_file()]
    plan.update(stage="SERVER_180_FIRST_EPOCH_PROFILE_ONLY", state="PREPARED_NOT_RUN_SERVER_APPROVAL_REQUIRED",
        expected_epoch_count=1, end_epoch=EPOCH, checkpoint_capture_epochs=[EPOCH],
        diagnostic_only=True, scientific_acceptance=False, local_host=platform.node(),
        runtime_supervision=dict(SERVER_LIMITS), svd_profile_deadline_seconds=15,
        debugger_stack_output_deadline_seconds=3, debugger_prefix=debugger_prefix(args.gdb_root),
        private_debugger_files=private_tools,
        fresh_input_equivalence={"passed": True, "files": len(verified),
            "total_bytes": sum(row["bytes"] for row in verified), "identities": verified},
        original_server_reference_plan=network.identity(args.reference_plan),
        producer_source_snapshot={"archive": network.identity(args.source_archive),
            "manifest": network.identity(args.source_manifest), "source_file_count": len(expected)},
        environment={key: val for key, val in env.items() if key.startswith("ZHANG_") or key.endswith("_NUM_THREADS")})
    plan["source_identity"].extend([network.identity(Path(profile.__file__)), network.identity(Path(__file__)),
        network.identity(args.source_archive), network.identity(args.source_manifest), *private_tools])
    plan["command"] = numerical_command(plan, network)
    stages.write_json(root / "stage_plan.json", plan)
    network.preflight(plan, env)
    stages.write_json(root / "prepared_server_profile_plan.json", plan)
    return plan


def execute(plan, approval_path, stages, network, profile):
    root = Path(plan["run_root"])
    original = json.loads((root / "prepared_server_profile_plan.json").read_text())
    current = json.loads((root / "stage_plan.json").read_text())
    if plan != original or current["state"] != "PREPARED_NOT_RUN_SERVER_APPROVAL_REQUIRED":
        raise RuntimeError("Only unused exact server profile plan may execute once")
    approval = json.loads(approval_path.read_text())
    required = ("approved_for_bounded_server_profile", "build_passed", "mandatory_tests_passed",
                "private_gdb_verified", "no_estimator_or_diagnostic_config_change")
    if any(approval.get(key) is not True for key in required):
        raise RuntimeError("Missing bounded server diagnostic approval")
    if (approval.get("host") != platform.node() or approval.get("binary_sha256") != plan["binary_sha256"]
            or approval.get("prepared_plan_sha256") != stages.file_sha(root / "prepared_server_profile_plan.json")
            or approval.get("runtime_supervision") != SERVER_LIMITS):
        raise RuntimeError("Server host/binary/plan/resource approval mismatch")
    network.verify_pins(plan)
    env = stages.environment("D")
    env.update(plan["environment"])
    env["DEBUGINFOD_URLS"] = ""
    network.preflight(plan, env)
    if profile.linux_memory()["MemAvailable"] < SERVER_LIMITS["minimum_host_available_memory_bytes"]:
        raise RuntimeError("Server memory reserve unavailable before launch")
    return supervise(plan, env, stages, profile)


def supervise(plan, env, stages, profile):
    root, started = Path(plan["run_root"]), time.monotonic()
    plan.update(state="SERVER_PROFILE_RUNNING", started_utc=utc_now())
    stages.write_json(root / "stage_plan.json", plan)
    child = identity = guard_error = None
    observed, events, stop_reasons = {}, [], []
    svd_at = capture_at = None
    svd_completed_before_sample = False
    remaining = []
    try:
        with (root / "runner.log").open("wb") as output:
            child = subprocess.Popen(plan["command"], cwd=plan["working_directory"], env=env,
                stdout=output, stderr=subprocess.STDOUT, start_new_session=True,
                preexec_fn=profile.hard_address_space_limit)
            identity = profile.proc_record(child.pid)
            if identity is None:
                raise RuntimeError("Owned server profiler initial identity unavailable")
            observed[(identity["pid"], identity["start_ticks"])] = identity
            stages.write_json(root / "profile_process.json", {"identity": identity, "command": plan["command"]})
            with (root / "runtime_resource_samples.jsonl").open("w", buffering=1) as samples:
                while child.poll() is None:
                    records = profile.owned_tree(child.pid, identity["start_ticks"])
                    if not records:
                        raise RuntimeError("Owned server resource monitoring unavailable")
                    observed.update({(row["pid"], row["start_ticks"]): row for row in records})
                    sample = sample_resources(root, started, records, profile)
                    samples.write(json.dumps(sample) + "\n")
                    stop_reasons = profile.violations(sample)
                    if sample["elapsed_seconds"] >= 900:
                        stop_reasons.append("PROFILE_WALLTIME_LIMIT")
                    if stop_reasons:
                        events.extend(profile.stop_owned(list(observed.values())))
                        break
                    text = (root / "runner.log").read_text(errors="replace")
                    if svd_at is None and "ZHANG_FLOAT_PHASE_BEGIN phase=FLOAT_RANK_JACOBI_SVD " in text:
                        svd_at = time.monotonic()
                    action = sampling_action(text, svd_at, capture_at, time.monotonic())
                    if action == "SVD_COMPLETED_BEFORE_STACK_SAMPLE":
                        svd_completed_before_sample = True
                        events.append({"event": action, "utc": utc_now(),
                            "note": "No later-stage stack may be labelled an SVD stack"})
                        events.extend(profile.stop_owned(list(observed.values())))
                        break
                    if action == "REQUEST_SVD_STACK":
                        events.append(interrupt_debugger(records, plan["debugger_prefix"], profile))
                        capture_at = time.monotonic()
                    if capture_at is not None and time.monotonic() - capture_at >= 3:
                        events.extend(profile.stop_owned(list(observed.values())))
                        break
                    # During the stack-output window the same hard guards run
                    # every 0.1 s; no unobserved three-second sleep.
                    time.sleep(0.1 if capture_at is not None else 1)
            child.wait(timeout=5)
    except BaseException as error:
        guard_error = type(error).__name__ + ": " + str(error)
    finally:
        if child is not None:
            if identity is None:
                try:
                    events.extend(profile.stop_unreaped_new_session(child))
                except BaseException as error:
                    guard_error = guard_error or "Initial identity cleanup: " + repr(error)
            else:
                try:
                    observed.update({(row["pid"], row["start_ticks"]): row for row in
                        profile.owned_tree(child.pid, identity["start_ticks"])})
                except BaseException as error:
                    guard_error = guard_error or "Cleanup discovery: " + repr(error)
            try:
                if profile.live_observed(list(observed.values())):
                    events.extend(profile.stop_owned(list(observed.values())))
            except BaseException as error:
                guard_error = guard_error or "Owned server cleanup: " + repr(error)
            try:
                child.wait(timeout=3)
            except BaseException as error:
                guard_error = guard_error or "Owned child not reaped: " + repr(error)
            try:
                remaining = profile.live_observed(list(observed.values()))
                if remaining or child.poll() is None:
                    guard_error = guard_error or "Owned server process remains live"
            except BaseException as error:
                guard_error = guard_error or "Cleanup verification: " + repr(error)
    try:
        text = (root / "runner.log").read_text(errors="replace") if (root / "runner.log").exists() else ""
    except BaseException as error:
        text = ""
        guard_error = guard_error or "Final log unavailable: " + repr(error)
    address_failure = any(marker in text for marker in ("std::bad_alloc", "Cannot allocate memory", "memory exhausted"))
    debugger_failure = "internal-error:" in text
    state = ("PROFILE_GUARD_FAILED" if guard_error else "PROFILE_RESOURCE_STOP" if stop_reasons or address_failure
        else "PROFILE_DEBUGGER_FAILED" if debugger_failure
        else "SVD_COMPLETED_BEFORE_STACK_SAMPLE" if svd_completed_before_sample
        else "PROFILE_SVD_STACK_CAPTURE_STOP" if capture_at is not None
        else "PROFILE_PROCESS_EXITED_NOT_SCIENTIFIC_ACCEPTANCE")
    plan.update(state=state, elapsed_seconds=time.monotonic() - started, finished_utc=utc_now(),
        exit_code=child.returncode if child else None, guard_failure_reason=guard_error,
        resource_stop_reasons=stop_reasons, svd_phase_observed=svd_at is not None,
        debugger_stack_requested=capture_at is not None, debugger_internal_error=debugger_failure,
        svd_completed_before_stack_sample=svd_completed_before_sample,
        svd_stack_attribution="Requires actual stack and phase-log review; request alone is not evidence",
        identity_verified_signal_events=events, owned_processes_still_live=remaining,
        child_reaped=child is not None and child.poll() is not None,
        possible_address_space_limit_failure=address_failure, scientific_acceptance=False,
        runner_sha256=stages.file_sha(root / "runner.log") if (root / "runner.log").exists() else None)
    stages.write_json(root / "stage_plan.json", plan)
    return plan


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("frozen-drivers", "binary", "audit", "run-root", "frozen-config", "input-root", "navigation",
            "coordinate-overlay", "coordinate-review", "reference-plan", "input-manifest", "source-archive",
            "source-manifest", "gdb-root", "prepared-plan", "approval"):
        parser.add_argument("--" + name, type=Path)
    parser.add_argument("--execute", action="store_true")
    args = parser.parse_args()
    if args.frozen_drivers is None:
        parser.error("Explicit frozen v2 drivers required")
    stages, network, profile = load_helpers(args.frozen_drivers)
    if args.execute:
        if args.prepared_plan is None or args.approval is None:
            parser.error("Execution requires exact prepared plan and host approval")
        plan = execute(json.loads(args.prepared_plan.read_text()), args.approval, stages, network, profile)
    else:
        if any(getattr(args, name) is None for name in ("binary", "audit", "run_root", "frozen_config",
                "input_root", "navigation", "coordinate_overlay", "coordinate_review", "reference_plan",
                "input_manifest", "source_archive", "source_manifest", "gdb_root")):
            parser.error("All explicit preparation inputs required")
        plan = prepare(args, stages, network, profile)
    print(json.dumps({key: plan.get(key) for key in ("state", "run_root", "binary_sha256",
        "expected_epoch_count", "svd_phase_observed", "resource_stop_reasons", "guard_failure_reason")}, indent=2))


if __name__ == "__main__":
    main()
