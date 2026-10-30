#!/usr/bin/env python3
"""Prepare or explicitly approve one owned local 180-station FLOAT profile.

No server operation, WSL restart, input download, or estimator modification.
Original physical inputs/coordinates/diagnostics are preserved. Resource stops
are diagnostic outcomes, never numerical QC failures or scientific acceptance.
"""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import json
import math
import os
from pathlib import Path
import platform
try:
    import resource
except ImportError:  # Module inspection/tests may run on Windows; execution may not.
    resource = None
import signal
import subprocess
import time

import run_network180_float as network
import run_stages as stages

GIB = 1024 ** 3
EPOCH = "2024-07-17 00:00:00"
LIMITS = {"maximum_address_space_bytes": 16 * GIB, "maximum_owned_tree_rss_bytes": 12 * GIB,
    "minimum_host_available_memory_bytes": 4 * GIB, "minimum_output_free_disk_bytes": 20 * GIB,
    "sample_interval_seconds": 1, "term_grace_seconds": 1, "kill_grace_seconds": 1}


def stamp():
    return datetime.now(timezone.utc).isoformat()


def profile_environment(values):
    result = {}
    for item in values:
        key, separator, value = item.partition("=")
        if not separator or not key.startswith("ZHANG_FLOAT_PROFILE") or not value or key in result:
            raise ValueError("Profile environment must be unique ZHANG_FLOAT_PROFILE*=VALUE flags")
        if any(ord(char) < 32 or ord(char) > 126 for char in item):
            raise ValueError("Profile environment must be printable ASCII")
        result[key] = value
    return result


def linux_memory():
    fields = {}
    for line in Path("/proc/meminfo").read_text().splitlines():
        key, value = line.split(":", 1)
        fields[key] = int(value.split()[0]) * 1024
    return {key: fields[key] for key in ("MemTotal", "MemAvailable", "SwapTotal", "SwapFree")}


def proc_record(pid):
    try:
        text = (Path("/proc") / str(pid) / "stat").read_text()
        fields = text[text.rfind(")") + 2:].split()
        return {"pid": int(pid), "parent_pid": int(fields[1]), "process_group": int(fields[2]),
            "session": int(fields[3]), "start_ticks": int(fields[19]), "virtual_bytes": int(fields[20]),
            "rss_bytes": int(fields[21]) * os.sysconf("SC_PAGE_SIZE"),
            "state": fields[0],
            "command": text[text.find("(") + 1:text.rfind(")")]}
    except (FileNotFoundError, ProcessLookupError, PermissionError, ValueError):
        return None


def owned_tree(pid, start_ticks):
    records = {int(path.name): record for path in Path("/proc").iterdir() if path.name.isdigit()
               and (record := proc_record(int(path.name))) is not None}
    if pid not in records or records[pid]["start_ticks"] != start_ticks:
        return []
    selected = {pid}
    while True:
        added = {key for key, value in records.items() if value["parent_pid"] in selected}
        if added <= selected:
            break
        selected |= added
    return [records[key] for key in sorted(selected)]


def violations(sample):
    result = []
    if sample["owned_tree_rss_bytes"] >= LIMITS["maximum_owned_tree_rss_bytes"]:
        result.append("OWNED_RSS_LIMIT")
    if sample["host_available_memory_bytes"] < LIMITS["minimum_host_available_memory_bytes"]:
        result.append("HOST_MEMAVAILABLE_RESERVE")
    if sample["output_free_disk_bytes"] < LIMITS["minimum_output_free_disk_bytes"]:
        result.append("OUTPUT_DISK_RESERVE")
    return result


def stop_owned(records):
    """Signal only identity-verified observed owned PIDs, never another job."""
    events = []
    for sig, wait in ((signal.SIGTERM, LIMITS["term_grace_seconds"]),
                      (signal.SIGKILL, LIMITS["kill_grace_seconds"])):
        for saved in reversed(records):
            actual = proc_record(saved["pid"])
            if actual is None or any(actual[key] != saved[key] for key in ("start_ticks", "session")):
                continue
            try:
                os.kill(saved["pid"], sig)
                events.append({"pid": saved["pid"], "start_ticks": saved["start_ticks"], "signal": sig.name})
            except ProcessLookupError:
                pass
            except OSError as error:
                events.append({"pid": saved["pid"], "start_ticks": saved["start_ticks"],
                               "signal": sig.name, "signal_error": str(error)})
        time.sleep(wait)
    return events


def hard_address_space_limit():
    if resource is None:
        raise RuntimeError("RLIMIT_AS requires native Linux; no unbounded fallback")
    resource.setrlimit(resource.RLIMIT_AS, (LIMITS["maximum_address_space_bytes"],) * 2)
    resource.setrlimit(resource.RLIMIT_CORE, (0, 0))


def live_observed(records):
    result = []
    for saved in records:
        current = proc_record(saved["pid"])
        if current and current["state"] != "Z" and all(
                current[key] == saved[key] for key in ("start_ticks", "session")):
            result.append(current)
    return result


def stop_unreaped_new_session(child):
    """Fallback only before this freshly spawned direct child was reaped.

    Popen created a new session. Its unreaped PID cannot have been reused;
    matching SID=PGID=PID therefore identifies this new owned session. GDB
    inferiors may use a different PGID, so enumerate the SID, not only PGID.
    Do not poll/wait between validation and both signal passes: the unreaped
    root PID reserves this SID. pidfds prevent member PID reuse after lookup.
    """
    if child.returncode is not None:
        return []
    try:
        if os.getsid(child.pid) != child.pid or os.getpgid(child.pid) != child.pid:
            return []
    except ProcessLookupError:
        return []
    if not hasattr(os, "pidfd_open") or not hasattr(signal, "pidfd_send_signal"):
        raise RuntimeError("Identity-safe new-session cleanup requires native Linux pidfds")
    events = []
    for sig, delay in ((signal.SIGTERM, 1), (signal.SIGKILL, 1)):
        for path in Path("/proc").iterdir():
            if not path.name.isdigit():
                continue
            pid, descriptor = int(path.name), None
            try:
                descriptor = os.pidfd_open(pid)
                if os.getsid(pid) != child.pid:
                    continue
                signal.pidfd_send_signal(descriptor, sig, None, 0)
                events.append({"pid": pid, "session": child.pid, "signal": sig.name,
                    "ownership_basis": "DIRECT_UNREAPED_POPEN_CHILD_NEW_SID_MEMBER_WITH_PIDFD"})
            except ProcessLookupError:
                pass
            except PermissionError as error:
                # Unrelated /proc members may be inaccessible; no signal was
                # authorised for them. The owned SID is checked before signal.
                events.append({"pid": pid, "signal": sig.name, "inspection_error": str(error)})
            finally:
                if descriptor is not None:
                    os.close(descriptor)
        time.sleep(delay)
    return events


def capture_rss_stop_stack(records, sample, stop_reasons):
    """At most one second of owned-parent stack evidence before RSS cleanup.

    No delay is authorised when the host memory reserve is breached. This is
    not an epoch grace period: the original TERM/KILL cleanup always follows.
    """
    if ("OWNED_RSS_LIMIT" not in stop_reasons or
            sample["host_available_memory_bytes"] < LIMITS["minimum_host_available_memory_bytes"]):
        return []
    for saved in records:
        if saved["command"] != "gdb":
            continue
        descriptor = None
        try:
            descriptor = os.pidfd_open(saved["pid"])
            actual = proc_record(saved["pid"])
            if actual is None or any(actual[key] != saved[key] for key in ("start_ticks", "session")):
                continue
            signal.pidfd_send_signal(descriptor, signal.SIGINT, None, 0)
            event = {"pid": saved["pid"], "start_ticks": saved["start_ticks"],
                "session": saved["session"], "signal": "SIGINT",
                "reason": "BOUNDED_RSS_STOP_STACK_CAPTURE", "maximum_wait_seconds": 1}
            time.sleep(1)
            return [event]
        except (OSError, AttributeError) as error:
            return [{"pid": saved["pid"], "reason": "RSS_STOP_STACK_CAPTURE_UNAVAILABLE",
                     "capture_error": str(error), "maximum_wait_seconds": 0}]
        finally:
            if descriptor is not None:
                os.close(descriptor)
    return []


def profiler_command(plan, gdb=False):
    pea = network.command(plan["binary"], plan["configs"], plan["pinned_input_contract"])
    if gdb:
        pea = ["/usr/bin/gdb", "--quiet", "--batch", "--return-child-result",
            "-ex", "set pagination off", "-ex", "set confirm off",
            "-ex", "set startup-with-shell off", "-ex", "set print thread-events off",
            "-ex", "handle SIGINT stop print nopass",
            "-ex", "run", "-ex", "thread apply all bt 16", "--args", *pea]
    return ["/usr/bin/time", "-v", *pea]


def prepare(args):
    profile_flags = profile_environment(args.profile_environment)
    reference = json.loads(args.server_plan.read_text())
    plan, env = network.prepare(args.binary, args.audit, args.run_root, args.frozen_config,
        args.input_root, args.navigation, args.coordinate_overlay, args.coordinate_review,
        timeout=args.timeout_seconds, float_threads=12, omp_threads=12)
    root = Path(plan["run_root"])
    overlay = Path(plan["configs"][-2])
    value = json.loads(overlay.read_text())
    value["processing_options"]["epoch_control"]["end_epoch"] = EPOCH
    value["processing_options"]["gnss_general"]["zhang_pppar"]["checkpoint_capture_epochs"] = [EPOCH]
    stages.write_json(overlay, value)
    names = network.read_station_selection(args.frozen_config / stages.CONFIG_NAMES[1],
                                          json.loads(args.audit.read_text()))
    # The original 180 observed files/nav are freshly hashed by this call.
    plan["pinned_input_contract"] = network.pin_inputs(args.input_root, args.navigation, names, plan["configs"])
    current = {Path(row["path"]).name: row for group in ("observations", "navigation")
               for row in plan["pinned_input_contract"][group]}
    for group in ("observations", "navigation"):
        for row in reference["pinned_input_contract"][group]:
            other = current.get(Path(row["path"]).name)
            if other is None or (other["bytes"], other["sha256"]) != (row["bytes"], row["sha256"]):
                raise RuntimeError("Local physical input differs from actual server plan: " + row["path"])
    model_reference = json.loads(args.server_input_manifest.read_text())
    verified = []
    for row in model_reference["files"]:
        actual = current.get(Path(row["path"]).name) or network.identity(Path(row["path"]))
        if (actual["bytes"], actual["sha256"]) != (row["bytes"], row["sha256"]):
            raise RuntimeError("Local model/observation differs from server SHA reference: " + row["path"])
        verified.append(actual)
    env.update(profile_flags)
    plan.update(stage="LOCAL_180_FIRST_EPOCH_PROFILE_ONLY", state="PREPARED_NOT_RUN_LOCAL_APPROVAL_REQUIRED",
        end_epoch=EPOCH, expected_epoch_count=1, checkpoint_capture_epochs=[EPOCH],
        command=profiler_command(plan, args.gdb), diagnostic_only=True, scientific_acceptance=False,
        gdb_parent=args.gdb, profile_environment=profile_flags,
        rss_stop_stack_capture={"enabled": args.gdb, "maximum_wait_seconds": 1,
            "requires_host_memory_reserve": True, "no_resource_limit_extension": True},
        svd_profile_deadline_seconds=args.svd_profile_deadline_seconds,
        server_reference_plan=network.identity(args.server_plan),
        server_input_reference_manifest=network.identity(args.server_input_manifest),
        fresh_input_equivalence={"passed": True, "files": len(verified),
            "total_bytes": sum(row["bytes"] for row in verified), "identities": verified},
        local_memory_snapshot=linux_memory(), local_host=platform.node(),
        runtime_supervision=dict(LIMITS),
        resource_feasibility={"original_server_peak_rss_bytes": 26177884160,
            "original_server_peak_exceeds_local_memtotal": 26177884160 > linux_memory()["MemTotal"],
            "full_original_unbounded_180_run_allowed": False,
            "stage_dimensions_available": False, "swap_is_not_ram": True,
            "scope": "Only root-approved bounded profiling; RLIMIT_AS may stop before hotspot is reached"},
        environment={key: val for key, val in env.items() if key.startswith("ZHANG_") or key.endswith("_NUM_THREADS")})
    plan["source_identity"].append(network.identity(Path(__file__)))
    source_manifest = json.loads(args.source_manifest.read_text())
    source_archive = network.identity(args.source_archive)
    if source_manifest.get("binary_sha256") != plan["binary_sha256"] or (
            source_manifest.get("archive_sha256") != source_archive["sha256"]):
        raise RuntimeError("Full source archive is not bound to the actual profile producer")
    plan["producer_source_snapshot"] = {"archive": source_archive,
        "manifest": network.identity(args.source_manifest),
        "source_file_count": len(source_manifest.get("source_files", []))}
    if not plan["producer_source_snapshot"]["source_file_count"]:
        raise RuntimeError("Empty full-source manifest is not a producer snapshot")
    plan["source_identity"].extend((source_archive, network.identity(args.source_manifest)))
    plan["effective_input_preflight"] = {"status": "NOT_RUN", "required_before_solver_launch": True}
    stages.write_json(root / "stage_plan.json", plan)
    if args.preflight:
        network.preflight(plan, env)
    stages.write_json(root / "prepared_profile_plan.json", plan)
    return plan, env


def execute(plan, env, approval_path):
    root = Path(plan["run_root"])
    prepared = json.loads((root / "prepared_profile_plan.json").read_text())
    if plan != prepared or plan.get("state") != "PREPARED_NOT_RUN_LOCAL_APPROVAL_REQUIRED":
        raise RuntimeError("Only the original unused immutable local profile plan may execute")
    current = json.loads((root / "stage_plan.json").read_text())
    if current.get("state") != "PREPARED_NOT_RUN_LOCAL_APPROVAL_REQUIRED":
        raise RuntimeError("Profile root already attempted; never duplicate a numerical launch")
    approval = json.loads(approval_path.read_text())
    required = ("approved_for_bounded_local_profile", "build_passed", "mandatory_tests_passed",
                "hard_address_space_limit_reviewed", "no_estimator_or_diagnostic_config_change")
    if any(approval.get(key) is not True for key in required):
        raise RuntimeError("Bounded local diagnostic approval missing; no numerical process launched")
    if approval.get("binary_sha256") != plan["binary_sha256"] or approval.get("host") != platform.node():
        raise RuntimeError("Local approval binary/host mismatch")
    if approval.get("prepared_plan_sha256") != stages.file_sha(root / "prepared_profile_plan.json"):
        raise RuntimeError("Local approval is not bound to this exact prepared plan")
    if approval.get("runtime_supervision") != LIMITS:
        raise RuntimeError("Local approval resource budget mismatch")
    network.verify_pins(plan)
    network.preflight(plan, env)
    memory = linux_memory()
    if memory["MemAvailable"] < LIMITS["minimum_host_available_memory_bytes"]:
        raise RuntimeError("No local memory reserve; numerical process not launched")
    if plan["gdb_parent"] and not Path("/usr/bin/gdb").is_file():
        raise RuntimeError("Requested local gdb parent unavailable")
    return supervise(plan, env)


def supervise(plan, env):
    """Single owned local job; every exception after spawn performs cleanup."""
    root = Path(plan["run_root"])
    plan.update(state="PROFILE_RUNNING", started_utc=stamp())
    stages.write_json(root / "stage_plan.json", plan)
    started = time.monotonic()
    child, identity, guard_error = None, None, None
    observed, stop_reasons, stop_events = {}, [], []
    svd_seen_at, gdb_interrupt_at, phase_stop = None, None, False
    remaining = []
    try:
        with (root / "runner.log").open("wb") as output:
            child = subprocess.Popen(plan["command"], cwd=plan["working_directory"], env=env,
                stdout=output, stderr=subprocess.STDOUT, start_new_session=True, preexec_fn=hard_address_space_limit)
            identity = proc_record(child.pid)
            if identity is None:
                raise RuntimeError("Owned local profiler process identity unavailable")
            observed[(identity["pid"], identity["start_ticks"])] = identity
            stages.write_json(root / "profile_process.json", {"pid": child.pid, "identity": identity,
                "command": plan["command"], "hard_address_space_limit_bytes": LIMITS["maximum_address_space_bytes"]})
            with (root / "runtime_resource_samples.jsonl").open("w", buffering=1) as samples:
                while child.poll() is None:
                    processes = owned_tree(child.pid, identity["start_ticks"])
                    if not processes:
                        raise RuntimeError("Owned process monitoring unavailable while child is live")
                    observed.update({(row["pid"], row["start_ticks"]): row for row in processes})
                    mem = linux_memory()
                    disk = os.statvfs(root)
                    sample = {"elapsed_seconds": time.monotonic() - started, "sampled_utc": stamp(),
                        "processes": processes, "owned_tree_rss_bytes": sum(row["rss_bytes"] for row in processes),
                        "host_available_memory_bytes": mem["MemAvailable"],
                        "output_free_disk_bytes": disk.f_bavail * disk.f_frsize,
                        "swap_used_bytes": mem["SwapTotal"] - mem["SwapFree"]}
                    samples.write(json.dumps(sample) + "\n")
                    stop_reasons = violations(sample)
                    if sample["elapsed_seconds"] > plan["timeout_seconds"]:
                        stop_reasons.append("PROFILE_WALLTIME_LIMIT")
                    if stop_reasons:
                        stop_events.extend(capture_rss_stop_stack(processes, sample, stop_reasons))
                        stop_events.extend(stop_owned(list(observed.values())))
                        break
                    # Resources have priority over the bounded GDB stack capture.
                    text = (root / "runner.log").read_text(errors="replace")
                    if svd_seen_at is None and "ZHANG_FLOAT_PHASE_BEGIN phase=FLOAT_RANK_JACOBI_SVD " in text:
                        svd_seen_at = time.monotonic()
                    if svd_seen_at is not None and not phase_stop and (
                            time.monotonic() - svd_seen_at >= plan["svd_profile_deadline_seconds"]):
                        phase_stop = True
                        for saved in processes:
                            actual = proc_record(saved["pid"])
                            if saved["command"] == "gdb" and actual and all(
                                    actual[key] == saved[key] for key in ("start_ticks", "session")):
                                os.kill(saved["pid"], signal.SIGINT)
                                gdb_interrupt_at = time.monotonic()
                                stop_events.append({"pid": saved["pid"], "start_ticks": saved["start_ticks"],
                                                    "signal": "SIGINT", "reason": "BOUNDED_SVD_STACK_CAPTURE"})
                                break
                        if gdb_interrupt_at is None:
                            stop_events.extend(stop_owned(list(observed.values())))
                            break
                    if gdb_interrupt_at is not None and time.monotonic() - gdb_interrupt_at > 5:
                        stop_events.extend(stop_owned(list(observed.values())))
                        break
                    time.sleep(LIMITS["sample_interval_seconds"])
            child.wait(timeout=5)
    except BaseException as error:
        guard_error = type(error).__name__ + ": " + str(error)
    finally:
        # Preserve all previously observed identities, including children that
        # may have been reparented. Never reinterpret an absent tree as RSS=0.
        if child is not None:
            if identity is None:
                try:
                    stop_events.extend(stop_unreaped_new_session(child))
                except BaseException as error:
                    guard_error = guard_error or "Initial identity fallback cleanup: " + repr(error)
            if identity is not None:
                try:
                    latest = owned_tree(child.pid, identity["start_ticks"])
                    observed.update({(row["pid"], row["start_ticks"]): row for row in latest})
                except BaseException as error:
                    guard_error = guard_error or "Cleanup discovery: " + repr(error)
            try:
                if live_observed(list(observed.values())):
                    stop_events.extend(stop_owned(list(observed.values())))
            except BaseException as error:
                guard_error = guard_error or "Cleanup signalling: " + repr(error)
            try:
                child.wait(timeout=3)
            except BaseException as error:
                guard_error = guard_error or "Child not reaped: " + repr(error)
            try:
                remaining = live_observed(list(observed.values()))
                if remaining or child.poll() is None:
                    guard_error = guard_error or "Owned processes still live after bounded cleanup"
            except BaseException as error:
                guard_error = guard_error or "Cleanup verification unavailable: " + repr(error)
    try:
        text = (root / "runner.log").read_text(errors="replace")
    except BaseException as error:
        text = ""
        guard_error = guard_error or "Final log unavailable: " + repr(error)
    address_failure = any(marker in text for marker in ("std::bad_alloc", "Cannot allocate memory", "memory exhausted"))
    if guard_error:
        state = "PROFILE_GUARD_FAILED"
    elif stop_reasons or address_failure:
        state = "PROFILE_RESOURCE_STOP"
    elif phase_stop:
        state = "PROFILE_HOTSPOT_CAPTURE_STOP"
    else:
        state = "PROFILE_PROCESS_EXITED_NOT_SCIENTIFIC_ACCEPTANCE"
    plan.update(state=state, exit_code=child.returncode if child else None, elapsed_seconds=time.monotonic() - started,
        finished_utc=stamp(), resource_stop_reasons=stop_reasons,
        svd_phase_observed=svd_seen_at is not None, bounded_svd_profile_stop=phase_stop,
        possible_address_space_limit_failure=address_failure, identity_verified_signal_events=stop_events,
        guard_failure_reason=guard_error, owned_processes_still_live=remaining,
        child_reaped=child is not None and child.poll() is not None,
        qc_failure_not_inferred_from_resource_stop=True,
        runner_sha256=stages.file_sha(root / "runner.log") if (root / "runner.log").is_file() else None)
    stages.write_json(root / "stage_plan.json", plan)
    return plan


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("binary", "audit", "run-root", "frozen-config", "input-root", "navigation",
                 "coordinate-overlay", "coordinate-review", "server-plan", "server-input-manifest",
                 "source-archive", "source-manifest"):
        parser.add_argument("--" + name, type=Path)
    parser.add_argument("--prepared-plan", type=Path)
    parser.add_argument("--profile-environment", action="append", default=[])
    parser.add_argument("--timeout-seconds", type=float, default=600)
    parser.add_argument("--svd-profile-deadline-seconds", type=float, default=15)
    parser.add_argument("--gdb", action="store_true")
    parser.add_argument("--preflight", action="store_true")
    parser.add_argument("--execute", action="store_true")
    parser.add_argument("--approval", type=Path)
    args = parser.parse_args()
    if not math.isfinite(args.timeout_seconds) or args.timeout_seconds <= 0:
        parser.error("Finite positive diagnostic timeout required")
    if not 1 <= args.svd_profile_deadline_seconds <= 60:
        parser.error("SVD profile deadline must be 1..60 seconds")
    if args.execute:
        if args.approval is None or args.prepared_plan is None:
            parser.error("--execute needs existing --prepared-plan and exact local --approval")
        plan = json.loads(args.prepared_plan.resolve(strict=True).read_text())
        env = stages.environment("D")
        env.update(plan["environment"])
        plan = execute(plan, env, args.approval)
    else:
        required = ("binary", "audit", "run_root", "frozen_config", "input_root", "navigation",
                    "coordinate_overlay", "coordinate_review", "server_plan", "server_input_manifest",
                    "source_archive", "source_manifest")
        if any(getattr(args, key) is None for key in required) or args.prepared_plan is not None:
            parser.error("Preparation requires all explicit input/review/reference paths")
        plan, env = prepare(args)
    print(json.dumps({key: plan[key] for key in ("state", "run_root", "binary_sha256", "expected_epoch_count",
                                               "fresh_input_equivalence", "resource_feasibility")}, indent=2))


if __name__ == "__main__":
    main()
