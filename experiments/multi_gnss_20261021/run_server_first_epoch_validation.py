#!/usr/bin/env python3
"""Prepare or explicitly execute one bounded original-180 FLOAT epoch.

This does not resume the two-hour experiment. The reviewed frozen-v2 ownership
helpers supervise one private producer, with a 900 s/40 GiB RSS guard. A default
off pure-observation rank diagnostic is the intended code change; all numerical
models, input identities and FLOAT thresholds remain those of the stopped run.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import platform
import re
import sys

import run_server_first_epoch_profile as server

EPOCH = server.EPOCH
PREPARED_STATE = "PREPARED_NOT_RUN_SERVER_VALIDATION_APPROVAL_REQUIRED"
ALLOWED_POLICY_ENV = {"ZHANG_FLOAT_PROFILE": "1", "ZHANG_FLOAT_PURE_OBS_RANK_DIAGNOSTIC": "off"}


def matrix_environment(matrix_threads, memory_mode):
    if matrix_threads is not None and (type(matrix_threads) is not int or not 1 <= matrix_threads <= 28):
        raise ValueError("Explicit matrix threads must be an integer in 1..28")
    if memory_mode not in ("serial", "parallel"):
        raise ValueError("Memory mode must be serial or parallel")
    result = {"ZHANG_FLOAT_MEMORY_MODE": memory_mode}
    if matrix_threads is not None:
        result["ZHANG_FLOAT_MATRIX_THREADS"] = str(matrix_threads)
    return result


def verify_source(args, stages, network):
    manifest = json.loads(args.source_manifest.read_text())
    if (manifest["binary_sha256"] != stages.file_sha(args.binary) or
            manifest["archive_sha256"] != stages.file_sha(args.source_archive)):
        raise RuntimeError("Actual producer/source archive SHA binding failed")
    expected = {row["path"]: row["sha256"] for row in manifest["source_files"]}
    for name in ("run_stages.py", "run_network180_float.py", "run_local_first_epoch_profile.py",
                 "audit_network180_float_completion.py", "run_static_checkpoint_replay.py"):
        relative = "experiments/multi_gnss_20261021/" + name
        if stages.file_sha(args.frozen_drivers / name) != expected.get(relative):
            raise RuntimeError("Frozen helper differs from the producer source archive: " + name)
    for path in (Path(__file__), Path(server.__file__)):
        relative = "experiments/multi_gnss_20261021/" + path.name
        if stages.file_sha(path) != expected.get(relative):
            raise RuntimeError("Server validation source is not in the producer source archive: " + path.name)
    runtime = [network.identity(path) for path in sorted((args.binary.parent / "lib").rglob("*"))
               if path.is_file()]
    if not runtime:
        raise RuntimeError("Private producer runtime library identities are required")
    return {"archive": network.identity(args.source_archive), "manifest": network.identity(args.source_manifest),
            "source_file_count": len(expected), "runtime_files": runtime}


def verify_reference_inputs(plan, reference, input_manifest, network):
    current = {Path(row["path"]).name: row for group in ("observations", "navigation")
               for row in plan["pinned_input_contract"][group]}
    for group in ("observations", "navigation"):
        old = reference["pinned_input_contract"][group]
        if len(old) != len(plan["pinned_input_contract"][group]):
            raise RuntimeError("Original and validation physical input counts differ")
        for row in old:
            actual = current.get(Path(row["path"]).name)
            if actual is None or (actual["bytes"], actual["sha256"]) != (row["bytes"], row["sha256"]):
                raise RuntimeError("Physical input differs from original server run: " + row["path"])
    verified = []
    for row in json.loads(input_manifest.read_text())["files"]:
        actual = current.get(Path(row["path"]).name) or network.identity(Path(row["path"]))
        if (actual["bytes"], actual["sha256"]) != (row["bytes"], row["sha256"]):
            raise RuntimeError("Model/input differs from server reference manifest: " + row["path"])
        verified.append(actual)
    return {"passed": True, "files": len(verified), "total_bytes": sum(row["bytes"] for row in verified),
            "identities": verified}


def prepare(args, stages, network, profile):
    matrix_policy = matrix_environment(args.matrix_threads, args.memory_mode)
    producer = verify_source(args, stages, network)
    reference = json.loads(args.reference_plan.read_text())
    if not (reference["ratio"] == 2 and reference["all_nis"] is False and
            reference["ar_execution"] is False and reference["epoch_only_ar"] is True):
        raise RuntimeError("Original reference is not the reviewed ratio2/NIS-off independent FLOAT run")
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
    plan["configuration_equivalence"] = server.configuration_equivalence(plan, reference, network)
    equivalent_inputs = verify_reference_inputs(plan, reference, args.input_manifest, network)
    for key, value in reference["environment"].items():
        if env.get(key) != value:
            raise RuntimeError("Original server numerical environment changed: " + key)
    # The sole non-timing environment change is the reviewed narrow off policy;
    # full_rank.output_diagnostics and estimator gates remain byte-identical.
    env.update(ALLOWED_POLICY_ENV)
    env.update(matrix_policy)
    plan.update(stage="SERVER_180_FIRST_EPOCH_VALIDATION_ONLY", state=PREPARED_STATE,
        expected_epoch_count=1, end_epoch=EPOCH, checkpoint_capture_epochs=[EPOCH],
        diagnostic_only=True, scientific_acceptance=False, local_host=platform.node(),
        runtime_supervision=dict(server.SERVER_LIMITS), timeout_seconds=900,
        svd_profile_deadline_seconds=15, gdb_parent=False,
        intended_rank_policy="NARROW_OFF_EXPLICIT_NO_ESTIMATOR_CHANGE",
        matrix_thread_budget=args.matrix_threads, memory_mode=args.memory_mode,
        fresh_input_equivalence=equivalent_inputs,
        original_server_reference_plan=network.identity(args.reference_plan),
        producer_source_snapshot=producer,
        environment={key: val for key, val in env.items() if key.startswith("ZHANG_") or key.endswith("_NUM_THREADS")})
    plan["source_identity"].extend([network.identity(args.frozen_drivers / name) for name in
        ("audit_network180_float_completion.py", "run_static_checkpoint_replay.py")])
    # verify_pins() rechecks source_identity immediately before and after the
    # numerical run. Pin every model file, not just observation/nav/config: the
    # fresh-equivalence report alone is a historical preparation-time statement.
    plan["source_identity"].extend([*equivalent_inputs["identities"],
        network.identity(args.reference_plan), network.identity(args.input_manifest)])
    plan["source_identity"].extend([network.identity(Path(profile.__file__)), network.identity(Path(__file__)),
        network.identity(Path(server.__file__)), network.identity(args.source_archive),
        network.identity(args.source_manifest), *producer["runtime_files"]])
    plan["command"] = ["/usr/bin/time", "-v", *network.command(plan["binary"], plan["configs"], plan["pinned_input_contract"])]
    stages.write_json(root / "stage_plan.json", plan)
    network.preflight(plan, env)
    stages.write_json(root / "prepared_server_validation_plan.json", plan)
    return plan


def output_audit(plan, stages):
    """Use complete product/full-Q validation; checkpoint decode remains separate."""
    # Import only when requested, and require the exact frozen audit driver.
    import audit_network180_float_completion as completion
    root = Path(plan["run_root"])
    report = {"scope": "Single real 180-station cold epoch only; not two-hour or accuracy acceptance",
              "passed": False, "checkpoint_full_kf_decode_claimed": False}
    try:
        if plan.get("exit_code") != 0 or plan.get("guard_failure_reason") or plan.get("resource_stop_reasons"):
            raise RuntimeError("Producer did not finish cleanly within the reviewed limits")
        if plan.get("owned_processes_still_live") or not plan.get("child_reaped"):
            raise RuntimeError("Owned process cleanup is not verified")
        text = (root / "runner.log").read_text(errors="replace")
        if any(marker in text for marker in (
                "ZHANG_FLOAT_PHASE_BEGIN phase=FLOAT_RANK_JACOBI_SVD ",
                "ZHANG_FLOAT_PHASE_BEGIN phase=FLOAT_PURE_OBSERVATION_COPY ")):
            raise RuntimeError("Narrow-off producer unexpectedly copied or decomposed the pure observation rank matrix")
        if "ZHANG_PURE_OBS_RANK_STATUS mode=off status=DISABLED" not in text or "rank=NOT_EVALUATED" not in text:
            raise RuntimeError("Explicit disabled/not-evaluated pure-rank status is absent")
        if "ZHANG_FLOAT_PHASE_END phase=FLOAT_TOTAL " not in text:
            raise RuntimeError("Completed FLOAT_TOTAL phase is absent")
        if plan.get("matrix_thread_budget") is not None or plan.get("memory_mode") is not None:
            policies = []
            for path in (root / "outputs").glob("Network-*.TRACE"):
                policies.extend(line for line in path.read_text(errors="replace").splitlines()
                                if line.startswith("ZHANG_FLOAT_THREADS "))
            if len(policies) != 1:
                raise RuntimeError("Exactly one observed FLOAT thread-policy record is required")
            observed = dict(re.findall(r"(\w+)=(\S+)", policies[0]))
            budget = plan.get("matrix_thread_budget")
            if ((budget is not None and observed.get("blas") != str(budget)) or
                    observed.get("requested") != "28" or observed.get("omp") != "28" or
                    observed.get("memory_mode") != plan.get("memory_mode")):
                raise RuntimeError("Observed matrix/outer/memory thread policy differs from the pinned plan")
            report["observed_thread_policy"] = observed
        import run_network180_float as network
        failures = network.CurrentRunLogTail(root).final_failures()
        if failures:
            raise RuntimeError("Authoritative FLOAT transaction failed: " + json.dumps(failures))
        directories = list((root / "outputs/zhang_internal_products.csv.epochs").iterdir())
        if len(directories) != 1 or not directories[0].is_dir():
            raise RuntimeError("Expected exactly one atomic product epoch")
        epoch = completion.audit_epoch(directories[0])
        if epoch["epoch"] != EPOCH:
            raise RuntimeError("Atomic product epoch is not the single cold first epoch")
        checkpoint_dirs = [path for path in (root / "checkpoints").iterdir() if path.is_dir()]
        if len(checkpoint_dirs) != 1:
            raise RuntimeError("Expected exactly one captured checkpoint directory")
        report.update(passed=True, epoch=epoch, checkpoint_directory=str(checkpoint_dirs[0]),
            checkpoint_files=[{"path": str(path), "bytes": path.stat().st_size,
                               "sha256": stages.file_sha(path)}
                              for path in sorted(checkpoint_dirs[0].rglob("*")) if path.is_file()],
            accuracy_acceptance=False)
    except Exception as error:
        report["failure_reason"] = type(error).__name__ + ": " + str(error)
    stages.write_json(root / "single_epoch_output_audit.json", report)
    return report


def execute(plan, approval_path, stages, network, profile):
    root = Path(plan["run_root"])
    prepared_path = root / "prepared_server_validation_plan.json"
    if (plan != json.loads(prepared_path.read_text()) or
            json.loads((root / "stage_plan.json").read_text())["state"] != PREPARED_STATE):
        raise RuntimeError("Only the exact unused server validation plan may execute once")
    approval = json.loads(approval_path.read_text())
    required = ("approved_for_bounded_server_validation", "build_passed", "mandatory_tests_passed",
                "diagnostic_default_off_reviewed", "no_estimator_or_model_config_change")
    if any(approval.get(key) is not True for key in required):
        raise RuntimeError("Missing bounded server validation approval")
    if (approval.get("host") != platform.node() or approval.get("binary_sha256") != plan["binary_sha256"] or
            approval.get("prepared_plan_sha256") != stages.file_sha(prepared_path) or
            approval.get("runtime_supervision") != server.SERVER_LIMITS):
        raise RuntimeError("Server host/binary/plan/resource approval mismatch")
    network.verify_pins(plan)
    env = stages.environment("D")
    env.update(plan["environment"])
    network.preflight(plan, env)
    if profile.linux_memory()["MemAvailable"] < server.SERVER_LIMITS["minimum_host_available_memory_bytes"]:
        raise RuntimeError("Reviewed host memory reserve unavailable before launch")
    result = profile.supervise(plan, env)
    network.verify_pins(result)
    report = output_audit(result, stages)
    result.update(state="SINGLE_EPOCH_NUMERICAL_SMOKE_PASSED_REQUIRES_CHECKPOINT_DECODE_AND_ACCURACY_REVIEW"
                  if report["passed"] else "SINGLE_EPOCH_VALIDATION_NOT_PASSED",
                  single_epoch_output_audit=report, scientific_acceptance=False)
    stages.write_json(root / "stage_plan.json", result)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("frozen-drivers", "binary", "audit", "run-root", "frozen-config", "input-root", "navigation",
                 "coordinate-overlay", "coordinate-review", "reference-plan", "input-manifest", "source-archive",
                 "source-manifest", "prepared-plan", "approval"):
        parser.add_argument("--" + name, type=Path)
    parser.add_argument("--execute", action="store_true")
    parser.add_argument("--matrix-threads", type=int)
    parser.add_argument("--memory-mode", choices=("serial", "parallel"), default="parallel")
    args = parser.parse_args()
    if args.frozen_drivers is None:
        parser.error("Explicit reviewed frozen-v2 helper directory required")
    stages, network, profile = server.load_helpers(args.frozen_drivers)
    if args.execute:
        if args.prepared_plan is None or args.approval is None:
            parser.error("Execution requires an exact prepared plan and root approval")
        plan = execute(json.loads(args.prepared_plan.read_text()), args.approval, stages, network, profile)
    else:
        if any(getattr(args, name) is None for name in ("binary", "audit", "run_root", "frozen_config",
                "input_root", "navigation", "coordinate_overlay", "coordinate_review", "reference_plan",
                "input_manifest", "source_archive", "source_manifest")):
            parser.error("Preparation requires all explicit producer/input/reference paths")
        plan = prepare(args, stages, network, profile)
    print(json.dumps({key: plan.get(key) for key in ("state", "run_root", "binary_sha256", "expected_epoch_count",
        "resource_stop_reasons", "guard_failure_reason", "single_epoch_output_audit")}, indent=2))


if __name__ == "__main__":
    main()
