#!/usr/bin/env python3
"""Real 12-station G/E/C FLOAT: legacy-full versus off versus unset.

Prepare-only by default. Execution requires an explicit same-binary approval.
Each independent cold process uses identical config/input/output paths and
runtime ID; only the narrow diagnostic environment variable changes. The
existing independent checkpoint reader compares x, the entire P, state order,
QC, transition/noise maps and every module section exactly. No scientific
quality or 180-station throughput acceptance is inferred from this experiment.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import re
import subprocess
import sys

import run_network180_float as network
import run_single_thread_cold_repeat as cold
import run_static_checkpoint_replay as replay
import run_stages as stages

VARIABLE = "ZHANG_FLOAT_PURE_OBS_RANK_DIAGNOSTIC"
ARMS = (("full", "full", "continuous"),
        ("off", "off", "interrupted_before_resume_cut"),
        ("unset", None, None))


def arm_environment(base, value):
    if value not in (None, "off", "full"):
        raise ValueError("Only reviewed full/off/unset comparison arms are permitted")
    result = dict(base)
    result.pop(VARIABLE, None)
    if value is not None:
        result[VARIABLE] = value
    return result


def diagnostic_evidence(log, value):
    text = Path(log).read_text(errors="replace")
    mode = "off" if value is None else value
    statuses = [line for line in text.splitlines()
                if "ZHANG_PURE_OBS_RANK_STATUS " in line and "mode=" + mode in line]
    starts = [line for line in text.splitlines()
              if "ZHANG_FLOAT_PHASE_BEGIN phase=FLOAT_RANK_JACOBI_SVD " in line]
    ends = [line for line in text.splitlines()
            if "ZHANG_FLOAT_PHASE_END phase=FLOAT_RANK_JACOBI_SVD " in line]
    rank_values = [line for line in text.splitlines() if "ZHANG_PURE_OBS_RANK time=" in line]
    if mode == "full":
        if len(starts) != 1 or len(ends) != 1 or len(rank_values) != 1 or not any(
                "status=COMPLETED " in line and "rank=EVALUATED" in line for line in statuses):
            raise RuntimeError("Full arm did not actually complete the original diagnostic SVD once")
    else:
        if starts or ends or rank_values or not any(
                "status=DISABLED " in line and "rank=NOT_EVALUATED" in line for line in statuses):
            raise RuntimeError("Off/default arm ran rank diagnostic or lacks explicit NOT_EVALUATED evidence")
    if not re.search(r"Processed epoch #1\b", text):
        raise RuntimeError("Arm lacks actual committed/processed first-epoch evidence")
    return {"mode": mode, "environment_value": value, "svd_starts": len(starts),
            "svd_ends": len(ends), "status_lines": statuses, "rank_lines": rank_values}


def verify_checkpoint_comparison(report, expected_receivers):
    required = {"core.x", "core.P", "core.dx", "core.kfIndexMap", "core.stateTransitionMap",
                "core.procNoiseMap", "core.initNoiseMap", "core.sigmaPass", "core.chiQCPass",
                "core.qc", "sections.inventory", "core.authoritative_spp_offset_all_GEC_present",
                "core.static_positions_all_XYZ_positive_variance", "core.no_receiver_position_rate"}
    checks = report.get("checks", {})
    if report.get("passed") is not True or not required.issubset(checks) or any(
            checks[key] is not True for key in required):
        raise RuntimeError("Independent reader did not pass all required state/full-P/order/QC/module checks")
    if set(report.get("coordinate_components_by_receiver", {})) != set(expected_receivers):
        raise RuntimeError("Committed coordinate states do not match all original selected receivers")
    numeric = report.get("numeric_P", {})
    if not numeric.get("all_finite") or numeric.get("exact_equal") is not True or \
            numeric.get("maximum_absolute_difference_native") != 0:
        raise RuntimeError("Full covariance comparison must be finite and exactly identical")
    sections = report.get("sections", [])
    if not sections or any(row.get("entire_payload_byte_identical") is not True for row in sections):
        raise RuntimeError("Missing or nonidentical checkpoint module payload")
    return True


def prepare(binary, audit, root, coordinate_overlay, coordinate_review, verifier, timeout=1800):
    # This reuses the reviewed static XYZ(1 m prior, Q=0) 12-site model,
    # not a fixed-coordinate network or a migrated GPS checkpoint.
    template = {"binary": str(binary.resolve(strict=True)), "verifier": str(verifier.resolve(strict=True)),
        "external_review_files": [{"path": str(coordinate_overlay.resolve(strict=True))},
                                  {"path": str(coordinate_review.resolve(strict=True))},
                                  {"path": str(audit.resolve(strict=True))}]}
    plan, env = cold.prepare(template, root)
    if not 0 < timeout <= 3600:
        raise ValueError("Bounded small-network diagnostic must have 0 < timeout <= 3600 seconds")
    env.update(ZHANG_FLOAT_PROFILE="1")
    env.pop(VARIABLE, None)
    plan.update(diagnostic="PURE_OBSERVATION_RANK_FULL_OFF_UNSET_EXACT_EQUIVALENCE",
        state="RANK_DIAGNOSTIC_PREPARED_NOT_RUN", timeout_seconds=timeout,
        arms=[{"label": label, "value": value} for label, value, _ in ARMS],
        environment={key: value for key, value in env.items() if key.startswith("ZHANG_") or
                     key.endswith("_NUM_THREADS") or key.startswith("OMP_")},
        verifier_sha256=stages.file_sha(verifier),
        comparison_policy="Strict entrywise exact x/full P/state order/QC/F/Q/priors/module payload; exact product/full-Q/certificate/COMMIT bytes; no normalization or tolerances",
        scientific_acceptance=False)
    plan["source_identity"].append(network.identity(Path(__file__).resolve()))
    stages.write_json(root / "replay_plan.json", plan)
    stages.write_json(root / "rank_equivalence_prepared_plan.json", plan)
    return plan, env


def execute(plan, env, approval):
    if sys.platform != "linux":
        raise RuntimeError("Actual controlled comparison must execute in native Linux")
    root = replay.owner_guard(Path(plan["run_root"]))
    if plan.get("state") != "RANK_DIAGNOSTIC_PREPARED_NOT_RUN" or (root / "rank_execution_started.json").exists():
        raise RuntimeError("This uniquely owned comparison has already been attempted")
    review = replay.execution_gate(plan, approval)
    if review.get("approved_for_small_network_legacy_full_rank_diagnostic") is not True:
        raise RuntimeError("Legacy full SVD needs explicit approval for this bounded 12-site experiment")
    if review.get("prepared_plan_sha256") != stages.file_sha(root / "rank_equivalence_prepared_plan.json"):
        raise RuntimeError("Approval must pin the exact prepared comparison plan")
    replay.verify_plan(plan)
    if stages.file_sha(plan["verifier"]) != plan["verifier_sha256"]:
        raise RuntimeError("Independent checkpoint reader changed")
    preflight = stages.configuration_preflight(Path(plan["binary"]), plan["configs"],
        plan["pinned_input_contract"], root, env)
    if not preflight["passed"]:
        raise RuntimeError("Actual production configuration parser rejected exact 12-site input")
    stages.write_json(root / "rank_execution_started.json", {"started_utc": replay.stamp(),
        "binary_sha256": plan["binary_sha256"], "approval": review})
    plan.update(state="RANK_DIAGNOSTIC_COMPARISON_RUNNING", configuration_preflight=preflight)
    stages.write_json(root / "replay_plan.json", plan)
    arms = {}
    for label, value, archive in ARMS:
        run = replay.run_phase(plan, arm_environment(env, value), "rank_" + label)
        evidence = diagnostic_evidence(run["log"], value)
        work = replay.archive_owned_outputs(root, archive) if archive else root / "work"
        threads = cold.thread_evidence(work / "outputs")
        products = cold.one_epoch_products(work / "outputs")
        checkpoint = replay.checkpoint_at(work / "checkpoints", cold.EPOCH)
        if checkpoint is None:
            raise RuntimeError("Arm did not publish a verified single-epoch checkpoint")
        arms[label] = {"run": run, "diagnostic_evidence": evidence, "thread_evidence": threads,
            "checkpoint": str(checkpoint), "products_sha256": products}
        if archive:
            for name in ("outputs", "checkpoints"):
                (root / "work" / name).mkdir()
    expected_receivers = {row["station"] for row in plan["pinned_input_contract"]["observations"]}
    comparisons = {}
    for label in ("off", "unset"):
        path = root / ("full_versus_" + label + "_full_checkpoint.json")
        with (root / ("full_versus_" + label + "_reader.log")).open("w") as stream:
            checked = subprocess.run([plan["verifier"], str(Path(arms["full"]["checkpoint"]) / "checkpoint.bundle"),
                str(Path(arms[label]["checkpoint"]) / "checkpoint.bundle"), plan["binary_sha256"],
                replay.MODE, str(path)], stdout=stream, stderr=subprocess.STDOUT, timeout=120)
        if checked.returncode:
            raise RuntimeError("Full versus " + label + " independent checkpoint comparison failed")
        report = json.loads(path.read_text())
        verify_checkpoint_comparison(report, expected_receivers)
        if arms["full"]["products_sha256"] != arms[label]["products_sha256"]:
            raise RuntimeError("Full versus " + label + " complete product/Q/certificate/COMMIT differs")
        comparisons[label] = {"checkpoint_full_wire_passed": True, "products_and_full_q_exact": True,
                              "report": str(path)}
    replay.verify_plan(plan)
    if stages.file_sha(plan["verifier"]) != plan["verifier_sha256"]:
        raise RuntimeError("Independent checkpoint reader changed during execution")
    result = {"passed": True, "experiment": plan["diagnostic"], "arms": arms, "comparisons": comparisons,
        "actual_epoch": cold.EPOCH, "original_selected_station_count": 12, "systems": ["GPS", "GAL", "BDS"],
        "ratio": 2, "all_nis": False, "independent_epoch_ar_configured": True, "ar_executed": False,
        "numeric_normalization": False, "scientific_acceptance": False,
        "scope": "Three independent cold single-thread FLOAT epochs only; not a full-hour, AR or scientific accuracy acceptance"}
    stages.write_json(root / "rank_equivalence_result.json", result)
    plan.update(state="RANK_DIAGNOSTIC_EXACT_EQUIVALENCE_PASSED", finished_utc=replay.stamp())
    stages.write_json(root / "replay_plan.json", plan)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("binary", "audit", "run-root", "coordinate-overlay", "coordinate-review", "verifier"):
        parser.add_argument("--" + name, type=Path)
    parser.add_argument("--prepared-plan", type=Path)
    parser.add_argument("--timeout-seconds", type=float, default=1800)
    parser.add_argument("--build-approval", type=Path)
    parser.add_argument("--execute", action="store_true")
    args = parser.parse_args()
    if args.execute and args.build_approval is None:
        parser.error("--execute requires --build-approval")
    if args.prepared_plan:
        plan = json.loads(args.prepared_plan.read_text())
        env = stages.environment("D")
        env.update(plan["environment"])
        env.pop(VARIABLE, None)
    else:
        names = ("binary", "audit", "run_root", "coordinate_overlay", "coordinate_review", "verifier")
        if any(getattr(args, key) is None for key in names):
            parser.error("Preparation requires --binary --audit --run-root --coordinate-overlay --coordinate-review --verifier")
        plan, env = prepare(*(getattr(args, key) for key in names), timeout=args.timeout_seconds)
    if args.execute:
        try:
            result = execute(plan, env, args.build_approval)
        except Exception as error:
            plan.update(state="RANK_DIAGNOSTIC_COMPARISON_FAILED", failure_reason=str(error),
                finished_utc=replay.stamp(), scientific_acceptance=False)
            stages.write_json(Path(plan["run_root"]) / "replay_plan.json", plan)
            raise
        print(json.dumps(result, indent=2), flush=True)
    else:
        print(json.dumps(plan, indent=2), flush=True)


if __name__ == "__main__":
    main()
