#!/usr/bin/env python3
"""Read-only final-binary cold FLOAT pilot replay evidence, not AR quality."""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import re

from audit_stage_results import StageAuditor, gpst, sha256, unavailable
from audit_candidate_replays import compare_cold_seeds
from run_clock_trials import resource_statistics

FINAL_SHA = "af149dfc28215bf25dbb932f691fb3a70dff4df0760e1490b27390a01770a1db"
START = 1405209600.0


def qc_transactions(root):
    commits, rollbacks, terminal, rejected, evidence = {}, [], [], [], []
    for path in sorted((root / "outputs").glob("Network*.TRACE")):
        with path.open(errors="replace") as stream:
            for number, line in enumerate(stream, 1):
                if "ZHANG_KF_MEASUREMENT_TRANSACTION" in line and "|/PPP|" in line:
                    match = re.search(r"time=(.*?) boundary=B13 status=(\S+)", line)
                    if match:
                        row = {"path": str(path), "line": number, "epoch": match[1], "text": line.strip()}
                        if match[2] == "COMMITTED":
                            commits[match[1]] = row
                        elif match[2] == "ROLLED_BACK":
                            rollbacks.append(row)
                        evidence.append(row)
                if "ZHANG_KF_FINAL_POSTFIT_QC" in line:
                    terminal.append({"line": number, "text": line.strip()})
                if "PPP AUTHORITATIVE MEASUREMENT TRANSACTION REJECTED" in line:
                    rejected.append({"line": number, "text": line.strip()})
    return {"committed_unique_epoch_count": len(commits), "committed_epochs": sorted(commits),
        "rollback_count": len(rollbacks), "rollbacks": rollbacks,
        "terminal_reconciliation_records": terminal, "authoritative_rejection_records": rejected,
        "commit_records": list(commits.values())}


def clock_lifecycle(root):
    """Stream lifecycle events, retaining small per-receiver offsets and cold seeds."""
    path = root / "clock_lifecycle_events.jsonl"
    if not path.is_file():
        return unavailable("No structured clock lifecycle events file")
    phases, sources, status_counts, per_epoch_source = Counter(), Counter(), Counter(), defaultdict(Counter)
    last_offsets, cold_receivers, cold_refs, cold_external = {}, {}, {}, {}
    adjustment_counts, adjustments, failures, samples = Counter(), defaultdict(list), [], []
    sat_residuals, cold_satellites = [], {}
    formula_max = offset_max = 0.0
    checked = matching_offsets = chart_matches = 0
    explicit_failed_spp, invalid = [], []
    with path.open() as stream:
        for number, line in enumerate(stream, 1):
            try:
                row = json.loads(line)
                event = row["event"]
                if event.get("schema") != "ZHANG_CLOCK_LIFECYCLE_V1":
                    raise ValueError("Unknown lifecycle schema")
            except (ValueError, TypeError, KeyError) as error:
                invalid.append({"jsonl_line": number, "reason": str(error)})
                continue
            phase = event.get("phase")
            phases[phase] += 1
            epoch, receiver = event.get("epoch"), event.get("receiver")
            system = event.get("system", event.get("spp_system"))
            identity = (receiver, system)
            first_epoch = gpst(epoch) == START
            if phase == "SPP_RESULT":
                status = event.get("spp_status", "UNAVAILABLE")
                status_counts[(epoch, event.get("spp_system"), status)] += 1
                if status != "SINGLE" and len(explicit_failed_spp) < 150:
                    explicit_failed_spp.append({"epoch": epoch, "receiver": receiver,
                        "selected_system": event.get("spp_system"), "status": status, "runner_line": row["source_line"]})
                for obs in event.get("observations", []):
                    if obs.get("spp_valid") is True:
                        source = obs.get("source", "UNAVAILABLE")
                        sources[source] += 1
                        per_epoch_source[(epoch, obs.get("system"))][source] += 1
            if phase == "SPP_OFFSET_SAVE_AFTER":
                offset = event.get("saved_system_offset_m", event.get("saved_offset_m"))
                last_offsets[identity] = {"offset": offset, "epoch": epoch,
                    "external_chart": event.get("external_chart"), "runner_line": row["source_line"]}
                if first_epoch and event.get("spp_system_clock_m") is not None:
                    cold_external[identity] = event["spp_system_clock_m"]
            if phase == "REC_CLOCK_ADJUST_AFTER":
                applied = event.get("applied_adjustment_m")
                if applied is not None:
                    adjustments[system].append(applied)
                    adjustment_counts[system] += 1
                spp = event.get("spp_system_clock_m", event.get("spp_clock_m"))
                offset = event.get("saved_system_offset_m", event.get("saved_offset_m"))
                prediction = event.get("predicted_clock_m")
                unrounded = event.get("unrounded_adjustment_m")
                if None not in (spp, offset, prediction, unrounded):
                    checked += 1
                    residual = spp + offset - prediction - unrounded
                    formula_max = max(formula_max, abs(residual))
                    old = last_offsets.get(identity)
                    if old is not None and old["offset"] is not None:
                        difference = offset - old["offset"]
                        offset_max = max(offset_max, abs(difference))
                        matching_offsets += abs(difference) <= 1e-9
                        chart_matches += event.get("external_chart") is not None and event.get("external_chart") == old["external_chart"]
                    else:
                        failures.append({"epoch": epoch, "receiver": receiver, "system": system,
                            "runner_line": row["source_line"], "reason": "Adjustment has no earlier observed committed offset save"})
                    if abs(residual) > 1e-8:
                        failures.append({"epoch": epoch, "receiver": receiver, "system": system,
                            "runner_line": row["source_line"], "formula_residual_m": residual})
                    if receiver == "NNOR" and (gpst(epoch) in (START + 30, START + 1800)):
                        samples.append({"epoch": epoch, "system": system, "receiver": receiver,
                            "runner_line": row["source_line"], "external_spp_clock_m": spp,
                            "saved_external_to_kf_offset_m": offset, "predicted_kf_clock_m": prediction,
                            "recomputed_adjustment_m": spp + offset - prediction, "applied_adjustment_m": applied})
            if first_epoch and phase == "COLD_REFERENCE_CLOCK_ZERO":
                cold_refs.setdefault(system, {"reference": event.get("model_reference"),
                    "model_clock_m": event.get("model_clock_m"), "selected_spp_system": event.get("spp_system"),
                    "selected_spp_clock_m": event.get("spp_clock_m"), "selected_spp_status": event.get("spp_status"),
                    "runner_line": row["source_line"]})
            if first_epoch and phase == "COLD_RECEIVER_CLOCK_SEED":
                cold_receivers.setdefault(identity, {"init_clock_m": event.get("init_clock_m"),
                    "selected_spp_system": event.get("spp_system"), "selected_spp_clock_m": event.get("spp_clock_m"),
                    "selected_spp_status": event.get("spp_status"),
                    "runner_line": row["source_line"]})
            if first_epoch and phase == "COLD_SATELLITE_CLOCK_SEED":
                cold_satellites.setdefault((event.get("system"), event.get("satellite")), {"seed_m": event.get("init_clock_m"),
                    "runner_line": row["source_line"]})
    seeds = []
    for (receiver, system), event in sorted(cold_receivers.items()):
        external = cold_external.get((receiver, system))
        if external is None and event["selected_spp_system"] == system and event["selected_spp_status"] == "SINGLE":
            external = event["selected_spp_clock_m"]
        reference = cold_refs.get(system, {})
        known_ref = reference.get("selected_spp_clock_m") if reference.get("selected_spp_system") == system and reference.get("selected_spp_status") == "SINGLE" else None
        translation = event["init_clock_m"] - external if external is not None else None
        seeds.append({"receiver": receiver, "system": system, "runner_line": event["runner_line"],
            "seed_m": event["init_clock_m"], "actual_external_system_spp_m": external,
            "external_coordinate_status": "AVAILABLE" if external is not None else "UNAVAILABLE_FAILED_OR_UNLOGGED_SPP_NOT_ZERO",
            "inferred_translation_m": translation, "directly_observed_reference_spp_m": known_ref,
            "translation_plus_reference_spp_m": translation + known_ref if None not in (translation, known_ref) else None})
    by_system = {}
    for system in {row["system"] for row in seeds}:
        values = [row["inferred_translation_m"] for row in seeds if row["system"] == system and row["inferred_translation_m"] is not None]
        residuals = [row["translation_plus_reference_spp_m"] for row in seeds if row["system"] == system and row["translation_plus_reference_spp_m"] is not None]
        by_system[system] = {"receiver_seed_count": len(values), "minimum_common_translation_m": min(values) if values else None,
            "maximum_common_translation_m": max(values) if values else None,
            "translation_spread_m": max(values) - min(values) if values else None,
            "max_abs_direct_reference_check_residual_m": max(map(abs, residuals)) if residuals else None,
            "direct_reference_check_scope": "UNAVAILABLE when reference's per-system external SPP coordinate was not logged; selected GPS is not substituted for GAL/BDS"}
    return {"status": "AVAILABLE", "events_sha256": sha256(path), "phase_counts": dict(phases),
        "invalid_records": invalid, "valid_spp_observation_clock_sources": dict(sources),
        "no_kalman_clock_used_for_valid_spp_observations": sources["KALMAN"] == 0 and bool(sources),
        "source_counts_by_epoch_observation_system": [{"epoch": key[0], "system": key[1], "sources": dict(value)} for key, value in sorted(per_epoch_source.items())],
        "explicit_spp_result_counts": [{"epoch": key[0], "selected_system": key[1], "status": key[2], "count": value} for key, value in sorted(status_counts.items())],
        "non_single_spp_samples": explicit_failed_spp,
        "offset_checks": {"formula_checked_adjustments": checked, "matching_previous_saved_offset_count": matching_offsets,
            "matching_external_descriptor_count": chart_matches, "maximum_formula_residual_m": formula_max,
            "maximum_previous_saved_offset_difference_m": offset_max, "failures": failures,
            "descriptor_scope": "Logged source descriptor equality, not an independent complete chart-transport certificate"},
        "actual_adjustment_summary_by_system": {system: {"sample_count": len(values), "maximum_abs_m": max(map(abs, values))} for system, values in adjustments.items()},
        "nnor_boundary_samples": samples, "cold_reference_zero_events": cold_refs,
        "cold_receiver_seeds": seeds, "cold_receiver_translation_by_system": by_system,
        "cold_satellite_seed_count": len(cold_satellites),
        "cold_satellite_seed_scope": "Direct exact old/new satellite translation compared separately for E/C; nominal satellite clock numeric seed before translation is not exported by this diagnostic schema"}


def audit_run(root, expected_epochs):
    plan = json.loads((root / "stage_plan.json").read_text())
    if plan.get("state") == "RUNNING" or "exit_code" not in plan:
        return {"run_root": str(root), "state": plan.get("state"), "completed": False,
            "scientific_acceptance": False, "scope": "Still running; no final-output acceptance is inferred"}
    products = StageAuditor(root, 30, .1, 1000, 2000).run()
    qc, lifecycle = qc_transactions(root), clock_lifecycle(root)
    covariance = [solution["full_covariance"] for epoch in products["products"]["epochs"] for solution in epoch["solutions"].values()]
    grid = [gpst(epoch) for epoch in qc["committed_epochs"]]
    expected_grid = {START + 30 * number for number in range(expected_epochs)}
    binary = Path(plan["binary"])
    binary_actual = sha256(binary) if binary.is_file() else None
    checks = {"final_binary_sha_matches": plan["binary_sha256"] == FINAL_SHA and binary_actual == FINAL_SHA,
        "exit_zero_not_timeout": plan.get("exit_code") == 0 and plan.get("timed_out") is False,
        "exact_12_plus_1_actual_preflight": plan["effective_input_preflight"].get("passed") is True and
            plan["effective_input_preflight"].get("actual_unique_observation_files") == 12 and plan["effective_input_preflight"].get("actual_navigation_files") == 1,
        "processed_expected_epochs": products["logs"]["processed_epochs"].get("count") == expected_epochs,
        "qc_committed_exact_expected_grid": len(grid) == expected_epochs and set(grid) == expected_grid,
        "no_qc_rollback_or_authoritative_rejection": not qc["rollbacks"] and not qc["authoritative_rejection_records"],
        "no_failed_terminal_reconciliation": not any("converged=0" in row["text"] for row in qc["terminal_reconciliation_records"]),
        "atomic_hash_verified_expected_epochs": products["atomic"]["hash_verified_bundles"] == expected_epochs and not products["atomic"]["pending_directories_not_published"],
        "full_q_all_expected_epochs_consistent": len(covariance) == expected_epochs and all(q.get("internally_consistent") is True for q in covariance),
        "finite_all_phase_product_rows": all(total.get("phase_rows") == total.get("finite_numeric_rows") for total in products["products"]["totals_by_solution_system"].values()),
        "no_kalman_spp_clock_feedback_observed": lifecycle.get("no_kalman_clock_used_for_valid_spp_observations") is True,
        "no_invalid_lifecycle_records": not lifecycle.get("invalid_records"),
        "original_pins_unchanged_after_run": plan.get("binary_unchanged_after_run") is True and plan.get("pinned_inputs_unchanged_after_run") is True,
        "offset_formula_and_save_match": lifecycle.get("offset_checks", {}).get("formula_checked_adjustments", 0) > 0 and
            not lifecycle["offset_checks"]["failures"] and lifecycle["offset_checks"]["maximum_formula_residual_m"] <= 1e-8 and
            lifecycle["offset_checks"]["maximum_previous_saved_offset_difference_m"] <= 1e-8}
    return {"run_root": str(root), "completed": True, "scientific_acceptance": False,
        "operational_evidence_checks": checks, "operational_evidence_checks_all_pass": all(checks.values()),
        "binary_sha256": plan["binary_sha256"], "actual_binary_file_sha256": binary_actual,
        "stage_plan_sha256": sha256(root / "stage_plan.json"), "exact_inputs": plan["pinned_input_contract"],
        "product_evidence": products, "qc_transactions": qc, "clock_lifecycle": lifecycle,
        "resources": resource_statistics(root / "runner.log"),
        "performance_scope": "Trials overlapped other tests; timing/RSS are actual observations, not controlled performance benchmarks"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-base", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    parser.add_argument("--epochs", default=61, type=int)
    args = parser.parse_args()
    runs = {label: audit_run(args.run_base / ("final03_" + label + "_61ep"), args.epochs) for label in ("G", "E", "C", "GEC")}
    for label in ("E", "C"):
        if runs[label]["completed"]:
            # Old instrumented single-system initialization is a separate
            # baseline. This comparison changes no original output files.
            runs[label]["cold_seed_old_new_exact_translation"] = compare_cold_seeds(
                args.run_base / ("diagnostic_" + label + "_2ep"), args.run_base / ("final03_" + label + "_61ep"),
                reference_receiver="MAL2" if label == "C" else "MATE")
    report = {"schema_version": 1, "generated_utc": datetime.now(timezone.utc).isoformat(),
        "required_binary_sha256": FINAL_SHA, "scientific_acceptance": False,
        "scope": "12 pilot input stations, cold FLOAT 2024-07-17 00:00--00:30 GPST; no AR, no convergence/continuity/product scientific acceptance",
        "performance_scope": "Parallel tests are not a controlled performance benchmark", "runs": runs}
    args.output_root.mkdir(parents=True, exist_ok=True)
    (args.output_root / "final_real_replay_audit.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    lines = ["FINAL03 REAL COLD FLOAT PILOT REPLAY EVIDENCE", report["scope"], report["performance_scope"], ""]
    for label, run in runs.items():
        lines.append(label + ": " + ("COMPLETED" if run["completed"] else "RUNNING / NOT ACCEPTED"))
        if not run["completed"]:
            continue
        product = run["product_evidence"]
        lines.extend(["  operational evidence checks " + json.dumps(run["operational_evidence_checks"]),
            "  processed " + json.dumps(product["logs"]["processed_epochs"]),
            "  atomic hash verified " + str(product["atomic"]["hash_verified_bundles"]),
            "  phase row totals " + json.dumps(product["products"]["totals_by_solution_system"]),
            "  epoch timing " + json.dumps(product["logs"]["epoch_timing"]),
            "  resources " + json.dumps(run["resources"]),
            "  valid SPP clock sources " + json.dumps(run["clock_lifecycle"].get("valid_spp_observation_clock_sources")),
            "  offset " + json.dumps(run["clock_lifecycle"].get("offset_checks")),
            "  cold receiver g " + json.dumps(run["clock_lifecycle"].get("cold_receiver_translation_by_system")),
            "  evidence issues " + json.dumps(product["evidence_issues"]), ""])
    (args.output_root / "final_real_replay_audit.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
