#!/usr/bin/env python3
"""Isolated cold FLOAT clock/SPP controls; preparation never launches PEA.

This diagnostic harness deliberately does not grant build/scientific approval
or change run_stages.py defaults. configured versus precise-only changes only
the requested satellite clock source selection, never QC/noise thresholds.
"""
from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timedelta, timezone
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
from audit_stage_results import StageAuditor, render_text, unavailable

SYSTEMS = {"G": ["gps"], "E": ["gal"], "C": ["bds"],
           "GEC": ["gps", "gal", "bds"]}


def end_time(epochs):
    if not 1 <= epochs <= 2880:
        raise ValueError("epochs must be between 1 and 2880 within the same GPST day")
    return (datetime(2024, 7, 17) + timedelta(seconds=30 * (epochs - 1))).strftime("%H:%M:%S")


def trial_overlay(systems, epochs, clock_sources, root, audit, spp_adjust="configured"):
    selected = SYSTEMS[systems]
    # D supplies the already-reviewed exact G/E/C signal and membership rules.
    # Select only enabled systems in a fresh local object, never mutate STAGES.
    value = stages.overlay("D", root, audit)
    gnss = value["processing_options"]["gnss_general"]
    for system, options in gnss["sys_options"].items():
        if system not in selected:
            options.clear()
            options.update(process=False, ambiguity_resolution=False)
    for name in ("zhang_full_rank", "zhang_pppar"):
        gnss[name]["sys_options"] = {system: options for system, options in
            gnss[name]["sys_options"].items() if system in selected}
    gnss["zhang_full_rank"]["bds3_only"] = "bds" in selected
    gnss["pivot_receiver"] = "MAL2" if selected == ["bds"] else "MATE"
    if spp_adjust == "off":
        gnss["adjust_rec_clocks_by_spp"] = False
    elif spp_adjust != "configured":
        raise ValueError("Unknown SPP receiver-clock adjustment control")
    for satellite, options in value["satellite_options"].items():
        if re.fullmatch(r"[GEC]\d\d", satellite):
            if satellite[0] not in {stages.SYSTEM_LETTERS[s] for s in selected}:
                options["exclude"] = True
    value["receiver_options"]["global"]["rec_reference_system"] = {
        "gps": "GPS", "gal": "GAL", "bds": "BDS"}[selected[0]]
    value["processing_options"]["epoch_control"]["end_epoch"] = "2024-07-17 " + end_time(epochs)
    if clock_sources == "precise-only":
        value["satellite_options"]["global"] = {"models": {"clock": {"sources": ["PRECISE"]}}}
    elif clock_sources != "configured":
        raise ValueError("Unknown satellite clock source control")
    value["outputs"]["metadata"]["config_description"] = "clock_trial_" + systems + "_" + clock_sources
    return value


def trial_environment(diagnostics=()):
    env = stages.environment("D")
    assert env["ZHANG_R51_AR_START_GPST_SECONDS"] == str(stages.epoch_seconds("23:59:30"))
    for item in diagnostics:
        key, separator, value = item.partition("=")
        if not separator or not re.fullmatch(r"ZHANG_[A-Z0-9_]*CLOCK[A-Z0-9_]*DIAGNOSTIC[A-Z0-9_]*", key):
            raise ValueError("--diagnostics-env accepts an explicit ZHANG_*CLOCK*DIAGNOSTIC*=value only")
        env[key] = value
    if any("RESTORE" in key or "REUSE_FLOAT" in key for key in env):
        raise RuntimeError("Cold clock trials cannot import checkpoints")
    return env


def prepare(binary, audit_path, root, systems, epochs=2, clock_sources="configured", diagnostics=(), timeout=300,
            spp_adjust="configured"):
    if not math.isfinite(timeout) or timeout <= 0:
        raise ValueError("timeout must be finite and positive")
    end = end_time(epochs)
    audit = json.loads(audit_path.read_text())
    if not audit["selected_complete"] or audit["missing_files"]:
        raise RuntimeError("Input audit incomplete")
    binary = binary.resolve(strict=True)
    root = root.resolve()
    env = trial_environment(diagnostics)
    # Refuse reuse even of an empty/prepared directory; preserve every trial.
    root.mkdir(parents=True, exist_ok=False)
    for directory in ("outputs", "checkpoints", "config"):
        (root / directory).mkdir()
    path = root / "config/clock_trial_overlay.yaml"
    stages.write_json(path, trial_overlay(systems, epochs, clock_sources, root, audit, spp_adjust))
    configs = [str(stages.FROZEN_CONFIG / name) for name in stages.CONFIG_NAMES] + [str(path)]
    pinned = stages.pin_stage_inputs(audit, configs)
    plan = {"schema_version": 1, "stage": "CLOCK_" + systems + "_" + clock_sources + "_spp_" + spp_adjust,
        "specification": {"systems": SYSTEMS[systems], "mode": "FLOAT_CLOCK_DIAGNOSTIC", "end": end},
        "state": "PREPARED_NOT_RUN", "run_root": str(root), "binary": str(binary),
        "binary_sha256": stages.file_sha(binary), "input_audit_path": str(audit_path.resolve()),
        "input_audit_sha256": stages.file_sha(audit_path), "configs": configs,
        "expected_epoch_count": epochs, "cold_start": True, "reuse_gps_checkpoint": False,
        "start_epoch": "2024-07-17 00:00:00", "end_epoch": "2024-07-17 " + end,
        "ratio": 2, "all_nis": False, "epoch_only_ar": True,
        "ar_start_gpst_seconds": stages.epoch_seconds("23:59:30"),
        "clock_source_control": clock_sources, "spp_receiver_adjustment_control": spp_adjust, "timeout_seconds": timeout,
        "control_interpretation": {
            "precise_only": "Global clock.sources also changes PPP clock components; not an SPP-only causal intervention",
            "spp_adjust_off": "Disables only SPP receiver-clock adjustment; diagnostic ablation, not an acceptable receiver-jump handling fix"},
        "diagnostic_control_not_a_fix": True, "scientific_acceptance": False,
        "qc_statistical_noise_models": "Original six YAMLs unchanged; no thresholds/noise/height cutoff overridden",
        "geometry_risk": "GAL AGGO local clock/ZTD rank mode is separate from common clock QC failure; station network unchanged",
        "pinned_input_contract": pinned,
        "driver_source_identity": [{"path": str(path.resolve()), "sha256": stages.file_sha(path)}
            for path in (Path(__file__), Path(stages.__file__), Path(__file__).with_name("audit_stage_results.py"))],
        "effective_input_preflight": {"status": "NOT_RUN", "required_before_solver_launch": True},
        "environment": {k: v for k, v in env.items() if k.startswith("ZHANG_") or k.endswith("_NUM_THREADS")}}
    plan["command"] = ["/usr/bin/time", "-v", *stages.independent_stage_command(binary, configs, pinned)]
    stages.write_json(root / "stage_plan.json", plan)
    return plan, env


def preflight(plan, env):
    binary = Path(plan["binary"])
    if stages.file_sha(binary) != plan["binary_sha256"]:
        raise RuntimeError("Pinned binary changed before dry-run")
    result = stages.configuration_preflight(binary, plan["configs"], plan["pinned_input_contract"],
                                            Path(plan["run_root"]), env)
    plan["effective_input_preflight"] = result
    if not result["passed"]:
        plan["state"] = "CONFIGURATION_PREFLIGHT_FAILED_NOT_RUN"
    stages.write_json(Path(plan["run_root"]) / "stage_plan.json", plan)
    if not result["passed"]:
        raise RuntimeError("Solver not launched: " + result["reason"])
    return result


def execute(plan, env):
    if sys.platform != "linux":
        raise RuntimeError("Execution requires the native Ubuntu-22.04-G environment")
    root = Path(plan["run_root"])
    preflight(plan, env)
    stages.verify_pinned_stage_inputs(plan["pinned_input_contract"])
    if stages.file_sha(plan["binary"]) != plan["binary_sha256"]:
        raise RuntimeError("Pinned binary changed before solver launch")
    plan.update(state="RUNNING", started_utc=datetime.now(timezone.utc).isoformat())
    stages.write_json(root / "stage_plan.json", plan)
    started = time.monotonic()
    timed_out = False
    with (root / "runner.log").open("wb") as log:
        process = subprocess.Popen(plan["command"], cwd=stages.NATIVE_INPUT, env=env,
            stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
        plan["pid"] = process.pid
        stages.write_json(root / "stage_plan.json", plan)
        try:
            process.wait(timeout=plan["timeout_seconds"])
        except subprocess.TimeoutExpired:
            timed_out = True
            os.killpg(process.pid, signal.SIGTERM)
            try:
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGKILL)
                process.wait()
    plan.update(state="COMPLETED_REQUIRES_SCIENTIFIC_AUDIT" if process.returncode == 0 and not timed_out else "EXECUTION_FAILED",
        exit_code=process.returncode, timed_out=timed_out, elapsed_seconds=time.monotonic() - started,
        finished_utc=datetime.now(timezone.utc).isoformat(), scientific_acceptance=False)
    plan["binary_unchanged_after_run"] = stages.file_sha(plan["binary"]) == plan["binary_sha256"]
    try:
        stages.verify_pinned_stage_inputs(plan["pinned_input_contract"])
        plan["pinned_inputs_unchanged_after_run"] = True
    except RuntimeError as error:
        plan["pinned_inputs_unchanged_after_run"] = False
        plan["input_integrity_error"] = str(error)
    stages.write_json(root / "stage_plan.json", plan)
    audit_trial(root)
    return process.returncode == 0 and not timed_out


def audit_trial(root):
    """Read-only product/atomic/full-Q audit plus explicit boundary availability.

    Old binaries do not print the clock lifecycle. Absence is UNAVAILABLE,
    not evidence that a reference/SPP transition did or did not occur.
    """
    report = {"schema_version": 1, "generated_utc": datetime.now(timezone.utc).isoformat(),
        "scientific_acceptance": False,
        "stages": [StageAuditor(root, 30, .1, 1000, 2000).run()]}
    report["clock_lifecycle"] = read_clock_lifecycle(root)
    report["resource_statistics"] = resource_statistics(root / "runner.log")
    stages.write_json(root / "clock_trial_audit.json", report)
    (root / "clock_trial_audit.txt").write_text(render_text(report) +
        "\nClock lifecycle capture: " + report["clock_lifecycle"]["status"] +
        "; causal proof and scientific acceptance remain false.\n", encoding="utf-8")
    return report


def resource_statistics(path):
    """Only /usr/bin/time's explicit records; a timeout can leave them absent."""
    fields = {"user_cpu_seconds": r"User time \(seconds\):\s*([\d.]+)",
        "system_cpu_seconds": r"System time \(seconds\):\s*([\d.]+)",
        "maximum_rss_kib": r"Maximum resident set size \(kbytes\):\s*(\d+)",
        "major_page_faults": r"Major \(requiring I/O\) page faults:\s*(\d+)",
        "minor_page_faults": r"Minor \(reclaiming a frame\) page faults:\s*(\d+)",
        "swaps": r"Swaps:\s*(\d+)"}
    found = {}
    if path.is_file():
        with path.open(errors="replace") as stream:
            for line in stream:
                for field, expression in fields.items():
                    match = re.search(expression, line)
                    if match:
                        found[field] = float(match[1]) if "seconds" in field else int(match[1])
    return {"status": "AVAILABLE" if found else "UNAVAILABLE", "values": found,
        "missing_fields": sorted(set(fields) - set(found)),
        "scope": "Explicit GNU time child-process statistics, not inferred memory/no-OOM acceptance"}


def read_clock_lifecycle(root):
    """Stream INFO JSON evidence; retain exact normalized events separately.

    Clock/source values are measurements of internal lifecycle, not proof of
    a physical receiver jump, correct gauge transport or full-P equivalence.
    Duplicate TRACE copies are deliberately not scanned with runner stdout.
    """
    names = ("satellite_clock_source_and_chart", "spp_clock", "saved_spp_ppp_offset_system_and_chart",
             "kf_receiver_clock_and_rate", "applied_adjustment_and_innovation", "reference_transport_full_covariance")
    counts, phases, systems, sources = Counter(), Counter(), Counter(), Counter()
    spp_results = {}
    invalid, records, samples = [], 0, []
    path = root / "runner.log"
    normalized = root / "clock_lifecycle_events.jsonl"
    def reject_constant(value):
        raise ValueError("Nonfinite JSON constant " + value)
    with normalized.open("w", encoding="utf-8") as output:
        if path.is_file():
            with path.open(errors="replace") as stream:
                for number, line in enumerate(stream, 1):
                    marker = "ZHANG_CLOCK_LIFECYCLE "
                    if marker not in line:
                        continue
                    try:
                        event = json.loads(line.split(marker, 1)[1], parse_constant=reject_constant)
                        if not isinstance(event, dict) or event.get("schema") != "ZHANG_CLOCK_LIFECYCLE_V1" or not event.get("phase"):
                            raise ValueError("Missing/unknown lifecycle schema or phase")
                        if any(not isinstance(event.get(key, []), list) for key in ("model_clocks", "observations")):
                            raise ValueError("Lifecycle model/observation arrays are malformed")
                        json.dumps(event, allow_nan=False)
                    except (ValueError, TypeError) as error:
                        if len(invalid) < 50:
                            invalid.append({"source_line": number, "reason": str(error)})
                        counts["invalid_records"] += 1
                        continue
                    records += 1
                    phases[event["phase"]] += 1
                    if event["phase"] == "SPP_RESULT":
                        identity = tuple(event.get(key) for key in ("epoch", "receiver", "spp_system"))
                        status = event.get("spp_status")
                        if all(identity) and isinstance(status, str) and status:
                            if identity in spp_results and spp_results[identity]["status"] != status:
                                counts["conflicting_spp_result_identities"] += 1
                            spp_results[identity] = {"epoch": identity[0], "receiver": identity[1],
                                "selected_spp_system": identity[2], "status": status, "source_line": number,
                                "valid_satellite_sources": dict(Counter(str(obs.get("source", "UNAVAILABLE"))
                                    for obs in event.get("observations", []) if isinstance(obs, dict) and obs.get("spp_valid") is True))}
                    wrapped = {"source": str(path), "source_line": number, "event": event}
                    output.write(json.dumps(wrapped, allow_nan=False) + "\n")
                    if len(samples) < 12:
                        samples.append(wrapped)
                    if event.get("spp_clock_m") is not None:
                        counts["spp_clock"] += 1
                    if event.get("saved_offset_m") is not None and event.get("saved_offset_system") is not None:
                        counts["saved_spp_ppp_offset_system_and_chart"] += 1
                    if event.get("applied_adjustment_m") is not None or event.get("unrounded_adjustment_m") is not None:
                        counts["applied_adjustment_and_innovation"] += 1
                    for model in event.get("model_clocks", []):
                        if isinstance(model, dict):
                            systems[str(model.get("system", "UNAVAILABLE"))] += 1
                            if model.get("kf_clock_m") is not None:
                                counts["kf_receiver_clock_and_rate"] += 1
                    for obs in event.get("observations", []):
                        if isinstance(obs, dict) and obs.get("source") is not None:
                            sources[str(obs["source"])] += 1
                            counts["satellite_clock_source_and_chart"] += 1
    # Existing schema records source selection, not a formally versioned chart
    # certificate or complete covariance. Availability never asserts that proof.
    boundaries = {name: {"status": "AVAILABLE", "nonnull_records": counts[name],
        "scope": "Explicit captured values only; chart identity/full covariance equivalence not certified"}
        if counts[name] else unavailable("No explicit non-null lifecycle value for this boundary") for name in names}
    return {"status": "AVAILABLE" if records else "UNAVAILABLE", "record_count": records,
        "phase_counts": dict(phases), "model_system_record_counts": dict(systems),
        "satellite_source_counts": dict(sources), "boundaries": boundaries,
        "invalid_record_count": counts["invalid_records"], "invalid_record_samples": invalid,
        "event_file": str(normalized), "event_file_sha256": stages.file_sha(normalized),
        "source_log_sha256": stages.file_sha(path) if path.is_file() else None,
        "spp_result_status": {"status": "AVAILABLE", "unique_epoch_receiver_selected_system_records": len(spp_results),
            "status_counts": dict(Counter(item["status"] for item in spp_results.values())),
            "conflicting_identity_count": counts["conflicting_spp_result_identities"],
            "records": sorted(spp_results.values(), key=lambda item: (item["epoch"], item["receiver"], item["selected_spp_system"])),
            "scope": "Explicit SPP_RESULT status only for its selected_spp_system; neither SINGLE_X nor FAILED is SINGLE; not a per-system marginal positioning proof"}
            if spp_results else unavailable("No explicit SPP_RESULT epoch/receiver/selected-system/status captures"),
        "samples": samples, "causal_proof": False, "clock_control_is_not_scientific_acceptance": True,
        "formal_chart_identity": unavailable("Lifecycle V1 exposes source selection and numeric values, not a certified complete chart-transport identity"),
        "missing_capture_warning": "Legacy binary has no lifecycle captures; no event is not a zero adjustment"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--systems", choices=SYSTEMS, required=True)
    parser.add_argument("--clock-sources", choices=("configured", "precise-only"), default="configured")
    parser.add_argument("--spp-adjust", choices=("configured", "off"), default="configured")
    parser.add_argument("--epochs", type=int, default=2)
    parser.add_argument("--timeout-seconds", type=float, default=300)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--audit", type=Path, required=True)
    parser.add_argument("--run-root", type=Path, required=True)
    parser.add_argument("--diagnostics-env", action="append", default=[], metavar="NAME=VALUE")
    parser.add_argument("--preflight", action="store_true", help="Run real parser dry-run only; no numerical solver")
    parser.add_argument("--execute", action="store_true", help="Launch this explicitly authorized FLOAT diagnostic, not formal AR")
    args = parser.parse_args()
    plan, env = prepare(args.binary, args.audit, args.run_root, args.systems,
        args.epochs, args.clock_sources, args.diagnostics_env, args.timeout_seconds, args.spp_adjust)
    if args.execute:
        raise SystemExit(0 if execute(plan, env) else 1)
    if args.preflight:
        preflight(plan, env)
    print(json.dumps(plan, indent=2))


if __name__ == "__main__":
    main()
