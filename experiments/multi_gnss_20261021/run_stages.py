#!/usr/bin/env python3
"""Prepare or explicitly execute isolated A--G multi-system experiments.

Preparation is not a run or an acceptance result. Execution requires a fresh
binary-pinned build/test approval JSON. B--G never import the GPS checkpoint.
Use Linux/WSL for --execute; original inputs/outputs are never overwritten.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import re
import signal
import subprocess
import sys
import time

NATIVE_INPUT = Path("/mnt/d/GINAN_R20/inputData")
DOWNLOADED_NAV3 = NATIVE_INPUT / "multi_gnss_2024199_20261021/navigation/BRDC00IGS_R_20241990000_01D_MN.rnx"
FROZEN_CONFIG = Path("/home/rx/GINAN/frozen-r51-satellite-fix-20260918/config")
CHECKPOINT = Path("/mnt/d/GINAN_R20/zhang_r51_ratio3_float2h_ar1h_2024199_180_20260917_checkpoints_source_before_ratioonly/E29-20240717015930-e240-rb06f1454a2e1")
BASELINE_BINARY = Path("/home/rx/GINAN/bridge_continuity_20261020/verified_v2/pea")
CONFIG_NAMES = ["zhang_global_2024199_180_base.yaml", "zhang_global_2024199_180_inputs.yaml",
                "zhang_global_2024199_180_product.yaml", "zhang_global_2024199_180_e29_180network_product_1h.yaml",
                "zhang_global_2024199_180_e29_service_30s_1h.yaml", "case.yaml"]
SIGNALS = {"gps": ["L1C", "L2W"], "gal": ["L1C", "L5Q"], "bds": ["L1P", "L5P"]}
NAV_TYPES = {"gps": "LNAV", "gal": "FNAV", "bds": "D1"}
SYSTEM_LETTERS = {"gps": "G", "gal": "E", "bds": "C"}
STAGES = {
    "A": {"systems": ["gps"], "mode": "PINNED_GPS_REGRESSION", "end": "02:00:00"},
    "B_E": {"systems": ["gal"], "mode": "FLOAT", "end": "02:05:00"},
    "B_C": {"systems": ["bds"], "mode": "FLOAT", "end": "02:05:00"},
    "C_GE": {"systems": ["gps", "gal"], "mode": "FLOAT", "end": "02:05:00"},
    "C_GC": {"systems": ["gps", "bds"], "mode": "FLOAT", "end": "02:05:00"},
    "D": {"systems": ["gps", "gal", "bds"], "mode": "FLOAT", "end": "02:05:00"},
    "E": {"systems": ["gps", "gal", "bds"], "mode": "AR", "end": "02:05:00"},
    "F": {"systems": ["gps", "gal", "bds"], "mode": "INDEPENDENT_USER_VALIDATION", "end": "02:05:00"},
    "G": {"systems": ["gps", "gal", "bds"], "mode": "AR", "end": "02:59:30"},
}


def file_sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(2**20), b""):
            h.update(chunk)
    return h.hexdigest()


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")


def pin_stage_inputs(audit, configs):
    """Bind an exact independent B--G subset; YAML sequence layers append.

    Do not flatten/reinterpret physical models. Production's --rnx_inputs
    clears the accumulated map; --nav_files overrides YAML after each layer.
    The real --dry-run parser must confirm both sets before PEA runs.
    """
    selected = audit["selected_12_station_candidates"]
    if len(selected) != 12 or len(set(selected)) != 12:
        raise RuntimeError("Independent stages require exactly 12 unique selected stations")
    headers = {row["filename"][:4]: row["filename"] for row in audit["station_headers"]}
    observations = []
    for station in selected:
        filename = headers[station]
        if Path(filename).name != filename or any(c in filename for c in "*?[]"):
            raise RuntimeError("Observation selection must name exact basenames, not globs or extra paths")
        observations.append((station, (NATIVE_INPUT / "data" / filename).resolve()))
    def identity(path, **extra):
        path = Path(path).resolve()
        if any(c.isspace() for c in str(path)):
            raise RuntimeError("Effective-summary path validation does not support whitespace paths")
        return dict(path=str(path), bytes=path.stat().st_size, sha256=file_sha(path), **extra)
    return {"selection_policy": "EXPLICIT_COMMAND_LINE_REPLACEMENT_REAL_PARSER_VERIFIED",
        "expected_unique_station_count": 12,
        "observations": [identity(path, station=station) for station, path in observations],
        "navigation": [identity(DOWNLOADED_NAV3)],
        "configurations": [identity(path) for path in configs],
        "cnv123_propagation_claim": False, "scientific_acceptance": False}


def exact_input_cli(pinned):
    observations = [row["path"] for row in pinned["observations"]]
    navigation = [row["path"] for row in pinned["navigation"]]
    if len(observations) != 12 or len(set(observations)) != 12 or len(navigation) != 1:
        raise RuntimeError("Pinned command-line input set is not exactly 12 unique RINEX plus one NAV")
    return ["--rnx_inputs", *observations, "--nav_files", *navigation]


def independent_stage_command(binary, configs, pinned, dry_run=False):
    # B--G need info-level processed-epoch/timing/end markers for auditing.
    # The isolated legacy A command is intentionally constructed separately.
    command = [str(binary), "-y", *configs, "-d", "case", *exact_input_cli(pinned)]
    return command + (["--dry-run"] if dry_run else [])


def verify_pinned_stage_inputs(pinned):
    exact_input_cli(pinned)
    for group in ("observations", "navigation", "configurations"):
        for row in pinned[group]:
            path = Path(row["path"])
            if not path.is_file() or path.stat().st_size != row["bytes"] or file_sha(path) != row["sha256"]:
                raise RuntimeError("Pinned input changed before solver launch: " + str(path))


def parse_effective_input_summary(text):
    # acsConfig.info prints these final collections, unlike FILE/RAW_CONFIG.
    text = re.sub(r"\x1b\[[0-?]*[ -/]*[@-~]", "", text)
    found = {}
    keys = {"rnx_inputs", "nav_files", "ubx_inputs", "sbf_inputs", "rtcm_inputs",
            "sp3_inputs", "snx_inputs", "custom_inputs", "qzl6_inputs", "sisnet_inputs"}
    for line in text.splitlines():
        match = re.fullmatch(r"\s*([a-z0-9_]+):\s*(.*?)\s*", line)
        if match and match[1] in keys:
            if match[1] in found:
                raise RuntimeError("Ambiguous repeated effective input summary: " + match[1])
            found[match[1]] = match[2].split()
    if not found.get("rnx_inputs") or not found.get("nav_files"):
        raise RuntimeError("Real parser did not emit effective rnx_inputs and nav_files summary")
    return found


def verify_effective_input_summary(text, pinned):
    found = parse_effective_input_summary(text)
    for key, group, count in (("rnx_inputs", "observations", 12), ("nav_files", "navigation", 1)):
        expected = {row["path"] for row in pinned[group]}
        actual = found[key]
        if len(actual) != count or len(set(actual)) != count or set(actual) != expected:
            raise RuntimeError("Real parser effective input mismatch for " + key +
                               f": expected {count} exact paths, got {len(actual)} / {len(set(actual))} unique")
    extra = {key: value for key, value in found.items() if key not in ("rnx_inputs", "nav_files") and value}
    if extra:
        raise RuntimeError("Unexpected additional observation/navigation streams: " + str(sorted(extra)))
    return {"passed": True, "actual_unique_observation_files": 12, "actual_navigation_files": 1,
        "actual_rnx_inputs": found["rnx_inputs"], "actual_nav_files": found["nav_files"],
        "scope": "Actual PEA configuration parsing only; no observations or numerical solver initialized"}


def configuration_preflight(binary, configs, pinned, root, env):
    command = independent_stage_command(binary, configs, pinned, dry_run=True)
    log = root / "config_preflight.log"
    try:
        verify_pinned_stage_inputs(pinned)
        completed = subprocess.run(command, cwd=NATIVE_INPUT, env=env, stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT, text=True, timeout=120)
        log.write_text(completed.stdout, encoding="utf-8")
        result = {"command": command, "exit_code": completed.returncode, "dump_sha256": file_sha(log), "log": str(log)}
        if completed.returncode != 0:
            raise RuntimeError("PEA --dry-run did not exit successfully")
        if "Dry-run successful: configuration and sanity checks passed" not in completed.stdout:
            raise RuntimeError("PEA emitted no explicit successful dry-run/configuration/sanity marker")
        result.update(verify_effective_input_summary(completed.stdout, pinned))
        # Detect mutations during the parse as well as before it.
        verify_pinned_stage_inputs(pinned)
        return result
    except (OSError, RuntimeError, subprocess.TimeoutExpired) as error:
        return {"passed": False, "command": command, "reason": str(error), "log": str(log),
                "dump_sha256": file_sha(log) if log.is_file() else None}


def epoch_seconds(hour):
    # GTime's solver schedule uses elapsed GPST seconds since 1980-01-06,
    # NOT Unix seconds. Input calendar is already GPST: do not apply UTC leaps.
    day = datetime.fromisoformat("2024-07-17T" + hour).replace(tzinfo=timezone.utc)
    origin = datetime(1980, 1, 6, tzinfo=timezone.utc)
    return int((day - origin).total_seconds())


def environment(stage):
    env = {k: v for k, v in os.environ.items() if not k.startswith("ZHANG_")}
    env.update(ZHANG_R51_ENABLE="1", ZHANG_R51_RATIO_ONLY="1",
               ZHANG_R51_RATIO_THRESHOLD="2", ZHANG_R51_EPOCH_ONLY_AR="1",
               ZHANG_R51_SINGLE_BASIS_HANDOFF="1", ZHANG_R51_PRODUCT_REPAIR="1",
               ZHANG_R51_INTERN_HANDOFF_PROOFS="1", ZHANG_R49_FUSION_WEIGHT="0.10",
               ZHANG_P0_RESOURCE_PROBE="1", OMP_NUM_THREADS="8", OPENBLAS_NUM_THREADS="1",
               MKL_NUM_THREADS="1", ZHANG_R51_FLOAT_THREADS="12")
    env["ZHANG_R51_AR_START_GPST_SECONDS"] = str(epoch_seconds("02:00:00" if STAGES[stage]["mode"] == "AR" or stage == "A" else "23:59:30"))
    if stage == "A":
        env.update(ZHANG_R51_REUSE_FLOAT_20240717="1", ZHANG_R51_REUSE_FLOAT_THREADS_8_1="1",
                   ZHANG_E29_RESTORE_BUNDLE=str(CHECKPOINT))
    else:
        assert not any("RESTORE" in k or "REUSE_FLOAT" in k for k in env)
    return env


def overlay(stage, root, audit):
    spec = STAGES[stage]
    systems = spec["systems"]
    by_station = {row["filename"][:4]: row for row in audit["station_headers"]}
    stations = audit["selected_12_station_candidates"]
    wanted_satellites = set(audit["product_membership"]["sp3_satellites"]) & set(audit["product_membership"]["clk_satellites"])
    c_white = [row["prn"] for row in audit["whitelist"]]
    full_rank = {"enable": True, "output_diagnostics": True,
                 "bds3_only": "bds" in systems, "bds3_satellite_whitelist": c_white,
                 "sys_options": {}}
    product_options = {"sys_options": {}, "deterministic_checkpoint": False,
        "checkpoint_restore_path": "", "product_solution": "FIXED", "product_relation_feedback": False,
        "product_filename": str(root / "outputs/zhang_internal_products.csv"),
        "product_covariance_filename": str(root / "outputs/zhang_internal_product_covariance.csv"),
        "output_products": True, "output_diagnostics": True}
    gnss = {"pivot_receiver": "MAL2" if systems == ["bds"] else "MATE",
            "use_tgd_bias": False, "minimise_ionosphere_offsets": False,
            "sys_options": {}, "zhang_full_rank": full_rank, "zhang_pppar": product_options}
    for system in ["gps", "gal", "bds", "glo", "qzs", "sbs", "irn", "leo"]:
        enabled = system in systems
        value = {"process": enabled, "ambiguity_resolution": enabled}
        if enabled:
            value.update(code_priorities=SIGNALS[system], reject_eclipse=False,
                         receiver_amb_pivot=False, network_amb_pivot=False,
                         used_nav_type=NAV_TYPES[system])
            ref = audit["reference_stations"][SYSTEM_LETTERS[system]]
            full_rank["sys_options"][system] = {
                "baseline_observables": SIGNALS[system], "reference_receiver": ref,
                "use_spanning_tree": True, "auto_reference_switch": False,
                "state_edge_grace_epochs": 10, "prefer_historical_edges": True,
                "core_skeleton": False, "product_core_min_satellite_support": 2}
            product_options["sys_options"][system] = {"baseline_observables": SIGNALS[system]}
        gnss["sys_options"][system] = value
    # Explicit day-specific satellite selection; BDS-3 provenance is the SINEX
    # whitelist, never PRN >= some threshold. Undo stale Galileo exclusions.
    satellite_options = {}
    for system, letter in SYSTEM_LETTERS.items():
        for prn in range(1, 64):
            satellite = letter + f"{prn:02d}"
            keep = system in systems and satellite in wanted_satellites
            if system == "gps" and satellite in ("G01", "G04", "G30"):
                keep = False  # Same pre-existing G-only analysis baseline.
            if system == "bds":
                keep &= satellite in c_white
            satellite_options[satellite] = {"exclude": not keep}
    # Explicit satellite clock-code definitions are tracked in the plan. Their
    # physical WUM code/phase datum compatibility is an execution review gate.
    satellite_options["GAL"] = {"clock_codes": SIGNALS["gal"]}
    satellite_options["BDS"] = {"clock_codes": SIGNALS["bds"]}
    return {"inputs": {"satellite_data": {"nav_files": [str(DOWNLOADED_NAV3)]},
            "gnss_observations": {"gnss_observations_root": str(NATIVE_INPUT / "data") + "/",
            "rnx_inputs": [by_station[s]["filename"] for s in stations]}},
        "processing_options": {"epoch_control": {"start_epoch": "2024-07-17 00:00:00",
            "end_epoch": "2024-07-17 " + spec["end"], "epoch_interval": 30, "wait_next_epoch": 86400},
            "gnss_general": gnss, "ambiguity_resolution": {"mode": "LAMBDA",
                "solution_ratio_threshold": 2, "success_rate_threshold": 0.999,
                "once_per_epoch": True, "fix_and_hold": False}},
        "satellite_options": satellite_options,
        "receiver_options": {"global": {"rec_reference_system": {"gps": "GPS", "gal": "GAL", "bds": "BDS"}[systems[0]]}},
        "outputs": {"outputs_root": str(root / "outputs"),
            "metadata": {"config_description": "multi_gnss_20261021_" + stage},
            "trace": {"level": 2, "output_residuals": True},
            "sp3": {"output": False}, "clocks": {"output": False}, "bias_sinex": {"output": False}}}


def approval(path, binary, stage):
    if path is None:
        raise RuntimeError("--execute requires --build-approval JSON; preparation is not approval")
    document = json.loads(path.read_text())
    required = ["build_passed", "mandatory_tests_passed", "approved_for_execution"]
    if stage != "A":
        required += ["multi_system_clock_model_review_passed", "precise_clock_datum_review_passed"]
    missing = [key for key in required if document.get(key) is not True]
    if missing:
        raise RuntimeError("Execution gates false or missing: " + ", ".join(missing))
    if document.get("binary_sha256") != file_sha(binary):
        raise RuntimeError("Build approval binary hash mismatch")
    if stage == "G" and document.get("short_ar_acceptance_passed") is not True:
        raise RuntimeError("G persistent trial requires actual short AR acceptance")
    return document


def compare_bundle(left, right):
    result = {}
    for filename in ("products.csv", "covariance.csv", "certificate.txt"):
        lpath, rpath = left / filename, right / filename
        item = {"left_sha256": file_sha(lpath), "right_sha256": file_sha(rpath)}
        item["byte_equal"] = item["left_sha256"] == item["right_sha256"]
        if filename.endswith(".csv"):
            with lpath.open(newline="") as stream:
                lrows = list(csv.DictReader(stream))
            with rpath.open(newline="") as stream:
                rrows = list(csv.DictReader(stream))
            item.update(left_rows=len(lrows), right_rows=len(rrows), exact_rows_equal=lrows == rrows)
            float_fields = ({"clock_m", "clock_sigma_m", "phase_m", "phase_sigma_m",
                "clock_phase_covariance_m2", "correction_m", "correction_sigma_m",
                "float_correction_shift_m"} if filename == "products.csv" else {"covariance_m2"})
            item["numeric_maximum_absolute_difference_by_field"] = {}
            item["identity_flag_or_nonfinite_mismatches"] = []
            for index, (lrow, rrow) in enumerate(zip(lrows, rrows)):
                for field in lrow.keys() | rrow.keys():
                    lvalue, rvalue = lrow.get(field), rrow.get(field)
                    if lvalue == rvalue:
                        continue
                    if field in float_fields:
                        try:
                            lnumber, rnumber = float(lvalue), float(rvalue)
                            if math.isfinite(lnumber) and math.isfinite(rnumber):
                                d = abs(lnumber - rnumber)
                                item["numeric_maximum_absolute_difference_by_field"][field] = max(
                                    d, item["numeric_maximum_absolute_difference_by_field"].get(field, 0.0))
                                continue
                        except (ValueError, TypeError):
                            pass
                    item["identity_flag_or_nonfinite_mismatches"].append({
                        "row": index + 2, "field": field, "baseline": lvalue, "candidate": rvalue})
        result[filename] = item
    result["all_exact_equal"] = all(value["byte_equal"] for value in result.values() if isinstance(value, dict))
    result["scope"] = "Full first committed products, covariance and integer certificate; no curve-only comparison"
    return result


def run_until_commit(command, root, env, timeout):
    log = root / "runner.log"
    with log.open("wb") as stream:
        process = subprocess.Popen(command, cwd=NATIVE_INPUT, env=env, stdout=stream,
                                   stderr=subprocess.STDOUT, start_new_session=True)
        start = time.monotonic()
        committed = None
        while process.poll() is None:
            candidates = list((root / "outputs/zhang_internal_products.csv.epochs").glob("*/COMMIT.json"))
            if candidates:
                committed = sorted(candidates)[0].parent
                os.killpg(process.pid, signal.SIGINT)
                break
            if time.monotonic() - start > timeout:
                os.killpg(process.pid, signal.SIGTERM)
                break
            time.sleep(0.25)
        try:
            process.wait(timeout=60)
        except subprocess.TimeoutExpired:
            os.killpg(process.pid, signal.SIGKILL)
            process.wait()
    return {"pid": process.pid, "exit_code": process.returncode,
            "controlled_stop_after_first_commit": committed is not None,
            "bundle": str(committed) if committed else None, "elapsed_seconds": time.monotonic() - start}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage", choices=STAGES, required=True)
    parser.add_argument("--audit", type=Path, required=True)
    parser.add_argument("--run-root", type=Path, required=True)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--baseline-binary", type=Path, default=BASELINE_BINARY)
    parser.add_argument("--build-approval", type=Path)
    parser.add_argument("--execute", action="store_true")
    parser.add_argument("--timeout-seconds", type=float, default=14400)
    parser.add_argument("--inside-a-namespace", action="store_true", help=argparse.SUPPRESS)
    args = parser.parse_args()
    audit = json.loads(args.audit.read_text())
    if not audit["selected_complete"] or audit["missing_files"]:
        raise RuntimeError("Input audit incomplete")
    if args.execute:
        if sys.platform != "linux":
            raise RuntimeError("Execute only inside Ubuntu-22.04-G")
        approval(args.build_approval, args.binary, args.stage)
    root = args.run_root.resolve()
    # Never overwrite an experiment, even a prepared-but-not-executed one.
    if not args.inside_a_namespace:
        root.mkdir(parents=True, exist_ok=False)
    for directory in ("outputs", "checkpoints", "config"):
        (root / directory).mkdir(exist_ok=True)
    plan = dict(stage=args.stage, specification=STAGES[args.stage], state="PREPARED_NOT_RUN",
                binary=str(args.binary), ratio=2, all_nis=False, epoch_only_ar=True,
                reuse_gps_checkpoint=args.stage == "A", cold_start=args.stage != "A",
                scientific_acceptance=False, run_root=str(root), input_audit_path=str(args.audit.resolve()),
                input_audit_sha256=file_sha(args.audit))
    plan["environment"] = {key: value for key, value in environment(args.stage).items()
                           if key.startswith("ZHANG_") or key.endswith("_NUM_THREADS")}
    if args.stage == "A":
        plan["configs"] = [str(FROZEN_CONFIG / name) for name in CONFIG_NAMES]
        plan["stop_policy"] = "Keep all six pinned YAMLs byte-identical; stop the launched child process group after first atomic COMMIT"
        plan["baseline_binary"] = str(args.baseline_binary)
    elif args.stage == "F":
        plan["state"] = "BLOCKED_UNTIL_VALIDATED_MULTI_SYSTEM_PRODUCTS_AND_INDEPENDENT_USER_CONFIGURATION"
        plan["required_evidence"] = ["Static independent user coordinates and exact target signal pairs",
            "FLOAT and AR solutions with matching orbit/clock/code datum and nonzero cross-system product covariance",
            "Elementwise full covariance reconstruction and application-sign check",
            "Convergence and fixed-output counts from actual solutions, not a generated configuration"]
    else:
        value = overlay(args.stage, root, audit)
        path = root / "config/overlay.yaml"
        write_json(path, value)  # JSON is a strict subset of YAML, accepted by yaml-cpp.
        plan["configs"] = [str(FROZEN_CONFIG / name) for name in CONFIG_NAMES] + [str(path)]
        plan["navigation_input"] = {"path": str(DOWNLOADED_NAV3),
            "sha256": file_sha(DOWNLOADED_NAV3),
            "cnv123_propagation_claim": False,
            "role": "Traditional mixed broadcast navigation for SPP initialization; final service uses precise orbit and estimated IF/GF clocks"}
        plan["pinned_input_contract"] = pin_stage_inputs(audit, plan["configs"])
        plan["effective_input_preflight"] = {"status": "NOT_RUN",
            "required_before_solver_launch": True,
            "reason": "YAML rnx_inputs sequences append; exact CLI replacement must be checked by the real parser"}
        plan["warmup"] = "Cold FLOAT from 00:00; AR stages begin only at 02:00; no GPS checkpoint migration"
        plan["expected_epoch_count"] = int((epoch_seconds(STAGES[args.stage]["end"]) - epoch_seconds("00:00:00")) / 30) + 1
        if args.stage == "G":
            plan["faults_must_be_verified_separately"] = ["Controlled E11 physical arc break must not mutate G/C integer authority",
                "Reference-tree switch with fixed physical observations preserves complete covariance",
                "Per-system missing-product recovery and bounded old-arc dependency lifecycle"]
            plan["required_runtime_metrics"] = ["per-system FLOAT/WL/L1/product/bridge components", "epoch total P50/P95/P99/max", "peak RSS and history-length trend"]
    write_json(root / "stage_plan.json", plan)
    if not args.execute:
        print(json.dumps(plan, indent=2))
        return
    if args.stage == "F":
        raise RuntimeError("F is not an implemented user experiment; supply and validate a real independent user setup before declaring completion")
    if args.stage == "A":
        if os.geteuid() != 0:
            raise RuntimeError("A requires root ONLY to create a private mount namespace; do not change original config output paths")
        if not args.inside_a_namespace:
            command = ["unshare", "--mount", "--propagation", "private", sys.executable, str(Path(__file__).resolve()),
                       *sys.argv[1:], "--inside-a-namespace"]
            returncode = subprocess.call(command)
            raise SystemExit(returncode)
        subprocess.run(["mount", "--make-rprivate", "/"], check=True)
        target_output = Path("/mnt/d/GINAN_R20/inputData/outputs/zhang_r51_satfix_ratio3_ar1h_2024199_180_20260918")
        target_checkpoint = Path("/mnt/d/GINAN_R20/zhang_r51_satfix_ratio3_ar1h_2024199_180_20260918_checkpoints")
        if not target_output.is_dir() or not target_checkpoint.is_dir():
            raise RuntimeError("Original explicit mount destinations missing; refusing to create/modify them")
        results = []
        for label, binary in (("baseline", args.baseline_binary), ("candidate", args.binary)):
            child = root / label
            (child / "outputs").mkdir(parents=True, exist_ok=False)
            (child / "checkpoints").mkdir()
            subprocess.run(["mount", "--bind", str(child / "outputs"), str(target_output)], check=True)
            subprocess.run(["mount", "--bind", str(child / "checkpoints"), str(target_checkpoint)], check=True)
            command = ["/usr/bin/time", "-v", str(binary), "-q", "-y", *plan["configs"], "-d", "case"]
            try:
                result = run_until_commit(command, child, environment("A"), args.timeout_seconds)
                results.append(dict(result, label=label, binary_sha256=file_sha(binary)))
            finally:
                subprocess.run(["umount", str(target_checkpoint)], check=True)
                subprocess.run(["umount", str(target_output)], check=True)
            if not result["bundle"]:
                break
        write_json(root / "run_results.json", results)
        if len(results) == 2 and all(r["bundle"] for r in results):
            result = compare_bundle(Path(results[0]["bundle"]), Path(results[1]["bundle"]))
            write_json(root / "comparison.json", result)
            raise SystemExit(0 if result["all_exact_equal"] else 1)
        raise SystemExit(1)
    env = environment(args.stage)
    preflight = configuration_preflight(args.binary, plan["configs"], plan["pinned_input_contract"], root, env)
    plan["effective_input_preflight"] = preflight
    if not preflight["passed"]:
        plan["state"] = "CONFIGURATION_PREFLIGHT_FAILED_NOT_RUN"
        write_json(root / "stage_plan.json", plan)
        raise RuntimeError("Solver not launched: " + preflight["reason"])
    command = ["/usr/bin/time", "-v", *independent_stage_command(
        args.binary, plan["configs"], plan["pinned_input_contract"])]
    plan["command"] = command
    plan.update(state="RUNNING", binary_sha256=file_sha(args.binary), started_utc=datetime.now(timezone.utc).isoformat())
    write_json(root / "stage_plan.json", plan)
    start = time.monotonic()
    with (root / "runner.log").open("wb") as stream:
        process = subprocess.Popen(command, cwd=NATIVE_INPUT, env=env, stdout=stream,
                                   stderr=subprocess.STDOUT, start_new_session=True)
        plan["pid"] = process.pid
        write_json(root / "stage_plan.json", plan)
        try:
            process.wait(timeout=args.timeout_seconds)
        except subprocess.TimeoutExpired:
            os.killpg(process.pid, signal.SIGTERM)
            try:
                process.wait(timeout=60)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGKILL)
                process.wait()
    plan.update(state="COMPLETED_REQUIRES_SCIENTIFIC_AUDIT" if process.returncode == 0 else "EXECUTION_FAILED",
                exit_code=process.returncode, elapsed_seconds=time.monotonic() - start,
                finished_utc=datetime.now(timezone.utc).isoformat(), scientific_acceptance=False)
    write_json(root / "stage_plan.json", plan)
    raise SystemExit(0 if process.returncode == 0 else 1)


if __name__ == "__main__":
    main()
