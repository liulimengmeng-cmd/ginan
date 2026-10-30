#!/usr/bin/env python3
"""Controlled real-data static G/E/C FLOAT checkpoint continuation contract.

A and B use one identical canonical configuration, runtime ID and output paths
in sequence. Every old output/checkpoint is preserved in an immutable archive
before reusing owned directories. Preparation never launches a solver.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import csv
import json
import math
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import time

import run_network180_float as network
import run_stages as stages

MODE = "E29_MULTI_GNSS_FLOAT_V2_GPS_GAL_BDS"
FIRST = "2024-07-17 00:00:30"
LAST = "2024-07-17 00:02:00"
RESUME = "2024-07-17 00:01:00"
RUNTIME = "E29-GEC-STATIC-RESTORE-20261021"
OWNERSHIP = "STATIC_CHECKPOINT_REPLAY_OWNED_V1"


def stamp():
    return datetime.now(timezone.utc).isoformat()


def inventory(root):
    rows = []
    for path in sorted(root.rglob("*")):
        if path.is_symlink():
            raise RuntimeError("Owned test outputs cannot contain symlinks: " + str(path))
        if path.is_file():
            rows.append({"relative": str(path.relative_to(root)), "bytes": path.stat().st_size,
                         "sha256": stages.file_sha(path)})
    return rows


def owner_guard(root):
    root = root.resolve(strict=True)
    marker = json.loads((root / "OWNERSHIP.json").read_text())
    if marker != {"schema": OWNERSHIP, "root": str(root)}:
        raise RuntimeError("Refuse output moves/cuts outside this uniquely owned test root")
    return root


def archive_owned_outputs(root, label):
    root = owner_guard(root)
    destination = root / label
    if destination.exists() or label not in ("continuous", "interrupted_before_resume_cut"):
        raise RuntimeError("Invalid or existing owned archive destination")
    source = root / "work"
    before = {name: inventory(source / name) for name in ("outputs", "checkpoints")}
    destination.mkdir()
    for name in ("outputs", "checkpoints"):
        path = (source / name).resolve(strict=True)
        if path.parent != source or path.name != name:
            raise RuntimeError("Resolved archive target escapes exact owned work directory")
        path.rename(destination / name)
    after = {name: inventory(destination / name) for name in ("outputs", "checkpoints")}
    if before != after:
        raise RuntimeError("Archived output identity changed during move")
    stages.write_json(destination / "ARCHIVE_SHA256.json", {"archived_utc": stamp(), "files": after,
        "policy": "All original bytes preserved before owned output reuse; no input/source/checkpoint deletion"})
    return destination


def checkpoint_at(directory, epoch):
    matches = []
    for manifest in directory.glob("*/checkpoint_manifest.json"):
        # Producer writes <final>.tmp.<steady-clock-nonce>, verifies its
        # bundle, then atomically renames the directory. Complete staging
        # files are not a published checkpoint and cannot authorize SIGINT.
        if ".tmp." in manifest.parent.name:
            continue
        doc = json.loads(manifest.read_text())
        if doc.get("epoch") == epoch and (manifest.parent / "checkpoint.bundle").is_file():
            matches.append(manifest.parent)
    if len(matches) > 1:
        raise RuntimeError("Ambiguous multiple checkpoints at " + epoch)
    return matches[0] if matches else None


def verify_atomic_bundle(path):
    expected = {}
    for line in (path / "SHA256SUMS").read_text().splitlines():
        if len(line) < 67 or line[64:66] != "  ":
            raise RuntimeError("Invalid atomic SHA manifest")
        digest, name = line[:64], line[66:]
        if name not in ("products.csv", "covariance.csv", "certificate.txt", "COMMIT.json") or name in expected:
            raise RuntimeError("Atomic hash manifest has unknown/duplicate filenames")
        expected[name] = digest
    if set(expected) != {"products.csv", "covariance.csv", "certificate.txt", "COMMIT.json"}:
        raise RuntimeError("Expected paired full product/Q atomic bundle")
    for name, digest in expected.items():
        if stages.file_sha(path / name) != digest:
            raise RuntimeError("Atomic hash mismatch: " + str(path / name))
    doc = json.loads((path / "COMMIT.json").read_text())
    if doc.get("ar") is not False:
        raise RuntimeError("AR output is forbidden in this FLOAT checkpoint test")
    return doc, expected


def epoch_seconds(epoch):
    value = datetime.fromisoformat(epoch)
    origin = datetime(1980, 1, 6)
    return (value - origin).total_seconds()


def cut_append_only_csv(source, destination, cutoff):
    # Preserve retained original text bytes rather than reformatting numbers.
    kept, removed = 0, 0
    with source.open("r", newline="") as inp, destination.open("w", newline="") as out:
        header = inp.readline()
        if not header or next(csv.reader([header]))[0] != "gpst_seconds":
            raise RuntimeError("Unrecognized append-only product/Q header")
        out.write(header)
        for line in inp:
            fields = next(csv.reader([line]))
            epoch = float(fields[0])
            if not math.isfinite(epoch):
                raise RuntimeError("Nonfinite append-only epoch")
            if epoch <= cutoff + 1e-9:
                out.write(line)
                kept += 1
            else:
                removed += 1
    return {"retained_rows": kept, "archived_future_rows_not_copied": removed,
            "source_sha256": stages.file_sha(source), "retained_sha256": stages.file_sha(destination)}


def restore_owned_prefix(root, archive, cutoff=FIRST):
    root = owner_guard(root)
    if archive.resolve() != root / "interrupted_before_resume_cut":
        raise RuntimeError("Resume prefix must originate in the exact owned interrupted archive")
    work = root / "work"
    for name in ("outputs", "checkpoints"):
        if (work / name).exists():
            raise RuntimeError("Refuse to overwrite existing active output/checkpoint directory")
        (work / name).mkdir()
    result = {"cutoff_epoch": cutoff, "archive": str(archive), "csv": {}, "retained_atomic_epochs": [],
              "future_atomic_epochs_preserved_only_in_archive": [], "numeric_values_reformatted": False}
    for source in (archive / "outputs").iterdir():
        target = work / "outputs" / source.name
        if source.name in ("zhang_internal_products.csv", "zhang_internal_product_covariance.csv") and source.is_file():
            result["csv"][source.name] = cut_append_only_csv(source, target, epoch_seconds(cutoff))
        elif source.name == "zhang_internal_products.csv.epochs" and source.is_dir():
            target.mkdir()
            for bundle in source.iterdir():
                if not bundle.is_dir() or not (bundle / "COMMIT.json").is_file():
                    raise RuntimeError("Uncommitted/unknown epoch directory in owned output archive")
                doc, hashes = verify_atomic_bundle(bundle)
                if epoch_seconds(doc["epoch"]) <= epoch_seconds(cutoff) + 1e-9:
                    shutil.copytree(bundle, target / bundle.name)
                    if inventory(bundle) != inventory(target / bundle.name):
                        raise RuntimeError("Retained atomic prefix changed")
                    result["retained_atomic_epochs"].append(doc["epoch"])
                else:
                    result["future_atomic_epochs_preserved_only_in_archive"].append(doc["epoch"])
        elif source.is_dir():
            shutil.copytree(source, target)
        elif source.is_file():
            shutil.copy2(source, target)
        else:
            raise RuntimeError("Unsupported owned output entity")
    checkpoint = checkpoint_at(archive / "checkpoints", cutoff)
    if checkpoint is None:
        raise RuntimeError("Archived interrupt checkpoint missing")
    restored = work / "checkpoints" / checkpoint.name
    shutil.copytree(checkpoint, restored)
    if inventory(checkpoint) != inventory(restored):
        raise RuntimeError("Retained source checkpoint bytes changed")
    result["source_checkpoint"] = str(restored)
    stages.write_json(root / "resume_prefix_audit.json", result)
    return restored


def prepare(binary, audit_path, root, coordinate_overlay, coordinate_review,
            verifier, timeout=1800, float_threads=12, omp_threads=8):
    if not math.isfinite(timeout) or timeout <= 0:
        raise ValueError("Timeout must be finite and positive")
    if not all(isinstance(value, int) and 1 <= value <= 1024 for value in (float_threads, omp_threads)):
        raise ValueError("Thread counts must be positive integers at most 1024")
    audit = json.loads(audit_path.read_text())
    if not audit.get("selected_complete") or audit.get("missing_files"):
        raise RuntimeError("Original 12-station input readiness failed")
    review, coordinate_value = network.coordinate_approval(coordinate_overlay, coordinate_review)
    root = root.resolve()
    root.mkdir(parents=True, exist_ok=False)
    stages.write_json(root / "OWNERSHIP.json", {"schema": OWNERSHIP, "root": str(root)})
    work = root / "work"
    for name in ("config", "outputs", "checkpoints"):
        (work / name).mkdir(parents=True)
    value = stages.overlay("D", work, audit)
    value["processing_options"]["epoch_control"]["end_epoch"] = LAST
    options = value["processing_options"]["gnss_general"]["zhang_pppar"]
    options.update(deterministic_checkpoint=True, checkpoint_runtime_id=RUNTIME,
        checkpoint_output_directory=str(work / "checkpoints"), checkpoint_capture_epochs=[FIRST, LAST],
        checkpoint_restore_path="")
    overlay = work / "config/static_checkpoint_overlay.yaml"
    stages.write_json(overlay, value)
    coord = work / "config/static_coordinate_overlay.yaml"
    stages.write_json(coord, coordinate_value)
    configs = [str(stages.FROZEN_CONFIG / name) for name in stages.CONFIG_NAMES] + [str(overlay), str(coord)]
    pinned = stages.pin_stage_inputs(audit, configs)
    binary, verifier = binary.resolve(strict=True), verifier.resolve()
    env = stages.environment("D")
    env.update(ZHANG_MULTI_GNSS_FLOAT_CHECKPOINT="1", ZHANG_R51_FLOAT_THREADS=str(float_threads),
               OMP_NUM_THREADS=str(omp_threads), ZHANG_CLOCK_LIFECYCLE_DIAGNOSTICS="1")
    plan = {"schema_version": 1, "state": "PREPARED_NOT_RUN", "run_root": str(root),
        "binary": str(binary), "binary_sha256": stages.file_sha(binary), "verifier": str(verifier),
        "expected_checkpoint_mode": MODE, "configs": configs, "pinned_input_contract": pinned,
        "expected_epoch_count": 5, "start_epoch": "2024-07-17 00:00:00", "end_epoch": LAST,
        "capture_epochs": [FIRST, LAST], "runtime_id": RUNTIME,
        "coordinate_review": review, "ratio": 2, "all_nis": False, "ar_execution": False,
        "cold_start": True, "timeout_seconds": timeout,
        "environment": {key: value for key, value in env.items() if key.startswith("ZHANG_") or key.endswith("_NUM_THREADS")},
        "command": stages.independent_stage_command(binary, configs, pinned),
        "scientific_acceptance": False,
        "comparison_policy": "Independent A continuous and fresh B resume use same paths/runtime/config; strict complete numeric/section/product equality; no value normalization",
        "preservation_policy": "Archive all owned outputs/checkpoints before directory reuse; retain exact <=00:00:30 output prefix before resume",
        "source_identity": [network.identity(Path(__file__)), network.identity(Path(network.__file__)),
                            network.identity(Path(stages.__file__)), network.identity(Path(__file__).resolve().parents[2] /
                                "src/tests/unit/multi_gnss_checkpoint_compare.cpp")],
        "external_review_files": [network.identity(coordinate_overlay), network.identity(coordinate_review), network.identity(audit_path)]}
    stages.write_json(root / "replay_plan.json", plan)
    return plan, env


def verify_plan(plan):
    owner_guard(Path(plan["run_root"]))
    stages.verify_pinned_stage_inputs(plan["pinned_input_contract"])
    if stages.file_sha(plan["binary"]) != plan["binary_sha256"]:
        raise RuntimeError("Producer binary changed after test preparation")
    for row in [*plan["source_identity"], *plan["external_review_files"]]:
        if stages.file_sha(row["path"]) != row["sha256"]:
            raise RuntimeError("Pinned test/review/input source changed")


def run_phase(plan, env, label, interrupt_at_first=False):
    root = owner_guard(Path(plan["run_root"]))
    verify_plan(plan)
    log = root / (label + ".log")
    start = time.monotonic()
    signalled, source = False, None
    with log.open("wb") as stream:
        child = subprocess.Popen(plan["command"], cwd=stages.NATIVE_INPUT, env=env,
                                 stdout=stream, stderr=subprocess.STDOUT, start_new_session=True)
        result = {"label": label, "pid": child.pid, "started_utc": stamp(), "log": str(log)}
        stages.write_json(root / (label + "_process.json"), result)
        while child.poll() is None:
            if interrupt_at_first and not signalled:
                source = checkpoint_at(root / "work/checkpoints", FIRST)
                if source:
                    os.killpg(child.pid, signal.SIGINT)
                    signalled = True
                    result["signal_sent_utc"] = stamp()
            if time.monotonic() - start > plan["timeout_seconds"]:
                os.killpg(child.pid, signal.SIGTERM)
                try:
                    child.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    os.killpg(child.pid, signal.SIGKILL)
                    child.wait()
                result["timed_out"] = True
                break
            time.sleep(.05)
        child.wait()
    result.update(exit_code=child.returncode, elapsed_seconds=time.monotonic() - start,
        finished_utc=stamp(), interrupt_sent_after_published_checkpoint=signalled,
        source_checkpoint=str(source) if source else None, log_sha256=stages.file_sha(log))
    stages.write_json(root / (label + "_process.json"), result)
    if interrupt_at_first and (not signalled or child.returncode != 130):
        raise RuntimeError("B must stop by real SIGINT at an epoch boundary with exit code 130")
    if not interrupt_at_first and child.returncode:
        raise RuntimeError(label + " numerical execution failed")
    return result


def final_atomic_comparison(left, right):
    def collect(root):
        result = {}
        for directory in (root / "zhang_internal_products.csv.epochs").iterdir():
            if not directory.is_dir():
                continue
            doc, hashes = verify_atomic_bundle(directory)
            if doc["epoch"] in result:
                raise RuntimeError("Repeated atomic product epoch")
            result[doc["epoch"]] = hashes
        return result
    first, second = collect(left), collect(right)
    expected = {"2024-07-17 00:" + minute + ":" + sec for minute, sec in
                (("00", "00"), ("00", "30"), ("01", "00"), ("01", "30"), ("02", "00"))}
    actual_first = {datetime.fromisoformat(value).strftime("%Y-%m-%d %H:%M:%S") for value in first}
    actual_second = {datetime.fromisoformat(value).strftime("%Y-%m-%d %H:%M:%S") for value in second}
    if actual_first != expected or actual_second != expected:
        raise RuntimeError("Continuous/resumed output must each commit all five physical epochs")
    return {"passed": first == second, "epoch_count": len(first),
        "left_epoch_sha256": first, "right_epoch_sha256": second,
        "scope": "Every complete products.csv/full covariance.csv/certificate/COMMIT, exact byte identity"}


def execution_gate(plan, path):
    review = json.loads(path.read_text())
    if any(review.get(key) is not True for key in
           ("build_passed", "mandatory_tests_passed", "approved_for_controlled_checkpoint_test")):
        raise RuntimeError("Controlled numerical checkpoint test requires explicit build/test approval")
    if review.get("binary_sha256") != plan["binary_sha256"]:
        raise RuntimeError("Controlled test approval binary SHA-256 mismatch")
    if not Path(plan["verifier"]).is_file():
        raise RuntimeError("Build the independent full checkpoint reader first")
    return review


def execute(plan, env, approval):
    if sys.platform != "linux":
        raise RuntimeError("Actual controlled replay runs only in native Linux")
    root = owner_guard(Path(plan["run_root"]))
    approved = execution_gate(plan, approval)
    verify_plan(plan)
    preflight = stages.configuration_preflight(Path(plan["binary"]), plan["configs"],
        plan["pinned_input_contract"], root, env)
    plan["configuration_preflight"] = preflight
    if not preflight["passed"]:
        stages.write_json(root / "replay_plan.json", plan)
        raise RuntimeError("Numerical solver not launched: " + preflight["reason"])
    plan.update(state="RUNNING", controlled_execution_approval=approved,
                verifier_sha256=stages.file_sha(plan["verifier"]))
    stages.write_json(root / "replay_plan.json", plan)
    run_phase(plan, env, "A_continuous")
    continuous = archive_owned_outputs(root, "continuous")
    if checkpoint_at(continuous / "checkpoints", LAST) is None:
        raise RuntimeError("Continuous final checkpoint missing")
    for name in ("outputs", "checkpoints"):
        (root / "work" / name).mkdir()
    run_phase(plan, env, "B_interrupted", True)
    interrupted = archive_owned_outputs(root, "interrupted_before_resume_cut")
    source = restore_owned_prefix(root, interrupted)
    resumed_env = dict(env, ZHANG_E29_RESTORE_BUNDLE=str(source))
    resumed = run_phase(plan, resumed_env, "B_fresh_resume")
    text = Path(resumed["log"]).read_text(errors="replace")
    required = ("Restored E29 checkpoint", "first resume epoch is " + RESUME,
                "E29 checkpoint resume boundary verified at epoch", "2024-07-17 00:01:00")
    if not all(marker in text for marker in required):
        raise RuntimeError("Fresh process restore/resume-next-epoch log evidence incomplete")
    final_left = checkpoint_at(continuous / "checkpoints", LAST)
    final_right = checkpoint_at(root / "work/checkpoints", LAST)
    if final_right is None:
        raise RuntimeError("Fresh resumed process failed to commit final checkpoint")
    compare_command = [plan["verifier"], str(final_left / "checkpoint.bundle"),
        str(final_right / "checkpoint.bundle"), plan["binary_sha256"], MODE, str(root / "full_checkpoint_comparison.json")]
    with (root / "checkpoint_comparison.log").open("w") as log:
        checked = subprocess.run(compare_command, stdout=log, stderr=subprocess.STDOUT)
    if checked.returncode:
        raise RuntimeError("Actual continuous versus resumed full checkpoint differs; inspect independent comparison report")
    comparison = json.loads((root / "full_checkpoint_comparison.json").read_text())
    expected_receivers = {row["station"] for row in plan["pinned_input_contract"]["observations"]}
    if set(comparison["coordinate_components_by_receiver"]) != expected_receivers:
        raise RuntimeError("Actual committed XYZ states do not belong to all selected service receivers")
    products = final_atomic_comparison(continuous / "outputs", root / "work/outputs")
    stages.write_json(root / "all_epoch_product_q_comparison.json", products)
    if not products["passed"]:
        raise RuntimeError("Actual continuous and resumed product/full-Q byte identity failed")
    if stages.file_sha(plan["verifier"]) != plan["verifier_sha256"]:
        raise RuntimeError("Test reader binary changed during controlled replay")
    verify_plan(plan)
    plan.update(state="CONTROLLED_RESTORE_CONTRACT_PASSED_NOT_180_STATION_SCIENTIFIC_ACCEPTANCE",
        controlled_restore_contract_passed=True, fresh_resume_pid=resumed["pid"],
        next_epoch_qc_commit_evidence="Real post-resume 00:01:00 atomic COMMIT/full Q plus final post-epoch checkpoint",
        finished_utc=stamp(), scientific_acceptance=False)
    stages.write_json(root / "replay_plan.json", plan)
    return True


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("binary", "audit", "run-root", "coordinate-overlay", "coordinate-review", "verifier"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--timeout-seconds", type=float, default=1800)
    parser.add_argument("--float-threads", type=int, default=12)
    parser.add_argument("--omp-threads", type=int, default=8)
    parser.add_argument("--build-approval", type=Path)
    parser.add_argument("--execute", action="store_true")
    args = parser.parse_args()
    if args.execute and args.build_approval is None:
        parser.error("--execute requires --build-approval")
    plan, env = prepare(args.binary, args.audit, args.run_root, args.coordinate_overlay,
        args.coordinate_review, args.verifier, args.timeout_seconds, args.float_threads, args.omp_threads)
    if args.execute:
        try:
            execute(plan, env, args.build_approval)
        except Exception as error:
            plan.update(state="CONTROLLED_RESTORE_CONTRACT_FAILED", failure_reason=str(error),
                        finished_utc=stamp(), controlled_restore_contract_passed=False)
            stages.write_json(Path(plan["run_root"]) / "replay_plan.json", plan)
            raise
    print(json.dumps(plan, indent=2))


if __name__ == "__main__":
    main()
