#!/usr/bin/env python3
"""Read-only minimum completion audit for the actual 180-station G/E/C FLOAT.

Exit zero is not enough. This checks all 241 committed epochs, paired complete
finite product covariance, and both real checkpoint envelopes/full-KF decode
through a separately built production-linked reader. It is not position/product
accuracy certification or a replacement for the coordinate quality audit.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
import hashlib
from datetime import datetime, timedelta, timezone
import json
import math
from pathlib import Path
import subprocess

import run_network180_float as network
import run_static_checkpoint_replay as replay
import run_stages as stages

FINITE_PRODUCT_FIELDS = ("clock_m", "clock_sigma_m", "phase_m", "phase_sigma_m",
    "clock_phase_covariance_m2", "correction_m", "correction_sigma_m")
REQUIRED_SECTIONS = {"pea.controller.v1", "pea.rinex_file_streams.v1", "pea.receiver_runtime.v1",
    "pea.satellite_runtime.v1", "zhang_graph_runtime", "zhang_ppp_ar_runtime",
    "zhang_ambres_runtime", "pea.rng_declaration.v1"}
SIGNALS = {"G": {"L1C", "L2W"}, "E": {"L1C", "L5Q"}, "C": {"L1P", "L5P"}}


def require(condition, why):
    if not condition:
        raise RuntimeError(why)


def expected_epochs():
    start = datetime(2024, 7, 17)
    return [(start + timedelta(seconds=30 * index)).strftime("%Y-%m-%d %H:%M:%S") for index in range(241)]


def canonical_epoch(value):
    return datetime.fromisoformat(value).strftime("%Y-%m-%d %H:%M:%S")


def validate_coordinate_inventory(actual, producer, epoch, dimension, bundle_sha):
    require(actual.get("schema") == "CHECKPOINT_COORDINATE_INVENTORY_V1" and
            actual.get("producer_sha256") == producer and actual.get("experiment_mode") == replay.MODE and
            actual.get("checkpoint_epoch") == epoch,
            "Checkpoint inventory schema/producer/mode/actual epoch mismatch")
    require(actual.get("full_x_finite") is True and actual.get("full_P_finite") is True and
            actual.get("state_dimension") == dimension and actual.get("P_entries") == dimension ** 2 and
            actual.get("bundle_sha256") == bundle_sha,
            "Full checkpoint x/P finite/dimension/hash inventory mismatch")


def audit_epoch(directory):
    doc, sums = replay.verify_atomic_bundle(directory)
    require(doc.get("products_bytes") == (directory / "products.csv").stat().st_size and
            doc.get("covariance_bytes") == (directory / "covariance.csv").stat().st_size,
            "Atomic manifest product/full-Q byte counts mismatch")
    epoch = canonical_epoch(doc["epoch"])
    target_seconds = replay.epoch_seconds(epoch)
    require(doc.get("schema_version") == 2 and set(doc.get("systems", {})) == {"GPS", "GAL", "BDS"},
            "Missing exact G/E/C per-system atomic manifest")
    for system, value in doc["systems"].items():
        require(value.get("fixed_phase_rows") == 0 and value.get("any_pppar_usable") is False,
                "Unexpected integer-authorized output in a FLOAT-only run: " + system)
    counts, satellites = Counter(), {system: set() for system in SIGNALS}
    product_labels = set()
    with (directory / "products.csv").open(newline="") as stream:
        for row in csv.DictReader(stream):
            require(row["solution"] == "FLOAT", "Non-FLOAT solution record")
            require(float(row["gpst_seconds"]) == target_seconds, "Product epoch differs from its COMMIT")
            satellite, code = row["satellite"], row["observable"]
            system = satellite[:1]
            require(system in SIGNALS and code in SIGNALS[system], "Unreviewed system/phase signal")
            for name in FINITE_PRODUCT_FIELDS:
                require(math.isfinite(float(row[name])), "Nonfinite physical product moment: " + name)
            require(all(float(row[name]) >= 0 for name in ("clock_sigma_m", "phase_sigma_m", "correction_sigma_m")),
                    "Negative product standard deviation")
            require(row.get("numeric_valid") == "1", "Published FLOAT moment is not numeric-valid")
            require(row.get("ar_valid", "0") == "0" and row.get("pppar_usable", "0") == "0",
                    "FLOAT row incorrectly authorizes AR")
            counts[system] += 1
            satellites[system].add(satellite)
            product_labels.add((satellite, "CLOCK", "NONE"))
            product_labels.add((satellite, "PHASE", code))
    require(set(counts) == set(SIGNALS), "At least one constellation has no FLOAT phase record")
    labels, pairs, covariance_rows, negative_diagonals = set(), set(), 0, 0
    with (directory / "covariance.csv").open(newline="") as stream:
        for row in csv.DictReader(stream):
            require(row["solution"] == "FLOAT" and float(row["gpst_seconds"]) == target_seconds,
                    "Covariance solution/epoch differs from COMMIT")
            value = float(row["covariance_m2"])
            require(math.isfinite(value), "Nonfinite product covariance entry")
            left = (row["row_satellite"], row["row_parameter"], row["row_observable"])
            right = (row["column_satellite"], row["column_parameter"], row["column_observable"])
            require(left[0][:1] in SIGNALS and right[0][:1] in SIGNALS, "Unknown covariance constellation")
            labels.update((left, right))
            pair = tuple(sorted((left, right)))
            require(pair not in pairs, "Duplicate unordered covariance entry")
            pairs.add(pair)
            covariance_rows += 1
            if left == right and value < -1e-12:
                negative_diagonals += 1
    require(product_labels <= labels, "Some physical product targets have no covariance variable")
    require(covariance_rows > 0 and covariance_rows == len(labels) * (len(labels) + 1) // 2,
            "Product covariance is not the complete unique upper triangle, including cross-system blocks")
    require(negative_diagonals == 0, "Negative physical product covariance diagonal")
    return {"epoch": epoch, "path": str(directory), "sha256": sums,
        "phase_rows_by_system": dict(counts), "phase_satellites_by_system": {key: sorted(value) for key, value in satellites.items()},
        "covariance_variable_count": len(labels), "full_upper_triangle_entries": covariance_rows,
        "physical_moments_finite": True,
        "continuity_metadata_policy": "Unavailable FLOAT warmup continuity fields may be nan; they are not physical moments or AR authorization"}


def audit_checkpoints(root, plan, verifier, output_dir, reader_kind="coordinate_inventory"):
    reports = []
    for epoch in network.CHECKPOINT_EPOCHS:
        directory = replay.checkpoint_at(root / "checkpoints", epoch)
        require(directory is not None, "Required checkpoint not published: " + epoch)
        manifest_path, bundle = directory / "checkpoint_manifest.json", directory / "checkpoint.bundle"
        doc = json.loads(manifest_path.read_text())
        require(doc.get("experiment_mode") == replay.MODE, "Wrong checkpoint experiment mode")
        require(doc.get("binary_sha256") == plan["binary_sha256"], "Checkpoint producer mismatch")
        require(doc.get("runtime_id") and isinstance(doc.get("state_dimension"), int) and doc["state_dimension"] > 0,
                "Checkpoint runtime/state identity missing")
        require(hashlib.sha256(doc["config_text"].encode()).hexdigest() == doc["config_sha256"], "Checkpoint config SHA mismatch")
        require(hashlib.sha256(doc["input_manifest_text"].encode()).hexdigest() == doc["input_manifest_sha256"], "Checkpoint input SHA mismatch")
        names = [row["name"] for row in doc.get("sections", [])]
        require(len(names) == len(set(names)) and set(names) == REQUIRED_SECTIONS,
                "Checkpoint module section inventory incomplete")
        require(all(row["payload_bytes"] > 0 and len(row["sha256"]) == 64 for row in doc["sections"]),
                "Checkpoint section payload metadata invalid")
        log = output_dir / ("checkpoint_" + epoch.replace(" ", "_").replace(":", "") + ".log")
        inventory_json = log.with_suffix(".inventory.json")
        command = [str(verifier), str(bundle), plan["binary_sha256"]]
        if reader_kind == "coordinate_inventory":
            command += [replay.MODE, str(inventory_json)]
        elif reader_kind != "legacy_probe":
            raise RuntimeError("Unknown explicit checkpoint reader protocol")
        with log.open("w") as stream:
            completed = subprocess.run(command, stdout=stream, stderr=subprocess.STDOUT)
        text = log.read_text()
        require(completed.returncode == 0, "Independent real checkpoint reader rejected " + epoch)
        if reader_kind == "coordinate_inventory":
            actual = json.loads(inventory_json.read_text())
            validate_coordinate_inventory(actual, plan["binary_sha256"], epoch, doc["state_dimension"], stages.file_sha(bundle))
        else:
            require("PASS checkpoint envelope" in text and "full KF decode" in text,
                    "Explicit legacy-protocol probe did not confirm full decode")
        reports.append({"epoch": epoch, "directory": str(directory), "bundle_sha256": stages.file_sha(bundle),
            "manifest_sha256": stages.file_sha(manifest_path), "state_dimension": doc["state_dimension"],
            "full_core_and_section_sha_verification": True, "verifier_log": str(log),
            "verifier_log_sha256": stages.file_sha(log), "reader_kind": reader_kind,
            "coordinate_inventory_json": str(inventory_json) if reader_kind == "coordinate_inventory" else None})
    return reports


def run(root, verifier, output_dir, reader_kind="coordinate_inventory"):
    root = root.resolve(strict=True)
    plan = json.loads((root / "stage_plan.json").read_text())
    require(plan.get("expected_epoch_count") == 241 and plan.get("exit_code") == 0,
            "Need actual 241-epoch job with exit0; process start/partial completion is insufficient")
    require(plan.get("cold_start") is True and plan.get("ar_execution") is False and plan.get("all_nis") is False,
            "Run policy is not the approved cold FLOAT/NIS-off experiment")
    require(len(plan.get("station_ids", [])) == 180 and len(set(plan["station_ids"])) == 180,
            "Not the 180-station service configuration")
    require(not plan.get("supervision_stop_reasons") and not plan.get("authoritative_float_failures"),
            "Runtime supervision or final FLOAT failure was recorded")
    network.verify_pins(plan)
    verifier = verifier.resolve(strict=True)
    output_dir.mkdir(parents=True, exist_ok=False)
    records = []
    for path in sorted((root / "outputs/zhang_internal_products.csv.epochs").iterdir()):
        require(path.is_dir() and (path / "COMMIT.json").is_file(), "Uncommitted/unknown epoch output entity")
        records.append(audit_epoch(path))
    epochs = [row["epoch"] for row in records]
    require(len(epochs) == 241 and len(set(epochs)) == 241 and sorted(epochs) == expected_epochs(),
            "Atomic COMMIT epoch sequence is not the complete 241 unique 30 s GPST schedule")
    failures = []
    logs = [root / "runner.log", *(root / "outputs").glob("Network-*.TRACE")]
    for path in logs:
        if not path.is_file():
            continue
        with path.open(errors="replace") as stream:
            for number, line in enumerate(stream, 1):
                if network.authoritative_failure(line):
                    failures.append({"path": str(path), "line_number": number, "line": line.strip()})
    require(not failures, "Terminal authoritative FLOAT rollback/rejection in this run's own logs")
    checkpoints = audit_checkpoints(root, plan, verifier, output_dir, reader_kind)
    return {"schema_version": 1, "completed_utc": datetime.now(timezone.utc).isoformat(),
        "minimum_completion_contract_passed": True, "scientific_quality_accepted": False,
        "scope": "241 real complete finite GEC FLOAT product/full-Q epochs and two independently decoded complete checkpoints; not coordinate/product accuracy acceptance",
        "run_root": str(root), "producer_sha256": plan["binary_sha256"],
        "checkpoint_verifier": str(verifier), "checkpoint_verifier_sha256": stages.file_sha(verifier),
        "epochs": records, "checkpoints": checkpoints, "terminal_authoritative_failures": failures,
        "external_accuracy_and_coordinate_review": "Separate required evidence; not inferred from file integrity, finite values or an exit code"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-root", type=Path, required=True)
    parser.add_argument("--checkpoint-verifier", type=Path, required=True)
    parser.add_argument("--checkpoint-reader-kind", choices=("coordinate_inventory", "legacy_probe"),
                        default="coordinate_inventory")
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if args.output_dir.exists():
        parser.error("Audit output directory already exists; previous evidence will not be overwritten")
    try:
        report = run(args.run_root, args.checkpoint_verifier, args.output_dir, args.checkpoint_reader_kind)
    except Exception as error:
        report = {"minimum_completion_contract_passed": False, "scientific_quality_accepted": False,
                  "failure_reason": str(error), "run_root": str(args.run_root)}
        args.output_dir.mkdir(parents=True, exist_ok=True)
    stages.write_json(args.output_dir / "minimum_completion_audit.json", report)
    print(json.dumps({key: value for key, value in report.items() if key not in ("epochs", "checkpoints")}))
    raise SystemExit(0 if report["minimum_completion_contract_passed"] else 1)


if __name__ == "__main__":
    main()
