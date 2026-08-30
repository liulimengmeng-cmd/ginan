#!/usr/bin/env python3
"""Audit frozen Experiment 0 RINEX epoch coverage before Ginan processing."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


EXPECTED_EPOCHS = 2880
MIN_COVERAGE = 0.95
MAX_GAP_SECONDS = 300.0
REQUIRED_GPS_OBSERVABLES = {"C1C", "L1C", "C2W", "L2W"}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_rinex(path: Path) -> dict[str, object]:
    observables: set[str] = set()
    epochs: list[datetime] = []
    in_header = True
    with path.open("rt", encoding="ascii", errors="replace") as stream:
        for line in stream:
            if in_header:
                if line.startswith("G") and "SYS / # / OBS TYPES" in line:
                    observables.update(line[7:60].split())
                if "END OF HEADER" in line:
                    in_header = False
                continue
            if not line.startswith(">"):
                continue
            fields = line.split()
            if len(fields) < 8 or int(fields[7]) not in (0, 1):
                continue
            second = float(fields[6])
            whole_second = int(second)
            microsecond = round((second - whole_second) * 1_000_000)
            epochs.append(
                datetime(
                    int(fields[1]),
                    int(fields[2]),
                    int(fields[3]),
                    int(fields[4]),
                    int(fields[5]),
                    whole_second,
                    microsecond,
                    tzinfo=timezone.utc,
                )
            )

    unique_epochs = sorted(set(epochs))
    intervals = [
        (right - left).total_seconds()
        for left, right in zip(unique_epochs, unique_epochs[1:])
    ]
    interval_counts = Counter(intervals)
    max_gap = max(intervals, default=0.0)
    missing_observables = sorted(REQUIRED_GPS_OBSERVABLES - observables)
    coverage = len(unique_epochs) / EXPECTED_EPOCHS
    passed = (
        coverage >= MIN_COVERAGE
        and max_gap <= MAX_GAP_SECONDS
        and not missing_observables
        and len(epochs) == len(unique_epochs)
    )
    return {
        "path": path.as_posix(),
        "epoch_count": len(epochs),
        "unique_epoch_count": len(unique_epochs),
        "duplicate_epoch_count": len(epochs) - len(unique_epochs),
        "coverage_fraction": coverage,
        "first_epoch": unique_epochs[0].isoformat() if unique_epochs else None,
        "last_epoch": unique_epochs[-1].isoformat() if unique_epochs else None,
        "max_gap_seconds": max_gap,
        "interval_histogram_seconds": {
            str(interval): count for interval, count in sorted(interval_counts.items())
        },
        "gps_observables": sorted(observables),
        "missing_required_gps_observables": missing_observables,
        "pass": passed,
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("data_root", type=Path)
    parser.add_argument("--output", type=Path)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    root = args.data_root.resolve()
    manifest_path = root / "input_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    observation_records = [
        record
        for record in manifest["files"]
        if record["role"] == "GNSS RINEX 3 observation, 30 s"
    ]
    audits: list[dict[str, object]] = []
    for record in observation_records:
        path = root / record["relative_path"]
        audit = parse_rinex(path)
        audit["relative_path"] = record["relative_path"]
        audit["station"] = str(record["source_metadata"]["siteId"]).upper()
        audit["experiment"] = record["experiment"]
        audit["sha256_matches_manifest"] = sha256_file(path) == record["sha256"]
        audits.append(audit)
        print(
            f"{audit['experiment']} {audit['station']}: "
            f"epochs={audit['unique_epoch_count']} max_gap={audit['max_gap_seconds']} "
            f"pass={audit['pass']}",
            flush=True,
        )

    expected_file_count = len(manifest["stations"]) * len(manifest["dates"])
    all_pass = (
        len(audits) == expected_file_count
        and all(audit["pass"] for audit in audits)
        and all(audit["sha256_matches_manifest"] for audit in audits)
    )
    payload = {
        "schema": "GINAN_EXPERIMENT0_2024_RINEX_QC_V1",
        "input_manifest": manifest_path.as_posix(),
        "input_manifest_sha256": sha256_file(manifest_path),
        "criteria": {
            "expected_epochs_per_day": EXPECTED_EPOCHS,
            "minimum_coverage_fraction": MIN_COVERAGE,
            "maximum_gap_seconds": MAX_GAP_SECONDS,
            "required_gps_observables": sorted(REQUIRED_GPS_OBSERVABLES),
            "duplicate_epochs_allowed": False,
        },
        "expected_file_count": expected_file_count,
        "audited_file_count": len(audits),
        "pass": all_pass,
        "files": audits,
    }
    output_path = args.output or root / "rinex_qc.json"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    partial = output_path.with_name(output_path.name + ".part")
    partial.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(partial, output_path)
    print(f"audit: {output_path}", flush=True)
    return 0 if all_pass else 1


if __name__ == "__main__":
    raise SystemExit(main())
