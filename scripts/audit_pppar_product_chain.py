#!/usr/bin/env python3
"""Audit one SP3/CLK/Bias-SINEX chain before a controlled PPP-AR run."""

from __future__ import annotations

import argparse
import hashlib
import math
from collections import defaultdict
from datetime import date
from pathlib import Path

from audit_experiment0_2024_products import (
    CLK_INTERVAL_SECONDS,
    REQUIRED_OBSERVABLES,
    SP3_INTERVAL_SECONDS,
    open_text,
    parse_bia,
    parse_clk,
    parse_sp3,
    summarize_bias_records,
    summarize_epoch_records,
    write_json_atomic,
)


SCHEMA = "GINAN_PPPAR_PRODUCT_CHAIN_AUDIT_V1"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def scalar_summary(values: list[float]) -> dict[str, float | int | None]:
    if not values:
        return {"count": 0, "min": None, "max": None, "mean": None}
    return {
        "count": len(values),
        "min": min(values),
        "max": max(values),
        "mean": sum(values) / len(values),
    }


def parse_osb_values(path: Path) -> dict[str, object]:
    estimates_ns: dict[str, list[float]] = defaultdict(list)
    sigmas_ns: dict[str, list[float]] = defaultdict(list)
    phase_pairs: dict[str, dict[str, float]] = defaultdict(dict)

    with open_text(path) as stream:
        for line in stream:
            if not line.startswith(" OSB"):
                continue
            tokens = line.split()
            if len(tokens) < 9 or not tokens[2].startswith("G"):
                continue
            observable = tokens[3]
            if observable not in REQUIRED_OBSERVABLES:
                continue
            try:
                estimate = float(tokens[7])
                sigma = float(tokens[8])
            except ValueError:
                continue
            estimates_ns[observable].append(estimate)
            sigmas_ns[observable].append(sigma)
            if observable in {"L1C", "L2W"}:
                phase_pairs[tokens[2]][observable] = estimate

    paired = [
        (signals["L1C"], signals["L2W"])
        for signals in phase_pairs.values()
        if "L1C" in signals and "L2W" in signals
    ]
    relation: dict[str, float | int | None] = {
        "pair_count": len(paired),
        "correlation": None,
        "centered_slope_l2w_on_l1c": None,
    }
    if len(paired) >= 2:
        mean_x = sum(x for x, _ in paired) / len(paired)
        mean_y = sum(y for _, y in paired) / len(paired)
        sxx = sum((x - mean_x) ** 2 for x, _ in paired)
        syy = sum((y - mean_y) ** 2 for _, y in paired)
        sxy = sum((x - mean_x) * (y - mean_y) for x, y in paired)
        if sxx > 0 and syy > 0:
            relation["correlation"] = sxy / math.sqrt(sxx * syy)
            relation["centered_slope_l2w_on_l1c"] = sxy / sxx

    return {
        "unit": "ns",
        "estimate_by_observable": {
            observable: scalar_summary(estimates_ns[observable])
            for observable in REQUIRED_OBSERVABLES
        },
        "formal_sigma_by_observable": {
            observable: {
                **scalar_summary(sigmas_ns[observable]),
                "zero_count": sum(value == 0 for value in sigmas_ns[observable]),
            }
            for observable in REQUIRED_OBSERVABLES
        },
        "phase_pair_relation": relation,
    }


def audit(day: date, sp3_path: Path, clk_path: Path, bia_path: Path) -> dict[str, object]:
    sp3_records = parse_sp3(sp3_path, day)
    clk_records = parse_clk(clk_path, day)
    bia_records = parse_bia(bia_path)
    sp3_summary, sp3_complete = summarize_epoch_records(
        sp3_records, day, SP3_INTERVAL_SECONDS
    )
    clk_summary, clk_complete = summarize_epoch_records(
        clk_records, day, CLK_INTERVAL_SECONDS
    )
    bia_summary, bia_complete = summarize_bias_records(bia_records, day)

    sp3_available = {
        satellite
        for satellite, summary in sp3_summary.items()
        if int(summary["unique_valid_epoch_count"]) > 0
    }
    clk_available = {
        satellite
        for satellite, summary in clk_summary.items()
        if int(summary["unique_valid_epoch_count"]) > 0
    }
    eligible = sp3_available & clk_available & bia_complete
    full_day_complete = sp3_complete & clk_complete & bia_complete
    incomplete = sorted(eligible - full_day_complete)

    return {
        "schema": SCHEMA,
        "date": day.isoformat(),
        "required_observables": list(REQUIRED_OBSERVABLES),
        "files": {
            "sp3": {"path": str(sp3_path), "sha256": sha256(sp3_path)},
            "clk": {"path": str(clk_path), "sha256": sha256(clk_path)},
            "bia": {"path": str(bia_path), "sha256": sha256(bia_path)},
        },
        "bia_metadata": bia_records.metadata,
        "bia_values": parse_osb_values(bia_path),
        "sp3_available_gps": sorted(sp3_available),
        "sp3_complete_gps": sorted(sp3_complete),
        "clk_available_gps": sorted(clk_available),
        "clk_complete_gps": sorted(clk_complete),
        "bia_complete_gps": sorted(bia_complete),
        "eligible_gps": sorted(eligible),
        "full_day_complete_gps": sorted(full_day_complete),
        "eligible_with_coverage_defects": incomplete,
        "satellite_details": {
            satellite: {
                "sp3": sp3_summary[satellite],
                "clk": clk_summary[satellite],
                "bia": bia_summary[satellite],
            }
            for satellite in sorted(set(sp3_summary) | set(clk_summary) | set(bia_summary))
        },
        "validation": {
            "has_eligible_satellites": bool(eligible),
            "eligible_are_full_day_complete": not incomplete,
            "bias_mode_is_absolute": str(
                bia_records.metadata.get("bias_mode", "")
            ).upper()
            == "ABSOLUTE",
            "passed": bool(eligible)
            and not incomplete
            and str(bia_records.metadata.get("bias_mode", "")).upper() == "ABSOLUTE",
        },
        "claim_limits": [
            "Product presence and exact temporal coverage do not prove CLK/OSB datum compatibility.",
            "Bias-SINEX formal sigmas do not provide cross-signal or cross-satellite covariance.",
            "A zero formal sigma is recorded as published and is not independently validated here.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--date", required=True, type=date.fromisoformat)
    parser.add_argument("--sp3", required=True, type=Path)
    parser.add_argument("--clk", required=True, type=Path)
    parser.add_argument("--bia", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    payload = audit(
        args.date,
        args.sp3.resolve(),
        args.clk.resolve(),
        args.bia.resolve(),
    )
    write_json_atomic(payload, args.output)
    print(
        "passed={} eligible={} complete={}".format(
            payload["validation"]["passed"],
            ",".join(payload["eligible_gps"]),
            ",".join(payload["full_day_complete_gps"]),
        )
    )
    return 0 if payload["validation"]["passed"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
