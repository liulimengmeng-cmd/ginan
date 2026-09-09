#!/usr/bin/env python3
"""Audit external satellite-bias application telemetry in a PEA trace.

An application key is (type, epoch, receiver, satellite, signal).  Duplicate
keys would show that the same external correction entered one observable more
than once.  The audit also verifies whether product formal variance was applied
per epoch or suppressed by a product-conditioned diagnostic configuration.
"""

from __future__ import annotations

import argparse
import json
import math
import re
from collections import Counter, defaultdict
from pathlib import Path


MARKER_RE = re.compile(
    r"PPP_EXTERNAL_BIAS_APPLICATION "
    r"type=(?P<type>\S+) time=(?P<time>.*?) receiver=(?P<receiver>\S+) "
    r"satellite=(?P<satellite>\S+) signal=(?P<signal>\S+) "
    r"bias_found=(?P<bias_found>[01]) bias_m=(?P<bias_m>\S+) "
    r"product_variance_m2=(?P<product_variance_m2>\S+) "
    r"applied_variance_m2=(?P<applied_variance_m2>\S+) "
    r"variance_mode=(?P<variance_mode>\S+)"
)


def _float(value: str) -> float:
    try:
        result = float(value)
    except ValueError:
        return math.nan
    return result if math.isfinite(result) else math.nan


def audit_trace(trace: Path) -> dict[str, object]:
    trace_paths = (
        sorted(
            path
            for path in trace.glob("*.trace")
            if not path.name.startswith("Network-")
        )
        if trace.is_dir()
        else [trace]
    )
    counts: Counter[str] = Counter()
    product_stats: dict[str, dict[str, float | int | None]] = defaultdict(
        lambda: {"count": 0, "min": None, "max": None, "sum": 0.0}
    )
    seen: set[tuple[str, str, str, str, str]] = set()
    duplicate_samples: list[dict[str, str]] = []

    marker_count = 0
    malformed_marker_count = 0
    duplicate_count = 0
    nonzero_applied_variance_count = 0
    missing_bias_count = 0

    for trace_path in trace_paths:
        with trace_path.open("r", encoding="utf-8", errors="replace") as stream:
            for line in stream:
                if "PPP_EXTERNAL_BIAS_APPLICATION" not in line:
                    continue

                match = MARKER_RE.search(line)
                if match is None:
                    malformed_marker_count += 1
                    continue

                record = match.groupdict()
                marker_count += 1
                key = (
                    record["type"],
                    record["time"],
                    record["receiver"],
                    record["satellite"],
                    record["signal"],
                )
                if key in seen:
                    duplicate_count += 1
                    if len(duplicate_samples) < 20:
                        duplicate_samples.append(
                            {
                                "type": record["type"],
                                "time": record["time"],
                                "receiver": record["receiver"],
                                "satellite": record["satellite"],
                                "signal": record["signal"],
                            }
                        )
                else:
                    seen.add(key)

                label = f'{record["type"]}/{record["signal"]}/{record["variance_mode"]}'
                counts[label] += 1

                if record["bias_found"] != "1":
                    missing_bias_count += 1

                applied_variance = _float(record["applied_variance_m2"])
                if math.isfinite(applied_variance) and applied_variance > 0:
                    nonzero_applied_variance_count += 1

                product_variance = _float(record["product_variance_m2"])
                if math.isfinite(product_variance):
                    stat = product_stats[f'{record["type"]}/{record["signal"]}']
                    stat["count"] = int(stat["count"]) + 1
                    stat["sum"] = float(stat["sum"]) + product_variance
                    stat["min"] = (
                        product_variance
                        if stat["min"] is None
                        else min(float(stat["min"]), product_variance)
                    )
                    stat["max"] = (
                        product_variance
                        if stat["max"] is None
                        else max(float(stat["max"]), product_variance)
                    )

    formatted_stats: dict[str, dict[str, float | int | None]] = {}
    for label, stat in sorted(product_stats.items()):
        count = int(stat["count"])
        formatted_stats[label] = {
            "count": count,
            "min": stat["min"],
            "max": stat["max"],
            "mean": float(stat["sum"]) / count if count else None,
        }

    return {
        "schema": "GINAN_EXTERNAL_BIAS_APPLICATION_AUDIT_V1",
        "trace": str(trace),
        "trace_files": [str(path) for path in trace_paths],
        "marker_count": marker_count,
        "unique_application_key_count": len(seen),
        "duplicate_application_key_count": duplicate_count,
        "duplicate_application_key_samples": duplicate_samples,
        "malformed_marker_count": malformed_marker_count,
        "missing_external_bias_count": missing_bias_count,
        "nonzero_applied_variance_count": nonzero_applied_variance_count,
        "counts_by_type_signal_mode": dict(sorted(counts.items())),
        "product_variance_m2_by_type_signal": formatted_stats,
        "claim_limits": [
            "zero duplicate keys checks exact-once source application, not product datum compatibility",
            "CONDITIONED suppresses per-epoch product variance and reports conditional precision only",
            "this audit does not certify integer truth or PPP-AR scientific benefit",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("trace", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    report = audit_trace(args.trace)
    encoded = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(encoded, encoding="utf-8")
    else:
        print(encoded, end="")


if __name__ == "__main__":
    main()
