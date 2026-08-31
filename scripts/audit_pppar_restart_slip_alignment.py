#!/usr/bin/env python3
"""Relate restart-disagreed PPP-AR candidates to ambiguity-arc reset history.

This is a diagnostic association audit, not a wrong-fix classifier.  It uses
the structured PDE-CS-DIAG/PDE-SCDIA-DIAG records already emitted by Ginan and
compares the most recent reset-capable event for every satellite appearing in
a shared candidate row whose integer right-hand side differs between runs.
"""

from __future__ import annotations

import argparse
import bisect
import json
import re
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

from audit_pppar_candidate_restart_consistency import (
    absolute_half_minute_index,
    parse_candidate_groups,
)


DIAG_RE = re.compile(
    r"(?P<kind>PDE-CS-DIAG|PDE-SCDIA-DIAG)\s+"
    r"(?:detector=(?P<detector>\S+)\s+)?action=(?P<action>\S+)\s+"
    r"epoch=(?P<epoch>.*?)\s+rec=(?P<receiver>\S+)\s+sat=(?P<satellite>\S+).*?"
    r"reason=(?P<reason>\S+)"
)
TERM_RE = re.compile(r"A\(([^,]+),([^,]+),([^)]+)\)")
RESET_ACTIONS = {
    "PDE-CS-DIAG": {"detected", "retracking", "flagged"},
    "PDE-SCDIA-DIAG": {"detected"},
}


def parse_utc(text: str) -> datetime:
    normalized = text.strip().replace("Z", "+00:00")
    try:
        value = datetime.fromisoformat(normalized)
    except ValueError:
        value = datetime.strptime(normalized, "%Y-%m-%d %H:%M:%S.%f")
    if value.tzinfo is None:
        value = value.replace(tzinfo=timezone.utc)
    return value.astimezone(timezone.utc)


def parse_reset_events(path: Path) -> dict[tuple[str, str], list[dict[str, object]]]:
    events: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    with path.open("r", encoding="utf-8", errors="replace") as stream:
        for line in stream:
            match = DIAG_RE.search(line)
            if not match or match.group("action") not in RESET_ACTIONS[match.group("kind")]:
                continue
            detector = match.group("detector") or "SCDIA"
            event = {
                "epoch": absolute_half_minute_index(match.group("epoch")),
                "detector": detector,
                "action": match.group("action"),
                "reason": match.group("reason"),
            }
            events[(match.group("receiver"), match.group("satellite"))].append(event)
    for rows in events.values():
        rows.sort(key=lambda row: int(row["epoch"]))
    return dict(events)


def absolute_groups(groups, day_start_index: int):
    result = {}
    for (epoch, receiver), rows in groups.items():
        absolute_epoch = epoch if epoch > 1_000_000 else day_start_index + epoch - 1
        result[(absolute_epoch, receiver)] = rows
    return result


def last_event(
    events: dict[tuple[str, str], list[dict[str, object]]],
    receiver: str,
    satellite: str,
    epoch: int,
    run_start: int,
) -> dict[str, object]:
    rows = events.get((receiver, satellite), [])
    positions = [int(row["epoch"]) for row in rows]
    index = bisect.bisect_right(positions, epoch) - 1
    if index < 0 or int(rows[index]["epoch"]) < run_start:
        return {"epoch": run_start, "detector": "RUN_START", "action": "initialised"}
    return rows[index]


def audit(
    primary_trace: Path,
    restart_trace: Path,
    day_start: datetime,
    restart_epoch_offset: int,
) -> dict[str, object]:
    day_start_index = round(day_start.timestamp() / 30)
    primary = absolute_groups(parse_candidate_groups(primary_trace), day_start_index)
    restart = absolute_groups(
        parse_candidate_groups(restart_trace, restart_epoch_offset), day_start_index
    )
    primary_events = parse_reset_events(primary_trace)
    restart_events = parse_reset_events(restart_trace)
    counts: Counter[str] = Counter()
    by_family: dict[str, Counter[str]] = defaultdict(Counter)
    samples: list[dict[str, object]] = []

    for key in sorted(primary.keys() & restart.keys()):
        primary_rows = {(row[0], row[1]): row for row in primary[key]}
        restart_rows = {(row[0], row[1]): row for row in restart[key]}
        for expression in primary_rows.keys() & restart_rows.keys():
            if abs(primary_rows[expression][2] - restart_rows[expression][2]) <= 1e-9:
                continue
            family, terms = expression
            satellites = sorted({match[1] for match in TERM_RE.findall(terms)})
            receiver = key[1]
            alignment = []
            for satellite in satellites:
                first = last_event(
                    primary_events, receiver, satellite, key[0], day_start_index
                )
                second = last_event(
                    restart_events,
                    receiver,
                    satellite,
                    key[0],
                    day_start_index + restart_epoch_offset,
                )
                alignment.append(
                    {
                        "satellite": satellite,
                        "primary_last_reset": first,
                        "restart_last_reset": second,
                        "same_reset_epoch": first["epoch"] == second["epoch"],
                    }
                )
            status = (
                "all_satellite_reset_epochs_match"
                if alignment and all(row["same_reset_epoch"] for row in alignment)
                else "one_or_more_satellite_reset_epochs_differ"
            )
            counts[status] += 1
            by_family[family][status] += 1
            if len(samples) < 50:
                samples.append(
                    {
                        "epoch": key[0],
                        "receiver": receiver,
                        "family": family,
                        "terms": terms,
                        "primary_rhs": primary_rows[expression][2],
                        "restart_rhs": restart_rows[expression][2],
                        "reset_alignment": alignment,
                        "status": status,
                    }
                )

    return {
        "schema": "GINAN_PPPAR_RESTART_SLIP_ALIGNMENT_V1",
        "primary_trace": str(primary_trace),
        "restart_trace": str(restart_trace),
        "day_start_utc": day_start.isoformat().replace("+00:00", "Z"),
        "restart_epoch_offset": restart_epoch_offset,
        "disagreed_integer_rhs_count": sum(counts.values()),
        "reset_alignment_counts": dict(sorted(counts.items())),
        "by_family": {
            family: dict(sorted(rows.items()))
            for family, rows in sorted(by_family.items())
        },
        "samples": samples,
        "claim_limits": [
            "different reset histories are an association and do not prove the integer is wrong",
            "matching reset epochs do not certify the candidate integer",
            "structured detector events do not replace independent ambiguity truth",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("primary_trace", type=Path)
    parser.add_argument("restart_trace", type=Path)
    parser.add_argument("--day-start", required=True, type=parse_utc)
    parser.add_argument("--restart-epoch-offset", required=True, type=int)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = audit(
        args.primary_trace,
        args.restart_trace,
        args.day_start,
        args.restart_epoch_offset,
    )
    payload = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(payload, encoding="utf-8")
    else:
        print(payload, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

