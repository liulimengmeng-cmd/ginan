#!/usr/bin/env python3
"""Compare full selected-subset PPP-AR candidates between two trace runs."""

from __future__ import annotations

import argparse
import json
import re
import statistics
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path


EPOCH_RE = re.compile(r"Epoch\s+(\d+)\s+=")
FIELD_RE = re.compile(r"([A-Za-z0-9_]+)=([^\s]+)")
RTS_STATE_TIME_RE = re.compile(r"^\*\t[^\t]*\t([^\t]+)\t")


CandidateRow = tuple[str, str, float, float, float, float]
GroupKey = tuple[int, str]


def absolute_half_minute_index(text: str) -> int:
    normalized = text.strip().replace("Z", "+00:00")
    try:
        epoch = datetime.fromisoformat(normalized)
    except ValueError:
        epoch = datetime.strptime(normalized, "%Y-%m-%d %H:%M:%S.%f")
    if epoch.tzinfo is None:
        epoch = epoch.replace(tzinfo=timezone.utc)
    return round(epoch.timestamp() / 30)


def parse_candidate_groups(path: Path, epoch_offset: int = 0) -> dict[GroupKey, set[CandidateRow]]:
    current_epoch: int | None = None
    groups: dict[GroupKey, set[CandidateRow]] = defaultdict(set)
    with path.open("r", encoding="utf-8", errors="replace") as stream:
        for line in stream:
            epoch_match = EPOCH_RE.search(line)
            if epoch_match:
                current_epoch = int(epoch_match.group(1)) + epoch_offset
                continue
            rts_time_match = RTS_STATE_TIME_RE.match(line)
            if rts_time_match:
                current_epoch = absolute_half_minute_index(rts_time_match.group(1))
            if (
                current_epoch is None
                or "PPP_AR DUAL_FREQUENCY_CANDIDATE_ROW" not in line
                or "candidate_scope=SELECTED_SUBSET" not in line
                or "status=INTEGER_MAPPED" not in line
            ):
                continue
            record = dict(FIELD_RE.findall(line))
            if "terms=" not in line or "action=" not in line:
                continue
            terms = line.split("terms=", 1)[1].rsplit("action=", 1)[0].strip()
            receiver = record.get("receiver", "MISSING")
            family = record.get("family", "MISSING")
            rhs = float(record["rhs"])
            groups[(current_epoch, receiver)].add(
                (
                    family,
                    terms,
                    rhs,
                    float(record.get("float_value", "nan")),
                    float(record.get("float_minus_integer", "nan")),
                    float(record.get("formal_sigma", "nan")),
                )
            )
    return dict(groups)


def compare_groups(
    primary: dict[GroupKey, set[CandidateRow]],
    restart: dict[GroupKey, set[CandidateRow]],
) -> dict[str, object]:
    shared_keys = sorted(primary.keys() & restart.keys())
    exact_group_count = 0
    differing_group_count = 0
    shared_row_count = 0
    agreed_rhs_count = 0
    disagreed_rhs_count = 0
    maximum_rhs_difference = 0.0
    by_receiver: dict[str, Counter[str]] = defaultdict(Counter)
    disagreement_samples: list[dict[str, object]] = []
    disagreement_epochs: set[int] = set()
    primary_disagreement_residuals: list[float] = []
    restart_disagreement_residuals: list[float] = []
    primary_disagreement_sigmas: list[float] = []
    restart_disagreement_sigmas: list[float] = []

    for key in shared_keys:
        receiver = key[1]
        primary_rows = primary[key]
        restart_rows = restart[key]
        primary_signatures = {row[:3] for row in primary_rows}
        restart_signatures = {row[:3] for row in restart_rows}
        if primary_signatures == restart_signatures:
            exact_group_count += 1
            by_receiver[receiver]["exact_group_count"] += 1
        else:
            differing_group_count += 1
            by_receiver[receiver]["differing_group_count"] += 1

        primary_by_expression = {
            (row[0], row[1]): row for row in primary_rows
        }
        restart_by_expression = {
            (row[0], row[1]): row for row in restart_rows
        }
        shared_expressions = primary_by_expression.keys() & restart_by_expression.keys()
        shared_row_count += len(shared_expressions)
        by_receiver[receiver]["shared_row_count"] += len(shared_expressions)
        for expression in shared_expressions:
            difference = abs(
                primary_by_expression[expression][2]
                - restart_by_expression[expression][2]
            )
            maximum_rhs_difference = max(maximum_rhs_difference, difference)
            if difference <= 1e-9:
                agreed_rhs_count += 1
                by_receiver[receiver]["agreed_rhs_count"] += 1
            else:
                disagreed_rhs_count += 1
                by_receiver[receiver]["disagreed_rhs_count"] += 1
                disagreement_epochs.add(key[0])
                primary_row = primary_by_expression[expression]
                restart_row = restart_by_expression[expression]
                if primary_row[4] == primary_row[4]:
                    primary_disagreement_residuals.append(abs(primary_row[4]))
                if restart_row[4] == restart_row[4]:
                    restart_disagreement_residuals.append(abs(restart_row[4]))
                if primary_row[5] == primary_row[5]:
                    primary_disagreement_sigmas.append(primary_row[5])
                if restart_row[5] == restart_row[5]:
                    restart_disagreement_sigmas.append(restart_row[5])
                if len(disagreement_samples) < 50:
                    disagreement_samples.append(
                        {
                            "epoch": key[0],
                            "receiver": receiver,
                            "family": expression[0],
                            "terms": expression[1],
                            "primary_rhs": primary_row[2],
                            "restart_rhs": restart_row[2],
                            "difference": (
                                restart_row[2] - primary_row[2]
                            ),
                            "primary_float_value": primary_row[3],
                            "restart_float_value": restart_row[3],
                            "primary_float_minus_integer": primary_row[4],
                            "restart_float_minus_integer": restart_row[4],
                            "primary_formal_sigma": primary_row[5],
                            "restart_formal_sigma": restart_row[5],
                        }
                    )

    primary_only = primary.keys() - restart.keys()
    restart_only = restart.keys() - primary.keys()
    for _, receiver in primary_only:
        by_receiver[receiver]["primary_only_group_count"] += 1
    for _, receiver in restart_only:
        by_receiver[receiver]["restart_only_group_count"] += 1

    return {
        "schema": "GINAN_PPPAR_CANDIDATE_RESTART_CONSISTENCY_V1",
        "primary_candidate_group_count": len(primary),
        "restart_candidate_group_count": len(restart),
        "shared_candidate_group_count": len(shared_keys),
        "exact_candidate_group_count": exact_group_count,
        "differing_candidate_group_count": differing_group_count,
        "primary_only_candidate_group_count": len(primary_only),
        "restart_only_candidate_group_count": len(restart_only),
        "shared_candidate_row_count": shared_row_count,
        "agreed_integer_rhs_count": agreed_rhs_count,
        "disagreed_integer_rhs_count": disagreed_rhs_count,
        "maximum_integer_rhs_difference": maximum_rhs_difference,
        "disagreement_epoch_count": len(disagreement_epochs),
        "first_disagreement_epoch": (
            min(disagreement_epochs) if disagreement_epochs else None
        ),
        "last_disagreement_epoch": (
            max(disagreement_epochs) if disagreement_epochs else None
        ),
        "disagreement_samples": disagreement_samples,
        "primary_disagreement_abs_residual_median": (
            statistics.median(primary_disagreement_residuals)
            if primary_disagreement_residuals
            else None
        ),
        "restart_disagreement_abs_residual_median": (
            statistics.median(restart_disagreement_residuals)
            if restart_disagreement_residuals
            else None
        ),
        "primary_disagreement_formal_sigma_median": (
            statistics.median(primary_disagreement_sigmas)
            if primary_disagreement_sigmas
            else None
        ),
        "restart_disagreement_formal_sigma_median": (
            statistics.median(restart_disagreement_sigmas)
            if restart_disagreement_sigmas
            else None
        ),
        "all_shared_integer_rows_agree": (
            shared_row_count > 0 and disagreed_rhs_count == 0
        ),
        "by_receiver": {
            receiver: dict(sorted(counts.items()))
            for receiver, counts in sorted(by_receiver.items())
        },
        "claim_limits": [
            "agreement is conditional on candidate rows shared by both independent runs",
            "restart consistency does not prove that the shared integer is true",
            "this comparison does not certify filter feedback or coordinate accuracy",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("primary_trace", type=Path)
    parser.add_argument("restart_trace", type=Path)
    parser.add_argument("--restart-epoch-offset", type=int, default=0)
    parser.add_argument("--output", type=Path)
    arguments = parser.parse_args()
    report = compare_groups(
        parse_candidate_groups(arguments.primary_trace),
        parse_candidate_groups(
            arguments.restart_trace,
            arguments.restart_epoch_offset,
        ),
    )
    report["primary_trace"] = str(arguments.primary_trace)
    report["restart_trace"] = str(arguments.restart_trace)
    report["restart_epoch_offset"] = arguments.restart_epoch_offset
    payload = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if arguments.output:
        arguments.output.write_text(payload, encoding="utf-8")
    else:
        print(payload, end="")
    return 0 if report["all_shared_integer_rows_agree"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
