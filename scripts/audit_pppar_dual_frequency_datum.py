#!/usr/bin/env python3
"""Audit explicit dual-frequency PPP-AR datum diagnostics in a PEA trace.

This audit deliberately separates a structurally complete integer basis from a
full integer candidate.  A full [wide-lane, d2] basis is necessary but does not
certify integer truth, filter acceptance, fixed STEC, or PPP-AR accuracy.
"""

from __future__ import annotations

import argparse
import json
import re
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path


EPOCH_RE = re.compile(r"Epoch\s+(\d+)\s+=")
FIELD_RE = re.compile(r"([A-Za-z0-9_]+)=([^\s]+)")
RTS_STATE_TIME_RE = re.compile(r"^\*\t[^\t]*\t([^\t]+)\t")


def absolute_half_minute_index(text: str) -> int:
    normalized = text.strip().replace("Z", "+00:00")
    try:
        epoch = datetime.fromisoformat(normalized)
    except ValueError:
        epoch = datetime.strptime(normalized, "%Y-%m-%d %H:%M:%S.%f")
    if epoch.tzinfo is None:
        epoch = epoch.replace(tzinfo=timezone.utc)
    return round(epoch.timestamp() / 30)


def fields(line: str) -> dict[str, str]:
    return dict(FIELD_RE.findall(line))


def as_int(record: dict[str, str], key: str) -> int:
    return int(record[key])


def longest_consecutive_epoch_run(epochs: list[int]) -> int:
    if not epochs:
        return 0
    longest = 1
    current = 1
    for previous, epoch in zip(epochs, epochs[1:]):
        if epoch == previous + 1:
            current += 1
            longest = max(longest, current)
        elif epoch != previous:
            current = 1
    return longest


def audit_trace(path: Path) -> dict[str, object]:
    current_epoch: int | None = None
    basis_by_epoch: dict[int, dict[str, str]] = {}
    summary_by_epoch: dict[int, dict[str, str]] = {}
    control_by_epoch: dict[int, dict[str, str]] = {}
    group_records: dict[int, list[dict[str, str]]] = defaultdict(list)
    stage_records: dict[int, list[dict[str, str]]] = defaultdict(list)
    candidate_records: list[dict[str, str]] = []
    action_violations: list[dict[str, object]] = []
    safety_field_violations: list[dict[str, object]] = []
    pseudoobs_submissions: list[dict[str, object]] = []

    with path.open("r", encoding="utf-8", errors="replace") as stream:
        for line_number, line in enumerate(stream, 1):
            epoch_match = EPOCH_RE.search(line)
            if epoch_match:
                current_epoch = int(epoch_match.group(1))
                continue
            rts_time_match = RTS_STATE_TIME_RE.match(line)
            if rts_time_match:
                current_epoch = absolute_half_minute_index(rts_time_match.group(1))
            if current_epoch is not None and "PPP_AR PSEUDOOBS_SUBMISSION" in line:
                submission = fields(line)
                pseudoobs_submissions.append(
                    {
                        "epoch": current_epoch,
                        "line": line_number,
                        "rows": int(submission.get("rows", "0")),
                        "status": submission.get("status"),
                    }
                )
            if "PPP_AR DUAL_FREQUENCY_" not in line or current_epoch is None:
                continue
            record = fields(line)
            if record.get("action") not in {
                "PROBE_ONLY_NOT_SUBMITTED",
                "INCOMPLETE_GRAPH_NOT_PROBED",
            }:
                action_violations.append(
                    {
                        "epoch": current_epoch,
                        "line": line_number,
                        "action": record.get("action"),
                    }
                )
            if "DUAL_FREQUENCY_CONTROL" in line:
                control_by_epoch[current_epoch] = record
            elif "DUAL_FREQUENCY_BASIS" in line:
                basis_by_epoch[current_epoch] = record
            elif "DUAL_FREQUENCY_GROUP" in line:
                group_records[current_epoch].append(record)
            elif "DUAL_FREQUENCY_STAGE" in line:
                stage_records[current_epoch].append(record)
            elif "DUAL_FREQUENCY_CANDIDATE_ROW" in line:
                candidate_records.append(record)
            elif "DUAL_FREQUENCY_DATUM_SUMMARY" in line:
                summary_by_epoch[current_epoch] = record
                if (
                    record.get("wrong_fix_certified") != "0"
                    or record.get("filter_feedback") != "0"
                ):
                    safety_field_violations.append(
                        {
                            "epoch": current_epoch,
                            "line": line_number,
                            "wrong_fix_certified": record.get(
                                "wrong_fix_certified"
                            ),
                            "filter_feedback": record.get("filter_feedback"),
                        }
                    )

    audited_epochs = sorted(set(basis_by_epoch) | set(summary_by_epoch))
    structural_pass_epochs: list[int] = []
    structural_failures: list[dict[str, object]] = []
    target_ranks: list[int] = []
    for epoch in audited_epochs:
        basis = basis_by_epoch.get(epoch)
        if not basis:
            continue
        expected_rank = as_int(basis, "expected_rank")
        actual_rank = as_int(basis, "actual_rank")
        target_ranks.append(expected_rank)
        if (
            basis.get("status") == "FULL_DUAL_FREQUENCY_INTEGER_BASIS"
            and actual_rank == expected_rank
            and basis.get("integer_valued") == "1"
            and basis.get("full_row_rank") == "1"
            and basis.get("covers_all") == "1"
            and as_int(basis, "unmatched_ambiguities") == 0
            and as_int(basis, "incomplete_groups") == 0
        ):
            structural_pass_epochs.append(epoch)
        else:
            structural_failures.append(
                {
                    "epoch": epoch,
                    "status": basis.get("status"),
                    "expected_rank": expected_rank,
                    "actual_rank": actual_rank,
                    "incomplete_groups": as_int(basis, "incomplete_groups"),
                    "unmatched_ambiguities": as_int(
                        basis,
                        "unmatched_ambiguities",
                    ),
                    "covers_all": basis.get("covers_all"),
                }
            )

    full_visible_candidate_epochs = [
        epoch
        for epoch, record in summary_by_epoch.items()
        if record.get("status") in {
            "FULL_INTEGER_DATUM_CANDIDATE_UNVERIFIED",
            "FULL_VISIBLE_INTEGER_DATUM_CANDIDATE_UNVERIFIED",
        }
        and as_int(record, "full_candidate_rows") == as_int(record, "target_rank")
        and as_int(record, "full_candidate_groups") == as_int(record, "target_groups")
    ]
    full_selected_candidate_epochs = [
        epoch
        for epoch, record in summary_by_epoch.items()
        if record.get("status")
        == "FULL_SELECTED_SUBSET_INTEGER_DATUM_CANDIDATE_UNVERIFIED"
        and int(record.get("selected_candidate_groups", "-1"))
        == as_int(record, "target_groups")
        and int(record.get("selected_candidate_rows", "0")) > 0
    ]
    full_candidate_epochs = sorted(
        set(full_visible_candidate_epochs) | set(full_selected_candidate_epochs)
    )
    second_family_reached = [
        (epoch, record)
        for epoch, records in stage_records.items()
        for record in records
        if record.get("second_d2_status") not in {None, "NOT_RUN"}
    ]
    full_group_candidates = [
        (epoch, record)
        for epoch, records in stage_records.items()
        for record in records
        if record.get("status") in {
            "FULL_GROUP_INTEGER_DATUM_CANDIDATE_UNVERIFIED",
            "FULL_VISIBLE_GROUP_INTEGER_DATUM_CANDIDATE_UNVERIFIED",
            "FULL_SELECTED_SUBSET_INTEGER_DATUM_CANDIDATE_UNVERIFIED",
        }
    ]
    selected_subset_candidates = [
        (epoch, record)
        for epoch, records in stage_records.items()
        for record in records
        if record.get("status")
        == "FULL_SELECTED_SUBSET_INTEGER_DATUM_CANDIDATE_UNVERIFIED"
    ]
    selected_subset_by_receiver: dict[str, dict[str, object]] = {}
    selected_candidates_by_receiver: dict[str, list[tuple[int, dict[str, str]]]] = (
        defaultdict(list)
    )
    for epoch, record in selected_subset_candidates:
        selected_candidates_by_receiver[record.get("receiver", "MISSING")].append(
            (epoch, record)
        )
    for receiver, candidates in sorted(selected_candidates_by_receiver.items()):
        epochs = sorted({epoch for epoch, _ in candidates})
        selected_counts = [
            int(record.get("selected_family_count", "0"))
            for _, record in candidates
        ]
        excluded_counts = [
            int(record.get("excluded_satellites", "0"))
            for _, record in candidates
        ]
        selected_subset_by_receiver[receiver] = {
            "candidate_epoch_count": len(epochs),
            "first_candidate_epoch": epochs[0],
            "last_candidate_epoch": epochs[-1],
            "longest_consecutive_epoch_run": longest_consecutive_epoch_run(epochs),
            "selected_family_count_min": min(selected_counts),
            "selected_family_count_max": max(selected_counts),
            "excluded_satellite_count_min": min(excluded_counts),
            "excluded_satellite_count_max": max(excluded_counts),
        }
    wide_lane_full_groups = [
        (epoch, record)
        for epoch, records in stage_records.items()
        for record in records
        if int(record.get("wide_lane_fixed", "-1"))
        == int(record.get("wide_lane_target", "-2"))
    ]
    second_family_fixed_groups = [
        (epoch, record)
        for epoch, records in stage_records.items()
        for record in records
        if int(record.get("second_d2_fixed", "0")) > 0
    ]
    second_family_fixed_row_count = sum(
        int(record["second_d2_fixed"])
        for _, record in second_family_fixed_groups
    )
    stage_status_counts = Counter(
        record.get("status", "MISSING")
        for records in stage_records.values()
        for record in records
    )
    wide_lane_status_counts = Counter(
        record.get("wide_lane_status", "MISSING")
        for records in stage_records.values()
        for record in records
    )
    second_family_status_counts = Counter(
        record.get("second_d2_status", "MISSING")
        for records in stage_records.values()
        for record in records
    )
    candidate_family_counts = Counter(
        record.get("family", "MISSING") for record in candidate_records
    )
    candidate_scope_counts = Counter(
        record.get("candidate_scope", "MISSING") for record in candidate_records
    )
    invalid_candidate_row_count = sum(
        record.get("status") != "INTEGER_MAPPED" for record in candidate_records
    )

    all_epochs_have_records = bool(audited_epochs) and all(
        epoch in basis_by_epoch and epoch in summary_by_epoch
        for epoch in audited_epochs
    )
    structural_basis_pass = (
        all_epochs_have_records
        and len(structural_pass_epochs) == len(audited_epochs)
        and not action_violations
        and not safety_field_violations
    )
    float_probe_control_epochs = [
        epoch
        for epoch, record in control_by_epoch.items()
        if record.get("legacy_feedback") == "0"
        and record.get("new_subset_feedback") == "0"
        and record.get("status") == "FLOAT_STATE_PROBE_ONLY"
    ]
    probe_isolated_from_filter_feedback = (
        bool(audited_epochs)
        and sorted(float_probe_control_epochs) == audited_epochs
        and not pseudoobs_submissions
    )

    return {
        "schema": "GINAN_PPPAR_DUAL_FREQUENCY_DATUM_AUDIT_V1",
        "trace": str(path),
        "audited_epoch_count": len(audited_epochs),
        "basis_record_count": len(basis_by_epoch),
        "summary_record_count": len(summary_by_epoch),
        "control_record_count": len(control_by_epoch),
        "group_record_count": sum(map(len, group_records.values())),
        "stage_record_count": sum(map(len, stage_records.values())),
        "structural_basis_pass": structural_basis_pass,
        "structural_pass_epoch_count": len(structural_pass_epochs),
        "structural_pass_epochs": structural_pass_epochs,
        "structural_failures": structural_failures,
        "target_rank_min": min(target_ranks) if target_ranks else None,
        "target_rank_max": max(target_ranks) if target_ranks else None,
        "second_family_reached_group_count": len(second_family_reached),
        "wide_lane_full_group_count": len(wide_lane_full_groups),
        "second_family_fixed_group_count": len(second_family_fixed_groups),
        "second_family_fixed_row_count": second_family_fixed_row_count,
        "full_group_candidate_count": len(full_group_candidates),
        "selected_subset_group_candidate_count": len(selected_subset_candidates),
        "selected_subset_by_receiver": selected_subset_by_receiver,
        "full_visible_epoch_candidate_count": len(full_visible_candidate_epochs),
        "full_visible_epoch_candidate_epochs": sorted(
            full_visible_candidate_epochs
        ),
        "full_selected_epoch_candidate_count": len(
            full_selected_candidate_epochs
        ),
        "full_selected_epoch_candidate_epochs": sorted(
            full_selected_candidate_epochs
        ),
        "full_epoch_candidate_count": len(full_candidate_epochs),
        "full_epoch_candidate_epochs": sorted(full_candidate_epochs),
        "candidate_row_count": len(candidate_records),
        "candidate_family_counts": dict(sorted(candidate_family_counts.items())),
        "candidate_scope_counts": dict(sorted(candidate_scope_counts.items())),
        "invalid_candidate_row_count": invalid_candidate_row_count,
        "stage_status_counts": dict(sorted(stage_status_counts.items())),
        "wide_lane_status_counts": dict(sorted(wide_lane_status_counts.items())),
        "second_family_status_counts": dict(
            sorted(second_family_status_counts.items())
        ),
        "action_violations": action_violations,
        "safety_field_violations": safety_field_violations,
        "pseudoobs_submission_count": len(pseudoobs_submissions),
        "pseudoobs_submitted_row_count": sum(
            int(record["rows"]) for record in pseudoobs_submissions
        ),
        "float_probe_control_epoch_count": len(float_probe_control_epochs),
        "probe_isolated_from_filter_feedback": probe_isolated_from_filter_feedback,
        "safety": {
            "diagnostic_only": (
                not action_violations
                and not safety_field_violations
                and not pseudoobs_submissions
            ),
            "filter_feedback_certified": False,
            "wrong_fix_certified": False,
        },
        "claim_limits": [
            "structural_basis_pass proves only an explicit integer-valued full-row-rank [wide-lane,d2] basis over the eligible dual-frequency graph",
            "a full visible or selected-subset integer datum candidate does not prove integer truth or filter acceptance",
            "this audit does not certify wrong-fix rejection, fixed STEC, coordinate accuracy, or scientific PPP-AR acceptance",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("trace", type=Path)
    parser.add_argument("--output", type=Path)
    arguments = parser.parse_args()
    report = audit_trace(arguments.trace)
    payload = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if arguments.output:
        arguments.output.write_text(payload, encoding="utf-8")
    else:
        print(payload, end="")
    return 0 if report["structural_basis_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
