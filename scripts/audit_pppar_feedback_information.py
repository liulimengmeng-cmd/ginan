#!/usr/bin/env python3
"""Audit guarded PPP-AR feedback availability and one-step information gain.

The trace statistics deliberately distinguish feedback incidence from a full
ambiguity fix rate.  Rank coverage is evaluated per receiver epoch.  The
analytic shadow is a one-step linear counterfactual from the exact pre-feedback
state; it is not a persistent no-feedback trajectory.
"""

from __future__ import annotations

import argparse
import json
import math
import re
from collections import Counter, defaultdict
from pathlib import Path
from statistics import median


EPOCH_RE = re.compile(r"Epoch\s+(\d+)\s+=")
FIELD_RE = re.compile(r"([A-Za-z0-9_]+)=([^\s]+)")


def fields(line: str) -> dict[str, str]:
    return dict(FIELD_RE.findall(line))


def number(record: dict[str, str], key: str, default: float = -1) -> float:
    try:
        return float(record.get(key, default))
    except (TypeError, ValueError):
        return default


def integer(record: dict[str, str], key: str, default: int = 0) -> int:
    try:
        return int(record.get(key, default))
    except (TypeError, ValueError):
        return default


def quantile(values: list[float], probability: float) -> float | None:
    clean = sorted(value for value in values if math.isfinite(value))
    if not clean:
        return None
    position = probability * (len(clean) - 1)
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return clean[lower]
    fraction = position - lower
    return clean[lower] * (1 - fraction) + clean[upper] * fraction


def summarize(values: list[float]) -> dict[str, float | int | None]:
    clean = [value for value in values if math.isfinite(value)]
    return {
        "count": len(clean),
        "min": min(clean) if clean else None,
        "median": median(clean) if clean else None,
        "p95": quantile(clean, 0.95),
        "max": max(clean) if clean else None,
        "mean": sum(clean) / len(clean) if clean else None,
    }


def stage_failure_gate(record: dict[str, str]) -> str:
    status = record.get("status", "MISSING")
    if status in {
        "FULL_VISIBLE_GROUP_INTEGER_DATUM_CANDIDATE_UNVERIFIED",
        "FULL_SELECTED_SUBSET_INTEGER_DATUM_CANDIDATE_UNVERIFIED",
        "FULL_GROUP_INTEGER_DATUM_CANDIDATE_UNVERIFIED",
    }:
        return "CANDIDATE_SEARCH_PASSED"
    if integer(record, "wide_lane_fixed") < integer(record, "wide_lane_target"):
        return "WIDE_LANE_INCOMPLETE_" + record.get("wide_lane_status", "MISSING")
    if integer(record, "second_d2_fixed") < integer(record, "second_d2_target"):
        return "SECOND_D2_" + record.get("second_d2_last_attempt_status", "MISSING")
    return "CANDIDATE_INTEGRITY_" + status


def sustained_epoch(binary: list[int], window: int = 120, duration: int = 240) -> int | None:
    if len(binary) < window + duration - 1:
        return None
    prefix = [0]
    for value in binary:
        prefix.append(prefix[-1] + value)
    rolling_pass = [False] * len(binary)
    for end in range(window - 1, len(binary)):
        count = prefix[end + 1] - prefix[end + 1 - window]
        rolling_pass[end] = count / window >= 0.5
    for start in range(window - 1, len(binary) - duration + 1):
        if all(rolling_pass[start : start + duration]):
            return start + 1
    return None


def audit_trace(
    path: Path,
    epoch_start: int | None = None,
    epoch_end: int | None = None,
) -> dict[str, object]:
    current_epoch: int | None = None
    all_epochs: set[int] = set()
    groups: dict[int, list[dict[str, str]]] = defaultdict(list)
    stages: dict[int, list[dict[str, str]]] = defaultdict(list)
    controls: dict[int, dict[str, str]] = {}
    submissions: dict[int, int] = {}
    shadow_summaries: list[dict[str, str]] = []
    shadow_blocks: dict[str, list[dict[str, str]]] = defaultdict(list)
    realised_blocks: dict[str, list[dict[str, str]]] = defaultdict(list)
    realisation_checks: list[dict[str, str]] = []

    with path.open("r", encoding="utf-8", errors="replace") as stream:
        for line in stream:
            match = EPOCH_RE.search(line)
            if match:
                current_epoch = int(match.group(1))
                if (
                    (epoch_start is None or current_epoch >= epoch_start)
                    and (epoch_end is None or current_epoch <= epoch_end)
                ):
                    all_epochs.add(current_epoch)
                continue
            if current_epoch is None or current_epoch not in all_epochs:
                continue
            record = fields(line)
            if "PPP_AR DUAL_FREQUENCY_GROUP" in line:
                groups[current_epoch].append(record)
            elif "PPP_AR DUAL_FREQUENCY_STAGE" in line:
                stages[current_epoch].append(record)
            elif "PPP_AR DUAL_FREQUENCY_CONTROL" in line:
                controls[current_epoch] = record
            elif "PPP_AR PSEUDOOBS_SUBMISSION" in line:
                submissions[current_epoch] = integer(record, "rows")
            elif "PPP_AR FEEDBACK_SHADOW_SUMMARY" in line:
                shadow_summaries.append(record)
            elif "PPP_AR FEEDBACK_SHADOW_BLOCK" in line:
                shadow_blocks[record.get("block", "MISSING")].append(record)
            elif "PPP_AR FEEDBACK_REALISED_BLOCK" in line:
                realised_blocks[record.get("block", "MISSING")].append(record)
            elif "PPP_AR FEEDBACK_SHADOW_REALISATION" in line:
                realisation_checks.append(record)

    stage_lookup: dict[tuple[int, str, str, str], dict[str, str]] = {}
    receivers: set[str] = set()
    for epoch, records in stages.items():
        for record in records:
            receiver = record.get("receiver", "MISSING")
            receivers.add(receiver)
            stage_lookup[(epoch, receiver, record.get("system", ""), record.get("reference", ""))] = record
    for records in groups.values():
        receivers.update(record.get("receiver", "MISSING") for record in records)

    receiver_epoch: dict[tuple[str, int], dict[str, object]] = defaultdict(
        lambda: {"eligible_rank": 0, "candidate_rank": 0, "gates": []}
    )
    success_margins: list[float] = []
    ratio_margins: list[float] = []
    success_threshold = 0.9999
    ratio_threshold = 3.0
    for epoch, records in groups.items():
        for group in records:
            receiver = group.get("receiver", "MISSING")
            key = (receiver, epoch)
            if group.get("complete_graph") != "1":
                receiver_epoch[key]["gates"].append("COMMON_DUAL_FREQUENCY_GRAPH")
                continue
            eligible = 2 * integer(group, "wide_lane_rows")
            receiver_epoch[key]["eligible_rank"] += eligible
            lookup = (
                epoch,
                receiver,
                group.get("system", ""),
                group.get("reference", ""),
            )
            stage = stage_lookup.get(lookup)
            if stage is None:
                receiver_epoch[key]["gates"].append("STAGE_RECORD_MISSING")
                continue
            gate = stage_failure_gate(stage)
            receiver_epoch[key]["gates"].append(gate)
            if gate == "CANDIDATE_SEARCH_PASSED":
                receiver_epoch[key]["candidate_rank"] += integer(stage, "combined_rank")
            for field_name, threshold, destination in (
                ("wide_lane_success_rate", success_threshold, success_margins),
                ("second_d2_success_rate", success_threshold, success_margins),
                ("wide_lane_ratio", ratio_threshold, ratio_margins),
                ("second_d2_ratio", ratio_threshold, ratio_margins),
            ):
                value = number(stage, field_name)
                if value >= 0:
                    destination.append(value - threshold)

    first_gate_counts: Counter[str] = Counter()
    submitted_receiver_epochs: set[tuple[str, int]] = set()
    candidate_rank_sum = 0
    submitted_rank_sum = 0
    eligible_rank_sum = 0
    per_receiver_first: dict[str, int] = {}
    per_receiver_rank_fraction: dict[str, list[float]] = defaultdict(list)
    for receiver in sorted(receivers):
        for epoch in sorted(all_epochs):
            record = receiver_epoch[(receiver, epoch)]
            eligible = int(record["eligible_rank"])
            candidate = int(record["candidate_rank"])
            eligible_rank_sum += eligible
            candidate_rank_sum += candidate
            submitted = epoch in submissions and candidate > 0
            if submitted:
                submitted_rank_sum += candidate
                submitted_receiver_epochs.add((receiver, epoch))
                per_receiver_first.setdefault(receiver, epoch)
                first_gate_counts["PASSED_ALL_IMPLEMENTED_GATES"] += 1
            elif candidate > 0:
                first_gate_counts[
                    "POST_CANDIDATE_" + controls.get(epoch, {}).get("status", "MISSING")
                ] += 1
            else:
                gates = list(record["gates"])
                first_gate_counts[gates[0] if gates else "NO_DUAL_FREQUENCY_GROUP"] += 1
            per_receiver_rank_fraction[receiver].append(
                candidate / eligible if eligible > 0 and submitted else 0.0
            )

    total_epochs = len(all_epochs)
    station_epoch_denominator = total_epochs * len(receivers)
    epoch_feedback_incidence = len(submissions) / total_epochs if total_epochs else 0
    station_epoch_feedback_incidence = (
        len(submitted_receiver_epochs) / station_epoch_denominator
        if station_epoch_denominator else 0
    )

    receiver_summary: dict[str, object] = {}
    sorted_epochs = sorted(all_epochs)
    for receiver in sorted(receivers):
        sequence = [
            int((receiver, epoch) in submitted_receiver_epochs)
            for epoch in sorted_epochs
        ]
        sustained_offset = sustained_epoch(sequence)
        sustained = (
            sorted_epochs[sustained_offset - 1]
            if sustained_offset is not None and sustained_offset <= len(sorted_epochs)
            else None
        )
        fractions = per_receiver_rank_fraction[receiver]
        receiver_summary[receiver] = {
            "first_feedback_epoch": per_receiver_first.get(receiver),
            "first_feedback_minutes_after_start": (
                (per_receiver_first[receiver] - min(all_epochs)) * 0.5
                if receiver in per_receiver_first and all_epochs else None
            ),
            "feedback_epoch_count": sum(sequence),
            "feedback_epoch_incidence": sum(sequence) / total_epochs if total_epochs else 0,
            "first_sustained_feedback_epoch": sustained,
            "first_rank_coverage_25_percent_epoch": next(
                (epoch for epoch, value in zip(sorted_epochs, fractions) if value >= 0.25),
                None,
            ),
            "first_rank_coverage_50_percent_epoch": next(
                (epoch for epoch, value in zip(sorted_epochs, fractions) if value >= 0.50),
                None,
            ),
        }

    block_summary: dict[str, object] = {}
    for block in sorted(set(shadow_blocks) | set(realised_blocks)):
        predicted = shadow_blocks.get(block, [])
        realised = realised_blocks.get(block, [])
        block_summary[block] = {
            "predicted_relative_trace_gain": summarize(
                [
                    number(record, "relative_trace_gain")
                    for record in predicted
                    if number(record, "relative_trace_gain") >= 0
                ]
            ),
            "coupled_predicted_relative_trace_gain": summarize(
                [
                    number(record, "coupled_relative_trace_gain")
                    for record in predicted
                    if number(record, "coupled_relative_trace_gain") >= 0
                ]
            ),
            "predicted_shift_norm": summarize(
                [number(record, "predicted_shift_norm") for record in predicted]
            ),
            "realised_relative_trace_gain": summarize(
                [
                    number(record, "relative_trace_gain")
                    for record in realised
                    if number(record, "relative_trace_gain") >= 0
                ]
            ),
            "coupled_realised_relative_trace_gain": summarize(
                [
                    number(record, "coupled_relative_trace_gain")
                    for record in realised
                    if number(record, "coupled_relative_trace_gain") >= 0
                ]
            ),
            "realised_shift_norm": summarize(
                [number(record, "realised_shift_norm") for record in realised]
            ),
        }

    valid_shadow_count = sum(
        record.get("status") == "VALID_LINEAR_SHADOW" for record in shadow_summaries
    )
    return {
        "schema": "GINAN_PPPAR_FEEDBACK_INFORMATION_AUDIT_V1",
        "trace": str(path),
        "epoch_window": {
            "requested_start": epoch_start,
            "requested_end": epoch_end,
            "observed_start": min(all_epochs) if all_epochs else None,
            "observed_end": max(all_epochs) if all_epochs else None,
        },
        "epoch_count": total_epochs,
        "receiver_count": len(receivers),
        "receivers": sorted(receivers),
        "feedback_epoch_count": len(submissions),
        "feedback_epoch_incidence": epoch_feedback_incidence,
        "station_epoch_denominator": station_epoch_denominator,
        "station_feedback_epoch_count": len(submitted_receiver_epochs),
        "station_epoch_feedback_incidence": station_epoch_feedback_incidence,
        "eligible_rank_epoch_sum": eligible_rank_sum,
        "candidate_rank_epoch_sum": candidate_rank_sum,
        "submitted_rank_epoch_sum": submitted_rank_sum,
        "candidate_rank_coverage": candidate_rank_sum / eligible_rank_sum if eligible_rank_sum else 0,
        "submitted_rank_coverage": submitted_rank_sum / eligible_rank_sum if eligible_rank_sum else 0,
        "submitted_rank_epoch_integral_hours": submitted_rank_sum * 30 / 3600,
        "first_failure_gate_counts": dict(sorted(first_gate_counts.items())),
        "success_rate_margin": summarize(success_margins),
        "ratio_margin": summarize(ratio_margins),
        "joint_nis": summarize([number(record, "joint_nis") for record in shadow_summaries]),
        "nis_per_row": summarize([number(record, "nis_per_row") for record in shadow_summaries]),
        "nis_gate_configured": False,
        "shadow_record_count": len(shadow_summaries),
        "valid_shadow_record_count": valid_shadow_count,
        "shadow_realisation_check_count": len(realisation_checks),
        "shadow_state_update_error_norm": summarize(
            [number(record, "state_update_error_norm") for record in realisation_checks]
        ),
        "shadow_covariance_error_norm": summarize(
            [number(record, "covariance_error_norm") for record in realisation_checks]
        ),
        "state_block_information": block_summary,
        "receiver_summary": receiver_summary,
        "claim_limits": [
            "feedback incidence is not a full ambiguity fix rate",
            "rank coverage is a per-epoch sum and does not certify independent rank accumulated across changing ambiguity arcs",
            "analytic one-step shadow quantifies the immediate linear update, not a persistent no-feedback trajectory",
            "joint NIS is diagnostic only because no NIS rejection gate is configured",
            "this audit does not certify integer truth, fixed STEC, or external accuracy",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("trace", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--epoch-start", type=int)
    parser.add_argument("--epoch-end", type=int)
    args = parser.parse_args()
    report = audit_trace(args.trace, args.epoch_start, args.epoch_end)
    payload = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(payload, encoding="utf-8")
    else:
        print(payload, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
