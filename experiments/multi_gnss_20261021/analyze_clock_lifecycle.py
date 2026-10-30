#!/usr/bin/env python3
"""Read-only numeric causal-boundary evidence from named cold FLOAT trials.

No PEA process is launched and no original trial files are rewritten.
Measured values, arithmetic identities and causal interpretations are separate.
"""
from __future__ import annotations

import argparse
from collections import Counter
import json
from pathlib import Path

from audit_stage_results import audit_logs, sha256, unavailable


def read_events(root):
    path = root / "clock_lifecycle_events.jsonl"
    rows = []
    with path.open() as stream:
        for number, line in enumerate(stream, 1):
            row = json.loads(line)
            row["jsonl_line"] = number
            rows.append(row)
    return rows


def compact(row):
    event = row["event"]
    return {"runner_source_line": row["source_line"], "jsonl_line": row["jsonl_line"],
        "event": {key: value for key, value in event.items() if key != "observations"},
        "spp_valid_source_counts": dict(Counter(obs.get("source", "UNAVAILABLE")
            for obs in event.get("observations", []) if obs.get("spp_valid") is True))}


def provenance(root):
    return {name: {"path": str(root / name), "sha256": sha256(root / name)}
        for name in ("stage_plan.json", "runner.log", "clock_lifecycle_events.jsonl") if (root / name).is_file()}


def spp_result_evidence(events):
    records = [{"epoch": row["event"].get("epoch"), "receiver": row["event"].get("receiver"),
        "selected_system": row["event"].get("spp_system"), "status": row["event"].get("spp_status"),
        "runner_source_line": row["source_line"]} for row in events if row["event"]["phase"] == "SPP_RESULT"]
    by_epoch = {}
    for row in records:
        key = str(row["epoch"]) + "|" + str(row["selected_system"])
        counts = by_epoch.setdefault(key, Counter())
        counts[row["status"]] += 1
    return {"status": "AVAILABLE" if records else "UNAVAILABLE", "records": records,
        "counts_by_epoch_selected_system": {key: dict(counts) for key, counts in by_epoch.items()},
        "scope": "Explicit SPP_RESULT events, not positions/absence-of-errors inference or a full per-system marginal solution"}


def galileo_chain(root):
    events = read_events(root)
    chains = []
    for receiver in ("NNOR", "MAL2", "MATE"):
        selected = [row for row in events if row["event"].get("receiver") == receiver]
        saves = [row for row in selected if row["event"]["phase"] == "SPP_OFFSET_SAVE_AFTER"]
        results = [row for row in selected if row["event"]["phase"] == "SPP_RESULT"]
        adjustments = [row for row in selected if row["event"]["phase"] == "REC_CLOCK_ADJUST_AFTER"]
        arithmetic = []
        for row in adjustments:
            event = row["event"]
            fields = ("spp_clock_m", "saved_offset_m", "predicted_clock_m", "unrounded_adjustment_m")
            if all(event.get(field) is not None for field in fields):
                independently_recomputed = event["spp_clock_m"] + event["saved_offset_m"] - event["predicted_clock_m"]
                arithmetic.append({"epoch": event["epoch"], "runner_source_line": row["source_line"],
                    "spp_plus_offset_minus_prediction_m": independently_recomputed,
                    "logged_unrounded_adjustment_m": event["unrounded_adjustment_m"],
                    "logged_applied_adjustment_m": event.get("applied_adjustment_m"),
                    "arithmetic_residual_m": independently_recomputed - event["unrounded_adjustment_m"],
                    "same_chart_no_saved_offset_counterfactual_m": event["spp_clock_m"] - event["predicted_clock_m"],
                    "counterfactual_scope": "Arithmetic diagnostic only, not proposed receiver-jump handling"})
        chains.append({"receiver": receiver, "observed_spp_results": [compact(row) for row in results],
            "observed_saved_offsets": [compact(row) for row in saves],
            "observed_adjustment_results": [compact(row) for row in adjustments],
            "independent_scalar_arithmetic": arithmetic})
    return {"run_root": str(root), "provenance": provenance(root), "receiver_chains": chains,
        "explicit_spp_results": spp_result_evidence(events),
        "measured_conclusion": "First successful SPP results use PRECISE clocks; next selected SPP results use KALMAN clocks while first-epoch saved offsets remain and are numerically added again",
        "causal_interpretation": "This locates a concrete auxiliary offset/source-origin mismatch consistent with the observed common receiver-clock QC innovations; explicit source labels alone are not a strict complete chart-transport proof",
        "not_proved": ["Full-state covariance equivalence of a proposed repair", "No physical receiver jump", "Scientific product acceptance"]}


def bds_initialization(root):
    events = read_events(root)
    def first(phase, receiver=None, satellite=None):
        return next(row for row in events if row["event"]["phase"] == phase and
            (receiver is None or row["event"].get("receiver") == receiver) and
            (satellite is None or row["event"].get("satellite") == satellite))
    ref_spp = first("SPP_RESULT", "MAL2")
    ref_zero = first("COLD_REFERENCE_CLOCK_ZERO", "MAL2")
    receiver_seed = first("COLD_RECEIVER_CLOCK_SEED", "NNOR")
    satellites = [first("COLD_SATELLITE_CLOCK_SEED", "MAL2", sat) for sat in ("C24", "C42")]
    satellite_details = []
    for row in satellites:
        event = row["event"]
        obs = next(obs for obs in event["observations"] if obs["satellite"] == event["satellite"])
        satellite_details.append({**compact(row), "same_endpoint_observed_clock": obs,
            "seed_minus_observed_clock_m": event["init_clock_m"] - obs["clock_m"],
            "scope": "Small transmit-time/nominal clock difference retained; it is not a 102717m chart translation"})
    reference = ref_spp["event"]["spp_clock_m"]
    rec = receiver_seed["event"]["init_clock_m"]
    sat = satellites[1]["event"]["init_clock_m"]
    return {"run_root": str(root), "provenance": provenance(root),
        "explicit_spp_results": spp_result_evidence(events),
        "reference_spp": compact(ref_spp), "reference_forced_zero": compact(ref_zero),
        "nonreference_receiver_seed": compact(receiver_seed), "satellite_seeds": satellite_details,
        "independent_scalar_example": {"reference_spp_minus_model_zero_m": reference - ref_zero["event"]["model_clock_m"],
            "nnor_seed_minus_own_spp_m": rec - receiver_seed["event"]["spp_clock_m"],
            "candidate_common_translation_m": -reference,
            "nonreference_rec_minus_sat_before_m": rec - sat,
            "nonreference_rec_minus_sat_after_both_translated_m": (rec - reference) - (sat - reference),
            "reference_zero_minus_current_sat_m": -sat,
            "reference_external_minus_current_sat_m": reference - sat,
            "reference_zero_minus_translated_sat_m": -(sat - reference)},
        "measured_conclusion": "MAL2 SPP clock is about +102717.4m, then its model clock is zero; NNOR receiver and satellite clock seeds remain at their external source origin",
        "causal_interpretation": "Forcing only the reference receiver to zero is not prediction-equivalent to a common rec/sat clock translation; this is a cold first-epoch boundary, not the second-epoch saved-offset reuse fault",
        "not_proved": ["Whether all required complete covariance blocks were transported", "The final corrected QC outcome", "Physical product quality"]}


def outcome(root):
    plan = json.loads((root / "stage_plan.json").read_text())
    trace_evidence = []
    patterns = ("ZHANG_KF_MEASUREMENT_TRANSACTION", "ZHANG_KF_FINAL_POSTFIT_QC", "PPP AUTHORITATIVE MEASUREMENT TRANSACTION REJECTED", "Bad state detected:")
    for path in sorted((root / "outputs").glob("Network*.TRACE")):
        with path.open(errors="replace") as stream:
            for number, line in enumerate(stream, 1):
                if any(pattern in line for pattern in patterns):
                    trace_evidence.append({"path": str(path), "line": number, "text": line.strip()})
    rolled_back = any("status=ROLLED_BACK" in row["text"] and "FINAL_POSTFIT_QC_NOT_CONVERGED" in row["text"] for row in trace_evidence)
    classification = "TIMEOUT" if plan.get("timed_out") is True else "QC_TRANSACTION_ABORT" if rolled_back else "OTHER_OR_UNPROVEN"
    return {"run_root": str(root), "plan_status": {key: plan.get(key) for key in
        ("state", "binary_sha256", "exit_code", "timed_out", "elapsed_seconds", "expected_epoch_count")},
        "classification": classification, "logs": audit_logs(root), "trace_evidence": trace_evidence,
        "spp_offset_chain": unavailable("Uninstrumented baseline contains no exact source/offset/adjustment lifecycle values; do not substitute a correlated QC clock magnitude"),
        "scientific_acceptance": False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-base", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    args = parser.parse_args()
    report = {"schema_version": 1, "scientific_acceptance": False,
        "galileo_saved_offset_boundary": galileo_chain(args.run_base / "diagnostic_E_2ep"),
        "bds_cold_initialization_boundary": bds_initialization(args.run_base / "diagnostic_C_2ep"),
        "baseline_outcomes": {name: outcome(args.run_base / name) for name in
            ("baseline_G_configured_2ep", "older_verified_v2_G_2ep", "baseline_GEC_configured_2ep")}}
    args.output_root.mkdir(parents=True, exist_ok=True)
    (args.output_root / "diagnostic_numeric_causes.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    lines = ["COLD FLOAT CLOCK DIAGNOSTIC NUMERIC EVIDENCE", "Scientific acceptance: false", "", "GAL: measured saved-offset/source boundary"]
    for chain in report["galileo_saved_offset_boundary"]["receiver_chains"]:
        lines.append(chain["receiver"] + ":")
        for row in chain["observed_spp_results"] + chain["observed_saved_offsets"] + chain["observed_adjustment_results"]:
            event = row["event"]
            lines.append("  " + json.dumps({"line": row["runner_source_line"], **{key: event.get(key) for key in
                ("phase", "epoch", "spp_clock_m", "saved_offset_m", "predicted_clock_m", "unrounded_adjustment_m", "applied_adjustment_m", "skip_reason")},
                "valid_sources": row["spp_valid_source_counts"]}))
        for item in chain["independent_scalar_arithmetic"]:
            lines.append("  arithmetic " + json.dumps(item))
    lines.extend(["", "BDS: measured reference-zero/external-seed boundary", json.dumps(report["bds_cold_initialization_boundary"]["independent_scalar_example"], indent=2)])
    for name, result in report["baseline_outcomes"].items():
        lines.extend(["", name + ": " + result["classification"], json.dumps(result["plan_status"]),
            "  " + json.dumps(result["logs"]["processed_epochs"]),
            "  Exact SPP offset lifecycle: UNAVAILABLE (baseline not instrumented)"])
        lines.extend("  line " + str(row["line"]) + ": " + row["text"] for row in result["trace_evidence"])
    lines.extend(["", "Interpretation limits:",
        "  Source tags and scalar arithmetic are measured boundary evidence, not full source-chart certificates.",
        "  A repair needs complete mean/P prediction-equivalence and real receiver-jump validation.",
        "  First-epoch cold initialization and later saved-offset reuse are different faults.",
        "  AGGO local geometry/rank risk remains separate; QC/noise/stations were not relaxed."])
    (args.output_root / "diagnostic_numeric_causes.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines[:4]))
    print("wrote", args.output_root)


if __name__ == "__main__":
    main()
