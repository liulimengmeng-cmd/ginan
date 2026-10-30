#!/usr/bin/env python3
"""Read-only intermediate candidate replay audit; never final acceptance."""
from __future__ import annotations

import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
import json
from pathlib import Path

from analyze_clock_lifecycle import compact, outcome, provenance, read_events, spp_result_evidence
from audit_stage_results import StageAuditor, finite, sha256, unavailable


def raw_phase_changes(root):
    previous, changes = {}, []
    base = root / "outputs/zhang_internal_products.csv.epochs"
    for directory in sorted(base.iterdir()) if base.is_dir() else []:
        if not (directory / "COMMIT.json").is_file() or directory.name.endswith(".pending"):
            continue
        with (directory / "products.csv").open(newline="") as stream:
            for row in csv.DictReader(stream):
                key = (row["solution"], row["satellite"], row["observable"])
                value = finite(row.get("phase_m"))
                if value is None:
                    continue
                if key in previous:
                    changes.append({"solution": key[0], "satellite": key[1], "observable": key[2],
                        "phase_difference_m": value - previous[key], "absolute_difference_m": abs(value - previous[key])})
                previous[key] = value
    return {"paired_endpoint_count": len(changes),
        "maximum_absolute_difference_m": max((r["absolute_difference_m"] for r in changes), default=None),
        "largest_ten": sorted(changes, key=lambda row: row["absolute_difference_m"], reverse=True)[:10],
        "interpretation": "Two cold-start FLOAT epochs only; raw changes are pending convergence/quality evaluation, not certified jumps, physical slips or continuity failures"}


def offset_chain(events, receiver="NNOR"):
    selected = [row for row in events if row["event"].get("receiver") == receiver and
        row["event"]["phase"] in ("SPP_RESULT", "SPP_OFFSET_SAVE_AFTER", "REC_CLOCK_ADJUST_AFTER")]
    formula = []
    for row in selected:
        event = row["event"]
        spp = event.get("spp_system_clock_m", event.get("spp_clock_m"))
        offset = event.get("saved_system_offset_m", event.get("saved_offset_m"))
        prediction = event.get("predicted_clock_m")
        if event["phase"] == "REC_CLOCK_ADJUST_AFTER" and None not in (spp, offset, prediction):
            formula.append({"runner_source_line": row["source_line"], "epoch": event["epoch"],
                "system": event.get("system", event.get("spp_system")), "external_spp_clock_m": spp,
                "saved_external_to_kf_offset_m": offset, "predicted_kf_clock_m": prediction,
                "independently_recomputed_adjustment_m": spp + offset - prediction,
                "observed_applied_adjustment_m": event.get("applied_adjustment_m"),
                "formula_residual_m": spp + offset - prediction - event["unrounded_adjustment_m"]})
    return {"receiver": receiver, "events": [compact(row) for row in selected], "formula_checks": formula}


def compare_cold_seeds(old_root, new_root, reference_receiver="MAL2"):
    def seeds(root):
        result = {}
        for row in read_events(root):
            event = row["event"]
            if event["phase"] in ("COLD_RECEIVER_CLOCK_SEED", "COLD_SATELLITE_CLOCK_SEED") and event.get("init_clock_m") is not None:
                key = tuple(event.get(k) for k in ("phase", "receiver", "satellite", "system"))
                result.setdefault(key, row)
        return result
    old, new = seeds(old_root), seeds(new_root)
    reference = next(row for row in read_events(new_root) if row["event"]["phase"] == "COLD_REFERENCE_CLOCK_ZERO" and row["event"].get("receiver") == reference_receiver)
    reference_spp = reference["event"]["spp_clock_m"]
    differences = []
    for key in sorted(old.keys() & new.keys()):
        before, after = old[key], new[key]
        difference = after["event"]["init_clock_m"] - before["event"]["init_clock_m"]
        differences.append({"phase": key[0], "receiver": key[1], "satellite": key[2], "system": key[3],
            "old_runner_line": before["source_line"], "new_runner_line": after["source_line"],
            "old_seed_m": before["event"]["init_clock_m"], "new_seed_m": after["event"]["init_clock_m"],
            "new_minus_old_m": difference, "difference_plus_reference_spp_m": difference + reference_spp})
    return {"reference_zero_event": compact(reference), "reference_spp_clock_m": reference_spp,
        "matched_seed_count": len(differences),
        "maximum_abs_translation_residual_m": max((abs(row["difference_plus_reference_spp_m"]) for row in differences), default=None),
        "seed_comparisons": differences,
        "scope": "Actual recorded cold scalar seeds; full state/covariance transport correctness needs independent production contracts, not just common scalar shifts"}


def audit(root):
    plan = json.loads((root / "stage_plan.json").read_text())
    if plan.get("state") == "RUNNING" or "exit_code" not in plan:
        return {"run_root": str(root), "state": plan.get("state"), "completion": "NOT_ESTABLISHED",
            "scientific_acceptance": False, "reason": "Snapshot not complete; mutable products not audited as final"}
    events = read_events(root)
    evidence = StageAuditor(root, 30, .1, 1000, 2000).run()
    execution = outcome(root)
    # outcome() also serves uninstrumented baseline runs. Do not repeat its
    # baseline-only unavailability disclaimer for an instrumented candidate.
    execution.pop("spp_offset_chain", None)
    return {"run_root": str(root), "scope": "Intermediate candidate replay, not final candidate03 source/build acceptance",
        "scientific_acceptance": False, "plan_sha256": sha256(root / "stage_plan.json"),
        "binary_sha256": plan["binary_sha256"], "effective_input_preflight": plan["effective_input_preflight"],
        "pinning_after_run": {key: plan.get(key) for key in ("binary_unchanged_after_run", "pinned_inputs_unchanged_after_run")},
        "provenance": provenance(root), "execution_and_qc": execution,
        "independently_recomputed_product_evidence": evidence,
        "explicit_spp_results": spp_result_evidence(events), "nnor_clock_chain": offset_chain(events),
        "raw_phase_changes": raw_phase_changes(root),
        "spp_valid_sources_all_result_events": dict(Counter(obs["source"] for row in events if row["event"]["phase"] == "SPP_RESULT"
            for obs in row["event"].get("observations", []) if obs.get("spp_valid") is True))}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-base", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    parser.add_argument("--runs", nargs="+", default=["candidate01_E_2ep", "candidate01_C_2ep", "candidate02_GEC_2ep"])
    args = parser.parse_args()
    results = {name: audit(args.run_base / name) for name in args.runs}
    if "candidate01_C_2ep" in results and results["candidate01_C_2ep"].get("completion") != "NOT_ESTABLISHED":
        results["candidate01_C_2ep"]["cold_seed_translation_against_instrumented_old"] = compare_cold_seeds(
            args.run_base / "diagnostic_C_2ep", args.run_base / "candidate01_C_2ep")
    report = {"schema_version": 1, "generated_utc": datetime.now(timezone.utc).isoformat(),
        "scientific_acceptance": False, "final_source_build_acceptance": False, "runs": results}
    args.output_root.mkdir(parents=True, exist_ok=True)
    (args.output_root / "candidate_replay_audit.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    lines = ["INTERMEDIATE CANDIDATE REAL REPLAY AUDIT", "Not final source/build or scientific acceptance", ""]
    for name, result in results.items():
        lines.append(name)
        if result.get("completion") == "NOT_ESTABLISHED":
            lines.append("  " + result["completion"])
            continue
        product = result["independently_recomputed_product_evidence"]
        lines.extend(["  binary " + result["binary_sha256"], "  processed " + json.dumps(product["logs"]["processed_epochs"]),
            "  timing " + json.dumps(product["logs"]["epoch_timing"]),
            "  hash-verified atomic bundles " + str(product["atomic"]["hash_verified_bundles"]),
            "  row totals " + json.dumps(product["products"]["totals_by_solution_system"]),
            "  SPP valid sources " + json.dumps(result["spp_valid_sources_all_result_events"]),
            "  NNOR actual formula " + json.dumps(result["nnor_clock_chain"]["formula_checks"]),
            "  evidence issues " + json.dumps(product["evidence_issues"]),
            "  raw cold phase max " + str(result["raw_phase_changes"]["maximum_absolute_difference_m"]) + "m; quality/continuity evaluation pending"])
        if "cold_seed_translation_against_instrumented_old" in result:
            translation = result["cold_seed_translation_against_instrumented_old"]
            lines.append("  cold common translation: reference=" + str(translation["reference_spp_clock_m"]) +
                "m, matched=" + str(translation["matched_seed_count"]) + ", max residual=" + str(translation["maximum_abs_translation_residual_m"]) + "m")
        lines.append("")
    (args.output_root / "candidate_replay_audit.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
