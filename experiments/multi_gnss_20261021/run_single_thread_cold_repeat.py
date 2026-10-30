#!/usr/bin/env python3
"""Diagnostic only: two identical cold single-epoch FLOAT executions.

Keeps original physical inputs/static priors and one identical canonical config
and output path. No checkpoint restore, AR, tolerance relaxation or numerical
normalization. Does not alter the failed multi-thread restart experiment.
"""
import argparse
import json
from pathlib import Path
import re
import subprocess

import run_network180_float as network
import run_static_checkpoint_replay as replay
import run_stages as stages

EPOCH = "2024-07-17 00:00:00"


def prepare(template, root):
    references = template["external_review_files"]
    plan, env = replay.prepare(Path(template["binary"]), Path(references[2]["path"]), root,
        Path(references[0]["path"]), Path(references[1]["path"]), Path(template["verifier"]),
        float_threads=1, omp_threads=1)
    overlay = Path(plan["configs"][-2])
    value = json.loads(overlay.read_text())
    value["processing_options"]["epoch_control"]["end_epoch"] = EPOCH
    options = value["processing_options"]["gnss_general"]["zhang_pppar"]
    options["checkpoint_capture_epochs"] = [EPOCH]
    stages.write_json(overlay, value)
    env.update(ZHANG_R51_AR_THREADS="1", OMP_THREAD_LIMIT="1", OMP_MAX_ACTIVE_LEVELS="1",
               OMP_DYNAMIC="FALSE", GOTO_NUM_THREADS="1", BLIS_NUM_THREADS="1")
    # Re-pin only after completing this unused diagnostic configuration.
    audit = json.loads(Path(references[2]["path"]).read_text())
    plan["pinned_input_contract"] = stages.pin_stage_inputs(audit, plan["configs"])
    plan["command"] = stages.independent_stage_command(Path(plan["binary"]), plan["configs"],
                                                       plan["pinned_input_contract"])
    plan.update(expected_epoch_count=1, start_epoch=EPOCH, end_epoch=EPOCH,
        capture_epochs=[EPOCH], diagnostic="COLD_REPEAT_SINGLE_THREAD_NOT_RESTORE_ACCEPTANCE",
        comparison_policy="Two independent cold runs; strict full core/module/product/Q equality",
        environment={key: val for key, val in env.items() if key.startswith("ZHANG_") or
                     key.endswith("_NUM_THREADS") or key.startswith("OMP_")})
    plan["source_identity"].append(network.identity(Path(__file__).resolve()))
    stages.write_json(root / "replay_plan.json", plan)
    return plan, env


def thread_evidence(outputs):
    text = "\n".join(path.read_text(errors="replace") for path in outputs.glob("Network-*.TRACE"))
    rows = {prefix: [line for line in text.splitlines() if prefix in line]
            for prefix in ("ZHANG_FLOAT_THREADS", "ZHANG_AR_THREADS")}
    for prefix, lines in rows.items():
        if not lines or any(not re.search(r"\brequested=1(?:\s|$)", line) for line in lines):
            raise RuntimeError("Actual single-thread scope not proved: " + prefix)
    for line in rows["ZHANG_FLOAT_THREADS"]:
        if any(not re.search(r"\b" + field + r"=1(?:\s|$)", line)
               for field in ("factor_blas", "omp", "blas", "eigen")):
            raise RuntimeError("Actual FLOAT OMP/BLAS/Eigen factor budget is not one")
    return rows


def one_epoch_products(outputs):
    bundles = [path for path in (outputs / "zhang_internal_products.csv.epochs").iterdir()
               if path.is_dir() and ".tmp." not in path.name]
    if len(bundles) != 1:
        raise RuntimeError("Diagnostic must publish exactly one physical epoch")
    doc, hashes = replay.verify_atomic_bundle(bundles[0])
    if doc["epoch"] != EPOCH:
        raise RuntimeError("Wrong diagnostic product epoch")
    return hashes


def execute(plan, env, approval):
    root = Path(plan["run_root"])
    review = replay.execution_gate(plan, approval)
    replay.verify_plan(plan)
    preflight = stages.configuration_preflight(Path(plan["binary"]), plan["configs"],
        plan["pinned_input_contract"], root, env)
    if not preflight["passed"]:
        raise RuntimeError("Actual production parser rejected single-thread diagnostic")
    plan.update(state="DIAGNOSTIC_RUNNING", configuration_preflight=preflight,
                controlled_execution_approval=review)
    stages.write_json(root / "replay_plan.json", plan)
    replay.run_phase(plan, env, "A_cold_one_thread")
    first = replay.archive_owned_outputs(root, "continuous")
    first_threads = thread_evidence(first / "outputs")
    for name in ("outputs", "checkpoints"):
        (root / "work" / name).mkdir()
    replay.run_phase(plan, env, "B_cold_one_thread")
    second = root / "work"
    second_threads = thread_evidence(second / "outputs")
    left = replay.checkpoint_at(first / "checkpoints", EPOCH)
    right = replay.checkpoint_at(second / "checkpoints", EPOCH)
    if left is None or right is None:
        raise RuntimeError("Actual diagnostic checkpoint publication missing")
    comparison_path = root / "cold_repeat_full_checkpoint_comparison.json"
    with (root / "checkpoint_comparison.log").open("w") as output:
        compared = subprocess.run([plan["verifier"], str(left / "checkpoint.bundle"),
            str(right / "checkpoint.bundle"), plan["binary_sha256"], replay.MODE, str(comparison_path)],
            stdout=output, stderr=subprocess.STDOUT)
    products_a, products_b = one_epoch_products(first / "outputs"), one_epoch_products(second / "outputs")
    report = {"diagnostic_only": True, "scientific_acceptance": False,
        "checkpoint_full_wire_passed": compared.returncode == 0,
        "products_and_full_q_exact": products_a == products_b,
        "first_product_sha256": products_a, "second_product_sha256": products_b,
        "actual_thread_scope_evidence": {"first": first_threads, "second": second_threads},
        "comparison_report": str(comparison_path), "numeric_normalization": False}
    stages.write_json(root / "cold_repeat_result.json", report)
    plan.update(state="DIAGNOSTIC_COLD_REPEAT_IDENTICAL" if compared.returncode == 0 and
                products_a == products_b else "DIAGNOSTIC_COLD_REPEAT_DIFFERENT",
                finished_utc=replay.stamp(), diagnostic_result=report)
    replay.verify_plan(plan)
    stages.write_json(root / "replay_plan.json", plan)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--template-plan", type=Path, required=True)
    parser.add_argument("--run-root", type=Path, required=True)
    parser.add_argument("--approval", type=Path, required=True)
    args = parser.parse_args()
    template = json.loads(args.template_plan.read_text())
    plan, env = prepare(template, args.run_root.resolve())
    try:
        result = execute(plan, env, args.approval)
    except Exception as error:
        plan.update(state="DIAGNOSTIC_FAILED", failure_reason=str(error), finished_utc=replay.stamp())
        stages.write_json(Path(plan["run_root"]) / "replay_plan.json", plan)
        raise
    print(json.dumps(result, indent=2), flush=True)


if __name__ == "__main__":
    main()
