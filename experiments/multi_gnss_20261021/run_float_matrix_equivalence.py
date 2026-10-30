#!/usr/bin/env python3
"""Prepare, or explicitly approve, two cold 12-site G/E/C FLOAT epochs.

Both arms use the same producer, paths, physical models, OMP=4 and BLAS=1.
Only ZHANG_FLOAT_MEMORY_MODE differs. No legacy rank SVD is executed. This is
an exact-equivalence contract, not AR, accuracy or 180-site throughput evidence.
--self-test runs guard tests only and never starts PEA.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import platform
import re
import subprocess
import sys

import run_float_rank_diagnostic_equivalence as rank
import run_single_thread_cold_repeat as cold
import run_static_checkpoint_replay as replay
import run_stages as stages
import run_network180_float as network
import audit_network180_float_completion as completion

STATE = "FLOAT_MATRIX_MEMORY_EQUIVALENCE_PREPARED_NOT_RUN"
VARIABLE = "ZHANG_FLOAT_MEMORY_MODE"
ARMS = (("serial", "continuous"), ("parallel", None))
POLICY = {"ZHANG_R51_FLOAT_THREADS": "4", "OMP_NUM_THREADS": "4",
          "OMP_THREAD_LIMIT": "4", "ZHANG_FLOAT_MATRIX_THREADS": "1",
          "ZHANG_FLOAT_PURE_OBS_RANK_DIAGNOSTIC": "off", "ZHANG_FLOAT_PROFILE": "1"}
CORE_FIELDS = ("x P dx prefitRatios postfitRatios time kfIndexMap stateTransitionMap "
    "gaussMarkovTauMap gaussMarkovMuMap procNoiseMap initNoiseMap sigmaMaxMap outageLimitMap "
    "exponentialNoiseMap pseudoStateMap pseudoParentMap errorCountMap filterChunkMap metaDataMap "
    "lsqRequired sigmaPass chiQCPass chi2 dof chi2PerDof qc id rtsBasename outputResiduals "
    "outputMongoMeasurements statisticsMap statisticsMapSum").split()


def timeout_value(value):
    if not math.isfinite(value) or not 0 < value <= 600:
        raise ValueError("A finite 0 < timeout <= 600 seconds per cold arm is required")
    return value


def arm_environment(base, mode):
    if mode not in ("serial", "parallel"):
        raise ValueError("Only serial and parallel memory arms are permitted")
    if any(base.get(key) != value for key, value in POLICY.items()):
        raise RuntimeError("Both arms must keep the reviewed outer4/matrix1/rank-off environment")
    result = dict(base)
    result[VARIABLE] = mode
    return result


def pinned_snapshot(binary, archive, manifest_path, input_manifest):
    """Bind actual producer, complete source snapshot, helpers, runtime and models."""
    manifest = json.loads(manifest_path.read_text())
    if (manifest.get("binary_sha256") != stages.file_sha(binary) or
            manifest.get("archive_sha256") != stages.file_sha(archive)):
        raise RuntimeError("Producer/source archive SHA binding failed")
    repo = Path(__file__).resolve().parents[2]
    identities = [network.identity(path) for path in (archive, manifest_path, input_manifest)]
    source_names = set()
    for row in manifest["source_files"]:
        path = (repo / row["path"]).resolve(strict=True)
        if repo not in path.parents:
            raise RuntimeError("Source manifest path escapes the reviewed repository")
        actual = network.identity(path)
        if actual["sha256"] != row["sha256"]:
            raise RuntimeError("Live source differs from frozen producer source: " + row["path"])
        source_names.add(row["path"])
        identities.append(actual)
    helpers = (Path(__file__), Path(rank.__file__), Path(cold.__file__), Path(replay.__file__),
               Path(stages.__file__), Path(network.__file__), Path(completion.__file__))
    if any(str(path.resolve().relative_to(repo)).replace("\\", "/") not in source_names for path in helpers):
        raise RuntimeError("All executed helpers must be in the frozen producer source archive")
    models = json.loads(input_manifest.read_text())["files"]
    if not models:
        raise RuntimeError("A nonempty physical/model input manifest is required")
    for row in models:
        actual = network.identity(Path(row["path"]))
        if (actual["bytes"], actual["sha256"]) != (row["bytes"], row["sha256"]):
            raise RuntimeError("Physical/model input changed: " + row["path"])
        identities.append(actual)
    runtime = [network.identity(path) for path in sorted((binary.parent / "lib").rglob("*")) if path.is_file()]
    if not runtime:
        raise RuntimeError("Private producer runtime libraries must be pinned")
    identities.extend(runtime)
    return identities


def prepare(binary, audit, root, coordinate_overlay, coordinate_review, verifier,
            source_archive, source_manifest, input_manifest, timeout=600):
    timeout_value(timeout)
    binary, verifier = binary.resolve(strict=True), verifier.resolve(strict=True)
    snapshot = pinned_snapshot(binary, source_archive, source_manifest, input_manifest)
    template = {"binary": str(binary), "verifier": str(verifier), "external_review_files": [
        {"path": str(coordinate_overlay.resolve(strict=True))},
        {"path": str(coordinate_review.resolve(strict=True))}, {"path": str(audit.resolve(strict=True))}]}
    plan, env = cold.prepare(template, root)
    baseline = dict(env)
    env.pop(VARIABLE, None)
    env.update(POLICY)
    changed = {key: {"old": baseline.get(key), "new": env.get(key)} for key in set(baseline) | set(env)
               if baseline.get(key) != env.get(key)}
    if set(changed) - set(POLICY) - {VARIABLE}:
        raise RuntimeError("Unreviewed environment difference from the original static cold model")
    plan["source_identity"].extend([*snapshot, network.identity(verifier)])
    plan.update(state=STATE, diagnostic="FLOAT_MEMORY_SERIAL_PARALLEL_EXACT_EQUIVALENCE",
        host=platform.node(), timeout_seconds=timeout, timeout_scope="PER_INDEPENDENT_COLD_ARM",
        verifier_sha256=stages.file_sha(verifier), arms=[mode for mode, _ in ARMS],
        environment={key: value for key, value in env.items() if key.startswith(("ZHANG_", "OMP_"))
                     or key.endswith("_NUM_THREADS")}, reviewed_environment_changes=changed,
        comparison_policy="Strict complete x/P/order/F/Q/priors/QC/module bytes and product/full-Q/certificate/COMMIT; no numerical normalization or tolerances",
        matrix_thread_policy="OMP/outer=4; every matrix BLAS=1 in both arms; only new memory-helper scheduling differs",
        memory_policy_scope="MEMORY_MODE governs new block/transpose/mirror helpers; existing copyDense is unchanged in both arms",
        scientific_acceptance=False)
    stages.write_json(root / "replay_plan.json", plan)
    stages.write_json(root / "matrix_equivalence_prepared_plan.json", plan)
    replay.verify_plan(plan)
    return plan, env


def thread_evidence(outputs, mode):
    text = "\n".join(path.read_text(errors="replace") for path in outputs.glob("Network-*.TRACE"))
    rows = [line for line in text.splitlines() if "ZHANG_FLOAT_THREADS " in line]
    expected = {"requested": 4, "omp": 4, "blas": 1, "factor_blas": 1, "eigen": 1}
    if not rows or any(not re.search(r"\b" + key + "=" + str(value) + r"(?:\s|$)", line)
                           for line in rows for key, value in expected.items()):
        raise RuntimeError("Actual outer4/matrix1/Eigen1 FLOAT thread budget was not demonstrated")
    kernels = [line for line in text.splitlines() if "ZHANG_FLOAT_KERNEL phase=GAIN_TRANSPOSE " in line]
    dimensions = [re.search(r"\bstates=(\d+) measurements=(\d+)\b", line) for line in kernels]
    if not kernels or any(match is None or int(match[1]) * int(match[2]) < 262144 for match in dimensions):
        raise RuntimeError("Actual gain transpose did not exercise the large-memory-helper size threshold")
    if "ZHANG_FLOAT_COVARIANCE_MODE symmetric=0 " not in text:
        raise RuntimeError("This memory comparison must observe the ordinary gain-transpose path, not infer it from a timer")
    if "boundary=B16 status=AUTHORITATIVE_COMMIT" not in text or "status=ROLLED_BACK" in text:
        raise RuntimeError("Authoritative FLOAT commit is missing or a rollback occurred")
    return {"float_thread_lines": rows, "large_memory_kernel_lines": kernels, "memory_mode": mode,
            "parallel_evidence_scope": "Reviewed helper policy plus actual OMP budget and matrix dimensions; not OS per-thread sampling"}


def strict_comparison(report, receivers):
    rank.verify_checkpoint_comparison(report, receivers)
    checks = report.get("checks", {})
    if any(checks.get("core." + key) is not True for key in CORE_FIELDS) or any(value is not True for value in checks.values()):
        raise RuntimeError("Complete serialized KF fields/QC/maps must all compare exactly")
    for key in ("x", "P", "dx", "prefitRatios", "postfitRatios"):
        metric = report.get("numeric_" + key, {})
        if (metric.get("all_finite") is not True or metric.get("exact_equal") is not True or
                metric.get("maximum_absolute_difference_native") != 0):
            raise RuntimeError("Nonfinite/nonexact checkpoint numeric field: " + key)


def execute(plan, env, approval):
    if sys.platform != "linux":
        raise RuntimeError("Numerical comparison requires native Linux")
    root = replay.owner_guard(Path(plan["run_root"]))
    prepared = root / "matrix_equivalence_prepared_plan.json"
    marker = root / "matrix_execution_started.json"
    if (plan.get("state") != STATE or marker.exists() or plan != json.loads(prepared.read_text()) or
            json.loads((root / "replay_plan.json").read_text())["state"] != STATE):
        raise RuntimeError("Only the exact unused prepared plan may run")
    timeout_value(plan["timeout_seconds"])
    for mode, _ in ARMS:
        arm_environment(env, mode)
    review = replay.execution_gate(plan, approval)
    if (review.get("approved_for_float_matrix_memory_equivalence") is not True or
            review.get("prepared_plan_sha256") != stages.file_sha(prepared) or
            review.get("host") != platform.node() or plan.get("host") != platform.node()):
        raise RuntimeError("Exact host/plan/root approval for memory equivalence is required")
    replay.verify_plan(plan)
    preflight = stages.configuration_preflight(Path(plan["binary"]), plan["configs"],
        plan["pinned_input_contract"], root, env)
    if not preflight["passed"]:
        raise RuntimeError("Actual production parser rejected the unchanged 12-site model")
    stages.write_json(marker, {"started_utc": replay.stamp(), "approval": review,
                              "binary_sha256": plan["binary_sha256"]})
    plan.update(state="FLOAT_MATRIX_MEMORY_EQUIVALENCE_RUNNING", configuration_preflight=preflight)
    stages.write_json(root / "replay_plan.json", plan)
    arms = {}
    for mode, archive in ARMS:
        replay.verify_plan(plan)
        run = replay.run_phase(plan, arm_environment(env, mode), "memory_" + mode)
        log = Path(run["log"]).read_text(errors="replace")
        if "phase=FLOAT_PURE_OBSERVATION_COPY " in log or any(network.authoritative_failure(line) for line in log.splitlines()):
            raise RuntimeError("Rank-off policy or authoritative FLOAT transaction failed")
        diagnostic = rank.diagnostic_evidence(run["log"], "off")
        work = replay.archive_owned_outputs(root, archive) if archive else root / "work"
        evidence = thread_evidence(work / "outputs", mode)
        products = cold.one_epoch_products(work / "outputs")
        bundles = [path for path in (work / "outputs/zhang_internal_products.csv.epochs").iterdir()
                   if path.is_dir() and ".tmp." not in path.name]
        product_audit = completion.audit_epoch(bundles[0])
        if product_audit["epoch"] != cold.EPOCH:
            raise RuntimeError("Finite complete G/E/C product/Q audit belongs to another epoch")
        checkpoint = replay.checkpoint_at(work / "checkpoints", cold.EPOCH)
        if checkpoint is None:
            raise RuntimeError("Actual single cold-epoch checkpoint is missing")
        arms[mode] = {"run": run, "diagnostic_evidence": diagnostic, "thread_evidence": evidence,
                      "checkpoint": str(checkpoint), "products_sha256": products,
                      "finite_complete_product_audit": product_audit}
        if archive:
            for name in ("outputs", "checkpoints"):
                (root / "work" / name).mkdir()
    path = root / "serial_versus_parallel_full_checkpoint.json"
    with (root / "serial_versus_parallel_reader.log").open("w") as output:
        checked = subprocess.run([plan["verifier"],
            str(Path(arms["serial"]["checkpoint"]) / "checkpoint.bundle"),
            str(Path(arms["parallel"]["checkpoint"]) / "checkpoint.bundle"),
            plan["binary_sha256"], replay.MODE, str(path)], stdout=output, stderr=subprocess.STDOUT, timeout=120)
    if checked.returncode:
        raise RuntimeError("Independent complete checkpoint comparison failed")
    strict_comparison(json.loads(path.read_text()), {row["station"] for row in plan["pinned_input_contract"]["observations"]})
    if arms["serial"]["products_sha256"] != arms["parallel"]["products_sha256"]:
        raise RuntimeError("Complete products/full-Q/certificate/COMMIT byte identities differ")
    replay.verify_plan(plan)
    result = {"passed": True, "arms": arms, "comparison_report": str(path),
        "checkpoint_full_wire_exact": True, "products_and_full_q_exact": True,
        "actual_epoch": cold.EPOCH, "station_count": 12, "systems": ["GPS", "GAL", "BDS"],
        "ratio": 2, "all_nis": False, "ar_executed": False, "rank_svd_executed": False,
        "numeric_normalization": False, "scientific_acceptance": False,
        "scope": "Two independent cold 12-site FLOAT epochs; OMP4 and BLAS1; not AR, 180-site performance, full-hour or accuracy acceptance"}
    stages.write_json(root / "matrix_equivalence_result.json", result)
    plan.update(state="FLOAT_MATRIX_MEMORY_EXACT_EQUIVALENCE_PASSED", finished_utc=replay.stamp())
    stages.write_json(root / "replay_plan.json", plan)
    return result


def guard_tests():
    import tempfile
    import unittest
    from unittest.mock import patch

    class Guards(unittest.TestCase):
        def test_timeout(self):
            self.assertEqual(timeout_value(600), 600)
            for value in (0, -1, 601, float("nan"), float("inf")):
                with self.assertRaises(ValueError): timeout_value(value)

        def test_arm_only_changes_memory_mode(self):
            base = dict(POLICY, SENTINEL="keep")
            first, second = arm_environment(base, "serial"), arm_environment(base, "parallel")
            self.assertEqual({key for key in first if first[key] != second[key]}, {VARIABLE})
            self.assertNotIn(VARIABLE, base)
            with self.assertRaises(ValueError): arm_environment(base, "full")
            base["ZHANG_FLOAT_MATRIX_THREADS"] = "4"
            with self.assertRaises(RuntimeError): arm_environment(base, "serial")

        def test_thread_and_large_kernel_evidence(self):
            with tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                path = root / "Network-case.TRACE"
                good = ("ZHANG_FLOAT_THREADS requested=4 omp=4 blas=1 factor_blas=1 eigen=1\n"
                        "ZHANG_FLOAT_KERNEL phase=GAIN_TRANSPOSE states=966 measurements=1079 elapsed_ms=1\n"
                        "ZHANG_FLOAT_COVARIANCE_MODE symmetric=0 reason=LEGACY_CANCELLATION_RISK\n"
                        "ZHANG_KF_MEASUREMENT_TRANSACTION boundary=B16 status=AUTHORITATIVE_COMMIT\n")
                path.write_text(good)
                self.assertEqual(thread_evidence(root, "parallel")["memory_mode"], "parallel")
                for changed in (good.replace("blas=1", "blas=4"), good.replace("966", "2"),
                                good.replace("symmetric=0", "symmetric=1"),
                                good.replace("AUTHORITATIVE_COMMIT", "ROLLED_BACK")):
                    path.write_text(changed)
                    with self.assertRaises(RuntimeError): thread_evidence(root, "parallel")

        def test_duplicate_plan_never_launches(self):
            with tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary).resolve()
                (root / "OWNERSHIP.json").write_text(json.dumps({"schema": replay.OWNERSHIP, "root": str(root)}))
                (root / "matrix_execution_started.json").write_text("{}")
                with patch.object(sys, "platform", "linux"), patch.object(replay, "execution_gate") as gate, \
                     patch.object(subprocess, "Popen") as launch:
                    with self.assertRaises(RuntimeError): execute({"state": STATE, "run_root": str(root)}, {}, root / "approval")
                    gate.assert_not_called()
                    launch.assert_not_called()

        def test_missing_exact_plan_approval_never_launches(self):
            with tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary).resolve()
                (root / "OWNERSHIP.json").write_text(json.dumps({"schema": replay.OWNERSHIP, "root": str(root)}))
                plan = {"state": STATE, "run_root": str(root), "timeout_seconds": 600, "host": platform.node()}
                for name in ("matrix_equivalence_prepared_plan.json", "replay_plan.json"):
                    (root / name).write_text(json.dumps(plan))
                with patch.object(sys, "platform", "linux"), \
                     patch.object(replay, "execution_gate", return_value={"approved_for_float_matrix_memory_equivalence": True}), \
                     patch.object(replay, "verify_plan") as verify, patch.object(subprocess, "Popen") as launch:
                    with self.assertRaises(RuntimeError): execute(plan, POLICY, root / "approval")
                    verify.assert_not_called()
                    launch.assert_not_called()

        def test_wrong_producer_archive_binding_fails_before_source_inventory(self):
            with tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                binary, archive, manifest = root / "pea", root / "source.zip", root / "manifest.json"
                binary.write_bytes(b"producer")
                archive.write_bytes(b"source")
                manifest.write_text(json.dumps({"binary_sha256": "wrong", "archive_sha256": stages.file_sha(archive)}))
                with self.assertRaisesRegex(RuntimeError, "Producer/source archive"):
                    pinned_snapshot(binary, archive, manifest, root / "never_read.json")

        def test_complete_scientific_fields_cannot_be_omitted(self):
            report = {"checks": {"core." + key: True for key in CORE_FIELDS}, **{
                "numeric_" + key: {"all_finite": True, "exact_equal": True, "maximum_absolute_difference_native": 0}
                for key in ("x", "P", "dx", "prefitRatios", "postfitRatios")}}
            with patch.object(rank, "verify_checkpoint_comparison", return_value=True):
                strict_comparison(report, set())
                report["checks"].pop("core.procNoiseMap")
                with self.assertRaises(RuntimeError): strict_comparison(report, set())
                report["checks"]["core.procNoiseMap"] = True
                report["numeric_P"]["maximum_absolute_difference_native"] = 1e-99
                with self.assertRaises(RuntimeError): strict_comparison(report, set())

    return unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(Guards)).wasSuccessful()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("binary", "audit", "run-root", "coordinate-overlay", "coordinate-review", "verifier",
                 "source-archive", "source-manifest", "input-manifest", "prepared-plan", "build-approval"):
        parser.add_argument("--" + name, type=Path)
    parser.add_argument("--timeout-seconds", type=float, default=600)
    parser.add_argument("--execute", action="store_true")
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        if args.execute: parser.error("--self-test never permits numerical execution")
        raise SystemExit(0 if guard_tests() else 1)
    if args.execute and (args.prepared_plan is None or args.build_approval is None):
        parser.error("Execution requires an exact existing prepared plan and root approval")
    if args.prepared_plan:
        plan = json.loads(args.prepared_plan.read_text())
        env = stages.environment("D")
        env.update(plan["environment"])
    else:
        names = ("binary", "audit", "run_root", "coordinate_overlay", "coordinate_review", "verifier",
                 "source_archive", "source_manifest", "input_manifest")
        if any(getattr(args, key) is None for key in names):
            parser.error("Preparation requires explicit producer/source/input/config/review/verifier paths")
        plan, env = prepare(*(getattr(args, key) for key in names), timeout=args.timeout_seconds)
    if args.execute:
        try:
            result = execute(plan, env, args.build_approval)
        except Exception as error:
            plan.update(state="FLOAT_MATRIX_MEMORY_EQUIVALENCE_FAILED", failure_reason=str(error),
                        finished_utc=replay.stamp(), scientific_acceptance=False)
            stages.write_json(Path(plan["run_root"]) / "replay_plan.json", plan)
            raise
        print(json.dumps(result, indent=2), flush=True)
    else:
        print(json.dumps(plan, indent=2), flush=True)


if __name__ == "__main__":
    main()
