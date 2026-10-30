"""Bounded server smoke guard contracts; no PEA or scientific run is started."""
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import run_server_first_epoch_validation as validation
import audit_network180_float_completion as completion
import run_network180_float as network
import run_stages as stages


class FirstEpochValidationContracts(unittest.TestCase):
    def test_one_epoch_limits_and_narrow_policy(self):
        self.assertEqual(validation.EPOCH, "2024-07-17 00:00:00")
        self.assertEqual(validation.server.SERVER_LIMITS["maximum_address_space_bytes"], 64 * 1024 ** 3)
        self.assertEqual(validation.ALLOWED_POLICY_ENV,
                         {"ZHANG_FLOAT_PROFILE": "1", "ZHANG_FLOAT_PURE_OBS_RANK_DIAGNOSTIC": "off"})

    def test_explicit_matrix_budget_and_memory_policy_are_strict(self):
        self.assertEqual(validation.matrix_environment(14, "parallel"),
            {"ZHANG_FLOAT_MATRIX_THREADS": "14", "ZHANG_FLOAT_MEMORY_MODE": "parallel"})
        self.assertEqual(validation.matrix_environment(None, "serial"), {"ZHANG_FLOAT_MEMORY_MODE": "serial"})
        for count in (0, -1, 29, True, "14"):
            with self.assertRaises(ValueError):
                validation.matrix_environment(count, "parallel")
        with self.assertRaises(ValueError):
            validation.matrix_environment(14, "unknown")

    def test_unapproved_or_used_plan_never_starts_numerics(self):
        for state in ("PROFILE_RUNNING", validation.PREPARED_STATE):
            with tempfile.TemporaryDirectory() as temp:
                root = Path(temp)
                plan = {"run_root": str(root)}
                (root / "prepared_server_validation_plan.json").write_text(json.dumps(plan))
                (root / "stage_plan.json").write_text(json.dumps({"state": state}))
                approval = root / "approval.json"
                approval.write_text("{}")
                with patch.object(validation.server.subprocess, "Popen") as launch:
                    with self.assertRaises(RuntimeError):
                        validation.execute(plan, approval, stages, None, None)
                    launch.assert_not_called()

    def test_reference_inputs_reject_count_difference(self):
        reference = {"pinned_input_contract": {"observations": [{"path": "/a"}], "navigation": []}}
        plan = {"pinned_input_contract": {"observations": [], "navigation": []}}
        with self.assertRaisesRegex(RuntimeError, "counts differ"):
            validation.verify_reference_inputs(plan, reference, None, None)

    def test_reference_inputs_reject_changed_same_name(self):
        old = {"path": "/old/a", "bytes": 10, "sha256": "old"}
        new = {"path": "/new/a", "bytes": 10, "sha256": "new"}
        reference = {"pinned_input_contract": {"observations": [old], "navigation": []}}
        plan = {"pinned_input_contract": {"observations": [new], "navigation": []}}
        with self.assertRaisesRegex(RuntimeError, "Physical input differs"):
            validation.verify_reference_inputs(plan, reference, None, None)

    def test_pinned_model_inputs_are_rechecked_at_execution_not_only_preparation(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            binary, audit, model = (root / name for name in ("pea", "audit.json", "orbit.sp3"))
            for path in (binary, audit, model):
                path.write_text("original")
            plan = {"binary": str(binary), "binary_sha256": stages.file_sha(binary),
                "pinned_input_contract": {"observations": [], "navigation": [], "configurations": []},
                "reviewed_coordinate_files": [], "input_audit": network.identity(audit),
                "source_identity": [network.identity(model)]}
            with patch.object(network, "exact_input_cli", return_value=[]):
                network.verify_pins(plan)
                model.write_text("changed precise orbit")
                with self.assertRaisesRegex(RuntimeError, "source changed"):
                    network.verify_pins(plan)

    def fixture(self, directory, text):
        root = Path(directory)
        (root / "runner.log").write_text(text)
        epoch = root / "outputs/zhang_internal_products.csv.epochs/e000"
        epoch.mkdir(parents=True)
        (epoch / "COMMIT.json").write_text("{}")
        checkpoint = root / "checkpoints/e000"
        checkpoint.mkdir(parents=True)
        (checkpoint / "manifest.json").write_text("{}")
        return {"run_root": str(root), "exit_code": 0, "guard_failure_reason": None,
                "resource_stop_reasons": [], "owned_processes_still_live": [], "child_reaped": True}

    @property
    def good_text(self):
        return ("ZHANG_PURE_OBS_RANK_STATUS mode=off status=DISABLED rank=NOT_EVALUATED\n"
                "ZHANG_FLOAT_PHASE_END phase=FLOAT_TOTAL epoch=2024-07-17_00:00:00\n")

    def test_no_exit_success_without_finite_atomic_epoch(self):
        with tempfile.TemporaryDirectory() as temp:
            plan = self.fixture(temp, self.good_text)
            with patch.object(completion, "audit_epoch", side_effect=RuntimeError("Nonfinite product covariance")):
                report = validation.output_audit(plan, stages)
            self.assertFalse(report["passed"])
            self.assertIn("Nonfinite", report["failure_reason"])

    def test_copy_or_svd_execution_fails_off_policy_even_with_exit_zero(self):
        for phase in ("FLOAT_RANK_JACOBI_SVD", "FLOAT_PURE_OBSERVATION_COPY"):
            with tempfile.TemporaryDirectory() as temp:
                plan = self.fixture(temp, self.good_text + "ZHANG_FLOAT_PHASE_BEGIN phase=" + phase + " epoch=x\n")
                report = validation.output_audit(plan, stages)
                self.assertFalse(report["passed"])
                self.assertIn("copied or decomposed", report["failure_reason"])

    def test_missing_explicit_off_status_fails(self):
        with tempfile.TemporaryDirectory() as temp:
            plan = self.fixture(temp, "ZHANG_FLOAT_PHASE_END phase=FLOAT_TOTAL epoch=x\n")
            self.assertFalse(validation.output_audit(plan, stages)["passed"])

    def test_actual_blas_budget_and_parallel_memory_must_match_pinned_plan(self):
        for blas, memory in ((14, "parallel"), (28, "parallel"), (14, "serial")):
            with tempfile.TemporaryDirectory() as temp:
                plan = self.fixture(temp, self.good_text)
                plan.update(matrix_thread_budget=14, memory_mode="parallel")
                (Path(temp) / "outputs/Network-case.TRACE").write_text(
                    "ZHANG_FLOAT_THREADS policy=phase_scoped_all_cpu_v2 requested=28 available=28 "
                    "factor_blas=" + str(blas) + " omp=28 blas=" + str(blas) +
                    " eigen=1 memory_mode=" + memory + " ar_threads=RESTORED_ON_EXIT\n")
                with patch.object(completion, "audit_epoch", return_value={"epoch": validation.EPOCH}):
                    report = validation.output_audit(plan, stages)
                self.assertEqual(report["passed"], blas == 14 and memory == "parallel")

    def test_final_float_rejection_fails_even_with_exit_zero(self):
        with tempfile.TemporaryDirectory() as temp:
            plan = self.fixture(temp, self.good_text + "PPP AUTHORITATIVE MEASUREMENT TRANSACTION REJECTED\n")
            report = validation.output_audit(plan, stages)
            self.assertFalse(report["passed"])
            self.assertIn("Authoritative FLOAT", report["failure_reason"])

    def test_single_finite_epoch_is_only_smoke_not_complete_checkpoint_decode(self):
        with tempfile.TemporaryDirectory() as temp:
            plan = self.fixture(temp, self.good_text)
            with patch.object(completion, "audit_epoch", return_value={"epoch": validation.EPOCH}):
                report = validation.output_audit(plan, stages)
            self.assertTrue(report["passed"])
            self.assertFalse(report["checkpoint_full_kf_decode_claimed"])
            self.assertFalse(report["accuracy_acceptance"])

    def test_wrong_epoch_or_duplicate_atomic_bundle_fails(self):
        for duplicate in (False, True):
            with tempfile.TemporaryDirectory() as temp:
                plan = self.fixture(temp, self.good_text)
                if duplicate:
                    (Path(temp) / "outputs/zhang_internal_products.csv.epochs/extra").mkdir()
                with patch.object(completion, "audit_epoch", return_value={"epoch": "2024-07-17 00:00:30"}):
                    self.assertFalse(validation.output_audit(plan, stages)["passed"])

    def test_resource_stop_and_unreaped_child_are_not_success(self):
        for change in ({"resource_stop_reasons": ["LIMIT"]}, {"child_reaped": False},
                       {"owned_processes_still_live": [{"pid": 123}]}, {"exit_code": -15}):
            with tempfile.TemporaryDirectory() as temp:
                plan = self.fixture(temp, self.good_text)
                plan.update(change)
                self.assertFalse(validation.output_audit(plan, stages)["passed"])


if __name__ == "__main__":
    unittest.main()
