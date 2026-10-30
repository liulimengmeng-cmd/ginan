"""Contract tests only: never launch a numerical PEA process."""
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import run_float_rank_diagnostic_equivalence as experiment


class FloatRankDiagnosticComparisonContracts(unittest.TestCase):
    def log(self, text):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        path = Path(temporary.name) / "runner.log"
        path.write_text(text)
        return path

    def test_arms_change_only_one_narrow_environment_variable(self):
        base = {"ZHANG_R51_RATIO_THRESHOLD": "2", "ZHANG_R51_RATIO_ONLY": "1",
                "ZHANG_R51_EPOCH_ONLY_AR": "1", "OMP_NUM_THREADS": "1"}
        before = dict(base)
        for value in ("full", "off", None):
            actual = experiment.arm_environment(base, value)
            stripped = dict(actual)
            stripped.pop(experiment.VARIABLE, None)
            self.assertEqual(stripped, before)
            if value is None:
                self.assertNotIn(experiment.VARIABLE, actual)
            else:
                self.assertEqual(actual[experiment.VARIABLE], value)
        self.assertEqual(base, before)
        with self.assertRaises(ValueError):
            experiment.arm_environment(base, "on")

    def test_unset_removes_inherited_diagnostic_optin(self):
        original = {experiment.VARIABLE: "full", "SENTINEL": "keep"}
        self.assertEqual(experiment.arm_environment(original, None), {"SENTINEL": "keep"})
        self.assertEqual(original[experiment.VARIABLE], "full")

    def test_full_requires_actual_svd_start_completion_and_evaluated_rank(self):
        text = ("ZHANG_FLOAT_PHASE_BEGIN phase=FLOAT_RANK_JACOBI_SVD epoch=x\n"
                "ZHANG_FLOAT_PHASE_END phase=FLOAT_RANK_JACOBI_SVD epoch=x inclusive_ms=1\n"
                "ZHANG_PURE_OBS_RANK time=x rows=2 rank=1\n"
                "ZHANG_PURE_OBS_RANK_STATUS mode=full status=COMPLETED epoch=x rank=EVALUATED\n"
                "Processed epoch #1 GPS time: 2323 259200 2024-07-17 00:00:00\n")
        result = experiment.diagnostic_evidence(self.log(text), "full")
        self.assertEqual(result["svd_starts"], 1)
        for remove in ("PHASE_BEGIN", "PHASE_END", "RANK time=", "status=COMPLETED", "Processed epoch"):
            broken = "\n".join(line for line in text.splitlines() if remove not in line)
            with self.assertRaises(RuntimeError):
                experiment.diagnostic_evidence(self.log(broken), "full")

    def test_off_and_default_are_not_a_fabricated_rank_success(self):
        text = ("ZHANG_PURE_OBS_RANK_STATUS mode=off status=DISABLED epoch=x rank=NOT_EVALUATED\n"
                "Processed epoch #1 GPS time: 2323 259200 2024-07-17 00:00:00\n")
        for mode in ("off", None):
            self.assertEqual(experiment.diagnostic_evidence(self.log(text), mode)["svd_starts"], 0)
            with self.assertRaises(RuntimeError):
                experiment.diagnostic_evidence(self.log(text.replace("NOT_EVALUATED", "EVALUATED")), mode)
            with self.assertRaises(RuntimeError):
                experiment.diagnostic_evidence(self.log(text +
                    "ZHANG_FLOAT_PHASE_BEGIN phase=FLOAT_RANK_JACOBI_SVD epoch=x\n"), mode)

    @staticmethod
    def valid_report():
        keys = ("core.x", "core.P", "core.dx", "core.kfIndexMap", "core.stateTransitionMap",
                "core.procNoiseMap", "core.initNoiseMap", "core.sigmaPass", "core.chiQCPass",
                "core.qc", "sections.inventory", "core.authoritative_spp_offset_all_GEC_present",
                "core.static_positions_all_XYZ_positive_variance", "core.no_receiver_position_rate")
        return {"passed": True, "checks": dict.fromkeys(keys, True),
            "coordinate_components_by_receiver": {"MAL2": [0, 1, 2], "MATE": [0, 1, 2]},
            "numeric_P": {"all_finite": True, "exact_equal": True, "maximum_absolute_difference_native": 0},
            "sections": [{"name": "private", "entire_payload_byte_identical": True}]}

    def test_full_state_order_qc_maps_and_complete_covariance_must_all_match(self):
        report = self.valid_report()
        self.assertTrue(experiment.verify_checkpoint_comparison(report, {"MAL2", "MATE"}))
        for key in report["checks"]:
            changed = self.valid_report()
            changed["checks"][key] = False
            with self.assertRaises(RuntimeError):
                experiment.verify_checkpoint_comparison(changed, {"MAL2", "MATE"})
        for field, value in (("all_finite", False), ("exact_equal", False),
                             ("maximum_absolute_difference_native", 1e-99)):
            changed = self.valid_report()
            changed["numeric_P"][field] = value
            with self.assertRaises(RuntimeError):
                experiment.verify_checkpoint_comparison(changed, {"MAL2", "MATE"})

    def test_missing_state_checks_sections_or_receiver_is_failure(self):
        changed = self.valid_report()
        changed["checks"].pop("core.kfIndexMap")
        with self.assertRaises(RuntimeError):
            experiment.verify_checkpoint_comparison(changed, {"MAL2", "MATE"})
        for sections in ([], [{"entire_payload_byte_identical": False}]):
            changed = self.valid_report()
            changed["sections"] = sections
            with self.assertRaises(RuntimeError):
                experiment.verify_checkpoint_comparison(changed, {"MAL2", "MATE"})
        with self.assertRaises(RuntimeError):
            experiment.verify_checkpoint_comparison(self.valid_report(), {"MAL2", "MATE", "EXTRA"})

    def test_duplicate_attempt_is_rejected_before_solver_or_approval(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp).resolve()
            (root / "OWNERSHIP.json").write_text(json.dumps({"schema": experiment.replay.OWNERSHIP,
                                                            "root": str(root)}))
            (root / "rank_execution_started.json").write_text("{}")
            plan = {"state": "RANK_DIAGNOSTIC_PREPARED_NOT_RUN", "run_root": str(root)}
            with patch.object(experiment.sys, "platform", "linux"), \
                 patch.object(experiment.replay, "execution_gate") as gate, \
                 patch.object(experiment.subprocess, "Popen") as popen:
                with self.assertRaisesRegex(RuntimeError, "already been attempted"):
                    experiment.execute(plan, {}, root / "not_read.json")
            gate.assert_not_called()
            popen.assert_not_called()


if __name__ == "__main__":
    unittest.main()
