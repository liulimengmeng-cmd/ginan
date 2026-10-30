"""Final replay reader contracts: never launch or accept a numerical solver."""
import json
from pathlib import Path
import tempfile
import unittest

import audit_final_real_replays as audit


class FinalReplayAuditTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def write_events(self, events):
        with (self.root / "clock_lifecycle_events.jsonl").open("w") as stream:
            for number, event in enumerate(events, 1):
                event.update(schema="ZHANG_CLOCK_LIFECYCLE_V1")
                stream.write(json.dumps({"source_line": number, "event": event}) + "\n")

    def test_failed_spp_cold_seed_is_not_called_an_actual_zero_clock(self):
        self.write_events([{"phase": "COLD_RECEIVER_CLOCK_SEED", "epoch": "2024-07-17 00:00:00",
            "receiver": "HOB2", "system": "GAL", "spp_system": "GAL", "spp_status": "FAILED",
            "spp_clock_m": 0, "init_clock_m": 650}])
        row = audit.clock_lifecycle(self.root)["cold_receiver_seeds"][0]
        self.assertIsNone(row["actual_external_system_spp_m"])
        self.assertIsNone(row["inferred_translation_m"])

    def test_offset_checks_recompute_formula_and_match_previous_committed_save(self):
        self.write_events([{"phase": "SPP_OFFSET_SAVE_AFTER", "epoch": "2024-07-17 00:00:00",
            "receiver": "NNOR", "system": "GAL", "saved_system_offset_m": 650, "spp_system_clock_m": 60, "external_chart": "external"},
            {"phase": "REC_CLOCK_ADJUST_AFTER", "epoch": "2024-07-17 00:00:30", "receiver": "NNOR", "system": "GAL",
             "spp_system_clock_m": 62, "saved_system_offset_m": 650, "predicted_clock_m": 710,
             "unrounded_adjustment_m": 2, "applied_adjustment_m": 2, "external_chart": "external"}])
        report = audit.clock_lifecycle(self.root)["offset_checks"]
        self.assertEqual(report["formula_checked_adjustments"], 1)
        self.assertEqual(report["matching_previous_saved_offset_count"], 1)
        self.assertEqual(report["maximum_formula_residual_m"], 0)
        self.assertEqual(report["matching_external_descriptor_count"], 1)

    def test_actual_qc_grid_and_terminal_rollback_are_distinct(self):
        (self.root / "outputs").mkdir()
        (self.root / "outputs/Network-test.TRACE").write_text(
            "ZHANG_KF_MEASUREMENT_TRANSACTION time=2024-07-17 00:00:00 boundary=B13 status=COMMITTED transaction_id=x|/PPP|1\n"
            "ZHANG_KF_MEASUREMENT_TRANSACTION time=2024-07-17 00:00:30 boundary=B13 status=ROLLED_BACK transaction_id=x|/PPP|2 reason=bad\n")
        report = audit.qc_transactions(self.root)
        self.assertEqual(report["committed_unique_epoch_count"], 1)
        self.assertEqual(report["rollback_count"], 1)

    def test_running_trial_is_not_a_completed_final_result(self):
        (self.root / "stage_plan.json").write_text(json.dumps({"state": "RUNNING"}))
        report = audit.audit_run(self.root, 61)
        self.assertFalse(report["completed"])
        self.assertFalse(report["scientific_acceptance"])


if __name__ == "__main__":
    unittest.main()
