"""Independent numeric reader contracts; no production execution."""
import json
from pathlib import Path
import tempfile
import unittest

import analyze_clock_lifecycle as analysis


class NumericClockEvidenceTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def write_events(self, events):
        with (self.root / "clock_lifecycle_events.jsonl").open("w") as stream:
            for number, event in enumerate(events, 1):
                stream.write(json.dumps({"source_line": number * 10, "event": event}) + "\n")

    def test_measured_offset_reapplication_is_independently_recomputed(self):
        self.write_events([{"phase": "SPP_OFFSET_SAVE_AFTER", "receiver": "NNOR", "epoch": "t0",
            "saved_offset_m": 650, "spp_clock_m": 60},
            {"phase": "REC_CLOCK_ADJUST_AFTER", "receiver": "NNOR", "epoch": "t1",
             "spp_clock_m": 714, "saved_offset_m": 650, "predicted_clock_m": 710,
             "unrounded_adjustment_m": 654, "applied_adjustment_m": 654}])
        result = analysis.galileo_chain(self.root)["receiver_chains"][0]
        arithmetic = result["independent_scalar_arithmetic"][0]
        self.assertEqual(arithmetic["spp_plus_offset_minus_prediction_m"], 654)
        self.assertEqual(arithmetic["arithmetic_residual_m"], 0)
        self.assertEqual(arithmetic["same_chart_no_saved_offset_counterfactual_m"], 4)
        self.assertEqual(result["observed_adjustment_results"][0]["runner_source_line"], 20)

    def test_full_qc_abort_is_not_timeout_and_missing_processed_is_unavailable(self):
        (self.root / "outputs").mkdir()
        (self.root / "stage_plan.json").write_text(json.dumps({"exit_code": 1, "timed_out": False, "elapsed_seconds": 40}))
        (self.root / "outputs/Network-test.TRACE").write_text("ZHANG_KF_MEASUREMENT_TRANSACTION time=t0 status=ROLLED_BACK reason=FINAL_POSTFIT_QC_NOT_CONVERGED\n")
        report = analysis.outcome(self.root)
        self.assertEqual(report["classification"], "QC_TRANSACTION_ABORT")
        self.assertEqual(report["logs"]["processed_epochs"]["status"], "UNAVAILABLE")
        self.assertEqual(report["spp_offset_chain"]["status"], "UNAVAILABLE")
        (self.root / "stage_plan.json").write_text(json.dumps({"exit_code": -15, "timed_out": True}))
        self.assertEqual(analysis.outcome(self.root)["classification"], "TIMEOUT")

    def test_explicit_single_and_failed_records_remain_distinct(self):
        rows = [{"source_line": i, "event": {"phase": "SPP_RESULT", "epoch": "t0", "receiver": receiver,
            "spp_system": "GAL", "spp_status": status}} for i, (receiver, status) in
            enumerate((("MATE", "SINGLE"), ("HOB2", "FAILED")), 1)]
        result = analysis.spp_result_evidence(rows)
        self.assertEqual(result["counts_by_epoch_selected_system"], {"t0|GAL": {"SINGLE": 1, "FAILED": 1}})


if __name__ == "__main__":
    unittest.main()
