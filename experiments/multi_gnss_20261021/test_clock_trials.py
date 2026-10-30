"""Driver and reader contracts, not PEA scientific acceptance."""
import copy
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import run_clock_trials as trials
import run_stages as stages


class ClockTrialContracts(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.audit = {"station_headers": [{"filename": f"A{i:03d}.rnx"} for i in range(12)],
            "selected_12_station_candidates": [f"A{i:03d}" for i in range(12)],
            "product_membership": {"sp3_satellites": ["G02", "E02", "C19"],
                                   "clk_satellites": ["G02", "E02", "C19"]},
            "whitelist": [{"prn": "C19"}], "reference_stations": {"G": "MATE", "E": "MATE", "C": "MAL2"}}

    def test_two_epoch_cold_schedule_and_long_controls(self):
        for count, end in ((2, "00:00:30"), (21, "00:10:00"), (251, "02:05:00")):
            self.assertEqual(trials.end_time(count), end)
        for count in (0, 2881):
            with self.assertRaises(ValueError):
                trials.end_time(count)

    def test_all_four_system_selections_and_exact_signals(self):
        old_stages = copy.deepcopy(stages.STAGES)
        old_audit = copy.deepcopy(self.audit)
        for label, selected in trials.SYSTEMS.items():
            value = trials.trial_overlay(label, 2, "configured", self.root, self.audit)
            gnss = value["processing_options"]["gnss_general"]
            self.assertEqual([s for s in stages.SIGNALS if gnss["sys_options"][s]["process"]], selected)
            self.assertEqual(set(gnss["zhang_full_rank"]["sys_options"]), set(selected))
            self.assertEqual(set(gnss["zhang_pppar"]["sys_options"]), set(selected))
            for system in selected:
                self.assertEqual(gnss["sys_options"][system]["code_priorities"], stages.SIGNALS[system])
            self.assertEqual(value["processing_options"]["epoch_control"]["start_epoch"], "2024-07-17 00:00:00")
            self.assertEqual(value["processing_options"]["epoch_control"]["end_epoch"], "2024-07-17 00:00:30")
            self.assertEqual(gnss["zhang_pppar"]["checkpoint_restore_path"], "")
            self.assertFalse(gnss["zhang_pppar"]["deterministic_checkpoint"])
            self.assertNotIn("global", value["satellite_options"])
        self.assertEqual(stages.STAGES, old_stages)
        self.assertEqual(self.audit, old_audit)

    def test_precise_control_changes_only_clock_source_and_description(self):
        configured = trials.trial_overlay("E", 2, "configured", self.root, self.audit)
        precise = trials.trial_overlay("E", 2, "precise-only", self.root, self.audit)
        self.assertEqual(precise["satellite_options"].pop("global"), {"models": {"clock": {"sources": ["PRECISE"]}}})
        precise["outputs"]["metadata"] = configured["outputs"]["metadata"]
        self.assertEqual(configured, precise)
        self.assertNotIn("statistics", precise["processing_options"])
        self.assertNotIn("models", precise["receiver_options"]["global"])

    def test_environment_rejects_checkpoint_and_policy_override(self):
        with patch.dict(trials.os.environ, {"ZHANG_E29_RESTORE_BUNDLE": "old", "ZHANG_R51_RATIO_THRESHOLD": "3"}):
            env = trials.trial_environment(["ZHANG_CLOCK_LIFECYCLE_DIAGNOSTICS=1"])
        self.assertEqual(env["ZHANG_CLOCK_LIFECYCLE_DIAGNOSTICS"], "1")
        self.assertEqual(env["ZHANG_R51_RATIO_THRESHOLD"], "2")
        self.assertEqual(env["ZHANG_R51_AR_START_GPST_SECONDS"], str(stages.epoch_seconds("23:59:30")))
        self.assertFalse(any("RESTORE" in k or "REUSE_FLOAT" in k for k in env))
        for item in ("ZHANG_R51_RATIO_THRESHOLD=3", "ZHANG_E29_RESTORE_BUNDLE=old", "CLOCK=1", "missing"):
            with self.assertRaises(ValueError):
                trials.trial_environment([item])

    def test_spp_adjust_ablation_changes_no_clock_source_or_qc(self):
        configured = trials.trial_overlay("E", 2, "configured", self.root, self.audit)
        disabled = trials.trial_overlay("E", 2, "configured", self.root, self.audit, "off")
        self.assertFalse(disabled["processing_options"]["gnss_general"].pop("adjust_rec_clocks_by_spp"))
        self.assertEqual(configured, disabled)
        self.assertNotIn("global", disabled["satellite_options"])

    def test_legacy_clock_boundaries_are_unavailable_not_zero(self):
        (self.root / "runner.log").write_text("receiver clock adjustment -650m\n")
        report = trials.read_clock_lifecycle(self.root)
        self.assertEqual(report["status"], "UNAVAILABLE")
        self.assertEqual(report["record_count"], 0)
        self.assertTrue(all(row["status"] == "UNAVAILABLE" for row in report["boundaries"].values()))
        self.assertFalse(report["causal_proof"])

    def test_preparation_is_new_directory_only_and_never_launches(self):
        self.audit.update(selected_complete=True, missing_files=[])
        audit = self.root / "audit.json"
        audit.write_text(json.dumps(self.audit))
        binary = self.root / "pea"
        binary.write_bytes(b"not an executable; synthetic preparation test")
        root = self.root / "new_trial"
        with patch.object(stages, "pin_stage_inputs", return_value={
            "observations": [{"path": f"/input/{i}.rnx"} for i in range(12)],
            "navigation": [{"path": "/input/nav.rnx"}]}), \
             patch.object(trials.subprocess, "Popen") as popen, \
             patch.object(stages, "configuration_preflight") as dryrun:
            plan, env = trials.prepare(binary, audit, root, "C")
            popen.assert_not_called()
            dryrun.assert_not_called()
            self.assertEqual(plan["state"], "PREPARED_NOT_RUN")
            self.assertEqual(plan["expected_epoch_count"], 2)
            self.assertFalse(plan["scientific_acceptance"])
            self.assertFalse(plan["reuse_gps_checkpoint"])
            self.assertEqual(len(plan["driver_source_identity"]), 3)
            with self.assertRaises(FileExistsError):
                trials.prepare(binary, audit, root, "C")

    def test_resource_statistics_missing_is_not_zero(self):
        path = self.root / "runner.log"
        self.assertEqual(trials.resource_statistics(path)["status"], "UNAVAILABLE")
        path.write_text("User time (seconds): 3.21\nMaximum resident set size (kbytes): 1234\nSwaps: 0\n")
        report = trials.resource_statistics(path)
        self.assertEqual(report["values"]["maximum_rss_kib"], 1234)
        self.assertEqual(report["values"]["user_cpu_seconds"], 3.21)
        self.assertEqual(report["values"]["swaps"], 0)
        self.assertNotIn("system_cpu_seconds", report["values"])

    def test_reader_real_info_json_shape_retains_null_and_source_line(self):
        event = {"schema": "ZHANG_CLOCK_LIFECYCLE_V1", "phase": "REC_CLOCK_ADJUST_BEFORE",
            "epoch": "2024-07-17 00:00:30.00", "receiver": "MATE", "spp_clock_m": 10,
            "saved_offset_m": 650, "saved_offset_system": "GAL", "applied_adjustment_m": 650,
            "model_clocks": [{"system": "GAL", "kf_clock_m": 10, "kf_rate_mps": None}],
            "observations": [{"satellite": "E02", "source": "KALMAN", "clock_m": 43781.5}]}
        (self.root / "runner.log").write_text("other info\n[INFO] ZHANG_CLOCK_LIFECYCLE " + json.dumps(event) + "\n")
        report = trials.read_clock_lifecycle(self.root)
        self.assertEqual(report["record_count"], 1)
        self.assertEqual(report["satellite_source_counts"], {"KALMAN": 1})
        self.assertEqual(report["samples"][0]["source_line"], 2)
        self.assertIsNone(report["samples"][0]["event"]["model_clocks"][0]["kf_rate_mps"])
        self.assertEqual(report["boundaries"]["reference_transport_full_covariance"]["status"], "UNAVAILABLE")
        self.assertFalse(report["causal_proof"])

    def test_nonfinite_or_unknown_schema_not_valid_boundary_evidence(self):
        (self.root / "runner.log").write_text('ZHANG_CLOCK_LIFECYCLE {"schema":"bad","phase":"SPP_ENTRY"}\n'
            'ZHANG_CLOCK_LIFECYCLE {"schema":"ZHANG_CLOCK_LIFECYCLE_V1","phase":"SPP_ENTRY","spp_clock_m":NaN}\n')
        report = trials.read_clock_lifecycle(self.root)
        self.assertEqual(report["record_count"], 0)
        self.assertEqual(report["invalid_record_count"], 2)

    def test_spp_single_requires_explicit_result_identity_and_retains_failures(self):
        events = [{"schema": "ZHANG_CLOCK_LIFECYCLE_V1", "phase": "SPP_RESULT",
            "epoch": "2024-07-17 00:00:00", "receiver": rec, "spp_system": "GAL", "spp_status": status}
            for rec, status in (("MATE", "SINGLE"), ("HOB2", "FAILED"), ("NNOR", "SINGLE_X"))]
        events.append({"schema": "ZHANG_CLOCK_LIFECYCLE_V1", "phase": "SPP_ENTRY",
            "epoch": "2024-07-17 00:00:00", "receiver": "AGGO", "spp_system": "GAL", "spp_status": "SINGLE"})
        (self.root / "runner.log").write_text("\n".join("ZHANG_CLOCK_LIFECYCLE " + json.dumps(event) for event in events))
        report = trials.read_clock_lifecycle(self.root)["spp_result_status"]
        self.assertEqual(report["unique_epoch_receiver_selected_system_records"], 3)
        self.assertEqual(report["status_counts"], {"SINGLE": 1, "FAILED": 1, "SINGLE_X": 1})
        self.assertEqual({r["selected_spp_system"] for r in report["records"]}, {"GAL"})

    def test_failed_preflight_never_launches_solver(self):
        plan = {"binary": "fake_pea", "binary_sha256": "pin", "configs": [],
            "pinned_input_contract": {}, "run_root": str(self.root)}
        with patch.object(trials.sys, "platform", "linux"), \
             patch.object(stages, "file_sha", return_value="pin"), \
             patch.object(stages, "configuration_preflight", return_value={"passed": False, "reason": "12+1 mismatch"}), \
             patch.object(trials.subprocess, "Popen") as popen:
            with self.assertRaisesRegex(RuntimeError, "12\\+1 mismatch"):
                trials.execute(plan, {})
        popen.assert_not_called()
        self.assertEqual(plan["state"], "CONFIGURATION_PREFLIGHT_FAILED_NOT_RUN")


if __name__ == "__main__":
    unittest.main()
