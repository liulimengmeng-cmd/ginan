"""180-station execution-safety contracts; not scientific result acceptance."""
import copy
import json
import os
import signal
import subprocess
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import run_network180_float as network
import run_stages as stages


class Network180Contracts(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.names = [f"A{i:03d}00XXX_R_20241990000_01D_30S_MO.rnx" for i in range(180)]
        self.audit = {"source_station_config": {"sha256": network.STATION_CONFIG_SHA256},
            "station_headers": [{"filename": name} for name in self.names], "missing_files": [],
            "selected_12_station_candidates": [name[:4] for name in self.names[:12]],
            "product_membership": {"sp3_satellites": ["G02", "E02", "C19"],
                                   "clk_satellites": ["G02", "E02", "C19"]},
            "whitelist": [{"prn": "C19"}], "reference_stations": {"G": "MATE", "E": "MATE", "C": "MAL2"}}
        self.station_yaml = self.root / "stations.yaml"
        self.station_yaml.write_text("inputs:\n" + "".join("    - " + name + "\n" for name in self.names))
        self.coordinate = {"estimation_parameters": {"receivers": {"global": {
            "pos": {"estimated": [True], "sigma": [1], "process_noise": [0]},
            "pos_rate": {"estimated": [False]}}}},
            "processing_options": {"minimum_constraints": {"enable": False}}}

    def test_original_service_list_not_smoke_or_new_selection(self):
        with patch.object(stages, "file_sha", return_value=network.STATION_CONFIG_SHA256):
            self.assertEqual(network.read_station_selection(self.station_yaml, self.audit), self.names)
            self.station_yaml.write_text("".join("    - " + name + "\n" for name in self.names[:12]))
            with self.assertRaisesRegex(RuntimeError, "180 unique"):
                network.read_station_selection(self.station_yaml, self.audit)

    def test_mismatched_station_manifest_and_header_list_rejected(self):
        with patch.object(stages, "file_sha", return_value="altered"):
            with self.assertRaisesRegex(RuntimeError, "SHA-256"):
                network.read_station_selection(self.station_yaml, self.audit)
        altered = copy.deepcopy(self.audit)
        altered["station_headers"][0]["filename"] = "HOLDOUT.rnx"
        with patch.object(stages, "file_sha", return_value=network.STATION_CONFIG_SHA256):
            with self.assertRaisesRegex(RuntimeError, "header audit"):
                network.read_station_selection(self.station_yaml, altered)

    def test_full_network_schedule_checkpoint_and_no_mutation(self):
        original = copy.deepcopy(self.audit)
        value = network.network_overlay(self.root, Path("/server/input"), Path("/server/nav.rnx"), self.audit, self.names)
        self.assertEqual(self.audit, original)
        epoch = value["processing_options"]["epoch_control"]
        self.assertEqual(epoch["start_epoch"], "2024-07-17 00:00:00")
        self.assertEqual(epoch["end_epoch"], "2024-07-17 02:00:00")
        self.assertEqual(epoch["epoch_interval"], 30)
        self.assertEqual(len(value["inputs"]["gnss_observations"]["rnx_inputs"]), 180)
        product = value["processing_options"]["gnss_general"]["zhang_pppar"]
        self.assertTrue(product["deterministic_checkpoint"])
        self.assertEqual(product["checkpoint_capture_epochs"], network.CHECKPOINT_EPOCHS)
        self.assertEqual(product["checkpoint_restore_path"], "")
        self.assertTrue(all(value["processing_options"]["gnss_general"]["sys_options"][system]["process"]
                            for system in ("gps", "gal", "bds")))
        self.assertEqual(value["inputs"]["satellite_data"]["nav_files"], [str(Path("/server/nav.rnx"))])

    def test_coordinate_review_sha_and_static_original_prior_are_mandatory(self):
        path, review_path = self.root / "coords.yaml", self.root / "coords_review.json"
        path.write_text(json.dumps(self.coordinate))
        review = {"approved_for_static_coordinate_estimation": True, "all_service_receivers_estimated": True,
            "no_pseudo_coordinate_constraints": True, "mode": "STATIC_ESTIMATED",
            "coordinate_overlay_sha256": stages.file_sha(path)}
        review_path.write_text(json.dumps(review))
        self.assertEqual(network.coordinate_approval(path, review_path)[1], self.coordinate)
        changed = copy.deepcopy(self.coordinate)
        changed["estimation_parameters"]["receivers"]["global"]["pos"]["process_noise"] = [100]
        path.write_text(json.dumps(changed))
        with self.assertRaisesRegex(RuntimeError, "SHA-256"):
            network.coordinate_approval(path, review_path)
        review["coordinate_overlay_sha256"] = stages.file_sha(path)
        review_path.write_text(json.dumps(review))
        with self.assertRaisesRegex(RuntimeError, "Q=0"):
            network.coordinate_approval(path, review_path)

    def test_cli_exact_replacement_180_plus_one_and_real_summary(self):
        pinned = {"observations": [{"path": f"/data/{name}", "station": name[:4]} for name in self.names],
                  "navigation": [{"path": "/nav/day.rnx"}]}
        cli = network.exact_input_cli(pinned)
        self.assertEqual(cli[0], "--rnx_inputs")
        self.assertEqual(cli[-2:], ["--nav_files", "/nav/day.rnx"])
        text = "rnx_inputs: " + " ".join(row["path"] for row in pinned["observations"]) + "\nnav_files: /nav/day.rnx\n"
        self.assertTrue(network.verify_effective_summary(text, pinned)["passed"])
        with self.assertRaisesRegex(RuntimeError, "differs"):
            network.verify_effective_summary(text.replace("/data/A000", "/data/USER"), pinned)
        with self.assertRaisesRegex(RuntimeError, "extra"):
            network.verify_effective_summary(text + "ubx_inputs: /extra.ubx\n", pinned)
        pinned["observations"] = pinned["observations"][:12]
        with self.assertRaisesRegex(RuntimeError, "180 unique"):
            network.exact_input_cli(pinned)

    def test_checkpoint_or_host_approval_failure_never_starts_solver(self):
        approval = self.root / "approval.json"
        approval.write_text(json.dumps({"approved_for_execution": True}))
        with patch.object(network.sys, "platform", "linux"), patch.object(network.subprocess, "Popen") as popen:
            with self.assertRaisesRegex(RuntimeError, "checkpoint"):
                network.execute({}, {}, approval)
        popen.assert_not_called()

    def test_real_parser_rejection_is_fail_closed_without_solver(self):
        pinned = {"observations": [{"path": f"/data/{name}", "station": name[:4]} for name in self.names],
                  "navigation": [{"path": "/nav/day.rnx"}], "configurations": []}
        plan = {"run_root": str(self.root), "binary": "/fake/pea", "configs": [],
                "working_directory": str(self.root), "pinned_input_contract": pinned}
        rejected = type("ParserResult", (), {"returncode": 1,
            "stdout": "Infra-0 v1 is restricted to GPS-only E29 processing"})()
        with patch.object(network, "verify_pins"), patch.object(network.subprocess, "run", return_value=rejected), \
             patch.object(network.subprocess, "Popen") as popen:
            with self.assertRaisesRegex(RuntimeError, "rejected"):
                network.preflight(plan, {})
        popen.assert_not_called()
        self.assertEqual(plan["state"], "CONFIGURATION_PREFLIGHT_FAILED_NOT_RUN")

    def test_preparation_does_not_launch_or_preflight_and_refuses_root_reuse(self):
        binary = self.root / "pea"
        binary.write_text("not executable")
        audit_path = self.root / "audit.json"
        audit_path.write_text(json.dumps(self.audit))
        coord, review = self.root / "coordinate.json", self.root / "review.json"
        coord.write_text(json.dumps(self.coordinate))
        review.write_text("{}")
        trial = self.root / "trial"
        pinned = {"observations": [{"path": f"/data/{name}", "station": name[:4]} for name in self.names],
                  "navigation": [{"path": "/nav/day.rnx"}], "configurations": []}
        with patch.object(network, "read_station_selection", return_value=self.names), \
             patch.object(network, "coordinate_approval", return_value=({}, self.coordinate)), \
             patch.object(network, "pin_inputs", return_value=pinned), \
             patch.object(network.subprocess, "Popen") as popen, \
             patch.object(network.subprocess, "run") as run:
            plan, env = network.prepare(binary, audit_path, trial, self.root, self.root, binary, coord, review)
            popen.assert_not_called()
            run.assert_not_called()
            self.assertEqual(plan["state"], "PREPARED_NOT_RUN")
            self.assertEqual(plan["expected_epoch_count"], 241)
            self.assertFalse(plan["ar_execution"])
            self.assertFalse(plan["scientific_acceptance"])
            self.assertFalse(plan["reuse_gps_checkpoint"])
            self.assertEqual(env["ZHANG_R51_AR_START_GPST_SECONDS"], str(stages.epoch_seconds("23:59:30")))
            self.assertEqual(env["ZHANG_MULTI_GNSS_FLOAT_CHECKPOINT"], "1")
            with self.assertRaises(FileExistsError):
                network.prepare(binary, audit_path, trial, self.root, self.root, binary, coord, review)

    def test_runtime_resource_thresholds_are_conservative_and_bounded(self):
        sample = {"owned_tree_rss_bytes": 84 * network.GIB,
                  "host_available_memory_bytes": 12 * network.GIB,
                  "output_free_disk_bytes": 20 * network.GIB}
        self.assertEqual(network.resource_violations(sample), [])
        sample.update(owned_tree_rss_bytes=85 * network.GIB,
                      host_available_memory_bytes=12 * network.GIB - 1,
                      output_free_disk_bytes=20 * network.GIB - 1)
        self.assertEqual(len(network.resource_violations(sample)), 3)
        self.assertLessEqual(network.RESOURCE_LIMITS["sample_interval_seconds"], 30)

    def test_only_final_current_authoritative_failure_is_network_fatal(self):
        self.assertFalse(network.authoritative_failure("QC failed; reweight and retry 1 of 20"))
        self.assertFalse(network.authoritative_failure("HOB2 SPP failed: insufficient satellites"))
        self.assertFalse(network.authoritative_failure("ZHANG_KF_MEASUREMENT_TRANSACTION status=ROLLED_BACK reason=SPP_RETRY"))
        outputs = self.root / "outputs"
        outputs.mkdir()
        (self.root / "runner.log").write_text("QC retry\n")
        trace = outputs / "Network-test.TRACE"
        trace.write_text("ordinary QC retries\n")
        reader = network.CurrentRunLogTail(self.root)
        self.assertEqual(reader.final_failures(), [])
        with trace.open("a") as stream:
            stream.write("PPP AUTHORITATIVE MEASUREMENT TRANSACTION REJECTED result=2 reason=FINAL_POSTFIT_QC_NOT_CONVERGED")
        failures = reader.final_failures()
        self.assertEqual(len(failures), 1)
        self.assertEqual(Path(failures[0]["path"]), trace)
        self.assertEqual(reader.final_failures(), [])

    @unittest.skipUnless(hasattr(os, "killpg"), "Native Linux-only process group supervisor")
    def test_supervisor_sigint_first_then_bounded_owned_escalation(self):
        child = type("Child", (), {"pid": 12345, "returncode": 137})()
        child.poll = unittest.mock.Mock(side_effect=[None, None, None, 137])
        child.wait = unittest.mock.Mock(side_effect=[subprocess.TimeoutExpired("test", 120),
            subprocess.TimeoutExpired("test", 10), 137])
        with patch.object(network, "proc_identity", return_value={"pid": 12345, "start_ticks": 99}), \
             patch.object(network.os, "getpgid", return_value=12345), patch.object(network.os, "killpg") as killed:
            result = network.stop_owned_child(child, {"start_ticks": 99})
        self.assertEqual([call.args for call in killed.call_args_list],
            [(12345, signal.SIGINT), (12345, signal.SIGTERM), (12345, signal.SIGKILL)])
        self.assertEqual([row["maximum_wait_seconds"] for row in result["signals"]], [120, 10, 10])

    @unittest.skipUnless(hasattr(os, "killpg"), "Native Linux-only process group supervisor")
    def test_supervisor_never_signals_recycled_or_nonowned_pid(self):
        child = type("Child", (), {"pid": 12345, "poll": lambda self: None})()
        with patch.object(network, "proc_identity", return_value={"pid": 12345, "start_ticks": 100}), \
             patch.object(network.os, "getpgid", return_value=12345), patch.object(network.os, "killpg") as killed:
            with self.assertRaisesRegex(RuntimeError, "changed/recycled/non-owned"):
                network.stop_owned_child(child, {"start_ticks": 99})
        killed.assert_not_called()


if __name__ == "__main__":
    unittest.main()
