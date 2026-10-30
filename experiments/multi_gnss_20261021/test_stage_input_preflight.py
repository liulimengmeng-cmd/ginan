"""Contracts for explicit input replacement; mocked parsing is not real PEA evidence."""
import copy
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

import run_stages as stages


class InputReplacementTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        (self.root / "data").mkdir()
        self.audit = {"selected_12_station_candidates": [], "station_headers": []}
        for number in range(12):
            station = f"A{number:03d}"
            filename = station + "00XXX_R_20241990000_01D_30S_MO.rnx"
            (self.root / "data" / filename).write_bytes(("SYNTHETIC_FILE_IDENTITY_" + station).encode())
            self.audit["selected_12_station_candidates"].append(station)
            self.audit["station_headers"].append({"filename": filename})
        self.nav = self.root / "mixed_nav3.rnx"
        self.nav.write_bytes(b"SYNTHETIC_NAV_FILE_IDENTITY_NOT_NAV_VALIDATION")
        self.config = self.root / "physical_model.yaml"
        self.config.write_bytes(b"keep_original_model: true\n")
        self.source_config_bytes = self.config.read_bytes()
        self.input_patch = patch.object(stages, "NATIVE_INPUT", self.root)
        self.nav_patch = patch.object(stages, "DOWNLOADED_NAV3", self.nav)
        self.input_patch.start()
        self.nav_patch.start()
        self.addCleanup(self.input_patch.stop)
        self.addCleanup(self.nav_patch.stop)
        self.pinned = stages.pin_stage_inputs(self.audit, [str(self.config)])

    def summary(self):
        return ("Final effective configuration\n\trnx_inputs: " +
            " ".join(v["path"] for v in self.pinned["observations"]) + "\n\tnav_files: " +
            self.pinned["navigation"][0]["path"] + "\n" +
            "Dry-run successful: configuration and sanity checks passed\n")

    def test_pin_exact_files_and_preserve_existing_physical_models(self):
        self.assertEqual(len(self.pinned["observations"]), 12)
        self.assertEqual(len({v["path"] for v in self.pinned["observations"]}), 12)
        self.assertEqual(len(self.pinned["navigation"]), 1)
        self.assertTrue(all(len(v["sha256"]) == 64 for v in self.pinned["observations"]))
        self.assertEqual(self.config.read_bytes(), self.source_config_bytes)
        self.assertFalse(self.pinned["scientific_acceptance"])
        self.assertFalse(self.pinned["cnv123_propagation_claim"])

    def test_command_line_has_explicit_absolute_replacement_flags(self):
        cli = stages.exact_input_cli(self.pinned)
        self.assertEqual(cli[0], "--rnx_inputs")
        self.assertEqual(cli[13], "--nav_files")
        self.assertTrue(all(Path(p).is_absolute() for p in cli[1:13]))
        self.assertEqual(cli[14], str(self.nav.resolve()))

    def test_b_to_g_execution_keeps_epoch_info_logs_and_preflight_has_no_solver_run(self):
        command = stages.independent_stage_command("pea", [str(self.config)], self.pinned)
        self.assertNotIn("-q", command)
        self.assertNotIn("--dry-run", command)
        self.assertIn("--rnx_inputs", command)
        preflight = stages.independent_stage_command("pea", [str(self.config)], self.pinned, dry_run=True)
        self.assertEqual(preflight, command + ["--dry-run"])

    def test_real_summary_shape_requires_exact_twelve_and_one(self):
        result = stages.verify_effective_input_summary(self.summary(), self.pinned)
        self.assertTrue(result["passed"])
        self.assertEqual(result["actual_unique_observation_files"], 12)
        self.assertEqual(result["actual_navigation_files"], 1)

    def test_180_station_merged_summary_is_rejected(self):
        text = self.summary().replace("\n\tnav_files:", " " + " ".join(f"/extra/station{i}.rnx" for i in range(168)) + "\n\tnav_files:")
        with self.assertRaisesRegex(RuntimeError, "effective input mismatch"):
            stages.verify_effective_input_summary(text, self.pinned)

    def test_duplicate_input_or_second_nav_is_rejected_even_if_unique_set_matches(self):
        duplicate = self.summary().replace("\n\tnav_files:", " " + self.pinned["observations"][0]["path"] + "\n\tnav_files:")
        with self.assertRaises(RuntimeError):
            stages.verify_effective_input_summary(duplicate, self.pinned)
        with self.assertRaises(RuntimeError):
            stages.verify_effective_input_summary(self.summary().replace(str(self.nav), str(self.nav) + " /extra/nav4.rnx"), self.pinned)

    def test_other_streams_cannot_bypass_station_or_nav_selection(self):
        with self.assertRaisesRegex(RuntimeError, "additional observation/navigation"):
            stages.verify_effective_input_summary(self.summary() + "\trtcm_inputs: ntrip://not_selected\n", self.pinned)

    def test_input_hash_mutation_fails_before_subprocess(self):
        Path(self.pinned["observations"][0]["path"]).write_bytes(b"CHANGED")
        with patch.object(stages.subprocess, "run") as run:
            result = stages.configuration_preflight("pea", [str(self.config)], self.pinned, self.root, {})
        self.assertFalse(result["passed"])
        self.assertIn("Pinned input changed", result["reason"])
        run.assert_not_called()

    def test_preflight_is_only_dry_run_and_requires_explicit_success_marker(self):
        completed = subprocess.CompletedProcess([], 0, self.summary())
        with patch.object(stages.subprocess, "run", return_value=completed) as run:
            result = stages.configuration_preflight("pea", [str(self.config)], self.pinned, self.root, {})
        self.assertTrue(result["passed"])
        command = run.call_args.args[0]
        self.assertIn("--dry-run", command)
        self.assertNotIn("--dump-config-only", command)
        self.assertNotIn("-q", command)
        self.assertEqual(run.call_args.kwargs["cwd"], self.root)
        self.assertTrue((self.root / "config_preflight.log").is_file())
        self.assertEqual(len(result["dump_sha256"]), 64)
        completed.stdout = self.summary().replace("Dry-run successful: configuration and sanity checks passed", "")
        with patch.object(stages.subprocess, "run", return_value=completed):
            result = stages.configuration_preflight("pea", [str(self.config)], self.pinned, self.root, {})
        self.assertFalse(result["passed"])

    def test_nonzero_dry_run_not_overridden_by_success_looking_summary(self):
        with patch.object(stages.subprocess, "run", return_value=subprocess.CompletedProcess([], 1, self.summary())):
            result = stages.configuration_preflight("pea", [str(self.config)], self.pinned, self.root, {})
        self.assertFalse(result["passed"])

    def test_raw_yaml_section_is_not_effective_summary(self):
        with self.assertRaises(RuntimeError):
            stages.verify_effective_input_summary("rnx_inputs:\n - station.rnx\nnav_files:\n - old.rnx\n", self.pinned)

    def test_selection_cannot_silently_duplicate_or_drop_station(self):
        malformed = copy.deepcopy(self.audit)
        malformed["selected_12_station_candidates"][-1] = malformed["selected_12_station_candidates"][0]
        with self.assertRaises(RuntimeError):
            stages.pin_stage_inputs(malformed, [str(self.config)])


if __name__ == "__main__":
    unittest.main()
