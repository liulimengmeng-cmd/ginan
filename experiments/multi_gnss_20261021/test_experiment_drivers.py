#!/usr/bin/env python3
"""Contract tests for experiment safety/configuration, not solver science."""
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import audit_inputs as audit
import run_stages as stages
import download_daily_inputs as daily


class StageContractTests(unittest.TestCase):
    def fixture(self):
        names = "MATE MAL2 AGGO HOB2 DARW DGAR KIRI NNOR IISC MAJU YEL2 KAT1".split()
        return {"station_headers": [{"filename": s + "00XXX_R_20241990000_01D_30S_MO.rnx"} for s in names],
                "selected_12_station_candidates": names,
                "product_membership": {"sp3_satellites": ["G02", "G04", "G30", "E05", "E10", "C02", "C19"],
                    "clk_satellites": ["G02", "G04", "G30", "E05", "E10", "C02", "C19"]},
                "whitelist": [{"prn": "C02"}, {"prn": "C19"}],
                "reference_stations": {"G": "MATE", "E": "MATE", "C": "MAL2"}}

    def test_exact_codes_no_pilot_x_substitution(self):
        value = stages.overlay("E", Path("/private/e"), self.fixture())
        gnss = value["processing_options"]["gnss_general"]
        for system in ("gps", "gal", "bds"):
            self.assertEqual(gnss["sys_options"][system]["code_priorities"], stages.SIGNALS[system])
            self.assertEqual(gnss["sys_options"][system]["used_nav_type"], stages.NAV_TYPES[system])
        self.assertEqual(gnss["zhang_full_rank"]["sys_options"]["bds"]["reference_receiver"], "MAL2")
        self.assertEqual(gnss["zhang_full_rank"]["sys_options"]["gal"]["reference_receiver"], "MATE")
        self.assertEqual(value["inputs"]["satellite_data"]["nav_files"], [str(stages.DOWNLOADED_NAV3)])
        self.assertFalse(value["outputs"]["bias_sinex"]["output"])

    def test_bds_metadata_not_prn_threshold(self):
        value = stages.overlay("B_C", Path("/private/b"), self.fixture())
        self.assertFalse(value["satellite_options"]["C02"]["exclude"])
        self.assertFalse(value["satellite_options"]["C19"]["exclude"])
        self.assertTrue(value["satellite_options"]["C20"]["exclude"])
        fullrank = value["processing_options"]["gnss_general"]["zhang_full_rank"]
        self.assertTrue(fullrank["bds3_only"])
        self.assertEqual(fullrank["bds3_satellite_whitelist"], ["C02", "C19"])

    def test_galileo_day_specific_exclusion_replaces_old_list(self):
        value = stages.overlay("B_E", Path("/private/b"), self.fixture())
        self.assertFalse(value["satellite_options"]["E05"]["exclude"])
        self.assertFalse(value["satellite_options"]["E10"]["exclude"])
        self.assertTrue(value["satellite_options"]["G02"]["exclude"])

    def test_other_constellations_disabled(self):
        for stage in ("B_E", "B_C", "C_GE", "C_GC", "D", "E", "G"):
            value = stages.overlay(stage, Path("/private/e"), self.fixture())
            gnss = value["processing_options"]["gnss_general"]
            enabled = {s for s, v in gnss["sys_options"].items() if v["process"]}
            self.assertEqual(enabled, set(stages.STAGES[stage]["systems"]))

    def test_original_gps_exclusions_preserved(self):
        value = stages.overlay("E", Path("/private/e"), self.fixture())
        self.assertTrue(value["satellite_options"]["G04"]["exclude"])
        self.assertTrue(value["satellite_options"]["G30"]["exclude"])

    def test_only_a_can_restore_gps_checkpoint(self):
        with patch.dict(os.environ, {"ZHANG_E29_RESTORE_BUNDLE": "contaminated", "ZHANG_R51_REUSE_FLOAT_20240717": "1"}):
            for stage in stages.STAGES:
                env = stages.environment(stage)
                self.assertEqual(env["ZHANG_R51_RATIO_THRESHOLD"], "2")
                self.assertEqual(env["ZHANG_R51_RATIO_ONLY"], "1")
                self.assertEqual(env["ZHANG_R51_EPOCH_ONLY_AR"], "1")
                self.assertEqual("ZHANG_E29_RESTORE_BUNDLE" in env, stage == "A")

    def test_float_does_not_start_ar_at_2h(self):
        for stage in ("B_E", "B_C", "C_GE", "C_GC", "D"):
            self.assertGreater(int(stages.environment(stage)["ZHANG_R51_AR_START_GPST_SECONDS"]), stages.epoch_seconds("02:05:00"))
        self.assertEqual(int(stages.environment("E")["ZHANG_R51_AR_START_GPST_SECONDS"]), 1405216800)

    def test_approval_cannot_be_missing_or_changed_binary(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            binary = root / "binary"
            binary.write_bytes(b"candidate-v1")
            with self.assertRaises(RuntimeError):
                stages.approval(None, binary, "E")
            approval = root / "approval.json"
            approval.write_text(json.dumps({"build_passed": True, "mandatory_tests_passed": True,
                "approved_for_execution": True, "binary_sha256": stages.file_sha(binary)}))
            with self.assertRaises(RuntimeError):
                stages.approval(approval, binary, "E")
            self.assertTrue(stages.approval(approval, binary, "A"))
            binary.write_bytes(b"candidate-v2")
            with self.assertRaises(RuntimeError):
                stages.approval(approval, binary, "A")

    def test_overlay_outputs_are_private_and_checkpoint_disabled(self):
        value = stages.overlay("E", Path("/private/e"), self.fixture())
        ppp = value["processing_options"]["gnss_general"]["zhang_pppar"]
        self.assertEqual(value["outputs"]["outputs_root"], str(Path("/private/e/outputs")))
        self.assertFalse(ppp["deterministic_checkpoint"])
        self.assertEqual(ppp["checkpoint_restore_path"], "")
        self.assertFalse(value["processing_options"]["ambiguity_resolution"]["fix_and_hold"])


class AuditContractTests(unittest.TestCase):
    def test_rinex_padded_satellite_and_undated_event_are_valid(self):
        def line(content, label):
            return content.ljust(60) + label + "\n"
        text = line("     3.04           OBSERVATION DATA    M", "RINEX VERSION / TYPE")
        text += line("G    4 C1C L1C C2W L2W", "SYS / # / OBS TYPES")
        text += line("", "END OF HEADER")
        text += "> 2024 07 17 00 00  0.0000000  0  1\n"
        text += "G 6" + "".join(f"{value:14.3f}  " for value in (20000000, 100000, 20000001, 100001)) + "\n"
        text += ">                              4  1\n"
        text += line("end statistics", "COMMENT")
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "test.rnx"
            path.write_text(text)
            result = daily.observation_audit(path, set())
            self.assertEqual(result["epoch_count"], 1)
            self.assertEqual(result["parse_error_count"], 0)
            self.assertEqual(result["actual_satellites_by_system"]["G"], ["G06"])
            self.assertEqual(result["target_nonzero_code_phase_pair_counts"]["G"], 1)
            self.assertEqual(result["event_flags"]["4"], 1)

    def test_metadata_generation_from_identity_and_effective_prn(self):
        lines = ["+SATELLITE/IDENTIFIER", " C201 2017-069A 43001 BDS-3M-CAST",
                 " C001 2007-011A 31115 BDS-2M", "+SATELLITE/PRN",
                 " C201 2018:319:00000 0000:000:00000 C02",
                 " C201 2017:309:00000 2018:164:00000 C19",
                 " C001 2010:001:00000 0000:000:00000 C58", "-SATELLITE/PRN"]
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "metadata.snx"
            path.write_text("\n".join(lines) + "\n")
            result = audit.metadata_audit(path)
            self.assertEqual([row["prn"] for row in result], ["C02"])
            self.assertEqual(result[0]["svn"], "C201")
            self.assertEqual(len(result[0]["identifier_line_sha256"]), 64)
            self.assertEqual(len(result[0]["prn_line_sha256"]), 64)

    def test_exact_signal_header_is_not_x_combined(self):
        def line(content, label):
            return content.ljust(60) + label + "\n"
        text = line("     3.04           OBSERVATION DATA    M", "RINEX VERSION / TYPE")
        text += line("E    4 C1X L1X C5X L5X", "SYS / # / OBS TYPES")
        text += line("C    4 C1P L1P C5P L5P", "SYS / # / OBS TYPES")
        text += line("", "END OF HEADER")
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "station.rnx"
            path.write_text(text)
            result = audit.read_header(path)
            self.assertFalse(result["header_target_support"]["E"])
            self.assertTrue(result["header_target_support"]["C"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
