#!/usr/bin/env python3

from __future__ import annotations

import io
import json
import sys
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from validate_stec_satellite_difference_covariance import main, validate_files


RAW_HEADER = """# GINAN_STEC_COVARIANCE_V3
# META,gps_week,gps_tow,status,state_count,upper_triangle_count,max_abs_asymmetry_tecu2,posterior_stage,ar_routine_invoked,ar_eligible_ambiguity_count,ar_integer_ambiguity_coordinate_count,ar_receiver_single_difference_applied,ar_receiver_datum_group_count,ar_dropped_singleton_group_count,ar_resolved_combination_count,ar_pseudoobservations_submitted,ar_mode,ar_configured_success_rate_threshold,ar_configured_solution_ratio_threshold,ar_diagnostic_status,ar_selected_decorrelated_ambiguity_count,ar_integer_candidate_count,ar_bootstrapped_success_rate,ar_best_squared_norm,ar_second_squared_norm,ar_solution_ratio
# STATE,gps_week,gps_tow,local_index,site,satellite,state_number,filter_index,estimate_tecu,variance_tecu2
# COV,gps_week,gps_tow,row_local_index,column_local_index,covariance_tecu2
"""

SD_HEADER = """# GINAN_STEC_SATELLITE_DIFFERENCE_COVARIANCE_V1
# META,gps_week,gps_tow,status,source_state_count,difference_state_count,datum_count,singleton_datum_count,upper_triangle_count,max_abs_asymmetry_tecu2,posterior_stage
# DATUM,gps_week,gps_tow,datum_index,site,constellation,state_number,reference_satellite,member_count,status
# TRANSFORM,gps_week,gps_tow,difference_local_index,source_local_index,coefficient
# SD_STATE,gps_week,gps_tow,local_index,datum_index,site,constellation,target_satellite,reference_satellite,state_number,target_source_local_index,reference_source_local_index,estimate_tecu,variance_tecu2
# SD_COV,gps_week,gps_tow,row_local_index,column_local_index,covariance_tecu2
"""

RAW_VALID = """META,2323,12345.5,OK,3,6,0,FILTER_POSTERIOR,0,0,0,0,0,0,0,0,OFF,0.9999,3,NOT_RUN,0,0,-1,-1,-1,-1
STATE,2323,12345.5,0,A,G01,0,7,10,4
STATE,2323,12345.5,1,A,G02,0,8,12,9
STATE,2323,12345.5,2,A,G03,0,9,15,16
COV,2323,12345.5,0,0,4
COV,2323,12345.5,0,1,1
COV,2323,12345.5,0,2,0.5
COV,2323,12345.5,1,1,9
COV,2323,12345.5,1,2,2
COV,2323,12345.5,2,2,16
"""

SD_VALID = """META,2323,12345.5,OK,3,2,1,0,3,0,FILTER_POSTERIOR
DATUM,2323,12345.5,0,A,G,0,G01,3,DIFFERENCES_CREATED
TRANSFORM,2323,12345.5,0,1,1
TRANSFORM,2323,12345.5,0,0,-1
TRANSFORM,2323,12345.5,1,2,1
TRANSFORM,2323,12345.5,1,0,-1
SD_STATE,2323,12345.5,0,0,A,G,G02,G01,0,1,0,2,11
SD_STATE,2323,12345.5,1,0,A,G,G03,G01,0,2,0,5,19
SD_COV,2323,12345.5,0,0,11
SD_COV,2323,12345.5,0,1,4.5
SD_COV,2323,12345.5,1,1,19
"""


class SatelliteDifferenceCovarianceTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary_directory.name)
        self.raw = self.root / "raw.STEC.COV"
        self.sd = self.root / "sd.STEC.SD.COV"
        self.output = self.root / "audit.json"
        self.write_case(RAW_VALID, SD_VALID)

    def tearDown(self) -> None:
        self.temporary_directory.cleanup()

    def write_case(self, raw_payload: str, sd_payload: str) -> None:
        self.raw.write_text(RAW_HEADER + raw_payload, encoding="utf-8")
        self.sd.write_text(SD_HEADER + sd_payload, encoding="utf-8")

    def test_accepts_complete_d_x_and_d_c_dt_epoch(self) -> None:
        report = validate_files(self.raw, self.sd)
        self.assertTrue(report["all_valid"])
        epoch = report["epochs"][0]
        self.assertEqual(epoch["datum_count"], 1)
        self.assertEqual(epoch["transform_term_count"], 4)
        self.assertEqual(epoch["sd_state_transform_max_abs_error_tecu"], 0)
        self.assertEqual(epoch["sd_covariance_transform_max_abs_error_tecu2"], 0)
        self.assertGreaterEqual(epoch["sd_minimum_eigenvalue_tecu2"], 0)

    def test_rejects_incorrect_sd_state_and_diagonal(self) -> None:
        corrupted = SD_VALID.replace(
            "SD_STATE,2323,12345.5,0,0,A,G,G02,G01,0,1,0,2,11",
            "SD_STATE,2323,12345.5,0,0,A,G,G02,G01,0,1,0,2.25,12",
        )
        self.write_case(RAW_VALID, corrupted)
        report = validate_files(self.raw, self.sd)
        self.assertFalse(report["all_valid"])
        errors = report["epochs"][0]["errors"]
        self.assertTrue(any("sd_state_not_equal_D_x" in error for error in errors))
        self.assertTrue(any("sd_state_variance_diagonal_mismatch" in error for error in errors))

    def test_rejects_wrong_full_off_diagonal_covariance(self) -> None:
        self.write_case(RAW_VALID, SD_VALID.replace("0,1,4.5", "0,1,0.0"))
        report = validate_files(self.raw, self.sd)
        self.assertFalse(report["all_valid"])
        self.assertTrue(
            any(
                "sd_covariance_not_equal_D_C_Dt" in error
                for error in report["epochs"][0]["errors"]
            )
        )

    def test_rejects_transform_without_two_distinct_plus_minus_one_terms(self) -> None:
        corrupted = SD_VALID.replace(
            "TRANSFORM,2323,12345.5,0,0,-1",
            "TRANSFORM,2323,12345.5,0,0,-2",
        )
        self.write_case(RAW_VALID, corrupted)
        report = validate_files(self.raw, self.sd)
        self.assertFalse(report["all_valid"])
        self.assertTrue(
            any(
                "must_have_two_distinct_plus_minus_one_terms" in error
                for error in report["epochs"][0]["errors"]
            )
        )

    def test_rejects_nondeterministic_datum_reference(self) -> None:
        corrupted = SD_VALID.replace(
            "DATUM,2323,12345.5,0,A,G,0,G01,3,DIFFERENCES_CREATED",
            "DATUM,2323,12345.5,0,A,G,0,G02,3,DIFFERENCES_CREATED",
        )
        self.write_case(RAW_VALID, corrupted)
        report = validate_files(self.raw, self.sd)
        self.assertFalse(report["all_valid"])
        self.assertTrue(
            any(
                "deterministic_definition_mismatch" in error
                for error in report["epochs"][0]["errors"]
            )
        )

    def test_rejects_matching_but_non_psd_covariance(self) -> None:
        raw = """META,2323,12345.5,OK,3,6,0,FILTER_POSTERIOR,0,0,0,0,0,0,0,0,OFF,0.9999,3,NOT_RUN,0,0,-1,-1,-1,-1
STATE,2323,12345.5,0,A,G01,0,7,10,1
STATE,2323,12345.5,1,A,G02,0,8,12,1
STATE,2323,12345.5,2,A,G03,0,9,15,1
COV,2323,12345.5,0,0,1
COV,2323,12345.5,0,1,0
COV,2323,12345.5,0,2,0
COV,2323,12345.5,1,1,1
COV,2323,12345.5,1,2,2
COV,2323,12345.5,2,2,1
"""
        sd = SD_VALID.replace(
            "SD_STATE,2323,12345.5,0,0,A,G,G02,G01,0,1,0,2,11",
            "SD_STATE,2323,12345.5,0,0,A,G,G02,G01,0,1,0,2,2",
        ).replace(
            "SD_STATE,2323,12345.5,1,0,A,G,G03,G01,0,2,0,5,19",
            "SD_STATE,2323,12345.5,1,0,A,G,G03,G01,0,2,0,5,2",
        ).replace("0,0,11", "0,0,2").replace("0,1,4.5", "0,1,3").replace("1,1,19", "1,1,2")
        self.write_case(raw, sd)
        report = validate_files(self.raw, self.sd)
        self.assertFalse(report["all_valid"])
        errors = report["epochs"][0]["errors"]
        self.assertTrue(any("raw_not_positive_semidefinite" in error for error in errors))
        self.assertTrue(any("sd_not_positive_semidefinite" in error for error in errors))

    def test_rejects_epoch_set_mismatch(self) -> None:
        self.write_case(RAW_VALID, SD_VALID.replace("12345.5", "12346.0"))
        report = validate_files(self.raw, self.sd)
        self.assertFalse(report["all_valid"])
        self.assertEqual(report["matched_epoch_count"], 0)
        self.assertEqual(len(report["file_errors"]), 2)

    def test_rejects_inconsistent_difference_and_upper_triangle_counts(self) -> None:
        corrupted = SD_VALID.replace(
            "META,2323,12345.5,OK,3,2,1,0,3,0,FILTER_POSTERIOR",
            "META,2323,12345.5,OK,3,1,1,0,1,0,FILTER_POSTERIOR",
        )
        self.write_case(RAW_VALID, corrupted)
        report = validate_files(self.raw, self.sd)
        self.assertFalse(report["all_valid"])
        errors = report["epochs"][0]["errors"]
        self.assertTrue(any("difference_state_count=1,expected=2" in error for error in errors))
        self.assertTrue(any("sd_invalid_upper_triangle_pairs" in error for error in errors))

    def test_cli_writes_atomic_json(self) -> None:
        with redirect_stdout(io.StringIO()):
            status = main([str(self.raw), str(self.sd), str(self.output)])
        self.assertEqual(status, 0)
        payload = json.loads(self.output.read_text(encoding="utf-8"))
        self.assertTrue(payload["all_valid"])
        self.assertFalse(list(self.root.glob("audit.json.*.tmp")))


if __name__ == "__main__":
    unittest.main()
