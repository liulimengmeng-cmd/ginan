#!/usr/bin/env python3

import tempfile
import unittest
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from validate_stec_covariance import validate_file


SCHEMA = """# GINAN_STEC_COVARIANCE_V3
# META,gps_week,gps_tow,status,state_count,upper_triangle_count,max_abs_asymmetry_tecu2,posterior_stage,ar_routine_invoked,ar_eligible_ambiguity_count,ar_integer_ambiguity_coordinate_count,ar_receiver_single_difference_applied,ar_receiver_datum_group_count,ar_dropped_singleton_group_count,ar_resolved_combination_count,ar_pseudoobservations_submitted,ar_mode,ar_configured_success_rate_threshold,ar_configured_solution_ratio_threshold,ar_diagnostic_status,ar_selected_decorrelated_ambiguity_count,ar_integer_candidate_count,ar_bootstrapped_success_rate,ar_best_squared_norm,ar_second_squared_norm,ar_solution_ratio
# STATE,gps_week,gps_tow,local_index,site,satellite,state_number,filter_index,estimate_tecu,variance_tecu2
# COV,gps_week,gps_tow,row_local_index,column_local_index,covariance_tecu2
"""

V2_SCHEMA = """# GINAN_STEC_COVARIANCE_V2
# META,gps_week,gps_tow,status,state_count,upper_triangle_count,max_abs_asymmetry_tecu2,posterior_stage,ar_routine_invoked,ar_eligible_ambiguity_count,ar_resolved_combination_count,ar_pseudoobservations_submitted,ar_mode,ar_configured_success_rate_threshold,ar_configured_solution_ratio_threshold,ar_diagnostic_status,ar_selected_decorrelated_ambiguity_count,ar_integer_candidate_count,ar_bootstrapped_success_rate,ar_best_squared_norm,ar_second_squared_norm,ar_solution_ratio
# STATE,gps_week,gps_tow,local_index,site,satellite,state_number,filter_index,estimate_tecu,variance_tecu2
# COV,gps_week,gps_tow,row_local_index,column_local_index,covariance_tecu2
"""


def write_case(payload: str) -> Path:
    directory = Path(tempfile.mkdtemp(prefix="ginan_stec_covariance_test_"))
    path = directory / "case.STEC.COV"
    path.write_text(SCHEMA + payload, encoding="utf-8")
    return path


class StecCovarianceValidatorTests(unittest.TestCase):
    def test_accepts_complete_positive_definite_epoch(self) -> None:
        path = write_case(
            """META,2323,12345.5,OK,2,3,0,FILTER_POSTERIOR_AFTER_AR_PSEUDOOBS_SUBMITTED_UNVERIFIED,1,5,3,1,2,0,2,1,LAMBDA_ALT,0.9999,3,RESOLVED_RATIO_ACCEPTED,2,2,0.99995,1.25,5,4
STATE,2323,12345.5,0,ALIC,G01,0,7,12.5,4
STATE,2323,12345.5,1,HOB2,G03,0,9,8.25,2.25
COV,2323,12345.5,0,0,4
COV,2323,12345.5,0,1,-0.75
COV,2323,12345.5,1,1,2.25
"""
        )
        report = validate_file(path)
        self.assertTrue(report["all_valid"])
        self.assertEqual(report["valid_epoch_count"], 1)

    def test_rejects_incomplete_or_indefinite_epoch(self) -> None:
        path = write_case(
            """META,2323,12345.5,OK,2,3,0,FILTER_POSTERIOR_NO_EPOCH_AR,0,0,0,0,0,0,0,0,OFF,0.9999,3,NOT_RUN,0,0,-1,-1,-1,-1
STATE,2323,12345.5,0,ALIC,G01,0,7,12.5,1
STATE,2323,12345.5,1,HOB2,G03,0,9,8.25,1
COV,2323,12345.5,0,0,1
COV,2323,12345.5,0,1,2
COV,2323,12345.5,1,1,1
"""
        )
        report = validate_file(path)
        self.assertFalse(report["all_valid"])
        errors = report["epochs"][0]["errors"]
        self.assertTrue(any("not_positive_semidefinite" in error for error in errors))

    def test_rejects_writer_skip_status(self) -> None:
        path = write_case(
            "META,2323,12345.5,STATE_LIMIT_EXCEEDED,600,0,0,FILTER_POSTERIOR_NO_EPOCH_AR,0,0,0,0,0,0,0,0,OFF,0.9999,3,NOT_RUN,0,0,-1,-1,-1,-1\n"
        )
        report = validate_file(path)
        self.assertFalse(report["all_valid"])
        self.assertIn("writer_status=STATE_LIMIT_EXCEEDED", report["epochs"][0]["errors"])

    def test_rejects_zero_state_epoch_without_crashing(self) -> None:
        path = write_case(
            "META,2323,12345.5,NO_STATES,0,0,0,FILTER_POSTERIOR_NO_EPOCH_AR,"
            "0,0,0,0,0,0,0,0,OFF,0.9999,3,NOT_RUN,0,0,-1,-1,-1,-1\n"
        )
        report = validate_file(path)
        self.assertFalse(report["all_valid"])
        errors = report["epochs"][0]["errors"]
        self.assertIn("writer_status=NO_STATES", errors)
        self.assertIn("state_count_not_positive", errors)
        self.assertIn("nonfinite_estimate_or_covariance", errors)

    def test_rejects_inconsistent_ar_evidence(self) -> None:
        path = write_case(
            "META,2323,12345.5,STATE_LIMIT_EXCEEDED,600,0,0,"
            "FILTER_POSTERIOR_NO_EPOCH_AR,0,5,3,1,2,0,2,1,LAMBDA_ALT,0.9999,3,"
            "RESOLVED_RATIO_ACCEPTED,2,2,0.99995,1.25,5,4\n"
        )
        report = validate_file(path)
        self.assertFalse(report["all_valid"])
        self.assertIn(
            "ar_evidence_present_without_routine_invocation",
            report["epochs"][0]["errors"],
        )

    def test_reports_malformed_short_record_without_index_error(self) -> None:
        path = write_case(
            "META,2323,12345.5,STATE_LIMIT_EXCEEDED,600,0,0,"
            "FILTER_POSTERIOR_NO_EPOCH_AR,0,0,0,0,0,0,0,0,OFF,0.9999,3,"
            "NOT_RUN,0,0,-1,-1,-1,-1\nSTATE\n"
        )
        with self.assertRaisesRegex(ValueError, "STATE expects 10 fields"):
            validate_file(path)

    def test_accepts_legacy_v2_metadata(self) -> None:
        directory = Path(tempfile.mkdtemp(prefix="ginan_stec_covariance_v2_test_"))
        path = directory / "legacy.STEC.COV"
        path.write_text(
            V2_SCHEMA
            + "META,2323,12345.5,OK,1,1,0,FILTER_POSTERIOR_NO_EPOCH_AR,0,0,0,0,"
            "OFF,0.9999,3,NOT_RUN,0,0,-1,-1,-1,-1\n"
            "STATE,2323,12345.5,0,ALIC,G01,0,7,12.5,4\n"
            "COV,2323,12345.5,0,0,4\n",
            encoding="utf-8",
        )
        report = validate_file(path)
        self.assertTrue(report["all_valid"])
        self.assertEqual(report["schema"], "GINAN_STEC_COVARIANCE_V2")


if __name__ == "__main__":
    unittest.main()
