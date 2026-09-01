from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from audit_pppar_dual_frequency_datum import audit_trace


class DualFrequencyDatumAuditTest(unittest.TestCase):
    def test_reads_rts_state_time_without_epoch_marker(self) -> None:
        payload = """
+STATES/PPP_RTS
*	-1	2024-07-17 12:00:00.00	AMBIGUITY	G02	A	L1C	1	0.1	0
-STATES/PPP_RTS
PPP_AR DUAL_FREQUENCY_BASIS original=8 rows=6 expected_rank=6 actual_rank=6 complete_groups=1 incomplete_groups=0 paired_ambiguities=8 unmatched_ambiguities=0 integer_valued=1 full_row_rank=1 covers_all=1 status=FULL_DUAL_FREQUENCY_INTEGER_BASIS action=PROBE_ONLY_NOT_SUBMITTED
PPP_AR DUAL_FREQUENCY_DATUM_SUMMARY full_candidate_groups=0 target_groups=1 full_candidate_rows=0 target_rank=6 status=NO_FULL_INTEGER_DATUM_CANDIDATE wrong_fix_certified=0 filter_feedback=0 action=PROBE_ONLY_NOT_SUBMITTED
"""
        with tempfile.TemporaryDirectory() as temporary:
            trace = Path(temporary) / "Network_smoothed.trace"
            trace.write_text(payload, encoding="utf-8")
            report = audit_trace(trace)
        self.assertEqual(report["audited_epoch_count"], 1)
        self.assertEqual(report["structural_pass_epoch_count"], 1)

    def test_separates_structural_basis_from_integer_candidate(self) -> None:
        payload = """
------=============== Epoch 1 =============-----------
PPP_AR DUAL_FREQUENCY_BASIS original=8 rows=6 expected_rank=6 actual_rank=6 complete_groups=1 incomplete_groups=0 paired_ambiguities=8 unmatched_ambiguities=0 integer_valued=1 full_row_rank=1 covers_all=1 status=FULL_DUAL_FREQUENCY_INTEGER_BASIS action=PROBE_ONLY_NOT_SUBMITTED
PPP_AR DUAL_FREQUENCY_GROUP receiver=A system=GPS first_signal=L1C second_signal=L2W reference=G01 common_satellites=4 unmatched_ambiguities=0 wide_lane_rows=3 second_d2_rows=3 complete_graph=1 action=PROBE_ONLY_NOT_SUBMITTED
PPP_AR DUAL_FREQUENCY_STAGE receiver=A system=GPS reference=G01 wide_lane_target=3 wide_lane_fixed=0 wide_lane_status=SUCCESS_RATE_BELOW_THRESHOLD wide_lane_success_rate=0.8 wide_lane_ratio=-1 second_d2_target=3 second_d2_fixed=0 second_d2_status=NOT_RUN combined_rank=0 target_rank=6 status=WIDE_LANE_FAMILY_NOT_FULL action=PROBE_ONLY_NOT_SUBMITTED
PPP_AR DUAL_FREQUENCY_DATUM_SUMMARY full_candidate_groups=0 target_groups=1 full_candidate_rows=0 target_rank=6 status=NO_FULL_INTEGER_DATUM_CANDIDATE wrong_fix_certified=0 filter_feedback=0 action=PROBE_ONLY_NOT_SUBMITTED
"""
        with tempfile.TemporaryDirectory() as temporary:
            trace = Path(temporary) / "Network.trace"
            trace.write_text(payload, encoding="utf-8")
            report = audit_trace(trace)
        self.assertTrue(report["structural_basis_pass"])
        self.assertEqual(report["structural_pass_epoch_count"], 1)
        self.assertEqual(report["full_epoch_candidate_count"], 0)
        self.assertEqual(report["second_family_reached_group_count"], 0)
        self.assertEqual(report["wide_lane_full_group_count"], 0)
        self.assertEqual(report["second_family_fixed_row_count"], 0)
        self.assertEqual(report["structural_failures"], [])

    def test_rejects_rank_or_submission_violation(self) -> None:
        payload = """
------=============== Epoch 1 =============-----------
PPP_AR DUAL_FREQUENCY_BASIS original=8 rows=6 expected_rank=6 actual_rank=5 complete_groups=1 incomplete_groups=0 paired_ambiguities=8 unmatched_ambiguities=0 integer_valued=1 full_row_rank=0 covers_all=1 status=INVALID_DUAL_FREQUENCY_INTEGER_BASIS action=SUBMITTED
PPP_AR DUAL_FREQUENCY_DATUM_SUMMARY full_candidate_groups=0 target_groups=1 full_candidate_rows=0 target_rank=6 status=NO_FULL_INTEGER_DATUM_CANDIDATE wrong_fix_certified=0 filter_feedback=0 action=SUBMITTED
"""
        with tempfile.TemporaryDirectory() as temporary:
            trace = Path(temporary) / "Network.trace"
            trace.write_text(payload, encoding="utf-8")
            report = audit_trace(trace)
        self.assertFalse(report["structural_basis_pass"])
        self.assertEqual(len(report["action_violations"]), 2)
        self.assertEqual(len(report["structural_failures"]), 1)

    def test_counts_mapped_second_family_rows(self) -> None:
        payload = """
------=============== Epoch 1 =============-----------
PPP_AR DUAL_FREQUENCY_BASIS original=8 rows=6 expected_rank=6 actual_rank=6 complete_groups=1 incomplete_groups=0 paired_ambiguities=8 unmatched_ambiguities=0 integer_valued=1 full_row_rank=1 covers_all=1 status=FULL_DUAL_FREQUENCY_INTEGER_BASIS action=PROBE_ONLY_NOT_SUBMITTED
PPP_AR DUAL_FREQUENCY_CANDIDATE_ROW receiver=A system=GPS reference=G01 candidate_scope=SECOND_D2_PARTIAL family=SECOND_D2 row=0 rhs=4 support=2 status=INTEGER_MAPPED terms=example action=PROBE_ONLY_NOT_SUBMITTED
PPP_AR DUAL_FREQUENCY_DATUM_SUMMARY full_candidate_groups=0 target_groups=1 full_candidate_rows=0 target_rank=6 status=NO_FULL_INTEGER_DATUM_CANDIDATE wrong_fix_certified=0 filter_feedback=0 action=PROBE_ONLY_NOT_SUBMITTED
"""
        with tempfile.TemporaryDirectory() as temporary:
            trace = Path(temporary) / "Network.trace"
            trace.write_text(payload, encoding="utf-8")
            report = audit_trace(trace)
        self.assertEqual(report["candidate_row_count"], 1)
        self.assertEqual(report["candidate_family_counts"], {"SECOND_D2": 1})
        self.assertEqual(
            report["candidate_scope_counts"],
            {"SECOND_D2_PARTIAL": 1},
        )
        self.assertEqual(report["invalid_candidate_row_count"], 0)

    def test_distinguishes_full_selected_subset_from_full_visible_group(self) -> None:
        payload = """
------=============== Epoch 1 =============-----------
PPP_AR DUAL_FREQUENCY_BASIS original=12 rows=10 expected_rank=10 actual_rank=10 complete_groups=1 incomplete_groups=0 paired_ambiguities=12 unmatched_ambiguities=0 integer_valued=1 full_row_rank=1 covers_all=1 status=FULL_DUAL_FREQUENCY_INTEGER_BASIS action=PROBE_ONLY_NOT_SUBMITTED
PPP_AR DUAL_FREQUENCY_STAGE receiver=A system=GPS reference=G01 visible_family_count=5 selected_family_count=3 excluded_satellites=2 wide_lane_target=3 wide_lane_fixed=3 wide_lane_status=RESOLVED_RATIO_TEST second_d2_target=3 second_d2_fixed=3 second_d2_status=FULL_INTEGER_FAMILY_RESOLVED combined_rank=6 target_rank=6 status=FULL_SELECTED_SUBSET_INTEGER_DATUM_CANDIDATE_UNVERIFIED action=PROBE_ONLY_NOT_SUBMITTED
PPP_AR DUAL_FREQUENCY_CANDIDATE_ROW receiver=A system=GPS reference=G01 candidate_scope=SELECTED_SUBSET family=WIDE_LANE row=0 rhs=1 support=4 status=INTEGER_MAPPED terms=example action=PROBE_ONLY_NOT_SUBMITTED
PPP_AR DUAL_FREQUENCY_DATUM_SUMMARY full_candidate_groups=0 target_groups=1 full_candidate_rows=0 target_rank=10 selected_candidate_groups=1 selected_candidate_rows=6 status=FULL_SELECTED_SUBSET_INTEGER_DATUM_CANDIDATE_UNVERIFIED wrong_fix_certified=0 filter_feedback=0 action=PROBE_ONLY_NOT_SUBMITTED
PPP_AR DUAL_FREQUENCY_CONTROL legacy_feedback=0 new_subset_feedback=0 status=FLOAT_STATE_PROBE_ONLY action=PROBE_ONLY_NOT_SUBMITTED
"""
        with tempfile.TemporaryDirectory() as temporary:
            trace = Path(temporary) / "Network.trace"
            trace.write_text(payload, encoding="utf-8")
            report = audit_trace(trace)
        self.assertTrue(report["structural_basis_pass"])
        self.assertEqual(report["selected_subset_group_candidate_count"], 1)
        self.assertEqual(report["full_selected_epoch_candidate_count"], 1)
        self.assertEqual(report["full_visible_epoch_candidate_count"], 0)
        self.assertEqual(report["full_epoch_candidate_count"], 1)
        self.assertEqual(report["candidate_scope_counts"], {"SELECTED_SUBSET": 1})
        self.assertTrue(report["probe_isolated_from_filter_feedback"])
        self.assertEqual(report["pseudoobs_submission_count"], 0)
        self.assertEqual(
            report["selected_subset_by_receiver"]["A"],
            {
                "candidate_epoch_count": 1,
                "first_candidate_epoch": 1,
                "last_candidate_epoch": 1,
                "longest_consecutive_epoch_run": 1,
                "selected_family_count_min": 3,
                "selected_family_count_max": 3,
                "excluded_satellite_count_min": 2,
                "excluded_satellite_count_max": 2,
            },
        )

    def test_detects_legacy_filter_feedback_in_same_run(self) -> None:
        payload = """
------=============== Epoch 1 =============-----------
PPP_AR DUAL_FREQUENCY_BASIS original=8 rows=6 expected_rank=6 actual_rank=6 complete_groups=1 incomplete_groups=0 paired_ambiguities=8 unmatched_ambiguities=0 integer_valued=1 full_row_rank=1 covers_all=1 status=FULL_DUAL_FREQUENCY_INTEGER_BASIS action=PROBE_ONLY_NOT_SUBMITTED
PPP_AR DUAL_FREQUENCY_DATUM_SUMMARY full_candidate_groups=0 target_groups=1 full_candidate_rows=0 target_rank=6 status=NO_FULL_INTEGER_DATUM_CANDIDATE wrong_fix_certified=0 filter_feedback=0 action=PROBE_ONLY_NOT_SUBMITTED
PPP_AR PSEUDOOBS_SUBMISSION rows=2 status=FILTER_CALL_RETURNED_SUBMITTED_UNVERIFIED
"""
        with tempfile.TemporaryDirectory() as temporary:
            trace = Path(temporary) / "Network.trace"
            trace.write_text(payload, encoding="utf-8")
            report = audit_trace(trace)
        self.assertFalse(report["probe_isolated_from_filter_feedback"])
        self.assertFalse(report["safety"]["diagnostic_only"])
        self.assertEqual(report["pseudoobs_submission_count"], 1)
        self.assertEqual(report["pseudoobs_submitted_row_count"], 2)

    def test_tracks_guarded_subset_feedback_without_calling_it_certified(self) -> None:
        payload = """
------=============== Epoch 1 =============-----------
PPP_AR DUAL_FREQUENCY_BASIS original=8 rows=6 expected_rank=6 actual_rank=6 complete_groups=1 incomplete_groups=0 paired_ambiguities=8 unmatched_ambiguities=0 integer_valued=1 full_row_rank=1 covers_all=1 status=FULL_DUAL_FREQUENCY_INTEGER_BASIS action=CANDIDATE_PENDING_SAFETY_GATES
PPP_AR DUAL_FREQUENCY_DATUM_SUMMARY full_candidate_groups=1 target_groups=1 full_candidate_rows=6 target_rank=6 selected_candidate_groups=1 selected_candidate_rows=6 status=FULL_VISIBLE_INTEGER_DATUM_CANDIDATE_UNVERIFIED feedback_requested=1 wrong_fix_certified=0 filter_feedback=0 action=CANDIDATE_PENDING_SAFETY_GATES
PPP_AR DUAL_FREQUENCY_PHASE_BIAS_GATE used_ambiguities=8 phase_bias_model_disabled=0 missing_phase_bias=0 invalid_map=0 status=COMPLETE_PRODUCT_COVERAGE
PPP_AR DUAL_FREQUENCY_FEEDBACK_MAP rows=6 ambiguities=8 status=CANONICAL_TO_FILTER_GAUGE
PPP_AR PSEUDOOBS_SUBMISSION rows=6 status=FILTER_CALL_RETURNED_SUBMITTED_UNVERIFIED
PPP_AR DUAL_FREQUENCY_CONTROL candidate_rows=6 legacy_feedback=0 new_subset_feedback=1 status=SUBMITTED_UNVERIFIED action=FILTER_CALL_RETURNED
"""
        with tempfile.TemporaryDirectory() as temporary:
            trace = Path(temporary) / "Network.trace"
            trace.write_text(payload, encoding="utf-8")
            report = audit_trace(trace)
        self.assertTrue(report["structural_basis_pass"])
        self.assertTrue(report["feedback_path_exercised"])
        self.assertTrue(report["feedback_submission_consistent"])
        self.assertTrue(report["feedback_safety_protocol_pass"])
        self.assertEqual(report["feedback_submitted_epoch_count"], 1)
        self.assertFalse(report["safety"]["filter_feedback_certified"])
        self.assertTrue(
            report["safety"]["filter_call_returned_submitted_unverified"]
        )

    def test_forward_state_rows_do_not_replace_epoch_number(self) -> None:
        payload = """
------=============== Epoch 7 =============-----------
*\t-1\t2024-07-17 00:03:00.00\tREC_POS\tA\tX\t1\t0.1\t0
PPP_AR DUAL_FREQUENCY_BASIS original=8 rows=6 expected_rank=6 actual_rank=6 complete_groups=1 incomplete_groups=0 paired_ambiguities=8 unmatched_ambiguities=0 integer_valued=1 full_row_rank=1 covers_all=1 status=FULL_DUAL_FREQUENCY_INTEGER_BASIS action=PROBE_ONLY_NOT_SUBMITTED
PPP_AR DUAL_FREQUENCY_DATUM_SUMMARY full_candidate_groups=0 target_groups=1 full_candidate_rows=0 target_rank=6 status=NO_FULL_INTEGER_DATUM_CANDIDATE wrong_fix_certified=0 filter_feedback=0 action=PROBE_ONLY_NOT_SUBMITTED
"""
        with tempfile.TemporaryDirectory() as temporary:
            trace = Path(temporary) / "Network.trace"
            trace.write_text(payload, encoding="utf-8")
            report = audit_trace(trace)
        self.assertEqual(report["structural_pass_epochs"], [7])

    def test_certifies_phase_bias_disabled_negative_control_protocol(self) -> None:
        payload = """
------=============== Epoch 1 =============-----------
PPP_AR DUAL_FREQUENCY_BASIS original=8 rows=6 expected_rank=6 actual_rank=6 complete_groups=1 incomplete_groups=0 paired_ambiguities=8 unmatched_ambiguities=0 integer_valued=1 full_row_rank=1 covers_all=1 status=FULL_DUAL_FREQUENCY_INTEGER_BASIS action=CANDIDATE_PENDING_SAFETY_GATES
PPP_AR DUAL_FREQUENCY_DATUM_SUMMARY full_candidate_groups=1 target_groups=1 full_candidate_rows=6 target_rank=6 selected_candidate_groups=1 selected_candidate_rows=6 status=FULL_VISIBLE_INTEGER_DATUM_CANDIDATE_UNVERIFIED feedback_requested=1 wrong_fix_certified=0 filter_feedback=0 action=CANDIDATE_PENDING_SAFETY_GATES
PPP_AR DUAL_FREQUENCY_PHASE_BIAS_GATE used_ambiguities=8 phase_bias_model_disabled=8 missing_phase_bias=0 invalid_map=0 status=REJECTED_INCOMPLETE_PRODUCT_COVERAGE
PPP_AR DUAL_FREQUENCY_CONTROL candidate_rows=6 legacy_feedback=0 new_subset_feedback=0 status=REJECTED_PHASE_BIAS_GATE action=NOT_SUBMITTED
"""
        with tempfile.TemporaryDirectory() as temporary:
            trace = Path(temporary) / "Network.trace"
            trace.write_text(payload, encoding="utf-8")
            report = audit_trace(trace)
        self.assertTrue(report["feedback_negative_control_pass"])
        self.assertEqual(report["candidate_epoch_count"], 1)
        self.assertEqual(report["rejected_product_gate_epoch_count"], 1)
        self.assertEqual(report["pseudoobs_submission_count"], 0)


if __name__ == "__main__":
    unittest.main()
