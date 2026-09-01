from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from audit_pppar_feedback_information import audit_trace


class FeedbackInformationAuditTest(unittest.TestCase):
    def test_incidence_rank_gates_and_information_blocks(self) -> None:
        payload = """
------=============== Epoch 1 =============-----------
PPP_AR DUAL_FREQUENCY_GROUP receiver=MATE system=GPS first_signal=L1C second_signal=L2W reference=G02 common_satellites=4 unmatched_ambiguities=0 wide_lane_rows=3 second_d2_rows=3 complete_graph=1 action=CANDIDATE_PENDING_SAFETY_GATES
PPP_AR DUAL_FREQUENCY_STAGE receiver=MATE system=GPS reference=G02 visible_family_count=3 selected_family_count=2 excluded_satellites=1 wide_lane_target=2 wide_lane_fixed=2 wide_lane_status=RESOLVED_RATIO_ACCEPTED wide_lane_success_rate=0.99999 wide_lane_ratio=4 second_d2_target=2 second_d2_fixed=2 second_d2_status=FULL_INTEGER_FAMILY_RESOLVED second_d2_minimum_decorrelated=1 second_d2_attempted_stages=1 second_d2_accepted_stages=1 second_d2_last_attempt_status=RESOLVED_RATIO_ACCEPTED second_d2_success_rate=0.99995 second_d2_ratio=5 second_d2_integer_valued=1 second_d2_row_rank=2 second_d2_unimodular=1 combined_integer_valued=1 combined_rank=4 target_rank=4 status=FULL_SELECTED_SUBSET_INTEGER_DATUM_CANDIDATE_UNVERIFIED action=CANDIDATE_PENDING_SAFETY_GATES
PPP_AR FEEDBACK_SHADOW_SUMMARY rows=4 innovation_norm=0.2 joint_nis=2 nis_per_row=0.5 nis_gate_configured=0 status=VALID_LINEAR_SHADOW action=ANALYTIC_ONE_STEP_DIAGNOSTIC_ONLY
PPP_AR FEEDBACK_SHADOW_BLOCK block=STEC states=4 prior_trace=10 predicted_trace_reduction=2 relative_trace_gain=0.2 predicted_shift_norm=0.1 coupled_states=2 coupled_prior_trace=4 coupled_predicted_trace_reduction=2 coupled_relative_trace_gain=0.5 status=VALID_LINEAR_SHADOW action=DIAGNOSTIC_ONLY
PPP_AR FEEDBACK_REALISED_BLOCK block=STEC states=4 realised_trace_reduction=2 relative_trace_gain=0.2 realised_shift_norm=0.1 coupled_states=2 coupled_realised_trace_reduction=2 coupled_relative_trace_gain=0.5 status=FILTER_CALL_RETURNED action=DIAGNOSTIC_ONLY
PPP_AR FEEDBACK_SHADOW_REALISATION state_update_error_norm=0 covariance_error_norm=1e-12 status=COMPARED_WITH_FILTER_RESULT action=DIAGNOSTIC_ONLY
PPP_AR PSEUDOOBS_SUBMISSION rows=4 status=FILTER_CALL_RETURNED_SUBMITTED_UNVERIFIED
PPP_AR DUAL_FREQUENCY_CONTROL candidate_rows=4 legacy_feedback=0 new_subset_feedback=1 status=SUBMITTED_UNVERIFIED action=FILTER_CALL_RETURNED
------=============== Epoch 2 =============-----------
PPP_AR DUAL_FREQUENCY_GROUP receiver=MATE system=GPS first_signal=L1C second_signal=L2W reference=G02 common_satellites=4 unmatched_ambiguities=0 wide_lane_rows=3 second_d2_rows=3 complete_graph=1 action=CANDIDATE_PENDING_SAFETY_GATES
PPP_AR DUAL_FREQUENCY_STAGE receiver=MATE system=GPS reference=G02 visible_family_count=3 selected_family_count=3 excluded_satellites=0 wide_lane_target=3 wide_lane_fixed=0 wide_lane_status=SUCCESS_RATE_BELOW_THRESHOLD wide_lane_success_rate=0.9 wide_lane_ratio=-1 second_d2_target=3 second_d2_fixed=0 second_d2_status=NOT_RUN second_d2_minimum_decorrelated=1 second_d2_attempted_stages=0 second_d2_accepted_stages=0 second_d2_last_attempt_status=NOT_RUN second_d2_success_rate=-1 second_d2_ratio=-1 second_d2_integer_valued=1 second_d2_row_rank=0 second_d2_unimodular=0 combined_integer_valued=0 combined_rank=0 target_rank=6 status=WIDE_LANE_FAMILY_NOT_FULL action=CANDIDATE_PENDING_SAFETY_GATES
PPP_AR DUAL_FREQUENCY_CONTROL candidate_rows=0 legacy_feedback=0 new_subset_feedback=0 status=REJECTED_INVALID_CANDIDATE_BLOCK action=NOT_SUBMITTED
"""
        with tempfile.TemporaryDirectory() as temporary:
            trace = Path(temporary) / "Network.trace"
            trace.write_text(payload, encoding="utf-8")
            report = audit_trace(trace)
        self.assertEqual(report["epoch_count"], 2)
        self.assertEqual(report["feedback_epoch_count"], 1)
        self.assertAlmostEqual(report["feedback_epoch_incidence"], 0.5)
        self.assertAlmostEqual(report["station_epoch_feedback_incidence"], 0.5)
        self.assertEqual(report["eligible_rank_epoch_sum"], 12)
        self.assertEqual(report["submitted_rank_epoch_sum"], 4)
        self.assertAlmostEqual(report["submitted_rank_coverage"], 1 / 3)
        self.assertEqual(report["valid_shadow_record_count"], 1)
        self.assertEqual(
            report["first_failure_gate_counts"]["WIDE_LANE_INCOMPLETE_SUCCESS_RATE_BELOW_THRESHOLD"],
            1,
        )
        self.assertAlmostEqual(
            report["state_block_information"]["STEC"]["realised_relative_trace_gain"]["median"],
            0.2,
        )
        self.assertAlmostEqual(
            report["state_block_information"]["STEC"]["coupled_realised_relative_trace_gain"]["median"],
            0.5,
        )

        with tempfile.TemporaryDirectory() as temporary:
            trace = Path(temporary) / "Network.trace"
            trace.write_text(payload, encoding="utf-8")
            sliced = audit_trace(trace, epoch_start=2, epoch_end=2)
        self.assertEqual(sliced["epoch_count"], 1)
        self.assertEqual(sliced["feedback_epoch_count"], 0)
        self.assertEqual(
            sliced["epoch_window"],
            {
                "requested_start": 2,
                "requested_end": 2,
                "observed_start": 2,
                "observed_end": 2,
            },
        )
        self.assertNotIn("STEC", sliced["state_block_information"])

    def test_sustained_epoch_is_reported_in_absolute_epoch_numbers(self) -> None:
        lines = []
        for epoch in range(1001, 1360):
            lines.extend(
                [
                    f"------=============== Epoch {epoch} =============-----------",
                    "PPP_AR DUAL_FREQUENCY_GROUP receiver=MATE system=GPS reference=G02 wide_lane_rows=1 complete_graph=1",
                    "PPP_AR DUAL_FREQUENCY_STAGE receiver=MATE system=GPS reference=G02 wide_lane_target=1 wide_lane_fixed=1 second_d2_target=1 second_d2_fixed=1 combined_rank=2 status=FULL_SELECTED_SUBSET_INTEGER_DATUM_CANDIDATE_UNVERIFIED",
                    "PPP_AR PSEUDOOBS_SUBMISSION rows=2 status=FILTER_CALL_RETURNED_SUBMITTED_UNVERIFIED",
                ]
            )
        with tempfile.TemporaryDirectory() as temporary:
            trace = Path(temporary) / "Network.trace"
            trace.write_text("\n".join(lines), encoding="utf-8")
            report = audit_trace(trace, epoch_start=1001, epoch_end=1359)
        self.assertEqual(
            report["receiver_summary"]["MATE"]["first_sustained_feedback_epoch"],
            1120,
        )


if __name__ == "__main__":
    unittest.main()
