from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from audit_pppar_dual_frequency_datum import audit_trace


class DualFrequencyDatumAuditTest(unittest.TestCase):
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


if __name__ == "__main__":
    unittest.main()
