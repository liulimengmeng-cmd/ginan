from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from audit_pppar_candidate_restart_consistency import (
    compare_groups,
    parse_candidate_groups,
)


def trace(epoch: int, rhs: int) -> str:
    return f"""
------=============== Epoch {epoch} =============-----------
PPP_AR DUAL_FREQUENCY_CANDIDATE_ROW receiver=A system=GPS reference=G01 candidate_scope=SELECTED_SUBSET family=SECOND_D2 row=3 rhs={rhs} float_value={rhs + 0.1} float_minus_integer=0.1 formal_sigma=0.05 support=2 status=INTEGER_MAPPED terms=+1 A(A,G02,L2W) -1 A(A,G01,L2W) action=PROBE_ONLY_NOT_SUBMITTED
"""


class CandidateRestartConsistencyTest(unittest.TestCase):
    def test_aligns_epoch_offset_and_detects_agreement(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            primary_path = Path(temporary) / "primary.trace"
            restart_path = Path(temporary) / "restart.trace"
            primary_path.write_text(trace(1441, 7), encoding="utf-8")
            restart_path.write_text(trace(1, 7), encoding="utf-8")
            report = compare_groups(
                parse_candidate_groups(primary_path),
                parse_candidate_groups(restart_path, 1440),
            )
        self.assertEqual(report["shared_candidate_group_count"], 1)
        self.assertEqual(report["exact_candidate_group_count"], 1)
        self.assertEqual(report["agreed_integer_rhs_count"], 1)
        self.assertTrue(report["all_shared_integer_rows_agree"])

    def test_detects_wrong_integer(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            primary_path = Path(temporary) / "primary.trace"
            restart_path = Path(temporary) / "restart.trace"
            primary_path.write_text(trace(1, 7), encoding="utf-8")
            restart_path.write_text(trace(1, 8), encoding="utf-8")
            report = compare_groups(
                parse_candidate_groups(primary_path),
                parse_candidate_groups(restart_path),
            )
        self.assertEqual(report["disagreed_integer_rhs_count"], 1)
        self.assertEqual(report["disagreement_epoch_count"], 1)
        self.assertEqual(report["disagreement_samples"][0]["difference"], 1)
        self.assertEqual(
            report["primary_disagreement_abs_residual_median"],
            0.1,
        )
        self.assertFalse(report["all_shared_integer_rows_agree"])


if __name__ == "__main__":
    unittest.main()
