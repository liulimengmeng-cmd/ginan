from __future__ import annotations

import tempfile
import unittest
from datetime import datetime, timezone
from pathlib import Path

from audit_pppar_restart_slip_alignment import audit


def candidate(epoch: int, rhs: int) -> str:
    return f"""
------=============== Epoch {epoch} =============-----------
PPP_AR DUAL_FREQUENCY_CANDIDATE_ROW receiver=A system=GPS reference=G01 candidate_scope=SELECTED_SUBSET family=SECOND_D2 row=3 rhs={rhs} float_value={rhs + 0.1} float_minus_integer=0.1 formal_sigma=0.05 support=2 status=INTEGER_MAPPED terms=+1 A(A,G02,L2W) -1 A(A,G01,L2W) action=PROBE_ONLY_NOT_SUBMITTED
"""


class RestartSlipAlignmentTest(unittest.TestCase):
    def test_reports_different_arc_start_without_common_reset(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            primary = Path(temporary) / "primary.trace"
            restart = Path(temporary) / "restart.trace"
            primary.write_text(candidate(11, 7), encoding="utf-8")
            restart.write_text(candidate(1, 8), encoding="utf-8")
            report = audit(
                primary,
                restart,
                datetime(2024, 7, 17, tzinfo=timezone.utc),
                10,
            )
        self.assertEqual(report["disagreed_integer_rhs_count"], 1)
        self.assertEqual(
            report["reset_alignment_counts"],
            {"one_or_more_satellite_reset_epochs_differ": 1},
        )

    def test_reports_matching_common_slip_epoch(self) -> None:
        event = (
            "PDE-CS-DIAG detector=GF action=detected "
            "epoch=2024-07-17 00:05:00.00 rec=A sat=G02 f1=F1 f2=F2 f3=NONE "
            "value1=1 value2=0 threshold=0.05 reason=phase_jump\n"
        )
        with tempfile.TemporaryDirectory() as temporary:
            primary = Path(temporary) / "primary.trace"
            restart = Path(temporary) / "restart.trace"
            primary.write_text(event + candidate(12, 7), encoding="utf-8")
            restart.write_text(event + candidate(12, 8), encoding="utf-8")
            report = audit(
                primary,
                restart,
                datetime(2024, 7, 17, tzinfo=timezone.utc),
                0,
            )
        self.assertEqual(
            report["reset_alignment_counts"],
            {"all_satellite_reset_epochs_match": 1},
        )


if __name__ == "__main__":
    unittest.main()

