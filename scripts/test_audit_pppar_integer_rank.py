from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from audit_pppar_integer_rank import audit_trace, exact_rank


def write_trace(text: str) -> Path:
    path = Path(tempfile.mkdtemp(prefix="ginan_pppar_rank_test_")) / "Network.trace"
    path.write_text(text, encoding="utf-8")
    return path


TRACE = """fixAndHoldAmbiguities: 2024-07-17 01:39:30.00
PPP_AR INTEGER_COORDINATES original=10 integer=8 receiver_sd_groups=4 identity_groups=0 dropped_singletons=0
PPP_AR RECEIVER_SD receiver=SOFI system=GPS signal=L1C reference=G12 members=3 status=INTEGER_COORDINATES_CREATED
PPP_AR RECEIVER_SD receiver=SOFI system=GPS signal=L2W reference=G12 members=3 status=INTEGER_COORDINATES_CREATED
      Applying:   -1 A(SOFI,G12,L1C) +1 A(SOFI,G12,L2W) +1 A(SOFI,G25,L1C) -1 A(SOFI,G25,L2W) = +24.00000
      Applying:   -1 A(SOFI,G12,L1C) +1 A(SOFI,G12,L2W) +1 A(SOFI,G32,L1C) -1 A(SOFI,G32,L2W) = +22.00000
PPP_AR PSEUDOOBS_SUBMISSION rows=2 status=FILTER_CALL_RETURNED_SUBMITTED_UNVERIFIED
"""


class IntegerRankAuditTests(unittest.TestCase):
    def test_exact_rank_detects_dependency(self) -> None:
        first = {("A", "G01", "L1C"): 1, ("A", "G02", "L1C"): -1}
        second = {key: 2 * value for key, value in first.items()}
        self.assertEqual(exact_rank([first, second]), 1)

    def test_reports_wide_lane_rank_and_full_target_deficit(self) -> None:
        report = audit_trace(write_trace(TRACE))
        self.assertEqual(report["submitted_epoch_count"], 1)
        self.assertEqual(report["wide_lane_row_count"], 2)
        self.assertEqual(report["single_signal_row_count"], 0)
        self.assertEqual(report["complete_rank_epoch_count"], 0)
        epoch = report["epochs"][0]
        self.assertEqual(epoch["submitted_rank"], 2)
        self.assertEqual(epoch["wide_lane_rank"], 2)
        self.assertEqual(epoch["rank_deficit"], 6)
        self.assertEqual(epoch["receiver_components"]["SOFI"]["target_rank"], 4)


if __name__ == "__main__":
    unittest.main()
