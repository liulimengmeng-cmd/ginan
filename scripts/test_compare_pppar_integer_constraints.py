from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from compare_pppar_integer_constraints import compare_constraints


def write_trace(text: str) -> Path:
    path = Path(tempfile.mkdtemp(prefix="ginan_pppar_constraint_test_")) / "Network.trace"
    path.write_text(text, encoding="utf-8")
    return path


TRACE_A = """fixAndHoldAmbiguities: 2024-07-17 03:08:30.00
      Applying:   -1 A(SOFI,G12,L1C)  +1 A(SOFI,G12,L2W)  +1 A(SOFI,G25,L1C)  -1 A(SOFI,G25,L2W) =  +24.00000
      Applying:   -1 A(SOFI,G12,L1C)  +1 A(SOFI,G12,L2W)  +1 A(SOFI,G32,L1C)  -1 A(SOFI,G32,L2W) =  +22.00000
"""

TRACE_B_EQUIVALENT = """fixAndHoldAmbiguities: 2024-07-17 03:08:30.00
      Applying:   +1 A(SOFI,G12,L1C)  -1 A(SOFI,G12,L2W)  -1 A(SOFI,G25,L1C)  +1 A(SOFI,G25,L2W) =  -24.00000
      Applying:   -1 A(SOFI,G25,L1C)  +1 A(SOFI,G25,L2W)  +1 A(SOFI,G32,L1C)  -1 A(SOFI,G32,L2W) =   -2.00000
"""

TRACE_B_MISMATCH = TRACE_B_EQUIVALENT.replace("-24.00000", "-23.00000")

TRACE_SINGLE_SIGNAL_A = """fixAndHoldAmbiguities: 2024-07-17 03:31:00.00
      Applying:   -1 A(MARS,G25,L2W)  +1 A(MARS,G29,L2W) =   -6.00000
"""

TRACE_SINGLE_SIGNAL_B = """fixAndHoldAmbiguities: 2024-07-17 03:31:00.00
      Applying:   +1 A(MARS,G25,L2W)  -1 A(MARS,G29,L2W) =   +6.00000
"""


class ConstraintComparisonTests(unittest.TestCase):
    def test_accepts_equivalent_constraints_with_different_pivot(self) -> None:
        report = compare_constraints(write_trace(TRACE_A), write_trace(TRACE_B_EQUIVALENT))
        self.assertTrue(report["integer_differences_consistent_on_shared_pairs"])
        self.assertEqual(report["integer_difference_mismatch_count"], 0)
        self.assertEqual(report["shared_implied_wide_lane_pair_comparisons"], 3)

    def test_detects_restart_integer_mismatch(self) -> None:
        report = compare_constraints(write_trace(TRACE_A), write_trace(TRACE_B_MISMATCH))
        self.assertFalse(report["integer_differences_consistent_on_shared_pairs"])
        self.assertGreater(report["integer_difference_mismatch_count"], 0)

    def test_accepts_equivalent_single_signal_satellite_difference(self) -> None:
        report = compare_constraints(
            write_trace(TRACE_SINGLE_SIGNAL_A), write_trace(TRACE_SINGLE_SIGNAL_B)
        )
        self.assertTrue(report["integer_differences_consistent_on_shared_pairs"])
        self.assertEqual(report["run_a_single_signal_rows"], 1)
        self.assertEqual(report["shared_implied_single_signal_pair_comparisons"], 1)


if __name__ == "__main__":
    unittest.main()
