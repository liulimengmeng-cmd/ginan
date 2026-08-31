from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from audit_pppar_integer_complement import audit_trace


def write_trace(text: str) -> Path:
    path = Path(tempfile.mkdtemp(prefix="ginan_pppar_complement_test_")) / "Network.trace"
    path.write_text(text, encoding="utf-8")
    return path


TRACE = """fixAndHoldAmbiguities: 2024-07-17 00:00:00.00
PPP_AR INTEGER_COMPLEMENT_DIAGNOSTIC stage1_rows=1 remaining_coordinates=3 stage2_probe_rows=1 combined_independent_rows=2 target_integer_rank=4 stage2_status=RESOLVED_RATIO_ACCEPTED stage2_success_rate=0.99995 stage2_selected_decorrelated=3 stage2_integer_candidates=2 stage2_best_squared_norm=0.1 stage2_second_squared_norm=0.5 stage2_ratio=5 action=PROBE_ONLY_NOT_SUBMITTED
PPP_AR INTEGER_COMPLEMENT_ROW stage=2 row=0 rhs=7 support=2 status=INTEGER_MAPPED terms=+1 A(SOFI,G02,L2W) -1 A(SOFI,G01,L2W) action=PROBE_ONLY_NOT_SUBMITTED
      Applying:  +1 A(SOFI,G02,L1C) -1 A(SOFI,G01,L1C) = +4.00000
fixAndHoldAmbiguities: 2024-07-17 00:00:30.00
PPP_AR INTEGER_COMPLEMENT_DIAGNOSTIC stage1_rows=1 remaining_coordinates=3 stage2_probe_rows=0 combined_independent_rows=1 target_integer_rank=4 stage2_status=INSUFFICIENT_DECORRELATED_AMBIGUITIES stage2_success_rate=0.99996 stage2_selected_decorrelated=2 stage2_integer_candidates=0 stage2_best_squared_norm=-1 stage2_second_squared_norm=-1 stage2_ratio=-1 action=PROBE_ONLY_NOT_SUBMITTED
      Applying:  +1 A(SOFI,G02,L1C) -1 A(SOFI,G01,L1C) = +4.00000
"""


class IntegerComplementAuditTests(unittest.TestCase):
    def test_exact_rank_and_probe_only_protocol(self) -> None:
        report = audit_trace(write_trace(TRACE))
        self.assertEqual(report["diagnostic_epoch_count"], 2)
        self.assertEqual(report["accepted_probe_epoch_count"], 1)
        self.assertEqual(report["insufficient_dimension_epoch_count"], 1)
        self.assertEqual(report["stage2_probe_row_count"], 1)
        self.assertEqual(report["stage2_single_signal_row_count"], 1)
        self.assertTrue(report["all_actions_probe_only_not_submitted"])
        self.assertEqual(report["epochs"][0]["combined_exact_rank"], 2)
        self.assertEqual(report["epochs"][1]["stage2_selected_decorrelated"], 2)

    def test_rejects_false_independence_claim(self) -> None:
        dependent = TRACE.replace(
            "+1 A(SOFI,G02,L2W) -1 A(SOFI,G01,L2W)",
            "+1 A(SOFI,G02,L1C) -1 A(SOFI,G01,L1C)",
            1,
        )
        with self.assertRaisesRegex(ValueError, "not independent"):
            audit_trace(write_trace(dependent))


if __name__ == "__main__":
    unittest.main()
