from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from audit_external_bias_application import audit_trace


class ExternalBiasApplicationAuditTest(unittest.TestCase):
    def test_counts_conditioning_and_duplicate_keys(self) -> None:
        phase = (
            "PPP_EXTERNAL_BIAS_APPLICATION type=PHASE time=2024-07-17 00:00:00.000 "
            "receiver=GRAZ satellite=G02 signal=L1C bias_found=1 bias_m=0.1 "
            "product_variance_m2=0.04 applied_variance_m2=0 variance_mode=CONDITIONED"
        )
        code = (
            "PPP_EXTERNAL_BIAS_APPLICATION type=CODE time=2024-07-17 00:00:00.000 "
            "receiver=GRAZ satellite=G02 signal=L1C bias_found=1 bias_m=0.2 "
            "product_variance_m2=0.0025 applied_variance_m2=0.0025 variance_mode=PER_EPOCH"
        )
        with tempfile.TemporaryDirectory() as temporary:
            trace = Path(temporary) / "Network.trace"
            trace.write_text("\n".join([phase, phase, code, "PPP_EXTERNAL_BIAS_APPLICATION bad"]), encoding="utf-8")
            report = audit_trace(trace)

        self.assertEqual(report["marker_count"], 3)
        self.assertEqual(report["unique_application_key_count"], 2)
        self.assertEqual(report["duplicate_application_key_count"], 1)
        self.assertEqual(report["malformed_marker_count"], 1)
        self.assertEqual(report["nonzero_applied_variance_count"], 1)
        self.assertEqual(
            report["counts_by_type_signal_mode"],
            {"CODE/L1C/PER_EPOCH": 1, "PHASE/L1C/CONDITIONED": 2},
        )
        self.assertAlmostEqual(
            report["product_variance_m2_by_type_signal"]["PHASE/L1C"]["mean"],
            0.04,
        )


if __name__ == "__main__":
    unittest.main()
