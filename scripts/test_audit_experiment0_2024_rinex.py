from __future__ import annotations

import hashlib
import json
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from audit_experiment0_2024_rinex import build_audit, parse_rinex


HEADER = """     3.04           OBSERVATION DATA    M                   RINEX VERSION / TYPE
G    4 C1C L1C C2W L2W                                      SYS / # / OBS TYPES
                                                            END OF HEADER
"""


class RinexAuditTests(unittest.TestCase):
    def test_build_audit_returns_deterministic_complete_payload(self) -> None:
        body = "> 2024 05 08 00 00  0.0000000  0  0\n"
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            rinex_path = root / "rinex" / "sample.rnx"
            rinex_path.parent.mkdir()
            rinex_path.write_text(HEADER + body, encoding="ascii")
            rinex_sha256 = hashlib.sha256(rinex_path.read_bytes()).hexdigest()
            manifest = {
                "stations": ["TEST"],
                "dates": ["2024-05-08"],
                "files": [
                    {
                        "relative_path": "rinex/sample.rnx",
                        "role": "GNSS RINEX 3 observation, 30 s",
                        "experiment": "quiet",
                        "sha256": rinex_sha256,
                        "source_metadata": {"siteId": "test"},
                    }
                ],
            }
            manifest_path = root / "input_manifest.json"
            manifest_path.write_text(json.dumps(manifest), encoding="utf-8")

            payload = build_audit(root)

        self.assertEqual(
            set(payload),
            {
                "schema",
                "input_manifest",
                "input_manifest_sha256",
                "criteria",
                "expected_file_count",
                "audited_file_count",
                "pass",
                "files",
            },
        )
        self.assertEqual(payload["schema"], "GINAN_EXPERIMENT0_2024_RINEX_QC_V1")
        self.assertEqual(payload["expected_file_count"], 1)
        self.assertEqual(payload["audited_file_count"], 1)
        self.assertFalse(payload["pass"])
        self.assertEqual(payload["files"][0]["station"], "TEST")
        self.assertTrue(payload["files"][0]["sha256_matches_manifest"])
        self.assertNotIn("generated_utc", payload)

    def test_parser_counts_epochs_and_gap(self) -> None:
        body = """> 2024 05 08 00 00  0.0000000  0  0
> 2024 05 08 00 00 30.0000000  0  0
> 2024 05 08 00 01 30.0000000  0  0
"""
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "sample.rnx"
            path.write_text(HEADER + body, encoding="ascii")
            result = parse_rinex(path)
        self.assertEqual(result["unique_epoch_count"], 3)
        self.assertEqual(result["max_gap_seconds"], 60.0)
        self.assertEqual(result["missing_required_gps_observables"], [])
        self.assertFalse(result["pass"])

    def test_parser_detects_missing_signal_and_duplicate(self) -> None:
        header = HEADER.replace(" C2W L2W", " C2W    ")
        epoch = "> 2024 05 08 00 00  0.0000000  0  0\n"
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "sample.rnx"
            path.write_text(header + epoch + epoch, encoding="ascii")
            result = parse_rinex(path)
        self.assertIn("L2W", result["missing_required_gps_observables"])
        self.assertEqual(result["duplicate_epoch_count"], 1)


if __name__ == "__main__":
    unittest.main()
