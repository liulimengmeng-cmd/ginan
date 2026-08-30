from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from test_validate_experiment0_2024_split import (
    EXPECTED_HELDOUT,
    EXPECTED_TRAIN,
    manifest_payload,
    model_yaml,
)
from validate_experiment0_2024_validation_split import build_validation_audit


class Experiment0ValidationSplitAuditTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary_directory.name)
        self.manifest = self.root / "input_manifest.json"
        self.quiet_model = self.root / "quiet_model.yaml"
        self.storm_model = self.root / "storm_model.yaml"
        self.quiet_heldout = self.root / "quiet_heldout.yaml"
        self.storm_heldout = self.root / "storm_heldout.yaml"
        self.manifest.write_text(json.dumps(manifest_payload()), encoding="utf-8")
        self.quiet_model.write_text(
            model_yaml("quiet", EXPECTED_TRAIN, 129), encoding="utf-8"
        )
        self.storm_model.write_text(
            model_yaml("storm", EXPECTED_TRAIN, 132), encoding="utf-8"
        )
        self.quiet_heldout.write_text(
            model_yaml("quiet", EXPECTED_HELDOUT, 129), encoding="utf-8"
        )
        self.storm_heldout.write_text(
            model_yaml("storm", EXPECTED_HELDOUT, 132), encoding="utf-8"
        )

    def tearDown(self) -> None:
        self.temporary_directory.cleanup()

    def test_accepts_exact_disjoint_validation_split(self) -> None:
        audit = build_validation_audit(
            self.manifest,
            self.quiet_model,
            self.storm_model,
            self.quiet_heldout,
            self.storm_heldout,
        )
        self.assertTrue(audit["valid"])
        self.assertTrue(audit["model_split_valid"])
        self.assertEqual(audit["intersection"], [])
        self.assertTrue(audit["roles_identical"])
        self.assertEqual(audit["per_day"]["quiet"]["model_input_count"], 0)
        self.assertTrue(audit["per_day"]["storm"]["paths_match_manifest"])

    def test_rejects_model_station_in_validation_yaml(self) -> None:
        contaminated = EXPECTED_HELDOUT[:-1] + (EXPECTED_TRAIN[0],)
        self.quiet_heldout.write_text(
            model_yaml("quiet", contaminated, 129), encoding="utf-8"
        )
        audit = build_validation_audit(
            self.manifest,
            self.quiet_model,
            self.storm_model,
            self.quiet_heldout,
            self.storm_heldout,
        )
        self.assertFalse(audit["valid"])
        self.assertEqual(audit["per_day"]["quiet"]["model_inputs"], ["ARMC"])

    def test_rejects_wrong_day_heldout_paths(self) -> None:
        self.quiet_heldout.write_text(
            model_yaml("storm", EXPECTED_HELDOUT, 132), encoding="utf-8"
        )
        audit = build_validation_audit(
            self.manifest,
            self.quiet_model,
            self.storm_model,
            self.quiet_heldout,
            self.storm_heldout,
        )
        self.assertFalse(audit["valid"])
        self.assertFalse(audit["per_day"]["quiet"]["paths_match_manifest"])


if __name__ == "__main__":
    unittest.main()
