from __future__ import annotations

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parent))

import experiment0_2024_blind_driver as driver
from experiment0_2024_blind_provenance import record_regular_file


SD_FORWARD = """# GINAN_STEC_SATELLITE_DIFFERENCE_COVARIANCE_V1
# META,gps_week,gps_tow,status,source_state_count,difference_state_count,datum_count,singleton_datum_count,upper_triangle_count,max_abs_asymmetry_tecu2,posterior_stage
META,2313,0,NO_STATES,0,0,0,0,0,0,FILTER_POSTERIOR_NO_EPOCH_AR
"""


class Experiment0BlindDriverTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.base = Path(self.temporary_directory.name)
        self.repo = self.base / "repo"
        self.repo.mkdir()
        self._git("init")
        self._git("config", "user.name", "Blind Driver Test")
        self._git("config", "user.email", "test@example.invalid")
        self._git("commit", "--allow-empty", "-m", "fixture")

        common = self.repo / "common.yaml"
        common.write_text("inputs: {}\n", encoding="utf-8")
        self.yamls: dict[str, Path] = {}
        for role in (
            "quiet_model",
            "storm_model",
            "quiet_heldout",
            "storm_heldout",
        ):
            path = self.repo / f"{role}.yaml"
            path.write_text(
                "inputs:\n  include_yamls:\n    - common.yaml\n", encoding="utf-8"
            )
            self.yamls[role] = path

        self.documents: dict[str, Path] = {}
        for role in (
            "manifest.input",
            "audit.rinex_qc",
            "audit.split",
            "audit.product",
            "audit.validation_split",
        ):
            path = self.repo / f"{role.replace('.', '_')}.json"
            path.write_text(json.dumps({"role": role}), encoding="utf-8")
            self.documents[role] = path

        self.models: dict[str, dict[str, Path]] = {}
        for day in driver.DAY_LABELS:
            day_paths = {
                "stec": self.repo / f"{day}_model.STEC",
                "sd": self.repo / f"{day}_model.STEC.SD.COV",
                "audit": self.repo / f"{day}_model.audit.json",
                "console": self.repo / f"{day}_model.console.log",
            }
            day_paths["stec"].write_text("model geometry\n", encoding="utf-8")
            day_paths["sd"].write_text(SD_FORWARD, encoding="utf-8")
            day_paths["audit"].write_text(
                json.dumps({"all_valid": False}), encoding="utf-8"
            )
            day_paths["console"].write_text("historical model log\n", encoding="utf-8")
            self.models[day] = day_paths

        self.pea = self.repo / "pea"
        self.pea.write_bytes(b"fake pea 1296966")
        self.freeze_path = self.repo / "freeze.json"
        self.freeze_receipt = self.repo / "freeze.receipt.json"
        self.report_path = self.repo / "report.json"
        self.evaluation_receipt = self.repo / "evaluation.receipt.json"
        self.heldout: dict[str, dict[str, Path]] = {}
        for day in driver.DAY_LABELS:
            root = self.base / f"{day}_heldout_output"
            self.heldout[day] = {
                "output_root": root,
                "stec": root / "ionstec" / f"{day}.STEC",
                "raw_cov": root / "ionstec" / f"{day}.STEC.COV",
                "sd": root / "ionstec" / f"{day}.STEC.SD.COV",
                "raw_audit": root / f"{day}.raw.audit.json",
                "sd_audit": root / f"{day}.sd.audit.json",
                "console": root / f"{day}.console.log",
            }
        self.plan_path = self.repo / "blind_plan.json"
        self.plan = self._plan_payload(common)
        self._write_plan()

    def tearDown(self) -> None:
        self.temporary_directory.cleanup()

    def _plan_payload(self, common: Path) -> dict[str, object]:
        static_inputs = {
            **{role: str(path) for role, path in driver.ACTUAL_CODE_PATHS.items()},
            **{role: str(path) for role, path in self.documents.items()},
            "yaml.quiet_model": str(self.yamls["quiet_model"]),
            "yaml.storm_model": str(self.yamls["storm_model"]),
            "yaml.quiet_heldout": str(self.yamls["quiet_heldout"]),
            "yaml.storm_heldout": str(self.yamls["storm_heldout"]),
            "yaml.include.common": str(common),
        }
        return {
            "schema": driver.PLAN_SCHEMA,
            "repository": str(self.repo),
            "static_inputs": static_inputs,
            "model": {
                day: {field: str(path) for field, path in paths.items()}
                for day, paths in self.models.items()
            },
            "pea": {
                "binary": str(self.pea),
                "commit": "1296966",
                "environment": {"LD_LIBRARY_PATH": str(self.repo / "lib")},
                **{
                    f"{day}_argv": [
                        str(self.pea),
                        "-y",
                        str(self.yamls[f"{day}_heldout"]),
                        "-a",
                        f"EXP0_OUTPUT_ROOT:{self.heldout[day]['output_root']}",
                    ]
                    for day in driver.DAY_LABELS
                },
            },
            "heldout": {
                day: {field: str(path) for field, path in paths.items()}
                for day, paths in self.heldout.items()
            },
            "freeze_output": str(self.freeze_path),
            "freeze_receipt_output": str(self.freeze_receipt),
        }

    def _write_plan(self) -> None:
        self.plan_path.write_text(
            json.dumps(self.plan, indent=2, sort_keys=True), encoding="utf-8"
        )

    def _freeze(self):
        payload = {
            "schema": spatial_freeze_schema(),
            "heldout_output_read": False,
            "scientific_scope": "test",
        }
        patcher = patch.object(driver.spatial, "build_freeze", return_value=payload)
        mocked = patcher.start()
        self.addCleanup(patcher.stop)
        result = driver.freeze_from_plan(self.plan_path)
        self.assertEqual(mocked.call_count, 1)
        return result, mocked

    def _materialise_heldout(self, *, all_valid: bool = True, stage=SD_FORWARD) -> None:
        for day in driver.DAY_LABELS:
            paths = self.heldout[day]
            paths["stec"].parent.mkdir(parents=True, exist_ok=True)
            paths["stec"].write_text("heldout geometry\n", encoding="utf-8")
            paths["raw_cov"].write_text("raw covariance\n", encoding="utf-8")
            paths["sd"].write_text(stage, encoding="utf-8")
            for field in ("raw_audit", "sd_audit"):
                paths[field].write_text(
                    json.dumps({"schema": field, "all_valid": all_valid}),
                    encoding="utf-8",
                )
            paths["console"].write_text(
                "PEA finished\nExit status: 0\n", encoding="utf-8"
            )

    def test_freeze_and_evaluate_closed_loop(self) -> None:
        freeze_result, freeze_mock = self._freeze()
        self.assertFalse(freeze_result["heldout_output_read"])
        self.assertFalse(
            any(paths["output_root"].exists() for paths in self.heldout.values())
        )
        freeze_payload = json.loads(self.freeze_path.read_text(encoding="utf-8"))
        receipt_payload = json.loads(self.freeze_receipt.read_text(encoding="utf-8"))
        self.assertEqual(receipt_payload["schema"], driver.FREEZE_RECEIPT_SCHEMA)
        self.assertEqual(
            freeze_payload["blind_provenance"], receipt_payload["provenance"]
        )
        self.assertTrue(receipt_payload["provenance"]["inputs_unchanged"])
        self.assertIn("executable.pea", receipt_payload["provenance"]["inputs_pre"])
        self.assertEqual(freeze_mock.call_count, 1)

        self._materialise_heldout()
        expected_sha = record_regular_file(self.freeze_path)["sha256"]
        report_payload = {"schema": "TEST_SPATIAL_REPORT"}
        with patch.object(
            driver.spatial, "build_report", return_value=report_payload
        ) as build_report:
            result = driver.evaluate_frozen(
                self.freeze_path,
                self.freeze_receipt,
                expected_sha,
                self.report_path,
                self.evaluation_receipt,
            )
        build_report.assert_called_once()
        self.assertTrue(self.report_path.is_file())
        self.assertTrue(self.evaluation_receipt.is_file())
        report = json.loads(self.report_path.read_text(encoding="utf-8"))
        self.assertEqual(report["blind_anchor"]["expected_freeze_sha256"], expected_sha)
        self.assertTrue(report["blind_provenance"]["inputs_unchanged"])
        self.assertEqual(
            report["blind_validation"]["quiet"]["sd_stage"]["posterior_stage"],
            driver.EXPECTED_POSTERIOR_STAGE,
        )
        self.assertEqual(result["report"]["path"], str(self.report_path.resolve()))

    def test_freeze_rejects_incomplete_yaml_closure_before_scientific_call(
        self,
    ) -> None:
        del self.plan["static_inputs"]["yaml.include.common"]
        self._write_plan()
        with patch.object(driver.spatial, "build_freeze") as build_freeze:
            with self.assertRaises(driver.BlindDriverError) as context:
                driver.freeze_from_plan(self.plan_path)
        build_freeze.assert_not_called()
        self.assertIn("complete include closure", str(context.exception))

    def test_freeze_rejects_existing_heldout_root_and_existing_result(self) -> None:
        self.heldout["quiet"]["output_root"].mkdir()
        with patch.object(driver.spatial, "build_freeze") as build_freeze:
            with self.assertRaises(FileExistsError):
                driver.freeze_from_plan(self.plan_path)
        build_freeze.assert_not_called()

        self.heldout["quiet"]["output_root"].rmdir()
        self._freeze()
        original = self.freeze_path.read_bytes()
        with patch.object(driver.spatial, "build_freeze") as build_freeze:
            with self.assertRaises(FileExistsError):
                driver.freeze_from_plan(self.plan_path)
        build_freeze.assert_not_called()
        self.assertEqual(self.freeze_path.read_bytes(), original)

    def test_wrong_external_hash_is_rejected_before_heldout_access(self) -> None:
        self._freeze()
        with patch.object(driver, "_heldout_inputs") as heldout_inputs:
            with self.assertRaises(driver.BlindDriverError) as context:
                driver.evaluate_frozen(
                    self.freeze_path,
                    self.freeze_receipt,
                    "0" * 64,
                    self.report_path,
                    self.evaluation_receipt,
                )
        heldout_inputs.assert_not_called()
        self.assertIn("SHA-256 mismatch", str(context.exception))

    def test_invalid_audit_or_stage_blocks_report_publication(self) -> None:
        self._freeze()
        self._materialise_heldout(all_valid=False)
        expected_sha = record_regular_file(self.freeze_path)["sha256"]
        with patch.object(driver.spatial, "build_report") as build_report:
            with self.assertRaises(driver.BlindDriverError) as context:
                driver.evaluate_frozen(
                    self.freeze_path,
                    self.freeze_receipt,
                    expected_sha,
                    self.report_path,
                    self.evaluation_receipt,
                )
        build_report.assert_not_called()
        self.assertIn("all_valid=true", str(context.exception))
        self.assertFalse(self.report_path.exists())

    def test_help_exposes_both_subcommands(self) -> None:
        help_text = driver.build_parser().format_help()
        self.assertIn("freeze", help_text)
        self.assertIn("evaluate", help_text)
        self.assertIn("GINAN_EXPERIMENT0_2024_BLIND_PLAN_V1", help_text)

    def _git(self, *arguments: str) -> None:
        subprocess.run(
            ["git", "-C", str(self.repo), *arguments],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            encoding="utf-8",
        )


def spatial_freeze_schema() -> str:
    return driver.spatial.FREEZE_SCHEMA


if __name__ == "__main__":
    unittest.main()
