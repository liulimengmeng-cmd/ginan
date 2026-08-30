from __future__ import annotations

import hashlib
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

        manifest_root = self.repo / "experiment0_2024"
        manifest_root.mkdir()
        manifest_entries = []
        self.manifest_files: list[Path] = []
        for index in range(driver.EXPECTED_MANIFEST_FILE_COUNT):
            relative = Path("selected") / f"input_{index:02d}.dat"
            path = manifest_root / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            payload = f"manifest fixture {index}\n".encode()
            path.write_bytes(payload)
            self.manifest_files.append(path)
            manifest_entries.append(
                {
                    "experiment": "shared",
                    "relative_path": relative.as_posix(),
                    "role": "test fixture",
                    "size": len(payload),
                    "sha256": hashlib.sha256(payload).hexdigest(),
                    "source": "local test fixture",
                    "source_metadata": {},
                }
            )
        manifest = manifest_root / "input_manifest.json"
        manifest.write_text(
            json.dumps(
                {
                    "schema": driver.EXPECTED_MANIFEST_SCHEMA,
                    "dates": driver.EXPECTED_DATES,
                    "model_stations": sorted(driver.EXPECTED_MODEL_STATIONS),
                    "heldout_stations": sorted(driver.EXPECTED_HELDOUT_STATIONS),
                    "stations": sorted(
                        driver.EXPECTED_MODEL_STATIONS
                        | driver.EXPECTED_HELDOUT_STATIONS
                    ),
                    "observation_policy": {
                        "constellation": "GPS",
                        "nominal_interval_seconds": 30,
                        "signals": ["C1C", "L1C", "C2W", "L2W"],
                    },
                    "precise_product_policy": {
                        "analysis_center": "Wuhan University",
                        "antenna_calibration": "shared/products/igs20_2303.atx",
                        "bias_apc_model": "IGS20_2303.ATX",
                        "series": "WUM0MGXFIN",
                    },
                    "files": manifest_entries,
                }
            ),
            encoding="utf-8",
        )
        self.documents: dict[str, Path] = {"manifest.input": manifest}
        for role in (
            "audit.rinex_qc",
            "audit.split",
            "audit.product",
            "audit.validation_split",
        ):
            path = manifest_root / f"{role.replace('.', '_')}.json"
            path.write_text(json.dumps({"role": role}), encoding="utf-8")
            self.documents[role] = path
        manifest_sha = hashlib.sha256(manifest.read_bytes()).hexdigest()
        self.documents["audit.rinex_qc"].write_text(
            json.dumps(
                {
                    "schema": driver.EXPECTED_RINEX_QC_SCHEMA,
                    "pass": True,
                    "expected_file_count": 26,
                    "audited_file_count": 26,
                    "files": [{"pass": True} for _ in range(26)],
                    "input_manifest": str(manifest.resolve()),
                    "input_manifest_sha256": manifest_sha,
                }
            ),
            encoding="utf-8",
        )
        self.documents["audit.split"].write_text(
            json.dumps(
                {
                    "schema": driver.EXPECTED_SPLIT_AUDIT_SCHEMA,
                    "valid": True,
                    "dates": driver.EXPECTED_DATES,
                    "train": sorted(driver.EXPECTED_MODEL_STATIONS),
                    "heldout": sorted(driver.EXPECTED_HELDOUT_STATIONS),
                    "intersection": [],
                    "roles_identical": True,
                    "manifest": str(manifest.resolve()),
                    "per_day": {
                        day: {
                            "yaml": str(self.yamls[f"{day}_model"].resolve()),
                            "rnx_input_count": 10,
                            "train": sorted(driver.EXPECTED_MODEL_STATIONS),
                            "heldout_input_count": 0,
                            "heldout_inputs": [],
                            "duplicate_stations": [],
                            "invalid_inputs": [],
                            "missing_train": [],
                            "unexpected_stations": [],
                        }
                        for day in driver.DAY_LABELS
                    },
                    "errors": [],
                }
            ),
            encoding="utf-8",
        )
        self.documents["audit.product"].write_text(
            json.dumps(
                {
                    "schema": driver.EXPECTED_PRODUCT_AUDIT_SCHEMA,
                    "validation_passed": False,
                    "experiments": {
                        day: {
                            "validation": {
                                "eligible_intersection_passed": True,
                                "eligible_satellites_with_coverage_defects": [
                                    f"G{index:02d}"
                                    for index in range(32 if day == "quiet" else 2, 33)
                                ],
                                "full_day_completeness_passed": False,
                                "metadata": {"passed": True},
                                "missing_expected_satellites": [],
                                "passed": False,
                                "unexpected_satellites": [],
                            }
                        }
                        for day in driver.DAY_LABELS
                    },
                }
            ),
            encoding="utf-8",
        )
        self.documents["audit.validation_split"].write_text(
            json.dumps(
                {
                    "schema": driver.EXPECTED_VALIDATION_SPLIT_SCHEMA,
                    "valid": True,
                    "manifest": str(manifest.resolve()),
                    "model_split_valid": True,
                    "model_stations": sorted(driver.EXPECTED_MODEL_STATIONS),
                    "heldout_stations": sorted(driver.EXPECTED_HELDOUT_STATIONS),
                    "intersection": [],
                    "roles_identical": True,
                    "per_day": {
                        day: {
                            "yaml": str(self.yamls[f"{day}_heldout"].resolve()),
                            "rnx_input_count": 3,
                            "stations": sorted(driver.EXPECTED_HELDOUT_STATIONS),
                            "model_input_count": 0,
                            "model_inputs": [],
                            "missing_heldout": [],
                            "unexpected_stations": [],
                            "paths_match_manifest": True,
                        }
                        for day in driver.DAY_LABELS
                    },
                    "errors": [],
                }
            ),
            encoding="utf-8",
        )

        self.models: dict[str, dict[str, Path]] = {}
        for day in driver.DAY_LABELS:
            day_paths = self._day_paths(self.repo, f"{day}_model")
            self._write_sidecars(day_paths)
            self._write_audits(day, day_paths)
            day_paths["console"].write_text("historical model log\n", encoding="utf-8")
            self.models[day] = day_paths

        self.pea = self.repo / "pea"
        self.pea.write_bytes(b"fake pea 1296966")
        self.freeze_path = self.repo / "freeze.json"
        self.freeze_receipt = self.repo / "freeze.receipt.json"
        self.report_path = self.repo / "report.json"
        self.evaluation_receipt = self.repo / "evaluation.receipt.json"
        self.evaluation_claim = self.repo / "evaluation.claim.json"
        self.heldout: dict[str, dict[str, Path]] = {}
        for day in driver.DAY_LABELS:
            root = self.base / f"{day}_heldout_output"
            paths = self._day_paths(root, day)
            paths["output_root"] = root
            paths["run_receipt"] = root / f"{day}.run.receipt.json"
            self.heldout[day] = paths
        self.plan_path = self.repo / "blind_plan.json"
        self.plan = self._plan_payload(common)
        self._write_plan()

    def tearDown(self) -> None:
        self.temporary_directory.cleanup()

    @staticmethod
    def _day_paths(root: Path, stem: str) -> dict[str, Path]:
        return {
            "stec": root / "ionstec" / f"{stem}.STEC",
            "raw_cov": root / "ionstec" / f"{stem}.STEC.COV",
            "sd": root / "ionstec" / f"{stem}.STEC.SD.COV",
            "raw_audit": root / "audit" / f"{stem}.raw.json",
            "sd_audit": root / "audit" / f"{stem}.sd.json",
            "console": root / f"{stem}.console.log",
        }

    @staticmethod
    def _write_sidecars(paths: dict[str, Path]) -> None:
        paths["stec"].parent.mkdir(parents=True, exist_ok=True)
        paths["stec"].write_text("geometry\n", encoding="utf-8")
        paths["raw_cov"].write_text("raw covariance\n", encoding="utf-8")
        paths["sd"].write_text(SD_FORWARD, encoding="utf-8")

    @staticmethod
    def _audit_payloads(
        day: str, paths: dict[str, Path]
    ) -> tuple[dict[str, object], dict[str, object]]:
        invalid_count = 0 if day == "quiet" else 7
        invalid_indices = set(range(120, 120 + invalid_count))
        raw_epochs = []
        sd_epochs = []
        start_tow = int(driver.EXPECTED_DAY_START_KEYS[day][1])
        for index in range(driver.EXPECTED_EPOCH_COUNT):
            key = {
                "gps_week": 2313,
                "gps_tow": start_tow + index * 30,
                "posterior_stage": driver.EXPECTED_POSTERIOR_STAGE,
            }
            if index in invalid_indices:
                raw_epochs.append(
                    {
                        **key,
                        "valid": False,
                        "writer_status": "NO_STATES",
                        "state_count": 0,
                        "errors": ["writer_status=NO_STATES"],
                    }
                )
                sd_epochs.append(
                    {
                        **key,
                        "valid": False,
                        "source_state_count": 0,
                        "difference_state_count": 0,
                        "errors": [
                            "raw_writer_status=NO_STATES",
                            "sd_writer_status=NO_STATES",
                        ],
                    }
                )
            else:
                raw_epochs.append({**key, "valid": True})
                sd_epochs.append({**key, "valid": True})
        valid_count = driver.EXPECTED_EPOCH_COUNT - invalid_count
        raw = {
            "schema": driver.sd_parser.RAW_SCHEMA,
            "input": str(paths["raw_cov"].resolve()),
            "epoch_count": driver.EXPECTED_EPOCH_COUNT,
            "valid_epoch_count": valid_count,
            "invalid_epoch_count": invalid_count,
            "all_valid": invalid_count == 0,
            "epochs": raw_epochs,
        }
        sd = {
            "schema": driver.sd_parser.REPORT_SCHEMA,
            "raw_schema": driver.sd_parser.RAW_SCHEMA,
            "satellite_difference_schema": driver.sd_parser.SD_SCHEMA,
            "raw_input": str(paths["raw_cov"].resolve()),
            "satellite_difference_input": str(paths["sd"].resolve()),
            "raw_epoch_count": driver.EXPECTED_EPOCH_COUNT,
            "satellite_difference_epoch_count": driver.EXPECTED_EPOCH_COUNT,
            "matched_epoch_count": driver.EXPECTED_EPOCH_COUNT,
            "valid_epoch_count": valid_count,
            "invalid_epoch_count": invalid_count,
            "all_valid": invalid_count == 0,
            "file_errors": [],
            "epochs": sd_epochs,
        }
        return raw, sd

    def _write_audits(self, day: str, paths: dict[str, Path]) -> None:
        raw, sd = self._audit_payloads(day, paths)
        paths["raw_audit"].parent.mkdir(parents=True, exist_ok=True)
        paths["raw_audit"].write_text(json.dumps(raw), encoding="utf-8")
        paths["sd_audit"].write_text(json.dumps(sd), encoding="utf-8")

    def _plan_payload(self, common: Path) -> dict[str, object]:
        static_inputs = {
            **{role: str(path) for role, path in driver.ACTUAL_CODE_PATHS.items()},
            **{role: str(path) for role, path in self.documents.items()},
            "yaml.quiet_model": str(self.yamls["quiet_model"]),
            "yaml.storm_model": str(self.yamls["storm_model"]),
            "yaml.quiet_heldout": str(self.yamls["quiet_heldout"]),
            "yaml.storm_heldout": str(self.yamls["storm_heldout"]),
            "yaml.quiet_model_inputs": str(self.yamls["quiet_model"]),
            "yaml.storm_model_inputs": str(self.yamls["storm_model"]),
            "yaml.include.common": str(common),
        }
        data_root = self.documents["manifest.input"].parent
        return {
            "schema": driver.PLAN_SCHEMA,
            "repository": str(self.repo),
            "static_inputs": static_inputs,
            "model": {
                day: {field: str(paths[field]) for field in driver.MODEL_FIELDS}
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
                        f"EXP0_DATA_ROOT:{data_root}",
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
            "audit_policy": {
                day: {
                    "expected_invalid_epochs": [
                        {
                            "gps_week": 2313,
                            "gps_tow": int(driver.EXPECTED_DAY_START_KEYS[day][1])
                            + index * 30,
                        }
                        for index in ([] if day == "quiet" else range(120, 127))
                    ]
                }
                for day in driver.DAY_LABELS
            },
            "freeze_output": str(self.freeze_path),
            "freeze_receipt_output": str(self.freeze_receipt),
            "report_output": str(self.report_path),
            "evaluation_receipt_output": str(self.evaluation_receipt),
            "evaluation_claim_output": str(self.evaluation_claim),
        }

    def _write_plan(self) -> None:
        self.plan_path.write_text(
            json.dumps(self.plan, indent=2, sort_keys=True), encoding="utf-8"
        )

    def _freeze(self) -> tuple[dict[str, object], object]:
        payload = {
            "schema": driver.spatial.FREEZE_SCHEMA,
            "heldout_output_read": False,
            "scientific_scope": "test",
        }
        patcher = patch.object(driver.spatial, "build_freeze", return_value=payload)
        mocked = patcher.start()
        self.addCleanup(patcher.stop)
        audit_patchers = (
            patch.object(
                driver.rinex_audit,
                "build_audit",
                return_value=json.loads(
                    self.documents["audit.rinex_qc"].read_text(encoding="utf-8")
                ),
            ),
            patch.object(
                driver.split_audit,
                "build_audit",
                return_value=json.loads(
                    self.documents["audit.split"].read_text(encoding="utf-8")
                ),
            ),
            patch.object(
                driver.product_audit,
                "audit_data_root",
                return_value=json.loads(
                    self.documents["audit.product"].read_text(encoding="utf-8")
                ),
            ),
            patch.object(
                driver.validation_split_audit,
                "build_validation_audit",
                return_value=json.loads(
                    self.documents["audit.validation_split"].read_text(encoding="utf-8")
                ),
            ),
        )
        for audit_patcher in audit_patchers:
            audit_patcher.start()
            self.addCleanup(audit_patcher.stop)
        result = driver.freeze_from_plan(self.plan_path)
        self.assertEqual(mocked.call_count, 1)
        return result, mocked

    def _fake_execution(
        self,
        command: list[str],
        working_directory: Path,
        _environment: dict[str, str],
        console: Path,
        *,
        exclusive: bool,
    ) -> dict[str, object]:
        is_pea = command[0] == "/usr/bin/time"
        text = "PEA finished\nExit status: 0\n" if is_pea else "validator\n"
        console.parent.mkdir(parents=True, exist_ok=True)
        with console.open("wb" if exclusive else "ab") as stream:
            stream.write(text.encode())
        if is_pea:
            day = "quiet" if "quiet" in str(console) else "storm"
            self.assertTrue(self.heldout[day]["output_root"].is_dir())
            self._write_sidecars(self.heldout[day])
            return {
                "command": command,
                "working_directory": str(working_directory.resolve()),
                "returncode": 0,
                "started_utc": "a",
                "ended_utc": "b",
            }
        if "--json-output" in command:
            audit_path = Path(command[command.index("--json-output") + 1])
            day = "quiet" if "quiet" in str(audit_path) else "storm"
            raw, _ = self._audit_payloads(day, self.heldout[day])
            audit_path.write_text(json.dumps(raw), encoding="utf-8")
            return {
                "command": command,
                "working_directory": str(working_directory.resolve()),
                "returncode": 0 if day == "quiet" else 1,
            }
        audit_path = Path(command[-1])
        day = "quiet" if "quiet" in str(audit_path) else "storm"
        _, sd = self._audit_payloads(day, self.heldout[day])
        audit_path.write_text(json.dumps(sd), encoding="utf-8")
        return {
            "command": command,
            "working_directory": str(working_directory.resolve()),
            "returncode": 0 if day == "quiet" else 1,
        }

    def _run_both(self, expected_sha: str) -> dict[str, str]:
        receipt_hashes: dict[str, str] = {}
        with patch.object(
            driver, "_execute_logged", side_effect=self._fake_execution
        ) as execute:
            for day in driver.DAY_LABELS:
                result = driver.run_frozen_day(
                    self.freeze_path,
                    self.freeze_receipt,
                    expected_sha,
                    day,
                )
                self.assertEqual(result["status"], "SUCCESS")
                receipt_hashes[day] = result["run_receipt"]["sha256"]
        self.assertEqual(execute.call_count, 6)
        return receipt_hashes

    def test_freeze_run_and_evaluate_closed_loop(self) -> None:
        freeze_result, freeze_mock = self._freeze()
        self.assertFalse(freeze_result["heldout_output_read"])
        receipt_payload = json.loads(self.freeze_receipt.read_text(encoding="utf-8"))
        manifest_roles = [
            role
            for role in receipt_payload["provenance"]["inputs_pre"]
            if role.startswith("manifest.file.")
        ]
        self.assertEqual(len(manifest_roles), driver.EXPECTED_MANIFEST_FILE_COUNT)
        self.assertEqual(freeze_mock.call_count, 1)

        expected_sha = record_regular_file(self.freeze_path)["sha256"]
        run_hashes = self._run_both(expected_sha)
        self.assertTrue(self.heldout["quiet"]["run_receipt"].is_file())
        with patch.object(
            driver.spatial,
            "build_report",
            return_value={"schema": "TEST_SPATIAL_REPORT"},
        ) as build_report:
            result = driver.evaluate_frozen(
                self.freeze_path,
                self.freeze_receipt,
                expected_sha,
                run_hashes["quiet"],
                run_hashes["storm"],
                self.report_path,
                self.evaluation_receipt,
            )
        build_report.assert_called_once()
        report = json.loads(self.report_path.read_text(encoding="utf-8"))
        self.assertTrue(report["blind_provenance"]["inputs_unchanged"])
        self.assertEqual(
            report["blind_validation"]["storm"]["audits"]["invalid_epoch_count"], 7
        )
        self.assertTrue(self.evaluation_claim.is_file())
        self.assertEqual(result["report"]["path"], str(self.report_path.resolve()))

    def test_manifest_content_mismatch_blocks_freeze(self) -> None:
        self.manifest_files[17].write_text("changed\n", encoding="utf-8")
        with patch.object(driver.spatial, "build_freeze") as build_freeze:
            with self.assertRaises(driver.BlindDriverError) as context:
                driver.freeze_from_plan(self.plan_path)
        build_freeze.assert_not_called()
        self.assertIn("manifest file identity mismatch", str(context.exception))

    def test_aliases_are_exact_and_unique(self) -> None:
        self.plan["pea"]["quiet_argv"].extend(
            ["-a", f"EXP0_DATA_ROOT:{self.documents['manifest.input'].parent}"]
        )
        self._write_plan()
        with self.assertRaises(driver.BlindDriverError) as context:
            driver._normalise_freeze_plan(self.plan_path)
        self.assertIn("repeats alias EXP0_DATA_ROOT", str(context.exception))

    def test_extra_yaml_or_alias_is_rejected(self) -> None:
        self.plan["pea"]["quiet_argv"].extend(
            ["-y", str(self.yamls["quiet_model"]), "-a", "EXTRA:value"]
        )
        self._write_plan()
        with self.assertRaises(driver.BlindDriverError):
            driver._normalise_freeze_plan(self.plan_path)

    def test_quiet_and_storm_audit_policy_and_key_match(self) -> None:
        model_strings = {
            day: {field: str(path) for field, path in paths.items()}
            for day, paths in self.models.items()
        }
        quiet_expected: set[tuple[int, driver.Decimal]] = set()
        storm_expected = {
            (
                2313,
                driver.EXPECTED_DAY_START_KEYS["storm"][1] + driver.Decimal(index * 30),
            )
            for index in range(120, 127)
        }
        quiet = driver._validate_day_audits(
            "quiet", model_strings["quiet"], "quiet", quiet_expected
        )
        storm = driver._validate_day_audits(
            "storm", model_strings["storm"], "storm", storm_expected
        )
        self.assertEqual(quiet["valid_epoch_count"], 2880)
        self.assertEqual(storm["invalid_epoch_count"], 7)

        sd = json.loads(self.models["storm"]["sd_audit"].read_text(encoding="utf-8"))
        first_invalid = next(epoch for epoch in sd["epochs"] if not epoch["valid"])
        first_invalid["gps_tow"] += 1
        self.models["storm"]["sd_audit"].write_text(json.dumps(sd), encoding="utf-8")
        with self.assertRaises(driver.BlindDriverError) as context:
            driver._validate_day_audits(
                "storm", model_strings["storm"], "storm", storm_expected
            )
        self.assertIn("complete epoch key order differs", str(context.exception))

    def test_artifact_hash_change_blocks_report(self) -> None:
        self._freeze()
        expected_sha = record_regular_file(self.freeze_path)["sha256"]
        run_hashes = self._run_both(expected_sha)
        self.heldout["quiet"]["stec"].write_text("changed", encoding="utf-8")
        with patch.object(driver.spatial, "build_report") as build_report:
            with self.assertRaises(driver.BlindDriverError) as context:
                driver.evaluate_frozen(
                    self.freeze_path,
                    self.freeze_receipt,
                    expected_sha,
                    run_hashes["quiet"],
                    run_hashes["storm"],
                    self.report_path,
                    self.evaluation_receipt,
                )
        build_report.assert_not_called()
        self.assertIn("artifact changed after run", str(context.exception))
        self.assertTrue(self.evaluation_claim.is_file())

    def test_run_root_is_an_atomic_claim(self) -> None:
        self._freeze()
        expected_sha = record_regular_file(self.freeze_path)["sha256"]
        self.heldout["quiet"]["output_root"].mkdir()
        with patch.object(driver, "_execute_logged") as execute:
            with self.assertRaises((driver.BlindDriverError, FileExistsError)):
                driver.run_frozen_day(
                    self.freeze_path,
                    self.freeze_receipt,
                    expected_sha,
                    "quiet",
                )
        execute.assert_not_called()

    def test_validator_failure_writes_terminal_receipt(self) -> None:
        self._freeze()
        expected_sha = record_regular_file(self.freeze_path)["sha256"]

        def fatal_raw_validator(*args, **kwargs):
            result = self._fake_execution(*args, **kwargs)
            command = args[0]
            if "--json-output" in command:
                result["returncode"] = 2
            return result

        with patch.object(driver, "_execute_logged", side_effect=fatal_raw_validator):
            with self.assertRaises(driver.BlindDriverError):
                driver.run_frozen_day(
                    self.freeze_path,
                    self.freeze_receipt,
                    expected_sha,
                    "quiet",
                )
        receipt = json.loads(
            self.heldout["quiet"]["run_receipt"].read_text(encoding="utf-8")
        )
        self.assertEqual(receipt["status"], "VALIDATION_FAILED")
        self.assertEqual(receipt["failure"]["stage"], "VALIDATOR_EXECUTION")

    def test_validator_exception_writes_terminal_receipt(self) -> None:
        self._freeze()
        expected_sha = record_regular_file(self.freeze_path)["sha256"]

        def validator_launch_error(*args, **kwargs):
            if "--json-output" in args[0]:
                raise OSError("validator launch failed")
            return self._fake_execution(*args, **kwargs)

        with patch.object(
            driver, "_execute_logged", side_effect=validator_launch_error
        ):
            with self.assertRaises(driver.BlindDriverError):
                driver.run_frozen_day(
                    self.freeze_path,
                    self.freeze_receipt,
                    expected_sha,
                    "quiet",
                )
        receipt = json.loads(
            self.heldout["quiet"]["run_receipt"].read_text(encoding="utf-8")
        )
        self.assertEqual(receipt["status"], "VALIDATION_FAILED")
        self.assertEqual(receipt["failure"]["stage"], "VALIDATOR_LAUNCH_OR_LOGGING")

    def test_console_publication_exception_writes_terminal_receipt(self) -> None:
        self._freeze()
        expected_sha = record_regular_file(self.freeze_path)["sha256"]
        with (
            patch.object(driver, "_execute_logged", side_effect=self._fake_execution),
            patch.object(
                driver,
                "_publish_existing_file_no_replace",
                side_effect=OSError("console publication failed"),
            ),
        ):
            with self.assertRaises(driver.BlindDriverError):
                driver.run_frozen_day(
                    self.freeze_path,
                    self.freeze_receipt,
                    expected_sha,
                    "quiet",
                )
        receipt = json.loads(
            self.heldout["quiet"]["run_receipt"].read_text(encoding="utf-8")
        )
        self.assertEqual(receipt["status"], "EXECUTION_FAILED")
        self.assertEqual(
            receipt["failure"]["stage"], "PEA_EXECUTION_OR_CONSOLE_PUBLICATION"
        )

    def test_receipt_command_mismatch_blocks_report(self) -> None:
        self._freeze()
        expected_sha = record_regular_file(self.freeze_path)["sha256"]
        run_hashes = self._run_both(expected_sha)
        receipt_path = self.heldout["quiet"]["run_receipt"]
        receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
        receipt["execution"]["command"] = ["not-the-frozen-command"]
        receipt_path.write_text(json.dumps(receipt), encoding="utf-8")
        run_hashes["quiet"] = record_regular_file(receipt_path)["sha256"]
        with patch.object(driver.spatial, "build_report") as build_report:
            with self.assertRaises(driver.BlindDriverError) as context:
                driver.evaluate_frozen(
                    self.freeze_path,
                    self.freeze_receipt,
                    expected_sha,
                    run_hashes["quiet"],
                    run_hashes["storm"],
                    self.report_path,
                    self.evaluation_receipt,
                )
        build_report.assert_not_called()
        self.assertIn("wrong PEA command", str(context.exception))

    def test_existing_evaluation_claim_blocks_second_evaluation(self) -> None:
        self._freeze()
        expected_sha = record_regular_file(self.freeze_path)["sha256"]
        self.evaluation_claim.write_text("claimed\n", encoding="utf-8")
        with patch.object(driver.spatial, "build_report") as build_report:
            with self.assertRaises(FileExistsError):
                driver.evaluate_frozen(
                    self.freeze_path,
                    self.freeze_receipt,
                    expected_sha,
                    "0" * 64,
                    "0" * 64,
                    self.report_path,
                    self.evaluation_receipt,
                )
        build_report.assert_not_called()

    def test_evaluate_rejects_nonfrozen_output_path(self) -> None:
        self._freeze()
        expected_sha = record_regular_file(self.freeze_path)["sha256"]
        with patch.object(driver.spatial, "build_report") as build_report:
            with self.assertRaises(driver.BlindDriverError) as context:
                driver.evaluate_frozen(
                    self.freeze_path,
                    self.freeze_receipt,
                    expected_sha,
                    "0" * 64,
                    "0" * 64,
                    self.repo / "alternate-report.json",
                    self.evaluation_receipt,
                )
        build_report.assert_not_called()
        self.assertIn("must exactly equal", str(context.exception))

    def test_wrong_external_hash_is_rejected_before_heldout_access(self) -> None:
        self._freeze()
        with patch.object(driver, "_heldout_inputs") as heldout_inputs:
            with self.assertRaises(driver.BlindDriverError):
                driver.evaluate_frozen(
                    self.freeze_path,
                    self.freeze_receipt,
                    "0" * 64,
                    "0" * 64,
                    "0" * 64,
                    self.report_path,
                    self.evaluation_receipt,
                )
        heldout_inputs.assert_not_called()

    def test_help_exposes_three_subcommands(self) -> None:
        help_text = driver.build_parser().format_help()
        for command in ("freeze", "run", "evaluate"):
            self.assertIn(command, help_text)
        self.assertIn(driver.PLAN_SCHEMA, help_text)

    def _git(self, *arguments: str) -> None:
        subprocess.run(
            ["git", "-C", str(self.repo), *arguments],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            encoding="utf-8",
        )


if __name__ == "__main__":
    unittest.main()
