#!/usr/bin/env python3
"""Lightweight tests for the corrected Experiment 1 driver.

No test launches PEA.  Dry-run execution is guarded with a mock that fails if
the condition executor is reached.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import experiment1_2024_corrected_driver as driver


def _write(path: Path, payload: bytes) -> dict[str, object]:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(payload)
    return {
        "size": len(payload),
        "sha256": hashlib.sha256(payload).hexdigest(),
    }


def _make_manifest(path: Path, data_root: Path) -> None:
    files = []
    for index in range(driver.EXPECTED_INPUT_COUNT):
        relative = Path("inputs") / f"file_{index:02d}.dat"
        payload = f"payload-{index}\n".encode("ascii")
        record = _write(data_root / relative, payload)
        files.append(
            {
                "role": f"role_{index:02d}",
                "relative_path": relative.as_posix(),
                "size": record["size"],
                "sha256": record["sha256"],
            }
        )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(
            {
                "schema": driver.MANIFEST_SCHEMA,
                "experiment": "unit-test",
                "files": files,
            }
        ),
        encoding="utf-8",
    )


def _fake_git(repository: Path) -> dict[str, object]:
    return {
        "repository": str(repository.resolve()),
        "head": "a" * 40,
        "branch": "unit-test",
        "dirty": False,
        "status_sha256": hashlib.sha256(b"").hexdigest(),
        "status_entries": [],
        "tracked_diff_sha256": hashlib.sha256(b"").hexdigest(),
    }


def _gate() -> dict[str, object]:
    return {
        "allowed": True,
        "concurrent_existing_pea": False,
        "active_pea_processes": [],
        "resource_snapshot": {
            "captured_utc": "2026-08-31T00:00:00+00:00",
            "mem_available_bytes": 16 * 1024**3,
            "swap_free_bytes": 16 * 1024**3,
        },
    }


def _write_fake_pea(proc: Path, pid: int, starttime_ticks: int = 1000) -> Path:
    process = proc / str(pid)
    process.mkdir()
    (process / "comm").write_text("pea\n", encoding="utf-8")
    (process / "cmdline").write_bytes(b"/path/pea\0-y\0config.yaml\0")
    stat_tail = ["S", *(["0"] * 18), str(starttime_ticks)]
    (process / "stat").write_text(
        f"{pid} (pea) {' '.join(stat_tail)}\n", encoding="utf-8"
    )
    return process


class DriverTests(unittest.TestCase):
    def test_windows_worktree_gitdir_is_mapped_for_wsl(self) -> None:
        repository = Path("/mnt/d/tec/ginan-main-stec-cov")
        mapped = driver._normalise_gitdir_path(
            "C:/Users/rx/Documents/GINAN/ginan/.git/worktrees/ginan-main-stec-cov",
            repository,
            platform_name="posix",
        )
        self.assertEqual(
            mapped,
            Path(
                "/mnt/c/Users/rx/Documents/GINAN/ginan/.git/worktrees/"
                "ginan-main-stec-cov"
            ),
        )

    def test_registered_manifest_has_exactly_23_valid_hashes(self) -> None:
        manifest_path = (
            Path(__file__).resolve().parents[1]
            / "Docs/stecCovarianceFeasibility/experiment_1_2024_input_manifest.json"
        )
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        self.assertEqual(manifest["schema"], driver.MANIFEST_SCHEMA)
        self.assertEqual(len(manifest["files"]), driver.EXPECTED_INPUT_COUNT)
        self.assertEqual(
            len({entry["relative_path"] for entry in manifest["files"]}),
            driver.EXPECTED_INPUT_COUNT,
        )
        for entry in manifest["files"]:
            self.assertRegex(entry["sha256"], r"^[0-9a-f]{64}$")
            self.assertGreater(entry["size"], 0)

    def test_manifest_verification_detects_changed_input(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            data_root = root / "data"
            manifest = root / "manifest.json"
            _make_manifest(manifest, data_root)
            snapshot = driver.verify_input_manifest(manifest, data_root)
            self.assertTrue(snapshot["all_match"])
            self.assertEqual(snapshot["file_count"], driver.EXPECTED_INPUT_COUNT)
            (data_root / "inputs/file_03.dat").write_bytes(b"changed\n")
            with self.assertRaises(driver.DriverError):
                driver.verify_input_manifest(manifest, data_root)

    def test_atomic_json_never_replaces_destination(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            destination = Path(temporary) / "receipt.json"
            first = driver.atomic_write_json_no_replace(destination, {"value": 1})
            self.assertEqual(first["sha256"], driver.record_file(destination)["sha256"])
            with self.assertRaises(FileExistsError):
                driver.atomic_write_json_no_replace(destination, {"value": 2})
            self.assertEqual(
                json.loads(destination.read_text(encoding="utf-8")), {"value": 1}
            )

    def test_freeze_binary_is_distinct_and_content_identical(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "bin/pea"
            destination = root / "batch/provenance/pea"
            _write(source, b"not-a-real-pea\n")
            source.chmod(0o755)
            frozen = driver._copy_binary_no_replace(source, destination)
            self.assertEqual(frozen["sha256"], driver.record_file(source)["sha256"])
            self.assertNotEqual(source.resolve(), destination.resolve())
            with self.assertRaises(FileExistsError):
                driver._copy_binary_no_replace(source, destination)

    def test_freeze_rejects_binary_replaced_after_plan(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "bin/pea"
            destination = root / "batch/provenance/pea"
            _write(source, b"planned binary\n")
            planned = driver.record_file(source)
            source.write_bytes(b"replaced binary\n")
            with self.assertRaises(driver.DriverError):
                driver.freeze_planned_binary(planned, source, destination)
            self.assertFalse(destination.exists())

    def test_recursive_config_closure_detects_included_yaml_drift(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            repository = Path(temporary) / "repo"
            leaf = repository / "config/leaf.yaml"
            base = repository / "config/base.yaml"
            overlay = repository / "config/overlay.yaml"
            _write(leaf, b"outputs:\n  metadata:\n    config_description: leaf\n")
            _write(
                base,
                b"inputs:\n  include_yamls:\n    - config/leaf.yaml\n",
            )
            _write(
                overlay,
                b"inputs:\n  include_yamls:\n    - config/base.yaml\n",
            )
            closure = driver.config_include_closure(overlay, repository)
            self.assertEqual(closure["file_count"], 3)
            self.assertEqual(
                [Path(record["path"]).name for record in closure["files"]],
                ["overlay.yaml", "base.yaml", "leaf.yaml"],
            )
            driver.verify_config_closure(closure)
            base.write_text("changed: true\n", encoding="utf-8")
            with self.assertRaises(driver.DriverError):
                driver.verify_config_closure(closure)

    def test_default_concurrency_gate_rejects_pea_and_override_is_strict(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            proc = root / "proc"
            meminfo = proc / "meminfo"
            proc.mkdir()
            meminfo.write_text(
                "MemAvailable: 4194304 kB\nSwapFree: 9437184 kB\n",
                encoding="utf-8",
            )
            _write_fake_pea(proc, 35124)
            with self.assertRaises(driver.DriverError):
                driver.concurrency_gate(None, proc_root=proc, meminfo_path=meminfo)
            allowed = driver.concurrency_gate(
                35124, concurrent_nice=5, proc_root=proc, meminfo_path=meminfo
            )
            self.assertTrue(allowed["concurrent_existing_pea"])
            self.assertEqual(allowed["execution_priority"]["nice_adjustment"], 5)
            self.assertEqual(allowed["execution_priority"]["ionice_class_number"], 3)
            _write_fake_pea(proc, 35125, 2000)
            with self.assertRaises(driver.DriverError):
                driver.concurrency_gate(
                    35124, proc_root=proc, meminfo_path=meminfo
                )

    def test_concurrent_override_rejects_when_named_pea_is_not_active(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            proc = Path(temporary) / "proc"
            proc.mkdir()
            meminfo = proc / "meminfo"
            meminfo.write_text(
                "MemAvailable: 4194304 kB\nSwapFree: 9437184 kB\n",
                encoding="utf-8",
            )
            with self.assertRaises(driver.DriverError):
                driver.concurrency_gate(35124, proc_root=proc, meminfo_path=meminfo)

    def test_concurrent_execution_policy_freezes_exact_wrapper_command(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            wrappers = {}
            for name in ("time", "nice", "ionice"):
                path = root / name
                _write(path, f"fake {name}\n".encode("ascii"))
                path.chmod(0o755)
                wrappers[name] = path
            gate = {
                "allowed": True,
                "concurrent_existing_pea": True,
                "allowed_existing_pea_pid": 35124,
                "active_pea_processes": [
                    {
                        "pid": 35124,
                        "starttime_ticks": 1000,
                        "comm": "pea",
                        "command": "/path/pea -y config.yaml",
                    }
                ],
                "resource_snapshot": {
                    "captured_utc": "2026-08-31T00:00:00+00:00",
                    "mem_available_bytes": 4 * 1024**3,
                    "swap_free_bytes": 9 * 1024**3,
                },
                "execution_priority": {
                    "nice_adjustment": 5,
                    "ionice_class_name": "idle",
                    "ionice_class_number": 3,
                },
            }
            policy = driver.build_execution_policy(
                gate,
                concurrent_nice=5,
                time_binary=wrappers["time"],
                nice_binary=wrappers["nice"],
                ionice_binary=wrappers["ionice"],
            )
            driver.verify_execution_policy(policy)
            driver.verify_gate_matches_policy(gate, policy)
            self.assertEqual(policy["mode"], "concurrent_low_priority")
            self.assertEqual(
                policy["command_prefix"],
                [
                    str(wrappers["time"].absolute()),
                    "-v",
                    "--",
                    str(wrappers["nice"].absolute()),
                    "-n",
                    "5",
                    str(wrappers["ionice"].absolute()),
                    "-c",
                    "3",
                ],
            )
            malformed_gates = []
            malformed_gates.append({**gate, "allowed": False})
            malformed_gates.append(
                {
                    **gate,
                    "active_pea_processes": [
                        *gate["active_pea_processes"],
                        {
                            "pid": 35125,
                            "starttime_ticks": 2000,
                            "comm": "pea",
                            "command": "/path/pea -y other.yaml",
                        },
                    ],
                }
            )
            malformed_gates.append(
                {
                    **gate,
                    "resource_snapshot": {
                        **gate["resource_snapshot"],
                        "mem_available_bytes": 2 * 1024**3,
                    },
                }
            )
            malformed_gates.append(
                {
                    **gate,
                    "resource_snapshot": {
                        **gate["resource_snapshot"],
                        "swap_free_bytes": 7 * 1024**3,
                    },
                }
            )
            for malformed in malformed_gates:
                with self.subTest(malformed=malformed):
                    with self.assertRaises(driver.DriverError):
                        driver.verify_gate_matches_policy(malformed, policy)
            for invalid_nice in (-20, 20, True):
                malformed_policy = json.loads(json.dumps(policy))
                malformed_policy["priority"]["nice_adjustment"] = invalid_nice
                with self.subTest(invalid_nice=invalid_nice):
                    with self.assertRaises(driver.DriverError):
                        driver.verify_execution_policy(malformed_policy)
            reused_pid_gate = json.loads(json.dumps(gate))
            reused_pid_gate["active_pea_processes"][0]["starttime_ticks"] = 2000
            driver.verify_gate_matches_policy(reused_pid_gate, policy)
            with self.assertRaises(driver.DriverError):
                driver.verify_allowed_process_identity(gate, reused_pid_gate)

    def test_effective_environment_rejects_dynamic_loader_injection(self) -> None:
        with patch.dict(os.environ, {"LD_PRELOAD": "/tmp/injected.so"}, clear=False):
            with self.assertRaises(driver.DriverError):
                driver.build_effective_environment("/fake/boost")

    def test_concurrent_override_rejects_low_memory_or_swap(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            proc = Path(temporary) / "proc"
            proc.mkdir()
            _write_fake_pea(proc, 7)
            meminfo = proc / "meminfo"
            meminfo.write_text(
                "MemAvailable: 2097152 kB\nSwapFree: 9437184 kB\n",
                encoding="utf-8",
            )
            with self.assertRaises(driver.DriverError):
                driver.concurrency_gate(7, proc_root=proc, meminfo_path=meminfo)
            meminfo.write_text(
                "MemAvailable: 4194304 kB\nSwapFree: 7340032 kB\n",
                encoding="utf-8",
            )
            with self.assertRaises(driver.DriverError):
                driver.concurrency_gate(7, proc_root=proc, meminfo_path=meminfo)

    def _dry_run_fixture(self, root: Path) -> argparse.Namespace:
        repository = root / "repo"
        data_root = root / "data"
        output_parent = root / "outputs"
        output_parent.mkdir(parents=True)
        for spec in driver.RUN_SPECS:
            _write(repository / spec.config, f"# {spec.condition}\n".encode())
        _write(repository / "scripts/validate_stec_covariance.py", b"# validator\n")
        _write(
            repository / "scripts/validate_stec_satellite_difference_covariance.py",
            b"# validator\n",
        )
        _write(
            repository / "scripts/compare_pppar_integer_constraints.py",
            b"# comparator\n",
        )
        pea = repository / "bin/pea"
        _write(pea, b"fake binary\n")
        pea.chmod(0o755)
        wrapper_root = root / "wrappers"
        time_binary = wrapper_root / "time"
        nice_binary = wrapper_root / "nice"
        ionice_binary = wrapper_root / "ionice"
        for wrapper in (time_binary, nice_binary, ionice_binary):
            _write(wrapper, f"fake {wrapper.name}\n".encode("ascii"))
            wrapper.chmod(0o755)
        manifest = root / "manifest.json"
        _make_manifest(manifest, data_root)
        return argparse.Namespace(
            repository=repository,
            data_root=data_root,
            output_parent=output_parent,
            input_manifest=manifest,
            pea_binary=pea,
            time_binary=time_binary,
            nice_binary=nice_binary,
            ionice_binary=ionice_binary,
            ld_library_path="/fake/boost",
            batch_id="unit_test_batch",
            allow_existing_pea_pid=None,
            concurrent_nice=19,
            dry_run=True,
        )

    def test_dry_run_plans_five_serial_runs_without_starting_pea(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            args = self._dry_run_fixture(Path(temporary))
            with patch.object(driver, "concurrency_gate", return_value=_gate()), patch.object(
                driver, "git_identity", side_effect=_fake_git
            ), patch.object(
                driver,
                "execute_condition",
                side_effect=AssertionError("dry-run attempted to start PEA"),
            ):
                result = driver.run_batch(args)
            self.assertEqual(result["status"], "DRY_RUN")
            self.assertFalse(result["pea_started"])
            self.assertEqual(
                [run["condition"] for run in result["plan"]["runs"]],
                [spec.condition for spec in driver.RUN_SPECS],
            )
            self.assertEqual(
                [run["epochs"] for run in result["plan"]["runs"]],
                [480, 480, 480, 480, 420],
            )
            self.assertFalse(Path(result["plan"]["batch_root"]).exists())
            self.assertEqual(result["plan"]["execution_policy"]["mode"], "normal")
            self.assertEqual(
                result["plan"]["execution_policy"]["command_prefix"],
                [str(args.time_binary.absolute()), "-v", "--"],
            )

    def test_existing_batch_directory_is_rejected_before_pea(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            args = self._dry_run_fixture(Path(temporary))
            args.dry_run = False
            batch_root = args.output_parent / args.batch_id
            batch_root.mkdir()
            with patch.object(driver, "concurrency_gate", return_value=_gate()), patch.object(
                driver, "git_identity", side_effect=_fake_git
            ), patch.object(
                driver,
                "execute_condition",
                side_effect=AssertionError("collision attempted to start PEA"),
            ):
                with self.assertRaises(driver.DriverError):
                    driver.run_batch(args)

    def test_audit_rehashes_receipts_without_consulting_live_pea_state(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            run_args = self._dry_run_fixture(root)
            policy = driver.build_execution_policy(
                _gate(),
                concurrent_nice=run_args.concurrent_nice,
                time_binary=run_args.time_binary,
                nice_binary=run_args.nice_binary,
                ionice_binary=run_args.ionice_binary,
            )
            with patch.object(driver, "git_identity", side_effect=_fake_git):
                plan = driver.build_plan(
                    repository=run_args.repository,
                    data_root=run_args.data_root,
                    output_parent=run_args.output_parent,
                    manifest_path=run_args.input_manifest,
                    pea_binary=run_args.pea_binary,
                    batch_id=run_args.batch_id,
                    allow_existing_pea_pid=None,
                    concurrency=_gate(),
                    execution_policy=policy,
                    ld_library_path=run_args.ld_library_path,
                )
            batch_root = Path(plan["batch_root"])
            (batch_root / "provenance").mkdir(parents=True)
            (batch_root / "audit").mkdir()
            plan_record = driver.atomic_write_json_no_replace(
                batch_root / "batch.plan.json", plan
            )
            frozen = driver.freeze_planned_binary(
                plan["build_binary"],
                run_args.pea_binary,
                Path(plan["frozen_binary"]),
            )
            run_results: list[dict[str, object]] = []
            registered: list[dict[str, object]] = []
            for planned_run in plan["runs"]:
                output_root = Path(planned_run["output_root"])
                output_root.mkdir()
                exit_record = driver.atomic_write_bytes_no_replace(
                    output_root / "run.exit_code", b"0\n"
                )
                condition = str(planned_run["condition"])
                trace_audit = {
                    "resolved_epoch_count": 1 if condition in {"primary", "restart_0030"} else 0,
                    "complete_rank_epoch_count": 0,
                    "maximum_resolved_combination_count": 2 if condition == "primary" else 0,
                    "maximum_integer_candidate_count": 2 if condition == "primary" else 0,
                    "multirow_runtime_coverage": True,
                }
                run_receipt = {
                    "schema": driver.RUN_RECEIPT_SCHEMA,
                    "status": "SUCCESS",
                    "condition": condition,
                    "config_include_closure": planned_run["config_include_closure"],
                    "pea": {
                        "binary": frozen,
                        "argv": planned_run["argv"],
                        "environment": plan["environment"],
                    },
                    "planned_execution_policy": policy,
                    "concurrency": _gate(),
                    "execution": {
                        "returncode": 0,
                        "command": [*policy["command_prefix"], *planned_run["argv"]],
                    },
                    "validators": {
                        name: {
                            "returncode": 0,
                            "command": [
                                plan["static_files"]["python_interpreter"]["path"],
                                plan["static_files"][tool_role]["path"],
                            ],
                            "tool_role": tool_role,
                            "tool": plan["static_files"][tool_role],
                            "python_interpreter": plan["static_files"][
                                "python_interpreter"
                            ],
                        }
                        for name, tool_role in {
                            "forward_raw": "raw_validator",
                            "forward_sd": "sd_validator",
                            "smoothed_raw": "raw_validator",
                            "smoothed_sd": "sd_validator",
                        }.items()
                    },
                    "validation": {"all_valid": True, "trace_audit": trace_audit},
                    "artifacts": {"run.exit_code": exit_record},
                }
                receipt_record = driver.atomic_write_json_no_replace(
                    output_root / "run.receipt.json", run_receipt
                )
                registered.append(
                    {
                        "condition": condition,
                        "status": "SUCCESS",
                        "receipt": receipt_record,
                    }
                )
                run_results.append({"payload": run_receipt})
            comparison_report = driver.atomic_write_json_no_replace(
                batch_root / "audit/restart_integer_comparison.json", {"shared": 1}
            )
            comparison_log = driver.atomic_write_bytes_no_replace(
                batch_root / "audit/restart_integer_comparison.console.log", b"ok\n"
            )
            comparison = {
                "status": "CONSISTENT_SHARED_FIXES",
                "execution": {
                    "returncode": 0,
                    "command": [
                        plan["static_files"]["python_interpreter"]["path"],
                        plan["static_files"]["restart_comparator"]["path"],
                    ],
                },
                "tool": plan["static_files"]["restart_comparator"],
                "python_interpreter": plan["static_files"]["python_interpreter"],
                "artifacts": {"report": comparison_report, "console": comparison_log},
            }
            batch_receipt = {
                "schema": driver.BATCH_RECEIPT_SCHEMA,
                "status": "SUCCESS",
                "batch_id": plan["batch_id"],
                "plan": plan_record,
                "binary": frozen,
                "concurrency_preflight": plan["concurrency_preflight"],
                "execution_policy": policy,
                "runs": registered,
                "restart_comparison": comparison,
                "scientific_assessment": driver._batch_scientific_assessment(
                    run_results, comparison
                ),
                "claim_limit": "unit-test",
            }
            driver.atomic_write_json_no_replace(
                batch_root / "batch.receipt.json", batch_receipt
            )
            audit_args = argparse.Namespace(
                batch_root=batch_root,
                output=None,
                dry_run=False,
                write_audit=False,
            )
            with patch.object(
                driver,
                "active_pea_processes",
                side_effect=AssertionError("audit consulted live PEA state"),
            ):
                audited = driver.audit_batch(audit_args)
            self.assertEqual(audited["status"], "PASS")
            self.assertFalse(audited["audit_written"])
            self.assertEqual(audited["report"]["status"], "PASS")
            self.assertEqual(
                audited["report"]["scientific_status"],
                "DEFECT_FIX_RUNTIME_EVIDENCE_COMPLETE",
            )
            self.assertFalse((batch_root / "audit/batch.audit.json").exists())

            output = batch_root / "audit/explicit.audit.json"
            write_args = argparse.Namespace(
                batch_root=batch_root,
                output=output,
                dry_run=False,
                write_audit=True,
            )
            written = driver.audit_batch(write_args)
            self.assertEqual(written["status"], "PASS")
            self.assertTrue(written["audit_written"])
            self.assertTrue(output.is_file())

            invalid_args = argparse.Namespace(
                batch_root=batch_root,
                output=batch_root / "audit/implicit.audit.json",
                dry_run=False,
                write_audit=False,
            )
            with self.assertRaises(driver.DriverError):
                driver.audit_batch(invalid_args)

    def test_output_discovery_uses_actual_smoothed_sd_suffix(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output_root = Path(temporary)
            ionstec = output_root / "ionstec"
            ionstec.mkdir()
            for relative in (
                "Network-test.trace",
                "Network-test_smoothed.trace",
                "ionstec/test.STEC.COV",
                "ionstec/test.STEC_smoothed.COV",
                "ionstec/test.STEC.SD.COV",
                "ionstec/test.STEC.SD_smoothed.COV",
            ):
                _write(output_root / relative, b"test\n")
            outputs = driver.discover_outputs(output_root)
            self.assertEqual(outputs["sd_cov_smoothed"].name, "test.STEC.SD_smoothed.COV")

    def test_trace_audit_links_resolved_rows_to_diagonal_noise_evidence(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            covariance = root / "test.STEC.COV"
            covariance.write_text(
                "# GINAN_STEC_COVARIANCE_V3\n"
                "# META,fields\n"
                "META,2323,259200,OK,1,1,0,FILTER_POSTERIOR_AFTER_AR,1,4,2,1,2,0,2,1,"
                "LAMBDA_ALT,0.9999,3,RESOLVED_RATIO_ACCEPTED,2,2,0.9999,1,4,4\n",
                encoding="utf-8",
            )
            trace = root / "Network.trace"
            trace.write_text(
                "fixAndHoldAmbiguities: 2024-07-17 00:00:00.00\n"
                "PPP_AR INTEGER_FEEDBACK rows=2 integer_coordinates=2 "
                "original_ambiguities=4 status=Z_TIMES_D_MAPPED\n"
                "PPP_AR PSEUDOOBS_DESIGN rows=2 original_ambiguities=4 "
                "z_rows=2 z_columns=4 status=MATCHES_Z_TIMES_D\n"
                "PPP_AR PSEUDOOBS_NOISE rows=2 variance=1.000e-08 "
                "status=INDEPENDENT_DIAGONAL\n"
                "PPP_AR PSEUDOOBS_SUBMISSION rows=2 "
                "status=FILTER_CALL_RETURNED_SUBMITTED_UNVERIFIED\n",
                encoding="utf-8",
            )
            report = driver.audit_trace(
                condition="primary",
                expected_epochs=1,
                receiver_sd_expected=True,
                raw_covariance=covariance,
                trace=trace,
            )
            self.assertTrue(report["all_valid"])
            self.assertEqual(report["resolved_epoch_count"], 1)
            self.assertEqual(report["complete_rank_epoch_count"], 1)
            self.assertEqual(report["maximum_resolved_combination_count"], 2)
            self.assertTrue(report["multirow_runtime_coverage"])
            self.assertEqual(report["maximum_integer_candidate_count"], 2)

            trace_without_submission = root / "Network-no-submission.trace"
            trace_without_submission.write_text(
                trace.read_text(encoding="utf-8").replace(
                    "PPP_AR PSEUDOOBS_SUBMISSION rows=2 "
                    "status=FILTER_CALL_RETURNED_SUBMITTED_UNVERIFIED\n",
                    "",
                ),
                encoding="utf-8",
            )
            missing_submission = driver.audit_trace(
                condition="primary",
                expected_epochs=1,
                receiver_sd_expected=True,
                raw_covariance=covariance,
                trace=trace_without_submission,
            )
            self.assertFalse(missing_submission["all_valid"])
            self.assertFalse(
                missing_submission["checks"][
                    "integer_feedback_design_noise_submission"
                ]
            )

            text = covariance.read_text(encoding="utf-8").replace(
                ",2,2,0.9999,1,4,4", ",2,3,0.9999,1,4,4"
            )
            covariance.write_text(text, encoding="utf-8")
            failed = driver.audit_trace(
                condition="primary",
                expected_epochs=1,
                receiver_sd_expected=True,
                raw_covariance=covariance,
                trace=trace,
            )
            self.assertFalse(failed["all_valid"])
            self.assertFalse(failed["checks"]["lambda_candidate_limit"])

            undifferenced_covariance = root / "undifferenced.STEC.COV"
            undifferenced_covariance.write_text(
                covariance.read_text(encoding="utf-8")
                .replace(",2,3,0.9999,1,4,4", ",2,2,0.9999,1,4,4")
                .replace(",4,2,1,2,0,2,1,", ",4,2,0,2,0,2,1,"),
                encoding="utf-8",
            )
            stopped = driver.audit_trace(
                condition="undifferenced",
                expected_epochs=1,
                receiver_sd_expected=False,
                raw_covariance=undifferenced_covariance,
                trace=trace,
            )
            self.assertFalse(stopped["all_valid"])
            self.assertFalse(stopped["checks"]["undifferentiated_no_complete_rank"])

    def test_phase_osb_disabled_submitted_feedback_is_scientifically_inconclusive(self) -> None:
        def result(
            condition: str, resolved: int, complete: int = 0, maximum: int = 0
        ) -> dict[str, object]:
            return {
                "payload": {
                    "condition": condition,
                    "validation": {
                        "trace_audit": {
                            "resolved_epoch_count": resolved,
                            "complete_rank_epoch_count": complete,
                            "maximum_resolved_combination_count": maximum,
                        }
                    },
                }
            }

        assessment = driver._batch_scientific_assessment(
            [
                result("primary", 20, maximum=4),
                result("phase_osb_disabled", 18),
                result("undifferenced", 0),
            ],
            {"status": "CONSISTENT_SHARED_FIXES"},
        )
        self.assertEqual(assessment["status"], "INCONCLUSIVE")
        self.assertIn(
            "PHASE_OSB_DISABLED_SUBMITTED_FEEDBACK_ROWS", assessment["flags"]
        )

        single_row = driver._batch_scientific_assessment(
            [
                result("primary", 20, maximum=1),
                result("phase_osb_disabled", 0),
                result("undifferenced", 0),
            ],
            {"status": "CONSISTENT_SHARED_FIXES"},
        )
        self.assertEqual(single_row["status"], "INCONCLUSIVE")
        self.assertIn(
            "PRIMARY_MULTIROW_RUNTIME_COVERAGE_NOT_ESTABLISHED",
            single_row["flags"],
        )


if __name__ == "__main__":
    unittest.main()
