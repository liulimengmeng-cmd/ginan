"""Safety/preservation contracts; synthetic fixtures are not solver evidence."""
import csv
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import run_static_checkpoint_replay as replay
import run_stages as stages


class StaticCheckpointContracts(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name).resolve()
        stages.write_json(self.root / "OWNERSHIP.json", {"schema": replay.OWNERSHIP, "root": str(self.root)})
        for directory in ("outputs", "checkpoints"):
            (self.root / "work" / directory).mkdir(parents=True)

    def atomic(self, root, epoch):
        name = str(int(replay.epoch_seconds(epoch)))
        bundle = root / "zhang_internal_products.csv.epochs" / name
        bundle.mkdir(parents=True)
        (bundle / "products.csv").write_text("epoch,phase\n" + epoch + ",1.12345678901234567\n")
        (bundle / "covariance.csv").write_text("covariance_m2\n0.0012345678901234567\n")
        (bundle / "certificate.txt").write_text("ar=0\n")
        stages.write_json(bundle / "COMMIT.json", {"epoch": epoch, "ar": False})
        names = ("products.csv", "covariance.csv", "certificate.txt", "COMMIT.json")
        (bundle / "SHA256SUMS").write_text("".join(stages.file_sha(bundle / name) + "  " + name + "\n" for name in names))
        return bundle

    def checkpoint(self, root, epoch):
        directory = root / ("checkpoint_" + str(int(replay.epoch_seconds(epoch))))
        directory.mkdir()
        stages.write_json(directory / "checkpoint_manifest.json", {"epoch": epoch})
        (directory / "checkpoint.bundle").write_bytes(b"UNIT_FIXTURE_ONLY_NOT_A_REAL_CHECKPOINT")
        return directory

    def test_owned_archive_keeps_all_bytes_and_rejects_root_escape(self):
        (self.root / "work/outputs/data.txt").write_text("do not lose existing output\n")
        self.checkpoint(self.root / "work/checkpoints", replay.FIRST)
        before = replay.inventory(self.root / "work/outputs")
        archived = replay.archive_owned_outputs(self.root, "interrupted_before_resume_cut")
        self.assertEqual(replay.inventory(archived / "outputs"), before)
        self.assertFalse((self.root / "work/outputs").exists())
        self.assertTrue((archived / "ARCHIVE_SHA256.json").is_file())
        with self.assertRaisesRegex(RuntimeError, "Invalid or existing"):
            replay.archive_owned_outputs(self.root, "../not_owned")
        (self.root / "OWNERSHIP.json").write_text("{}")
        with self.assertRaisesRegex(RuntimeError, "owned test root"):
            replay.owner_guard(self.root)

    def test_late_stop_cut_preserves_original_future_epochs_without_rounding(self):
        outputs, checkpoints = self.root / "work/outputs", self.root / "work/checkpoints"
        for epoch in ("2024-07-17 00:00:00", replay.FIRST, replay.RESUME, replay.LAST):
            self.atomic(outputs, epoch)
        self.checkpoint(checkpoints, replay.FIRST)
        self.checkpoint(checkpoints, replay.LAST)
        flat = outputs / "zhang_internal_products.csv"
        original = "gpst_seconds,phase_m\r\n" + "".join(
            f"{replay.epoch_seconds(epoch):.1f},1.12345678901234567\r\n" for epoch in
            ("2024-07-17 00:00:00", replay.FIRST, replay.RESUME, replay.LAST))
        flat.write_bytes(original.encode())
        sha = stages.file_sha(flat)
        archived = replay.archive_owned_outputs(self.root, "interrupted_before_resume_cut")
        source = replay.restore_owned_prefix(self.root, archived)
        self.assertTrue(source.is_dir())
        self.assertEqual(stages.file_sha(archived / "outputs" / flat.name), sha)
        retained = (self.root / "work/outputs" / flat.name).read_bytes()
        self.assertEqual(retained, "".join(original.splitlines(keepends=True)[:3]).encode())
        active_epochs = list((self.root / "work/outputs/zhang_internal_products.csv.epochs").iterdir())
        self.assertEqual(len(active_epochs), 2)
        self.assertEqual(len(list((archived / "outputs/zhang_internal_products.csv.epochs").iterdir())), 4)
        self.assertEqual(len(list((self.root / "work/checkpoints").iterdir())), 1)
        audit = json.loads((self.root / "resume_prefix_audit.json").read_text())
        self.assertFalse(audit["numeric_values_reformatted"])
        self.assertEqual(audit["csv"][flat.name]["archived_future_rows_not_copied"], 2)

    def test_atomic_hash_failure_or_ar_output_is_rejected(self):
        bundle = self.atomic(self.root / "work/outputs", replay.FIRST)
        replay.verify_atomic_bundle(bundle)
        (bundle / "covariance.csv").write_text("changed")
        with self.assertRaisesRegex(RuntimeError, "hash mismatch"):
            replay.verify_atomic_bundle(bundle)

    def test_all_five_epoch_paired_products_and_full_q_must_match(self):
        left, right = self.root / "left", self.root / "right"
        epochs = ["2024-07-17 00:00:00", replay.FIRST, replay.RESUME, "2024-07-17 00:01:30", replay.LAST]
        for epoch in epochs:
            self.atomic(left, epoch)
            self.atomic(right, epoch)
        self.assertTrue(replay.final_atomic_comparison(left, right)["passed"])
        extra = self.atomic(right, "2024-07-17 00:02:30")
        with self.assertRaisesRegex(RuntimeError, "all five"):
            replay.final_atomic_comparison(left, right)

    def test_actual_sigint130_required_not_successful_premature_exit(self):
        fake = type("Child", (), {"pid": 999, "returncode": 0, "poll": lambda self: 0, "wait": lambda self: 0})()
        plan = {"run_root": str(self.root), "command": ["not_executed"], "timeout_seconds": 1}
        with patch.object(replay, "verify_plan"), patch.object(replay.subprocess, "Popen", return_value=fake):
            with self.assertRaisesRegex(RuntimeError, "real SIGINT"):
                replay.run_phase(plan, {}, "B_interrupted", True)

    def test_staging_checkpoint_is_not_published_until_atomic_rename(self):
        checkpoints = self.root / "work/checkpoints"
        final = self.checkpoint(checkpoints, replay.FIRST)
        staging = final.with_name(final.name + ".tmp.123456789")
        final.rename(staging)
        self.assertIsNone(replay.checkpoint_at(checkpoints, replay.FIRST))
        staging.rename(final)
        self.assertEqual(replay.checkpoint_at(checkpoints, replay.FIRST), final)
        # A second complete staging copy must not be reported as ambiguity.
        staging.mkdir()
        stages.write_json(staging / "checkpoint_manifest.json", {"epoch": replay.FIRST})
        (staging / "checkpoint.bundle").write_bytes(b"UNPUBLISHED_STAGING_FIXTURE")
        self.assertEqual(replay.checkpoint_at(checkpoints, replay.FIRST), final)

    def test_partial_staging_manifest_cannot_authorize_interrupt(self):
        staging = self.root / "work/checkpoints/checkpoint.tmp.98765"
        staging.mkdir()
        (staging / "checkpoint_manifest.json").write_text("{incomplete")
        (staging / "checkpoint.bundle").write_bytes(b"UNPUBLISHED_STAGING_FIXTURE")
        self.assertIsNone(replay.checkpoint_at(staging.parent, replay.FIRST))
        polls = iter((None, 0))
        fake = type("Child", (), {"pid": 999, "returncode": 0,
                    "poll": lambda self: next(polls), "wait": lambda self: 0})()
        plan = {"run_root": str(self.root), "command": ["not_executed"], "timeout_seconds": 1}
        with patch.object(replay, "verify_plan"), \
             patch.object(replay.subprocess, "Popen", return_value=fake), \
             patch.object(replay.os, "killpg", create=True) as kill:
            with self.assertRaisesRegex(RuntimeError, "real SIGINT"):
                replay.run_phase(plan, {}, "B_unpublished_staging", True)
        kill.assert_not_called()

    def test_missing_controlled_build_gate_never_launches(self):
        path = self.root / "approval.json"
        path.write_text(json.dumps({"build_passed": True, "mandatory_tests_passed": True}))
        with patch.object(replay.sys, "platform", "linux"), patch.object(replay.subprocess, "Popen") as popen:
            with self.assertRaisesRegex(RuntimeError, "explicit build/test approval"):
                replay.execute({"run_root": str(self.root)}, {}, path)
        popen.assert_not_called()

    def test_preparation_five_epochs_static_coordinates_and_no_launch(self):
        selected = [f"A{i:03d}" for i in range(12)]
        audit = {"selected_complete": True, "missing_files": [],
            "selected_12_station_candidates": selected,
            "station_headers": [{"filename": station + ".rnx"} for station in selected],
            "product_membership": {"sp3_satellites": ["G02", "E02", "C19"],
                                   "clk_satellites": ["G02", "E02", "C19"]},
            "whitelist": [{"prn": "C19"}], "reference_stations": {"G": "MATE", "E": "MATE", "C": "MAL2"}}
        audit_path = self.root / "audit.json"
        audit_path.write_text(json.dumps(audit))
        coordinate = self.root / "coordinate.json"
        coordinate.write_text("{}")
        review = self.root / "review.json"
        review.write_text("{}")
        binary = self.root / "fake_pea"
        binary.write_text("NOT_EXECUTABLE_FIXTURE")
        pinned = {"observations": [{"station": station, "path": "/fake/" + station + ".rnx"} for station in selected],
                  "navigation": [{"path": "/fake/nav.rnx"}]}
        target = self.root / "new_controlled_replay"
        with patch.object(replay.network, "coordinate_approval", return_value=({"mode": "STATIC_ESTIMATED"}, {})), \
             patch.object(stages, "pin_stage_inputs", return_value=pinned), \
             patch.object(replay.subprocess, "Popen") as popen, \
             patch.object(stages, "configuration_preflight") as parser:
            plan, env = replay.prepare(binary, audit_path, target, coordinate, review, self.root / "future_reader")
            popen.assert_not_called()
            parser.assert_not_called()
            self.assertEqual(plan["expected_epoch_count"], 5)
            self.assertEqual(plan["capture_epochs"], [replay.FIRST, replay.LAST])
            self.assertEqual(plan["expected_checkpoint_mode"], "E29_MULTI_GNSS_FLOAT_V2_GPS_GAL_BDS")
            self.assertEqual(plan["state"], "PREPARED_NOT_RUN")
            self.assertFalse(plan["ar_execution"])
            self.assertFalse(plan["scientific_acceptance"])
            self.assertEqual(env["ZHANG_MULTI_GNSS_FLOAT_CHECKPOINT"], "1")
            self.assertNotIn("ZHANG_E29_RESTORE_BUNDLE", env)
            overlay = json.loads(Path(plan["configs"][-2]).read_text())
            self.assertEqual(overlay["processing_options"]["epoch_control"]["end_epoch"], replay.LAST)
            self.assertTrue(overlay["processing_options"]["gnss_general"]["zhang_pppar"]["deterministic_checkpoint"])


if __name__ == "__main__":
    unittest.main()
