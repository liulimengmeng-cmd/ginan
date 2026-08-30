from __future__ import annotations

import hashlib
import json
import os
import subprocess
import sys
import tempfile
import threading
import unittest
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parent))

from experiment0_2024_blind_provenance import (
    FILE_RECORD_FIELDS,
    InputChangedError,
    ProvenanceError,
    UnsafePathError,
    assert_input_snapshot_unchanged,
    assert_planned_outputs_absent,
    atomic_write_json_no_replace,
    begin_provenance,
    capture_input_snapshot,
    complete_provenance,
    git_identity,
    record_regular_file,
    runtime_identity,
)


class Experiment0BlindProvenanceTests(unittest.TestCase):
    def test_file_record_and_unchanged_snapshot(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            path = Path(temporary_directory) / "input.bin"
            payload = b"stable input\x00\xff"
            path.write_bytes(payload)

            record = record_regular_file(path)
            self.assertEqual(set(record), set(FILE_RECORD_FIELDS))
            self.assertEqual(record["path"], str(path.resolve()))
            self.assertEqual(record["resolved_path"], str(path.resolve()))
            self.assertEqual(record["size"], len(payload))
            self.assertEqual(record["sha256"], hashlib.sha256(payload).hexdigest())
            self.assertIsInstance(record["mtime_ns"], int)
            self.assertIsInstance(record["dev"], int)
            self.assertIsInstance(record["inode"], int)

            before = capture_input_snapshot({"quiet_model": path})
            after = assert_input_snapshot_unchanged(before)
            self.assertEqual(before, after)

    def test_changed_content_and_roles_are_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            path = Path(temporary_directory) / "input.txt"
            path.write_text("before", encoding="utf-8")
            before = capture_input_snapshot({"input": path})
            path.write_text("after!", encoding="utf-8")
            with self.assertRaises(InputChangedError):
                assert_input_snapshot_unchanged(before, {"input": path})
            with self.assertRaises(InputChangedError):
                assert_input_snapshot_unchanged(before, {"renamed": path})

    def test_missing_directory_and_symlink_inputs_are_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            with self.assertRaises(UnsafePathError):
                record_regular_file(root / "missing")
            with self.assertRaises(UnsafePathError):
                record_regular_file(root)

            target = root / "target.txt"
            link = root / "link.txt"
            target.write_text("target", encoding="utf-8")
            try:
                link.symlink_to(target)
            except OSError as exc:
                self.skipTest(f"symlink creation unavailable: {exc}")
            with self.assertRaises(UnsafePathError):
                record_regular_file(link)

    def test_planned_outputs_must_be_absent_and_distinct(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            first = root / "first.json"
            second = root / "second.json"
            planned = assert_planned_outputs_absent(
                {"report": first, "provenance": second}
            )
            self.assertEqual(planned["report"], str(first.resolve()))
            first.write_text("occupied", encoding="utf-8")
            with self.assertRaises(FileExistsError):
                assert_planned_outputs_absent([first])
            with self.assertRaises(ProvenanceError) as context:
                assert_planned_outputs_absent({"one": second, "two": second})
            self.assertIn("resolve to", str(context.exception))

    def test_atomic_json_never_replaces_existing_destination(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            output = Path(temporary_directory) / "result.json"
            payload = {"finite": 1.25, "items": [3, 2, 1]}
            with patch("os.replace", side_effect=AssertionError("must not be used")):
                record = atomic_write_json_no_replace(output, payload)
            self.assertEqual(json.loads(output.read_text(encoding="utf-8")), payload)
            self.assertEqual(
                record["sha256"], hashlib.sha256(output.read_bytes()).hexdigest()
            )

            original = output.read_bytes()
            with self.assertRaises(FileExistsError):
                atomic_write_json_no_replace(output, {"replacement": True})
            self.assertEqual(output.read_bytes(), original)
            self.assertEqual(list(output.parent.glob(f".{output.name}.*.tmp")), [])

    def test_post_publish_temporary_cleanup_failure_is_not_ambiguous(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            output = Path(temporary_directory) / "result.json"
            with patch.object(Path, "unlink", side_effect=PermissionError("injected")):
                record = atomic_write_json_no_replace(output, {"committed": True})
            self.assertEqual(
                json.loads(output.read_text(encoding="utf-8")), {"committed": True}
            )
            self.assertEqual(record["path"], str(output.resolve()))
            temporary_aliases = list(output.parent.glob(f".{output.name}.*.tmp"))
            self.assertEqual(len(temporary_aliases), 1)
            temporary_aliases[0].unlink()

    def test_concurrent_publish_has_exactly_one_winner(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            output = Path(temporary_directory) / "result.json"
            barrier = threading.Barrier(2)

            def publish(value: int) -> str:
                barrier.wait()
                try:
                    atomic_write_json_no_replace(output, {"writer": value})
                except FileExistsError:
                    return "exists"
                return "published"

            with ThreadPoolExecutor(max_workers=2) as executor:
                outcomes = list(executor.map(publish, (1, 2)))
            self.assertEqual(sorted(outcomes), ["exists", "published"])
            self.assertIn(
                json.loads(output.read_text(encoding="utf-8"))["writer"], (1, 2)
            )
            self.assertEqual(list(output.parent.glob(f".{output.name}.*.tmp")), [])

    def test_git_and_runtime_identity(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            repo = Path(temporary_directory) / "repo"
            repo.mkdir()
            self._git(repo, "init")
            self._git(repo, "config", "user.name", "Provenance Test")
            self._git(repo, "config", "user.email", "test@example.invalid")
            tracked = repo / "tracked.txt"
            tracked.write_text("committed\n", encoding="utf-8")
            self._git(repo, "add", "tracked.txt")
            self._git(repo, "commit", "-m", "fixture")

            identity = git_identity(repo)
            self.assertTrue(identity["tracked_clean"])
            self.assertEqual(len(identity["head"]), 40)
            self.assertTrue(identity["branch"])

            (repo / "untracked.txt").write_text("ignored\n", encoding="utf-8")
            self.assertTrue(git_identity(repo)["tracked_clean"])
            tracked.write_text("modified\n", encoding="utf-8")
            dirty = git_identity(repo)
            self.assertFalse(dirty["tracked_clean"])
            self.assertTrue(dirty["tracked_status_porcelain_v1"])

        runtime = runtime_identity()
        self.assertEqual(
            runtime["python"]["executable"], str(Path(sys.executable).resolve())
        )
        self.assertTrue(runtime["python"]["version"])
        self.assertTrue(runtime["numpy"]["version"])
        self.assertTrue(runtime["numpy"]["module_path"])

    def test_begin_and_complete_provenance(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            repo = root / "repo"
            repo.mkdir()
            self._git(repo, "init")
            self._git(repo, "config", "user.name", "Provenance Test")
            self._git(repo, "config", "user.email", "test@example.invalid")
            self._git(repo, "commit", "--allow-empty", "-m", "fixture")
            input_path = root / "input.dat"
            input_path.write_bytes(b"input")
            output_path = root / "planned.json"

            provenance = begin_provenance(
                repo, {"model": input_path}, {"report": output_path}
            )
            completed = complete_provenance(provenance)
            self.assertTrue(completed["inputs_unchanged"])
            self.assertEqual(completed["inputs_pre"], completed["inputs_post"])
            self.assertFalse(output_path.exists())

    @staticmethod
    def _git(repo: Path, *arguments: str) -> None:
        subprocess.run(
            ["git", "-C", str(repo), *arguments],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            encoding="utf-8",
        )


if __name__ == "__main__":
    unittest.main()
