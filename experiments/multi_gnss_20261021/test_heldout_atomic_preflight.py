"""Structural preflight fixtures only, not scientific integer certificates."""
import csv
import json
from pathlib import Path
import tempfile
import unittest

from check_heldout_ready import audit_atomic_bundles
from prepare_heldout_users import sha
from run_stages import epoch_seconds


class AtomicPreflightTests(unittest.TestCase):
    def fixture(self, base):
        columns = ["gpst_seconds", "solution", "satellite", "same_segment_continuous",
                   "requires_user_ambiguity_reset", "public_phase_segment_id", "joint_covariance_complete"]
        columns += ["test_column_" + str(i) for i in range(65 - len(columns))]
        for i in range(120):
            path = base / ("epoch_" + str(i))
            path.mkdir()
            epoch = epoch_seconds("02:00:00") + 30 * i
            minute, second = divmod(30 * i, 60)
            with (path / "products.csv").open("w", newline="") as stream:
                writer = csv.DictWriter(stream, fieldnames=columns)
                writer.writeheader()
                for sat in ("G02", "E11", "C19"):
                    writer.writerow(dict(gpst_seconds=epoch, solution="PRODUCT_FIXED", satellite=sat))
            (path / "covariance.csv").write_text("STRUCTURAL_TEST_NOT_REAL_COVARIANCE\n")
            (path / "certificate.txt").write_text("STRUCTURAL_TEST_NOT_AN_INTEGER_CERTIFICATE\n")
            metadata = dict(epoch=f"2024-07-17 02:{minute:02d}:{second:02d}", schema_version=2,
                            systems={s: {} for s in ("GPS", "GAL", "BDS")},
                            products_bytes=(path / "products.csv").stat().st_size,
                            covariance_bytes=(path / "covariance.csv").stat().st_size)
            (path / "COMMIT.json").write_text(json.dumps(metadata))
            names = ("COMMIT.json", "products.csv", "covariance.csv", "certificate.txt")
            (path / "SHA256SUMS").write_text("".join(sha(path / name) + "  " + name + "\n" for name in names))

    def test_complete_hash_structure_and_corruption(self):
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory)
            self.fixture(base)
            self.assertEqual(len(audit_atomic_bundles(base)), 64)
            (base / "epoch_7" / "covariance.csv").write_text("CHANGED\n")
            with self.assertRaisesRegex(RuntimeError, "HASH_MISMATCH"):
                audit_atomic_bundles(base)

    def test_pending_and_duplicate_are_not_authority(self):
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory)
            self.fixture(base)
            (base / "epoch_119").rename(base / "epoch_119.pending")
            with self.assertRaisesRegex(RuntimeError, "120_EPOCH"):
                audit_atomic_bundles(base)
            (base / "epoch_119.pending").rename(base / "epoch_119")
            import shutil
            shutil.copytree(base / "epoch_4", base / "duplicate_epoch")
            with self.assertRaisesRegex(RuntimeError, "AMBIGUOUS"):
                audit_atomic_bundles(base)

    def test_flat_file_is_not_atomic_authority(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "products.csv"
            path.write_text("no atomic proof\n")
            with self.assertRaisesRegex(RuntimeError, "NO_REAL_MULTI_SYSTEM_ATOMIC"):
                audit_atomic_bundles(path)


if __name__ == "__main__":
    unittest.main()
