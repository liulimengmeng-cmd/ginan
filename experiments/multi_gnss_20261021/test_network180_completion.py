"""Minimum completion audit contracts, not a scientific acceptance fixture."""
import csv
from itertools import combinations_with_replacement
import json
from pathlib import Path
import tempfile
import unittest

import audit_network180_float_completion as audit
import run_stages as stages


class CompletionContracts(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def fixture(self):
        path = self.root / "atomic_epoch"
        path.mkdir()
        epoch = "2024-07-17 00:00:00"
        seconds = audit.replay.epoch_seconds(epoch)
        fields = ["gpst_seconds", "solution", "satellite", "observable", *audit.FINITE_PRODUCT_FIELDS,
                  "numeric_valid", "ar_valid", "pppar_usable", "previous_epoch_gpst_seconds"]
        products, labels = [], set()
        for system, codes in audit.SIGNALS.items():
            satellite = system + "02"
            for code in sorted(codes):
                row = dict.fromkeys(fields, "1")
                row.update(gpst_seconds=str(seconds), solution="FLOAT", satellite=satellite, observable=code,
                           numeric_valid="1", ar_valid="0", pppar_usable="0", previous_epoch_gpst_seconds="nan")
                products.append(row)
                labels.update(((satellite, "CLOCK", "NONE"), (satellite, "PHASE", code)))
        with (path / "products.csv").open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fields)
            writer.writeheader()
            writer.writerows(products)
        cov_fields = ["gpst_seconds", "solution", "row_satellite", "row_parameter", "row_observable",
            "column_satellite", "column_parameter", "column_observable", "covariance_m2"]
        with (path / "covariance.csv").open("w", newline="") as stream:
            writer = csv.writer(stream)
            writer.writerow(cov_fields)
            for left, right in combinations_with_replacement(sorted(labels), 2):
                writer.writerow([seconds, "FLOAT", *left, *right, 1 if left == right else .05])
        (path / "certificate.txt").write_text("FLOAT; NO_INTEGER_AUTHORITY\n")
        stages.write_json(path / "COMMIT.json", {"epoch": epoch, "ar": False, "schema_version": 2,
            "systems": {system: {"fixed_phase_rows": 0, "any_pppar_usable": False} for system in ("GPS", "GAL", "BDS")},
            "products_bytes": (path / "products.csv").stat().st_size,
            "covariance_bytes": (path / "covariance.csv").stat().st_size})
        self.rehash(path)
        return path

    def rehash(self, path):
        (path / "SHA256SUMS").write_text("".join(stages.file_sha(path / name) + "  " + name + "\n"
            for name in ("products.csv", "covariance.csv", "certificate.txt", "COMMIT.json")))

    def test_241_schedule_includes_both_two_hour_boundary_epochs(self):
        values = audit.expected_epochs()
        self.assertEqual(len(values), 241)
        self.assertEqual(values[0], "2024-07-17 00:00:00")
        self.assertEqual(values[-2:], ["2024-07-17 01:59:30", "2024-07-17 02:00:00"])
        self.assertEqual(len(set(values)), 241)

    def test_finite_physical_full_cross_system_q_with_optional_nan_metadata(self):
        result = audit.audit_epoch(self.fixture())
        self.assertTrue(result["physical_moments_finite"])
        self.assertEqual(result["phase_rows_by_system"], {"G": 2, "E": 2, "C": 2})
        self.assertEqual(result["covariance_variable_count"], 9)
        self.assertEqual(result["full_upper_triangle_entries"], 45)

    def test_correctly_rehashed_truncated_q_is_not_complete(self):
        path = self.fixture()
        covariance = path / "covariance.csv"
        lines = covariance.read_bytes().splitlines(keepends=True)
        covariance.write_bytes(b"".join(lines[:-1]))
        manifest = json.loads((path / "COMMIT.json").read_text())
        manifest["covariance_bytes"] = covariance.stat().st_size
        stages.write_json(path / "COMMIT.json", manifest)
        self.rehash(path)
        with self.assertRaisesRegex(RuntimeError, "complete unique upper triangle"):
            audit.audit_epoch(path)

    def test_forged_ar_authority_or_nonfinite_physical_value_is_rejected(self):
        path = self.fixture()
        manifest = json.loads((path / "COMMIT.json").read_text())
        manifest["ar"] = True
        stages.write_json(path / "COMMIT.json", manifest)
        self.rehash(path)
        with self.assertRaisesRegex(RuntimeError, "AR output is forbidden"):
            audit.audit_epoch(path)
        manifest["ar"] = False
        with (path / "products.csv").open(newline="") as stream:
            reader = csv.DictReader(stream)
            fields, rows = reader.fieldnames, list(reader)
        rows[0]["phase_m"] = "nan"
        with (path / "products.csv").open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fields)
            writer.writeheader()
            writer.writerows(rows)
        manifest["products_bytes"] = (path / "products.csv").stat().st_size
        stages.write_json(path / "COMMIT.json", manifest)
        self.rehash(path)
        with self.assertRaisesRegex(RuntimeError, "Nonfinite physical product moment"):
            audit.audit_epoch(path)

    def test_job_exit_zero_does_not_replace_expected_epoch_contract(self):
        stages.write_json(self.root / "stage_plan.json", {"exit_code": 0, "expected_epoch_count": 2})
        with self.assertRaisesRegex(RuntimeError, "actual 241-epoch"):
            audit.run(self.root, self.root / "not_a_verifier", self.root / "audit")

    def test_generic_coordinate_reader_protocol_full_core_hash_and_mode(self):
        actual = {"schema": "CHECKPOINT_COORDINATE_INVENTORY_V1", "producer_sha256": "producer",
            "experiment_mode": audit.replay.MODE, "checkpoint_epoch": "2024-07-17 02:00:00",
            "full_x_finite": True, "full_P_finite": True, "state_dimension": 6,
            "P_entries": 36, "bundle_sha256": "bundle"}
        audit.validate_coordinate_inventory(actual, "producer", "2024-07-17 02:00:00", 6, "bundle")
        actual["P_entries"] = 35
        with self.assertRaisesRegex(RuntimeError, "finite/dimension/hash"):
            audit.validate_coordinate_inventory(actual, "producer", "2024-07-17 02:00:00", 6, "bundle")
        actual["P_entries"] = 36
        actual["experiment_mode"] = "E29_GPS_L1C_L2W_ZHANG_FULL_RANK"
        with self.assertRaisesRegex(RuntimeError, "schema/producer/mode"):
            audit.validate_coordinate_inventory(actual, "producer", "2024-07-17 02:00:00", 6, "bundle")


if __name__ == "__main__":
    unittest.main()
