"""Synthetic audit-wire contracts, not claims about a GNSS experiment result."""
import csv
import hashlib
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

import audit_stage_results as audit


class ResultAuditTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        (self.root / "outputs").mkdir()

    def plan(self, systems=("gps", "gal", "bds"), count=2, state="COMPLETED_REQUIRES_SCIENTIFIC_AUDIT", code=0):
        plan = {"stage": "E", "state": state, "specification": {"systems": list(systems), "mode": "AR"},
                "expected_epoch_count": count}
        if code is not None:
            plan["exit_code"] = code
        (self.root / "stage_plan.json").write_text(json.dumps(plan))

    def bundle(self, epoch=1405209600, systems=("G", "E", "C"), matrix=None, omit_pair=False,
               reverse_conflict=False, missing_endpoint=None, nonfinite=False, segment="S1"):
        directory = self.root / "outputs/zhang_internal_products.csv.epochs" / str(epoch)
        directory.mkdir(parents=True)
        parameters = []
        for system in systems:
            sat = system + ("19" if system == "C" else "01")
            parameters.extend([(sat, "CLOCK", "NONE"), *[(sat, "PHASE", c) for c in audit.SIGNALS[system]]])
        indices = {key: i for i, key in enumerate(parameters)}
        if matrix is None:
            matrix = np.eye(len(parameters)) * 2 + np.ones((len(parameters), len(parameters))) * .02
        products = []
        for satellite, kind, code in parameters:
            if kind != "PHASE" or (satellite, code) == missing_endpoint:
                continue
            phase, clock = indices[(satellite, kind, code)], indices[(satellite, "CLOCK", "NONE")]
            row = {"gpst_seconds": epoch, "solution": "PRODUCT_FIXED", "satellite": satellite, "observable": code,
                "clock_m": 1, "clock_sigma_m": float(np.sqrt(abs(matrix[clock, clock]))), "phase_m": .2,
                "phase_sigma_m": float(np.sqrt(abs(matrix[phase, phase]))),
                "clock_phase_covariance_m2": float(matrix[clock, phase]), "correction_m": .8,
                "correction_sigma_m": float(np.sqrt(abs(matrix[clock, clock] + matrix[phase, phase] - 2 * matrix[clock, phase]))),
                "discontinuity_counter": "0", "datum_version": "0", "public_phase_segment_id": segment,
                "reset_reason": "TEST_ONLY", "invalid_reason": "NONE"}
            for name in audit.FLAGS:
                row[name] = "0" if name == "requires_user_ambiguity_reset" else "1"
            products.append(row)
        if nonfinite:
            products[0]["phase_m"] = "nan"
        with (directory / "products.csv").open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(products[0]))
            writer.writeheader()
            writer.writerows(products)
        cov = []
        for i, left in enumerate(parameters):
            for j in range(i, len(parameters)):
                if omit_pair and i == 0 and j == len(parameters) - 1:
                    continue
                right = parameters[j]
                cov.append(dict(zip(audit.COV_COLUMNS,
                    [epoch, "PRODUCT_FIXED", *left, *right, float(matrix[i, j])], strict=True)))
        if reverse_conflict:
            left, right = parameters[0], parameters[1]
            cov.append(dict(zip(audit.COV_COLUMNS,
                [epoch, "PRODUCT_FIXED", *right, *left, float(matrix[0, 1] + 1)], strict=True)))
        with (directory / "covariance.csv").open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=audit.COV_COLUMNS)
            writer.writeheader()
            writer.writerows(cov)
        certificate = "".join("multi_system_certificate_begin system=" + audit.NAMES[s] + "\n" +
            "SYNTHETIC_TEST_NOT_INTEGER_PROOF\n" + "multi_system_certificate_end system=" + audit.NAMES[s] + "\n"
            for s in systems)
        (directory / "certificate.txt").write_bytes(certificate.encode())
        manifest = {"epoch": str(epoch), "route": "SYNTHETIC_TEST", "root": "TEST_ONLY", "ar": True,
            "products_bytes": (directory / "products.csv").stat().st_size,
            "covariance_bytes": (directory / "covariance.csv").stat().st_size,
            "schema_version": 2, "systems": {}}
        for system in systems:
            name = audit.NAMES[system]
            section = "multi_system_certificate_begin system=" + name + "\n" + "SYNTHETIC_TEST_NOT_INTEGER_PROOF\n" + \
                "multi_system_certificate_end system=" + name + "\n"
            manifest["systems"][name] = {"signals": list(audit.SIGNALS[system]),
                "any_pppar_usable": True, "fixed_phase_rows": sum(r["satellite"][0] == system for r in products),
                "pppar_usable_phase_rows": sum(r["satellite"][0] == system for r in products),
                "epoch_certificate_sha256": hashlib.sha256(certificate.encode()).hexdigest(),
                "system_certificate_sha256": hashlib.sha256(section.encode()).hexdigest()}
        (directory / "COMMIT.json").write_text(json.dumps(manifest))
        self.rehash(directory)
        return directory

    def rehash(self, directory):
        names = ("products.csv", "covariance.csv", "certificate.txt", "COMMIT.json")
        (directory / "SHA256SUMS").write_text("".join(audit.sha256(directory / name) + "  " + name + "\n" for name in names))

    def log(self, count=2, explicit_spp=False):
        rows = []
        for index in range(count):
            seconds = index * 30
            date = f"2024-07-17 00:{seconds // 60:02d}:{seconds % 60:02d}"
            rows.append(f"Processed epoch #{index + 1} - GPS time: 2323 {259200 + seconds} - {date} "
                        f"(data handling took 0.01, processing took 1, epoch took {index + 1}.0)")
            if explicit_spp:
                rows.append(f"SPP time={date} receiver=TEST status=SINGLE")
        (self.root / "runner.log").write_text("\n".join(rows))

    def run_audit(self, **options):
        return audit.StageAuditor(self.root, **options).run()

    def test_prepared_is_not_run_and_missing_metrics_unavailable(self):
        self.plan(state="PREPARED_NOT_RUN", code=None)
        report = self.run_audit()
        self.assertEqual(report["execution_classification"], "PREPARED_NOT_RUN")
        self.assertFalse(report["scientific_acceptance"])
        for key in ("processed_epochs", "spp_input_observations", "spp_single"):
            self.assertEqual(report["logs"][key]["status"], "UNAVAILABLE")

    def test_exit_zero_cannot_establish_scientific_or_output_success(self):
        self.plan()
        report = self.run_audit()
        self.assertEqual(report["execution_classification"], "COMPLETED_WITH_EVIDENCE_FAILURES")
        self.assertIn("COMPLETED_WITHOUT_PRODUCT_EPOCHS", report["evidence_issues"])
        self.assertFalse(report["scientific_acceptance"])

    def test_complete_three_system_full_covariance_and_hashes(self):
        self.plan()
        self.bundle()
        self.bundle(1405209630)
        self.log(explicit_spp=True)
        report = self.run_audit()
        self.assertEqual(report["execution_classification"], "COMPLETED_REQUIRES_SCIENTIFIC_AUDIT")
        self.assertFalse(report["scientific_acceptance"])
        self.assertEqual(report["atomic"]["hash_verified_bundles"], 2)
        self.assertEqual(report["products"]["epoch_count"], 2)
        q = report["products"]["epochs"][0]["solutions"]["PRODUCT_FIXED"]["full_covariance"]
        self.assertTrue(q["internally_consistent"])
        self.assertEqual(q["parameter_count"], 9)
        self.assertEqual(q["unique_pairs"], 45)
        self.assertEqual(q["missing_all_product_parameters"], [])
        self.assertTrue(q["psd"]["is_psd_with_tolerance"])
        self.assertTrue(all(v["exact_nonzero_entries"] > 0 for v in q["cross_system_blocks"].values()))
        self.assertEqual(report["logs"]["spp_single"]["explicit_status_records"]["SINGLE"], 2)
        json.dumps(report, allow_nan=False)

    def test_tamper_is_quarantined_not_counted_as_verified_product_epoch(self):
        self.plan(count=1)
        path = self.bundle()
        with (path / "products.csv").open("a") as stream:
            stream.write("TAMPER\n")
        report = self.run_audit()
        self.assertEqual(report["atomic"]["hash_verified_bundles"], 0)
        self.assertEqual(report["products"]["epoch_count"], 0)
        self.assertIn("ATOMIC_BUNDLE_INTEGRITY_FAILED", report["evidence_issues"])

    def test_missing_q_pair_not_filled_zero_for_psd(self):
        self.plan(count=1)
        self.bundle(omit_pair=True)
        report = self.run_audit()
        q = report["products"]["epochs"][0]["solutions"]["PRODUCT_FIXED"]["full_covariance"]
        self.assertFalse(q["complete_pair_coverage"])
        self.assertEqual(q["psd"]["status"], "UNAVAILABLE")

    def test_negative_eigenvalue_fails_full_q(self):
        self.plan(count=1)
        matrix = np.eye(9)
        matrix[0, 1] = matrix[1, 0] = 2
        self.bundle(matrix=matrix)
        report = self.run_audit()
        q = report["products"]["epochs"][0]["solutions"]["PRODUCT_FIXED"]["full_covariance"]
        self.assertFalse(q["psd"]["is_psd_with_tolerance"])
        self.assertLess(q["psd"]["minimum_eigenvalue_m2"], 0)

    def test_reverse_covariance_conflict_is_not_hidden_by_symmetrizing(self):
        self.plan(count=1)
        self.bundle(reverse_conflict=True)
        report = self.run_audit()
        q = report["products"]["epochs"][0]["solutions"]["PRODUCT_FIXED"]["full_covariance"]
        self.assertEqual(q["contradictory_duplicate_pairs"], 1)
        self.assertGreater(q["maximum_observed_directional_asymmetry_m2"], .9)
        self.assertFalse(q["internally_consistent"])

    def test_nonfinite_usable_phase_cannot_pass(self):
        self.plan(count=1)
        self.bundle(nonfinite=True)
        report = self.run_audit()
        self.assertEqual(report["evidence_issues"]["USABLE_AR_WITH_INVALID_NUMERIC_OR_AUTHORITY_FLAGS"], 1)
        self.assertEqual(report["products"]["totals_by_solution_system"]["PRODUCT_FIXED:G"]["nonfinite:phase_m"], 1)
        json.dumps(report, allow_nan=False)

    def test_old_gps_output_is_not_relabelled_new_three_system_completion(self):
        self.plan(count=1)
        self.bundle(systems=("G",))
        report = self.run_audit()
        self.assertEqual(report["actual_product_systems"], ["G"])
        self.assertEqual(report["systems_without_any_product_rows"], ["C", "E"])
        self.assertIn("ACTUAL_PRODUCT_SYSTEMS_DIFFER_FROM_STAGE_NOT_NEW_MULTI_RESULT", report["evidence_issues"])
        self.assertEqual(report["execution_classification"], "COMPLETED_WITH_EVIDENCE_FAILURES")

    def test_spp_position_and_no_failure_do_not_prove_single(self):
        self.plan(state="PREPARED_NOT_RUN", code=None)
        (self.root / "outputs/TEST.TRACE").write_text("sppos pos: 6378137 0 0\n"
            "sppos: time=2024-07-17 00:00:00, nObs=14\n")
        report = self.run_audit()
        self.assertEqual(report["logs"]["spp_single"]["status"], "UNAVAILABLE")
        self.assertEqual(report["logs"]["spp_input_observations"]["sum_nObs"], 14)

    def test_product_gap_and_declared_segmentation_are_separate(self):
        self.plan(count=3)
        self.bundle()
        self.bundle(1405209660, segment="S2", missing_endpoint=("E01", "L5Q"))
        report = self.run_audit()
        self.assertEqual(report["products"]["missing_grid_epochs_between_observed_products"], 1)
        events = report["continuity"]["event_counts"]
        self.assertEqual(events["PRODUCT_ENDPOINT_GAP"], 5)
        self.assertEqual(events["DECLARED_SEGMENT_DATUM_OR_COUNTER_CHANGE"], 5)
        self.assertEqual(events["PREVIOUSLY_PUBLISHED_ENDPOINT_ABSENT_THIS_EPOCH"], 1)
        second = report["products"]["epochs"][1]["solutions"]["PRODUCT_FIXED"]["systems"]["E"]
        self.assertEqual(second["satellites_missing_published_frequency_partner"], 1)

    def test_failed_run_keeps_real_partial_output_without_declaring_complete(self):
        self.plan(state="EXECUTION_FAILED", code=1)
        self.bundle()
        report = self.run_audit()
        self.assertEqual(report["execution_classification"], "EXECUTION_FAILED")
        self.assertEqual(report["products"]["epoch_count"], 1)

    def test_251_epochs_keep_only_current_covariance_in_memory(self):
        self.plan(systems=("gps",), count=251)
        for epoch in range(251):
            self.bundle(1405209600 + epoch * 30, systems=("G",))
        report = self.run_audit(max_events=5)
        self.assertEqual(report["products"]["epoch_count"], 251)
        self.assertEqual(report["streaming_memory"]["maximum_covariance_rows_in_memory"], 6)
        self.assertEqual(report["streaming_memory"]["maximum_covariance_parameters_in_memory"], 3)
        self.assertEqual(len(report["products"]["epochs"]), 251)

    def test_explicit_memory_limit_yields_unavailable_not_partial_psd_claim(self):
        self.plan(count=1)
        self.bundle()
        report = self.run_audit(max_parameters=4)
        q = report["products"]["epochs"][0]["solutions"]["PRODUCT_FIXED"]["full_covariance"]
        self.assertEqual(q["psd"]["status"], "UNAVAILABLE")
        self.assertFalse(q["internally_consistent"])

    def test_nonmonotonic_aggregate_refused_instead_of_loading_all_history(self):
        self.plan(systems=("gps",), count=2)
        with (self.root / "outputs/zhang_internal_products.csv").open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=["gpst_seconds", "solution"])
            writer.writeheader()
            writer.writerows([{"gpst_seconds": 30, "solution": "FLOAT"}, {"gpst_seconds": 0, "solution": "FLOAT"}])
        report = self.run_audit()
        self.assertIn("AGGREGATE_PARSE_FAILED", report["evidence_issues"])

    def test_pending_files_are_execution_evidence_not_prepared_or_published(self):
        self.plan(state="PREPARED_NOT_RUN", code=None)
        (self.root / "outputs/zhang_internal_products.csv.epochs/epoch.pending").mkdir(parents=True)
        report = self.run_audit()
        self.assertEqual(report["execution_classification"], "RUN_EVIDENCE_PRESENT_COMPLETION_UNCONFIRMED")
        self.assertEqual(report["atomic"]["pending_directories_not_published"], 1)
        self.assertEqual(report["atomic"]["hash_verified_bundles"], 0)

    def test_full_q_must_cover_clock_as_well_as_phase(self):
        self.plan(count=1)
        path = self.bundle()
        with (path / "covariance.csv").open(newline="") as stream:
            rows = [r for r in csv.DictReader(stream) if not
                    (r["row_satellite"] == "G01" and r["row_parameter"] == "CLOCK") and not
                    (r["column_satellite"] == "G01" and r["column_parameter"] == "CLOCK")]
        with (path / "covariance.csv").open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=audit.COV_COLUMNS)
            writer.writeheader()
            writer.writerows(rows)
        manifest = json.loads((path / "COMMIT.json").read_text())
        manifest["covariance_bytes"] = (path / "covariance.csv").stat().st_size
        (path / "COMMIT.json").write_text(json.dumps(manifest))
        self.rehash(path)
        report = self.run_audit()
        q = report["products"]["epochs"][0]["solutions"]["PRODUCT_FIXED"]["full_covariance"]
        self.assertEqual(q["missing_all_product_parameters"], [["G01", "CLOCK", "NONE"]])
        self.assertFalse(q["internally_consistent"])

    def test_manifest_per_system_claims_checked_against_actual_rows(self):
        self.plan(count=1)
        path = self.bundle()
        manifest = json.loads((path / "COMMIT.json").read_text())
        manifest["systems"]["GAL"]["pppar_usable_phase_rows"] = 100
        (path / "COMMIT.json").write_text(json.dumps(manifest))
        self.rehash(path)
        report = self.run_audit()
        self.assertFalse(report["atomic"]["bundles"][0]["integrity_verified"])
        self.assertIn("SYSTEM_PRODUCT_COUNT_OR_FLAG_MISMATCH:GAL", report["atomic"]["bundles"][0]["errors"])

    def test_production_hh_mm_ss_epoch_duration_is_seconds_not_leading_zero(self):
        self.plan(count=1)
        (self.root / "runner.log").write_text(
            "Processed epoch #1 - GPS time: 2323 259200 - 2024-07-17 00:00:00.00 "
            "(data handling took 00:00:03.07, processing took 00:00:01.46, epoch took 00:00:04.53)\n")
        timing = self.run_audit()["logs"]["epoch_timing"]
        self.assertEqual(timing["status"], "AVAILABLE")
        self.assertEqual(timing["sample_count"], 1)
        self.assertAlmostEqual(timing["p50"], 4.53)
        self.assertAlmostEqual(timing["maximum"], 4.53)
        parsed = audit.audit_logs(self.root)["processed_epochs"]
        self.assertEqual(parsed["status"], "AVAILABLE")
        self.assertEqual(parsed["count"], 1)
        self.assertEqual(parsed["first_gpst_seconds"], 1405209600)
        self.assertEqual(parsed["last_gpst_seconds"], 2323 * 604800 + 259200)

    def test_pea_two_digit_calendar_fraction_works_on_python_310(self):
        for value, expected in (("2024-07-17 00:00:00.00", 1405209600),
                                ("2024-07-17 00:00:00.50", 1405209600.5),
                                ("2024-07-17 00:00:00", 1405209600),
                                ("2024-07-17T00:00:00.050Z", 1405209600.05),
                                ("1405209600.5", 1405209600.5)):
            self.assertAlmostEqual(audit.gpst(value), expected)

    def test_invalid_processed_calendar_cannot_be_available_with_null_time(self):
        (self.root / "runner.log").write_text(
            "Processed epoch #1 - GPS time: 2323 259200 - 2024-02-30 00:00:00.00 (epoch took 00:00:04.53)\n")
        report = audit.audit_logs(self.root)
        self.assertEqual(report["processed_epochs"]["status"], "UNAVAILABLE")
        self.assertEqual(report["invalid_processed_epoch_calendar_records"], 1)
        self.assertEqual(report["epoch_timing"]["status"], "UNAVAILABLE")
        self.assertIsNone(audit.gpst("2024-02-30 00:00:00.00"))

    def test_duration_parser_retains_numeric_seconds_and_arbitrary_hours(self):
        for text, seconds in (("4.53", 4.53), (".5", .5), ("1.2e2", 120),
                              ("01:02:03.4", 3723.4), ("123:59:59.99", 446399.99), ("0", 0)):
            self.assertAlmostEqual(audit.duration_seconds(text), seconds)

    def test_malformed_or_missing_duration_is_unavailable_not_zero(self):
        for text in (None, "nan", "inf", "garbage", "00:99:03", "00:00:99", "-1"):
            self.assertIsNone(audit.duration_seconds(text))
        (self.root / "runner.log").write_text(
            "Processed epoch #1 - GPS time: 2323 259200 - 2024-07-17 00:00:00.00 (epoch took nan)\n")
        self.assertEqual(audit.audit_logs(self.root)["epoch_timing"]["status"], "UNAVAILABLE")


if __name__ == "__main__":
    unittest.main()
