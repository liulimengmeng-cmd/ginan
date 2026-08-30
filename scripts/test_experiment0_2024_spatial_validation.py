from __future__ import annotations

import math
import io
import json
import sys
import tempfile
import unittest
from decimal import Decimal
from pathlib import Path
from unittest.mock import patch

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

import experiment0_2024_spatial_validation as spatial_validation
from experiment0_2024_spatial_validation import (
    EXPECTED_FORWARD_STAGE,
    CalibrationRecord,
    METHODS,
    PredictionRow,
    SpatialEpoch,
    SpatialInputError,
    WGS84_A,
    _bootstrap_comparisons,
    _valid_spatial_epoch,
    build_freeze,
    calibration_nll,
    calibration_variances,
    collect_loso_calibration,
    comparison_metrics,
    fit_spatial_plane,
    optimise_kappa,
    prediction_target_covariance,
    predict_spatial_plane,
    ray_geometry,
    iter_geometry_epochs,
)


class Experiment0SpatialValidationTests(unittest.TestCase):
    @staticmethod
    def _joint_spatial_epoch() -> SpatialEpoch:
        design = np.array(
            [
                [1.0, -2.0, -1.0],
                [1.0, -1.0, 1.0],
                [1.0, 0.0, -2.0],
                [1.0, 0.0, 2.0],
                [1.0, 1.0, -1.0],
                [1.0, 1.0, 1.0],
                [1.0, 2.0, -2.0],
                [1.0, 2.0, 2.0],
            ]
        )
        sites = np.array(["A", "A", "B", "B", "C", "C", "D", "D"], dtype=object)
        covariance = 0.8 * np.eye(8) + 0.2 * np.ones((8, 8))
        return SpatialEpoch(
            gps_week=2313,
            gps_tow=Decimal(259200),
            status="OK",
            posterior_stage="FILTER_POSTERIOR_NO_EPOCH_AR",
            sites=sites,
            target_satellites=[f"G{index + 2:02d}" for index in range(8)],
            reference_satellites=["G01"] * 8,
            design=design,
            values=design @ np.array([5.0, 2.0, -1.0]),
            covariance=covariance,
            declared_difference_count=8,
            dropped_geometry_count=0,
        )

    def test_zenith_ray_has_unit_mapping_and_receiver_ipp(self) -> None:
        receiver = np.array([WGS84_A, 0.0, 0.0])
        satellite = np.array([26_000_000.0, 0.0, 0.0])
        geometry = ray_geometry("TEST", "G02", 0, receiver, satellite)
        self.assertAlmostEqual(geometry.elevation_rad, math.pi / 2, places=12)
        self.assertAlmostEqual(geometry.ipp_lat_rad, 0.0, places=12)
        self.assertAlmostEqual(geometry.ipp_lon_rad, 0.0, places=12)
        self.assertAlmostEqual(geometry.mapping_factor, 1.0, places=12)

    def test_all_covariance_arms_recover_exact_plane(self) -> None:
        design = np.array(
            [
                [1.0, -1.0, -0.5],
                [0.5, 0.2, 1.0],
                [-0.3, 1.2, 0.1],
                [1.3, -0.4, 0.8],
                [-0.7, -0.5, 1.5],
                [0.2, 0.8, -1.1],
            ]
        )
        truth = np.array([8.0, -2.0, 3.5])
        values = design @ truth
        base = np.array(
            [
                [2.0, 0.2, 0.0, 0.0, 0.0, 0.0],
                [0.2, 1.5, 0.1, 0.0, 0.0, 0.0],
                [0.0, 0.1, 1.8, 0.2, 0.0, 0.0],
                [0.0, 0.0, 0.2, 2.2, 0.1, 0.0],
                [0.0, 0.0, 0.0, 0.1, 1.7, 0.2],
                [0.0, 0.0, 0.0, 0.0, 0.2, 1.9],
            ]
        )
        for method in METHODS:
            fit = fit_spatial_plane(design, values, base, method)
            prediction, variance = predict_spatial_plane(fit, design[:2])
            np.testing.assert_allclose(fit.beta, truth, atol=1e-10)
            np.testing.assert_allclose(prediction, values[:2], atol=1e-10)
            self.assertTrue(np.all(variance >= 0))

    def test_all_covariance_arms_reject_rank_deficient_plane(self) -> None:
        design = np.array(
            [
                [1.0, 0.0, 0.0],
                [1.0, 1.0, 2.0],
                [1.0, 2.0, 4.0],
                [1.0, 3.0, 6.0],
            ]
        )
        values = np.arange(4, dtype=float)
        covariance = np.eye(4)
        for method in METHODS:
            with self.subTest(method=method):
                with self.assertRaises(SpatialInputError):
                    fit_spatial_plane(design, values, covariance, method)

    def test_wrong_forward_posterior_stage_is_rejected(self) -> None:
        epoch = self._joint_spatial_epoch()
        self.assertEqual(epoch.posterior_stage, EXPECTED_FORWARD_STAGE)
        self.assertTrue(_valid_spatial_epoch(epoch, set()))
        epoch.posterior_stage = "SMOOTHED_POSTERIOR"
        self.assertFalse(_valid_spatial_epoch(epoch, set()))

    def test_nll_calibration_finds_expected_scale(self) -> None:
        records = [
            CalibrationRecord(
                gps_week=2313,
                gps_tow=Decimal(index * 30),
                site="TEST",
                residual=2.0,
                prediction_variance=1.0,
                target_variance=0.1,
            )
            for index in range(120)
        ]
        kappa, objective = optimise_kappa(records)
        self.assertTrue(math.isfinite(objective))
        self.assertAlmostEqual(kappa, math.sqrt(3.9), places=3)
        self.assertLess(objective, calibration_nll(records, 1.0))

    def test_kappa_search_refines_below_smallest_positive_grid_point(self) -> None:
        expected_kappa = 8e-5
        target_variance = 1e-8
        records = [
            CalibrationRecord(
                gps_week=2313,
                gps_tow=Decimal(0),
                site="TEST",
                residual=math.sqrt(target_variance + expected_kappa**2),
                prediction_variance=1.0,
                target_variance=target_variance,
            )
        ]
        with patch.object(
            spatial_validation,
            "_calibration_arrays",
            wraps=spatial_validation._calibration_arrays,
        ) as extractor:
            kappa, objective = optimise_kappa(records)

        self.assertEqual(extractor.call_count, 1)
        self.assertAlmostEqual(kappa, expected_kappa, places=10)
        self.assertLess(objective, calibration_nll(records, 0.0))
        self.assertLess(objective, calibration_nll(records, 1e-4))

    def test_kappa_search_accepts_physical_zero_boundary(self) -> None:
        records = [
            CalibrationRecord(
                gps_week=2313,
                gps_tow=Decimal(0),
                site="TEST",
                residual=1.0,
                prediction_variance=1.0,
                target_variance=1.0,
            )
        ]
        kappa, objective = optimise_kappa(records)
        self.assertEqual(kappa, 0.0)
        self.assertTrue(math.isfinite(objective))

    def test_kappa_search_rejects_upper_grid_boundary(self) -> None:
        records = [
            CalibrationRecord(
                gps_week=2313,
                gps_tow=Decimal(0),
                site="TEST",
                residual=20_000.0,
                prediction_variance=1.0,
                target_variance=1.0,
            )
        ]
        with self.assertRaisesRegex(SpatialInputError, "upper search boundary 10000"):
            optimise_kappa(records)

    def test_freeze_records_no_kappa_boundary_hit(self) -> None:
        record = CalibrationRecord(
            gps_week=2313,
            gps_tow=Decimal(0),
            site="TEST",
            residual=1.0,
            prediction_variance=0.5,
            target_variance=0.5,
        )
        records = {method: [record] for method in METHODS}
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            quiet_stec = root / "quiet.STEC"
            quiet_sd = root / "quiet.STEC.SD.COV"
            product_audit = root / "products.json"
            split_audit = root / "split.json"
            quiet_stec.write_text("fixture\n", encoding="utf-8")
            quiet_sd.write_text("fixture\n", encoding="utf-8")
            product_audit.write_text("{}\n", encoding="utf-8")
            split_audit.write_text(
                json.dumps(
                    {
                        "valid": True,
                        "model_split_valid": True,
                        "model_stations": ["TEST"],
                    }
                ),
                encoding="utf-8",
            )
            with patch(
                "experiment0_2024_spatial_validation.collect_loso_calibration",
                return_value=(records, {"fixture": True}),
            ):
                freeze = build_freeze(
                    quiet_stec, quiet_sd, product_audit, split_audit
                )

        for method in METHODS:
            self.assertIs(freeze["calibration"][method]["boundary_hit"], False)

    def test_full_cross_covariance_matches_independent_linear_transform(
        self,
    ) -> None:
        generator = np.random.default_rng(20240830)
        training_count = 7
        target_count = 3
        training_design = np.column_stack(
            [np.ones(training_count), generator.normal(size=(training_count, 2))]
        )
        target_design = np.column_stack(
            [np.ones(target_count), generator.normal(size=(target_count, 2))]
        )
        raw = generator.normal(
            size=(training_count + target_count, training_count + target_count)
        )
        factor = np.geomspace(0.4, 2.5, training_count + target_count)[:, None] * raw
        joint_covariance = factor @ factor.T + np.diag(
            np.linspace(0.3, 1.3, training_count + target_count)
        )
        training_covariance = joint_covariance[:training_count, :training_count]
        training_target_covariance = joint_covariance[
            :training_count, training_count:
        ]

        fit = fit_spatial_plane(
            training_design,
            generator.normal(size=training_count),
            training_covariance,
            "full",
        )
        _, actual_prediction_variance = predict_spatial_plane(fit, target_design)
        actual_cross = prediction_target_covariance(
            fit, target_design, training_target_covariance
        )

        covariance_inverse_design = np.linalg.solve(
            training_covariance, training_design
        )
        direct_beta_covariance = np.linalg.solve(
            training_design.T @ covariance_inverse_design, np.eye(3)
        )
        linear_map = (
            target_design
            @ direct_beta_covariance
            @ covariance_inverse_design.T
        )
        expected_prediction_covariance = (
            linear_map @ training_covariance @ linear_map.T
        )
        expected_cross_matrix = linear_map @ training_target_covariance
        residual_transform = np.hstack((-linear_map, np.eye(target_count)))
        direct_residual_covariance = (
            residual_transform @ joint_covariance @ residual_transform.T
        )

        np.testing.assert_allclose(
            actual_prediction_variance,
            np.diag(expected_prediction_covariance),
            rtol=1e-10,
            atol=1e-11,
        )
        np.testing.assert_allclose(
            actual_cross,
            np.diag(expected_cross_matrix),
            rtol=1e-10,
            atol=1e-11,
        )
        np.testing.assert_allclose(
            np.diag(direct_residual_covariance),
            np.diag(joint_covariance[training_count:, training_count:])
            + actual_prediction_variance
            - 2 * actual_cross,
            rtol=1e-10,
            atol=1e-11,
        )

    def test_full_loso_uses_known_training_target_cross_covariance(self) -> None:
        epoch = self._joint_spatial_epoch()
        with (
            patch(
                "experiment0_2024_spatial_validation.product_gap_keys",
                return_value=set(),
            ),
            patch(
                "experiment0_2024_spatial_validation.iter_spatial_epochs",
                return_value=iter([epoch]),
            ),
        ):
            records, summary = collect_loso_calibration(
                Path("model.STEC"),
                Path("model.STEC.SD.COV"),
                Path("product-audit.json"),
                {"A", "B", "C", "D"},
            )

        self.assertEqual(
            {method: len(records[method]) for method in METHODS},
            {method: 8 for method in METHODS},
        )
        self.assertTrue(summary["strict_common_method_rows"])
        self.assertTrue(summary["identical_method_record_key_sequence"])
        key_sequences = {
            method: [record.key for record in records[method]] for method in METHODS
        }
        for method in METHODS[1:]:
            self.assertListEqual(key_sequences[method], key_sequences[METHODS[0]])
        self.assertTrue(all(record.target_satellite for record in records["full"]))
        self.assertTrue(all(record.reference_satellite for record in records["full"]))
        self.assertTrue(
            all(record.prediction_target_covariance == 0 for record in records["none"])
        )
        self.assertTrue(
            all(
                record.prediction_target_covariance == 0
                for record in records["diagonal"]
            )
        )
        self.assertTrue(
            any(
                abs(record.prediction_target_covariance) > 1e-12
                for record in records["full"]
            )
        )
        record = records["full"][0]
        expected_variance = (
            record.prediction_variance
            + record.target_variance
            - 2 * record.prediction_target_covariance
        )
        self.assertAlmostEqual(
            calibration_variances([record], 1.0)[0], expected_variance
        )
        expected_nll = 0.5 * (
            math.log(2 * math.pi * expected_variance)
            + record.residual * record.residual / expected_variance
        )
        self.assertAlmostEqual(calibration_nll([record], 1.0), expected_nll)

    def test_nonpositive_cross_covariance_variance_has_infinite_nll(self) -> None:
        record = CalibrationRecord(
            gps_week=2313,
            gps_tow=Decimal(0),
            site="TEST",
            residual=1.0,
            prediction_variance=0.1,
            target_variance=0.1,
            prediction_target_covariance=0.2,
        )
        self.assertLess(calibration_variances([record], 1.0)[0], 0)
        self.assertTrue(math.isinf(calibration_nll([record], 1.0)))

    def test_loso_aborts_if_any_covariance_arm_fit_fails(self) -> None:
        epoch = self._joint_spatial_epoch()
        real_fit = fit_spatial_plane

        def fail_full(design, values, covariance, method):
            if method == "full":
                raise SpatialInputError("injected full-arm failure")
            return real_fit(design, values, covariance, method)

        with (
            patch(
                "experiment0_2024_spatial_validation.product_gap_keys",
                return_value=set(),
            ),
            patch(
                "experiment0_2024_spatial_validation.iter_spatial_epochs",
                return_value=iter([epoch]),
            ),
            patch(
                "experiment0_2024_spatial_validation.fit_spatial_plane",
                side_effect=fail_full,
            ),
        ):
            with self.assertRaisesRegex(
                SpatialInputError, "refusing method-specific rows"
            ):
                collect_loso_calibration(
                    Path("model.STEC"),
                    Path("model.STEC.SD.COV"),
                    Path("product-audit.json"),
                    {"A", "B", "C", "D"},
                )

    def test_full_vs_diagonal_metrics_are_directional_log_score_differences(
        self,
    ) -> None:
        scores = {
            "none": {
                "calibrated": {
                    "coverage_percent": 80.0,
                    "mean_interval_score_tecu": 10.0,
                    "mean_width_tecu": 2.0,
                    "nll": 3.0,
                }
            },
            "diagonal": {
                "calibrated": {
                    "coverage_percent": 94.0,
                    "mean_interval_score_tecu": 8.0,
                    "mean_width_tecu": 2.0,
                    "nll": 2.0,
                }
            },
            "full": {
                "calibrated": {
                    "coverage_percent": 91.0,
                    "mean_interval_score_tecu": 7.0,
                    "mean_width_tecu": 2.1,
                    "nll": 1.7,
                }
            },
        }
        metrics = comparison_metrics(scores)
        self.assertEqual(
            metrics["full_vs_diagonal_coverage_error_reduction_points"], -3.0
        )
        self.assertAlmostEqual(
            metrics["full_vs_diagonal_mean_log_score_improvement_nats_per_prediction"],
            0.3,
        )
        self.assertFalse(metrics["full_vs_diagonal_coverage_material_benefit"])

    def test_bootstrap_log_score_interval_preserves_improvement_direction(
        self,
    ) -> None:
        rows: list[PredictionRow] = []
        for hour, diagonal_prediction in enumerate((2.0, 4.0)):
            for index in range(2):
                methods = {
                    method: {
                        "predicted_tecu": (
                            0.0 if method == "full" else diagonal_prediction
                        ),
                        "prediction_variance_tecu2": 0.5,
                        "target_variance_tecu2": 0.5,
                        "uncalibrated_total_variance_tecu2": 1.0,
                        "calibrated_total_variance_tecu2": 1.0,
                    }
                    for method in METHODS
                }
                rows.append(
                    PredictionRow(
                        gps_week=2313,
                        gps_tow=Decimal(hour * 3600 + index * 30),
                        site="TEST",
                        target_satellite=f"G{index + 2:02d}",
                        reference_satellite="G01",
                        observed_tecu=0.0,
                        methods=methods,
                    )
                )
        with patch.object(spatial_validation, "BOOTSTRAP_REPLICATES", 64):
            bootstrap = _bootstrap_comparisons(rows)
        interval = bootstrap["percentile_95_intervals"][
            "full_vs_diagonal_mean_log_score_improvement_nats_per_prediction"
        ]
        self.assertGreater(interval["lower"], 0.0)

    def test_blank_observation_geometry_is_dropped_without_parse_failure(self) -> None:
        payload = (
            "IONO_MEA,2313,259200.000,TEST,G02,1.0,2.0e+00, , ,             ,"
            "             ,             ,    ,  6378137.000,        0.000,        0.000\n"
        )
        epochs = list(iter_geometry_epochs(io.StringIO(payload)))
        self.assertEqual(len(epochs), 1)
        self.assertEqual(epochs[0][0], (2313, Decimal("259200.000")))
        self.assertEqual(epochs[0][1], {})


if __name__ == "__main__":
    unittest.main()
