from __future__ import annotations

import math
import io
import sys
import unittest
from decimal import Decimal
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from experiment0_2024_spatial_validation import (
    CalibrationRecord,
    METHODS,
    WGS84_A,
    calibration_nll,
    fit_spatial_plane,
    optimise_kappa,
    predict_spatial_plane,
    ray_geometry,
    iter_geometry_epochs,
)


class Experiment0SpatialValidationTests(unittest.TestCase):
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
