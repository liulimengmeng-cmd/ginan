from __future__ import annotations

import sys
import unittest
from datetime import date
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from fetch_experiment0_2024_data import (
    HELDOUT_STATIONS,
    ELIGIBLE_GPS_SATELLITES,
    MODEL_STATIONS,
    day_of_year,
    ga_output_filename,
    gps_day_of_week,
    gps_week,
)


class Experiment0DownloadTests(unittest.TestCase):
    def test_gps_week_and_doy(self) -> None:
        self.assertEqual(gps_week(date(2024, 5, 8)), 2313)
        self.assertEqual(gps_week(date(2024, 5, 11)), 2313)
        self.assertEqual(gps_week(date(2024, 5, 12)), 2314)
        self.assertEqual(gps_day_of_week(date(2024, 5, 8)), 3)
        self.assertEqual(gps_day_of_week(date(2024, 5, 11)), 6)
        self.assertEqual(day_of_year(date(2024, 5, 8)), 129)
        self.assertEqual(day_of_year(date(2024, 5, 11)), 132)

    def test_model_and_heldout_stations_are_frozen_and_disjoint(self) -> None:
        self.assertEqual(len(MODEL_STATIONS), 10)
        self.assertEqual(len(HELDOUT_STATIONS), 3)
        self.assertFalse(set(MODEL_STATIONS) & set(HELDOUT_STATIONS))
        self.assertEqual(HELDOUT_STATIONS, ("STR2", "BALL", "PARK"))

    def test_wum_eligibility_is_frozen_by_day(self) -> None:
        self.assertNotIn("G01", ELIGIBLE_GPS_SATELLITES["quiet"])
        self.assertNotIn("G10", ELIGIBLE_GPS_SATELLITES["quiet"])
        self.assertNotIn("G01", ELIGIBLE_GPS_SATELLITES["storm"])
        self.assertIn("G10", ELIGIBLE_GPS_SATELLITES["storm"])

    def test_ga_observation_filename_is_plain_rinex(self) -> None:
        entry = {
            "fileType": "obs",
            "fileLocation": (
                "https://example.invalid/default/"
                "MOBS00AUS_R_20241290000_01D_30S_MO.crx.gz?signature=discard"
            ),
        }
        self.assertEqual(
            ga_output_filename(entry),
            "MOBS00AUS_R_20241290000_01D_30S_MO.rnx",
        )

    def test_ga_navigation_filename_keeps_rinex_extension(self) -> None:
        entry = {
            "fileType": "nav",
            "fileLocation": (
                "https://example.invalid/default/"
                "BRDC00IGS_R_20241290000_01D_MN.rnx.gz?signature=discard"
            ),
        }
        self.assertEqual(
            ga_output_filename(entry),
            "BRDC00IGS_R_20241290000_01D_MN.rnx",
        )


if __name__ == "__main__":
    unittest.main()
