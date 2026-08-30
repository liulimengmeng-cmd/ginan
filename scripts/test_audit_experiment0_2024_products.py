from __future__ import annotations

import gzip
import json
import sys
import tempfile
import unittest
from collections import Counter
from datetime import date, timedelta
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from audit_experiment0_2024_products import (
    CLK_INTERVAL_SECONDS,
    EXPERIMENTS,
    REQUIRED_OBSERVABLES,
    SP3_INTERVAL_SECONDS,
    EpochRecords,
    audit_data_root,
    expected_epochs,
    summarize_epoch_records,
    write_json_atomic,
)


def write_sp3(
    path: Path,
    day: date,
    satellites: tuple[str, ...],
    missing: set[tuple[str, int]] | None = None,
) -> None:
    missing = missing or set()
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wt", encoding="ascii", newline="\n") as stream:
        stream.write("#cP synthetic Experiment 0 SP3\n")
        for index, epoch in enumerate(expected_epochs(day, SP3_INTERVAL_SECONDS)):
            stream.write(
                f"*  {epoch.year:4d} {epoch.month:2d} {epoch.day:2d} "
                f"{epoch.hour:2d} {epoch.minute:2d} {epoch.second:2d}.00000000\n"
            )
            for satellite in satellites:
                if (satellite, index) in missing:
                    continue
                stream.write(
                    f"P{satellite}  15600.000000  20100.000000  21700.000000"
                    "      0.000000\n"
                )


def write_clk(
    path: Path,
    day: date,
    satellites: tuple[str, ...],
    missing: set[tuple[str, int]] | None = None,
) -> None:
    missing = missing or set()
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wt", encoding="ascii", newline="\n") as stream:
        stream.write("     3.04           C                                       RINEX VERSION / TYPE\n")
        stream.write("                                                            END OF HEADER\n")
        for index, epoch in enumerate(expected_epochs(day, CLK_INTERVAL_SECONDS)):
            for satellite in satellites:
                if (satellite, index) in missing:
                    continue
                stream.write(
                    f"AS {satellite}  {epoch.year:4d} {epoch.month:02d} {epoch.day:02d} "
                    f"{epoch.hour:02d} {epoch.minute:02d} {epoch.second:2d}.000000  "
                    "1  1.234567890123e-04\n"
                )


def write_bia(
    path: Path,
    day: date,
    complete_satellites: tuple[str, ...],
    partial_satellites: dict[str, tuple[str, ...]] | None = None,
) -> None:
    partial_satellites = partial_satellites or {}
    path.parent.mkdir(parents=True, exist_ok=True)
    next_day = day + timedelta(days=1)
    start = f"{day.year}:{day.timetuple().tm_yday:03d}:00000"
    end = f"{next_day.year}:{next_day.timetuple().tm_yday:03d}:00000"
    lines = [
        "%=BIA 1.00 TST synthetic\n",
        "+BIAS/DESCRIPTION\n",
        " OBSERVATION_SAMPLING                             30\n",
        " PARAMETER_SPACING                             86400\n",
        " DETERMINATION_METHOD                    PSEUDO-ABSOLUTE_BIAS_ESTIMATION\n",
        " BIAS_MODE                               ABSOLUTE\n",
        " TIME_SYSTEM                             G\n",
        " SATELLITE_CLOCK_REFERENCE_OBSERVABLES   G  C1W  C2W\n",
        " APC_MODEL                               IGS20_2303.ATX\n",
        "-BIAS/DESCRIPTION\n",
        "+BIAS/SOLUTION\n",
    ]
    for satellite in complete_satellites:
        for observable in REQUIRED_OBSERVABLES:
            lines.append(
                f" OSB  G000 {satellite}           {observable}       {start} {end} "
                "ns       0.125000000000000    0.010000\n"
            )
    for satellite, observables in partial_satellites.items():
        for observable in observables:
            lines.append(
                f" OSB  G000 {satellite}           {observable}       {start} {end} "
                "ns       0.125000000000000    0.010000\n"
            )
    lines.extend(["-BIAS/SOLUTION\n", "%=ENDBIA\n"])

    if path.suffix.lower() == ".gz":
        with gzip.open(path, "wt", encoding="ascii", newline="\n") as stream:
            stream.writelines(lines)
    else:
        path.write_text("".join(lines), encoding="ascii", newline="\n")


def build_complete_tree(
    root: Path,
    clock_missing_by_label: dict[str, set[tuple[str, int]]] | None = None,
) -> None:
    clock_missing_by_label = clock_missing_by_label or {}
    for label, specification in EXPERIMENTS.items():
        day = specification["day"]
        doy = specification["doy"]
        eligible = tuple(specification["eligible"])
        assert isinstance(day, date)
        assert isinstance(doy, int)
        directory = root / label / "products"
        prefix = f"WUM0MGXFIN_2024{doy:03d}0000_"
        write_sp3(directory / (prefix + "01D_05M_ORB.SP3"), day, eligible)
        write_clk(
            directory / (prefix + "01D_30S_CLK.CLK"),
            day,
            eligible,
            clock_missing_by_label.get(label),
        )
        bia_path = directory / (prefix + "01D_01D_OSB.BIA")
        partial: dict[str, tuple[str, ...]] = {}
        if label == "quiet":
            partial = {"G10": ("C1C", "C2W")}
        else:
            bia_path = bia_path.with_name(bia_path.name + ".gz")
        write_bia(bia_path, day, eligible, partial)


class Experiment0ProductAuditTests(unittest.TestCase):
    def test_expected_intersections_metadata_and_atomic_json(self) -> None:
        with tempfile.TemporaryDirectory() as directory_name:
            root = Path(directory_name)
            build_complete_tree(root)

            payload = audit_data_root(root)

            self.assertTrue(payload["validation_passed"])
            quiet = payload["experiments"]["quiet"]
            storm = payload["experiments"]["storm"]
            self.assertEqual(
                quiet["eligible_gps_satellites"],
                list(EXPERIMENTS["quiet"]["eligible"]),
            )
            self.assertEqual(
                storm["eligible_gps_satellites"],
                list(EXPERIMENTS["storm"]["eligible"]),
            )
            self.assertEqual(quiet["sp3"]["expected_epoch_count"], 288)
            self.assertEqual(quiet["clk"]["expected_epoch_count"], 2880)
            self.assertEqual(
                quiet["bia"]["metadata"]["apc_model"], "IGS20_2303.ATX"
            )
            self.assertEqual(
                quiet["bia"]["metadata"][
                    "satellite_clock_reference_observables"
                ]["G"],
                ["C1W", "C2W"],
            )
            self.assertTrue(
                quiet["bia"]["satellites"]["G10"]["observables"]["C1C"][
                    "covers_full_day"
                ]
            )
            self.assertFalse(
                quiet["bia"]["satellites"]["G10"]["observables"]["L1C"][
                    "covers_full_day"
                ]
            )
            self.assertEqual(
                quiet["clk"]["satellites"]["G01"]["missing_epoch_count"],
                2880,
            )

            output = root / "audit.json"
            output.write_text("stale\n", encoding="utf-8")
            write_json_atomic(payload, output)
            written = json.loads(output.read_text(encoding="utf-8"))
            self.assertTrue(written["validation_passed"])
            self.assertFalse(list(root.glob(".audit.json.*.tmp")))

    def test_eligible_roster_is_distinct_from_strict_clock_completeness(self) -> None:
        with tempfile.TemporaryDirectory() as directory_name:
            root = Path(directory_name)
            quiet_last_eight = {("G32", index) for index in range(2872, 2880)}
            storm_shared_ten = {
                (satellite, index)
                for satellite in EXPERIMENTS["storm"]["eligible"]
                for index in range(430, 440)
            }
            build_complete_tree(
                root,
                {
                    "quiet": quiet_last_eight,
                    "storm": storm_shared_ten,
                },
            )

            payload = audit_data_root(root)

            self.assertFalse(payload["validation_passed"])
            quiet = payload["experiments"]["quiet"]
            storm = payload["experiments"]["storm"]
            self.assertEqual(
                quiet["eligible_gps_satellites"],
                list(EXPERIMENTS["quiet"]["eligible"]),
            )
            self.assertEqual(
                storm["eligible_gps_satellites"],
                list(EXPERIMENTS["storm"]["eligible"]),
            )
            self.assertTrue(quiet["validation"]["eligible_intersection_passed"])
            self.assertTrue(storm["validation"]["eligible_intersection_passed"])
            self.assertFalse(quiet["validation"]["full_day_completeness_passed"])
            self.assertFalse(storm["validation"]["full_day_completeness_passed"])
            self.assertEqual(
                quiet["validation"]["eligible_satellites_with_coverage_defects"],
                ["G32"],
            )
            self.assertEqual(
                storm["validation"]["eligible_satellites_with_coverage_defects"],
                list(EXPERIMENTS["storm"]["eligible"]),
            )
            self.assertEqual(
                quiet["clk"]["satellites"]["G32"]["gaps"],
                [
                    {
                        "first_epoch": "2024-05-08T23:56:00+00:00",
                        "last_epoch": "2024-05-08T23:59:30+00:00",
                        "count": 8,
                    }
                ],
            )
            self.assertEqual(
                storm["clk"]["satellites"]["G02"]["gaps"],
                [
                    {
                        "first_epoch": "2024-05-11T03:35:00+00:00",
                        "last_epoch": "2024-05-11T03:39:30+00:00",
                        "count": 10,
                    }
                ],
            )

    def test_one_clock_gap_is_reported_and_removes_completeness(self) -> None:
        day = date(2024, 5, 8)
        epochs = expected_epochs(day, CLK_INTERVAL_SECONDS)
        missing_index = 137
        records = EpochRecords()
        records.valid["G02"] = Counter(
            epoch for index, epoch in enumerate(epochs) if index != missing_index
        )
        records.raw_count["G02"] = len(epochs) - 1

        summaries, complete = summarize_epoch_records(
            records, day, CLK_INTERVAL_SECONDS
        )

        self.assertNotIn("G02", complete)
        self.assertEqual(summaries["G02"]["missing_epoch_count"], 1)
        self.assertEqual(summaries["G02"]["gaps"][0]["count"], 1)
        self.assertEqual(
            summaries["G02"]["gaps"][0]["first_epoch"],
            epochs[missing_index].isoformat(),
        )


if __name__ == "__main__":
    unittest.main()
