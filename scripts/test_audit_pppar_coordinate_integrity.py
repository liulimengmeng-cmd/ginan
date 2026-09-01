from __future__ import annotations

import json
import math
import sys
import tempfile
import unittest
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from audit_pppar_coordinate_integrity import (
    Thresholds,
    audit_coordinate_integrity,
    main,
    parse_sinex,
    parse_sinex_epoch,
    parse_trace,
    parse_trace_epoch,
    select_sinex_reference,
)


def rec_pos(epoch: str, receiver: str, x: float, y: float, z: float) -> str:
    rows = []
    for axis, value in zip("XYZ", (x, y, z)):
        rows.append(
            f"*\t-1\t{epoch}\t     REC_POS\t    \t{receiver}\t"
            f"      {axis}\t {value:.7f}\t 1.0\t0.0\n"
        )
    return "".join(rows)


def state_block(name: str, rows: str) -> str:
    return f"+STATES/{name}\n{rows}-STATES/{name}\n"


def sinex(stations: dict[str, tuple[float, float, float]]) -> str:
    rows = [
        "% synthetic coordinate reference\n",
        "+SOLUTION/EPOCHS\n",
        "*CODE PT SOLN T _DATA_START_ __DATA_END__ _MEAN_EPOCH_\n",
    ]
    for receiver in stations:
        rows.append(
            f" {receiver:<4s} A    1 C 24:196:00000 24:202:86399 24:199:00000\n"
        )
    rows.extend(
        [
            "-SOLUTION/EPOCHS\n",
        "+SOLUTION/ESTIMATE\n",
        "*INDEX TYPE__ CODE PT SOLN _REF_EPOCH__ UNIT S __ESTIMATED_VALUE____ _STD_DEV___\n",
        ]
    )
    index = 1
    for receiver, values in stations.items():
        for estimate_type, value in zip(("STAX", "STAY", "STAZ"), values):
            rows.append(
                f" {index:5d} {estimate_type:<6s} {receiver:<4s} A    1 "
                f"24:199:00000 m    2 {value:21.14E} 1.00000E-03\n"
            )
            index += 1
    rows.append("-SOLUTION/ESTIMATE\n")
    return "".join(rows)


class CoordinateIntegrityAuditTests(unittest.TestCase):
    def test_trace_epoch_accepts_ginan_two_digit_fraction_on_all_supported_pythons(self) -> None:
        parsed = parse_trace_epoch("2024-07-17 00:00:00.00")
        self.assertEqual(parsed, datetime(2024, 7, 17, tzinfo=timezone.utc))

    def test_trace_keeps_last_complete_duplicate_without_mixing_partial_group(self) -> None:
        epoch = "2024-07-17 00:00:00.00"
        rows = rec_pos(epoch, "DYNG", 1.0, 2.0, 3.0)
        rows += rec_pos(epoch, "DYNG", 4.0, 5.0, 6.0)
        # An incomplete later fragment must not replace the last complete group.
        rows += (
            f"*\t-1\t{epoch}\t     REC_POS\t    \tDYNG\t      X\t 9.0\t1\t0\n"
            f"*\t-1\t{epoch}\t     REC_POS\t    \tDYNG\t      Y\t 9.0\t1\t0\n"
        )
        text = state_block("AR", rows)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "primary.trace"
            path.write_text(text, encoding="utf-8")
            parsed = parse_trace(path, "primary")

        sample = next(iter(parsed.samples.values()))
        self.assertEqual(sample.ecef_m, (4.0, 5.0, 6.0))
        self.assertEqual(parsed.metadata["complete_vector_occurrence_count"], 2)
        self.assertEqual(parsed.metadata["duplicate_complete_vector_count"], 1)
        self.assertEqual(parsed.metadata["discarded_incomplete_fragment_count"], 1)

    def test_constraint_attempt_marker_without_applying_row_is_not_an_event(self) -> None:
        epoch = "2024-07-17 00:00:00.00"
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "primary.trace"
            path.write_text(
                f"fixAndHoldAmbiguities: {epoch}\n"
                + state_block("AR", rec_pos(epoch, "DYNG", 1.0, 2.0, 3.0)),
                encoding="utf-8",
            )
            parsed = parse_trace(path, "primary")

        self.assertEqual(parsed.pseudoobs_attempt_rows_by_epoch, {})
        self.assertEqual(parsed.metadata["constraint_attempt_marker_count"], 1)
        self.assertEqual(
            parsed.metadata["visible_integer_pseudoobservation_attempt_row_count"], 0
        )

    def test_primary_selects_ar_block_and_reports_ambiguous_ar_blocks(self) -> None:
        epoch = "2024-07-17 00:00:00.00"
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "primary.trace"
            path.write_text(
                state_block("AR", rec_pos(epoch, "DYNG", 1.0, 2.0, 3.0))
                + state_block("PPP", rec_pos(epoch, "DYNG", 9.0, 9.0, 9.0)),
                encoding="utf-8",
            )
            parsed = parse_trace(path, "primary")
            self.assertEqual(next(iter(parsed.samples.values())).ecef_m, (1.0, 2.0, 3.0))
            self.assertEqual(parsed.metadata["selected_state_block"], "STATES/AR")

            path.write_text(
                state_block("AR", rec_pos(epoch, "DYNG", 1.0, 2.0, 3.0))
                + state_block("AR_RTS", rec_pos(epoch, "DYNG", 4.0, 5.0, 6.0)),
                encoding="utf-8",
            )
            ambiguous = parse_trace(path, "primary")

        self.assertEqual(ambiguous.samples, {})
        self.assertEqual(
            ambiguous.metadata["state_block_selection_status"],
            "ambiguous_multiple_role_appropriate_state_blocks",
        )

    def test_sinex_year_pivot_and_solution_validity_selection(self) -> None:
        self.assertEqual(parse_sinex_epoch("51:001:00000").year, 1951)
        self.assertEqual(parse_sinex_epoch("50:001:00000").year, 2050)
        content = """+SOLUTION/EPOCHS
*CODE PT SOLN T _DATA_START_ __DATA_END__ _MEAN_EPOCH_
 DYNG A    1 C 24:196:00000 24:198:86399 24:197:00000
 DYNG A    2 C 24:199:00000 24:202:86399 24:200:00000
-SOLUTION/EPOCHS
+SOLUTION/ESTIMATE
*INDEX TYPE__ CODE PT SOLN _REF_EPOCH__ UNIT S __ESTIMATED_VALUE____ _STD_DEV___
     1 STAX   DYNG A 1 24:197:00000 m 2 1.0 0.001
     2 STAY   DYNG A 1 24:197:00000 m 2 2.0 0.001
     3 STAZ   DYNG A 1 24:197:00000 m 2 3.0 0.001
     4 STAX   DYNG A 2 24:200:00000 m 2 4.0 0.001
     5 STAY   DYNG A 2 24:200:00000 m 2 5.0 0.001
     6 STAZ   DYNG A 2 24:200:00000 m 2 6.0 0.001
-SOLUTION/ESTIMATE
"""
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "reference.snx"
            path.write_text(content, encoding="ascii")
            parsed = parse_sinex(path)
            selection = select_sinex_reference(
                parsed.references_by_receiver["DYNG"],
                parse_sinex_epoch("24:199:43200"),
            )

        self.assertEqual(selection.status, "unique_time_valid_solution")
        self.assertEqual(selection.reference.solution, "2")
        self.assertEqual(selection.reference.ecef_m, (4.0, 5.0, 6.0))

    def test_reports_ecef_enu_distances_and_receiver_statistics(self) -> None:
        epoch = "2024-07-17 00:00:00.00"
        reference = (6378137.0, 0.0, 0.0)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            primary = root / "primary.trace"
            floating = root / "float.trace"
            crd = root / "reference.snx"
            primary.write_text(
                state_block("AR", rec_pos(epoch, "DYNG", 6378138.0, 2.0, 3.0))
            )
            floating.write_text(
                state_block("PPP", rec_pos(epoch, "DYNG", 6378137.5, 1.0, 1.0))
            )
            crd.write_text(sinex({"DYNG": reference}), encoding="ascii")

            report = audit_coordinate_integrity(primary, floating, crd)

        comparison = report["station_epoch_comparisons"][0]
        self.assertAlmostEqual(
            comparison["ar_vs_float"]["distance_3d_m"], math.sqrt(5.25)
        )
        self.assertEqual(
            comparison["primary_vs_sinex"]["ecef_delta_m"],
            {"x": 1.0, "y": 2.0, "z": 3.0},
        )
        self.assertEqual(
            comparison["primary_vs_sinex"]["enu_delta_m"],
            {"east": 2.0, "north": 3.0, "up": 1.0},
        )
        stats = report["statistics"]["by_receiver"]["DYNG"]
        self.assertEqual(stats["ar_vs_float_3d"]["count"], 1)
        self.assertAlmostEqual(stats["primary_vs_sinex_3d"]["maximum_m"], math.sqrt(14))
        self.assertEqual(report["screening"]["status"], "not_evaluated_no_thresholds_configured")
        self.assertNotIn("pass", report)

    def test_flags_jump_with_submitted_unverified_event_in_interval(self) -> None:
        epoch_0 = "2024-07-17 00:00:00.00"
        epoch_1 = "2024-07-17 00:00:30.00"
        reference = (6378137.0, 0.0, 0.0)
        primary_text = state_block("AR", rec_pos(epoch_0, "DYNG", *reference))
        primary_text += f"fixAndHoldAmbiguities: {epoch_1}\n"
        primary_text += (
            "      Applying:   -1 A(DYNG,G12,L1C) +1 A(DYNG,G25,L1C) = +4.00000\n"
        )
        primary_text += (
            "PPP_AR PSEUDOOBS_DESIGN rows=1 original_ambiguities=2 "
            "status=MATCHES_Z_TIMES_D\n"
            "PPP_AR PSEUDOOBS_NOISE rows=1 variance=1e-08 "
            "status=INDEPENDENT_DIAGONAL\n"
            "PPP_AR PSEUDOOBS_SUBMISSION rows=1 "
            "status=FILTER_CALL_RETURNED_SUBMITTED_UNVERIFIED\n"
        )
        primary_text += state_block(
            "AR", rec_pos(epoch_1, "DYNG", 6378137.3, 0.0, 0.0)
        )
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            primary = root / "primary.trace"
            floating = root / "float.trace"
            crd = root / "reference.snx"
            primary.write_text(primary_text, encoding="utf-8")
            floating.write_text(
                state_block(
                    "PPP",
                    rec_pos(epoch_0, "DYNG", *reference)
                    + rec_pos(epoch_1, "DYNG", *reference),
                ),
                encoding="utf-8",
            )
            crd.write_text(sinex({"DYNG": reference}), encoding="ascii")

            report = audit_coordinate_integrity(
                primary,
                floating,
                crd,
                Thresholds(ar_float_m=0.2, sinex_m=0.2, jump_m=0.2),
            )

        jump = report["primary_coordinate_jumps"][0]
        self.assertTrue(jump["exceeds_jump_threshold"])
        self.assertEqual(len(jump["integer_pseudoobservation_attempts_in_interval"]), 1)
        self.assertEqual(len(jump["receiver_attempts_in_interval"]), 1)
        self.assertEqual(
            len(jump["receiver_submitted_unverified_events_in_interval"]), 1
        )
        self.assertEqual(
            report["screening"]["exceedance_counts"][
                "submitted_unverified_interval_primary_epoch_jump_3d"
            ],
            1,
        )
        self.assertEqual(
            report["screening"]["status"],
            "configured_threshold_exceedance_observed",
        )

    def test_reports_missing_trace_and_sinex_data(self) -> None:
        epoch_0 = "2024-07-17 00:00:00.00"
        epoch_1 = "2024-07-17 00:00:30.00"
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            primary = root / "primary.trace"
            floating = root / "float.trace"
            crd = root / "reference.snx"
            primary.write_text(
                state_block("AR", rec_pos(epoch_0, "NONE", 1.0, 2.0, 3.0))
            )
            floating.write_text(
                state_block("PPP", rec_pos(epoch_1, "NONE", 1.0, 2.0, 3.0))
            )
            crd.write_text(sinex({"DYNG": (6378137.0, 0.0, 0.0)}), encoding="ascii")

            report = audit_coordinate_integrity(primary, floating, crd)

        missing = report["missing_data"]
        self.assertEqual(len(missing["missing_from_primary_trace"]), 1)
        self.assertEqual(len(missing["missing_from_float_trace"]), 1)
        self.assertEqual(missing["receivers_without_complete_sinex_reference"], ["NONE"])
        self.assertTrue(all(row["missing"] for row in report["station_epoch_comparisons"]))

    def test_cli_writes_valid_json_without_turning_screening_into_acceptance(self) -> None:
        epoch = "2024-07-17 00:00:00.00"
        reference = (6378137.0, 0.0, 0.0)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            primary = root / "primary.trace"
            floating = root / "float.trace"
            crd = root / "reference.snx"
            output = root / "audit.json"
            summary_output = root / "audit-summary.json"
            primary.write_text(state_block("AR", rec_pos(epoch, "DYNG", *reference)))
            floating.write_text(state_block("PPP", rec_pos(epoch, "DYNG", *reference)))
            crd.write_text(sinex({"DYNG": reference}), encoding="ascii")

            return_code = main(
                [str(primary), str(floating), str(crd), "--output", str(output)]
            )
            payload = json.loads(output.read_text(encoding="utf-8"))
            summary_return_code = main(
                [
                    str(primary),
                    str(floating),
                    str(crd),
                    "--output",
                    str(summary_output),
                    "--summary-only",
                ]
            )
            summary_payload = json.loads(
                summary_output.read_text(encoding="utf-8")
            )

        self.assertEqual(return_code, 0)
        self.assertEqual(payload["schema"], "GINAN_PPPAR_COORDINATE_INTEGRITY_AUDIT_V1")
        self.assertEqual(
            payload["screening"]["status"], "not_evaluated_no_thresholds_configured"
        )
        self.assertIn("not PPP-AR acceptance", payload["screening"]["scientific_interpretation"])
        self.assertEqual(summary_return_code, 0)
        self.assertTrue(summary_payload["schema"].endswith("_SUMMARY"))
        self.assertNotIn("station_epoch_comparisons", summary_payload)
        self.assertEqual(
            summary_payload["summary_only"]["omitted_record_counts"][
                "station_epoch_comparisons"
            ],
            1,
        )


if __name__ == "__main__":
    unittest.main()
