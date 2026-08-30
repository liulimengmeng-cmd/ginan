from __future__ import annotations

import io
import json
import sys
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from validate_experiment0_2024_split import (
    EXPECTED_HELDOUT,
    EXPECTED_TRAIN,
    build_audit,
    explicit_rnx_inputs,
    main,
)


def rinex_path(label: str, station: str, doy: int) -> str:
    source = "S" if station in {"BALL", "BATH", "WGGA"} else "R"
    return f"{label}/data/{station}00AUS_{source}_2024{doy:03d}0000_01D_30S_MO.rnx"


def model_yaml(label: str, stations: tuple[str, ...], doy: int) -> str:
    rows = "\n".join(f"      - {rinex_path(label, station, doy)}" for station in stations)
    return (
        "inputs:\n"
        "  include_yamls:\n"
        "    - common.yaml\n"
        "  gnss_observations:\n"
        "    rnx_inputs:\n"
        f"{rows}\n"
        "outputs:\n"
        "  metadata:\n"
        f"    config_description: {label}\n"
    )


def manifest_payload() -> dict[str, object]:
    files: list[dict[str, object]] = []
    for label, doy in (("quiet", 129), ("storm", 132)):
        for station in EXPECTED_TRAIN + EXPECTED_HELDOUT:
            files.append(
                {
                    "experiment": label,
                    "role": "GNSS RINEX 3 observation, 30 s",
                    "relative_path": rinex_path(label, station, doy),
                }
            )
    return {
        "schema": "GINAN_EXPERIMENT0_2024_INPUT_MANIFEST_V1",
        "dates": {"quiet": "2024-05-08", "storm": "2024-05-11"},
        "model_stations": list(EXPECTED_TRAIN),
        "heldout_stations": list(EXPECTED_HELDOUT),
        "stations": list(EXPECTED_TRAIN + EXPECTED_HELDOUT),
        "files": files,
    }


class Experiment0SplitAuditTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary_directory.name)
        self.manifest = self.root / "input_manifest.json"
        self.quiet = self.root / "quiet.yaml"
        self.storm = self.root / "storm.yaml"
        self.output = self.root / "audit.json"
        self.manifest.write_text(json.dumps(manifest_payload()), encoding="utf-8")
        self.quiet.write_text(model_yaml("quiet", EXPECTED_TRAIN, 129), encoding="utf-8")
        self.storm.write_text(model_yaml("storm", EXPECTED_TRAIN, 132), encoding="utf-8")

    def tearDown(self) -> None:
        self.temporary_directory.cleanup()

    def test_valid_split_and_required_summary_fields(self) -> None:
        audit = build_audit(self.manifest, self.quiet, self.storm)
        self.assertTrue(audit["valid"])
        self.assertEqual(audit["train"], sorted(EXPECTED_TRAIN))
        self.assertEqual(audit["heldout"], sorted(EXPECTED_HELDOUT))
        self.assertEqual(audit["intersection"], [])
        self.assertEqual(audit["heldout_input_count"], 0)
        self.assertTrue(audit["roles_identical"])
        self.assertEqual(audit["per_day"]["quiet"]["rnx_input_count"], 10)
        self.assertEqual(audit["per_day"]["storm"]["rnx_input_count"], 10)

    def test_heldout_observation_is_rejected_and_counted(self) -> None:
        contaminated = EXPECTED_TRAIN[:-1] + ("STR2",)
        self.quiet.write_text(model_yaml("quiet", contaminated, 129), encoding="utf-8")
        audit = build_audit(self.manifest, self.quiet, self.storm)
        self.assertFalse(audit["valid"])
        self.assertEqual(audit["heldout_input_count"], 1)
        self.assertEqual(audit["per_day"]["quiet"]["heldout_input_count"], 1)
        self.assertIn("STR2", audit["per_day"]["quiet"]["unexpected_stations"])
        self.assertFalse(audit["roles_identical"])

    def test_manifest_overlap_and_nonfrozen_role_are_rejected(self) -> None:
        payload = manifest_payload()
        payload["model_stations"][-1] = "STR2"
        self.manifest.write_text(json.dumps(payload), encoding="utf-8")
        audit = build_audit(self.manifest, self.quiet, self.storm)
        self.assertFalse(audit["valid"])
        self.assertEqual(audit["intersection"], ["STR2"])

    def test_quiet_and_storm_roles_must_match(self) -> None:
        changed = EXPECTED_TRAIN[:-1] + ("PARK",)
        self.storm.write_text(model_yaml("storm", changed, 132), encoding="utf-8")
        audit = build_audit(self.manifest, self.quiet, self.storm)
        self.assertFalse(audit["valid"])
        self.assertFalse(audit["roles_identical"])
        self.assertEqual(audit["heldout_input_count"], 1)

    def test_only_explicit_day_level_rnx_inputs_are_read(self) -> None:
        self.quiet.write_text(
            self.quiet.read_text(encoding="utf-8")
            + "unrelated:\n"
            + "  rnx_inputs:\n"
            + f"    - {rinex_path('quiet', 'STR2', 129)}\n",
            encoding="utf-8",
        )
        values = explicit_rnx_inputs(self.quiet)
        self.assertEqual(len(values), 10)
        audit = build_audit(self.manifest, self.quiet, self.storm)
        self.assertTrue(audit["valid"])

    def test_cli_writes_json_atomically_and_returns_failure_for_contamination(self) -> None:
        self.quiet.write_text(
            model_yaml("quiet", EXPECTED_TRAIN[:-1] + ("BALL",), 129),
            encoding="utf-8",
        )
        with redirect_stdout(io.StringIO()):
            status = main(
                [
                    str(self.manifest),
                    str(self.quiet),
                    str(self.storm),
                    str(self.output),
                ]
            )
        self.assertEqual(status, 1)
        payload = json.loads(self.output.read_text(encoding="utf-8"))
        self.assertFalse(payload["valid"])
        self.assertEqual(payload["heldout_input_count"], 1)
        self.assertFalse(list(self.root.glob("audit.json.*.tmp")))


if __name__ == "__main__":
    unittest.main()
