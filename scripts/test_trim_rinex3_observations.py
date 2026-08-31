from __future__ import annotations

import tempfile
import unittest
from datetime import datetime
from pathlib import Path

from trim_rinex3_observations import parse_epoch, trim_file


HEADER = """     3.05           OBSERVATION DATA    M                   RINEX VERSION / TYPE
  2024    07    17    00    00    0.0000000     GPS         TIME OF FIRST OBS
                                                            END OF HEADER
"""


class TrimRinex3ObservationsTest(unittest.TestCase):
    def test_parses_epoch(self) -> None:
        self.assertEqual(
            parse_epoch("> 2024 07 17 12 00 30.0000000  0  1\n"),
            datetime(2024, 7, 17, 12, 0, 30),
        )

    def test_trims_without_modifying_source(self) -> None:
        payload = (
            HEADER
            + "> 2024 07 17 11 59 30.0000000  0  1\nG01 before\n"
            + "> 2024 07 17 12 00  0.0000000  0  1\nG01 first\n"
            + "> 2024 07 17 12 00 30.0000000  0  1\nG01 second\n"
        )
        with tempfile.TemporaryDirectory() as temporary:
            source = Path(temporary) / "source.rnx"
            destination = Path(temporary) / "subset" / source.name
            source.write_text(payload, encoding="ascii")
            record = trim_file(
                source,
                destination,
                datetime(2024, 7, 17, 12, 0),
            )
            output = destination.read_text(encoding="ascii")
            self.assertEqual(source.read_text(encoding="ascii"), payload)
        self.assertNotIn("before", output)
        self.assertIn("G01 first", output)
        self.assertIn("G01 second", output)
        self.assertEqual(record["first_epoch"], "2024-07-17T12:00:00")


if __name__ == "__main__":
    unittest.main()
