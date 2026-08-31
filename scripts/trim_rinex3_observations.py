#!/usr/bin/env python3
"""Create non-destructive RINEX 3 observation subsets from a start epoch."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
from datetime import datetime
from pathlib import Path


EPOCH_RE = re.compile(
    r"^>\s+(\d{4})\s+(\d{1,2})\s+(\d{1,2})\s+"
    r"(\d{1,2})\s+(\d{1,2})\s+([0-9.]+)"
)


def parse_epoch(line: str) -> datetime | None:
    match = EPOCH_RE.match(line)
    if not match:
        return None
    year, month, day, hour, minute = map(int, match.groups()[:5])
    seconds = float(match.group(6))
    whole_seconds = int(seconds)
    microseconds = round((seconds - whole_seconds) * 1_000_000)
    return datetime(year, month, day, hour, minute, whole_seconds, microseconds)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def first_observation_header(epoch: datetime, original: str) -> str:
    time_system = original[48:51] if len(original) >= 51 else "GPS"
    return (
        f"{epoch.year:6d}{epoch.month:6d}{epoch.day:6d}"
        f"{epoch.hour:6d}{epoch.minute:6d}"
        f"{epoch.second + epoch.microsecond / 1_000_000:13.7f}"
        f"     {time_system:<3}         TIME OF FIRST OBS\n"
    )


def trim_file(source: Path, destination: Path, start: datetime) -> dict[str, object]:
    if destination.exists():
        raise FileExistsError(f"refusing to overwrite {destination}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(destination.name + ".partial")
    if temporary.exists():
        raise FileExistsError(f"refusing to overwrite {temporary}")

    header: list[str] = []
    copied_lines = 0
    first_epoch: datetime | None = None
    try:
        with source.open("r", encoding="ascii", errors="strict", newline="") as input_stream:
            for line in input_stream:
                header.append(line)
                if line[60:].strip() == "END OF HEADER":
                    break
            else:
                raise ValueError(f"missing END OF HEADER in {source}")
            if not header or "OBSERVATION DATA" not in header[0]:
                raise ValueError(f"not a RINEX observation file: {source}")

            for index, line in enumerate(header):
                if line[60:].strip() == "TIME OF FIRST OBS":
                    header[index] = first_observation_header(start, line)
                    break

            with temporary.open(
                "w", encoding="ascii", errors="strict", newline=""
            ) as output_stream:
                output_stream.writelines(header)
                for line in input_stream:
                    epoch = parse_epoch(line)
                    if epoch is None or epoch < start:
                        continue
                    first_epoch = epoch
                    output_stream.write(line)
                    copied_lines += 1
                    break
                if first_epoch is None:
                    raise ValueError(f"no observation epoch at or after {start} in {source}")
                for line in input_stream:
                    output_stream.write(line)
                    copied_lines += 1
        os.replace(temporary, destination)
    except Exception:
        if temporary.exists():
            temporary.unlink()
        raise

    return {
        "source": str(source),
        "destination": str(destination),
        "requested_start": start.isoformat(),
        "first_epoch": first_epoch.isoformat(),
        "copied_data_line_count": copied_lines,
        "size_bytes": destination.stat().st_size,
        "sha256": sha256(destination),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("inputs", nargs="+", type=Path)
    parser.add_argument("--start", required=True, type=datetime.fromisoformat)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--manifest", type=Path)
    arguments = parser.parse_args()

    records = [
        trim_file(source, arguments.output_dir / source.name, arguments.start)
        for source in arguments.inputs
    ]
    report = {
        "schema": "RINEX3_OBSERVATION_TRIM_MANIFEST_V1",
        "file_count": len(records),
        "start": arguments.start.isoformat(),
        "files": records,
    }
    payload = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if arguments.manifest:
        arguments.manifest.write_text(payload, encoding="utf-8")
    print(payload, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
