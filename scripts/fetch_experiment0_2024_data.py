#!/usr/bin/env python3
"""Fetch and hash the frozen 2024 Australian Experiment 0 inputs.

The Geoscience Australia API is requested with ``decompress=true`` and returns
plain RINEX text even though its transient object-lambda URL contains the
original ``.crx.gz`` object name.  Transient signed URLs are never persisted.

Dynamic precise products are taken from Wuhan University's anonymous IGS MGEX
archive using one internally consistent WUM final series.  Only orbit files are
fetched for the adjacent days; clock, ERP and OSB products are fetched for the
experiment day itself.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import shutil
import time
import urllib.parse
import urllib.request
from dataclasses import asdict, dataclass
from datetime import date, datetime, time as datetime_time, timedelta, timezone
from pathlib import Path
from typing import Iterable


GA_RINEX_API = "https://data.gnss.ga.gov.au/api/rinexFiles"
WUM_PRODUCTS = "ftp://igs.gnsswhu.cn/pub/gnss/products/mgex"
GA_PRODUCTS = "https://ga-gnss-products-v1.s3.amazonaws.com/public"
GFZ_KP_API = "https://kp.gfz.de/app/json/"
KYOTO_DST = "https://wdc.kugi.kyoto-u.ac.jp/dst_provisional/202405/dst2405.for.request"

MODEL_STATIONS = (
    "ARMC",
    "BATH",
    "HOB2",
    "MCHL",
    "MOBS",
    "STR1",
    "SYDN",
    "TID1",
    "WGGA",
    "YARR",
)
HELDOUT_STATIONS = ("STR2", "BALL", "PARK")
STATIONS = MODEL_STATIONS + HELDOUT_STATIONS

REQUIRED_GPS_OBSERVABLES = ("C1C", "L1C", "C2W", "L2W")

STATIC_PRODUCTS = (
    (
        "https://files.igs.org/pub/station/general/pcv_archive/igs20_2303.atx",
        "igs20_2303.atx",
        "IGS20_2303 antenna calibration used by WUM final products",
        False,
    ),
    (
        "https://files.igs.org/pub/station/general/igs_satellite_metadata.snx",
        "tables/igs_satellite_metadata.snx",
        "IGS satellite metadata SINEX",
        False,
    ),
    (
        "https://peanpod.s3.ap-southeast-2.amazonaws.com/aux/products/tables/igrf14coeffs.txt.gz",
        "tables/igrf14coeffs.txt",
        "IGRF14 geomagnetic coefficients",
        True,
    ),
    (
        "https://peanpod.s3.ap-southeast-2.amazonaws.com/aux/products/tables/DE436.1950.2050.gz",
        "tables/DE436.1950.2050",
        "JPL DE436 planetary ephemeris",
        True,
    ),
    (
        "https://peanpod.s3.ap-southeast-2.amazonaws.com/aux/products/tables/gpt_25.grd.gz",
        "tables/gpt_25.grd",
        "GPT2 5-degree troposphere grid",
        True,
    ),
    (
        "https://peanpod.s3.ap-southeast-2.amazonaws.com/aux/products/tables/OLOAD_GO.BLQ.gz",
        "tables/OLOAD_GO.BLQ",
        "ocean tide loading BLQ",
        True,
    ),
    (
        "https://peanpod.s3.ap-southeast-2.amazonaws.com/aux/products/tables/ALOAD_GO.BLQ.gz",
        "tables/ALOAD_GO.BLQ",
        "atmospheric tide loading BLQ",
        True,
    ),
    (
        "https://peanpod.s3.ap-southeast-2.amazonaws.com/aux/products/tables/opoleloadcoefcmcor.txt.gz",
        "tables/opoleloadcoefcmcor.txt",
        "ocean pole tide loading coefficients",
        True,
    ),
    (
        "https://peanpod.s3.ap-southeast-2.amazonaws.com/aux/products/tables/sat_yaw_bias_rate.snx.gz",
        "tables/sat_yaw_bias_rate.snx",
        "satellite yaw bias and rate model",
        True,
    ),
)

EXPERIMENT_DATES = {
    "quiet": date(2024, 5, 8),
    "storm": date(2024, 5, 11),
}

ELIGIBLE_GPS_SATELLITES = {
    "quiet": tuple(f"G{prn:02d}" for prn in range(2, 33) if prn != 10),
    "storm": tuple(f"G{prn:02d}" for prn in range(2, 33)),
}


@dataclass(frozen=True)
class FrozenFile:
    role: str
    experiment: str
    source: str
    relative_path: str
    size: int
    sha256: str
    source_metadata: dict[str, object]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def gps_week(day: date) -> int:
    return (day - date(1980, 1, 6)).days // 7


def day_of_year(day: date) -> int:
    return int(day.strftime("%j"))


def gps_day_of_week(day: date) -> int:
    return (day - date(1980, 1, 6)).days % 7


def query_ga_rinex(stations: Iterable[str], day: date, file_type: str) -> list[dict[str, object]]:
    start = datetime.combine(day, datetime_time.min, tzinfo=timezone.utc)
    end = start + timedelta(days=1) - timedelta(seconds=1)
    params = {
        "metadataStatus": "valid",
        "stationId": ",".join(stations),
        "fileType": file_type,
        "rinexVersion": 3,
        "filePeriod": "01D",
        "decompress": "true",
        "startDate": start.strftime("%Y-%m-%dT%H:%M:%SZ"),
        "endDate": end.strftime("%Y-%m-%dT%H:%M:%SZ"),
        "tenantId": "default",
    }
    url = GA_RINEX_API + "?" + urllib.parse.urlencode(params)
    with urllib.request.urlopen(url, timeout=120) as response:
        payload = json.load(response)
    if not isinstance(payload, list):
        raise RuntimeError(f"unexpected GA API response for {day}: {type(payload)}")
    return payload


def ga_output_filename(entry: dict[str, object]) -> str:
    path_name = Path(urllib.parse.urlsplit(str(entry["fileLocation"])).path).name
    for suffix in (".crx.gz", ".rnx.gz", ".gz"):
        if path_name.lower().endswith(suffix):
            path_name = path_name[: -len(suffix)]
            break
    if str(entry["fileType"]) == "obs":
        return path_name + ".rnx"
    return path_name if path_name.lower().endswith(".rnx") else path_name + ".rnx"


def source_metadata_without_url(entry: dict[str, object]) -> dict[str, object]:
    return {key: value for key, value in entry.items() if key != "fileLocation"}


def download_atomic(url: str, destination: Path, attempts: int = 4) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    partial = destination.with_name(destination.name + ".part")
    for attempt in range(1, attempts + 1):
        try:
            request = urllib.request.Request(url, headers={"User-Agent": "Ginan-Experiment0/1.0"})
            with urllib.request.urlopen(request, timeout=180) as response, partial.open("wb") as output:
                shutil.copyfileobj(response, output, length=1024 * 1024)
            os.replace(partial, destination)
            return
        except Exception:
            partial.unlink(missing_ok=True)
            if attempt == attempts:
                raise
            time.sleep(2**attempt)


def download_gzip_atomic(url: str, destination: Path, attempts: int = 4) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    compressed = destination.with_name(destination.name + ".gz.part")
    partial = destination.with_name(destination.name + ".part")
    for attempt in range(1, attempts + 1):
        try:
            request = urllib.request.Request(url, headers={"User-Agent": "Ginan-Experiment0/1.0"})
            with urllib.request.urlopen(request, timeout=180) as response, compressed.open("wb") as output:
                shutil.copyfileobj(response, output, length=1024 * 1024)
            with gzip.open(compressed, "rb") as source, partial.open("wb") as output:
                shutil.copyfileobj(source, output, length=1024 * 1024)
            os.replace(partial, destination)
            compressed.unlink(missing_ok=True)
            return
        except Exception:
            compressed.unlink(missing_ok=True)
            partial.unlink(missing_ok=True)
            if attempt == attempts:
                raise
            time.sleep(2**attempt)


def validate_rinex(path: Path) -> None:
    with path.open("rb") as stream:
        prefix = stream.read(256)
    if b"RINEX VERSION / TYPE" not in prefix:
        raise RuntimeError(f"GA response is not decompressed RINEX text: {path}")


def validate_gps_observables(path: Path) -> None:
    gps_observables: set[str] = set()
    with path.open("rt", encoding="ascii", errors="replace") as stream:
        for line in stream:
            if "SYS / # / OBS TYPES" in line and line.startswith("G"):
                gps_observables.update(line[7:60].split())
            if "END OF HEADER" in line:
                break
    missing = sorted(set(REQUIRED_GPS_OBSERVABLES) - gps_observables)
    if missing:
        raise RuntimeError(f"missing required GPS observables {missing}: {path}")


def validate_sinex(path: Path) -> None:
    with path.open("rb") as stream:
        prefix = stream.read(16)
    if not prefix.startswith(b"%=SNX"):
        raise RuntimeError(f"APREF response is not plain SINEX text: {path}")


def frozen_record(
    root: Path,
    path: Path,
    role: str,
    experiment: str,
    source: str,
    source_metadata: dict[str, object],
) -> FrozenFile:
    return FrozenFile(
        role=role,
        experiment=experiment,
        source=source,
        relative_path=path.relative_to(root).as_posix(),
        size=path.stat().st_size,
        sha256=sha256_file(path),
        source_metadata=source_metadata,
    )


def fetch_observations(root: Path, label: str, day: date) -> list[FrozenFile]:
    entries = query_ga_rinex(STATIONS, day, "obs")
    by_station: dict[str, list[dict[str, object]]] = {}
    for entry in entries:
        by_station.setdefault(str(entry["siteId"]).upper(), []).append(entry)
    missing = sorted(set(STATIONS) - set(by_station))
    duplicates = {station: values for station, values in by_station.items() if len(values) != 1}
    if missing or duplicates:
        raise RuntimeError(f"GA observation query is not one-to-one: missing={missing}, duplicates={list(duplicates)}")

    records: list[FrozenFile] = []
    for station in STATIONS:
        entry = by_station[station][0]
        destination = root / label / "data" / ga_output_filename(entry)
        if not destination.exists():
            download_atomic(str(entry["fileLocation"]), destination)
        validate_rinex(destination)
        validate_gps_observables(destination)
        records.append(
            frozen_record(
                root,
                destination,
                "GNSS RINEX 3 observation, 30 s",
                label,
                "Geoscience Australia GNSS Data Centre API",
                source_metadata_without_url(entry),
            )
        )
        print(f"{label}: observation {station}: {destination.name}", flush=True)
    return records


def fetch_broadcast_navigation(root: Path, label: str, day: date) -> FrozenFile:
    entries = query_ga_rinex(("BRDC",), day, "nav")
    if len(entries) != 1:
        raise RuntimeError(f"expected one BRDC file for {day}, received {len(entries)}")
    entry = entries[0]
    destination = root / label / "products" / ga_output_filename(entry)
    if not destination.exists():
        download_atomic(str(entry["fileLocation"]), destination)
    validate_rinex(destination)
    print(f"{label}: broadcast navigation: {destination.name}", flush=True)
    return frozen_record(
        root,
        destination,
        "broadcast navigation",
        label,
        "Geoscience Australia GNSS Data Centre API",
        source_metadata_without_url(entry),
    )


def fetch_wum_products(root: Path, label: str, day: date) -> list[FrozenFile]:
    target_doy = day_of_year(day)
    records: list[FrozenFile] = []
    specifications: list[tuple[int, str, str]] = []
    for offset in (-1, 0, 1):
        orbit_day = day + timedelta(days=offset)
        specifications.append((day_of_year(orbit_day), "01D_05M_ORB.SP3", "precise orbit"))
    specifications.extend(
        [
            (target_doy, "01D_30S_CLK.CLK", "precise clock, 30 s"),
            (target_doy, "01D_01D_ERP.ERP", "Earth rotation parameters"),
            (target_doy, "01D_01D_OSB.BIA", "WUM final Bias-SINEX OSB"),
        ]
    )

    for doy, suffix, role in specifications:
        product_day = date(day.year, 1, 1) + timedelta(days=doy - 1)
        product_week = gps_week(product_day)
        name = f"WUM0MGXFIN_{day.year}{doy:03d}0000_{suffix}"
        url = f"{WUM_PRODUCTS}/{product_week}/{name}.gz"
        destination = root / label / "products" / name
        if not destination.exists():
            download_gzip_atomic(url, destination)
        records.append(
            frozen_record(
                root,
                destination,
                role,
                label,
                "Wuhan University IGS MGEX anonymous FTP archive",
                {"url": url, "analysis_center": "WUM", "series": "MGXFIN"},
            )
        )
        print(f"{label}: product {role}: {name}", flush=True)
    return records


def fetch_static_products(root: Path) -> list[FrozenFile]:
    records: list[FrozenFile] = []
    for url, relative_name, role, is_gzip in STATIC_PRODUCTS:
        destination = root / "shared" / "products" / relative_name
        if not destination.exists():
            if is_gzip:
                download_gzip_atomic(url, destination)
            else:
                download_atomic(url, destination)
        records.append(
            frozen_record(
                root,
                destination,
                role,
                "shared",
                "official IGS archive" if "files.igs.org" in url else "Ginan auxiliary-data archive",
                {"url": url},
            )
        )
        print(f"shared: static product: {relative_name}", flush=True)
    return records


def fetch_apref_sinex(root: Path, label: str, day: date) -> FrozenFile:
    week = gps_week(day)
    day_of_week = gps_day_of_week(day)
    name = f"AUT{week}{day_of_week}.SNX"
    url = f"{GA_PRODUCTS}/{week}/{name}.gz"
    destination = root / label / "products" / name
    if not destination.exists():
        download_gzip_atomic(url, destination)
    validate_sinex(destination)
    print(f"{label}: APREF daily SINEX: {name}", flush=True)
    return frozen_record(
        root,
        destination,
        "APREF daily station coordinates and metadata",
        label,
        "Geoscience Australia GNSS Products S3",
        {"url": url, "gps_week": week, "gps_day_of_week": day_of_week},
    )


def fetch_space_weather(root: Path) -> list[FrozenFile]:
    directory = root / "selection" / "space_weather"
    kp_params = {
        "start": "2024-05-08T00:00:00Z",
        "end": "2024-05-12T00:00:00Z",
        "index": "Kp",
    }
    kp_url = GFZ_KP_API + "?" + urllib.parse.urlencode(kp_params)
    kp_path = directory / "gfz_kp_2024-05-08_2024-05-12.json"
    if not kp_path.exists():
        download_atomic(kp_url, kp_path)
    kp_payload = json.loads(kp_path.read_text(encoding="utf-8"))
    kp_by_time = dict(zip(kp_payload["datetime"], kp_payload["Kp"]))
    quiet_kp = [value for timestamp, value in kp_by_time.items() if timestamp.startswith("2024-05-08")]
    storm_kp = [value for timestamp, value in kp_by_time.items() if timestamp.startswith("2024-05-11")]
    if len(quiet_kp) != 8 or max(quiet_kp) > 2.0:
        raise RuntimeError(f"quiet-day Kp criterion failed: {quiet_kp}")
    if len(storm_kp) != 8 or max(storm_kp) < 8.0:
        raise RuntimeError(f"storm-day Kp criterion failed: {storm_kp}")

    dst_path = directory / "kyoto_dst_provisional_2024-05.txt"
    if not dst_path.exists():
        download_atomic(KYOTO_DST, dst_path)
    dst_lines = dst_path.read_text(encoding="ascii").splitlines()
    dst_by_day: dict[int, list[int]] = {}
    for line in dst_lines:
        if not line.startswith("DST2405*"):
            continue
        day_number = int(line[8:10])
        values = [int(line[position : position + 4]) for position in range(20, len(line), 4)]
        if len(values) < 24:
            raise RuntimeError(f"unexpected Kyoto Dst row: {line}")
        dst_by_day[day_number] = values[:24]
    quiet_dst = dst_by_day[8]
    storm_dst = dst_by_day[11]
    if min(quiet_dst) < -30 or max(storm_dst) > -100:
        raise RuntimeError(
            f"Dst selection criteria failed: quiet={quiet_dst}, storm={storm_dst}"
        )

    return [
        frozen_record(
            root,
            kp_path,
            "3-hour planetary Kp index used to freeze quiet/storm dates",
            "selection",
            "GFZ German Research Centre for Geosciences",
            {
                "url": kp_url,
                "quiet_day_max_kp": max(quiet_kp),
                "storm_day_max_kp": max(storm_kp),
            },
        ),
        frozen_record(
            root,
            dst_path,
            "hourly provisional Dst index used to freeze quiet/storm dates",
            "selection",
            "WDC for Geomagnetism, Kyoto",
            {
                "url": KYOTO_DST,
                "quiet_day_min_dst_nt": min(quiet_dst),
                "storm_day_min_dst_nt": min(storm_dst),
            },
        ),
    ]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output_root", type=Path)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    root = args.output_root.resolve()
    root.mkdir(parents=True, exist_ok=True)

    records: list[FrozenFile] = []
    records.extend(fetch_space_weather(root))
    records.extend(fetch_static_products(root))
    for label, day in EXPERIMENT_DATES.items():
        records.extend(fetch_observations(root, label, day))
        records.append(fetch_broadcast_navigation(root, label, day))
        records.extend(fetch_wum_products(root, label, day))
        records.append(fetch_apref_sinex(root, label, day))

    manifest = {
        "schema": "GINAN_EXPERIMENT0_2024_INPUT_MANIFEST_V1",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "dates": {label: day.isoformat() for label, day in EXPERIMENT_DATES.items()},
        "model_stations": list(MODEL_STATIONS),
        "heldout_stations": list(HELDOUT_STATIONS),
        "stations": list(STATIONS),
        "observation_policy": {
            "constellation": "GPS",
            "signals": list(REQUIRED_GPS_OBSERVABLES),
            "nominal_interval_seconds": 30,
        },
        "eligible_gps_satellites": {
            label: list(satellites) for label, satellites in ELIGIBLE_GPS_SATELLITES.items()
        },
        "excluded_gps_satellites": {
            "quiet": ["G01", "G10"],
            "storm": ["G01"],
        },
        "precise_product_policy": {
            "analysis_center": "Wuhan University",
            "series": "WUM0MGXFIN",
            "bias_apc_model": "IGS20_2303.ATX",
            "antenna_calibration": "shared/products/igs20_2303.atx",
        },
        "files": [asdict(record) for record in records],
    }
    manifest_path = root / "input_manifest.json"
    partial = manifest_path.with_name(manifest_path.name + ".part")
    partial.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(partial, manifest_path)
    print(f"manifest: {manifest_path}", flush=True)
    print(f"frozen files: {len(records)}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
