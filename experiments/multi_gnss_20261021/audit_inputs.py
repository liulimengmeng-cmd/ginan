#!/usr/bin/env python3
"""Read-only GNSS input audit. Outputs are evidence, never a quality certificate.

Read the configured 180 station headers and only the first five minutes of
selected candidate stations. Never rewrite/download raw RINEX or products.
"""
from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re

SIGNALS = {"G": ("L1C", "L2W"), "E": ("L1C", "L5Q"), "C": ("L1P", "L5P")}
PRODUCT_NAMES = {
    "sp3": "WUM0MGXFIN_20241990000_01D_05M_ORB.SP3",
    "clk": "WUM0MGXFIN_20241990000_01D_30S_CLK.CLK",
    "bias": "WUM0MGXFIN_20241990000_01D_01D_OSB.BIA",
    "atx": "igs20_2303.atx",
    "metadata": "igs_satellite_metadata_2203_plus.snx",
}
CANDIDATES = ("MATE MAL2 AGGO HOB2 DARW DGAR FAIR KIRI NNOR IISC KARR REUN "
              "MAJU YEL2 KAT1 ABMF PIMO CUSV WUH2 HRAO SEYG HARB").split()


def sha256(path: Path) -> str:
    result = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            result.update(chunk)
    return result.hexdigest()


def epoch_value(value: str, open_end=False):
    if value == "0000:000:00000":
        return (9999, 999, 99999) if open_end else (0, 0, 0)
    return tuple(map(int, value.split(":")))


def read_header(path: Path):
    result = {"filename": path.name, "exists": path.is_file(), "signals": {}}
    if not result["exists"]:
        return result
    digest = hashlib.sha256()
    system = None
    with path.open("rb") as stream:
        for number, raw in enumerate(stream, 1):
            digest.update(raw)
            line = raw.decode("ascii", "replace").rstrip("\r\n")
            label = line[60:].strip()
            if label == "RINEX VERSION / TYPE":
                result["version"] = line[:9].strip()
            elif label == "MARKER NAME":
                result["marker_name"] = line[:60].strip()
            elif label == "ANT # / TYPE":
                result["antenna_type"] = line[20:40].rstrip()
            elif label == "APPROX POSITION XYZ":
                result["approx_xyz"] = list(map(float, line[:60].split()))
            elif label == "SYS / # / OBS TYPES":
                if line[0] != " ":
                    system = line[0]
                    result["signals"][system] = []
                result["signals"][system].extend(line[6:60].split())
            elif label == "INTERVAL":
                result["header_interval_seconds"] = float(line[:10])
            elif label == "END OF HEADER":
                result["header_lines"] = number
                break
        else:
            raise ValueError(f"No END OF HEADER: {path}")
    result["header_sha256"] = digest.hexdigest()
    result["header_target_support"] = {
        system: all(code in result["signals"].get(system, [])
                    for phase in pair for code in (phase, "C" + phase[1:]))
        for system, pair in SIGNALS.items()
    }
    result["size_bytes"] = path.stat().st_size
    return result


def window_audit(path: Path, header, whitelist):
    """Parse the first 5 minutes only; headers do not prove usable observations."""
    counts = Counter()
    satellites = {s: set() for s in SIGNALS}
    exact_pair_epochs = {s: Counter() for s in SIGNALS}
    epoch_times = []
    with path.open("r", encoding="ascii", errors="replace") as stream:
        for line in stream:
            if line[60:].strip() == "END OF HEADER":
                break
        for line in stream:
            if not line.startswith(">"):
                continue
            fields = line[1:].split()
            stamp = tuple(map(int, fields[:5])) + (float(fields[5]),)
            seconds = stamp[3] * 3600 + stamp[4] * 60 + stamp[5]
            if seconds > 300:
                break
            flag, nrecords = int(fields[6]), int(fields[7])
            if flag > 1:
                for _ in range(nrecords):
                    next(stream)
                continue
            epoch_times.append(stamp)
            for _ in range(nrecords):
                record = next(stream).rstrip("\r\n")
                sat, system = record[:3], record[0]
                if system not in SIGNALS or (system == "C" and sat not in whitelist):
                    continue
                codes = header["signals"].get(system, [])
                desired = [code for phase in SIGNALS[system]
                           for code in (phase, "C" + phase[1:])]
                if any(code not in codes for code in desired):
                    continue
                complete = True
                for code in desired:
                    index = codes.index(code)
                    field = record[3 + 16 * index:3 + 16 * index + 14].strip()
                    try:
                        complete &= bool(field) and float(field) != 0
                    except ValueError:
                        complete = False
                if complete:
                    counts[system] += 1
                    satellites[system].add(sat)
                    exact_pair_epochs[system][seconds] += 1
    return {"scope": "00:00:00 through 00:05:00 GPST, not the full day",
            "parsed_epochs": len(epoch_times), "first_epoch": epoch_times[0] if epoch_times else None,
            "last_epoch": epoch_times[-1] if epoch_times else None,
            "nonzero_code_and_phase_pair_observations": dict(counts),
            "satellites": {s: sorted(v) for s, v in satellites.items()},
            "minimum_pair_satellites_per_epoch": {
                s: min((exact_pair_epochs[s].get(t[3] * 3600 + t[4] * 60 + t[5], 0)
                        for t in epoch_times), default=0) for s in SIGNALS},
            "warning": "Raw signal presence only; no elevation, residual, antenna or cycle-slip acceptance applied."}


def metadata_audit(path: Path):
    block = ""
    identifiers = {}
    mappings = []
    with path.open("rb") as stream:
        for number, raw in enumerate(stream, 1):
            line = raw.decode("ascii", "replace").rstrip("\r\n")
            if line.startswith("+"):
                block = line.strip()
            if block == "+SATELLITE/IDENTIFIER":
                match = re.match(r" (C\d{3})\s+(\S+)\s+(\d+)\s+(\S+)", line)
                if match:
                    svn, cospar, norad, satellite_block = match.groups()
                    identifiers[svn] = {"svn": svn, "cospar": cospar, "norad": norad,
                        "block": satellite_block, "identifier_line": number,
                        "identifier_line_sha256": hashlib.sha256(raw).hexdigest()}
            elif block == "+SATELLITE/PRN":
                match = re.match(r" (C\d{3})\s+(\S+)\s+(\S+)\s+(C\d{2})", line)
                if match:
                    svn, start, stop, prn = match.groups()
                    identity = identifiers.get(svn, {})
                    if (identity.get("block", "").startswith("BDS-3")
                            and epoch_value(start) <= (2024, 199, 0)
                            and epoch_value(stop, True) >= (2024, 200, 0)):
                        mappings.append(dict(identity, prn=prn, valid_from=start, valid_until=stop,
                            prn_line=number, prn_line_sha256=hashlib.sha256(raw).hexdigest()))
    return sorted(mappings, key=lambda row: row["prn"])


def product_membership(root: Path):
    sp3 = set()
    with (root / PRODUCT_NAMES["sp3"]).open("r", encoding="ascii", errors="replace") as stream:
        for line in stream:
            if line.startswith("+ "):
                sp3.update(re.findall(r"[GREC]\d{2}", line))
            if line.startswith("* "):
                break
    clocks = set()
    clock_counts = Counter()
    with (root / PRODUCT_NAMES["clk"]).open("r", encoding="ascii", errors="replace") as stream:
        for line in stream:
            if line.startswith("AS "):
                sat = line[3:6]
                clocks.add(sat)
                clock_counts[sat] += 1
    bias = {}
    with (root / PRODUCT_NAMES["bias"]).open("r", encoding="ascii", errors="replace") as stream:
        for line in stream:
            match = re.match(r" OSB\s+\S+\s+([GREC]\d{2})\s+([CL]\d\w)\s+", line)
            if match:
                sat, signal = match.groups()
                bias.setdefault(sat, set()).add(signal)
    return {"sp3_satellites": sorted(sp3), "clk_satellites": sorted(clocks),
            "clk_record_count_by_satellite": dict(sorted(clock_counts.items())),
            "bias_codes_by_satellite": {s: sorted(v) for s, v in sorted(bias.items())},
            "note": "WUM OSB is available for later user comparisons; server base bsx_files is empty."}


def antex_audit(path: Path):
    blocks, antenna = [], None
    with path.open("r", encoding="ascii", errors="replace") as stream:
        for line in stream:
            label = line[60:].strip()
            if label == "START OF ANTENNA":
                antenna = {"frequencies": []}
            if antenna is None:
                continue
            if label == "TYPE / SERIAL NO":
                antenna.update(type=line[:20].rstrip(), serial=line[20:40].strip(), svn=line[40:50].strip())
            elif label in ("VALID FROM", "VALID UNTIL"):
                parts = line[:50].split()
                antenna["from" if label == "VALID FROM" else "until"] = tuple(map(int, parts[:5])) + (float(parts[5]),)
            elif label == "START OF FREQUENCY":
                antenna["frequencies"].append(line[:10].strip())
            elif label == "END OF ANTENNA":
                if (antenna.get("from", (0,)) <= (2024, 7, 17, 0, 0, 0)
                        and antenna.get("until", (9999,)) >= (2024, 7, 18, 0, 0, 0)):
                    blocks.append(antenna)
                antenna = None
    return blocks


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-root", type=Path, required=True)
    parser.add_argument("--station-config", type=Path, required=True)
    parser.add_argument("--checkpoint-manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    names = re.findall(r"^    - (.+\.rnx)$", args.station_config.read_text(), re.MULTILINE)
    if len(names) != len(set(names)) or len(names) != 180:
        raise ValueError("Expected original unique 180 service stations")
    products_root = args.input_root / "products"
    product_files = {key: dict(path=str(products_root / value), sha256=sha256(products_root / value),
                              size_bytes=(products_root / value).stat().st_size)
                     for key, value in PRODUCT_NAMES.items()}
    products = product_membership(products_root)
    meta = metadata_audit(products_root / PRODUCT_NAMES["metadata"])
    whitelist = [row for row in meta if row["prn"] in products["sp3_satellites"]
                 and row["prn"] in products["clk_satellites"]]
    atx = antex_audit(products_root / PRODUCT_NAMES["atx"])
    for row in whitelist:
        satellite = [b for b in atx if b.get("serial") == row["prn"] and b.get("svn") == row["svn"]]
        row["antex_L1P_L5P_supported"] = any(all(f in b["frequencies"] for f in ("C01", "C05")) for b in satellite)
    headers = [read_header(args.input_root / "data" / name) for name in names]
    antenna_models = {b["type"]: b for b in atx if not re.match(r"^[GRECJIS]\d{2}$", b.get("serial", ""))}
    for row in headers:
        model = antenna_models.get(row.get("antenna_type"), {})
        row["receiver_antex_target_frequency_support"] = {
            s: all(freq in model.get("frequencies", []) for freq in frequencies)
            for s, frequencies in {"G": ("G01", "G02"), "E": ("E01", "E05"), "C": ("C01", "C05")}.items()}
    by_station = {row["filename"][:4]: row for row in headers}
    satellite_set = {row["prn"] for row in whitelist}
    windows = {}
    for station in CANDIDATES:
        if station in by_station:
            header = by_station[station]
            windows[station] = window_audit(args.input_root / "data" / header["filename"], header, satellite_set)
    selected = ["MATE"]
    for station in CANDIDATES:
        if station == "MATE" or station not in by_station:
            continue
        if (all(by_station[station].get("header_target_support", {}).values())
                and all(by_station[station]["receiver_antex_target_frequency_support"].values())
                and all(windows[station]["minimum_pair_satellites_per_epoch"][s] >= 3 for s in SIGNALS)):
            selected.append(station)
        if len(selected) == 12:
            break
    checkpoint = json.loads(args.checkpoint_manifest.read_text())
    references = {"G": "MATE", "E": "MATE", "C": "MAL2"}
    reference_usable = {s: (by_station[ref]["header_target_support"][s]
                       and windows[ref]["minimum_pair_satellites_per_epoch"][s] > 0)
                       for s, ref in references.items()}
    summary = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "Input/header/first-5-minute readiness only, not full-day integrity or scientific acceptance",
        "source_station_config": {"path": str(args.station_config), "sha256": sha256(args.station_config)},
        "products": product_files, "product_membership": products,
        "configured_count": len(names), "missing_files": [h["filename"] for h in headers if not h["exists"]],
        "exact_header_support_station_counts": {s: sum(h.get("header_target_support", {}).get(s, False) for h in headers) for s in SIGNALS},
        "common_exact_header_support_station_count": sum(all(h.get("header_target_support", {}).values()) for h in headers),
        "common_header_and_receiver_antex_support_station_count": sum(
            all(h.get("header_target_support", {}).values())
            and all(h["receiver_antex_target_frequency_support"].values()) for h in headers),
        "signal_pairs": SIGNALS, "reference_stations": references,
        "reference_first_5min_nonzero_supported": reference_usable,
        "selected_12_station_candidates": selected,
        "candidate_rejected_antenna_support": [s for s in CANDIDATES if s in by_station
            and all(by_station[s].get("header_target_support", {}).values())
            and not all(by_station[s]["receiver_antex_target_frequency_support"].values())],
        "selected_complete": len(selected) == 12,
        "input_checkpoint": {k: checkpoint.get(k) for k in ("checkpoint_id", "epoch", "experiment_mode", "binary_sha256")},
        "input_checkpoint_valid_for_multi_system_hot_start": False,
        "multi_system_start_requirement": "Cold FLOAT from 00:00, all selected systems present; create a new G/E/C checkpoint only after validated FLOAT warmup.",
        "solver_started": False, "whitelist": whitelist, "station_headers": headers,
        "selected_window_observations": windows,
        "known_limits": ["Declared header signals do not establish all-day usable observations.",
                         "Exact L1C/L5Q excludes L1X/L5X-only Galileo receivers; no silent substitution.",
                         "WUM precise clock datum must still be checked in the observation model.",
                         "SINEX release 2203 is the locally supplied identity source, not a fresh global-generation catalogue.",
                         "12-station smoke tests cannot substitute for 180-station or all-day acceptance."]}
    (args.output / "input_readiness.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    (args.output / "bds3_whitelist.json").write_text(json.dumps({
        "epoch_day": "2024-07-17", "method": "SINEX IDENTIFIER block + effective PRN assignment + WUM SP3/CLK intersection",
        "source": product_files["metadata"], "satellites": whitelist}, indent=2), encoding="utf-8")
    lines = ["GPS/Galileo/BDS-3 read-only input audit", summary["scope"],
             f"Configured files: {len(names)}; missing: {len(summary['missing_files'])}",
             f"Exact code+phase header support: {summary['exact_header_support_station_counts']}",
             f"Common exact G/E/C header support: {summary['common_exact_header_support_station_count']}; with exact receiver ANTEX support: {summary['common_header_and_receiver_antex_support_station_count']}",
             f"References: {references}; first-5-minute raw pairs: {reference_usable}",
             f"12-station candidate network: {' '.join(selected)}", "BDS3 whitelist: " + " ".join(sorted(satellite_set)),
             "Original checkpoint experiment_mode: " + str(checkpoint.get("experiment_mode")),
             "NOT a multi-system hot start. No solver run was started.", *summary["known_limits"]]
    (args.output / "INPUT_READINESS_REPORT.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(json.dumps({key: summary[key] for key in ("exact_header_support_station_counts", "reference_first_5min_nonzero_supported", "selected_12_station_candidates", "selected_complete")}, indent=2))


if __name__ == "__main__":
    main()
