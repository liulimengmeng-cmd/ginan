#!/usr/bin/env python3
"""Download verified official daily links and audit existing/new GNSS inputs.

No existing source is overwritten. Download workers <= 3. No solver is run.
Only explicit directory-listed links are downloaded; no guessed endless retry.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import time
import urllib.error
import urllib.request
from urllib.parse import urljoin
import zipfile

import audit_inputs

BRDC_URL = "https://igs.bkg.bund.de/root_ftp/IGS/BRDC/2024/199/"
OBS_URL = "https://igs.bkg.bund.de/root_ftp/IGS/obs/2024/199/"
TOOL_URL = "https://terras.gsi.go.jp/ja/crx2rnx/RNXCMP_4.2.0_Windows_mingw_64bit.zip"
TOOL_LICENSE = "https://terras.gsi.go.jp/ja/crx2rnx/LICENSE.txt"
TOOL_README = "https://terras.gsi.go.jp/ja/crx2rnx/README.txt"
TARGET_NAMES = ["BRDC00IGS_R_20241990000_01D_MN.rnx.gz", "BRD400DLR_S_20241990000_01D_MN.rnx.gz"]
TARGET_OBS = ["MAL200KEN_R_20241990000_01D_30S_MO.crx.gz",
              "HOB200AUS_R_20241990000_01D_30S_MO.crx.gz",
              "AGGO00ARG_R_20241990000_01D_30S_MO.crx.gz"]
SIGNALS = audit_inputs.SIGNALS


def dump(path, value):
    path.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")


def sha(path):
    return audit_inputs.sha256(path)


def get(url, destination):
    request = urllib.request.Request(url, headers={"User-Agent": "GNSS-research-input-audit/1.0"})
    result = {"url": url, "destination": str(destination), "started_utc": datetime.now(timezone.utc).isoformat()}
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists():
        result.update(status="EXISTING_NEW_DOWNLOAD_REUSED", size_bytes=destination.stat().st_size, sha256=sha(destination))
        return result
    part = destination.with_name(destination.name + ".part")
    if part.exists():
        raise RuntimeError(f"Unverified partial exists; no overwrite/resume assumption: {part}")
    started = time.monotonic()
    try:
        with urllib.request.urlopen(request, timeout=60) as response, part.open("xb") as output:
            result.update(http_status=response.status, effective_url=response.url,
                          headers=dict(response.headers.items()))
            while data := response.read(1024 * 1024):
                output.write(data)
        result.update(size_bytes=part.stat().st_size, sha256=sha(part))
        if part.name.endswith(".gz.part"):
            with gzip.open(part, "rb") as stream:
                for _ in iter(lambda: stream.read(1024 * 1024), b""):
                    pass
            result["gzip_crc_passed"] = True
        part.rename(destination)
        result.update(status="DOWNLOADED", elapsed_seconds=time.monotonic() - started)
    except Exception as error:
        result.update(status="FAILED", error=f"{type(error).__name__}: {error}", elapsed_seconds=time.monotonic() - started)
    return result


def epoch_seconds(fields):
    y, m, d, hour, minute = map(int, fields[:5])
    seconds = float(fields[5])
    if (y, m, d) != (2024, 7, 17):
        return None
    return hour * 3600 + minute * 60 + seconds


def observation_audit(path, whitelist):
    header = audit_inputs.read_header(path)
    systems = defaultdict(Counter)
    epoch_set = set()
    duplicates = 0
    first, last = None, None
    event_flags = Counter()
    body_hash = hashlib.sha256()
    invalid = []
    invalid_count = 0
    padded_satellite_ids = 0
    desired = {s: [code for phase in pair for code in (phase, "C" + phase[1:])]
               for s, pair in SIGNALS.items()}
    indexes = {s: [header["signals"].get(s, []).index(code) for code in wanted]
               for s, wanted in desired.items() if all(code in header["signals"].get(s, []) for code in wanted)}
    target_pair_epochs = {s: set() for s in SIGNALS}
    target_pair_satellites = {s: set() for s in SIGNALS}
    records = Counter()
    with path.open("rb") as raw:
        for line in raw:
            if line.decode("ascii", "replace")[60:].strip() == "END OF HEADER":
                break
        for rawline in raw:
            body_hash.update(rawline.rstrip(b"\r\n") + b"\n")
            line = rawline.decode("ascii", "replace").rstrip("\r\n")
            if not line.startswith(">"):
                invalid_count += 1
                if len(invalid) < 20:
                    invalid.append("Expected epoch: " + line[:40])
                if invalid_count > 20:
                    break
                continue
            try:
                # RINEX3 epoch fields are fixed-width. Special event records
                # can legally have no date/time (e.g. HKWS trailing flag 4).
                flag, nrecords = int(line[31:32]), int(line[32:35])
                stamp = epoch_seconds(line[2:29].split()) if line[2:29].strip() else None
            except (ValueError, IndexError):
                invalid_count += 1
                if len(invalid) < 20:
                    invalid.append("Malformed epoch: " + line[:80])
                continue
            event_flags[str(flag)] += 1
            valid_epoch = flag <= 1 and stamp is not None and 0 <= stamp < 86400
            if valid_epoch:
                duplicates += stamp in epoch_set
                epoch_set.add(stamp)
                first = stamp if first is None else min(first, stamp)
                last = stamp if last is None else max(last, stamp)
            for _ in range(nrecords):
                rawrecord = next(raw, b"")
                body_hash.update(rawrecord.rstrip(b"\r\n") + b"\n")
                if not rawrecord:
                    invalid_count += 1
                    if len(invalid) < 20:
                        invalid.append("Truncated record block")
                    break
                if flag > 1:
                    continue
                record = rawrecord.decode("ascii", "replace").rstrip("\r\n")
                system = record[0:1]
                # A1,I2 permits blank padding for one-digit satellite numbers.
                # Do not turn valid 'G 6' into a corrupt record or drop G06.
                try:
                    prn = int(record[1:3])
                    satellite = system + f"{prn:02d}"
                    valid_satellite = system in "GRECJIS" and 0 < prn < 100
                except ValueError:
                    valid_satellite = False
                if not valid_satellite:
                    invalid_count += 1
                    if len(invalid) < 20:
                        invalid.append("Unexpected satellite/continuation record: " + record[:20])
                    continue
                padded_satellite_ids += record[1:3].startswith(" ")
                records[system] += 1
                systems[system][satellite] += 1
                if system not in indexes or (system == "C" and satellite not in whitelist):
                    continue
                values = [record[3 + 16 * i:3 + 16 * i + 14].strip() for i in indexes[system]]
                try:
                    complete = all(v and math.isfinite(float(v)) and float(v) != 0 for v in values)
                except ValueError:
                    complete = False
                if complete:
                    records[system + "_target_code_phase_pair"] += 1
                    target_pair_satellites[system].add(satellite)
                    if valid_epoch:
                        target_pair_epochs[system].add(stamp)
    expected = set(range(0, 86400, 30))
    result = {"path": str(path), "sha256": sha(path), "size_bytes": path.stat().st_size,
        "header": header, "epoch_count": len(epoch_set), "duplicate_epochs": duplicates,
        "first_seconds_of_day": first, "last_seconds_of_day": last,
        "missing_30s_epochs": sorted(expected - epoch_set),
        "off_grid_epoch_count": len(epoch_set - expected), "event_flags": dict(event_flags),
        "full_2880_epoch_grid": epoch_set == expected and duplicates == 0,
        "all_observation_records_parseable": invalid_count == 0,
        "blank_padded_satellite_id_records": padded_satellite_ids,
        "actual_observation_records_by_system": {s: records[s] for s in sorted(systems)},
        "actual_satellites_by_system": {s: sorted(v) for s, v in sorted(systems.items())},
        "target_nonzero_code_phase_pair_counts": {s: records[s + "_target_code_phase_pair"] for s in SIGNALS},
        "target_nonzero_pair_epoch_counts": {s: len(v) for s, v in target_pair_epochs.items()},
        "target_nonzero_pair_satellites": {s: sorted(v) for s, v in target_pair_satellites.items()},
        "normalized_observation_body_sha256": body_hash.hexdigest(), "parse_errors": invalid,
        "parse_error_count": invalid_count,
        "scientific_quality_accepted": False}
    return result


def nav_audit(path):
    records = Counter()
    satellites = defaultdict(set)
    epochs = defaultdict(list)
    message_types = Counter()
    with path.open("r", encoding="ascii", errors="replace") as stream:
        version = next(stream)[:9].strip()
        for line in stream:
            if line[60:].strip() == "END OF HEADER":
                break
        for line in stream:
            if line.startswith(">"):
                fields = line.split()
                if len(fields) >= 4:
                    message_types[" ".join(fields[1:4])] += 1
            if re.match(r"^[GRECJIS]\d{2} \d{4}", line):
                satellite = line[:3]
                system = satellite[0]
                records[system] += 1
                satellites[system].add(satellite)
                epochs[system].append(line[4:23].strip())
    return {"path": str(path), "sha256": sha(path), "size_bytes": path.stat().st_size,
        "rinex_version": version, "ephemeris_records_by_system": dict(records),
        "satellites_by_system": {s: sorted(v) for s, v in sorted(satellites.items())},
        "native_epoch_range_by_system": {s: [min(v), max(v)] for s, v in sorted(epochs.items())},
        "rinex4_message_type_counts": dict(message_types),
        "BDS_CNAV_explicit_message_records": sum(n for key, n in message_types.items() if key.startswith("EPH C") and any(t in key for t in ("CNV1", "CNV2", "CNV3"))),
        "G_E_C_present": all(records[s] > 0 for s in "GEC"),
        "warning": "Mixed RINEX3 presence does not certify BDS CNAV. Native nav time fields are not normalized to GPST here.",
        "scientific_quality_accepted": False}


def precise_products_audit(root, whitelist):
    files = audit_inputs.PRODUCT_NAMES
    result = {"scope": "Record identities/grids/finite values, not external accuracy or datum acceptance"}
    for kind, expected_step in (("sp3", 300), ("clk", 30)):
        path = root / files[kind]
        times = defaultdict(set)
        duplicates = Counter()
        invalid = Counter()
        current = None
        with path.open("r", encoding="ascii", errors="replace") as stream:
            for line in stream:
                if kind == "sp3":
                    if line.startswith("* "):
                        current = epoch_seconds(line[1:].split())
                    if not line.startswith("P") or not re.match(r"^[GREC]\d{2}$", line[1:4]):
                        continue
                    satellite = line[1:4]
                    numbers = [line[4 + i * 14:18 + i * 14].strip() for i in range(3)]
                    stamp = current
                else:
                    if not line.startswith("AS "):
                        continue
                    fields = line.split()
                    satellite = fields[1]
                    stamp = epoch_seconds(fields[2:8])
                    numbers = fields[9:10]
                if stamp in times[satellite]:
                    duplicates[satellite] += 1
                times[satellite].add(stamp)
                try:
                    valid = stamp is not None and all(math.isfinite(float(v)) and abs(float(v)) < 999999 for v in numbers)
                except ValueError:
                    valid = False
                if not valid:
                    invalid[satellite] += 1
        expected = set(range(0, 86400, expected_step))
        by_system = Counter(s[0] for s in times)
        result[kind] = {"path": str(path), "sha256": sha(path), "size_bytes": path.stat().st_size,
            "satellite_counts": dict(by_system), "total_satellites": len(times),
            "expected_records_per_satellite": len(expected),
            "all_satellites_complete_grid": all(t == expected for t in times.values()),
            "duplicate_records": dict(duplicates), "invalid_records": dict(invalid),
            "per_satellite": {s: {"records": len(t), "missing_count": len(expected - t),
                                  "off_grid_count": len(t - expected)} for s, t in sorted(times.items())}}
    path = root / files["bias"]
    intervals = defaultdict(list)
    invalid = []
    with path.open("r", encoding="ascii", errors="replace") as stream:
        for number, line in enumerate(stream, 1):
            match = re.match(r" OSB\s+\S+\s+([GRECJ]\d{2})\s+([CL]\d\w)\s+(\d{4}:\d{3}:\d{5})\s+(\d{4}:\d{3}:\d{5})\s+\w+\s+([\d.+Ee-]+)", line)
            if not match:
                continue
            satellite, code, start, stop, value = match.groups()
            if not math.isfinite(float(value)):
                invalid.append(number)
            intervals[(satellite, code)].append((start, stop))
    coverage = {}
    for (satellite, code), spans in intervals.items():
        system = satellite[0]
        if system not in SIGNALS or code not in SIGNALS[system] or (system == "C" and satellite not in whitelist):
            continue
        covered = set()
        for start, stop in spans:
            if start.startswith("2024:199:"):
                a = int(start[-5:])
                b = int(stop[-5:]) if stop.startswith("2024:199:") else 86400 if stop.startswith("2024:200:") else a
                covered.update(range(a, b, 30))
        coverage[satellite + ":" + code] = {"covered_30s_ticks": len(covered), "complete_2880": len(covered) == 2880}
    result["bias"] = {"path": str(path), "sha256": sha(path), "size_bytes": path.stat().st_size,
        "satellite_counts": dict(Counter(s[0] for s in {key[0] for key in intervals})),
        "invalid_nonfinite_line_numbers": invalid, "target_phase_coverage": coverage,
        "warning": "Availability/coverage of OSB is not proof of its accuracy or service/user datum compatibility"}
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--destination", type=Path, required=True)
    parser.add_argument("--existing-input-root", type=Path, required=True)
    parser.add_argument("--audit", type=Path, required=True)
    parser.add_argument("--download-only", action="store_true")
    parser.add_argument("--audit-only", action="store_true")
    args = parser.parse_args()
    root = args.destination.resolve()
    root.mkdir(parents=True, exist_ok=True)
    audit = json.loads(args.audit.read_text())
    whitelist = {r["prn"] for r in audit["whitelist"]}
    manifest_path = root / "download_manifest.json"
    manifest = json.loads(manifest_path.read_text()) if manifest_path.exists() else {
        "created_utc": datetime.now(timezone.utc).isoformat(), "date": "2024-07-17", "concurrency": 3,
        "source_policy": "Official public directory-listed files only; no auth bypass; preserve originals",
        "downloads": [], "solver_started": False}
    if not args.audit_only:
        jobs = []
        for base, names, directory in ((BRDC_URL, TARGET_NAMES, "navigation"), (OBS_URL, TARGET_OBS, "observations")):
            index_path = root / "source_indexes" / (directory + ".html")
            response = get(base, index_path)
            manifest["downloads"].append(response)
            dump(manifest_path, manifest)
            if response["status"] == "FAILED":
                continue
            links = set(re.findall(r'href="([^"]+)"', index_path.read_text(errors="replace")))
            for name in names:
                if name not in links:
                    manifest["downloads"].append({"url": urljoin(base, name), "status": "NOT_IN_OFFICIAL_DIRECTORY"})
                    continue
                jobs.append((urljoin(base, name), root / directory / name))
        jobs += [(TOOL_URL, root / "tools" / TOOL_URL.rsplit("/", 1)[1]),
                 (TOOL_LICENSE, root / "tools/LICENSE.txt"), (TOOL_README, root / "tools/README.txt")]
        with ThreadPoolExecutor(max_workers=3) as pool:
            futures = [pool.submit(get, *job) for job in jobs]
            for future in as_completed(futures):
                response = future.result()
                manifest["downloads"].append(response)
                dump(manifest_path, manifest)
                print(json.dumps({k: response.get(k) for k in ("status", "url", "size_bytes", "sha256", "error")}), flush=True)
    nav_paths = []
    for compressed in sorted((root / "navigation").glob("*.gz")):
        output = compressed.with_suffix("")
        if not output.exists():
            with gzip.open(compressed, "rb") as source, output.open("xb") as destination:
                shutil.copyfileobj(source, destination)
        nav_paths.append(output)
    tool_zip = root / "tools" / TOOL_URL.rsplit("/", 1)[1]
    converter = None
    if tool_zip.exists():
        unpacked = root / "tools/RNXCMP_4.2.0"
        unpacked.mkdir(exist_ok=True)
        with zipfile.ZipFile(tool_zip) as archive:
            for member in archive.infolist():
                candidate = (unpacked / member.filename).resolve()
                if not candidate.is_relative_to(unpacked.resolve()):
                    raise RuntimeError("Unsafe RNXCMP zip member")
            archive.extractall(unpacked)
        converters = list(unpacked.rglob("crx2rnx.exe")) + list(unpacked.rglob("CRX2RNX.exe"))
        converter = converters[0] if converters else None
    converted = []
    for compressed in sorted((root / "observations").glob("*.crx.gz")):
        crx = compressed.with_suffix("")
        rnx = crx.with_suffix(".rnx")
        if not crx.exists():
            with gzip.open(compressed, "rb") as source, crx.open("xb") as destination:
                shutil.copyfileobj(source, destination)
        receipt = {"compressed": str(compressed), "crx_sha256": sha(crx), "rnx": str(rnx), "converter": str(converter)}
        if converter and not rnx.exists():
            part = rnx.with_name(rnx.name + ".part")
            with crx.open("rb") as source, part.open("xb") as output:
                process = subprocess.run([str(converter), "-"], stdin=source, stdout=output, stderr=subprocess.PIPE)
            receipt.update(exit_code=process.returncode, stderr=process.stderr.decode(errors="replace"))
            if process.returncode == 0:
                part.rename(rnx)
        if rnx.exists():
            receipt.update(rnx_sha256=sha(rnx), size_bytes=rnx.stat().st_size)
        converted.append(receipt)
    dump(root / "conversion_manifest.json", converted)
    manifest["navigation_audit"] = [nav_audit(path) for path in nav_paths]
    dump(manifest_path, manifest)
    if args.download_only:
        return
    evidence = {"scanner_contract": "rinex_fixed_width_space_padded_prn_and_undated_event_v2",
                "auditor_source_sha256": sha(Path(__file__)),
                "scientific_quality_accepted": False, "precise_products": precise_products_audit(args.existing_input_root / "products", whitelist),
                "existing_observations": [], "new_observations": []}
    dump(root / "validation_progress.json", {"state": "AUDITING_EXISTING_OBSERVATIONS", "total": len(audit["station_headers"]), "completed": 0})
    jobs = [(args.existing_input_root / "data" / row["filename"], "existing_observations") for row in audit["station_headers"]]
    jobs += [(Path(r["rnx"]), "new_observations") for r in converted if Path(r["rnx"]).exists()]
    with ThreadPoolExecutor(max_workers=3) as pool:
        futures = {pool.submit(observation_audit, path, whitelist): (path, bucket) for path, bucket in jobs}
        for future in as_completed(futures):
            path, bucket = futures[future]
            try:
                result = future.result()
            except Exception as error:
                result = {"path": str(path), "error": f"{type(error).__name__}: {error}", "full_2880_epoch_grid": False}
            evidence[bucket].append(result)
            complete = sum(map(len, (evidence["existing_observations"], evidence["new_observations"])))
            dump(root / "validation_progress.json", {"state": "AUDITING_OBSERVATIONS", "completed": complete, "total": len(jobs), "latest": str(path)})
            dump(root / "validation_results.json", evidence)
            if complete % 10 == 0 or bucket == "new_observations":
                print(f"Audited {complete}/{len(jobs)}; latest={path.name} epochs={result.get('epoch_count')} full={result.get('full_2880_epoch_grid')}", flush=True)
    old_by_name = {Path(row["path"]).name: row for row in evidence["existing_observations"]}
    comparisons = []
    for row in evidence["new_observations"]:
        old = old_by_name.get(Path(row["path"]).name)
        comparisons.append({"station": Path(row["path"]).name[:9], "new": row["path"], "old": old.get("path") if old else None,
            "complete_file_sha_equal": old is not None and row.get("sha256") == old.get("sha256"),
            "normalized_observation_body_sha_equal": old is not None and row.get("normalized_observation_body_sha256") == old.get("normalized_observation_body_sha256"),
            "actual_target_pair_counts_equal": old is not None and row.get("target_nonzero_code_phase_pair_counts") == old.get("target_nonzero_code_phase_pair_counts")})
    evidence["downloaded_vs_existing_comparison"] = comparisons
    evidence["summary"] = {
        "existing_total": len(evidence["existing_observations"]),
        "existing_complete_2880_grid": sum(r["full_2880_epoch_grid"] for r in evidence["existing_observations"]),
        "existing_all_records_parseable": sum(r.get("all_observation_records_parseable", False) for r in evidence["existing_observations"]),
        "existing_with_actual_GEC_observations": sum(all(r.get("actual_observation_records_by_system", {}).get(s, 0) > 0 for s in "GEC") for r in evidence["existing_observations"]),
        "existing_with_exact_GEC_target_pairs": sum(all(r.get("target_nonzero_code_phase_pair_counts", {}).get(s, 0) > 0 for s in "GEC") for r in evidence["existing_observations"]),
        "downloaded_observation_count": len(evidence["new_observations"]),
        "downloaded_complete_2880_grid": sum(r["full_2880_epoch_grid"] for r in evidence["new_observations"]),
        "unresolved_partial_or_invalid_original_files": [r["path"] for r in evidence["existing_observations"]
            if not r["full_2880_epoch_grid"] or not r.get("all_observation_records_parseable", False)],
        "scientific_quality_accepted": False, "solver_started": False}
    dump(root / "validation_results.json", evidence)
    dump(root / "validation_progress.json", {"state": "COMPLETED_RECORD_AUDIT_NOT_SCIENTIFIC_ACCEPTANCE", **evidence["summary"]})
    print(json.dumps(evidence["summary"], indent=2), flush=True)


if __name__ == "__main__":
    main()
