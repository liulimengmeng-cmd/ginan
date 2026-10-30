#!/usr/bin/env python3
"""Read-only evidence audit of real A--G stage outputs; never launches PEA.

File integrity, finite products and a PSD product covariance are necessary
diagnostics, not scientific validation. In particular this script never emits
scientific_acceptance=true or an execution/build approval. Only one epoch's
product rows and covariance entries are retained at a time.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import re

import numpy as np

LETTERS = {"gps": "G", "gal": "E", "bds": "C", "GPS": "G", "GAL": "E", "BDS": "C"}
NAMES = {"G": "GPS", "E": "GAL", "C": "BDS"}
SIGNALS = {"G": ("L1C", "L2W"), "E": ("L1C", "L5Q"), "C": ("L1P", "L5P")}
NUMBERS = ("clock_m", "clock_sigma_m", "phase_m", "phase_sigma_m",
           "clock_phase_covariance_m2", "correction_m", "correction_sigma_m")
FLAGS = ("numeric_valid", "branch_valid", "continuity_valid", "ppp_usable", "pppar_usable",
         "integer_valid", "dual_frequency_ar_valid", "ar_valid", "same_segment_continuous",
         "requires_user_ambiguity_reset")
REQUIRED_AR = ("numeric_valid", "branch_valid", "integer_valid", "dual_frequency_ar_valid", "ar_valid")
COV_COLUMNS = ("gpst_seconds", "solution", "row_satellite", "row_parameter", "row_observable",
               "column_satellite", "column_parameter", "column_observable", "covariance_m2")
DATE = r"\d{4}-\d\d-\d\d[ T]\d\d:\d\d:\d\d(?:\.\d+)?"
EPOCH_RE = re.compile(r"Processed epoch #(\d+).*?GPS time:\s*(\d+)\s+([\d.]+).*?(" + DATE + r")")
SPP_OBS_RE = re.compile(r"sppos:\s*time=(" + DATE + r"),\s*nObs=(\d+)")
SPP_STATUS_RE = re.compile(r"(?:SPP|sppos).*?\bstatus=(SINGLE_X|SINGLE|NONE)\b")


def unavailable(reason):
    return {"status": "UNAVAILABLE", "reason": reason}


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def finite(value):
    try:
        result = float(value)
        return result if math.isfinite(result) else None
    except (ValueError, TypeError):
        return None


def gpst(value):
    number = finite(value)
    if number is not None:
        return number
    try:
        text = str(value).replace("Z", "+00:00")
        # Python3.10 (Ubuntu22) accepts 3/6 fractional digits, whereas newer
        # Python also accepts PEA's .00/.50. Normalize supported precision
        # explicitly instead of changing the timestamp or treating it as0.
        text = re.sub(r"\.(\d{1,6})(?=(?:[+-]\d\d:\d\d)?$)",
                      lambda match: "." + match[1].ljust(6, "0"), text)
        date = datetime.fromisoformat(text)
        if date.tzinfo is None:
            date = date.replace(tzinfo=timezone.utc)
        # Calendar strings from PEA are already GPST, not UTC observations.
        return (date - datetime(1980, 1, 6, tzinfo=timezone.utc)).total_seconds()
    except (ValueError, TypeError):
        return None


def duration_seconds(value):
    """PEA GTime durations use HH:MM:SS, not a leading decimal hour.

    Older numeric-second records remain supported. Malformed/missing tokens
    return None, never a fabricated zero-second sample.
    """
    match = re.fullmatch(r"(\d+):([0-5]\d):([0-5]\d(?:\.\d+)?)", str(value))
    if match:
        return 3600 * int(match[1]) + 60 * int(match[2]) + float(match[3])
    if re.fullmatch(r"[+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?", str(value)):
        result = finite(value)
        return result if result is not None and result >= 0 else None
    return None


def csv_epochs(path):
    """Yield one complete epoch. Non-monotonic input is rejected, not sorted in RAM."""
    with Path(path).open(newline="", encoding="utf-8-sig") as stream:
        reader = csv.DictReader(stream)
        if not reader.fieldnames or len(set(reader.fieldnames)) != len(reader.fieldnames):
            raise ValueError("MISSING_OR_DUPLICATE_CSV_HEADER")
        if "gpst_seconds" not in reader.fieldnames:
            raise ValueError("CSV_EPOCH_COLUMN_MISSING")
        current = None
        rows = []
        for line, row in enumerate(reader, 2):
            if None in row or any(value is None for value in row.values()):
                raise ValueError(f"CSV_WIDTH_MISMATCH:{line}")
            epoch = finite(row["gpst_seconds"])
            if epoch is None:
                raise ValueError(f"CSV_NONFINITE_EPOCH:{line}")
            if current is not None and epoch < current:
                raise ValueError(f"CSV_EPOCH_ORDER_REVERSED:{line}")
            if current is not None and epoch != current:
                yield current, rows
                rows = []
            current = epoch
            rows.append(row)
        if current is not None:
            yield current, rows


def audit_logs(root):
    paths = []
    if (root / "runner.log").is_file():
        paths.append(root / "runner.log")
    outputs = root / "outputs"
    if outputs.is_dir():
        paths.extend(sorted(p for p in outputs.rglob("*") if p.is_file() and p.suffix.lower() == ".trace"))
    epochs, obs, statuses, failures = {}, {}, {}, set()
    times, endings = {}, []
    invalid_epoch_timestamps = 0
    for path in paths:
        with path.open(errors="replace") as stream:
            for line in stream:
                match = EPOCH_RE.search(line)
                if match:
                    timestamp = gpst(match[4])
                    if timestamp is None:
                        invalid_epoch_timestamps += 1
                        continue
                    key = (match[2], match[3])
                    epochs[key] = {"epoch_number": int(match[1]), "gpst_seconds": timestamp}
                    duration = re.search(r"epoch took\s+([^\s)]+)", line)
                    seconds = duration_seconds(duration[1]) if duration else None
                    if seconds is not None:
                        times[key] = seconds
                match = SPP_OBS_RE.search(line)
                if match:
                    # Receiver TRACE filename is part of this source identity.
                    # Do not call these raw RINEX records or accepted PPP rows.
                    obs[(str(path), match[1])] = int(match[2])
                match = SPP_STATUS_RE.search(line)
                if match:
                    statuses[(str(path), line.strip())] = match[1]
                if "SPP failed for " in line:
                    failures.add(line[line.index("SPP failed for "):].strip())
                if "as end epoch" in line and "has been reached" in line:
                    endings.append(line.strip())
    timing = unavailable("No parseable explicit processed-epoch duration records; missing/malformed values are not zero")
    if times:
        values = list(times.values())
        timing = {"status": "AVAILABLE", "sample_count": len(times), "unit": "seconds",
                  **{name: float(np.percentile(values, quantile)) for name, quantile in
                     (("p50", 50), ("p95", 95), ("p99", 99))}, "maximum": max(values)}
    return {"files_scanned": [str(p) for p in paths],
            "invalid_processed_epoch_calendar_records": invalid_epoch_timestamps,
            "processed_epochs": {"status": "AVAILABLE", "count": len(epochs),
                "first_gpst_seconds": min(e["gpst_seconds"] for e in epochs.values()),
                "last_gpst_seconds": max(e["gpst_seconds"] for e in epochs.values())} if epochs else
                unavailable("No explicit Processed epoch records; product epochs are reported separately"),
            "spp_input_observations": {"status": "AVAILABLE", "records": len(obs), "sum_nObs": sum(obs.values()),
                "scope": "nObs reported at sppos entry, not accepted/independent measurement count"} if obs else
                unavailable("Trace has no sppos time/nObs evidence"),
            "spp_single": {"status": "AVAILABLE", "explicit_status_records": dict(Counter(statuses.values())),
                "scope": "Explicit SPP status records only; SINGLE_X is not SINGLE"} if statuses else
                unavailable("Current spp.cpp normally does not print a successful status; positions or absence of failures do not prove SINGLE"),
            "spp_failure_records": len(failures), "end_epoch_markers": endings[-3:], "epoch_timing": timing}


def verify_bundle(directory):
    errors, hashes = [], {}
    manifest = {}
    try:
        if directory.name.endswith(".pending"):
            raise ValueError("PENDING_NOT_COMMITTED")
        manifest = json.loads((directory / "COMMIT.json").read_text())
        declared = {}
        for line in (directory / "SHA256SUMS").read_text().splitlines():
            match = re.fullmatch(r"([0-9a-f]{64})  (products\.csv|covariance\.csv|certificate\.txt|COMMIT\.json)", line)
            if not match or match[2] in declared:
                raise ValueError("INVALID_OR_DUPLICATE_HASH_MANIFEST_ENTRY")
            declared[match[2]] = match[1]
        required = {"products.csv", "certificate.txt", "COMMIT.json"}
        if (directory / "covariance.csv").is_file():
            required.add("covariance.csv")
        if set(declared) != required:
            errors.append("HASH_MANIFEST_FILE_COVERAGE_MISMATCH")
        for name in sorted(required):
            path = directory / name
            if not path.is_file():
                errors.append("MISSING_FILE:" + name)
                continue
            hashes[name] = sha256(path)
            if hashes[name] != declared.get(name):
                errors.append("SHA256_MISMATCH:" + name)
        for name, field in (("products.csv", "products_bytes"), ("covariance.csv", "covariance_bytes")):
            actual = (directory / name).stat().st_size if (directory / name).is_file() else 0
            if manifest.get(field) != actual:
                errors.append("BYTE_COUNT_MISMATCH:" + name)
        if manifest.get("ar") is True and "covariance.csv" not in required:
            errors.append("AR_COVARIANCE_MISSING")
        if gpst(manifest.get("epoch")) is None or not manifest.get("root") or not manifest.get("route"):
            errors.append("INVALID_EPOCH_ROOT_OR_ROUTE")
        if "schema_version" in manifest:
            if manifest["schema_version"] != 2 or not isinstance(manifest.get("systems"), dict):
                errors.append("INVALID_MULTI_SYSTEM_SCHEMA")
            else:
                # Hash exact bytes, retaining original CR/LF distinction.
                certificate = (directory / "certificate.txt").read_bytes()
                actual_fixed = defaultdict(Counter)
                with (directory / "products.csv").open(newline="", encoding="utf-8-sig") as stream:
                    for row in csv.DictReader(stream):
                        if row.get("solution") != "PRODUCT_FIXED":
                            continue
                        system = NAMES.get(row.get("satellite", "")[:1], "UNDECLARED")
                        actual_fixed[system]["rows"] += 1
                        actual_fixed[system]["usable"] += row.get("pppar_usable") == "1"
                for system, item in manifest["systems"].items():
                    letter = LETTERS.get(system)
                    if letter not in SIGNALS or item.get("signals") != list(SIGNALS[letter]):
                        errors.append("SYSTEM_SIGNAL_DECLARATION_MISMATCH:" + system)
                    if item.get("fixed_phase_rows") != actual_fixed[system]["rows"] or \
                            item.get("pppar_usable_phase_rows") != actual_fixed[system]["usable"] or \
                            item.get("any_pppar_usable") != (actual_fixed[system]["usable"] > 0):
                        errors.append("SYSTEM_PRODUCT_COUNT_OR_FLAG_MISMATCH:" + system)
                    if item.get("epoch_certificate_sha256") != hashlib.sha256(certificate).hexdigest():
                        errors.append("EPOCH_CERTIFICATE_HASH_MISMATCH:" + system)
                    begin = ("multi_system_certificate_begin system=" + system + "\n").encode()
                    end = ("multi_system_certificate_end system=" + system + "\n").encode()
                    declared_hash = item.get("system_certificate_sha256")
                    if certificate.count(begin) == 1 and certificate.count(end) == 1:
                        first, last = certificate.find(begin), certificate.find(end)
                        section = certificate[first:last + len(end)]
                        if last < first or declared_hash != hashlib.sha256(section).hexdigest():
                            errors.append("SYSTEM_CERTIFICATE_HASH_MISMATCH:" + system)
                    elif declared_hash is not None or item.get("any_pppar_usable") is True:
                        errors.append("SYSTEM_CERTIFICATE_SECTION_MISSING_OR_AMBIGUOUS:" + system)
    except (OSError, ValueError, TypeError, KeyError) as error:
        errors.append(str(error))
    return {"directory": str(directory), "integrity_verified": not errors,
            "errors": errors, "computed_sha256": hashes, "manifest": manifest,
            "scope": "File identity and metadata only; integer proofs/scientific accuracy not certified"}


def parameter(row, prefix):
    return (row[prefix + "satellite"], row[prefix + "parameter"], row[prefix + "observable"])


def covariance_audit(product_rows, covariance_rows, max_parameters=2000):
    result = {"status": "AVAILABLE", "entry_rows": len(covariance_rows)}
    if not covariance_rows:
        return unavailable("No covariance records for this product epoch/solution")
    pairs, directional, parameters, nonfinite = {}, {}, set(), 0
    duplicate, contradiction, asymmetric = 0, 0, 0.0
    invalid_parameters = set()
    for row in covariance_rows:
        if not all(name in row for name in COV_COLUMNS):
            return {"status": "INVALID", "reason": "COVARIANCE_HEADER_MISSING_COLUMNS"}
        left, right = parameter(row, "row_"), parameter(row, "column_")
        for value in (left, right):
            if not re.fullmatch(r"[GEC]\d{2}", value[0]) or value[1] not in ("CLOCK", "PHASE") or \
                    (value[1] == "CLOCK" and value[2] != "NONE") or \
                    (value[1] == "PHASE" and value[2] not in SIGNALS.get(value[0][:1], ())):
                invalid_parameters.add(value)
        parameters.update((left, right))
        value = finite(row["covariance_m2"])
        if value is None:
            nonfinite += 1
        key = tuple(sorted((left, right)))
        if key in pairs:
            duplicate += 1
            if value is None or pairs[key] is None or abs(value - pairs[key]) > 1e-12:
                contradiction += 1
        else:
            pairs[key] = value
        if (right, left) in directional and value is not None and directional[(right, left)] is not None:
            asymmetric = max(asymmetric, abs(value - directional[(right, left)]))
        directional[(left, right)] = value
    required_all, required_usable = set(), set()
    for row in product_rows:
        wanted = {(row["satellite"], "CLOCK", "NONE"), (row["satellite"], "PHASE", row["observable"])}
        required_all.update(wanted)
        if row.get("ppp_usable") == "1" or row.get("pppar_usable") == "1":
            required_usable.update(wanted)
    count = len(parameters)
    complete = len(pairs) == count * (count + 1) // 2
    by_system = defaultdict(Counter)
    for satellite, kind, code in parameters:
        by_system[satellite[0]][kind] += 1
        if kind == "PHASE":
            by_system[satellite[0]]["PHASE:" + code] += 1
    cross = {}
    systems = sorted({p[0][0] for p in parameters})
    for i, system in enumerate(systems):
        for other in systems[i + 1:]:
            values = [value for (left, right), value in pairs.items()
                      if {left[0][0], right[0][0]} == {system, other} and value is not None]
            cross[system + "_" + other] = {"finite_entries": len(values),
                "exact_nonzero_entries": sum(value != 0.0 for value in values),
                "maximum_absolute_covariance_m2": max((abs(v) for v in values), default=None),
                "status": "NONZERO_BLOCK_PRESENT" if any(v != 0.0 for v in values) else "NO_NONZERO_ENTRY_OBSERVED"}
    result.update(parameter_count=count, expected_unique_triangle_pairs=count * (count + 1) // 2,
        unique_pairs=len(pairs), complete_pair_coverage=complete, nonfinite_entries=nonfinite,
        duplicate_or_reverse_rows=duplicate, contradictory_duplicate_pairs=contradiction,
        maximum_observed_directional_asymmetry_m2=asymmetric,
        symmetry_scope="Upper triangular wire format implies reconstructed symmetry; independently observed only where both directions exist",
        missing_all_product_parameters=[list(v) for v in sorted(required_all - parameters)],
        missing_usable_product_parameters=[list(v) for v in sorted(required_usable - parameters)],
        invalid_parameters=[list(v) for v in sorted(invalid_parameters)],
        parameter_coverage_by_system={s: dict(v) for s, v in by_system.items()}, cross_system_blocks=cross)
    result["psd"] = unavailable("Incomplete/nonfinite/contradictory covariance cannot be completed with fabricated zeros")
    if count > max_parameters:
        result["psd"] = unavailable("Parameter count exceeds explicit current-epoch memory limit")
    elif complete and not nonfinite and not contradiction and not invalid_parameters:
        ordered = sorted(parameters)
        indices = {key: i for i, key in enumerate(ordered)}
        matrix = np.empty((count, count), dtype=float)
        for (left, right), value in pairs.items():
            i, j = indices[left], indices[right]
            matrix[i, j] = matrix[j, i] = value
        spectrum = np.linalg.eigvalsh(matrix)
        tolerance = max(1.0, float(np.max(np.abs(matrix)))) * 1e-10
        result["psd"] = {"status": "AVAILABLE", "minimum_eigenvalue_m2": float(spectrum[0]),
            "tolerance_m2": tolerance, "is_psd_with_tolerance": bool(spectrum[0] >= -tolerance),
            "method": "numpy.linalg.eigvalsh on complete current-epoch full Q, no jitter/clipping"}
    result["internally_consistent"] = bool(complete and not nonfinite and not contradiction and not invalid_parameters
        and not (required_all - parameters) and asymmetric <= 1e-12
        and result["psd"].get("is_psd_with_tolerance") is True)
    return result


class StageAuditor:
    def __init__(self, root, step=30.0, jump_threshold=0.1, max_events=1000, max_parameters=2000):
        self.root = Path(root)
        self.step, self.jump_threshold = step, jump_threshold
        self.max_events, self.max_parameters = max_events, max_parameters
        self.issues, self.event_counts, self.events = Counter(), Counter(), []
        self.last = {}
        self.last_endpoints = {}
        self.expected_systems = set()
        self.totals = defaultdict(Counter)
        self.reason_counts = defaultdict(Counter)
        self.epoch_reports, self.epochs, self.systems = [], set(), set()
        self.max_rows, self.max_cov_rows, self.max_cov_parameters = 0, 0, 0

    def event(self, kind, **details):
        self.event_counts[kind] += 1
        if len(self.events) < self.max_events:
            self.events.append({"kind": kind, **details})

    def process_epoch(self, epoch, products, covariance, authority):
        self.max_rows = max(self.max_rows, len(products))
        self.max_cov_rows = max(self.max_cov_rows, len(covariance))
        self.epochs.add(epoch)
        epoch_result = {"gpst_seconds": epoch, "authority": authority, "solutions": {}}
        current_solutions = {r.get("solution", "UNAVAILABLE") for r in products}
        for solution in sorted(current_solutions):
            rows = [r for r in products if r.get("solution") == solution]
            qrows = [r for r in covariance if r.get("solution") == solution]
            counters, observed = defaultdict(Counter), defaultdict(set)
            seen = set()
            for row in rows:
                if not all(name in row for name in ("satellite", "observable", "solution", *NUMBERS)):
                    self.issues["PRODUCT_REQUIRED_COLUMNS_MISSING"] += 1
                    continue
                satellite, code = row["satellite"], row["observable"]
                system = satellite[:1]
                self.systems.add(system)
                current = counters[system]
                current["phase_rows"] += 1
                observed[system].add((satellite, code))
                key = (solution, satellite, code)
                if key in seen:
                    self.issues["DUPLICATE_PRODUCT_ENDPOINT"] += 1
                seen.add(key)
                numbers = {name: finite(row.get(name)) for name in NUMBERS}
                numeric = all(value is not None for value in numbers.values())
                current["finite_numeric_rows"] += numeric
                if any(numbers.get(name) is not None and numbers[name] < 0 for name in
                       ("clock_sigma_m", "phase_sigma_m", "correction_sigma_m")):
                    self.issues["NEGATIVE_PRODUCT_SIGMA"] += 1
                for name in NUMBERS:
                    current["nonfinite:" + name] += numbers[name] is None
                if not re.fullmatch(r"[GEC]\d{2}", satellite) or code not in SIGNALS.get(system, ()):
                    self.issues["WRONG_SYSTEM_OR_EXACT_SIGNAL_PRODUCT"] += 1
                for flag in FLAGS:
                    if flag not in row:
                        current["unavailable_flag:" + flag] += 1
                    elif row[flag] not in ("0", "1"):
                        self.issues["INVALID_BOOLEAN_FLAG:" + flag] += 1
                    else:
                        current[flag + "_rows"] += row[flag] == "1"
                if row.get("pppar_usable") == "1" and (not numeric or any(row.get(f) != "1" for f in REQUIRED_AR)):
                    self.issues["USABLE_AR_WITH_INVALID_NUMERIC_OR_AUTHORITY_FLAGS"] += 1
                if row.get("ppp_usable") == "1" and (not numeric or row.get("numeric_valid") != "1" or row.get("branch_valid") != "1"):
                    self.issues["USABLE_PPP_WITH_INVALID_NUMERIC_OR_BRANCH_FLAGS"] += 1
                if row.get("pppar_usable") != "1":
                    self.reason_counts[solution + ":" + system][row.get("invalid_reason", "UNAVAILABLE")] += 1
                previous = self.last.get(key)
                continuity = ("public_phase_segment_id", "discontinuity_counter", "datum_version")
                if previous:
                    gap = epoch - previous["epoch"]
                    if gap > self.step + 1e-6:
                        self.event("PRODUCT_ENDPOINT_GAP", epoch=epoch, solution=solution, satellite=satellite,
                            observable=code, absent_grid_epochs=max(0, round(gap / self.step) - 1))
                    if all(row.get(v) not in (None, "", "UNAVAILABLE") and previous.get(v) not in (None, "", "UNAVAILABLE") for v in continuity):
                        changed = [v for v in continuity if row[v] != previous[v]]
                        if changed:
                            self.event("DECLARED_SEGMENT_DATUM_OR_COUNTER_CHANGE", epoch=epoch, solution=solution,
                                satellite=satellite, observable=code, changed=changed, reason=row.get("reset_reason", "UNAVAILABLE"))
                    else:
                        current["continuity_comparison_unavailable_rows"] += 1
                    for name in ("phase_m", "correction_m"):
                        first, second = previous.get(name), numbers[name]
                        if first is not None and second is not None and abs(second - first) > self.jump_threshold:
                            self.event("RAW_ADJACENT_STEP_OVER_THRESHOLD_NOT_ERROR_PROOF", epoch=epoch,
                                solution=solution, satellite=satellite, observable=code, field=name,
                                difference_m=second - first, same_segment_flag=row.get("same_segment_continuous", "UNAVAILABLE"))
                if row.get("requires_user_ambiguity_reset") == "1":
                    self.event("EXPLICIT_USER_RESET_FLAG", epoch=epoch, solution=solution, satellite=satellite, observable=code)
                self.last[key] = {"epoch": epoch, **numbers, **{v: row.get(v) for v in continuity}}
            absent_endpoints = self.last_endpoints.get(solution, set()) - seen
            for unused, satellite, code in sorted(absent_endpoints):
                self.event("PREVIOUSLY_PUBLISHED_ENDPOINT_ABSENT_THIS_EPOCH", epoch=epoch,
                           solution=solution, satellite=satellite, observable=code)
            self.last_endpoints[solution] = seen
            systems_report = {}
            for system, current in counters.items():
                # Count an actually emitted one-band endpoint as a missing
                # partner, not as proof that the other band was observable.
                partners = defaultdict(set)
                for satellite, code in observed[system]:
                    partners[satellite].add(code)
                current["satellites"] = len(partners)
                current["satellites_missing_published_frequency_partner"] = sum(len(v) != 2 for v in partners.values())
                for satellite, codes in partners.items():
                    if len(codes) != 2 and any(r.get("satellite") == satellite and
                            (r.get("pppar_usable") == "1" or r.get("dual_frequency_ar_valid") == "1") for r in rows):
                        self.issues["DUAL_AR_AUTHORITY_WITH_MISSING_PUBLISHED_PARTNER"] += 1
                self.totals[solution + ":" + system].update(current)
                systems_report[system] = dict(current)
            q = covariance_audit(rows, qrows, self.max_parameters)
            self.max_cov_parameters = max(self.max_cov_parameters, q.get("parameter_count", 0))
            if not q.get("internally_consistent", False):
                self.issues["FULL_Q_UNAVAILABLE_OR_INTERNAL_CHECK_FAILED"] += 1
            epoch_result["solutions"][solution] = {"systems": systems_report,
                "systems_without_rows": sorted(self.expected_systems - set(counters)),
                "previously_published_endpoints_absent": [list(v[1:]) for v in sorted(absent_endpoints)],
                "full_covariance": q}
        for solution in sorted(set(self.last_endpoints) - current_solutions):
            absent = self.last_endpoints[solution]
            for unused, satellite, code in sorted(absent):
                self.event("PREVIOUSLY_PUBLISHED_ENDPOINT_ABSENT_THIS_EPOCH", epoch=epoch,
                           solution=solution, satellite=satellite, observable=code)
            if absent:
                epoch_result["solutions"][solution] = {"systems": {}, "no_product_rows": True,
                    "systems_without_rows": sorted(self.expected_systems),
                    "previously_published_endpoints_absent": [list(v[1:]) for v in sorted(absent)],
                    "full_covariance": unavailable("Previously emitted solution has no current product rows")}
            self.last_endpoints[solution] = set()
        self.epoch_reports.append(epoch_result)

    def paired_files(self, products, covariance, authority, manifest_epoch=None):
        qstream = iter(csv_epochs(covariance)) if covariance and covariance.is_file() else iter(())
        next_q = next(qstream, None)
        for epoch, rows in csv_epochs(products):
            if manifest_epoch is not None and abs(epoch - manifest_epoch) > 1e-6:
                self.issues["PRODUCT_EPOCH_DIFFERS_FROM_COMMIT"] += 1
            while next_q and next_q[0] < epoch:
                self.issues["COVARIANCE_EPOCH_WITHOUT_PRODUCTS"] += 1
                next_q = next(qstream, None)
            qrows = []
            consumed_q = False
            if next_q and next_q[0] == epoch:
                qrows = next_q[1]
                consumed_q = True
            self.process_epoch(epoch, rows, qrows, authority)
            # Release the current Q before advancing. Advancing before the
            # audit would retain both current and following epoch matrices.
            qrows = []
            if consumed_q:
                next_q = None
                next_q = next(qstream, None)
        while next_q:
            self.issues["COVARIANCE_EPOCH_WITHOUT_PRODUCTS"] += 1
            next_q = next(qstream, None)

    def run(self):
        try:
            plan = json.loads((self.root / "stage_plan.json").read_text())
        except (OSError, ValueError) as error:
            plan = {}
            self.issues["STAGE_PLAN_MISSING_OR_INVALID"] += 1
        expected_systems = {LETTERS[v] for v in plan.get("specification", {}).get("systems", []) if v in LETTERS}
        self.expected_systems = expected_systems
        logs = audit_logs(self.root)
        base = self.root / "outputs/zhang_internal_products.csv.epochs"
        bundles, epoch_directories = [], Counter()
        pending_count = 0
        if base.is_dir():
            for directory in sorted(base.iterdir()):
                if not directory.is_dir():
                    continue
                if directory.name.endswith(".pending"):
                    pending_count += 1
                    continue
                if not (directory / "COMMIT.json").is_file():
                    self.issues["UNCOMMITTED_NONPENDING_DIRECTORY"] += 1
                    continue
                item = verify_bundle(directory)
                bundles.append(item)
                if not item["integrity_verified"]:
                    self.issues["ATOMIC_BUNDLE_INTEGRITY_FAILED"] += 1
                    continue
                manifest = item["manifest"]
                epoch_directories[manifest.get("epoch")] += 1
                if len(expected_systems) > 1 and manifest.get("schema_version") != 2:
                    self.issues["MULTI_SYSTEM_COMMIT_SCHEMA_MISSING"] += 1
                if "systems" in manifest:
                    actual = {LETTERS.get(v, v) for v in manifest["systems"]}
                    if expected_systems and actual != expected_systems:
                        self.issues["COMMIT_SYSTEMS_DIFFER_FROM_STAGE"] += 1
                try:
                    self.paired_files(directory / "products.csv", directory / "covariance.csv",
                        "HASH_VERIFIED_ATOMIC_BUNDLE", gpst(manifest["epoch"]))
                except (OSError, ValueError, KeyError, TypeError, np.linalg.LinAlgError) as error:
                    self.issues["ATOMIC_PRODUCT_OR_Q_PARSE_FAILED"] += 1
                    item["content_error"] = str(error)
        aggregate = self.root / "outputs/zhang_internal_products.csv"
        if not bundles and aggregate.is_file():
            self.issues["ONLY_UNCOMMITTED_AGGREGATE_AVAILABLE"] += 1
            try:
                self.paired_files(aggregate, self.root / "outputs/zhang_internal_product_covariance.csv",
                                  "UNCOMMITTED_AGGREGATE_NOT_PUBLICATION_AUTHORITY")
            except (OSError, ValueError, KeyError, TypeError, np.linalg.LinAlgError) as error:
                self.issues["AGGREGATE_PARSE_FAILED"] += 1
        if any(count > 1 for count in epoch_directories.values()):
            self.issues["DUPLICATE_COMMITTED_EPOCH"] += sum(count - 1 for count in epoch_directories.values() if count > 1)
        if expected_systems and self.systems and self.systems != expected_systems:
            self.issues["ACTUAL_PRODUCT_SYSTEMS_DIFFER_FROM_STAGE_NOT_NEW_MULTI_RESULT"] += 1
        expected = plan.get("expected_epoch_count")
        processed = logs["processed_epochs"].get("count")
        if expected is not None and processed is not None and processed != expected:
            self.issues["ACTUAL_PROCESSED_EPOCH_COUNT_DIFFERS_FROM_PLAN"] += 1
        state = plan.get("state", "UNAVAILABLE")
        evidence = bool(bundles or pending_count or aggregate.is_file() or processed or logs["files_scanned"])
        if plan.get("exit_code") not in (None, 0) or state == "EXECUTION_FAILED":
            execution = "EXECUTION_FAILED"
        elif plan.get("exit_code") == 0 and state.startswith("COMPLETED"):
            if not self.epochs:
                self.issues["COMPLETED_WITHOUT_PRODUCT_EPOCHS"] += 1
            execution = "COMPLETED_WITH_EVIDENCE_FAILURES" if self.issues else "COMPLETED_REQUIRES_SCIENTIFIC_AUDIT"
        elif not evidence and state == "PREPARED_NOT_RUN":
            execution = "PREPARED_NOT_RUN"
        elif not evidence and state.startswith("BLOCKED"):
            execution = state
        elif evidence:
            execution = "RUN_EVIDENCE_PRESENT_COMPLETION_UNCONFIRMED"
        else:
            execution = "NO_VERIFIABLE_RUN_EVIDENCE"
        absent = expected_systems - self.systems
        actual_epochs = sorted(self.epochs)
        gaps = [int(round((b - a) / self.step)) - 1 for a, b in zip(actual_epochs, actual_epochs[1:]) if b - a > self.step + 1e-6]
        return {"stage": plan.get("stage", "UNAVAILABLE"), "run_root": str(self.root),
            "declared_state": state, "execution_classification": execution,
            "scientific_acceptance": False, "scientific_acceptance_status": "NOT_ASSESSED_REQUIRES_INDEPENDENT_VALIDATION",
            "scope": "Diagnostics of these exact output files only; no old GPS run is relabelled multi-GNSS",
            "expected_systems": sorted(expected_systems), "actual_product_systems": sorted(self.systems),
            "systems_without_any_product_rows": sorted(absent), "plan_expected_epoch_count": expected,
            "logs": logs, "atomic": {"bundle_directories": len(bundles),
                "hash_verified_bundles": sum(v["integrity_verified"] for v in bundles),
                "pending_directories_not_published": pending_count, "bundles": bundles},
            "products": {"epoch_count": len(self.epochs), "first_gpst_seconds": actual_epochs[0] if actual_epochs else None,
                "last_gpst_seconds": actual_epochs[-1] if actual_epochs else None,
                "missing_grid_epochs_between_observed_products": sum(gaps),
                "totals_by_solution_system": {k: dict(v) for k, v in sorted(self.totals.items())},
                "unusable_AR_reason_counts": {k: dict(v) for k, v in sorted(self.reason_counts.items())},
                "visibility_based_expected_satellites": unavailable("No current retained-observation visibility table is exported; absent product is not proof of an observable satellite"),
                "epochs": self.epoch_reports},
            "continuity": {"event_counts": dict(self.event_counts), "event_samples": self.events,
                "event_sample_limit": self.max_events, "omitted_event_samples": sum(self.event_counts.values()) - len(self.events),
                "raw_step_threshold_m": self.jump_threshold,
                "warning": "Raw adjacent value changes are not by themselves physical slips, datum errors or certified continuity"},
            "evidence_issues": dict(self.issues),
            "streaming_memory": {"policy": "Current epoch rows/Q only; endpoint last values and compact epoch reports retained",
                "maximum_product_rows_in_memory": self.max_rows, "maximum_covariance_rows_in_memory": self.max_cov_rows,
                "maximum_covariance_parameters_in_memory": self.max_cov_parameters, "matrix_parameter_limit": self.max_parameters}}


def render_text(report):
    lines = ["MULTI-GNSS REAL OUTPUT EVIDENCE AUDIT", "Scientific acceptance: NOT ASSESSED (always false)", ""]
    for stage in report["stages"]:
        lines.extend([f"{stage['stage']}: {stage['execution_classification']}", f"  root: {stage['run_root']}",
            f"  systems expected/actual: {stage['expected_systems']} / {stage['actual_product_systems']}",
            f"  processed log epochs: {stage['logs']['processed_epochs']}",
            f"  product epochs: {stage['products']['epoch_count']}; verified atomic bundles: {stage['atomic']['hash_verified_bundles']}",
            f"  SPP SINGLE evidence: {stage['logs']['spp_single']}",
            f"  issues: {stage['evidence_issues']}", f"  continuity events: {stage['continuity']['event_counts']}"])
        for key, value in stage["products"]["totals_by_solution_system"].items():
            lines.append(f"  {key}: rows={value.get('phase_rows', 0)} finite={value.get('finite_numeric_rows', 0)} "
                         f"PPP={value.get('ppp_usable_rows', 0)} AR={value.get('pppar_usable_rows', 0)}")
        lines.append("")
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_roots", nargs="+", type=Path)
    parser.add_argument("--json-output", required=True, type=Path)
    parser.add_argument("--text-output", type=Path)
    parser.add_argument("--epoch-interval", type=float, default=30)
    parser.add_argument("--raw-step-threshold-m", type=float, default=.1)
    parser.add_argument("--max-events", type=int, default=1000)
    parser.add_argument("--max-covariance-parameters", type=int, default=2000)
    args = parser.parse_args()
    if args.epoch_interval <= 0 or args.max_events < 0 or args.max_covariance_parameters < 1:
        parser.error("Invalid positive interval/memory limit or negative event limit")
    report = {"schema_version": 1, "generated_utc": datetime.now(timezone.utc).isoformat(),
        "scientific_acceptance": False,
        "stages": [StageAuditor(root, args.epoch_interval, args.raw_step_threshold_m, args.max_events,
                                args.max_covariance_parameters).run() for root in args.run_roots]}
    args.json_output.parent.mkdir(parents=True, exist_ok=True)
    args.json_output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    text = render_text(report)
    if args.text_output:
        args.text_output.parent.mkdir(parents=True, exist_ok=True)
        args.text_output.write_text(text, encoding="utf-8")
    print(text)


if __name__ == "__main__":
    main()
