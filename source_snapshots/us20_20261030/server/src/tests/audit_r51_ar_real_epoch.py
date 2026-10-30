"""Validate committed probe bundles and extract real (nested) phase timings.
Usage: python3 audit_r51_ar_real_epoch.py PROBE_DIRECTORY EVIDENCE_DIRECTORY CHECKPOINT_VERIFIER
"""
from collections import Counter, defaultdict
from pathlib import Path
import csv
import hashlib
import json
import math
import re
import shutil
import subprocess
import sys

probe_dir, evidence = map(Path, sys.argv[1:3])
probe = json.loads((probe_dir / "probe.json").read_text())
assert probe["returncode"] == probe["expected_returncode"] == 130 and \
    not probe["solver_still_running"] and not probe["timed_out"] and probe["checkpoint_bundles"]
run = Path(probe["run_root"])
assert run.parent == Path("/home/rx/GINAN/local_timing")
evidence.mkdir(parents=True, exist_ok=True)
trace = run / "outputs/Network-case-202419900.TRACE"
phases = defaultdict(lambda: dict(calls=0, inclusive_seconds=0., exclusive_seconds=0.))
floats = []
records = defaultdict(list)
tags = ("ZHANG_AR_THREADS", "ZHANG_FLOAT_THREADS", "R49_KF_ROOT_WORKSPACE",
        "ZHANG_R51_ROOT_NIS_POLICY", "ZHANG_R51_EPOCH_ONLY_AR", "ZHANG_AR_SUMMARY",
        "ZHANG_R51_BLOCK_", "ZHANG_R51_PAR_SEARCH", "ZHANG_R51_RATIO",
        "ZHANG_RATIO_GATE",
        "ZHANG_R51_FLOAT_AUTHORITY", "E29_RESTORE", "ZHANG_E29_RESTORE")
ar_active=False
with trace.open() as stream:
    for line in stream:
        if line.startswith("ZHANG_AR_THREADS "):
            ar_active=True
        if line.startswith("ZHANG_AR_PHASE_TIMER "):
            if not ar_active:
                continue
            fields = dict(re.findall(r"(\w+)=(\S+)", line))
            p = phases[fields["phase"]]
            p["calls"] += 1
            p["inclusive_seconds"] += float(fields["inclusive_ms"]) / 1000
            p["exclusive_seconds"] += float(fields["exclusive_ms"]) / 1000
            if fields["phase"]=="AR_TOTAL":
                ar_active=False
        elif line.startswith("ZHANG_FLOAT_PHASE_TIMER phase=FLOAT_TOTAL "):
            floats.append(line.rstrip())
        else:
            for tag in tags:
                if line.startswith(tag):
                    records[tag].append(line.rstrip())
                    break
assert phases["AR_TOTAL"]["calls"] == 1 and len(floats)==1, "Probe did not stop after one real epoch"
assert any("requested=12" in s for s in records["ZHANG_AR_THREADS"])
assert any("policy=hotspot_all_cpu_v2" in s and "eigen_tridiagonal_blas=1" in s
           for s in records["ZHANG_AR_THREADS"]), "Adaptive policy was not used"
assert records["ZHANG_R51_ROOT_NIS_POLICY"] and all("evaluation=SKIPPED_RATIO_ONLY" in s
    for s in records["ZHANG_R51_ROOT_NIS_POLICY"]), "Joint/root NIS policy changed"
assert any("prior_integer_reuse=0" in s for s in records["ZHANG_R51_EPOCH_ONLY_AR"])
assert records["ZHANG_RATIO_GATE"] and all(float(re.search(r"\bthreshold=(\S+)", s)[1]) == 2.
    for s in records["ZHANG_RATIO_GATE"]), "Effective ratio is not uniformly 2"
assert any(int(re.search(r"\bprotected_float_hits=(\d+)", s)[1]) > 0
           for s in records["R49_KF_ROOT_WORKSPACE"]), \
    "Protected FLOAT root was not reused in production"

bundles = []
for commit_path in map(Path, probe["committed_epochs"]):
    commit = json.loads(commit_path.read_text())
    directory = commit_path.parent
    for line in (directory / "SHA256SUMS").read_text().splitlines():
        expected, name = line.split(maxsplit=1)
        assert Path(name).name == name and name in (
            "products.csv", "covariance.csv", "certificate.txt", "COMMIT.json")
        value = hashlib.sha256((directory / name).read_bytes()).hexdigest()
        assert value == expected, f"Bundle hash mismatch: {directory / name}"
    assert (directory / "products.csv").stat().st_size == commit["products_bytes"]
    assert (directory / "covariance.csv").stat().st_size == commit["covariance_bytes"]
    with (directory / "products.csv").open() as stream:
        rows = list(csv.DictReader(stream))
    fixed = [r for r in rows if r["ar_valid"] == "1" and r["pppar_usable"] == "1"]
    for row in fixed:
        for field in ("clock_m", "clock_sigma_m", "phase_m", "phase_sigma_m",
                      "clock_phase_covariance_m2", "correction_m", "correction_sigma_m"):
            assert math.isfinite(float(row[field])), f"Nonfinite published AR {field}"
        clock_sigma, phase_sigma, correction_sigma = (
            float(row[f + "_sigma_m"]) for f in ("clock", "phase", "correction"))
        cross = float(row["clock_phase_covariance_m2"])
        assert min(clock_sigma, phase_sigma, correction_sigma) >= 0
        assert abs(cross) <= clock_sigma * phase_sigma + 1e-7
        expected = clock_sigma ** 2 + phase_sigma ** 2 - 2 * cross
        assert abs(correction_sigma ** 2 - expected) < 1e-7 * max(1., abs(expected)), \
            "Published clock/phase/correction covariance identity failed"
    satellites = defaultdict(set)
    for row in fixed:
        satellites[row["satellite"]].add(row["observable"])
    summary = dict(commit=commit, path=str(directory), rows=len(rows),
                   solution_counts=dict(Counter(r["solution"] for r in rows)),
                   product_state_counts=dict(Counter(r["product_state"] for r in rows)),
                   ar_usable_rows=len(fixed),
                   dual_frequency_ar_satellites=sorted(s for s, f in satellites.items() if len(f) >= 2))
    bundles.append(summary)
    target = evidence / "bundles" / directory.name
    target.mkdir(parents=True, exist_ok=True)
    for name in ("COMMIT.json", "SHA256SUMS", "products.csv", "covariance.csv", "certificate.txt"):
        shutil.copy2(directory / name, target / name)

assert any(b["commit"]["ar"] and b["dual_frequency_ar_satellites"] for b in bundles), \
    "No dual-frequency AR-usable published products"
runner = (run / "runner.log").read_text()
assert "E29_CAPTURE_PROVENANCE_CHANGED_SINCE_STARTUP" not in runner
assert "E29_CAPTURE_SCHEDULE_INTERRUPTED" in runner, "Incomplete hour was not explicitly marked interrupted"
checkpoints = []
assert len(sys.argv)>3, "Provide the production checkpoint verifier executable"
for stored in probe["checkpoint_bundles"]:
    path=Path(stored)
    assert path.resolve().is_relative_to((run/"checkpoints").resolve())
    metadata=path.parent/"checkpoint_manifest.json"
    m=json.loads(metadata.read_text())
    assert m["binary_sha256"]==probe["binary_sha256"]
    for field in ("config","input_manifest"):
        assert hashlib.sha256(m[field+"_text"].encode()).hexdigest()==m[field+"_sha256"]
    verified=subprocess.run([sys.argv[3],str(path),probe["binary_sha256"]],
        text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=90)
    (evidence/(path.parent.name+"_checkpoint_verify.log")).write_text(verified.stdout)
    assert verified.returncode==0, verified.stdout
    value=hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda:stream.read(4*1024*1024),b""):
            value.update(block)
    checkpoints.append(dict(path=str(path),sha256=value.hexdigest(),bytes=path.stat().st_size,
        epoch=m["epoch"],state_dimension=m["state_dimension"],production_decode_passed=True))
    shutil.copy2(metadata,evidence/(path.parent.name+"_manifest.json"))
rss = re.search(r"Maximum resident set size \(kbytes\):\s*(\d+)", runner)
audit = dict(scope="Bounded real-epoch acceptance, not a paired full-hour speed test",
             probe=probe, committed_epoch_times=sorted(set(b["commit"]["epoch"] for b in bundles)),
             bundles=bundles, phases=dict(phases), float_total_records=floats,
             verified_checkpoints=checkpoints, intentional_interruption=True, hour_complete=False,
             trace_records=dict(records), peak_rss_mib=int(rss[1]) / 1024 if rss else None,
             warning="Inclusive timings are nested. One epoch cannot establish continuity or precision.")
for name, source in (("probe.json", probe_dir / "probe.json"),
                     ("launch.log", probe_dir / "launch.log"),
                     ("runner.log", run / "runner.log"), (trace.name, trace)):
    shutil.copy2(source, evidence / name)
(evidence / "real_epoch_audit.json").write_text(json.dumps(audit, indent=2))
print(json.dumps({"epochs": audit["committed_epoch_times"], "peak_rss_mib": audit["peak_rss_mib"],
                  "ar_total": phases["AR_TOTAL"], "bundles": bundles}, indent=2))
