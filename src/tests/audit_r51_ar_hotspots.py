"""Sequential paired production-kernel audit; not a full-epoch benchmark.

Usage: python3 audit_r51_ar_hotspots.py BASELINE OPTIMIZED OUTPUT [REPEATS]
Both inputs are directories produced by run_r51_ar_hotspots.py.
"""
from pathlib import Path
import hashlib
import json
import os
import re
import statistics
import subprocess
import sys

# Keep large comparison arrays out of the timing parent: fork/exec children
# can inherit a pre-exec RSS high-water mark from a memory-heavy parent.
if len(sys.argv) > 1 and sys.argv[1] == "--compare-root":
    import numpy as np
    left = np.fromfile(sys.argv[2], dtype=np.float64)
    right = np.fromfile(sys.argv[3], dtype=np.float64)
    assert left.size == right.size == 1800 * 1800
    assert np.isfinite(left).all() and np.isfinite(right).all()
    delta = np.abs(left - right)
    maximum = float(delta.max())
    assert maximum < 2e-11 * max(1., float(np.abs(left).max()))
    print(json.dumps(dict(entries=left.size, maximum_absolute_difference=maximum,
                          mean_absolute_difference=float(delta.mean()))))
    sys.exit(0)

baseline, optimized, output = map(Path, sys.argv[1:4])
repeats = int(sys.argv[4]) if len(sys.argv) > 4 else 5
output.mkdir(parents=True, exist_ok=False)
env = dict(os.environ, OMP_NUM_THREADS="4", OPENBLAS_NUM_THREADS="4",
           MKL_NUM_THREADS="1", ZHANG_R51_AR_THREADS="all")
records = []
checks = []
sizes = {"root": ["1800"], "compat": ["768"], "closure": ["32", "20"]}
labels = {"root": "ROOT", "compat": "COMPATIBILITY", "closure": "CLOSURE_STRUCTURAL"}

for repeat in range(repeats):
    for mode in sizes:
        # Alternate which version runs first to reduce thermal/order bias.
        versions = [("before", baseline), ("after", optimized)]
        if repeat % 2:
            versions.reverse()
        paths = {}
        for version, directory in versions:
            prefix = output / f"{repeat + 1}_{mode}_{version}"
            result = subprocess.run(
                [str(directory / "ar_hotspot_benchmark"), mode, str(prefix), *sizes[mode]],
                env=env, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                timeout=180)
            prefix.with_suffix(".log").write_text(result.stdout)
            if result.returncode:
                raise RuntimeError(f"{prefix}: exit={result.returncode}\n{result.stdout}")
            seconds = float(re.search(labels[mode] + r" seconds=([0-9.eE+-]+)", result.stdout)[1])
            rss = int(re.search(r"peak_rss_kib=(\d+)", result.stdout)[1])
            record = dict(repeat=repeat + 1, mode=mode, version=version,
                          seconds=seconds, peak_rss_kib=rss, log=prefix.name + ".log")
            records.append(record)
            paths[version] = prefix
            print(json.dumps(record), flush=True)
        if mode == "root":
            comparison = json.loads(subprocess.check_output([
                sys.executable, str(Path(__file__).resolve()), "--compare-root",
                str(paths["before"]) + ".root", str(paths["after"]) + ".root"], text=True))
            checks.append(dict(repeat=repeat + 1, mode=mode, **comparison))
        elif mode == "closure":
            left = Path(str(paths["before"]) + ".basis").read_bytes()
            right = Path(str(paths["after"]) + ".basis").read_bytes()
            assert left == right, "Product basis, witness, search ordering or offsets changed"
            checks.append(dict(repeat=repeat + 1, mode=mode, byte_identical=True,
                               bytes=len(left), sha256=hashlib.sha256(left).hexdigest()))

summary = {}
for mode in sizes:
    summary[mode] = {}
    for version in ("before", "after"):
        selected = [r for r in records if r["mode"] == mode and r["version"] == version]
        summary[mode][version] = dict(
            median_seconds=statistics.median(r["seconds"] for r in selected),
            minimum_seconds=min(r["seconds"] for r in selected),
            maximum_seconds=max(r["seconds"] for r in selected),
            maximum_peak_rss_mib=max(r["peak_rss_kib"] for r in selected) / 1024)
    old, new = (summary[mode][v]["median_seconds"] for v in ("before", "after"))
    summary[mode]["speedup"] = old / new
    summary[mode]["reduction_percent"] = (1 - new / old) * 100

report = dict(scope="Controlled production kernels only; not full product closure or AR epoch",
              repeats=repeats, baseline_omp_blas=[4, 4], optimized_omp_blas=[12, 12],
              optimized_eigensolver_blas_cap=6, sizes=sizes, records=records,
              comparisons=checks, summary=summary)
(output / "results.json").write_text(json.dumps(report, indent=2))
print(json.dumps(summary, indent=2), flush=True)
