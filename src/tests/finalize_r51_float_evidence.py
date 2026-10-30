"""Freeze the exact changed sources and generate an audit manifest/report."""
from datetime import datetime, timezone
import difflib
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys

root = Path(__file__).resolve().parents[2]
output = Path(sys.argv[1]).resolve()
results = json.loads((output / 'paired_final/results.json').read_text())
gates = {
    'compile_pea': '[100%] Built target pea' in (output / 'build_final3.log').read_text(),
    'ctest_35': '100% tests passed, 0 tests failed out of 35' in (output / 'ctest_final.log').read_text(),
    'production_float_qc': 'PASS production FLOAT QC' in (output / 'float_qc_final.log').read_text(),
    'production_graph': 'PASS production graph runtime' in (output / 'graph_final.log').read_text(),
    'production_slip_qc': 'production_adapter=PASS quarantines=0' in (output / 'slip_qc_final.log').read_text(),
    'paired_full_posterior': results['status'] == 'PASS',
    'nis_policy_unchanged': hashlib.sha256((root / 'src/cpp/common/zhangRatioOnly.hpp').read_bytes()).hexdigest()
        == '11360b0f0c15751b3af6f40d521f1999671dc45bd5f66910a414f384f4c61ee8'}
if not all(gates.values()):
    raise RuntimeError(f'failed evidence gate: {gates}')
files = ['src/cpp/common/algebra.cpp', 'src/cpp/common/zhangFloatThreading.hpp',
         'src/cpp/common/zhangRatioOnly.hpp',
         'src/cpp/pea/main.cpp', 'src/cpp/pea/ppppp.cpp', 'src/tests/CMakeLists.txt',
         'src/tests/unit/r51_float_parallel_benchmark.cpp', 'src/tests/run_r51_float_parallel.py',
         'src/tests/audit_r51_float_parallel.py', 'src/tests/finalize_r51_float_evidence.py',
         'experiments/float_all_cpu_20261019/README.txt']
sources = {}
for rel in files:
    destination = output / 'tested_sources' / rel
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(root / rel, destination)
    sources[rel] = hashlib.sha256(destination.read_bytes()).hexdigest()
for rel in ('src/cpp/common/algebra.cpp', 'src/cpp/pea/main.cpp', 'src/cpp/pea/ppppp.cpp'):
    before = output / 'before' / Path(rel).name
    patch = ''.join(difflib.unified_diff(before.read_text().splitlines(keepends=True),
        (root / rel).read_text().splitlines(keepends=True), fromfile='before/' + rel,
        tofile='after/' + rel))
    (output / (Path(rel).name + '.patch')).write_text(patch)
summary = results['summary']
base = summary['baseline_4']
current = summary['optimized_12']
speed = 1 - current['median_seconds'] / base['median_seconds']
memory = 1 - current['max_peak_rss_mib'] / base['max_peak_rss_mib']
manifest = {'status': results['status'], 'created_utc': datetime.now(timezone.utc).isoformat(),
            'test_gates': gates,
            'sources_sha256': sources, 'pea_sha256': hashlib.sha256((root / 'bin/pea').read_bytes()).hexdigest(),
            'benchmark': 'paired_final/results.json', 'thread_policy': 'phase_scoped_all_cpu_v2',
            'scope_limit': results['scope'], 'nis_policy_changed': False, 'ar_started': False,
            'median_kernel_time_reduction_fraction': speed,
            'max_rss_reduction_fraction': memory,
            'hardware': subprocess.check_output(['lscpu'], text=True)}
(output / 'MANIFEST.json').write_text(json.dumps(manifest, indent=2))
lines = ['FLOAT optimization verification (2026-10-19)', '',
         'Scope: deterministic synthetic input; actual production kFilter. NOT a full GNSS epoch.',
         'Host: Ubuntu-22.04-G, Ryzen 9600X, 6 physical cores / 12 logical CPUs.',
         'WSL: processors=12; original 20 GB memory / 16 GB swap retained.',
         'Default: 12 station/GEMM threads, factorisation capped at 6 physical cores.',
         'Inherited AR threading restored at FLOAT exit; ratio/NIS algorithms unchanged.', '',
         'Five sequential repetitions per mode (rotated run order):']
for name, row in summary.items():
    lines.append(f'{name}: median={row["median_seconds"]:.6f}s; '
                 f'min/max={row["min_seconds"]:.6f}/{row["max_seconds"]:.6f}s; '
                 f'peak RSS={row["max_peak_rss_mib"]:.2f} MiB')
lines += ['', f'Median kernel time reduction: {speed:.2%}', f'Peak RSS reduction: {memory:.2%}',
          'Complete mean/covariance differences:']
for mode, difference in results['full_posterior_comparison'].items():
    lines.append(f'{mode}: mean max abs={difference["mean_max_abs"]:.4e}, '
                 f'covariance max abs={difference["full_covariance_max_abs"]:.4e}')
lines += ['', 'See CTest and production-linked QC logs for exact test scope/results.',
          'The initial all-suite CTest attempt was deliberately interrupted to remove CPU contention.',
          'Long slip Monte Carlo tests were not rerun; slip detection code is unchanged.',
          'Earlier v1 timing logs include a deliberately untuned 12-thread factorisation; '
          'only paired_final is the final tuned comparison.',
          'No AR job launched. Old frozen binaries, observations, checkpoints and products preserved.',
          'Existing WSL/Windows-worktree Git version-generation warnings remain in build logs; '
          'the binary and tested source content are identified by SHA256 in MANIFEST.json.']
(output / 'FINAL_REPORT.txt').write_text('\n'.join(lines) + '\n')
print('\n'.join(lines))
