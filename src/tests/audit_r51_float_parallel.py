"""Paired, sequential production-kernel benchmarks and full-posterior audit.

No GNSS jobs or AR search are started. The dense inputs are synthetic; the
test does not establish complete real-network epoch timing/product quality.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import statistics
import subprocess

import numpy as np


def digest(path):
    value = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            value.update(block)
    return value.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline', type=Path, required=True)
    parser.add_argument('--optimized', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--repeats', type=int, default=5)
    parser.add_argument('--states', type=int, default=4200)
    parser.add_argument('--measurements', type=int, default=5600)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    modes = {'baseline_4': (args.baseline, 4), 'optimized_4': (args.optimized, 4),
             'optimized_6': (args.optimized, 6), 'optimized_12': (args.optimized, 12)}
    observations = {mode: [] for mode in modes}
    environment = {**os.environ, 'OMP_NUM_THREADS': '4', 'OPENBLAS_NUM_THREADS': '4',
                   'MKL_NUM_THREADS': '1'}
    for repeat in range(args.repeats):
        # Rotate order to avoid assigning warm/cold CPU state to only one mode.
        order = list(modes)
        order = order[repeat % len(order):] + order[:repeat % len(order)]
        for mode in order:
            binary, threads = modes[mode]
            command = [str(binary), str(threads), str(args.states), str(args.measurements)]
            if repeat == 0:
                command.append(str(args.output / (mode + '.posterior')))
            run = subprocess.run(command, text=True, stdout=subprocess.PIPE,
                                 stderr=subprocess.STDOUT, env=environment)
            (args.output / f'{mode}_{repeat + 1}.log').write_text(run.stdout)
            if run.returncode:
                raise RuntimeError(f'{mode} failed; see preserved log')
            match = re.search(r'^BENCH (.*)$', run.stdout, re.M)
            if not match:
                raise RuntimeError(f'{mode} produced no benchmark record')
            record = {key: float(value) for key, value in
                      re.findall(r'(\w+)=([0-9.eE+-]+)', match.group(1))}
            record['phases_ms'] = {key: float(value) for key, value in re.findall(
                r'ZHANG_FLOAT_KERNEL phase=(\w+).*?elapsed_ms=([0-9.eE+-]+)', run.stdout)}
            observations[mode].append(record)
            print(f'{mode} repeat={repeat + 1} wall={record["wall_s"]:.6f}s', flush=True)
    baseline = np.fromfile(args.output / 'baseline_4.posterior', dtype=np.float64)
    if len(baseline) != args.states + args.states ** 2 or not np.isfinite(baseline).all():
        raise RuntimeError('invalid baseline posterior')
    comparisons = {}
    for mode in modes:
        values = np.fromfile(args.output / (mode + '.posterior'), dtype=np.float64)
        if values.shape != baseline.shape or not np.isfinite(values).all():
            raise RuntimeError(f'invalid posterior: {mode}')
        error = np.abs(values - baseline)
        comparisons[mode] = {'mean_max_abs': float(error[:args.states].max()),
                             'full_covariance_max_abs': float(error[args.states:].max()),
                             'tolerance': 2e-11}
        if error.max() > 2e-11:
            raise RuntimeError(f'full posterior mismatch: {mode}')
    summary = {}
    for mode, rows in observations.items():
        walls = [r['wall_s'] for r in rows]
        summary[mode] = {'median_seconds': statistics.median(walls),
                         'mean_seconds': statistics.mean(walls),
                         'min_seconds': min(walls), 'max_seconds': max(walls),
                         'max_peak_rss_mib': max(r['peak_rss_kib'] for r in rows) / 1024,
                         'mean_average_cores': statistics.mean(r['average_cores'] for r in rows),
                         'actual_omp_teams': sorted(set(r['omp_team'] for r in rows)),
                         'actual_blas_threads': sorted(set(r['blas'] for r in rows))}
    result = {'status': 'PASS', 'scope': 'synthetic inputs, real production kFilter, not a full GNSS epoch',
              'states': args.states, 'measurements': args.measurements, 'repeats': args.repeats,
              'binaries': {str(p): digest(p) for p in (args.baseline, args.optimized)},
              'summary': summary, 'full_posterior_comparison': comparisons, 'observations': observations}
    (args.output / 'results.json').write_text(json.dumps(result, indent=2))
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == '__main__':
    main()
