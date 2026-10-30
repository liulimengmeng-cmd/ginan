#!/usr/bin/env python3
"""Isolated deterministic LP64 OpenBLAS benchmark; never starts PEA.

No numpy, packages, global thread settings, or GNSS files are required. The
output and its live sample sidecar must be new files. Resource guards apply to
this one process only. Benchmark timing excludes memcpy and residual checking.
"""
import argparse
from array import array
import ctypes as c
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import resource
import signal
import statistics
import threading
import time

GIB = 1024 ** 3
THREADS = (1, 4, 8, 14, 28)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def proc_ticks(path):
    text = path.read_text()
    fields = text[text.rfind(')') + 2:].split()
    return int(fields[11]) + int(fields[12])


def cpu_snapshot():
    tasks = {}
    for task in Path('/proc/self/task').iterdir():
        try:
            tasks[task.name] = proc_ticks(task / 'stat')
        except FileNotFoundError:
            pass
    return proc_ticks(Path('/proc/self/stat')), tasks


def cgroup_info():
    suffix = next(line.split('::', 1)[1] for line in Path('/proc/self/cgroup').read_text().splitlines()
                  if line.startswith('0::'))
    path = Path('/sys/fs/cgroup') / suffix.lstrip('/')
    rows = []
    while True:
        row = {'path': str(path)}
        for name in ('cpu.max', 'cpu.weight', 'cpuset.cpus.effective', 'cpu.stat'):
            file = path / name
            row[name] = file.read_text().strip() if file.exists() else None
        rows.append(row)
        if path == Path('/sys/fs/cgroup'):
            return rows
        path = path.parent


class Blas:
    def __init__(self, path, n):
        self.lib = c.CDLL(str(path))
        self.lib.openblas_get_config.restype = c.c_char_p
        self.lib.openblas_get_corename.restype = c.c_char_p
        self.config = self.lib.openblas_get_config().decode()
        if 'USE64BITINT' in self.config or 'INTERFACE64' in self.config:
            raise RuntimeError('LP64 benchmark cannot call an ILP64 library')
        for name in ('dgemm_', 'dpotrf_', 'dpotrs_', 'dpotri_', 'dnrm2_', 'daxpy_', 'dcopy_'):
            getattr(self.lib, name)
        self.lib.dnrm2_.restype = c.c_double
        self.n, self.total, self.one = c.c_int(n), c.c_int(n * n), c.c_int(1)
        self.upper, self.lower = c.c_char(b'U'), c.c_char(b'L')
        self.normal, self.transposed = c.c_char(b'N'), c.c_char(b'T')
        self.alpha, self.zero, self.minus = c.c_double(1), c.c_double(0), c.c_double(-1)
        self.bytes = n * n * 8

    def gemm(self, a, b, out, transpose_b=False, scale=1.):
        alpha = c.c_double(scale)
        self.lib.dgemm_(c.byref(self.normal), c.byref(self.transposed if transpose_b else self.normal),
            c.byref(self.n), c.byref(self.n), c.byref(self.n), c.byref(alpha), a, c.byref(self.n),
            b, c.byref(self.n), c.byref(self.zero), out, c.byref(self.n))

    def lapack(self, name, matrix, rhs=None):
        info = c.c_int(-999)
        args = [c.byref(self.lower), c.byref(self.n)]
        if rhs is not None:
            args += [c.byref(self.n), matrix, c.byref(self.n), rhs, c.byref(self.n)]
        else:
            args += [matrix, c.byref(self.n)]
        getattr(self.lib, name)(*args, c.byref(info))
        if info.value != 0:
            raise RuntimeError(name + ' info=' + str(info.value))

    def norm(self, value):
        return self.lib.dnrm2_(c.byref(self.total), value, c.byref(self.one))

    def relative_difference(self, actual, expected, scratch):
        c.memmove(scratch, actual, self.bytes)
        self.lib.daxpy_(c.byref(self.total), c.byref(self.minus), expected, c.byref(self.one), scratch, c.byref(self.one))
        return self.norm(scratch) / max(self.norm(expected), 1e-300)

    def symmetrize_lower(self, matrix):
        n = self.n.value
        base = c.addressof(matrix)
        for column in range(1, n):
            count = c.c_int(column)
            self.lib.dcopy_(c.byref(count), c.c_void_p(base + column * 8), c.byref(self.n),
                            c.c_void_p(base + column * n * 8), c.byref(self.one))

    def clear_upper(self, matrix):
        for column in range(self.n.value):
            c.memset(c.addressof(matrix) + column * self.n.value * 8, 0, column * 8)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--library', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--size', type=int, default=2560)
    parser.add_argument('--repetitions', type=int, default=2)
    parser.add_argument('--wall-seconds', type=int, default=360)
    args = parser.parse_args()
    if not 2500 <= args.size <= 4000 or not 1 <= args.repetitions <= 3 or not 30 <= args.wall_seconds <= 360:
        parser.error('Require n=2500..4000, repetitions=1..3, wall=30..360 s')
    sample_path = args.output.with_suffix('.samples.jsonl')
    if args.output.exists() or sample_path.exists():
        raise RuntimeError('Refuse to overwrite a benchmark output')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    # Virtual address reservations by pthread/BLAS arenas are not resident RAM.
    # Keep a distinct 8 GiB address ceiling and monitor actual RSS at 3 GiB.
    resource.setrlimit(resource.RLIMIT_AS, (8 * GIB,) * 2)
    resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
    os.environ['OPENBLAS_NUM_THREADS'] = '1'
    os.environ['OMP_NUM_THREADS'] = '1'
    started = time.monotonic()
    hz = os.sysconf('SC_CLK_TCK')
    live = {'phase': 'SETUP', 'threads': 1}
    done = threading.Event()
    report = {'schema': 'ISOLATED_DETERMINISTIC_OPENBLAS_KERNELS_V1', 'state': 'RUNNING',
        'pid': os.getpid(), 'source_sha256': sha(__file__), 'library': str(args.library.resolve()),
        'library_sha256': sha(args.library), 'affinity': sorted(os.sched_getaffinity(0)),
        'cgroup_ancestors': cgroup_info(), 'matrix_size': args.size, 'repetitions': args.repetitions,
        'thread_counts': list(THREADS), 'hard_address_space_bytes': 8 * GIB,
        'rss_limit_bytes': 3 * GIB, 'wall_limit_seconds': args.wall_seconds,
        'input_definition': 'Column-major A[k]=(((k*1664525+1013904223)&0xffffffff)/2**32)-0.5; '
                            'B[k]=(((k*22695477+1)&0xffffffff)/2**32)-0.5; S=I+A*A.T/n',
        'records': [], 'validation': [], 'scientific_gnss_acceptance': False}

    def save():
        args.output.write_text(json.dumps(report, indent=2))

    def monitor():
        with sample_path.open('x', buffering=1) as stream:
            while not done.is_set():
                try:
                    stat = Path('/proc/self/stat').read_text()
                    fields = stat[stat.rfind(')') + 2:].split()
                    rss = int(fields[21]) * os.sysconf('SC_PAGE_SIZE')
                    ticks, tids = cpu_snapshot()
                    row = {'elapsed_seconds': time.monotonic() - started, 'phase': live['phase'],
                        'threads_requested': live['threads'], 'rss_bytes': rss,
                        'process_cpu_ticks': ticks, 'thread_cpu_ticks': tids}
                    stream.write(json.dumps(row) + '\n')
                    if rss >= 3 * GIB or row['elapsed_seconds'] >= args.wall_seconds:
                        stream.write(json.dumps({'RESOURCE_STOP': row}) + '\n')
                        os.kill(os.getpid(), signal.SIGTERM)
                except BaseException as error:
                    stream.write(json.dumps({'MONITOR_FAILED': repr(error)}) + '\n')
                    os.kill(os.getpid(), signal.SIGTERM)
                done.wait(.1)

    watchdog = threading.Thread(target=monitor, daemon=True)
    watchdog.start()
    try:
        blas = Blas(args.library, args.size)
        blas.lib.openblas_set_num_threads(c.c_int(1))
        report.update(openblas_config=blas.config, openblas_corename=blas.lib.openblas_get_corename().decode(),
                      openblas_parallel_kind=int(blas.lib.openblas_get_parallel()))
        n, count = args.size, args.size ** 2
        matrix_type = c.c_double * count
        # Building arrays uses no random library or wall-clock-derived seed.
        aa = array('d', ((((k * 1664525 + 1013904223) & 0xffffffff) / 2 ** 32) - .5 for k in range(count)))
        bb = array('d', ((((k * 22695477 + 1) & 0xffffffff) / 2 ** 32) - .5 for k in range(count)))
        a, b = matrix_type.from_buffer(aa), matrix_type.from_buffer(bb)
        report['input_A_sha256'] = hashlib.sha256(aa).hexdigest()
        report['input_B_sha256'] = hashlib.sha256(bb).hexdigest()
        s, work, result, scratch, residual = (matrix_type() for _ in range(5))
        blas.gemm(a, a, s, transpose_b=True, scale=1. / n)
        for diagonal in range(n):
            s[diagonal + n * diagonal] += 1.
        report['input_S_sha256'] = hashlib.sha256(s).hexdigest()
        baselines = {}
        save()

        def kernel(name, fn, threads, repetition):
            live.update(phase=name, threads=threads)
            before, before_tids = cpu_snapshot()
            cpu_begin = time.process_time()
            begin = time.perf_counter()
            fn()
            wall = time.perf_counter() - begin
            cpu_elapsed = time.process_time() - cpu_begin
            after, after_tids = cpu_snapshot()
            thread_cpu = {tid: (value - before_tids.get(tid, value)) / hz for tid, value in after_tids.items()}
            report['records'].append({'kernel': name, 'threads_requested': threads,
                'threads_actual': int(blas.lib.openblas_get_num_threads()), 'repetition': repetition,
                'wall_seconds': wall, 'cpu_seconds': cpu_elapsed,
                'core_equivalent': cpu_elapsed / wall,
                'CPU_metric': 'CLOCK_PROCESS_CPUTIME high resolution across timed call, including monitor thread',
                'proc_tick_cpu_seconds': (after - before) / hz,
                'per_thread_cpu_seconds': thread_cpu,
                'per_thread_scope': 'Wider before/after proc snapshots; SC_CLK_TCK resolution; includes snapshot overhead'})
            live['phase'] = 'VALIDATION'

        for threads in THREADS:
            blas.lib.openblas_set_num_threads(c.c_int(threads))
            if blas.lib.openblas_get_num_threads() != threads:
                raise RuntimeError('Requested OpenBLAS thread count not applied')
            for name in ('DGEMM_NN', 'DGEMM_NT', 'POTRF', 'POTRS_MULTI_RHS', 'POTRI'):
                for repetition in range(args.repetitions):
                    if name.startswith('DGEMM'):
                        kernel(name, lambda: blas.gemm(a, b, result, name == 'DGEMM_NT'), threads, repetition)
                        target = result
                    elif name == 'POTRF':
                        c.memmove(work, s, blas.bytes)
                        kernel(name, lambda: blas.lapack('dpotrf_', work), threads, repetition)
                        target = work
                    elif name == 'POTRS_MULTI_RHS':
                        c.memmove(work, s, blas.bytes)
                        blas.lapack('dpotrf_', work)
                        c.memmove(result, b, blas.bytes)
                        kernel(name, lambda: blas.lapack('dpotrs_', work, result), threads, repetition)
                        target = result
                    else:
                        c.memmove(work, s, blas.bytes)
                        blas.lapack('dpotrf_', work)
                        kernel(name, lambda: blas.lapack('dpotri_', work), threads, repetition)
                        target = work
                    if name not in baselines:
                        baselines[name] = matrix_type()
                        c.memmove(baselines[name], target, blas.bytes)
                    difference = blas.relative_difference(target, baselines[name], scratch)
                    if not math.isfinite(difference) or difference > 5e-11:
                        raise RuntimeError(name + ' differs from the 1-thread baseline: ' + str(difference))
                    report['records'][-1]['relative_Frobenius_difference_vs_1_thread'] = difference
                # Independent operation residual, outside all timed calls.
                if name == 'POTRF':
                    c.memmove(result, work, blas.bytes)
                    blas.clear_upper(result)
                    blas.gemm(result, result, residual, transpose_b=True)
                    error = blas.relative_difference(residual, s, scratch)
                elif name == 'POTRS_MULTI_RHS':
                    blas.gemm(s, result, residual)
                    error = blas.relative_difference(residual, b, scratch)
                elif name == 'POTRI':
                    blas.symmetrize_lower(work)
                    blas.gemm(s, work, residual)
                    for diagonal in range(n):
                        residual[diagonal + n * diagonal] -= 1.
                    error = blas.norm(residual) / math.sqrt(n)
                else:
                    error = None
                if error is not None and (not math.isfinite(error) or error > 5e-10):
                    raise RuntimeError(name + ' residual failed: ' + str(error))
                report['validation'].append({'kernel': name, 'threads': threads,
                                              'relative_residual': error, 'passed': True})
                save()
        report['summary'] = []
        for name in ('DGEMM_NN', 'DGEMM_NT', 'POTRF', 'POTRS_MULTI_RHS', 'POTRI'):
            rows = [row for row in report['records'] if row['kernel'] == name]
            baseline = statistics.median(row['wall_seconds'] for row in rows if row['threads_requested'] == 1)
            for threads in THREADS:
                selected = [row for row in rows if row['threads_requested'] == threads]
                median = statistics.median(row['wall_seconds'] for row in selected)
                report['summary'].append({'kernel': name, 'threads': threads, 'median_wall_seconds': median,
                    'speedup_vs_1': baseline / median,
                    'median_core_equivalent': statistics.median(row['core_equivalent'] for row in selected)})
        report.update(state='PASSED', all_numerical_comparisons_passed=True)
    except BaseException as error:
        report.update(state='FAILED', failure=type(error).__name__ + ': ' + str(error))
        raise
    finally:
        done.set()
        watchdog.join(timeout=2)
        report.update(elapsed_seconds=time.monotonic() - started,
            peak_rss_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            finished_utc=datetime.now(timezone.utc).isoformat())
        save()
    print(json.dumps({'state': report['state'], 'elapsed_seconds': report['elapsed_seconds'],
                      'peak_rss_KiB': report['peak_rss_KiB'], 'output': str(args.output),
                      'summary': report['summary']}, indent=2))


if __name__ == '__main__':
    main()
