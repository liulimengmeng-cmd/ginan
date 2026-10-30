#!/usr/bin/env python3
"""Link a real-section COMPONENT test to all finished PEA production objects.

The exact rejected e135 core and source GEC bytes are immutable inputs. This
does not read a checkpoint bundle through a relaxed provenance path, execute
PEA main, claim nine-section restoration, or modify any producer object.
Execution must be wrapped by the caller's unchanged host/WSL resource guard.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import os
from pathlib import Path
import shlex
import subprocess

REJECTED = Path('/home/rx/GINAN/multi_gnss_float_profile_20261022/'
    'server134_fixedcoord_local_ar_e122_v6_pilot_20261025_01/restore_audit/'
    'ar_roundtrip_E29-20240717010700-e135-r32f626de9463-validation-only.rejected')
CORE = REJECTED / 'core.payload'
GEC = REJECTED / 'source_GEC_AR_RUNTIME_V1.payload'
CORE_BYTES = 425556082
CORE_SHA = '6163ae25308fec882ad2f63c967e53cce7a2e52b6d90d0ea5d3c1901bcb43cdf'
GEC_BYTES = 504018523
GEC_SHA = 'fb07adc9f827bd5d311dd50b35629d5cbdfaa16f76ef5b00df7c98530369cbde'
RUNTIME = 'E29-GEC134-FIXEDCOORD-20261023-one_hour'
SCOPE = 'REAL_GEC_SECTION_WITH_ACTUAL_KF_CORE_NOT_FULL_NINE_SECTION_RESTORE'


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def sha(path):
    value = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            value.update(chunk)
    return value.hexdigest()


def identity(path):
    return dict(path=str(path), bytes=Path(path).stat().st_size, sha256=sha(path))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('build', type=Path)
    parser.add_argument('--production-source-root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--dimension', type=int, required=True,
                        help='Exact e135 dimension from independently authenticated checkpoint manifest')
    parser.add_argument('--blas-library', type=Path)
    parser.add_argument('--runtime-timeout-seconds', type=int, default=1800)
    parser.add_argument('--syntax-only', action='store_true')
    parser.add_argument('--baseline-imprecise-red-only', '--legacy-gec-imprecise-json-red-only',
                        action='store_true', dest='baseline_imprecise_red_only')
    args = parser.parse_args()
    require(args.dimension == 7262 and 0 < args.runtime_timeout_seconds <= 3600,
            'Invalid bounded real component test dimension/deadline')
    repo = Path(__file__).resolve().parents[2]
    production = args.production_source_root.resolve(strict=True)
    build = args.build.resolve(strict=True)
    output = args.output.resolve()
    require(not output.exists(), 'Fresh real-component output directory required')
    output.mkdir(parents=True, exist_ok=False)
    cache = (build / 'CMakeCache.txt').read_text()
    require('CMAKE_HOME_DIRECTORY:INTERNAL=' + str(production / 'src') in cache.splitlines(),
            'Production objects belong to another immutable source tree')
    driver = repo / 'src/tests/unit/gec_ar_real_checkpoint_roundtrip.cpp'
    sources = [driver, Path(__file__).resolve(),
        *[production / relative for relative in ('src/cpp/common/zhangCheckpoint.hpp',
            'src/cpp/common/zhangCheckpoint.cpp', 'src/cpp/common/zhangCurrentEpochEvidence.hpp',
            'src/cpp/common/zhangIntegerDecisionProof.hpp', 'src/cpp/common/zhangMultiGnssInteger.hpp',
            'src/cpp/pea/ppp_ambres.cpp', 'src/cpp/pea/ppp_ambres_checkpoint.inc',
            'src/cpp/pea/zhangPppAr.hpp', 'src/cpp/pea/zhangPppAr.cpp')]]
    source_before = [identity(path) for path in sources]
    expected_inputs = [dict(path=str(CORE), bytes=CORE_BYTES, sha256=CORE_SHA),
                       dict(path=str(GEC), bytes=GEC_BYTES, sha256=GEC_SHA)]
    require([identity(path) for path in (CORE, GEC)] == expected_inputs,
            'Real rejected e135 core/GEC source identity differs')
    settings = {}
    for line in (build / 'cpp/CMakeFiles/pea.dir/flags.make').read_text().splitlines():
        key, separator, value = line.partition(' = ')
        if separator:
            settings[key] = shlex.split(value)
    link = shlex.split((build / 'cpp/CMakeFiles/pea.dir/link.txt').read_text())
    object_tokens = [token for token in link if token.endswith('.o')]
    require(len(object_tokens) == 158 and len(set(object_tokens)) == 158,
            'Complete real component must link all 158 distinct production objects')
    objects = [build / 'cpp' / token for token in object_tokens]
    object_before = [{**identity(path), 'mtime_ns': path.stat().st_mtime_ns} for path in objects]
    if args.blas_library:
        blas = args.blas_library.resolve(strict=True)
        indices = [index for index, item in enumerate(link) if 'openblas' in Path(item).name.lower()]
        require(indices, 'No explicit OpenBLAS production binding; no fallback allowed')
        for index in indices:
            link[index] = str(blas)
        link.append('-Wl,-rpath,' + str(blas.parent))
    prefix = [link[0], *settings['CXX_DEFINES'], *settings['CXX_INCLUDES'],
              *settings['CXX_FLAGS'], '-O0', '-std=c++20']
    if args.syntax_only:
        return subprocess.call([*prefix, '-fsyntax-only', str(driver)])
    driver_object = output / 'gec_ar_real_checkpoint_roundtrip.o'
    renamed_main = output / 'original_pea_main.o'
    executable = output / 'gec_ar_real_checkpoint_roundtrip'
    with (output / 'build.log').open('x') as stream:
        subprocess.run([*prefix, '-c', str(driver), '-o', str(driver_object)], stdout=stream,
                       stderr=subprocess.STDOUT, check=True)
        main_token = 'CMakeFiles/pea.dir/pea/main.cpp.o'
        require(link.count(main_token) == 1, 'Exactly one production PEA main object is required')
        subprocess.run(['objcopy', '--redefine-sym', 'main=pea_original_main',
                        str(build / 'cpp' / main_token), str(renamed_main)],
                       stdout=stream, stderr=subprocess.STDOUT, check=True)
        link[link.index(main_token)] = str(renamed_main)
        link.insert(1, str(driver_object)); link[link.index('-o') + 1] = str(executable)
        subprocess.run(link, cwd=build / 'cpp', stdout=stream, stderr=subprocess.STDOUT, check=True)
    require(source_before == [identity(path) for path in sources], 'Contract source changed during compile/link')
    require(object_before == [{**identity(path), 'mtime_ns': path.stat().st_mtime_ns} for path in objects],
            'Production objects changed during real-component link')
    raw_after = output / 'after_GEC_AR_RUNTIME_V1.payload'
    command = [str(executable), str(CORE), CORE_SHA, str(GEC), GEC_SHA, str(raw_after), str(args.dimension), RUNTIME]
    env = {key: value for key, value in os.environ.items() if not key.startswith('ZHANG_')}
    env.update(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
    with (output / 'runtime.log').open('x') as stream:
        result = subprocess.run(command, cwd=output, env=env, stdout=stream,
                                stderr=subprocess.STDOUT, timeout=args.runtime_timeout_seconds)
    marker = 'REAL_GEC_ROUNDTRIP_RESULT '
    records = [json.loads(line.partition(marker)[2]) for line in (output / 'runtime.log').read_text().splitlines()
               if line.startswith(marker)]
    native = records[0] if len(records) == 1 else None
    source_after = [identity(path) for path in sources]
    object_after = [{**identity(path), 'mtime_ns': path.stat().st_mtime_ns} for path in objects]
    inputs_after = [identity(path) for path in (CORE, GEC)]
    immutable = source_before == source_after and object_before == object_after and inputs_after == expected_inputs
    after_identity = identity(raw_after) if raw_after.is_file() else None
    expected_red = (args.baseline_imprecise_red_only and result.returncode == 1 and native is not None
        and native.get('scope') == SCOPE and native.get('stage') == 'COMPLETE_WIRE_COMPARE'
        and native.get('error') == 'GEC_AR_REAL_COMPONENT_WIRE_CHANGED_NO_TOLERANCE'
        and native.get('first_different_byte_offset') == 269945724 and immutable
        and native.get('actual_core_unchanged_after_import_and_export') is True
        and after_identity is not None and after_identity['sha256'] != GEC_SHA)
    passed = (not args.baseline_imprecise_red_only and result.returncode == 0 and native is not None
        and native.get('scope') == SCOPE and native.get('passed') is True and immutable
        and native.get('actual_production_roundtrips') == 2 and after_identity is not None
        and after_identity['bytes'] == GEC_BYTES and after_identity['sha256'] == GEC_SHA)
    report = dict(schema='REAL_GEC_COMPONENT_NATIVE_TEST_V1', scope=SCOPE, passed=passed,
        expected_imprecise_JSON_RED_observed=expected_red, baseline_imprecise_RED_only=args.baseline_imprecise_red_only,
        native_exit_code=result.returncode, native_result=native, command=command,
        linked_production_objects=158, production_source_root=str(production), production_cmake_home_exact=True,
        tested_files=source_before, tested_files_after=source_after, sources_unchanged=source_before == source_after,
        production_objects_before=object_before, production_objects_unchanged=object_before == object_after,
        actual_source_inputs=expected_inputs, actual_source_inputs_after=inputs_after,
        actual_source_inputs_unchanged=inputs_after == expected_inputs, raw_reexport=after_identity,
        executable=identity(executable), PEA_main_executed=False, numerical_epochs_consumed=0,
        full_nine_section_restore=False, checkpoint_binary_sha_gate_bypassed=False,
        scientific_precision_accepted=False)
    with (output / 'result.json').open('x') as stream:
        json.dump(report, stream, indent=2, allow_nan=False); stream.write('\n')
    print(json.dumps({key: report[key] for key in ('scope', 'passed', 'expected_imprecise_JSON_RED_observed',
        'native_exit_code', 'native_result', 'raw_reexport')}, indent=2))
    # Expected RED remains a deliberately failing contract, never acceptance.
    return 0 if passed else 1


if __name__ == '__main__':
    raise SystemExit(main())
