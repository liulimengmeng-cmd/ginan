"""Link a read-only boundary-comparator test to completed production objects.

No production source/object is compiled or modified here and PEA is not run.
The raw source and optional actual post-restore payload are hash-bound.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shlex
import subprocess

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('build', type=Path)
parser.add_argument('--output', required=True, type=Path)
parser.add_argument('--source-payload', required=True, type=Path)
parser.add_argument('--after-payload', type=Path)
parser.add_argument('--blas-library', type=Path)
args = parser.parse_args()
repo = Path(__file__).resolve().parents[2]
build = args.build.resolve(strict=True)
output = args.output.resolve()
output.mkdir(parents=True, exist_ok=False)
settings = {}
for line in (build / 'cpp/CMakeFiles/pea.dir/flags.make').read_text().splitlines():
    key, separator, value = line.partition(' = ')
    if separator:
        settings[key] = shlex.split(value)
link = shlex.split((build / 'cpp/CMakeFiles/pea.dir/link.txt').read_text())
if args.blas_library:
    library = args.blas_library.resolve(strict=True)
    selected = [i for i, part in enumerate(link) if 'openblas' in Path(part).name.lower()]
    if not selected:
        raise RuntimeError('No explicit OpenBLAS link argument; no fallback is permitted')
    for index in selected:
        link[index] = str(library)
    link.append('-Wl,-rpath,' + str(library.parent))
driver = repo / 'src/tests/unit/receiver_boundary_reset_equivalence.cpp'
sources = [driver, Path(__file__).resolve(), repo / 'src/cpp/pea/zhangReceiverCheckpoint.cpp',
           repo / 'src/cpp/pea/zhangReceiverCheckpoint.hpp',
           repo / 'src/tests/unit/multi_gnss_receiver_runtime_audit_schema.hpp',
           repo / 'src/cpp/common/zhangCheckpoint.hpp', repo / 'src/cpp/common/zhangCheckpoint.cpp']
objects = [build / 'cpp' / part for part in link if part.endswith('.o')]


def identity(path):
    path = Path(path)
    return {'path': str(path), 'bytes': path.stat().st_size,
            'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}


def stamps():
    return {str(path): [path.stat().st_size, path.stat().st_mtime_ns] for path in objects}


source_bindings = [identity(path) for path in sources]
payloads = [args.source_payload.resolve(strict=True)]
if args.after_payload:
    payloads.append(args.after_payload.resolve(strict=True))
payload_bindings = [identity(path) for path in payloads]
before = stamps()
object_path = output / 'receiver_boundary_reset_equivalence.o'
renamed_main = output / 'original_pea_main.o'
executable = output / 'receiver_boundary_reset_equivalence'
with (output / 'build.log').open('w') as log:
    subprocess.run([link[0], *settings['CXX_DEFINES'], *settings['CXX_INCLUDES'],
        *settings['CXX_FLAGS'], '-O0', '-std=c++20', '-c', str(driver), '-o', str(object_path)],
        stdout=log, stderr=subprocess.STDOUT, check=True)
    subprocess.run(['objcopy', '--redefine-sym', 'main=pea_original_main',
        str(build / 'cpp/CMakeFiles/pea.dir/pea/main.cpp.o'), str(renamed_main)],
        stdout=log, stderr=subprocess.STDOUT, check=True)
    link[link.index('CMakeFiles/pea.dir/pea/main.cpp.o')] = str(renamed_main)
    link.insert(1, str(object_path))
    link[link.index('-o') + 1] = str(executable)
    subprocess.run(link, cwd=build / 'cpp', stdout=log, stderr=subprocess.STDOUT, check=True)
if stamps() != before:
    raise RuntimeError('Production object stamps changed during test-only linking')
env = {key: value for key, value in os.environ.items() if not key.startswith('ZHANG_')}
env.update(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
with (output / 'runtime.log').open('w') as log:
    native = subprocess.run([str(executable), *map(str, payloads)], cwd=output,
        env=env, stdout=log, stderr=subprocess.STDOUT, timeout=180)
marker = 'RECEIVER_BOUNDARY_RESET_TEST_RESULT '
records = [json.loads(line[len(marker):]) for line in (output / 'runtime.log').read_text().splitlines()
           if line.startswith(marker)]
unchanged = stamps() == before and source_bindings == [identity(path) for path in sources] and \
    payload_bindings == [identity(path) for path in payloads]
report = {'schema': 'RECEIVER_BOUNDARY_RESET_NATIVE_TEST_V1', 'native_exit_code': native.returncode,
          'production_objects_sources_and_payloads_unchanged': unchanged,
          'source_bindings': source_bindings, 'payload_bindings': payload_bindings,
          'executable': identity(executable), 'result': records[-1] if records else None,
          'passed': native.returncode == 0 and len(records) == 1 and records[0]['status'] == 'PASS' and unchanged,
          'PEA_epoch_executed': False}
(output / 'result.json').write_text(json.dumps(report, indent=2) + '\n')
print(json.dumps(report, indent=2))
raise SystemExit(0 if report['passed'] else 1)
