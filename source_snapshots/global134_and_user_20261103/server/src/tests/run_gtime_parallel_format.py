"""Actual-flags production-linked shared GTime formatting contract.

The serial oracle is a separate executable linked to the preserved old
production gTime.cpp.o, not another call to the new implementation.
Never launches PPP or AR and never replaces production objects or binaries.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shlex
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument('build', type=Path)
parser.add_argument('--output', required=True, type=Path)
parser.add_argument('--legacy-object', required=True, type=Path)
parser.add_argument('--generic-blas', required=True, type=Path)
args = parser.parse_args()
root = Path(__file__).resolve().parents[2]
build = args.build.resolve(strict=True)
output = args.output.resolve()
output.mkdir(exist_ok=False, parents=True)
legacy = args.legacy_object.resolve(strict=True)
generic = args.generic_blas.resolve(strict=True)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


require(sha(legacy) == '5d74b4f029f0220b6503e9f244e578d8e9bb9b71d804738ccff95fea8999fea8',
        'Wrong preserved production legacy gTime object')
require('CMAKE_HOME_DIRECTORY:INTERNAL=' + str(root / 'src') in
        (build / 'CMakeCache.txt').read_text(), 'Wrong production checkout')
flags = {}
for line in (build / 'cpp/CMakeFiles/pea.dir/flags.make').read_text().splitlines():
    if ' = ' in line:
        key, value = line.split(' = ', 1)
        flags[key] = shlex.split(value)
link = shlex.split((build / 'cpp/CMakeFiles/pea.dir/link.txt').read_text())
source = root / 'src/tests/unit/gtime_parallel_format.cpp'
sources = [source, Path(__file__).resolve(), root / 'src/cpp/common/gTime.cpp',
           root / 'src/cpp/common/gTime.hpp']
source_pins = {str(path): sha(path) for path in sources}
objects = [(build / 'cpp' / token).resolve(strict=True) for token in link if token.endswith('.o')]
object_pins = {str(path): sha(path) for path in objects}
libraries = [(build / 'cpp' / token).resolve(strict=True) for token in link if token.endswith('.a')]
library_pins = {str(path): sha(path) for path in libraries}
replacements = sum(token.endswith('/libopenblas.a') for token in link)
require(replacements == 2, 'Unexpected original static BLAS link inputs')
link = [str(generic) if token.endswith('/libopenblas.a') else token for token in link]
link.insert(1, '-Wl,-rpath,' + str(generic.parent))
prefix = [link[0], *flags['CXX_DEFINES'], *flags['CXX_INCLUDES'],
          *flags['CXX_FLAGS'], '-std=c++20']
with (output / 'build.log').open('x') as log:
    subprocess.run([*prefix, '-c', str(source), '-o', str(output / 'driver.o')],
                   stdout=log, stderr=subprocess.STDOUT, check=True)
    subprocess.run(['objcopy', '--redefine-sym', 'main=pea_original_main',
                    str(build / 'cpp/CMakeFiles/pea.dir/pea/main.cpp.o'), str(output / 'original_main.o')],
                   stdout=log, stderr=subprocess.STDOUT, check=True)
    link[link.index('CMakeFiles/pea.dir/pea/main.cpp.o')] = str(output / 'original_main.o')
    link.insert(1, str(output / 'driver.o'))
    new_exe = output / 'gtime_parallel_format'
    old_exe = output / 'gtime_legacy_serial_oracle'
    link[link.index('-o') + 1] = str(new_exe)
    subprocess.run(link, cwd=build / 'cpp', stdout=log, stderr=subprocess.STDOUT, check=True)
    old_link = list(link)
    old_link[old_link.index('CMakeFiles/pea.dir/common/gTime.cpp.o')] = str(legacy)
    old_link[old_link.index('-o') + 1] = str(old_exe)
    subprocess.run(old_link, cwd=build / 'cpp', stdout=log, stderr=subprocess.STDOUT, check=True)

env = {**os.environ, 'OMP_NUM_THREADS': '1', 'OPENBLAS_NUM_THREADS': '1'}
reports = {}
raw_outputs = {}
for name, command in [('legacy_serial', [str(old_exe), '--serial-only']),
                      ('new_serial', [str(new_exe), '--serial-only']),
                      ('legacy_serial_bytes', [str(old_exe), '--dump-serial']),
                      ('new_serial_bytes', [str(new_exe), '--dump-serial']),
                      ('new_parallel_poison', [str(new_exe)])]:
    result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            env=env, timeout=120, text=True)
    (output / (name + '.log')).write_text(result.stdout)
    require(result.returncode == 0, name + ' failed: ' + result.stdout)
    raw_outputs[name] = result.stdout
    reports[name] = json.loads(result.stdout.strip().splitlines()[-1])
require(reports['legacy_serial']['serial_fingerprint'] ==
        reports['new_serial']['serial_fingerprint'] ==
        reports['new_parallel_poison']['serial_fingerprint'], 'Old/new serial formatting drift')
require(raw_outputs['legacy_serial_bytes'] == raw_outputs['new_serial_bytes'],
        'Entrywise old/new serial output byte mismatch')
require(object_pins == {str(path): sha(path) for path in objects}, 'Producer objects changed during contract')
require(source_pins == {str(path): sha(path) for path in sources}, 'Sources changed during contract')
require(library_pins == {str(path): sha(path) for path in libraries}, 'Static libraries changed during contract')
report = {'passed': True, 'scope': 'Actual production object link; independent preserved old-object serial oracle; no solver',
          'source_sha256': source_pins, 'production_object_sha256': object_pins,
          'static_library_sha256': library_pins, 'legacy_object_sha256': sha(legacy),
          'new_executable_sha256': sha(new_exe), 'legacy_oracle_executable_sha256': sha(old_exe),
          'generic_blas_sha256': sha(generic), 'static_blas_replacements': replacements,
          'production_inputs_unchanged': True, 'results': reports,
          'serial_entrywise_byte_equal': True,
          'serial_output_bytes': len(raw_outputs['legacy_serial_bytes'].encode()),
          'serial_output_sha256': hashlib.sha256(raw_outputs['legacy_serial_bytes'].encode()).hexdigest(),
          'limitations': ['Fresh legacy time=-1 precision0 cache-sentinel bug is deliberately not preserved.',
              'Existing long-int rounding overflow for modern epochs at precision>=10 is outside this minimal fix; high-precision equality is empirical, not an undefined-behavior safety claim.',
              'No concurrent writes to bigTime are authorized; the contract uses shared immutable epoch values.'],
          'link_command': link}
(output / 'result.json').write_text(json.dumps(report, indent=2) + '\n')
print(json.dumps({key: value for key, value in report.items()
                  if key not in ('production_object_sha256', 'static_library_sha256', 'link_command')}, indent=2))
