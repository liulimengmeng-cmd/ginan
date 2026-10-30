"""Production-flag and production-object-linked checkpoint policy contract.

Never launches a network/PPP/AR solver. Use --syntax-only while PEA compiles;
run the full link only after the producer build and policy header are stable.
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
parser.add_argument('--output', type=Path)
parser.add_argument('--syntax-only', action='store_true')
parser.add_argument('--build-only', action='store_true')
parser.add_argument('--generic-blas', type=Path)
args = parser.parse_args()
root = Path(__file__).resolve().parents[2]
build = args.build.resolve()
assert 'CMAKE_HOME_DIRECTORY:INTERNAL=' + str(root / 'src') in (build / 'CMakeCache.txt').read_text(), 'Wrong checkout objects'
flags = {}
for line in (build / 'cpp/CMakeFiles/pea.dir/flags.make').read_text().splitlines():
    if ' = ' in line:
        key, value = line.split(' = ', 1)
        flags[key] = shlex.split(value)
link = shlex.split((build / 'cpp/CMakeFiles/pea.dir/link.txt').read_text())
prefix = [link[0], *flags['CXX_DEFINES'], *flags['CXX_INCLUDES'], *flags['CXX_FLAGS'], '-O0', '-std=c++20']
source = root / 'src/tests/unit/multi_gnss_float_checkpoint_policy.cpp'
if args.syntax_only:
    raise SystemExit(subprocess.call([*prefix, '-fsyntax-only', str(source)]))
output = args.output.resolve() if args.output else build / 'multi_gnss_float_checkpoint_policy_contract'
output.mkdir(parents=True, exist_ok=True)

def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()

sources = [source, Path(__file__).resolve(), root / 'src/cpp/common/zhangMultiGnssFloatCheckpoint.hpp',
           root / 'src/cpp/common/zhangR51Warmup.hpp', root / 'src/cpp/common/acsConfig.hpp',
           root / 'src/cpp/pea/zhangCheckpointRuntime.cpp', root / 'src/cpp/pea/zhangCheckpointRuntime.hpp']
source_pins = {str(path): sha(path) for path in sources}
objects = [(build / 'cpp' / token).resolve() for token in link if token.endswith('.o')]
object_pins = {str(path): sha(path) for path in objects}
generic = args.generic_blas.resolve() if args.generic_blas else None
replacements = 0
if generic:
    replacements = sum(token.endswith('/libopenblas.a') for token in link)
    assert replacements, 'Original static BLAS library not found in link command'
    link = [str(generic) if token.endswith('/libopenblas.a') else token for token in link]
    link.insert(1, '-Wl,-rpath,' + str(generic.parent))
executable = output / 'multi_gnss_float_checkpoint_policy'
with (output / 'build.log').open('w') as log:
    subprocess.run([*prefix, '-c', str(source), '-o', str(output / 'driver.o')],
                   stdout=log, stderr=subprocess.STDOUT, check=True)
    subprocess.run(['objcopy', '--redefine-sym', 'main=pea_original_main',
                    str(build / 'cpp/CMakeFiles/pea.dir/pea/main.cpp.o'), str(output / 'original_main.o')],
                   stdout=log, stderr=subprocess.STDOUT, check=True)
    link[link.index('CMakeFiles/pea.dir/pea/main.cpp.o')] = str(output / 'original_main.o')
    link.insert(1, str(output / 'driver.o'))
    link[link.index('-o') + 1] = str(executable)
    subprocess.run(link, cwd=build / 'cpp', stdout=log, stderr=subprocess.STDOUT, check=True)
assert object_pins == {str(path): sha(path) for path in objects}, 'Producer objects changed during policy link; retry'
assert source_pins == {str(path): sha(path) for path in sources}, 'Policy/ABI sources changed during contract build'
report = {'scope': 'Actual-flags production-linked policy only; no solver launch',
          'executable': str(executable), 'binary_sha256': sha(executable), 'sources_sha256': source_pins,
          'production_object_sha256': object_pins, 'inputs_unchanged': True, 'link_command': link,
          'generic_blas': str(generic) if generic else None,
          'generic_blas_sha256': sha(generic) if generic else None,
          'static_blas_replacements': replacements}
if args.build_only:
    (output / 'result.json').write_text(json.dumps(report, indent=2) + '\n')
    print(executable)
    raise SystemExit(0)
with (output / 'runtime.log').open('w') as log:
    result = subprocess.run([str(executable)], stdout=log, stderr=subprocess.STDOUT,
                            env={**os.environ, 'OMP_NUM_THREADS': '1', 'OPENBLAS_NUM_THREADS': '1'})
report['exit_code'] = result.returncode
report['passed'] = result.returncode == 0
(output / 'result.json').write_text(json.dumps(report, indent=2) + '\n')
print((output / 'runtime.log').read_text())
print(json.dumps({key: value for key, value in report.items() if key not in ('production_object_sha256','link_command')}, indent=2))
raise SystemExit(result.returncode)
