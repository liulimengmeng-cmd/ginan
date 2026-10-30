"""Link only a test driver against stable production PEA objects; no solver."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shlex
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument('build', type=Path)
parser.add_argument('--output', type=Path, required=True)
parser.add_argument('--syntax-only', action='store_true')
parser.add_argument('--blas-library', type=Path)
args = parser.parse_args()
repo = Path(__file__).resolve().parents[2]
build, output = args.build.resolve(), args.output.resolve()
output.mkdir(parents=True, exist_ok=True)
flags = {}
for line in (build / 'cpp/CMakeFiles/pea.dir/flags.make').read_text().splitlines():
    key, separator, value = line.partition(' = ')
    if separator:
        flags[key] = shlex.split(value)
link = shlex.split((build / 'cpp/CMakeFiles/pea.dir/link.txt').read_text())
if args.blas_library:
    library = args.blas_library.resolve(strict=True)
    selected = [i for i, value in enumerate(link) if 'openblas' in Path(value).name.lower()]
    if not selected:
        raise RuntimeError('No explicit OpenBLAS link argument')
    for index in selected:
        link[index] = str(library)
    link.append('-Wl,-rpath,' + str(library.parent))
driver = repo / 'src/tests/unit/multi_gnss_linear_combo_contract.cpp'
compile_prefix = [link[0], *flags['CXX_DEFINES'], *flags['CXX_INCLUDES'], *flags['CXX_FLAGS'], '-O0']
if args.syntax_only:
    raise SystemExit(subprocess.call([*compile_prefix, '-fsyntax-only', str(driver)]))
objects = [build / 'cpp' / value for value in link if value.endswith('.o')]
def stamps():
    return {str(path.relative_to(build)): [path.stat().st_size, path.stat().st_mtime_ns] for path in objects}
before = stamps()
with (output / 'build.log').open('w') as log:
    obj, renamed = output / 'driver.o', output / 'original_pea_main.o'
    subprocess.run([*compile_prefix, '-c', str(driver), '-o', str(obj)], check=True,
                   stdout=log, stderr=subprocess.STDOUT)
    subprocess.run(['objcopy', '--redefine-sym', 'main=pea_original_main',
                    str(build / 'cpp/CMakeFiles/pea.dir/pea/main.cpp.o'), str(renamed)],
                   check=True, stdout=log, stderr=subprocess.STDOUT)
    link[link.index('CMakeFiles/pea.dir/pea/main.cpp.o')] = str(renamed)
    link.insert(1, str(obj))
    executable = output / 'linear_combo_contract'
    link[link.index('-o') + 1] = str(executable)
    subprocess.run(link, cwd=build / 'cpp', stdout=log, stderr=subprocess.STDOUT, check=True)
assert before == stamps(), 'Production objects changed during test link'
env = {key: value for key, value in os.environ.items() if not key.startswith('ZHANG_')}
env.update(OMP_NUM_THREADS='8', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
with (output / 'runtime.log').open('w') as log:
    completed = subprocess.run([str(executable)], cwd=output, env=env, timeout=90,
                               stdout=log, stderr=subprocess.STDOUT)
marker = 'LINEAR_COMBO_CONTRACT_RESULT '
results = [json.loads(line.partition(marker)[2]) for line in (output / 'runtime.log').read_text().splitlines() if marker in line]
report = {'passed': completed.returncode == 0 and bool(results), 'exit_code': completed.returncode,
          'result': results[-1] if results else None,
          'producer_objects_unchanged': before == stamps(), 'production_object_stamps': before,
          'driver_sha256': hashlib.sha256(driver.read_bytes()).hexdigest(),
          'selected_blas_library': str(args.blas_library.resolve()) if args.blas_library else None,
          'executable_sha256': hashlib.sha256(executable.read_bytes()).hexdigest()}
(output / 'result.json').write_text(json.dumps(report, indent=2) + '\n')
compact = dict(report['result']) if report['result'] else None
if compact:
    compact.pop('oracle_rows', None)
    compact.pop('field_order', None)
print(json.dumps({'passed': report['passed'], 'exit_code': report['exit_code'], 'result': compact}, indent=2))
if not report['passed']:
    print((output / 'runtime.log').read_text())
raise SystemExit(0 if report['passed'] else 1)
