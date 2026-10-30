"""Link controlled G/E/C runtime tests to the actual built PEA objects.

Does not launch PEA, mutate existing objects, or read any RINEX observations.
Usage: python3 src/tests/run_multi_gnss_runtime.py BUILD [OUTPUT_DIRECTORY]
Optional: --nav3 RINEX3 --nav4 RINEX4 [--nav-only], or --syntax-only.
"""
from pathlib import Path
import json
import os
import shlex
import subprocess
import sys
import time

root = Path(__file__).resolve().parents[2]
syntax_only = '--syntax-only' in sys.argv
if syntax_only:
    sys.argv.remove('--syntax-only')
nav_only = '--nav-only' in sys.argv
if nav_only:
    sys.argv.remove('--nav-only')
nav_paths = []
for option in ('--nav3', '--nav4'):
    if option in sys.argv:
        index = sys.argv.index(option)
        nav_paths.append(str(Path(sys.argv[index + 1]).resolve()))
        del sys.argv[index:index + 2]
if nav_paths and len(nav_paths) != 2:
    raise SystemExit('Both --nav3 and --nav4 paths are required')
if nav_only and len(nav_paths) != 2:
    raise SystemExit('--nav-only requires --nav3 and --nav4')
build = Path(sys.argv[1]).resolve()
output = Path(sys.argv[2]).resolve() if len(sys.argv) > 2 else build / 'multi_gnss_runtime'
output.mkdir(parents=True, exist_ok=True)
flags = {}
for line in (build / 'cpp/CMakeFiles/pea.dir/flags.make').read_text().splitlines():
    if ' = ' in line:
        key, value = line.split(' = ', 1)
        flags[key] = shlex.split(value)
command = shlex.split((build / 'cpp/CMakeFiles/pea.dir/link.txt').read_text())
started = time.monotonic()
if syntax_only:
    result = subprocess.run([
        command[0], *flags['CXX_DEFINES'], *flags['CXX_INCLUDES'], *flags['CXX_FLAGS'],
        '-O1', '-fopenmp', '-std=c++20', '-fsyntax-only',
        str(root / 'src/tests/unit/multi_gnss_runtime.cpp')])
    raise SystemExit(result.returncode)
with (output / 'build.log').open('w') as log:
    subprocess.run([
        command[0], *flags['CXX_DEFINES'], *flags['CXX_INCLUDES'], *flags['CXX_FLAGS'],
        '-O1', '-fopenmp', '-std=c++20', '-c',
        str(root / 'src/tests/unit/multi_gnss_runtime.cpp'),
        '-o', str(output / 'driver.o')], check=True, stdout=log, stderr=subprocess.STDOUT)
    subprocess.run([
        'objcopy', '--redefine-sym', 'main=pea_original_main',
        str(build / 'cpp/CMakeFiles/pea.dir/pea/main.cpp.o'),
        str(output / 'original_main.o')], check=True, stdout=log, stderr=subprocess.STDOUT)
    command[command.index('CMakeFiles/pea.dir/pea/main.cpp.o')] = str(output / 'original_main.o')
    command.insert(1, str(output / 'driver.o'))
    command[command.index('-o') + 1] = str(output / 'multi_gnss_runtime')
    subprocess.run(command, cwd=build / 'cpp', check=True, stdout=log, stderr=subprocess.STDOUT)
with (output / 'runtime.log').open('w') as log:
    arguments = [str(output / 'multi_gnss_runtime'),
                 '--nav-only' if nav_only else str(output / 'runtime.TRACE'), *nav_paths]
    result = subprocess.run(arguments,
        stdout=log, stderr=subprocess.STDOUT,
        env={**os.environ, 'ZHANG_R51_ENABLE': '1', 'ZHANG_R51_RATIO_ONLY': '1',
             'ZHANG_R51_RATIO_THRESHOLD': '2', 'ZHANG_R51_EPOCH_ONLY_AR': '1',
             'OPENBLAS_NUM_THREADS': '1', 'OMP_NUM_THREADS': '1'})
summary = {
    'production_objects': str(build), 'exit_code': result.returncode,
    'elapsed_seconds': time.monotonic() - started,
    'synthetic_controlled_runtime_only': not nav_only,
    'navigation_parser_inputs': nav_paths,
    'navigation_parser_requested': bool(nav_paths),
    'scientific_ppp_run': False, 'full_hour_acceptance': False,
    'trace': str(output / 'runtime.TRACE'), 'runtime_log': str(output / 'runtime.log')}
(output / 'result.json').write_text(json.dumps(summary, indent=2) + '\n')
print((output / 'runtime.log').read_text())
print(json.dumps(summary))
raise SystemExit(result.returncode)
