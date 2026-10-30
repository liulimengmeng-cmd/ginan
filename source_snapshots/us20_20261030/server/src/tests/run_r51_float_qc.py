"""Exercise actual production FLOAT transaction rollback after a PEA build.

Usage: python3 src/tests/run_r51_float_qc.py <Unix-Makefiles-build-directory>
Only a copy of the main object is renamed; production objects and pea are intact.
"""
from pathlib import Path
import os
import shlex
import subprocess
import sys

root = Path(__file__).resolve().parents[2]
build = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else root / 'build'
output = build / 'float_qc_contract'
output.mkdir(exist_ok=True)
flags = {}
for line in (build / 'cpp/CMakeFiles/pea.dir/flags.make').read_text().splitlines():
    if ' = ' in line:
        key, value = line.split(' = ', 1)
        flags[key] = shlex.split(value)
command = shlex.split((build / 'cpp/CMakeFiles/pea.dir/link.txt').read_text())
subprocess.run([
    command[0], *flags['CXX_DEFINES'], *flags['CXX_INCLUDES'], *flags['CXX_FLAGS'],
    '-O1', '-fopenmp', '-std=c++20', '-c',
    str(root / 'src/tests/unit/r51_float_qc_transaction.cpp'),
    '-o', str(output / 'driver.o')], check=True)
subprocess.run([
    'objcopy', '--redefine-sym', 'main=pea_original_main',
    str(build / 'cpp/CMakeFiles/pea.dir/pea/main.cpp.o'),
    str(output / 'original_main.o')], check=True)
command[command.index('CMakeFiles/pea.dir/pea/main.cpp.o')] = str(output / 'original_main.o')
command.insert(1, str(output / 'driver.o'))
command[command.index('-o') + 1] = str(output / 'float_qc_contract')
subprocess.run(command, cwd=build / 'cpp', check=True)
subprocess.run([str(output / 'float_qc_contract')], check=True,
    env={**os.environ, 'OPENBLAS_NUM_THREADS': '1', 'OMP_NUM_THREADS': '1'})
