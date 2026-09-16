"""Exercise production KF transactions after a Unix Makefiles PEA build.

Usage: python3 src/tests/run_r51_algebra_transport.py [build-directory]
The main object is copied and renamed; production objects and pea are unchanged.
"""
from pathlib import Path
import shlex
import os
import subprocess
import sys

root = Path(__file__).resolve().parents[2]
build = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else root / 'build'
output = build / 'transport_contract'
output.mkdir(exist_ok=True)
flags = {}
for line in (build / 'cpp/CMakeFiles/pea.dir/flags.make').read_text().splitlines():
    if ' = ' in line:
        key, value = line.split(' = ', 1)
        flags[key] = shlex.split(value)
subprocess.run([
    'c++', *flags['CXX_DEFINES'], *flags['CXX_INCLUDES'],
    '-O1', '-fopenmp', '-std=c++20', '-c',
    str(root / 'src/tests/unit/r51_algebra_transport_contract.cpp'),
    '-o', str(output / 'driver.o')], check=True)
subprocess.run([
    'objcopy', '--redefine-sym', 'main=pea_original_main',
    str(build / 'cpp/CMakeFiles/pea.dir/pea/main.cpp.o'),
    str(output / 'original_main.o')], check=True)
command = shlex.split((build / 'cpp/CMakeFiles/pea.dir/link.txt').read_text())
command[command.index('CMakeFiles/pea.dir/pea/main.cpp.o')] = str(output / 'original_main.o')
command.insert(1, str(output / 'driver.o'))
command[command.index('-o') + 1] = str(output / 'algebra_contract')
subprocess.run(command, cwd=build / 'cpp', check=True)
subprocess.run([str(output / 'algebra_contract')], check=True,
               env={**os.environ, 'ZHANG_R51_ENABLE': '1', 'OPENBLAS_NUM_THREADS': '1'})
