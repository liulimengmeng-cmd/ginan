"""Link active dependency regression to the real candidate PEA objects."""
from pathlib import Path
import os
import shlex
import subprocess
import sys

root = Path(__file__).resolve().parents[2]
build = Path(sys.argv[1]).resolve()
output = build / 'handoff_float_active_dependencies'
output.mkdir(exist_ok=True)
flags = {}
for line in (build / 'cpp/CMakeFiles/pea.dir/flags.make').read_text().splitlines():
    if ' = ' in line:
        key, value = line.split(' = ', 1)
        flags[key] = shlex.split(value)
command = shlex.split((build / 'cpp/CMakeFiles/pea.dir/link.txt').read_text())
subprocess.run([command[0], *flags['CXX_DEFINES'], *flags['CXX_INCLUDES'],
    *flags['CXX_FLAGS'], '-O1', '-std=c++20', '-c',
    str(root / 'src/tests/unit/r51_handoff_float_active_dependencies.cpp'),
    '-o', str(output / 'driver.o')], check=True)
subprocess.run(['objcopy', '--redefine-sym', 'main=pea_original_main',
    str(build / 'cpp/CMakeFiles/pea.dir/pea/main.cpp.o'),
    str(output / 'original_main.o')], check=True)
command[command.index('CMakeFiles/pea.dir/pea/main.cpp.o')] = str(output / 'original_main.o')
command.insert(1, str(output / 'driver.o'))
command[command.index('-o') + 1] = str(output / 'handoff_float_active_dependencies')
subprocess.run(command, cwd=build / 'cpp', check=True)
subprocess.run([str(output / 'handoff_float_active_dependencies')], check=True,
    env={**os.environ, 'ZHANG_R51_ENABLE': '1', 'OPENBLAS_NUM_THREADS': '1', 'OMP_NUM_THREADS': '1'})
