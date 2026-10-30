"""Compile header kernels with the production flags/libraries, no mock core.

BUILD OUTPUT [BEFORE_DIRECTORY] builds the current or frozen baseline driver.
Baseline is an explicit snapshot from the start of THIS task, not git HEAD.
"""
from pathlib import Path
import shlex
import shutil
import subprocess
import sys

root=Path(__file__).resolve().parents[2]
build=Path(sys.argv[1]).resolve();output=Path(sys.argv[2]).resolve()
output.mkdir(parents=True,exist_ok=True)
flags={}
for line in (build/'cpp/CMakeFiles/pea.dir/flags.make').read_text().splitlines():
    if ' = ' in line:
        key,value=line.split(' = ',1);flags[key]=shlex.split(value)
command=shlex.split((build/'cpp/CMakeFiles/pea.dir/link.txt').read_text())
extra=[]
if len(sys.argv)>3:
    before=Path(sys.argv[3]).resolve();include=output/'include/common';include.mkdir(parents=True,exist_ok=True)
    for path in before.glob('*.hpp'):shutil.copy2(path,include/path.name)
    extra=['-DCLOSURE_BASELINE=1','-I'+str(include.parent)]
executable=output/'product_closure_contract'
subprocess.run([command[0],*extra,*flags['CXX_DEFINES'],*flags['CXX_INCLUDES'],*flags['CXX_FLAGS'],
    '-O2','-std=c++20','-c',str(root/'src/tests/unit/r51_product_closure_contract.cpp'),
    '-o',str(output/'driver.o')],check=True)
command=[part for part in command if not part.endswith('.o') and not part.startswith('-Wl,--wrap=')]
command.insert(1,str(output/'driver.o'));command[command.index('-o')+1]=str(executable)
subprocess.run(command,cwd=build/'cpp',check=True)
print(executable)
