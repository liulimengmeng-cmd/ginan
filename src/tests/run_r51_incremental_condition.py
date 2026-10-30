"""Compile current/frozen conditioner independently with production flags.
BUILD OUTPUT [FROZEN_COMMON_DIRECTORY]. No PEA launch or original output writes.
"""
from pathlib import Path
import shlex
import subprocess
import sys

root=Path(__file__).resolve().parents[2]
build=Path(sys.argv[1]).resolve();out=Path(sys.argv[2]).resolve()
out.mkdir(parents=True,exist_ok=True)
flags={}
for line in (build/'cpp/CMakeFiles/pea.dir/flags.make').read_text().splitlines():
    if ' = ' in line:
        key,value=line.split(' = ',1);flags[key]=shlex.split(value)
command=shlex.split((build/'cpp/CMakeFiles/pea.dir/link.txt').read_text())
extra=[]
if len(sys.argv)>3:extra=['-DCONDITION_BASELINE=1','-I'+str(Path(sys.argv[3]).resolve().parent)]
exe=out/'incremental_condition_contract'
subprocess.run([command[0],*extra,*flags['CXX_DEFINES'],*flags['CXX_INCLUDES'],*flags['CXX_FLAGS'],
    '-O2','-std=c++20','-c',str(root/'src/tests/unit/r51_incremental_condition_contract.cpp'),'-o',str(out/'driver.o')],check=True)
command=[part for part in command if not part.endswith('.o') and not part.startswith('-Wl,--wrap=')]
command.insert(1,str(out/'driver.o'));command[command.index('-o')+1]=str(exe)
subprocess.run(command,cwd=build/'cpp',check=True)
print(exe)
