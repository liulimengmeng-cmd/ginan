"""Link production kernels with actual PEA build flags and libraries/objects.
Usage: python3 src/tests/run_r51_ar_hotspots.py BUILD OUTPUT [baseline|contracts|math|condition|checkpoint]
"""
from pathlib import Path
import shlex
import subprocess
import sys
import shutil

root=Path(__file__).resolve().parents[2]
build=Path(sys.argv[1]).resolve();output=Path(sys.argv[2]).resolve()
output.mkdir(parents=True,exist_ok=True)
flags={}
for line in (build/'cpp/CMakeFiles/pea.dir/flags.make').read_text().splitlines():
    if ' = ' in line:
        key,value=line.split(' = ',1);flags[key]=shlex.split(value)
command=shlex.split((build/'cpp/CMakeFiles/pea.dir/link.txt').read_text())
baseline=len(sys.argv)>3 and sys.argv[3]=='baseline'
condition_profile=len(sys.argv)>3 and sys.argv[3]=='condition'
checkpoint_probe=len(sys.argv)>3 and sys.argv[3]=='checkpoint'
math_only=condition_profile or (len(sys.argv)>3 and sys.argv[3]=='math')
contracts=math_only or (len(sys.argv)>3 and sys.argv[3]=='contracts')
profile=[] if baseline else ['-DAR_OPTIMIZATION_PROFILE=1']
if baseline:
    before=root.parent/'artifacts/ar_optimization_20261019/before'
    include=output/'baseline_include/common';include.mkdir(parents=True,exist_ok=True)
    for path in before.glob('*.hpp'):shutil.copy2(path,include/path.name)
    profile.insert(0,'-I'+str(include.parent))
subprocess.run([command[0],*profile,*flags['CXX_DEFINES'],*flags['CXX_INCLUDES'],*flags['CXX_FLAGS'],
    '-O2','-std=c++20','-c',str(root/'src/tests/unit'/
        ('r51_checkpoint_probe_verify.cpp' if checkpoint_probe else
         'r51_ar_condition_profile.cpp' if condition_profile else
         'r51_ar_optimization_contract.cpp' if contracts else 'r51_ar_hotspot_benchmark.cpp')),
    '-o',str(output/'driver.o')],check=True)
if math_only:
    command=[part for part in command if not part.endswith('.o') and not part.startswith('-Wl,--wrap=')]
    command.insert(1,str(output/'driver.o'))
    executable=output/('ar_condition_profile' if condition_profile else 'ar_optimization_contract')
    command[command.index('-o')+1]=str(executable)
    subprocess.run(command,cwd=build/'cpp',check=True)
    if not condition_profile:subprocess.run([str(executable)],check=True)
    print(executable)
    sys.exit(0)
if not baseline or not (output/'original_main.o').exists():
    subprocess.run(['objcopy','--redefine-sym','main=pea_original_main',
        str(build/'cpp/CMakeFiles/pea.dir/pea/main.cpp.o'),str(output/'original_main.o')],check=True)
if not baseline or not (output/'product_compiler.o').exists():
    shutil.copy2(build/'cpp/CMakeFiles/pea.dir/pea/ppp_ambres.cpp.o',output/'product_compiler.o')
command[command.index('CMakeFiles/pea.dir/pea/main.cpp.o')]=str(output/'original_main.o')
command[command.index('CMakeFiles/pea.dir/pea/ppp_ambres.cpp.o')]=str(output/'product_compiler.o')
command.insert(1,str(output/'driver.o'))
executable=output/('ar_checkpoint_verify' if checkpoint_probe else
                  'ar_optimization_contract' if contracts else 'ar_hotspot_benchmark')
command[command.index('-o')+1]=str(executable)
subprocess.run(command,cwd=build/'cpp',check=True)
print(executable)
