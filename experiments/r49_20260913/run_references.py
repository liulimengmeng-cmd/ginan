from pathlib import Path
import subprocess,os,json
w=Path(__file__).resolve().parent;r=Path('/home/rx/GINAN/r49-condition-domain-20260913')
prefix=Path('/home/rx/GINAN/vcpkg_installed/linux/x64-linux')
env={**os.environ,'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1'}
results=[]
for name in ['selection','posterior','domain']:
 target='r49_'+name+'_reference'
 command=['g++','-std=c++20','-O1','-g0','-DEIGEN_USE_BLAS=1','-I'+str(r/'src/cpp'),'-I'+str(prefix/'include/eigen3'),'-I'+str(prefix/'include'),str(r/'src/tests/unit'/f'{target}.cpp'),str(r/'src/tests/unit/zhangBlasArgumentAudit.cpp'),str(prefix/'lib/libopenblas.a'),'-Wl,--wrap=dgemv_','-lpthread','-lgfortran','-o',str(r/target)]
 with (w/(target+'_build.log')).open('w') as f:rc=subprocess.run(command,stdout=f,stderr=subprocess.STDOUT,env=env).returncode
 assert rc==0,target
 with (w/(target+'.log')).open('w') as f:rc=subprocess.run([str(r/target)],stdout=f,stderr=subprocess.STDOUT,env=env).returncode
 assert rc==0,target
 results.append({'test':target,'exit':rc,'backend':'EIGEN_USE_BLAS=1 + actual LP64 OpenBLAS + DGEMV wrapper'})
 print(target+' PASS',flush=True)
(w/'reference_validation.json').write_text(json.dumps(results,indent=2)+'\n')
