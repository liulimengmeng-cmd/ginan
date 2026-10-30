"""Strict sequential frozen-baseline/current kernel comparisons.

Usage: BUILD EVIDENCE_DIRECTORY WORK_DIRECTORY [--replay INPUT_DIRECTORY]
Do not run while production compilation or a solver is competing for CPUs.
"""
from pathlib import Path
import hashlib
import json
import os
import re
import statistics
import subprocess
import sys
import time

root=Path(__file__).resolve().parents[2]
build,evidence,work=map(lambda x:Path(x).resolve(),sys.argv[1:4])
work.mkdir(parents=True,exist_ok=True)
env=dict(os.environ,OPENBLAS_NUM_THREADS='4',OMP_NUM_THREADS='4',MKL_NUM_THREADS='1',
         ZHANG_R51_RATIO_ONLY='1',ZHANG_R51_RATIO_THRESHOLD='2',ZHANG_R51_EPOCH_ONLY='1',ZHANG_R51_AR_THREADS='all')
env.pop('ZHANG_R51_EXACT_REPLAY_DIRECTORY',None)
results={'scope':'production header kernels; not full AR epoch timing','runs':[]}
result_name='replay_comparison.json' if len(sys.argv)>4 else 'kernel_comparison.json'
def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def save():
    (evidence/result_name).write_text(json.dumps(results,indent=2))
if len(sys.argv)>4:
    assert sys.argv[4]=='--replay'
    cases=[('replay',str(p)) for p in sorted(Path(sys.argv[5]).glob('frame_*.txt'))]
else:
    cases=[('tests',None)]+[('bench','800')]*5+[('marginal','2350')]
for name in ('baseline','current'):
    directory=work/name
    if len(sys.argv)<=4:
        command=[sys.executable,str(root/'src/tests/run_r51_product_closure_tests.py'),str(build),str(directory)]
        if name=='baseline':command.append(str(evidence/'before'))
        with (evidence/f'kernel_build_{name}.log').open('w') as log:
            subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True,env=env)
    results[name+'_binary_sha256']=digest(directory/'product_closure_contract')
for index,(mode,argument) in enumerate(cases):
    pair={};record={'mode':mode,'argument':argument,'index':index}
    prefix='replay' if mode=='replay' else 'kernel'
    for name in ('baseline','current'):
        output=work/f'{prefix}_{index}_{name}.txt'
        command=[str(work/name/'product_closure_contract'),mode,str(output)]
        if argument is not None:command.append(argument)
        start=time.monotonic()
        completed=subprocess.run(command,env=env,text=True,capture_output=True,timeout=360)
        (evidence/f'{prefix}_{index}_{name}.log').write_text(completed.stdout+completed.stderr)
        completed.check_returncode()
        record[name]={'wall_seconds':time.monotonic()-start,'output_sha256':digest(output),
                      'output_bytes':output.stat().st_size,'stdout':completed.stdout}
        key=re.search(r'(?:frame|marginal)_seconds=([0-9.e+-]+)',completed.stdout)
        if key:record[name]['kernel_seconds']=float(key[1])
        rss=re.search(r'peak_rss_kib=(\d+)',completed.stdout)
        if rss:record[name]['peak_rss_mib']=int(rss[1])/1024
        pair[name]=output
    if mode=='marginal':
        import numpy as np
        old=np.loadtxt(pair['baseline']);new=np.loadtxt(pair['current'])
        record['maximum_numeric_difference']=float(np.abs(old-new).max())
        record['relative_frobenius_difference']=float(np.linalg.norm(old-new)/np.linalg.norm(old))
        assert np.allclose(old,new,atol=1e-10,rtol=1e-10)
    else:
        assert digest(pair['baseline'])==digest(pair['current']),f'exact coordinate/witness mismatch: {argument}'
        record['byte_identical']=True
    results['runs'].append(record);save();print(json.dumps(record),flush=True)
if len(sys.argv)<=4:
    samples=[r for r in results['runs'] if r['mode']=='bench']
    results['frame_median_seconds']={n:statistics.median(r[n]['kernel_seconds'] for r in samples) for n in ('baseline','current')}
results['all_passed']=True;save()
