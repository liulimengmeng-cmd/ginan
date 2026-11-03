"""Alternating frozen/current production-kernel runs. No GNSS solver launch.
ROOT contains current/{incremental_condition_contract,float_parallel_benchmark},
baseline_condition/incremental_condition_contract, frozen_float/float_parallel_benchmark.
"""
from pathlib import Path
import json
import os
import re
import statistics
import subprocess
import sys
import numpy as np

root=Path(sys.argv[1]).resolve()
if len(sys.argv)>2 and sys.argv[2]=='normalize-missing-metrics':
    path=root/'paired_kernels/results.json'
    result=json.loads(path.read_text())
    for name,records in result['observations'].items():
        if not any('peak_rss_kib' in row for row in records):
            result['summary'][name]['peak_rss_mib']=None
    result['missing_metric_policy']='RSS absent from a kernel record is unavailable, not zero.'
    path.write_text(json.dumps(result,indent=2))
    print('Normalized only unavailable RSS; timing observations and comparisons unchanged.')
    sys.exit(0)
out=root/(sys.argv[2] if len(sys.argv)>2 else 'paired_kernels');out.mkdir(exist_ok=False)
environment={**os.environ,'ZHANG_R51_AR_THREADS':'all','OMP_NUM_THREADS':'12',
             'OPENBLAS_NUM_THREADS':'12','MKL_NUM_THREADS':'1',
             'ZHANG_R51_BLOCK_CONDITIONER':'incremental','ZHANG_R51_FLOAT_MATRIX_MODE':'symmetric'}
observations={}
def run(name,command):
    result=subprocess.run(list(map(str,command)),env=environment,text=True,
        stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    number=len(observations.get(name,[]))+1
    (out/f'{name}_{number}.log').write_text(result.stdout)
    assert result.returncode==0,f'{name} failed; see preserved log'
    match=re.search(r'^(?:BENCH|CONDITION) (.+)$',result.stdout,re.M)
    assert match
    values={k:float(v) for k,v in re.findall(r'(\w+)=([0-9.eE+-]+)',match[1])}
    observations.setdefault(name,[]).append(values)
    print(name,number,values.get('wall_s',values.get('seconds')),flush=True)

for repeat in range(5):
    order=['old_condition','new_condition']
    if repeat%2:order.reverse()
    for name in order:
        binary=root/('baseline_condition' if name=='old_condition' else 'current')/'incremental_condition_contract'
        run(name,[binary,1175,30,32])

float_old=root/'frozen_float/float_parallel_benchmark'
float_new=root/'current/float_parallel_benchmark'
for repeat in range(3):
    order=[1,2,4,6,8,12];order=order[repeat:]+order[:repeat]
    for threads in order:
        run(f'float_grid_{threads}',[float_new,threads,2350,2800])

for repeat in range(5):
    order=['old_float_12','new_float_12']
    if repeat%2:order.reverse()
    for name in order:
        command=[float_old if name=='old_float_12' else float_new,12,4755,5600]
        if repeat==0:command.append(out/(name+'.posterior'))
        run(name,command)

a=np.fromfile(out/'old_float_12.posterior',dtype=np.float64)
b=np.fromfile(out/'new_float_12.posterior',dtype=np.float64)
assert len(a)==len(b)==4755+4755**2 and np.isfinite(a).all() and np.isfinite(b).all()
error=np.abs(a-b)
comparison=dict(entries=len(a),mean_maximum_difference=float(error[:4755].max()),
    covariance_maximum_difference=float(error[4755:].max()),tolerance=2e-11)
assert error.max()<comparison['tolerance'],'full FLOAT kernel posterior mismatch'
summary={}
for name,records in observations.items():
    times=[r.get('wall_s',r.get('seconds')) for r in records]
    summary[name]=dict(repetitions=len(times),median_seconds=statistics.median(times),
        minimum_seconds=min(times),maximum_seconds=max(times),
        peak_rss_mib=max(r['peak_rss_kib'] for r in records)/1024 if
            all('peak_rss_kib' in r for r in records) else None)
result=dict(scope='Synthetic inputs, real production kernels; not a full-network throughput test',
    settings=dict(ratio=2,joint_nis=False,epoch_only_ar=True),
    summary=summary,full_float_comparison=comparison,observations=observations,passed=True)
(out/'results.json').write_text(json.dumps(result,indent=2))
print(json.dumps(summary,indent=2),flush=True)
