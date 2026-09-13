from pathlib import Path
import os, json, hashlib, datetime, subprocess, time, sys
w=Path(sys.argv[1] if len(sys.argv)>1 else '/mnt/c/Users/rx/Documents/GINAN/r48_user_model_work_20260914')
spec=json.loads((w/'experiment.json').read_text()); frozen=json.loads((w/'frozen_user.json').read_text())
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def now(): return datetime.datetime.now().astimezone().isoformat()
assert sha(frozen['binary'])==frozen['binary_sha256']
for name,item in json.loads((w/'product_manifest.json').read_text()).items(): assert sha(w/'products'/name)==item['sha256']
results=[]
for case in spec['cases']:
 out=Path(case['output']); cfg=w/'frozen/config'/Path(case['config']).name
 assert sha(cfg)==case['config_sha256']
 assert not out.exists(), str(out)
 out.mkdir()
 env=dict(os.environ)
 for key in ['ZHANG_R49_SNAPSHOT_DIRECTORY','ZHANG_R49_PARTIAL_SHADOW','ZHANG_BLAS_INJECT_INVALID','ZHANG_USER_MEASUREMENT_DUMP','ZHANG_USER_NAMED_ILS']: env.pop(key,None)
 env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
 if case['joint_ils']: env['ZHANG_USER_NAMED_ILS']='1'
 started=now(); clock=time.monotonic()
 with (out/'run.log').open('w') as log:
  rc=subprocess.run(['nice','-n','15',frozen['binary'],'-y',str(cfg),'-d',case['case'],'-r',case['observation']],cwd='/mnt/d/GINAN_R20/inputData',env=env,stdout=log,stderr=subprocess.STDOUT).returncode
 result={**case,'executed_config':str(cfg),'binary_sha256':frozen['binary_sha256'],'start':started,'end':now(),'wall_seconds':time.monotonic()-clock,'exit_code':rc,'threads':1}
 (w/(case['case']+'.result.json')).write_text(json.dumps(result,indent=2)+'\n')
 results.append(result)
 (w/'suite_progress.json').write_text(json.dumps({'completed_cases':len(results),'total_cases':len(spec['cases']),'results':results},indent=2)+'\n')
 print(case['case'],rc,round(result['wall_seconds'],2),flush=True)
assert sha(frozen['binary'])==frozen['binary_sha256']
for name,item in json.loads((w/'product_manifest.json').read_text()).items(): assert sha(w/'products'/name)==item['sha256']
(w/'suite_complete.json').write_text(json.dumps({'end':now(),'cases':len(results),'exit_zero':sum(x['exit_code']==0 for x in results)},indent=2)+'\n')
