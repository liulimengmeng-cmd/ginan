#!/usr/bin/env python3
"""Run frozen R48 user cases serially with exact provenance and no overwrite."""
import argparse,datetime,hashlib,json,os,signal,subprocess,time
from pathlib import Path
WORK=Path('/mnt/c/Users/rx/Documents/GINAN/r48_uduc_work_20260913')
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--case');ap.add_argument('--spec',default='experiment.json');args=ap.parse_args()
 s=json.loads((WORK/args.spec).read_text())
 assert sha(s['binary'])==s['binary_sha256']
 for x in s['snapshot'].values():assert sha(x['path'])==x['sha256']
 for st in s['stations']:assert sha(st['path'])==st['sha256']
 assert sha(s['reference_sinex']['path'])==s['reference_sinex']['sha256']
 selected=[c for c in s['cases'] if not args.case or c['case']==args.case];assert selected
 for case in selected:
  output=Path(case['output']);receipt=WORK/(case['case']+'.result.json')
  if receipt.exists():
   previous=json.loads(receipt.read_text());assert previous['exit_code']==0 and previous['binary_sha256']==s['binary_sha256'] and previous['config_sha256']==case['config_sha256']
   print('already completed',case['case'],flush=True);continue
  assert not output.exists(),f'Refusing overwrite: {output}'
  assert sha(case['config'])==case['config_sha256']
  available=int(next(x.split()[1] for x in Path('/proc/meminfo').read_text().splitlines() if x.startswith('MemAvailable:')))
  assert available>1200000,f'Insufficient spare memory: {available} KiB'
  output.mkdir()
  env=os.environ.copy();env.update(s['threads'])
  command=[s['binary'],'-q','-y',case['config'],'-d',case['case'],'-r',case['observation']]
  result=dict(case,command=command,binary_sha256=s['binary_sha256'],thread_environment=s['threads'],started_at=datetime.datetime.now().isoformat(),runner_commit=subprocess.check_output(['git','-C',s['source_root'],'rev-parse','HEAD'],text=True).strip())
  started=time.monotonic()
  with (output/'run.log').open('w') as f:
   proc=subprocess.Popen(command,cwd=s['working_directory'],env=env,stdout=f,stderr=subprocess.STDOUT,start_new_session=True)
   os.setpriority(os.PRIO_PROCESS,proc.pid,10)
   result['pid']=proc.pid
   (WORK/'current_case.json').write_text(json.dumps(result,indent=2))
   print('started',case['case'],'pid',proc.pid,flush=True)
   try:status=proc.wait(timeout=1200)
   except subprocess.TimeoutExpired:
    result['timed_out']=True;os.killpg(proc.pid,signal.SIGTERM)
    try:status=proc.wait(timeout=15)
    except subprocess.TimeoutExpired:os.killpg(proc.pid,signal.SIGKILL);status=proc.wait()
  result.update(exit_code=status,wall_seconds=time.monotonic()-started,finished_at=datetime.datetime.now().isoformat())
  receipt.write_text(json.dumps(result,indent=2));print('finished',case['case'],status,round(result['wall_seconds'],2),flush=True)
  if status!=0:print('FAILED CASE RETAINED; continuing independent cases',case['case'],flush=True)
 (WORK/'runner_completed.json').write_text(json.dumps({'selection':[c['case'] for c in selected],'time':datetime.datetime.now().isoformat()},indent=2))
if __name__=='__main__':main()
