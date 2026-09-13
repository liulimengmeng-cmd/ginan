#!/usr/bin/env python3
"""Reproduce first UDUC failure with info logs while retaining all numerical gates."""
import hashlib,json,os,subprocess,time
from pathlib import Path
import yaml
W=Path('/mnt/c/Users/rx/Documents/GINAN/r48_uduc_work_20260913')
s=json.loads((W/'experiment.json').read_text());c=next(c for c in s['cases'] if c['station']=='MARS' and c['mode']=='float')
name='r48_uduc_20260913_MARS_float_diagnostic';out=Path(s['working_directory'])/'outputs'/name;out.mkdir(exist_ok=False)
cfg=yaml.safe_load(Path(c['config']).read_text());cfg['outputs']['metadata']['config_description']=name;cfg['processing_options']['epoch_control']['end_epoch']='2024-07-17 00:03:30';cfg['processing_options']['gnss_general']['zhang_pppar']['checkpoint_runtime_id']=name
cp=W/'config'/(name+'.yaml');assert not cp.exists();cp.write_text(yaml.safe_dump(cfg,sort_keys=False))
env=os.environ.copy();env.update(s['threads']);command=[s['binary'],'-y',str(cp),'-d',name,'-r',c['observation']]
start=time.monotonic()
with (out/'run.log').open('w') as f:
 proc=subprocess.Popen(command,cwd=s['working_directory'],env=env,stdout=f,stderr=subprocess.STDOUT);os.setpriority(os.PRIO_PROCESS,proc.pid,10);status=proc.wait(timeout=120)
r={'case':name,'exit_code':status,'wall_seconds':time.monotonic()-start,'command':command,'binary_sha256':s['binary_sha256'],'config_sha256':hashlib.sha256(cp.read_bytes()).hexdigest(),'log':str(out/'run.log')};(W/'diagnostic.json').write_text(json.dumps(r,indent=2));print(json.dumps(r,indent=2))
