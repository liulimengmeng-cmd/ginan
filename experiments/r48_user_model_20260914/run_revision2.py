from pathlib import Path
import time,subprocess,json,shutil,hashlib,datetime,yaml,os
w=Path('/mnt/c/Users/rx/Documents/GINAN/r48_user_model_work_20260914');r=Path('/home/rx/GINAN/r48-user-model-20260914')
v=w/'revision2';assert not v.exists();v.mkdir();(v/'config').mkdir();(v/'products').symlink_to(w/'products',target_is_directory=True)
shutil.copy2(w/'product_manifest.json',v/'product_manifest.json')
spec=json.loads((w/'experiment.json').read_text());spec['work']=str(v);spec['previous_attempt']=str(w)
for case in spec['cases']:
 old=case['case'];new=old+'_fix1';cfg=yaml.safe_load(Path(case['config']).read_text())
 cfg['outputs']['metadata']['config_description']=new;cfg['processing_options']['gnss_general']['zhang_pppar']['checkpoint_runtime_id']=new
 path=v/'config'/f'{new}.yaml';path.write_text(yaml.safe_dump(cfg,sort_keys=False))
 case.update(case=new,config=str(path),output=case['output']+'_fix1',config_sha256=hashlib.sha256(path.read_bytes()).hexdigest())
(v/'experiment.json').write_text(json.dumps(spec,indent=2)+'\n')
while not (w/'suite_complete.json').exists():time.sleep(5)
print('Initial suite ended; compiling revision2',flush=True)
source=(w/'freeze_user.py').read_text().replace("w=Path('/mnt/c/Users/rx/Documents/GINAN/r48_user_model_work_20260914')", "w=Path('/mnt/c/Users/rx/Documents/GINAN/r48_user_model_work_20260914/revision2')")
exec(compile(source,'freeze_revision2','exec'),{})
print('Revision2 frozen; starting matched suite',flush=True)
rc=subprocess.run(['python3',str(w/'run_suite.py'),str(v)]).returncode
assert rc==0
subprocess.run(['python3',str(w/'analyse_suite.py'),str(v)],check=True)
