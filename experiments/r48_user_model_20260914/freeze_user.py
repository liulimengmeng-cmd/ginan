from pathlib import Path
import subprocess,hashlib,shutil,json,datetime,time
while Path('/mnt/c/Users/rx/Documents/GINAN/r48_user_model_work_20260914/hold_build').exists():time.sleep(5)
w=Path('/mnt/c/Users/rx/Documents/GINAN/r48_user_model_work_20260914');r=Path('/home/rx/GINAN/r48-user-model-20260914')
def source_hash():
 h=hashlib.sha256()
 for p in sorted((r/'src').rglob('*')):
  if p.is_file() and p.name!='peaCommitVersion.h':h.update(str(p.relative_to(r)).encode());h.update(p.read_bytes())
 return h.hexdigest()
before=source_hash()
with (w/'build_final.log').open('w') as f:rc=subprocess.run(['nice','-n','19','cmake','--build',str(r/'build'),'--target','pea','-j1'],stdout=f,stderr=subprocess.STDOUT).returncode
assert rc==0 and before==source_hash()
frozen=w/'frozen';assert not frozen.exists();frozen.mkdir();shutil.copy2(r/'bin/pea',frozen/'pea');shutil.copytree(w/'config',frozen/'config')
sha=hashlib.sha256((frozen/'pea').read_bytes()).hexdigest();commit=subprocess.check_output(['git','-C',str(r),'rev-parse','HEAD'],text=True).strip()
result={'source_commit':commit,'source_sha256':before,'binary_sha256':sha,'binary':str(frozen/'pea'),'build_exit':rc,'build_log':str(w/'build_final.log'),'created':datetime.datetime.now().astimezone().isoformat()}
(w/'frozen_user.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
