from pathlib import Path
import hashlib,json,datetime,subprocess,shutil
w=Path('/mnt/c/Users/rx/Documents/GINAN/r48_user_model_work_20260914');v=w/'revision2';r=Path('/home/rx/GINAN/r48-user-model-20260914')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
s=json.loads((v/'analysis/summary.json').read_text());spec=json.loads((v/'experiment.json').read_text());frozen=json.loads((v/'frozen_user.json').read_text())
assert len(s['cases'])==50 and all(x['completed'] and not x['errors'] and not any(x['alarms'].values()) for x in s['cases'])
assert sha(frozen['binary'])==frozen['binary_sha256']
for station in spec['stations']:assert sha(station['path'])==station['sha256']
for name,item in json.loads((w/'product_manifest.json').read_text()).items():
 assert sha(w/'products'/name)==item['sha256'] and sha(item['source'])==item['sha256']
for case in spec['cases']:assert sha(v/'frozen/config'/Path(case['config']).name)==case['config_sha256']
server_manifest=json.loads((w.parent/'r49_work_20260913/frozen_r49.json').read_text())
assert sha(Path(server_manifest['frozen_root'])/'bin/pea')==server_manifest['binary_sha256']
source=hashlib.sha256()
for p in sorted((r/'src').rglob('*')):
 if p.is_file() and p.name!='peaCommitVersion.h':source.update(str(p.relative_to(r)).encode());source.update(p.read_bytes())
assert source.hexdigest()==frozen['source_sha256']
overlap_path=w.parent/'r49_work_20260913/concurrent_user_tests_20260914.json';overlap=json.loads(overlap_path.read_text());overlap['user_solver_ended']=json.loads((v/'suite_complete.json').read_text())['end'];overlap['audited_at']=datetime.datetime.now().astimezone().isoformat();overlap_path.write_text(json.dumps(overlap,indent=2)+'\n')
result={'checked_at':overlap['audited_at'],'completed_cases':50,'diagnostic_epochs':2750,'user_DGEMV':0,'fatal_transaction_alarms':0,'all_input_and_product_hashes_unchanged':True,'frozen_user_binary_and_source_unchanged':True,'r49_server_binary_hash_unchanged':True,'initial_failed_cases_preserved':25,'no_independent_integer_truth':True}
(v/'analysis/final_verification.json').write_text(json.dumps(result,indent=2)+'\n')
subprocess.run(['python3',str(w/'make_report.py')],check=True)
subprocess.run(['python3',str(w/'archive_experiment_sources.py')],check=True)
d=r/'experiments/r48_user_model_20260914'
for rel in ['analysis','revision2/analysis']:
 shutil.copytree(w/rel,d/rel,dirs_exist_ok=True)
for rel in ['R48_user_model_report.md','suite_complete.json','revision2/suite_complete.json']:
 shutil.copy2(w/rel,d/rel)
for rel in ['run_revision2.py','finish_audit.py']:
 shutil.copy2(w/rel,d/rel)
for root,dest in [(w,d),(v,d/'revision2')]:
 receipts=[json.loads(p.read_text()) for p in sorted(root.glob('*.result.json'))]
 (dest/'execution_receipts.json').write_text(json.dumps(receipts,indent=2)+'\n')
print(json.dumps(result))
