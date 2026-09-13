from pathlib import Path
import hashlib,json,subprocess,time
w=Path('/mnt/c/Users/rx/Documents/GINAN/r49_work_20260913')
root=Path('/home/rx/GINAN/r49-condition-domain-20260913')
def source_hash():
    h=hashlib.sha256()
    for p in sorted((root/'src').rglob('*')):
        if not p.is_file() or p.name=='peaCommitVersion.h':continue
        h.update(str(p.relative_to(root)).encode());h.update(p.read_bytes())
    return h.hexdigest()
before=source_hash();start=time.time()
with (w/'build_final_consistency.log').open('w') as log:
    rc=subprocess.run(['cmake','--build',str(root/'build'),'--target','pea','zhang_full_rank_tests','zhang_checkpoint_infra0_tests','-j2'],stdout=log,stderr=subprocess.STDOUT).returncode
after=source_hash()
result={'up_to_date':rc==0 and before==after,'build_exit':rc,'source_hash_before':before,'source_hash_after':after,
        'generated_version_sha256':hashlib.sha256((root/'src/cpp/pea/peaCommitVersion.h').read_bytes()).hexdigest(),'elapsed_seconds':time.time()-start}
(w/'build_current.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result));assert result['up_to_date']
with (w/'test_results_final.log').open('w') as log:
    rc=subprocess.run([str(root/'bin/zhang_full_rank_tests'),'--report_level=detailed'],stdout=log,stderr=subprocess.STDOUT).returncode
print('final_test_exit='+str(rc));assert rc==0

with (w/'test_checkpoint_final.log').open('w') as log:
    rc=subprocess.run([str(root/'bin/zhang_checkpoint_infra0_tests'),'--report_level=detailed'],stdout=log,stderr=subprocess.STDOUT).returncode
assert rc==0
