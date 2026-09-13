from pathlib import Path
import json,hashlib,subprocess,shutil,re,datetime
w=Path(__file__).resolve().parent
spec=json.loads((w/'experiment.json').read_text());source=Path(spec['source_root']);frozen=Path(spec['frozen_root'])
assert not frozen.exists(),'Never overwrite a frozen run'
current=json.loads((w/'build_current.json').read_text());assert current['up_to_date']
tests=(w/'test_results_final.log').read_text();match=re.search(r'(\d+) test cases out of (\d+) passed',tests)
assert match and match[1]==match[2] and int(match[1])>=433 and 'has failed' not in tests
assert '9 test cases out of 9 passed' in (w/'test_checkpoint_final.log').read_text()
assert json.loads((w/'blas_injection_result.json').read_text())['audit_intercepted']
smoke=json.loads((w/'smoke_result.json').read_text());assert smoke['status']==0
assert 'R49_INVALID_DGEMV_ARGUMENT' not in (w/'smoke.log').read_text()
assert 'illegal value' not in (w/'smoke.log').read_text()
manifests=list(Path(smoke['checkpoint_directory']).glob('*/checkpoint_manifest.json'));assert len(manifests)==1
manifest=json.loads(manifests[0].read_text());assert manifest['runtime_id']=='R49-NETWORK' and manifest['state_dimension']>0
for label in ['selection','posterior','domain']:
 assert 'PASS' in (w/f'r49_{label}_reference.log').read_text()
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def source_hash():
 h=hashlib.sha256()
 for p in sorted((source/'src').rglob('*')):
  if p.is_file() and p.name!='peaCommitVersion.h':h.update(str(p.relative_to(source)).encode());h.update(p.read_bytes())
 return h.hexdigest()
assert source_hash()==current['source_hash_after']
git=lambda *args:subprocess.check_output(['git','-C',str(source),*args],text=True)
assert set(git('diff','--name-only','HEAD').splitlines())<= {'src/cpp/pea/peaCommitVersion.h'}
commit=git('rev-parse','HEAD').strip();assert git('rev-parse','refs/remotes/fork/codex/r49-condition-domain').strip()==commit
frozen.mkdir();shutil.copytree(source/'src',frozen/'source/src');shutil.copytree(w/'config',frozen/'config')
(frozen/'bin').mkdir();(frozen/'evidence').mkdir()
for n in ['pea','zhang_full_rank_tests','zhang_checkpoint_infra0_tests']:shutil.copy2(source/'bin'/n,frozen/'bin'/n)
for n in ['r49_selection_reference','r49_posterior_reference','r49_domain_reference','r49_selection_benchmark']:
 shutil.copy2(source/n,frozen/'bin'/n)
for p in w.iterdir():
 if p.is_file() and p.name not in ['source_conversation.md'] and p.suffix in ['.py','.json','.log','.md','.sh','.txt']:
  shutil.copy2(p,frozen/'evidence'/p.name)
shutil.copy2(w/'INPUTS_SHA256SUMS',frozen/'evidence/INPUTS_SHA256SUMS')
shutil.copy2(source/'R49_IMPLEMENTATION.md',frozen/'evidence/R49_IMPLEMENTATION.md')
shutil.copy2(manifests[0],frozen/'evidence/smoke_checkpoint_manifest.json')
(frozen/'evidence/source_diff_from_r48fix1.patch').write_text(git('diff','820c2cb','--','src'))
runtime={'source_commit':commit,'uname':subprocess.check_output(['uname','-a'],text=True).strip(),'ldd':subprocess.check_output(['ldd',str(frozen/'bin/pea')],text=True),'generated_version':(source/'src/cpp/pea/peaCommitVersion.h').read_text(),'policy':spec['r49_policy']}
(frozen/'evidence/runtime.json').write_text(json.dumps(runtime,indent=2)+'\n')
entries=[f'{sha(p)}  {p.relative_to(frozen).as_posix()}' for p in sorted(frozen.rglob('*')) if p.is_file()]
(frozen/'SHA256SUMS').write_text('\n'.join(entries)+'\n')
result={'frozen_root':str(frozen),'source_commit':commit,'binary_sha256':sha(frozen/'bin/pea'),'tests_passed':True,'test_count':int(match[1]),'source_matches_build':True,'manifest_files':len(entries),'input_data_copied':False,'frozen_at':datetime.datetime.now().astimezone().isoformat()}
(w/'frozen_r49.json').write_text(json.dumps(result,indent=2)+'\n')
for p in frozen.rglob('*'):
 if p.is_file():p.chmod(p.stat().st_mode & ~0o222)
print(json.dumps(result,indent=2))
