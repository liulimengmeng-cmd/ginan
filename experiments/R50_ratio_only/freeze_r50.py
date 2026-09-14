from pathlib import Path
import subprocess,json,shutil,hashlib,datetime,re
w=Path(__file__).resolve().parent;spec=json.loads((w/'experiment.json').read_text());r=Path(spec['source_root']);f=Path(spec['frozen_root']);assert not f.exists()
test=json.loads((w/'test_results.json').read_text());assert all(v['status']==0 for v in test.values())
smoke=json.loads((w/'smoke_result.json').read_text());assert smoke['status']==0
log=(w/'smoke.log').read_text();assert 'illegal value' not in log and 'R49_INVALID_DGEMV_ARGUMENT' not in log
trace='\n'.join(p.read_text(errors='replace') for p in Path(smoke['output']).glob('*.TRACE'))
assert 'R50_VALIDATION_POLICY ratio_only=1 ratio_threshold=3' in trace
assert 'R50_RATIO_TEST' in trace
assert not re.search(r'R49_SEARCH_DIAGNOSTIC[^\n]*mode=LAMBDA[^\n]*ratio_test_executed=0[^\n]*reason=ACCEPTED',trace)
assert all(x['baseline_identical'] for x in json.loads((w/'protected_algorithms.json').read_text()))
assert not json.loads((w/'r50_input_preflight.json').read_text())['failures']
git=lambda *a:subprocess.check_output(['git','-C',str(r),*a],text=True)
assert set(git('diff','--name-only','HEAD').splitlines())<={'src/cpp/pea/peaCommitVersion.h'}
records_commit=git('rev-parse','HEAD').strip();assert git('rev-parse','refs/remotes/fork/codex/r50-ratio-only').strip()==records_commit
version=(r/'src/cpp/pea/peaCommitVersion.h').read_text()
embedded=re.search(r'GINAN_COMMIT_HASH\s+"([a-f0-9]+)"',version).group(1)
commit=git('rev-parse',embedded).strip()
assert set(git('diff','--name-only',commit,'HEAD','--','src').splitlines())<= {'src/cpp/pea/peaCommitVersion.h'}
current=json.loads((w/'build_current.json').read_text());assert current['up_to_date']
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def source_hash():
 h=hashlib.sha256()
 for p in sorted((r/'src').rglob('*')):
  if p.is_file() and p.name!='peaCommitVersion.h':h.update(str(p.relative_to(r)).encode());h.update(p.read_bytes())
 return h.hexdigest()
assert source_hash()==current['source_hash_after']
f.mkdir();shutil.copytree(r/'src',f/'source/src');shutil.copytree(w/'config',f/'config');(f/'bin').mkdir();(f/'evidence').mkdir()
for n in ['pea','zhang_full_rank_tests','zhang_checkpoint_infra0_tests']:shutil.copy2(r/'bin'/n,f/'bin'/n)
for p in w.iterdir():
 if p.is_file() and p.suffix in ['.py','.json','.log','.md','.sh','.txt']:shutil.copy2(p,f/'evidence'/p.name)
shutil.copy2(w/'INPUTS_SHA256SUMS',f/'evidence/INPUTS_SHA256SUMS');shutil.copy2(r/'R50_IMPLEMENTATION.md',f/'evidence/R50_IMPLEMENTATION.md')
(f/'evidence/source_diff_from_r49.patch').write_text(git('diff',spec['baseline_snapshot_commit'],'--','src'))
(f/'evidence/runtime.json').write_text(json.dumps({'source_commit':commit,'records_commit':records_commit,'ldd':subprocess.check_output(['ldd',str(f/'bin/pea')],text=True),'policy':spec['r50_policy']},indent=2)+'\n')
(f/'SHA256SUMS').write_text('\n'.join(f'{sha(p)}  {p.relative_to(f).as_posix()}' for p in sorted(f.rglob('*')) if p.is_file())+'\n')
result={'frozen_root':str(f),'source_commit':commit,'records_commit':records_commit,'binary_sha256':sha(f/'bin/pea'),'tests_passed':True,'source_matches_build':True,'input_data_copied':False,'frozen_at':datetime.datetime.now().astimezone().isoformat()}
(w/'frozen_r50.json').write_text(json.dumps(result,indent=2)+'\n')
for p in f.rglob('*'):
 if p.is_file():p.chmod(p.stat().st_mode & ~0o222)
print(json.dumps(result,indent=2))
