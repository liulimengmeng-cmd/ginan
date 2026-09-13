"""One-shot validation of frozen inputs and final R49 outputs. Never polls."""
from pathlib import Path
import argparse,json,csv,hashlib,subprocess,re,collections,itertools
from metrics import quality,metrics
W=Path(__file__).resolve().parent
SPEC=json.loads((W/'experiment.json').read_text())
OUT=Path('/mnt/d/GINAN_R20/inputData/outputs')/SPEC['case']
def sha(path):
 h=hashlib.sha256()
 with path.open('rb') as f:
  for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
 return h.hexdigest()
def preflight():
 out={'files':0,'bytes':0,'failures':[]}
 options=subprocess.check_output(['findmnt','-nro','OPTIONS','/'],text=True).strip()
 if 'rw' not in options.split(','):out['failures'].append('ROOT_READ_ONLY')
 for line in (W/'INPUTS_SHA256SUMS').read_text().splitlines():
  digest,name=line.split('  ',1);p=Path('/mnt/d/GINAN_R20')/name
  if not p.is_file():out['failures'].append('MISSING:'+name);continue
  out['files']+=1;out['bytes']+=p.stat().st_size
  if sha(p)!=digest:out['failures'].append('HASH:'+name)
 if out['files']!=202:out['failures'].append('INPUT_COUNT')
 (W/'r49_input_preflight.json').write_text(json.dumps(out,indent=2)+'\n')
 print(json.dumps(out));return 0 if not out['failures'] else 98
def enrich(epochs):
 previous=set()
 for row in epochs:
  pairs={tuple(pair) for group in row['components'].values() for pair in itertools.combinations(sorted(group),2)}
  row['certified_pair_count']=len(pairs)
  row['pair_loss_vs_previous']=len(previous-pairs)
  row['pair_gain_vs_previous']=len(pairs-previous)
  row['component_count']=len(row['components'])
  previous=pairs
 return epochs
def final(status):
 result={'case':SPEC['case'],'pea_status':status,'failures':[],'alarms':{},'claim_boundary':SPEC['claims']}
 if status!=0:result['failures'].append('PEA_EXIT_NONZERO')
 try:
  for name in ['zhang_internal_products.csv','zhang_internal_product_covariance.csv']:
   with (OUT/name).open(newline='') as f:
    for row in csv.DictReader(f):
     if None in row or any(v is None for v in row.values()):raise RuntimeError('INCOMPLETE_CSV:'+name)
  epochs=enrich(quality(OUT/'zhang_internal_products.csv'))
  result['epochs']=epochs;result['metrics']=metrics(epochs)
  result['metrics']['certified_pair_epoch_integral']=sum(r['certified_pair_count'] for r in epochs)
  if [e['epoch'] for e in epochs]!=list(range(1405209600,1405216801,30)):result['failures'].append('EPOCH_COVERAGE')
  if any(e['duplicate_product_satellites'] for e in epochs):result['failures'].append('DUPLICATE_PRODUCT_ROWS')
  result['checkpoint_bundles']=len(list(Path(SPEC['checkpoints']['directory']).glob('*/checkpoint_manifest.json')))
  if result['checkpoint_bundles']!=241:result['failures'].append('CHECKPOINT_COVERAGE')
  alarms=collections.Counter();events=collections.Counter();timers=collections.defaultdict(float)
  records=[]
  for path in list(OUT.glob('*.TRACE'))+[W/(SPEC['case']+'.runner.log')]:
   if not path.exists():continue
   with path.open(errors='replace') as f:
    for line in f:
     for token in ['R49_INVALID_DGEMV_ARGUMENT','std::bad_alloc','INVALID KALMAN FILTER BLOCK','PRODUCT_PAIR_AFFINE_PULLBACK_MISMATCH','TRUE_PHYSICAL_EXPANSION_FAILED','EXACT_PHYSICAL_ROUNDTRIP_FAILED','R49_SEARCH_FAMILY_BUDGET_EXCEEDED','R49_BLAS_ABI_MISMATCH','R49_CONDITION_RECORD_COVERAGE_MISMATCH']:
      if token in line:alarms[token]+=1
     if 'DGEMV' in line and 'illegal value' in line:alarms['DGEMV_ILLEGAL_VALUE']+=1
     if line.startswith('R49_'):
      tag=line.split()[0];events[tag]+=1;records.append(line.rstrip())
     if line.startswith('ZHANG_AR_PHASE_TIMER'):
      fields=dict(re.findall(r'(\w+)=([^\s]+)',line))
      if 'exclusive_ms' in fields:timers[fields['phase']]+=float(fields['exclusive_ms'])/1000
  result['alarms']=dict(alarms);result['r49_event_counts']=dict(events);result['exclusive_timer_seconds']=dict(timers)
  if alarms:result['failures'].append('CONTRACT_OR_NUMERICAL_ALARM')
  (W/'r49_event_records.txt').write_text('\n'.join(records)+'\n')
  result['matched_comparison']={}
  for key in ['r47','r48fix1']:
   baseline=enrich(quality(Path(SPEC['baseline_outputs'][key])/'zhang_internal_products.csv'))
   b={e['epoch']:e for e in baseline};c={e['epoch']:e for e in epochs};common=sorted(b.keys()&c.keys())
   bm=metrics([b[e] for e in common]);cm=metrics([c[e] for e in common])
   bm['certified_pair_epoch_integral']=sum(b[e]['certified_pair_count'] for e in common)
   cm['certified_pair_epoch_integral']=sum(c[e]['certified_pair_count'] for e in common)
   result['matched_comparison'][key]={'epoch_count':len(common),'baseline':bm,'r49':cm}
  lines=['# R49 completion audit','',f'PEA exit: {status}. Completed epochs: {len(epochs)}/241.','',
    '|Epoch|Strict dual AR satellites|Largest component|Components|Certified satellite pairs|','|---|---:|---:|---:|---:|']
  lines += [f"|{e['time']}|{e['strict_ar']}|{e['largest_component']}|{e['component_count']}|{e['certified_pair_count']}|" for e in epochs]
  lines += ['', 'Failures: '+(', '.join(result['failures']) or 'NONE'),'',SPEC['claims']]
  (W/'R49_result_summary.md').write_text('\n'.join(lines)+'\n')
 except Exception as e:result['failures'].append(type(e).__name__+': '+str(e))
 result['runtime_pass']=not result['failures']
 print(json.dumps(result,ensure_ascii=False,indent=2));return 0 if result['runtime_pass'] else 109
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--preflight',action='store_true');p.add_argument('--pea-status',type=int,default=-1);a=p.parse_args()
 raise SystemExit(preflight() if a.preflight else final(a.pea_status))
