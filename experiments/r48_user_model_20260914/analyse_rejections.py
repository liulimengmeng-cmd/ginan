from pathlib import Path
import json,re,csv,math,collections
w=Path('/mnt/c/Users/rx/Documents/GINAN/r48_user_model_work_20260914/revision2');a=w/'analysis'
spec=json.loads((w/'experiment.json').read_text());rows=[];search_rows=[]
def fields(line):return dict(re.findall(r'([A-Za-z0-9_]+)=([^\s]+)',re.sub(r'time=(\d{4}-\d{2}-\d{2}) (\d{2}:\d{2}:\d{2})',r'time=\1T\2',line)))
for case in spec['cases']:
 if not (w/(case['case']+'.result.json')).exists():continue
 trace=next(Path(case['output']).glob('Network*.TRACE'));groups=collections.defaultdict(list);targets=collections.defaultdict(list)
 for line in trace.read_text().splitlines():
  if line.startswith('R49_SEARCH_DIAGNOSTIC'):
   r=fields(line);groups[r['label']].append(r);search_rows.append({'case':case['case'],'station':case['station'],'model':case['model'],'mode':case['mode'],**r})
  if line.startswith(('ZHANG_E26_USER_TARGET','ZHANG_E27_USER_TARGET')):
   r=fields(line);targets[r['stage']].append(r)
 for stage,items in groups.items():
  perr=[float(x['perr']) for x in targets[stage] if math.isfinite(float(x['perr']))];sigmas=[math.sqrt(float(x['variance'])) for x in targets[stage] if float(x['variance'])>=0]
  rows.append({'station':case['station'],'model':case['model'],'mode':case['mode'],'stage':stage,'attempts':len(items),'search_exit_reasons':dict(collections.Counter(x['search_exit_reason'] for x in items)),'ratio_executed':sum(x['ratio_test_executed']=='1' for x in items),'maximum_best_1d_bootstrap':max(float(x['bootstrap_best_1d']) for x in items),'minimum_named_perr':min(perr,default=None),'minimum_named_sigma_cycles':min(sigmas,default=None)})
(a/'rejection_summary.json').write_text(json.dumps(rows,indent=2)+'\n')
if search_rows:
 with (a/'search_diagnostics.csv').open('w',newline='') as f:
  wr=csv.DictWriter(f,fieldnames=list(dict.fromkeys(k for r in search_rows for k in r)));wr.writeheader();wr.writerows(search_rows)
print(json.dumps([r for r in rows if r['mode']=='ils_feedback'],indent=2))
