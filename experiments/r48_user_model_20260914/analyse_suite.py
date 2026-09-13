from pathlib import Path
import csv,datetime,hashlib,json,math,re,statistics,sys
from collections import defaultdict,Counter
w=Path(sys.argv[1] if len(sys.argv)>1 else '/mnt/c/Users/rx/Documents/GINAN/r48_user_model_work_20260914')
spec=json.loads((w/'experiment.json').read_text());out=w/'analysis';out.mkdir(exist_ok=True)
def fields(line):return dict(re.findall(r'([A-Za-z0-9_]+)=([^\s]+)',re.sub(r'time=(\d{4}-\d{2}-\d{2}) (\d{2}:\d{2}:\d{2})',r'time=\1T\2',line)))
def rms(v):return math.sqrt(statistics.fmean(x*x for x in v)) if v else None
def csvwrite(name,rows):
 if not rows:return
 keys=list(dict.fromkeys(k for r in rows for k in r))
 with (out/name).open('w',newline='') as f:
  wr=csv.DictWriter(f,fieldnames=keys);wr.writeheader();wr.writerows(rows)
expected=[(datetime.datetime(2024,7,17)+datetime.timedelta(seconds=e-1405209600)).isoformat() for e in range(spec['first_epoch'],spec['last_epoch']+1,30)]
summaries=[];epochs=[];stages=[];targets=[];relations=[];positions={};hashes=[];commit_rows=[]
for case in spec['cases']:
 receipt=w/(case['case']+'.result.json')
 if not receipt.exists():continue
 receipt=json.loads(receipt.read_text());root=Path(case['output']);trace=next(root.glob('Network-*.TRACE'));txt=trace.read_text();log=(root/'run.log').read_text();rec=defaultdict(list)
 relation_context={}
 for line in txt.splitlines():
  if line.startswith('ZHANG_'):
   tag=line.split()[0];parsed=fields(line)
   if tag=='ZHANG_USER_ILS_RELATIONS':relation_context={k:parsed.get(k,'') for k in ['time','receiver','stage']}
   if tag=='ZHANG_USER_ILS_ROW':parsed={**relation_context,**parsed,'coefficients':line.split('coefficients=',1)[1].strip()}
   rec[tag].append(parsed)
 ident={k:case[k] for k in ['case','station','model','mode']}
 diag=rec['ZHANG_USER_DIAGNOSTIC'];times=[r['time'] for r in diag];assert times==expected[:len(times)],case['case']
 processing=re.findall(r'Time  (2024-07-17 \d\d:\d\d:\d\d)',txt)
 assert [t.replace(' ','T') for t in processing]==expected[:len(processing)]
 if receipt['exit_code']==0:assert times==expected and len(processing)==len(expected),case['case']
 ar={r['time']:r for r in rec['ZHANG_USER_AR_SUMMARY']};st=defaultdict(dict)
 for r in rec['ZHANG_E27_USER_STAGE']+rec['ZHANG_E26_USER_STAGE']:
  passed=r.get('reliability_gate')=='1' and int(r.get('selected',0))>0 and float(r.get('maximum_perr','inf'))<=.001 and float(r.get('joint_nis','nan'))<=float(r.get('joint_threshold','nan'))
  st[r['time']][r['stage']]={**r,'passed':passed};stages.append({**ident,**r,'passed':passed})
 for r in rec['ZHANG_E27_USER_TARGET']+rec['ZHANG_E26_USER_TARGET']:targets.append({**ident,**r})
 for tag in ['ZHANG_USER_ILS_RELATIONS','ZHANG_USER_ILS_ROW','ZHANG_USER_ILS_NAMED_RECOVERY']:
  for r in rec[tag]:relations.append({**ident,'tag':tag,**r})
 pos=next(root.glob('*.POS'));ref=re.search(r'XYZ Reference position\s*:\s*([-\d.e+]+)\s+([-\d.e+]+)\s+([-\d.e+]+)',pos.read_text());assert ref
 xyz=list(map(float,ref.groups()));snx=next(r['coordinates_m'] for r in spec['reference_sinex']['reference'] if r['station']==case['station'])
 refdiff=math.sqrt(sum((a-snx[b])**2 for a,b in zip(xyz,['STAX','STAY','STAZ'])));assert refdiff<.01
 pp={};wlpassed=0;l1passed=0
 prefix='USER_IF_' if case['model']=='IF' else 'USER_'
 for r in diag:
  t=r['time'];enu=[float(r[k]) for k in ['east_error_m','north_error_m','up_error_m']];assert all(math.isfinite(x) for x in enu);pp[t]=enu
  wl=st[t].get(prefix+'WL_SD',{});l1=st[t].get(prefix+'L1_SD_CONDITIONAL',{});a=ar.get(t,{})
  wlpassed+=bool(wl.get('passed'));l1passed+=bool(l1.get('passed'))
  epochs.append({**ident,'time':t,'east_m':enu[0],'north_m':enu[1],'up_m':enu[2],'error_3d_m':math.sqrt(sum(x*x for x in enu)),'ambiguities':r.get('ambiguities'),'integer_valid_ambiguities':r.get('integer_valid_ambiguities'),'wl_selected':wl.get('selected',0),'wl_passed':bool(wl.get('passed')),'l1_selected':l1.get('selected',0),'l1_passed':bool(l1.get('passed')),'newly_fixed':int(a.get('newly_fixed',0)),'held_integer_rank':int(a.get('held_integer_rank',0))})
 positions[(case['station'],case['model'],case['mode'])]=pp
 newly=sum(int(r.get('newly_fixed',0)) for r in ar.values());held=max([int(r.get('held_integer_rank',0)) for r in ar.values()]+[0])
 if case['mode']=='float':assert not ar and newly==0
 effective=sorted({r.get('mode') for r in rec['ZHANG_USER_EFFECTIVE_SEARCH']})
 if effective:assert effective==(['LAMBDA_ALT'] if case['joint_ils'] else ['ROUND']),effective
 alarms={k:(log+'\n'+txt).count(k) for k in ['DGEMV','AUTHORITATIVE MEASUREMENT TRANSACTION REJECTED','AUTHORITATIVE PREDICTION TRANSACTION REJECTED']}
 errors=[l[l.index('Error:'):] for l in log.splitlines() if 'Error:' in l]
 summaries.append({**ident,'completed':receipt['exit_code']==0 and times==expected,'epochs':len(times),'exit_code':receipt['exit_code'],'wall_seconds':receipt['wall_seconds'],'ar_attempt_epochs':len(ar),'wl_passed_epochs':wlpassed,'l1_passed_epochs':l1passed,'reported_newly_fixed_sum':newly,'held_integer_rank_peak':held,'newly_fixed_epochs':sum(int(r.get('newly_fixed',0))>0 for r in ar.values()),'rms_3d_m':rms([math.sqrt(sum(x*x for x in p)) for p in pp.values()]),'final_enu_m':list(pp.values())[-1] if pp else None,'effective_search':effective,'general_rows_total':sum(int(r.get('general_rows',0)) for r in rec['ZHANG_USER_ILS_RELATIONS']),'named_recovered_total':sum(int(r.get('named_recovered',0)) for r in rec['ZHANG_USER_ILS_NAMED_RECOVERY']),'reference_difference_m':refdiff,'errors':errors,'alarms':alarms,'trace':str(trace)})
 commits=rec['ZHANG_USER_INTEGER_FACTOR_COMMIT']
 commit_rows.extend({**ident,**r} for r in commits)
 summaries[-1].update(integer_factor_commits=len(commits),integer_factor_rows=sum(int(r['rows']) for r in commits),integer_factor_l1_rows=sum(int(r['rows']) for r in commits if 'L1' in r.get('provenance','')),held_positive_epochs=sum(int(r.get('held_integer_rank',0))>0 for r in ar.values()))
 for f in [trace,pos,root/'run.log']:hashes.append({'path':str(f),'bytes':f.stat().st_size,'sha256':hashlib.sha256(f.read_bytes()).hexdigest()})
pairs=[]
for s in summaries:
 if s['mode']=='float':continue
 base=positions.get((s['station'],s['model'],'float'));p=positions[(s['station'],s['model'],s['mode'])]
 if not base:continue
 common=sorted(set(base)&set(p));diffs=[math.sqrt(sum((a-b)**2 for a,b in zip(base[t],p[t]))) for t in common]
 if s['completed'] and 'shadow' in s['mode']:assert max(diffs,default=0)==0,(s['case'],max(diffs))
 pairs.append({k:s[k] for k in ['station','model','mode','completed','epochs','wl_passed_epochs','l1_passed_epochs','reported_newly_fixed_sum','held_integer_rank_peak','general_rows_total','named_recovered_total','integer_factor_commits','integer_factor_rows','integer_factor_l1_rows']}|{'matched_epochs':len(common),'float_matched_3d_rms_m':rms([math.sqrt(sum(x*x for x in base[t])) for t in common]),'test_matched_3d_rms_m':rms([math.sqrt(sum(x*x for x in p[t])) for t in common]),'max_position_difference_m_trace_precision':max(diffs,default=None)})
result={'expected_epochs':len(expected),'cases':summaries,'paired':pairs,'independent_integer_truth_available':False,'product_temporal_covariance_available':False,'all_cases_finished':len(summaries)==len(spec['cases'])}
(out/'summary.json').write_text(json.dumps(result,indent=2)+'\n');(out/'output_sha256.json').write_text(json.dumps(hashes,indent=2)+'\n')
for name,rows in [('case_summary.csv',summaries),('paired_summary.csv',pairs),('per_epoch.csv',epochs),('stage_gates.csv',stages),('named_targets.csv',targets),('ils_relations.csv',relations),('integer_commits.csv',commit_rows)]:csvwrite(name,rows)
print(json.dumps({'cases':len(summaries),'completed':sum(s['completed'] for s in summaries),'reported_newly_fixed_sum':sum(s['reported_newly_fixed_sum'] for s in summaries),'failed':[{k:s[k] for k in ['station','model','mode','epochs','errors']} for s in summaries if not s['completed']]},indent=2))
