"""Count committed network WL/L1 constraints, separate from product publication."""
from pathlib import Path
import json,re,collections
W=Path(__file__).resolve().parent
spec=json.loads((W/'experiment.json').read_text())
root=Path('/mnt/d/GINAN_R20/inputData/outputs')/spec['case']
epochs=collections.defaultdict(lambda: {'layers':{},'commits':{},'canonical_counts':collections.Counter(),'sources':{}})
tags={'ZHANG_AR_SUMMARY','ZHANG_LAYERED_AR_RESULT','ZHANG_FIXED_SUBTRANSACTION','ZHANG_R47_WRITER_ATOMIC_COMMIT','ZHANG_CANONICAL_INTEGER_COORDINATE'}
for path in root.glob('Network*.TRACE'):
 with path.open(errors='replace') as stream:
  for line_no,line in enumerate(stream,1):
   tag=line.split(' ',1)[0]
   if tag not in tags:continue
   match=re.search(r'time=2024-07-17 (\d\d:\d\d:\d\d)',line)
   if not match:continue
   fields=dict(re.findall(r'\b(\w+)=([^\s]+)',line));e=epochs[match[1]]
   if tag=='ZHANG_CANONICAL_INTEGER_COORDINATE':
    e['canonical_counts'][fields.get('type','UNKNOWN')]+=1
   elif tag=='ZHANG_AR_SUMMARY':
    assert 'summary' not in e,match[1]
    e['summary']=fields;e['sources']['summary']=line_no
   elif tag=='ZHANG_LAYERED_AR_RESULT' and 'candidates' in fields:
    stage=fields['stage'];assert stage not in e['layers'],(match[1],stage)
    e['layers'][stage]=fields;e['sources'][stage]=line_no
   elif tag=='ZHANG_FIXED_SUBTRANSACTION' and fields['status']=='COMMITTED':
    stage=fields['stage'];assert stage not in e['commits'],(match[1],stage)
    e['commits'][stage]=int(fields['committed_rows'])
   elif tag=='ZHANG_R47_WRITER_ATOMIC_COMMIT':e['product_committed']=fields['committed']=='1'
rows=[]
for time,e in sorted(epochs.items()):
 if 'summary' not in e:
  rows.append({'gpst':time,'missing_summary':True,'candidates':None,'canonical_counts':dict(e['canonical_counts'])})
  continue
 s=e['summary'];wl=e['layers']['WL'];l1=e['layers']['L1']
 n=int(s['candidates']);nw=int(wl['candidates']);nl=int(l1['candidates'])
 fw=e['commits'].get('LAYERED_WIDE_LANE',0);fl=e['commits'].get('LAYERED_FIRST_SIGNAL',0)
 assert fw==int(wl['fixed']) and fl==int(l1['fixed']),(time,'stage search/commit mismatch')
 assert fw+fl==int(s['newly_fixed']),(time,'combined fixed mismatch')
 assert 0<=fw<=nw and 0<=fl<=nl and fw+fl<=n,time
 rows.append({'gpst':time,'candidates':n,'fixed_constraints':fw+fl,'total_rate_pct':100*(fw+fl)/n,
 'wl_candidates':nw,'wl_fixed':fw,'wl_rate_pct':100*fw/nw if nw else None,
 'l1_candidates':nl,'l1_fixed':fl,'l1_rate_pct':100*fl/nl if nl else None,
 'product_committed':e.get('product_committed',False),'canonical_counts':dict(e['canonical_counts']),
 'held_rank_diagnostic':int(s['held_integer_rank']),'source_lines':e['sources']})
assert len(rows)==61,len(rows)
result={'definition':'Committed network WL plus conditional L1 integer constraints divided by AR_SUMMARY candidates. Constraint-dimension coverage, not correctness rate, raw per-link fully-fixed fraction, or product availability.',
 'l1_no_search':'NA when candidates=0; no L1 search was performed.',
 'rows':rows,'checks':'61 epoch entries; where AR summary exists, stage fixed counts equal committed rows and WL+L1 equals summary newly_fixed. Missing diagnostics remain NA.'}
(W/'R49_cycle_fix_rates.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
lines=['# R49逐历元循环模糊度固定比例','','总比例＝（WL已提交固定约束数＋WL条件下L1已提交固定约束数）／AR_SUMMARY候选模糊度维数。它是网络整数约束维数覆盖率，不是整数正确率，也不是卫星产品合格率。固定的是整数线性组合，不能按此推算每条原始链路的L1/L2均已完全固定。L1候选数为0时标记未搜索，不能用0/0计算。','','|GPST|总固定数/候选维数|总比例|WL固定/候选|WL比例|L1固定/候选|L1比例|产品批次提交|','|---|---:|---:|---:|---:|---:|---:|---|']
for r in rows:
 if r.get('missing_summary'):
  lines.append(f"|{r['gpst']}|NA|NA|NA|NA|NA|NA|无AR汇总记录|")
  continue
 l1pct=f"{r['l1_rate_pct']:.2f}%" if r['l1_rate_pct'] is not None else '未搜索'
 lines.append(f"|{r['gpst']}|{r['fixed_constraints']}/{r['candidates']}|{r['total_rate_pct']:.2f}%|{r['wl_fixed']}/{r['wl_candidates']}|{r['wl_rate_pct']:.2f}%|{r['l1_fixed']}/{r['l1_candidates']}|{l1pct}|{'是' if r['product_committed'] else '否'}|")
lines+=['','有AR汇总记录的历元均已交叉核对阶段固定数、子事务COMMITTED行数与AR_SUMMARY.newly_fixed，一致。缺记录历元记NA，不补零。held_integer_rank只保留为协方差诊断交叉参考，不作为固定数的独立依据。','数据源：'+str(next(root.glob('Network*.TRACE'))),'代码口径：冻结ppp_ambres.cpp的ZHANG_LAYERED_AR_RESULT、ZHANG_FIXED_SUBTRANSACTION及traceZhangAmbiguityAndFixedProducts。']
(W/'R49_cycle_fix_rates.md').write_text('\n'.join(lines)+'\n')
print('\n'.join(lines))
print('canonical_samples',json.dumps([{k:r[k] for k in ['gpst','candidates','canonical_counts']} for r in [rows[0],rows[11],rows[-1]]]))
