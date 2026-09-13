from pathlib import Path
import csv,json,re,datetime,statistics,collections,hashlib
w=Path('/mnt/c/Users/rx/Documents/GINAN/r48_user_model_work_20260914');v=w/'revision2';out=v/'analysis/station_ar_product_usage';out.mkdir(exist_ok=True)
spec=json.loads((v/'experiment.json').read_text())
def fields(line):return dict(re.findall(r'([A-Za-z0-9_]+)=([^\s]+)',re.sub(r'time=(\d{4}-\d{2}-\d{2}) (\d{2}:\d{2}:\d{2})',r'time=\1T\2',line)))
def stamp(e):return (datetime.datetime(2024,7,17)+datetime.timedelta(seconds=int(float(e))-1405209600)).isoformat()
times=[stamp(e) for e in range(spec['first_epoch'],spec['last_epoch']+1,30)]
products=collections.defaultdict(dict)
for row in csv.DictReader((w/'products/zhang_internal_products.csv').open()):
 if row['solution']=='PRODUCT_FIXED' and stamp(row['gpst_seconds']) in times:products[stamp(row['gpst_seconds'])][row['satellite'],row['observable']]=row
flags=['numeric_valid','branch_valid','continuity_valid','pppar_usable','integer_valid','integer_structure_valid','integer_datum_continuous','integer_precision_valid','dual_frequency_ar_valid','ar_valid']
eligible={};components={}
for t,p in products.items():
 eligible[t]=set();components[t]={}
 for sat in {s for s,c in p}:
  pair=[p.get((sat,c),{}) for c in ['L1C','L2W']]
  if all(r.get('product_state')=='AR_VALID' and all(r.get(f)=='1' for f in flags) for r in pair) and pair[0]['integer_component_id']==pair[1]['integer_component_id'] and pair[0]['integer_component_id'] not in ['', 'NONE']:
   eligible[t].add(sat);components[t][sat]=pair[0]['integer_component_id']
allrows=[];summary=[];hashes=[]
for case in spec['cases']:
 root=Path(case['output']);receiver=next(root.glob(case['station']+'*.TRACE'));network=next(root.glob('Network*.TRACE'))
 used=collections.defaultdict(lambda:collections.defaultdict(set));user=collections.defaultdict(lambda:collections.defaultdict(set));lattice={}
 for line in receiver.read_text().splitlines():
  if line.startswith('ZHANG_E27_PRODUCT_COVARIANCE'):
   r=fields(line)
   if r.get('valid')=='1' and r.get('measurement')=='PHASE':used[r['time']][r['satellite']].add(r['observable'])
 for line in network.read_text().splitlines():
  if line.startswith('ZHANG_USER_AMBIGUITY '):
   r=fields(line)
   if r.get('integer_valid')=='1':user[r['time']][r['satellite']].update(r['observable'].split('-'))
  if line.startswith(('ZHANG_E26_USER_LATTICE','ZHANG_E27_USER_LATTICE')):
   r=fields(line);lattice[r['time']]=r
 assert sorted(used)==times,case['case']
 rows=[]
 for t in times:
  dual={sat for sat,codes in used[t].items() if {'L1C','L2W'}<=codes};ar=dual&eligible[t];local={sat for sat,codes in user[t].items() if {'L1C','L2W'}<=codes}
  groups=collections.defaultdict(list)
  for sat in sorted(ar):groups[components[t][sat]].append(sat)
  rr={'case':case['case'],'station':case['station'],'model':case['model'],'mode':case['mode'],'time':t,'dual_product_satellites':len(dual),'ar_qualified_product_satellites':len(ar),'ar_satellites':';'.join(sorted(ar)),'ar_components':json.dumps(dict(groups),sort_keys=True),'component_count':len(groups),'maximum_within_component_sd_rank':sum(max(0,len(s)-1) for s in groups.values()),'runtime_integer_valid_satellites_excluding_reference':len(local),'runtime_satellites':';'.join(sorted(local)),'search_common_satellites':int(lattice.get(t,{}).get('common_satellites',0)),'search_named_rank':int(lattice.get(t,{}).get('named_rank',0))}
  rows.append(rr)
 allrows.extend(rows);counts=[r['ar_qualified_product_satellites'] for r in rows];positive=[r for r in rows if r['ar_qualified_product_satellites']]
 summary.append({k:case[k] for k in ['case','station','model','mode']}|{'epochs':len(rows),'dual_products_min':min(r['dual_product_satellites'] for r in rows),'dual_products_max':max(r['dual_product_satellites'] for r in rows),'ar_min':min(counts),'ar_max':max(counts),'ar_mean':statistics.mean(counts),'ar_satellite_epochs':sum(counts),'ar_unique_satellites':len(set(sat for r in rows for sat in r['ar_satellites'].split(';') if sat)),'ar_positive_epochs':len(positive),'first_ar_epoch':positive[0]['time'] if positive else '', 'last_ar_count':counts[-1],'last_ar_satellites':rows[-1]['ar_satellites'],'runtime_ar_mean_excluding_reference':statistics.mean(r['runtime_integer_valid_satellites_excluding_reference'] for r in rows),'runtime_ar_max_excluding_reference':max(r['runtime_integer_valid_satellites_excluding_reference'] for r in rows),'search_common_max':max(r['search_common_satellites'] for r in rows)})
 for p in [receiver,network]:hashes.append({'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
def save(name,rows):
 with (out/name).open('w',newline='') as f:
  wr=csv.DictWriter(f,fieldnames=list(rows[0]));wr.writeheader();wr.writerows(rows)
save('all_cases_per_epoch.csv',allrows);save('all_cases_summary.csv',summary)
base=[r for r in summary if r['mode']=='float'];save('station_summary.csv',base)
cross=[]
for t in times:
 r={'time':t}
 for station in ['MARS','DYNG','NICO','BREW','JPLM']:
  found=[x for x in allrows if x['time']==t and x['station']==station and x['mode']=='float'];assert len(found)==2
  for x in found:r[station+'_'+x['model']]=x['ar_qualified_product_satellites']
 cross.append(r)
save('station_per_epoch.csv',cross)
same=[]
for station in ['MARS','DYNG','NICO','BREW','JPLM']:
 groups=collections.defaultdict(list)
 for x in allrows:
  if x['station']==station:groups[x['case']].append(x['ar_satellites'])
 same.append({'station':station,'same_product_eligibility_across_all_ten_cases':len({tuple(a) for a in groups.values()})==1})
(out/'summary.json').write_text(json.dumps({'definition':'Dual-frequency products applied during phase observation construction, intersected with strict PRODUCT_FIXED AR_VALID product flags and common L1C/L2W component; deduplicated across QC retries. This is product eligibility, not final measurement acceptance or integer fixing.','flags':flags,'summary':base,'cross_case_consistency':same,'hashes':hashes},indent=2)+'\n')
assert all(x['same_product_eligibility_across_all_ten_cases'] for x in same)
compact=[{'time':r['time'],**{station:r[station+'_UDUC'] for station in ['MARS','DYNG','NICO','BREW','JPLM']}} for r in cross]
save('station_per_epoch_counts.csv',compact)
report=['# 各站使用的 R48 AR 合格产品统计','','时段：2024-07-17 00:03:00—00:30:00 GPST，55历元。统计单位为卫星，两频合并计1颗。','', '口径：站端相位观测构建时实际引用 L1C、L2W 产品协方差的卫星，与同历元 PRODUCT_FIXED / AR_VALID 双频产品严格资格取交集，两频分量编号须相同且非空。QC多次构建去重；不是最终观测验收、整数搜索或固定成功数。统计读取冻结产品和已完成的50个用例TRACE，没有重跑解算。','', '五站各自十个IF/UDUC及FLOAT/shadow/feedback用例的产品资格集合完全一致，下面每站只列一次。','', '|站点|使用双频产品卫星范围|AR合格范围|AR平均颗数|有合格产品历元|合格卫星历元累计|累计不同合格卫星|00:30合格颗数|','|---|---:|---:|---:|---:|---:|---:|---:|']
for station in ['MARS','DYNG','NICO','BREW','JPLM']:
 x=next(x for x in base if x['station']==station and x['model']=='UDUC')
 report.append(f"|{station}|{x['dual_products_min']}–{x['dual_products_max']}|{x['ar_min']}–{x['ar_max']}|{x['ar_mean']:.2f}|{x['ar_positive_epochs']}/55|{x['ar_satellite_epochs']}|{x['ar_unique_satellites']}|{x['last_ar_count']}|")
report+=['', '平均值包含没有AR合格产品的历元；累计卫星历元为逐历元合格颗数之和，不是独立卫星数。若按L1C/L2W两条产品记录计数，对应卫星历元数乘2。', '', 'BREW、JPLM首次出现合格产品在00:10:00，其余三站为00:10:30。产品资格仍需与用户参考星、连通分量、连续性、状态存在性及观测筛选相结合。all_cases_per_epoch.csv 另列运行时可固定模糊度卫星数和实际搜索共同卫星数；前者排除接收机参考星，不能与产品资格颗数直接混用。', '', '逐历元五站表：station_per_epoch_counts.csv；含卫星编号、分量及搜索计数的全部50用例明细：all_cases_per_epoch.csv。']
(out/'README.md').write_text('\n'.join(report)+'\n')
print(json.dumps({'summary':base,'cross_case_consistency':same},indent=2))
