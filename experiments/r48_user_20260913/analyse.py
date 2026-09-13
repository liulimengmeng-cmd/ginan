#!/usr/bin/env python3
"""Audit actual user constraints and matched positioning errors; no truth proxy."""
import csv,datetime,hashlib,json,math,re,statistics
from collections import Counter,defaultdict
from pathlib import Path
W=Path('/mnt/c/Users/rx/Documents/GINAN/r48_user_work_20260913_prefix')
KV=re.compile(r'([A-Za-z0-9_]+)=([^\s]+)')
def fields(line):
 return dict(KV.findall(re.sub(r'time=(\d{4}-\d{2}-\d{2}) (\d{2}:\d{2}:\d{2})',r'time=\1T\2',line)))
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def rms(v):return math.sqrt(statistics.fmean(x*x for x in v)) if v else None
def writecsv(p,rows):
 with p.open('w',newline='') as f:
  wr=csv.DictWriter(f,fieldnames=list(rows[0]));wr.writeheader();wr.writerows(rows)
def main():
 s=json.loads((W/'experiment_if.json').read_text());out=W/'analysis';out.mkdir(exist_ok=True)
 expected=[(datetime.datetime(2024,7,17)+datetime.timedelta(seconds=e-1405209600)).isoformat() for e in range(s['first_epoch'],s['last_epoch']+1,30)]
 summaries=[];epoch_rows=[];stage_rows=[];target_rows=[];allcases={};artifact_hashes=[]
 for case in s['cases']:
  receipt=json.loads((W/(case['case']+'.result.json')).read_text())
  p=Path(case['output']);trace=next(p.glob('Network-*.TRACE'));lines=trace.read_text().splitlines();log=(p/'run.log').read_text()
  rec=defaultdict(list)
  for l in lines:
   if l.startswith('ZHANG_'):rec[l.split()[0]].append(fields(l))
  diag=rec['ZHANG_USER_DIAGNOSTIC'];times=[r['time'] for r in diag];assert times==expected[:len(times)],(case['case'],len(times))
  if receipt['exit_code']==0:assert times==expected
  processing=[re.search(r'Time  (2024-07-17 \d\d:\d\d:\d\d)',l).group(1).replace(' ','T') for l in lines if re.search(r'Time  (2024-07-17 \d\d:\d\d:\d\d)',l)]
  assert processing==expected[:len(processing)]
  if receipt['exit_code']==0:assert processing==expected
  ar={r['time']:r for r in rec['ZHANG_USER_AR_SUMMARY']};stage=defaultdict(dict);lattice={r['time']:r for r in rec['ZHANG_E27_USER_LATTICE']}
  for r in rec['ZHANG_E27_USER_STAGE']:
   passed=r.get('reliability_gate')=='1' and int(r.get('selected','0'))>0 and float(r.get('maximum_perr','inf'))<=1e-3 and float(r.get('joint_nis','nan'))<=float(r.get('joint_threshold','nan'))
   stage[r['time']][r['stage']]=dict(r,passed=passed)
   stage_rows.append({'station':case['station'],'mode':case['mode'],'time':r['time'],'stage':r['stage'],'selected':r.get('selected'),'provisional_selected':r.get('provisional_selected'),'joint_nis':r.get('joint_nis'),'joint_threshold':r.get('joint_threshold'),'maximum_perr':r.get('maximum_perr'),'reliability_gate':r.get('reliability_gate'),'passed':passed,'feedback_flag_is_mode_only':r.get('feedback')})
  for r in rec['ZHANG_E27_USER_TARGET']:
   target_rows.append({'station':case['station'],'mode':case['mode'],**{k:r.get(k,'') for k in ['time','stage','reference','satellite','mean','variance','fractional','perr','selected','provisional_selected']}})
  pos=next(p.glob('*.POS'));pl=pos.read_text();ref=re.search(r'XYZ Reference position\s*:\s*([-\d.e+]+)\s+([-\d.e+]+)\s+([-\d.e+]+)',pl)
  assert ref
  xyz=list(map(float,ref.groups()));snx=next(x['coordinates_m'] for x in s['reference_sinex']['reference'] if x['station']==case['station']);refdiff=math.sqrt(sum((a-snx[b])**2 for a,b in zip(xyz,['STAX','STAY','STAZ'])))
  assert refdiff<0.01,(case['station'],refdiff)
  positions=[];wlpassed=[];l1passed=[]
  for r in diag:
   t=r['time'];a=ar.get(t,{});st=stage.get(t,{})
   wl=st.get('USER_IF_WL_SD',{});l1=st.get('USER_IF_L1_SD_CONDITIONAL',{})
   enu=[float(r[k]) for k in ['east_error_m','north_error_m','up_error_m']];assert all(math.isfinite(x) for x in enu)
   positions.append(enu)
   if wl.get('passed'):wlpassed.append(t)
   if l1.get('passed'):l1passed.append(t)
   epoch_rows.append({'station':case['station'],'mode':case['mode'],'time':t,'east_m':enu[0],'north_m':enu[1],'up_m':enu[2],'horizontal_m':math.hypot(*enu[:2]),'error_3d_m':math.sqrt(sum(v*v for v in enu)),'ambiguity_states':r['ambiguities'],'product_eligible_ambiguity_states':r['integer_valid_ambiguities'],'common_component_satellites':lattice.get(t,{}).get('common_satellites',''),'named_rank':lattice.get(t,{}).get('named_rank',''),'wl_stage_selected':wl.get('selected',0),'wl_stage_passed':bool(wl.get('passed')),'l1_stage_selected':l1.get('selected',0),'l1_stage_passed':bool(l1.get('passed')),'newly_fixed':int(a.get('newly_fixed',0)),'held_integer_rank':a.get('held_integer_rank',''),'ar_summary_present':bool(a)})
  newly=sum(int(r.get('newly_fixed',0)) for r in ar.values());held=max([int(r.get('held_integer_rank',0)) for r in ar.values()]+[0])
  if case['mode']=='float':assert not ar and newly==0
  errors=[l for l in log.splitlines() if l.startswith('Error:')]
  warnings=Counter(l for l in log.splitlines() if l.startswith('Warning:'))
  alarms={k:log.count(k)+trace.read_text().count(k) for k in ['DGEMV','ROLLED_BACK','AUTHORITATIVE MEASUREMENT TRANSACTION REJECTED','AUTHORITATIVE PREDICTION TRANSACTION REJECTED']}
  summary={'station':case['station'],'mode':case['mode'],'completed':receipt['exit_code']==0 and times==expected,'epochs':len(diag),'processing_epochs':len(processing),'wall_seconds':receipt['wall_seconds'],'exit_code':receipt['exit_code'],'ar_attempt_epochs':len(ar),'wl_passed_epochs':len(wlpassed),'first_wl_passed':wlpassed[0] if wlpassed else None,'wl_selected_peak':max([int(r['selected']) for st in stage.values() for name,r in st.items() if name=='USER_IF_WL_SD']+[0]),'l1_passed_epochs':len(l1passed),'newly_fixed_total':newly,'held_integer_rank_peak':held,'fixed_epochs':sum(int(r.get('newly_fixed',0))>0 for r in ar.values()),'east_rms_m':rms([v[0] for v in positions]),'north_rms_m':rms([v[1] for v in positions]),'up_rms_m':rms([v[2] for v in positions]),'rms_3d_m':rms([math.sqrt(sum(x*x for x in v)) for v in positions]),'final_enu_m':positions[-1],'reference_matches_weekly_sinex_within_m':refdiff,'errors':errors,'warnings':dict(warnings),'alarms':alarms,'trace':str(trace)}
  summaries.append(summary);allcases[(case['station'],case['mode'])]=positions
  for f in p.iterdir():
   if f.is_file():artifact_hashes.append({'path':str(f),'bytes':f.stat().st_size,'sha256':sha(f)})
 pairs=[]
 for station in s['stations']:
  ss=station['station'];f=next(r for r in summaries if r['station']==ss and r['mode']=='float');a=next(r for r in summaries if r['station']==ss and r['mode']=='ar')
  diff=max(math.sqrt(sum((x-y)**2 for x,y in zip(v,w))) for v,w in zip(allcases[(ss,'float')],allcases[(ss,'ar')]))
  pairs.append({'station':ss,'complete_pair':f['completed'] and a['completed'],'matched_epochs':min(f['epochs'],a['epochs']),'float_3d_rms_m':f['rms_3d_m'],'ar_3d_rms_m':a['rms_3d_m'],'max_paired_position_difference_m_at_trace_precision':diff,'wl_passed_epochs':a['wl_passed_epochs'],'first_wl_passed':a['first_wl_passed'],'wl_selected_peak':a['wl_selected_peak'],'l1_passed_epochs':a['l1_passed_epochs'],'newly_fixed_total':a['newly_fixed_total'],'held_integer_rank_peak':a['held_integer_rank_peak']})
 assert all(not r['errors'] and not any(r['alarms'].values()) for r in summaries if r['completed'])
 for x in s['snapshot'].values():assert sha(x['path'])==x['sha256']
 result={'expected_epochs':44,'start':'2024-07-17 00:03:00','end':s['end'],'cases':summaries,'paired':pairs,'process_and_epoch_checks_passed':all(r['completed'] for r in summaries),'upstream_numerical_alarm_resolved':False,'independent_integer_truth_available':False,'successful_full_user_ar_demonstrated':any(r['newly_fixed_total']>0 and r['l1_passed_epochs']>0 for r in pairs),'position_reference':'IGS weekly SINEX (POS header matched within 1 cm); independent of R48 estimation, not withheld from initialization or IGS network','note':'feedback=1 on E27 stage logs indicates configured mode; only selected gates and actually committed constraints are counted. No integer wrong-fix probability is inferred from coordinate error or rounding.'}
 (out/'summary.json').write_text(json.dumps(result,indent=2));writecsv(out/'per_epoch.csv',epoch_rows);writecsv(out/'paired_summary.csv',pairs)
 if stage_rows:writecsv(out/'stage_gates.csv',stage_rows)
 if target_rows:writecsv(out/'named_targets.csv',target_rows)
 (out/'output_sha256.json').write_text(json.dumps(artifact_hashes,indent=2))
 report=['# R48 用户端 PPP-AR 配对试验','',f"窗口：2024-07-17 00:03:00–00:24:30，44历元。五个独立用户站、10个用例；8个完成44历元，BREW两组在00:06:30被事务检查拒绝（完成7历元）。",'', '采用自由位置估计、CANONICAL_USER_IF_WL_L1、同一 PRODUCT_FIXED 前缀与完整协方差。两组仅改变 AR 模式和整数反馈开关；成功率0.99、ratio 3、目标 perr≤1e-3、NIS alpha=1e-6 保持原值。','', '|站点|完整配对|完成历元|WL通过历元|WL最大选中数|L1通过历元|实际新增固定|FLOAT 3D RMS(m)|AR 3D RMS(m)|','|---|---|---:|---:|---:|---:|---:|---:|---:|']
 for p in pairs:report.append(f"|{p['station']}|{p['complete_pair']}|{p['matched_epochs']}|{p['wl_passed_epochs']}|{p['wl_selected_peak']}|{p['l1_passed_epochs']}|{p['newly_fixed_total']}|{p['float_3d_rms_m']:.4f}|{p['ar_3d_rms_m']:.4f}|")
 report += ['', '定位误差相对 IGS 周解 SINEX 的先验坐标，已核对 POS 参考坐标一致。用户站未进入 R48 180站产品估计，但参考坐标参与了用户解的初始化；位置自由估计，先验标准差100 m。全窗 RMS 包含冷启动，不能替代长期稳态精度。', '', '上游仍有5次 DGEMV 告警；本次10个用户用例均无 DGEMV；BREW两组有测量事务回滚和致命退出，其余8个用例正常。上游产品的最终数值验收仍未通过。', '', 'WL选中或 feedback=1 不代表用户 IF/L1 整数约束已成功提交；以实际新增固定、保持秩和 L1 门控共同判定。没有独立整数真值，不报告“误固定率为零”。', '', '预检记录：首次冻结请求被拒绝，原因是00:00–00:02:30没有 PRODUCT_FIXED；初次用户运行缺少IF联合协方差回调需要的因子缓冲区；第二次因关闭 STEC 设计列导致 IF 组合未构造。正式用例已启用因子记录并恢复 STEC 设计列，两组一致。所有失败产物保留。', '', '产物：per_epoch.csv（包含完成与失败用例实际已发布历元）、stage_gates.csv、named_targets.csv、paired_summary.csv、summary.json、output_sha256.json。输入数据只引用与哈希，不复制。']
 if not result['successful_full_user_ar_demonstrated']:report[2:2]=['结论：四站完成用户产品接入，BREW两组失败。本次未实现完整用户端 PPP-AR 固定，也未证明定位改善。','']
 report += ['', '## 失败定位与后续工作', '', 'BREW 在00:06:30的独立详细日志复现显示 RAW_SQUARE_ROOT_POSTERIOR_MISMATCH。平方根后验协方差相对误差0.000236645，超过代码固定容差0.0001；均值相对误差约4.23e-9。原滤波后验与平方根重算结果不一致，触发测量事务回滚。当前证据尚不能区分数值病态与实现误差，应保留该历元作为回归用例；不提高容差掩盖失败。', '', 'MARS、DYNG、NICO 的条件 L1 目标最小 perr 分别为0.817209、0.812519、0.674789，均远高于1e-3。三站虽有部分 WL 通过，L1 不确定性仍不足以支持固定。JPLM 没有 WL 通过，未进入条件 L1 阶段。', '', '下一步应先修复 IF 联合协方差与因子记录的配置耦合及 BREW 后验不一致，再检查用户 IF 产品协方差、基准变化和 MARS 米级误差的来源。延长已成熟产品覆盖的用户冷启动窗口，以相同 FLOAT 对照重新检验 L1 收敛；不能用降低成功率、ratio 或 perr 门限获得表面固定。', '', '详细拒绝证据：/mnt/d/GINAN_R20/inputData/outputs/r48_user_20260913_BREW_float_diagnostic/run.log:269。实验源代码与配置步骤均已逐步提交并同步到 fork/codex/r48-continuity-bridge。']
 (out/'R48_user_report.md').write_text('\n'.join(report)+'\n')
 print(json.dumps(pairs,indent=2))
if __name__=='__main__':main()
