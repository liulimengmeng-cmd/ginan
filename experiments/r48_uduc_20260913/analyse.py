#!/usr/bin/env python3
"""Close out all UDUC attempts, preserving failures and rejecting false AR claims."""
import csv,hashlib,json,re,shutil,subprocess
from pathlib import Path
import yaml
W=Path('/mnt/c/Users/rx/Documents/GINAN/r48_uduc_work_20260913')
HERE=Path(__file__).resolve().parent
KV=re.compile(r'([A-Za-z0-9_]+)=([^\s]+)')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def flat(x,p=''):
 out={}
 for k,v in x.items():
  key=(p+'.' if p else '')+k
  if isinstance(v,dict):out.update(flat(v,key))
  else:out[key]=v
 return out
def fields(l):return dict(KV.findall(re.sub(r'time=(\d{4}-\d{2}-\d{2}) (\d\d:\d\d:\d\d)',r'time=\1T\2',l)))
def main():
 s=json.loads((W/'experiment.json').read_text());a=W/'analysis';a.mkdir(exist_ok=True)
 summaries=[];epochs=[];hashes=[];comparison=[];configs={}
 metadata={'outputs.metadata.config_description','processing_options.gnss_general.zhang_pppar.checkpoint_runtime_id'}
 science={'processing_options.gnss_general.zhang_pppar.integer_strategy','processing_options.ppp_filter.ionospheric_components.use_if_combo','processing_options.ppp_filter.ionospheric_components.corr_mode'}
 for c in s['cases']:
  cfg=flat(yaml.safe_load(Path(c['config']).read_text()));old=flat(yaml.safe_load(Path(c['if_config']).read_text()));diff={k:[old.get(k),cfg.get(k)] for k in set(cfg)|set(old) if cfg.get(k)!=old.get(k)}
  assert set(diff)==metadata|science;configs[(c['station'],c['mode'])]=cfg
  assert sha(c['config'])==c['config_sha256'];comparison.append({'case':c['case'],'differences_from_if':diff})
  receipt=json.loads((W/(c['case']+'.result.json')).read_text());p=Path(c['output']);tf=next(p.glob('Network-*.TRACE'));t=tf.read_text();log=re.sub(r'\x1b\[[0-9;]*m','',(p/'run.log').read_text());lines=t.splitlines()
  d=[fields(l) for l in lines if l.startswith('ZHANG_USER_DIAGNOSTIC ')];ar=[fields(l) for l in lines if l.startswith('ZHANG_USER_AR_SUMMARY ')];stage=[fields(l) for l in lines if l.startswith('ZHANG_E26_USER_STAGE ')]
  processing=re.findall(r'Time  (2024-07-17 \d\d:\d\d:\d\d)',t)
  failure=re.findall(r'failed closed at (2024-07-17 \d\d:\d\d:\d\d)',log)
  assert d and all(x['time']=='2024-07-17T00:03:00' for x in d)
  assert processing==['2024-07-17 00:03:00','2024-07-17 00:03:30']
  assert 'IONO_STEC' in t and 'L1C-L2W' not in '\n'.join(l for l in lines if l.startswith('%'))
  for r in d:epochs.append({'station':c['station'],'mode':c['mode'],'time':r['time'],'east_m':r['east_error_m'],'north_m':r['north_error_m'],'up_m':r['up_error_m'],'ambiguity_states':r['ambiguities'],'product_eligible_states':r['integer_valid_ambiguities'],'interpretation':'single committed startup epoch only; no convergence or performance claim'})
  item={'station':c['station'],'mode':c['mode'],'exit_code':receipt['exit_code'],'expected_epochs':44,'committed_diagnostic_epochs':len(d),'attempted_epochs':len(processing),'last_committed_epoch':d[-1]['time'],'failed_epoch':failure[-1] if failure else None,'wall_seconds':receipt['wall_seconds'],'ar_summary_count':len(ar),'ar_stage_count':len(stage),'newly_fixed_total':sum(int(r.get('newly_fixed',0)) for r in ar),'wl_selected_total':sum(int(r.get('selected',0)) for r in stage if r.get('stage')=='USER_WL_SD'),'l1_selected_total':sum(int(r.get('selected',0)) for r in stage if r.get('stage')=='USER_L1_SD_CONDITIONAL'),'DGEMV_count':log.count('DGEMV'),'transaction_rollback_count':t.count('status=ROLLED_BACK'),'last_postfit_passed':t.rfind('Postfit check passed')>t.rfind('Postfit check failed'),'raw_dual_frequency_measurements_verified':True,'stec_states_present':True,'trace':str(tf)}
  summaries.append(item)
  for f in p.iterdir():
   if f.is_file():hashes.append({'path':str(f),'size':f.stat().st_size,'sha256':sha(f)})
 for station in s['stations']:
  st=station['station'];f=configs[(st,'float')];r=configs[(st,'ar')];diff={k for k in set(f)|set(r) if f.get(k)!=r.get(k)}
  assert diff==metadata|{'processing_options.ambiguity_resolution.mode','processing_options.gnss_general.zhang_pppar.canonical_user_target_feedback'}
 for v in s['snapshot'].values():assert sha(v['path'])==v['sha256']
 assert sha(s['binary'])==s['binary_sha256']
 dl=Path(json.loads((W/'diagnostic.json').read_text())['log']);reject=[l for l in dl.read_text().splitlines() if 'event=MEASUREMENT status=REJECTED' in l];assert len(reject)==1
 rejection=fields(reject[0]);assert rejection['failure_reason']=='RAW_SQUARE_ROOT_POSTERIOR_MISMATCH'
 result={'cases':summaries,'attempts_completed':10,'full_window_completed':sum(r['committed_diagnostic_epochs']==44 and r['exit_code']==0 for r in summaries),'scientific_user_ar_comparison_available':False,'products_unchanged':True,'binary_unchanged':True,'if_configuration_differences':comparison,'detailed_mars_rejection':rejection,'replay_covariance_tolerance':1e-4,'upstream_DGEMV_unresolved':True,'claim':'All UDUC cases failed the second measurement transaction before useful AR. This is not evidence that UDUC positioning or fixing is intrinsically worse than IF.'}
 (a/'summary.json').write_text(json.dumps(result,indent=2));(a/'output_sha256.json').write_text(json.dumps(hashes,indent=2))
 for name,rows in [('case_summary.csv',summaries),('per_epoch.csv',epochs)]:
  with (a/name).open('w',newline='') as f:wr=csv.DictWriter(f,fieldnames=list(rows[0]));wr.writeheader();wr.writerows(rows)
 report=['# R48 UDUC 用户端对照试验','', '结论：五站 FLOAT／AR 共10个用例均已执行，但全部在00:03:30测量事务检查失败；仅完成00:03:00一个历元。尚未得到有效 UDUC 固定性能或定位精度对照。','', '计划窗口00:03:00–00:24:30，共44历元。与前一轮 IF 使用相同五个独立用户站、同一 PRODUCT_FIXED 产品与完整协方差、同一冻结二进制；只改变 UDUC 整数策略、IF组合开关与电离层处理模式。保留因子重算检查和全部门限。', '', '|站点|FLOAT完成历元|AR完成历元|两组退出时刻|实际新增固定|','|---|---:|---:|---|---:|']
 for st in s['stations']:report.append('|'+st['station']+'|1/44|1/44|00:03:30|0|')
 report += ['', '## 具体失败证据', '', 'MARS详细日志复现：RAW_SQUARE_ROOT_POSTERIOR_MISMATCH；后验协方差相对差异0.0061983，固定容差0.0001（约62倍）；均值相对差异1.27534e-8。后验残差检查已经通过，拒绝发生在因子重算与滤波后验的一致性审计阶段。只有MARS的详细重算根因已通过信息级日志确认；其他站已确认相同事务拒绝边界，未将其内部原因未经验证地全部等同。', '', '全部10个用例无 DGEMV 日志，均尚未进入 E26 用户整数阶段。首历元真实存在 IONO_STEC 状态与独立 L1C/L2W 观测，确认 UDUC 模型已生效。不能将这次失败解释为宽巷或 L1 门限过严，也不能用单个冷启动历元报告收敛精度。', '', '## 与IF的区别及下一步', '', 'IF轮次有四站完成44历元、BREW于00:06:30失败；UDUC本轮所有站在第二历元失败。两个试验都暴露了后验重算一致性问题，但尚不能据此得出模型优劣结论。下一步应针对00:03:30保存的可复现输入，检查滤波先验、过程噪声传播、相关产品观测噪声及平方根后验更新的一致性，修复后再重跑同条件UDUC对照。不关闭因子检查，不提高容差，不降低固定门限。', '', '用户位置仍自由估计；STEC沿用现有先验与过程噪声，未引入额外电离层约束。UDUC的逐频持续参考星选择策略未修改；候选提取会检查相对产品分量、版本、对齐代次及双频资格。后续性能试验仍需单独评估参考星覆盖限制。', '', '原始日志、实际配置、每个失败用例收据、详细复现日志均保留。per_epoch.csv仅含10条已完成首历元记录。输入数据只引用，不复制。服务端R48运行不受配置或二进制修改。']
 (a/'R48_UDUC_report.md').write_text('\n'.join(report)+'\n')
 for name in ['results','executed_config']:(HERE/name).mkdir(exist_ok=True)
 for f in a.iterdir():
  if f.is_file():shutil.copy2(f,HERE/'results'/f.name)
 for c in s['cases']:
  shutil.copy2(c['config'],HERE/'executed_config'/Path(c['config']).name)
  shutil.copy2(W/(c['case']+'.result.json'),HERE/'results'/(c['case']+'.result.json'))
 for name in ['experiment.json','diagnostic.json']:shutil.copy2(W/name,HERE/'results'/name)
 shutil.copy2(dl,HERE/'results/MARS_diagnostic.log')
 print(json.dumps({'attempted':10,'complete':result['full_window_completed'],'epochs_recorded':len(epochs),'detailed_rejection':rejection},indent=2))
if __name__=='__main__':main()
