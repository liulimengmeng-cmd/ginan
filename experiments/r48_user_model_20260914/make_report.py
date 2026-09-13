from pathlib import Path
import json,csv,math,datetime,hashlib,collections
w=Path('/mnt/c/Users/rx/Documents/GINAN/r48_user_model_work_20260914');v=w/'revision2';a=v/'analysis'
s=json.loads((a/'summary.json').read_text());old=json.loads((w/'analysis/summary.json').read_text());f=json.loads((v/'frozen_user.json').read_text())
assert s['all_cases_finished']
by={(x['station'],x['model'],x['mode']):x for x in s['cases']}
oldby={(x['station'],x['model'],x['mode']):x for x in old['cases']}
checks=[]
for station in ['DYNG','NICO','BREW','MARS','JPLM']:
 for model in ['IF','UDUC']:
  x=by[station,model,'float'];y=oldby[station,model,'float']
  assert x['rms_3d_m']==y['rms_3d_m'] and x['final_enu_m']==y['final_enu_m']
  checks.append({'station':station,'model':model,'float_unchanged_at_trace_precision':True})
rows=list(csv.DictReader((a/'per_epoch.csv').open()));oldrows=list(csv.DictReader((w/'analysis/per_epoch.csv').open()))
oldfloat={(r['station'],r['model'],r['time']):r for r in oldrows if r['mode']=='float'}
for r in rows:
 if r['mode']=='float':assert all(r[k]==oldfloat[r['station'],r['model'],r['time']][k] for k in ['east_m','north_m','up_m'])
complete=sum(x['completed'] for x in s['cases']);full=[x for x in s['cases'] if x['mode'].endswith('feedback') and x['integer_factor_l1_rows']>0]
lines=['# R48 产品用户端 IF / UDUC PPP-AR 对照试验','',f'正式修复版：{complete}/50 例完整运行至 00:30:00。每例预期 55 历元，2024-07-17 00:03:00—00:30:00 GPST，30 s 间隔。R48 的 00:00—00:02:30 无 PRODUCT_FIXED，因此不将这些历元计入本次有效用户窗口。','',f'具有实际 L1 因子提交的反馈用例数：{len(full)}。宽巷提交、一般整数关系、具名 L1 提交及位置收益分别统计；本试验没有独立整数真值，不能报告已验证的误固定率。','', '## 试验设置','', '已阅读 OSB 项目“分析解算速度并测试”的最新用户模型分析，来源会话 6aa68697-6fac-83e9-be18-acaf20169431，消息 42ee9bb3-8212-43dc-b6db-3a194b280b38。原文保留在本地 source_user_model_conversation.md。','', '五个用户站 MARS、DYNG、NICO、BREW、JPLM 均未进入 R48 的 180 站产品估计。使用 R48fix1 最终产品及完整产品协方差的冻结副本；输入观测只引用路径和 SHA256，不复制。每种模型设置 FLOAT、ROUND shadow、ROUND feedback、联合 ILS shadow、联合 ILS feedback，共 50 例。单线程、nice 15，编译 nice 19；R49 服务端二进制和配置保持原样。重叠时段存在资源竞争，不用于评价 R49 无干扰运行速度。','', '产品模型按会话核对：码使用 P+C，相位使用 λL+C−B_j；保留完整产品协方差，不重复添加另一套钟差或整数平移。UDUC 估计位置、接收机钟、ZTD、可估 STEC 及模糊度；该 STEC 含码偏差基准，不能直接解释为绝对物理电离层。IF 保留组合前的 STEC 设计列，最终消去电离层；N_IF 本身不是整数，先 WL、再条件 L1。','', '最终具名目标 perr≤0.001、NIS alpha=1e-6 均保持。ROUND 的 success=0.99、ratio=3 沿用基线，其中 ratio 是逐目标距离门控；联合 ILS 明确选 LAMBDA_ALT，success 前筛收紧至0.999，ratio=3 执行候选比值检验。记录实际模式及门控执行状态。一般整数关系先记录，再做精确具名恢复，不能直接冒充完整双频固定。','', '## 正式结果','', '|站点|模型|FLOAT 3D RMS m|ROUND反馈 RMS m|ILS反馈 RMS m|ROUND WL/L1通过历元|ILS WL/L1通过历元|ROUND/ILS L1提交行数|','|---|---|---:|---:|---:|---:|---:|---:|']
for station in ['DYNG','NICO','BREW','MARS','JPLM']:
 for model in ['IF','UDUC']:
  x=by[station,model,'float'];r=by[station,model,'round_feedback'];i=by[station,model,'ils_feedback']
  def val(q):return f"{q['rms_3d_m']:.4f}"+('（未完成）' if not q['completed'] else '')
  lines.append(f"|{station}|{model}|{val(x)}|{val(r)}|{val(i)}|{r['wl_passed_epochs']}/{r['l1_passed_epochs']}|{i['wl_passed_epochs']}/{i['l1_passed_epochs']}|{r['integer_factor_l1_rows']}/{i['integer_factor_l1_rows']}|")
lines+=['', '表内 RMS 为各用例已发布历元；若未完成，不直接据此作全窗优劣比较，匹配历元结果见 paired_summary.csv。位置参考为 IGS 周解 SINEX，已核对 POS 头坐标差小于 1 cm；该参考参与初始化，不能称为完全盲测。全窗包含冷启动，不是长期稳态精度。影子用例与各自 FLOAT 的所有匹配坐标在 TRACE 输出精度下一致。','', '## 已修复及验证','', '1. 数值更新：原 MARS UDUC 第二历元及 BREW IF 00:06:30 失败均先复现。冻结最终测量包 x/P/H/v/R/xp/Pp，以独立白化 QR 和 Joseph 公式重算。首次创新协方差条件数分别约 1.24e11 和 5.51e10；普通更新对 QR 的协方差相对差约 6.09e-7 和 5.38e-7，Joseph 更新分别约 2.62e-14 和 4.37e-14。后续失败历元局部更新与其自身先验较一致，因此不能把累计后验不一致直接说成该历元 H/R 构造错误。启用 Joseph 后，十个 FLOAT 用例全部完成；检查容差未改变。','', '2. 联合 ILS 拒绝路径：首次对照的 20 个 ILS 用例因新增诊断读取无对应整数值的残留变换矩阵而段错误。拒绝搜索后的变换基或拒绝候选不具备整数授权，现将其清空后再记录和恢复具名目标，防止越界及错误晋升。最初失败产物与首版二进制均保留。','', '3. 用户反馈提交链：首次对照的五个 UDUC ROUND feedback 在首次 WL 条件更新的下一历元失败。原路径只更新 KF x/P，未同步因子记录后验。修复后，在发布用户后验前提交精确条件矩映射 F=I−KA、b=Kz、Q=0，并推进同一提交序号；标签 USER_INTEGER_CONDITION 表示条件更新，不算新的独立观测。映射通过原有后验一致性检查后才提交。此修改限于用户反馈分支，不修改正在运行的 R49。','', '最初 50 例中 25 例完成、20 例 ILS 段错误、5 例 UDUC 反馈链失败。正式修复版重新执行全部 50 例；十个 FLOAT 每个历元的坐标与首版完全一致（TRACE 精度）。','', '## 仍不能据此宣称的问题已解决','', 'MARS IF 的 00:12:30 坐标放松仍复现：3D误差约0.397→13.423 m，没有AR反馈。它是 FLOAT/QC 链的问题，不能归因于错误整数固定。UDUC 基线在该窗口更好，并不构成对所有数据与长期性能的普遍结论。','', '本次完成数值回归、FLOAT 配对、显式算法对照及反馈链修复与回放；原始 MW 累积和产品整数资格生命周期的完全分离、多连通分量各自参考星及联合跨分量搜索尚未实现。跨历元产品协方差仍不可得，也未加入区域大气改正；本次是第一层钟差/相位产品驱动的用户测试，不能作为快速 PPP-RTK 性能证明。','', '## 来源与产物','',f"源提交：`{f['source_commit']}`；源目录 SHA256：`{f['source_sha256']}`；二进制 SHA256：`{f['binary_sha256']}`。",'', '正式数据：revision2/analysis/per_epoch.csv（逐历元）、stage_gates.csv（逐阶段门控）、named_targets.csv（目标均值/方差/perr）、ils_relations.csv（一般关系及具名恢复）、paired_summary.csv（匹配位置对照）、summary.json 和 output_sha256.json。原失败版位于上级 analysis/；数值独立重算见 independent_all_measurements.json。']
lines += ['', '## 对结果的具体解释', '', 'BREW UDUC 的 ROUND 在00:23:00提交1条条件L1（G10−G32），联合ILS在00:25:00提交1条条件L1；两组均完成后续历元。其全窗3D RMS分别为0.2793 m和0.2808 m，FLOAT为0.2826 m，改善幅度很小。这是有限整数方向上的实际反馈，不能表述为全部卫星双频固定或可靠误固定率已获验证。', '', 'NICO UDUC 的联合ILS在10个历元生成通过搜索和关系门控的一般L1关系（逐历元行数累计12，不能解释为12条独立新增信息），但没有恢复出通过最终具名L1门控的目标。因而不能笼统地说所有联合整数方向都没有信息；当前具名恢复与反馈范围仍限制了结果。', '', 'reported_newly_fixed_sum 只保留原日志计数；原实现会重复计入已经满足的约束，不代表独立新增固定。实际提交另见 integer_commits.csv。保持秩是条件协方差的秩诊断，也不等同于整数真值验证。所有50例正式运行的DGEMV及致命测量/预测事务告警均为0。']
failed=[x for x in s['cases'] if not x['completed']]
if failed:
 lines+=['','## 正式版未完成用例','']
 for x in failed:lines.append(f"- {x['case']}：{x['epochs']}历元，退出{x['exit_code']}；{x['errors']}")
(w/'R48_user_model_report.md').write_text('\n'.join(lines)+'\n')
(a/'regression_checks.json').write_text(json.dumps({'matched_float_epochs_checked':550,'float_controls':checks,'completed':complete,'actual_l1_feedback_cases':[x['case'] for x in full]},indent=2)+'\n')
print(json.dumps({'completed':complete,'actual_l1_feedback_cases':[x['case'] for x in full]},indent=2))
