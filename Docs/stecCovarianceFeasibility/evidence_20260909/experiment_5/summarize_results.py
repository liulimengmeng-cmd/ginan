import json
from pathlib import Path
import sys
sys.path.insert(0,'/mnt/d/tec/ginan-main-code-products/scripts')
from experiment5_code_osb_sigma import sha, write_json

r=Path('/mnt/d/tec/code-exp5-sigma-20260909-v2')
d=Path('/mnt/d/tec/ginan-main-code-products/Docs/stecCovarianceFeasibility')
a=r/'audits'
def read(p): return json.loads(p.read_text())
models=[]
for mm in (0,3,10):
    if mm==0:
        f=read(d/'experiment_4c_feedback_information.json')
        c=read(d/'experiment_4cd_coordinate_integrity_summary.json')
        s=read(a/'baseline_zero_screening.json')
        b=read(d/'experiment_4c_external_bias_application.json')
        t=read(d/'experiment_4c_dual_frequency_datum.json')
    else:
        f=read(a/f'sigma_{mm}mm_ar_feedback.json')
        c=read(a/f'sigma_{mm}mm_coordinate.json')
        s=read(a/f'sigma_{mm}mm_ar_screening.json')
        b=read(a/f'sigma_{mm}mm_ar_bias.json')
        t=read(a/f'sigma_{mm}mm_ar_datum.json')
        assert read(a/f'sigma_{mm}mm_float_feedback.json')['feedback_epoch_count']==0
    assert c['data_quality']['primary_float_exact_epoch_coverage_fraction']==1
    assert c['post_burn_in_statistics']['overall']['primary_vs_sinex_3d']['count']==7200
    models.append({'sigma_mm':mm,'feedback':f,'coordinate':c['post_burn_in_statistics'],
                   'screening':s,'bias':b,
                   'datum':{k:v for k,v in t.items() if not isinstance(v,list)},
                   'verified_fixed_rate':None,'certified_fixed_STEC':False,
                   'scientific_PPPAR_acceptance':False})

inputs=read(r/'all_inputs.json')
preservation=read(a/'original_input_preservation.json')
assert preservation['all_unchanged']
receipt={'schema':'GINAN_CODE_PHASE_OSB_SIGMA_SENSITIVITY_V1','baseline_commit':'91c1ceae9ecf1192419d1db83841d1751ef7cb37',
         'plan':read(r/'plan.json'),'original_inputs':inputs,'preservation':preservation,
         'models':models,'run_receipts':{p.name:read(p) for p in r.glob('*.finished.json')},
         'tests':{'selected_unittests':27,'passed':27,'log':str(r/'tests.log')},
         'source_hashes':read(r/'source_file_hashes.json'),
         'sigma_10mm_nis_extremes':read(r/'sigma_10mm_ar_nis_extremes.json'),
         'audit_hashes':{p.name:sha(p) for p in a.glob('*.json')},
         'claim_limits':['Assumed sigma, not calibrated CODE accuracy','Feedback is SUBMITTED_UNVERIFIED','Full-state shadow norm has mixed units','One day, six sites; external weekly SINEX without velocity propagation','No independent integer truth or certified fixed STEC']}
write_json(d/'experiment_5_sigma_run_receipt.json',receipt)

lines=['# CODE phase-OSB 非零先验不确定度受控实验（2024-07-17）','',
'## 结论','',
'完成 3 mm、10 mm 两档日持久 phase-OSB 先验的全天 AR/FLOAT 配对实验。非零先验改变了候选选择、反馈活动与 NIS，但本实验仍未建立可信固定率，不能认证 PPP-AR 或 fixed STEC。下表分别列出统计诊断与外部坐标变化；人为 sigma 不等于经过校准的 CODE 精度。','',
'源码和原始 trace 还纠正了旧报告的一项解释：shadow 全状态差范数混合多种单位，不能标成米；大的差异与实际滤波中的后验筛选失败同时出现。没有筛选失败时，零方差及非零方差模型均接近数值舍入量级。因此不能把该范数直接解释成非线性误差，也不能把它作为零 OSB 方差导致物理误差的独立证据。','',
'## 反馈活动与统计诊断','',
'| 先验 sigma | 网络反馈历元/2880 | 六站反馈历元/17280 | 提交秩覆盖 | NIS/row 中位数 | P95 | 最大值 |',
'|---|---:|---:|---:|---:|---:|---:|']
for m in models:
    f=m['feedback']; n=f['nis_per_row']
    lines.append(f"| {m['sigma_mm']} mm | {f['feedback_epoch_count']} ({100*f['feedback_epoch_incidence']:.2f}%) | {f['station_feedback_epoch_count']} ({100*f['station_epoch_feedback_incidence']:.2f}%) | {100*f['submitted_rank_coverage']:.2f}% | {n['median']:.4f} | {n['p95']:.4f} | {n['max']:.4f} |")
lines += ['', '10 mm 档 NIS/row 最大值 592.727 出现在第 2503 历元（20:51:00 UTC），对应 20 行约束，处于预设分析窗口内，不能简单归为启动瞬态。其 NIS 均值为 6.092，甚至高于 3 mm 档的 5.921：仅看中位数会掩盖恶化的尾部。详见原始 trace 和 `sigma_10mm_ar_nis_extremes.json`。', '', '上述反馈仅表示 `FILTER_CALL_RETURNED_SUBMITTED_UNVERIFIED`。候选、提交、经过滤波筛选后的实际更新和整数真值是不同层次。可信固定率目前不可估计，不能把反馈活动率改称固定率。NIS 为筛选前的联合线性诊断，没有额外 NIS 拒绝门；候选选择及重复约束使其不能直接充当独立卡方检验样本。两档 FLOAT 的反馈提交均为 0。', '',
'## 14:00–24:00 UTC 正向坐标精度','',
'单位为 mm；参考为同一 IGS 周解 SINEX，不推断站速传播。每档 AR/FLOAT 精确时间匹配覆盖率为 100%，窗口内各 7200 个站历元。最后一列正数表示 AR 比自身 FLOAT 更差。', '',
'| sigma | AR 3D RMS | FLOAT 3D RMS | AR 水平 RMS | FLOAT 水平 RMS | AR 垂直 RMS | FLOAT 垂直 RMS | AR−FLOAT 3D RMS |',
'|---|---:|---:|---:|---:|---:|---:|---:|']
for m in models:
    o=m['coordinate']['overall']; val=lambda k:o[k]['rms_m']*1000
    values=[val(k) for k in ['primary_vs_sinex_3d','float_vs_sinex_3d','primary_vs_sinex_horizontal','float_vs_sinex_horizontal','primary_vs_sinex_vertical','float_vs_sinex_vertical']]
    lines.append('| '+str(m['sigma_mm'])+' mm | '+' | '.join(f'{v:.4f}' for v in values)+f' | {values[0]-values[1]:+.4f} |')
lines += ['', '### 分站三维 RMS（mm）','', '| sigma | 站 | AR | FLOAT | AR−FLOAT | 首次反馈（分钟） |', '|---|---|---:|---:|---:|---:|']
for m in models:
    for rec,c in sorted(m['coordinate']['by_receiver'].items()):
        ar=c['primary_vs_sinex_3d']['rms_m']*1000; fl=c['float_vs_sinex_3d']['rms_m']*1000
        first=m['feedback']['receiver_summary'][rec]['first_feedback_minutes_after_start']
        lines.append(f"| {m['sigma_mm']} mm | {rec} | {ar:.4f} | {fl:.4f} | {ar-fl:+.4f} | {first} |")
lines += ['', '该单日、六站实验不能证明跨日稳定收益；MARS 垂向参考差对网络 RMS 的影响须保留。网络平均值与分站符号必须同时阅读；若 sigma 改变 FLOAT 本身，不能把这一变化计为 AR 收益。', '',
'## Shadow 与实际筛选','',
'| sigma | 反馈调用 | 有后验筛选失败的调用 | 后验失败事件 | 无筛选失败时最大范数 | 有筛选失败时最大范数 |',
'|---|---:|---:|---:|---:|---:|']
for m in models:
    s=m['screening']
    lines.append(f"| {m['sigma_mm']} mm | {s['feedback_calls']} | {s['calls_with_screening_failure']} | {s['postfit_failure_events']} | {s['mixed_unit_norm_without_screening']['max']:.4g} | {s['mixed_unit_norm_with_screening']['max']:.6g} |")
lines += ['', '两列范数均无统一物理单位，不是坐标位移。`ppp_ambres.cpp:1158` 的 shadow 使用筛选前协方差和创新；`:1248` 调用实际滤波；`algebra.cpp:3060` 执行后验检查；`ppp_callbacks.cpp:66` 可改变测量协方差。原始 trace 中的条件分组支持这一解释，但没有输出逐状态 shadow 误差向量，不能追溯每个状态的误差分量。', '',
'## 候选与门控','', '| sigma | 候选行 | 候选网络历元 | 全日结构审计通过 |', '|---|---:|---:|---|']
for m in models:
    t=m['datum']
    lines.append(f"| {m['sigma_mm']} mm | {t.get('candidate_row_count')} | {t.get('candidate_epoch_count')} | {t.get('structural_basis_pass')} |")
lines += ['', '完整失败门计数如下；未降低成功率 0.9999 或 ratio 3。结构审计并非整数真值认证。', '', '| 首个失败门/通过状态 | 0 mm | 3 mm | 10 mm |', '|---|---:|---:|---:|']
keys=sorted(set().union(*(m['feedback']['first_failure_gate_counts'] for m in models)))
for k in keys:
    lines.append('| '+k+' | '+' | '.join(str(m['feedback']['first_failure_gate_counts'].get(k,0)) for m in models)+' |')
lines += ['', '## 不确定度注入、数据及保全','',
'同一 CODE rapid SP3/CLK/ERP/OSB 产品链，保持全部偏差均值、观测、门限和分析窗口。只在独立命名的实验 BIA 副本中改变 64 条 GPS L1C/L2W 标准差字段；原始文件未改。每个卫星/信号先验独立，全天和六站共享，Q=0。G01 排除后实际使用 62 个状态，四个新解中均只初始化一次；运行方差与副本编码一致，没有把 phase 先验重复注入为逐历元观测噪声。', '',
'四个新解均正常结束，原始 23 个输入的前后 SHA-256 一致；外部偏差审计无缺失、无重复。完整逐文件清单（6 个 RINEX、3 个 SP3、CLK、BIA、ERP、BRDC、ANTEX、周解 SINEX、IGRF14、DE436、GPT2.5、OLOAD、ALOAD、极潮负荷、卫星 metadata、yaw 表）见 `EXPERIMENT_5_CODE_PHASE_OSB_DATA_MANIFEST.md`，不含 WUM/CAS 混用或 ORBEX。', '',
'| 新解 | 外部偏差记录 | 初始化状态数 | 运行秒数 |', '|---|---:|---:|---:|']
for label in ['sigma_3mm_ar','sigma_3mm_float','sigma_10mm_ar','sigma_10mm_float']:
    b=read(a/f'{label}_bias.json'); run=read(r/f'{label}.finished.json')
    lines.append(f"| {label} | {b['marker_count']} | {b['state_prior_initialisation_count']} | {run['elapsed_seconds']:.3f} |")
lines += ['', '## 交付与复核','',
'- 设计与可重复命令：`EXPERIMENT_5_CODE_PHASE_OSB_SIGMA_DESIGN_20260909.md`。',
'- 汇总运行收据：`experiment_5_sigma_run_receipt.json`，含全部输入哈希、审计哈希、候选/反馈/精度统计和声明边界。',
'- 原始 trace、stdout/stderr、逐运行开始/结束收据、冻结 PEA、实验 BIA、配置与审计：`D:/tec/code-exp5-sigma-20260909-v2`。',
'- 源码新增：`scripts/experiment5_code_osb_sigma.py`、`scripts/audit_experiment5_code_osb_sigma.py`、`scripts/test_experiment5_code_osb_sigma.py`；核心估计器未改。',
'- 测试：27 项选定单元测试通过；5 历元冒烟运行通过；四个全天运行退出码 0。通用 unittest discovery 曾触发无关 GUI 的 PySide6 缺失，改为显式运行上述相关测试模块后全部通过；原缓存保留。',
'- 基线分支/提交仍为 `codex/code-products-pppar` / `91c1cea`。本轮新文件留在工作区供审阅，未提交或推送；其他工作树、旧日志和未跟踪实验成果未删除、归档或覆盖。', '',
'后续应先诊断 NIS 尾部和候选一致性，再以外部证据约束跨信号/跨卫星相关性及跨日表现。不能继续按坐标结果挑选 sigma、降低门限，或将本轮敏感性改善写成已实现可信固定。', '']
(d/'PPPAR_CODE_PHASE_OSB_NONZERO_SIGMA_RESULTS_20260909.md').open('x',encoding='utf-8').write('\n'.join(lines))
print(json.dumps([{ 'sigma_mm':m['sigma_mm'], 'feedback_station_rate':m['feedback']['station_epoch_feedback_incidence'], 'nis':m['feedback']['nis_per_row'], 'coordinates':m['coordinate']['overall']} for m in models],indent=2))
