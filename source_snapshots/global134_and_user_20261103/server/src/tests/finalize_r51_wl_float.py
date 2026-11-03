"""Finalize only completed, passing evidence. Does not launch/restart a solver.
WORK_ROOT EVIDENCE_DIRECTORY. All comparisons must already exist.
"""
from pathlib import Path
import difflib
import hashlib
import json
import re
import shutil
import sys

work,evidence=map(lambda p:Path(p).resolve(),sys.argv[1:3])
repo=Path(__file__).resolve().parents[2]
previous=repo.parent/'artifacts/product_closure_optimization_20261019'
def load(path):return json.loads(path.read_text())
selected=load(work/'candidate_selection.json') if (work/'candidate_selection.json').exists() else {}
probe_name=selected.get('probe_directory','one_epoch_probe')
probe_work=work/probe_name
probe_evidence=evidence/probe_name
def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda:stream.read(4*1024*1024),b''):h.update(block)
    return h.hexdigest()

suite=load(work/'suite_comparison_accepted.json')
kernels=load(work/selected.get('kernel_directory','paired_kernels')/'results.json')
comparison=load(work/selected.get('epoch_comparison','epoch_comparison.json'))
real_replays=load(work/selected.get('real_replay_report','real_conditioner_replays.json'))
old=load(previous/'one_epoch_probe/real_epoch_audit.json')
new=load(probe_evidence/'real_epoch_audit.json')
ctest=(work/'ctest.log').read_text()
full=(work/selected.get('full_float_comparison_log','full_float_compare.log')).read_text()
suffix='_'+selected['tag'] if selected.get('tag') else ''
full_check=load(work/('full_checkpoint_comparison'+suffix+'.json'))
assert all(x['passed'] for x in (suite,kernels,comparison))
assert all(x['passed'] for x in real_replays.values())
assert full_check['passed'] and full_check['new_binary_sha256']==new['probe']['binary_sha256']
assert full_check['metrics']['mean_component_scaled']<=2e-10
assert full_check['metrics']['covariance_entry_scaled']<=2e-11
assert Path(comparison['new_bundle'])==Path(new['probe']['committed_epochs'][0]).parent
assert '100% tests passed, 0 tests failed out of 40' in ctest
assert 'REAL_FLOAT_FULL_COMPARE' in full and 'state_order_identical=1' in full
assert new['probe']['config_sha256']==old['probe']['config_sha256']
assert not new['probe']['solver_still_running'] and new['probe']['returncode']==130
assert old['probe']['binary_sha256']==digest(work/'baseline/pea')
assert new['probe']['binary_sha256']==digest(work/selected.get('binary_directory','verified')/'pea')
exact=[]
for path in sorted((probe_work/'exact_inputs').glob('frame_*.txt')):
    reference=previous/'exact_inputs'/path.name
    assert reference.is_file() and digest(reference)==digest(path),'exact integer frame changed'
    exact.append(dict(name=path.name,sha256=digest(path)))
assert len(exact)==8,'missing exact integer frames'

def phase(a,name):return a['phases'][name]['inclusive_seconds']
timings={name:dict(old_seconds=phase(old,name),new_seconds=phase(new,name)) for name in (
    'AR_TOTAL','R51_WL_BLOCK_SEARCH','R51_L1_BLOCK_SEARCH','R51_WL_CONDITION','R51_L1_CONDITION',
    'PRODUCT_RELATION_SEARCH_TOTAL') if name in old['phases'] and name in new['phases']}
trace=(probe_evidence/'Network-case-202419900.TRACE').read_text()
old_trace=(previous/'one_epoch_probe/Network-case-202419900.TRACE').read_text()
float_timings={}
for name in ('FLOAT_TOTAL','FLOAT_KF_UPDATE'):
    values=[]
    for text in (old_trace,trace):
        match=re.search(r'ZHANG_FLOAT_PHASE_TIMER phase='+name+r' .*?inclusive_ms=([0-9.eE+-]+)',text)
        assert match,'missing FLOAT phase '+name
        values.append(float(match[1])/1000)
    float_timings[name]=dict(old_seconds=values[0],new_seconds=values[1])
conditioner=[line for line in trace.splitlines() if line.startswith('ZHANG_R51_CONDITIONER ')]
assert len(conditioner)==2
assert 'evaluation=SKIPPED_RATIO_ONLY' in trace and 'prior_integer_reuse=0' in trace
records=re.findall(r'ZHANG_FLOAT_KERNEL phase=(\w+) states=(\d+) measurements=(\d+) elapsed_ms=([0-9.eE+-]+)',trace)
assert records,'missing detailed FLOAT kernel timing'
summary=dict(passed=True,scope='One authenticated real epoch plus independent/frozen kernel comparisons; not an hour or real-time claim',
    settings=dict(ratio=2,joint_nis=False,epoch_only_ar=True),
    baseline_binary_sha256=old['probe']['binary_sha256'],candidate_binary_sha256=new['probe']['binary_sha256'],
    ctest_passed=40,long_cycle_slip_monte_carlo_not_rerun=10,suite=suite,
    exact_integer_frames=exact,kernel_comparisons=kernels['summary'],
    full_float_kernel_comparison=kernels['full_float_comparison'],real_epoch_products=comparison,
    real_conditioner_replays=real_replays,
    real_epoch_timings=timings,real_epoch_float_timings=float_timings,conditioner_records=conditioner,
    float_kernel_timings=[dict(phase=p,states=int(n),measurements=int(m),milliseconds=float(t)) for p,n,m,t in records],
    full_float_checkpoint_comparison=full,peak_rss_mib=dict(old=old['peak_rss_mib'],new=new['peak_rss_mib']),
    full_hour_started=False,solver_left_running=False,
    known_limits=['Cycle-slip scientific readiness unchanged','No multi-epoch continuity or memory-growth acceptance',
                  'Real epoch runtime uses an authenticated historical baseline, not repeated paired whole-epoch timing'])
failed_log=work/'full_float_compare.log'
if failed_log.exists() and 'equivalence failed' in failed_log.read_text():
    summary['rejected_unconditional_symmetric_candidate']=dict(
        binary_sha256=load(work/'one_epoch_probe/probe.json')['binary_sha256'],
        reason='Full FLOAT covariance equivalence failed despite close published products',
        evidence_log=str(failed_log),source_snapshot='failed_symmetric_candidate_source')
summary['selected_candidate']=selected
summary['full_checkpoint_comparison_metrics']=full_check['metrics']
summary['float_covariance_mode_records']=[line for line in trace.splitlines()
    if line.startswith('ZHANG_FLOAT_COVARIANCE_MODE ')]

evidence.mkdir(parents=True,exist_ok=True)
for path in work.glob('*.log'):shutil.copy2(path,evidence/path.name)
for path in work.glob('*.json'):shutil.copy2(path,evidence/path.name)
shutil.copytree(work/selected.get('kernel_directory','paired_kernels'),evidence/'paired_kernels',dirs_exist_ok=True,
                ignore=shutil.ignore_patterns('*.posterior'))
shutil.copytree(probe_work/'condition_inputs',evidence/'condition_inputs',dirs_exist_ok=True)
shutil.copytree(probe_work/'exact_inputs',evidence/'exact_inputs',dirs_exist_ok=True)
for mode,name in (('baseline','pea_before'),(selected.get('binary_directory','verified'),'pea_after')):
    shutil.copy2(work/mode/'pea',evidence/name)
changes=[]
for before in (evidence/'before').rglob('*'):
    if not before.is_file():continue
    relative=before.relative_to(evidence/'before');current=repo/relative
    if current.is_file():
        changes.extend(difflib.unified_diff(before.read_text().splitlines(True),current.read_text().splitlines(True),
            fromfile='before/'+relative.as_posix(),tofile='after/'+relative.as_posix()))
(evidence/'task_only_existing_files.patch').write_text(''.join(changes))
(evidence/'summary.json').write_text(json.dumps(summary,indent=2))
lines=['WL/L1 条件化与 FLOAT KF 优化验证','',
       '验收通过。ratio=2；整网联合 NIS 关闭；每历元独立 AR。',
       '范围：一次 02:00:00 真实历元；不是完整一小时，也不是实时性或多历元连续性验收。','',
       f'回归：40/40；10 项未修改的长时周跳 Monte Carlo 未重跑。',
       f'条件化：{suite["records"]} 次；{suite["compared_entries"]} 个数值；状态差异 {len(suite["status_changes"])}。',
       f'真实 WL/L1 后验重放：{sum(r["compared_entries"] for r in real_replays.values())} 个均值/协方差数值，逐次接受/拒绝和秩一致。',
       f'完整 FLOAT 核心：{kernels["full_float_comparison"]["entries"]} 个均值/协方差数值。',
       '完整真实 FLOAT 检查点：4755 状态，22610025 协方差元素，含逐元素尺度检查。',
       f'产品：{comparison["product_rows"]} 行；联合协方差 {comparison["covariance_entries"]} 项；发布/分段标志不变。',
       '8 份精确整数系统输入逐字节一致；接受分块决策一致。','',
       '历史基线 / 本次受控历元耗时（父子时间不可重复相加）：']
for name,value in timings.items():lines.append(f'{name}: {value["old_seconds"]:.6f} -> {value["new_seconds"]:.6f} s')
for name,value in float_timings.items():lines.append(f'{name}: {value["old_seconds"]:.6f} -> {value["new_seconds"]:.6f} s')
lines+=['','条件化分解记录：',*conditioner,'','同机交替内核测试：']
for name,value in kernels['summary'].items():lines.append(f'{name}: median={value["median_seconds"]:.6f} s; repeats={value["repetitions"]}')
lines+=['',f'峰值 RSS (MiB): {old["peak_rss_mib"]} -> {new["peak_rss_mib"]}',
        '保留了病态矩阵中间失败记录；最终通过报告为 suite_comparison_accepted.json。',
        '无条件对称更新候选未通过完整 FLOAT 协方差检查，未被作为最终版本。',
        '最终版在条件数/噪声尺度不适合时回退旧公式，不更改统计验收门槛。',
        *summary['float_covariance_mode_records'],
        '秩亏/病态初始商空间、秩或 Gram 歧义均保留旧版严格检查；并非所有输入都会走缓存快速路径。',
        '缓存仅属同一历元同一搜索阶段；不复用跨历元整数。',
        '可用 ZHANG_R51_BLOCK_CONDITIONER=legacy / ZHANG_R51_FLOAT_MATRIX_MODE=legacy 做受控回退。',
        '不声称修复周跳检出延迟，不声称通过多历元连续性、长期内存增长或完整一小时实验。',
        '解算已退出，没有遗留运行中的 AR。','',full]
(evidence/'REPORT.txt').write_text('\n'.join(lines),encoding='utf-8')
manifest=[]
for path in sorted(evidence.rglob('*')):
    if path.is_file() and path.name!='MANIFEST.json':
        manifest.append(dict(path=str(path.relative_to(evidence)),bytes=path.stat().st_size,sha256=digest(path)))
(evidence/'MANIFEST.json').write_text(json.dumps(manifest,indent=2))
assert all((evidence/item['path']).stat().st_size==item['bytes'] and digest(evidence/item['path'])==item['sha256'] for item in manifest)
print(json.dumps(dict(passed=True,evidence=str(evidence),files_verified=len(manifest),timings=timings),indent=2))
