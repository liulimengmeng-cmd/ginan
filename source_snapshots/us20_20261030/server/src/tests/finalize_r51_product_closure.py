"""Collect verified results, scoped source delta, hashes and readable report."""
from pathlib import Path
import difflib
import hashlib
import json
import re
import shutil
import sys

root=Path(__file__).resolve().parents[2]
evidence=Path(sys.argv[1]).resolve()
def digest(path):
    value=hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda:stream.read(4*1024*1024),b''):value.update(block)
    return value.hexdigest()
old=json.loads((evidence.parent/'ar_optimization_20261019/one_epoch_probe/real_epoch_audit.json').read_text())
new=json.loads((evidence/'one_epoch_probe/real_epoch_audit.json').read_text())
kernels=json.loads((evidence/'kernel_comparison.json').read_text())
replay=json.loads((evidence/'replay_comparison.json').read_text())
products=json.loads((evidence/'product_epoch_comparison.json').read_text())
probe=json.loads((evidence/'one_epoch_probe/probe.json').read_text())
inputs=Path(probe['binary']).parent/'exact_inputs'
if inputs.is_dir():shutil.copytree(inputs,evidence/'exact_inputs',dirs_exist_ok=True)
assert kernels['all_passed'] and replay['all_passed']
assert new['intentional_interruption'] and not new['hour_complete']
assert '100% tests passed, 0 tests failed out of 39' in (evidence/'ctest.log').read_text()
files=[
    'cpp/common/'+name for name in (
        'zhangIntegerAudit.hpp','zhangProductRelationSolver.hpp','zhangR47ProductDomain.hpp',
        'zhangR48ExactAudit.hpp','zhangR48Marginal.hpp','zhangIntegerCandidateNis.hpp',
        'zhangPhaseTimer.hpp','zhangProductFrameReplay.hpp','zhangR48ValidationCache.hpp')]
files+=['cpp/pea/ppp_ambres.cpp','tests/CMakeLists.txt','tests/run_r51_ar_one_epoch_probe.py',
        'tests/run_r51_product_closure_tests.py','tests/unit/r51_product_closure_contract.cpp',
        'tests/audit_r51_product_closure.py','tests/finalize_r51_product_closure.py']
patch=[];sources={}
for name in files:
    source=root/'src'/name;target=evidence/'source_current/src'/name
    target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(source,target)
    before=evidence/'before'/('tests_CMakeLists.txt' if name=='tests/CMakeLists.txt' else source.name)
    a=before.read_text().splitlines(True) if before.is_file() else []
    b=source.read_text().splitlines(True)
    patch.extend(difflib.unified_diff(a,b,fromfile='before/src/'+name,tofile='after/src/'+name))
    sources['src/'+name]={'before_sha256':digest(before) if before.is_file() else None,'after_sha256':digest(source)}
(evidence/'product_closure.patch').write_text(''.join(patch))
shutil.copy2(root/'bin/pea',evidence/'pea_after')
build=Path('/home/rx/GINAN/build-r51-handoff-20260930')
for source,label in [(build/'cpp/CMakeFiles/pea.dir/flags.make','production_flags.make'),
                     (build/'cpp/CMakeFiles/pea.dir/link.txt','production_link.txt')]:
    shutil.copy2(source,evidence/label)
def phase(a,name):return a['phases'][name]['inclusive_seconds']
names=['R51_PRODUCT_CLOSURE_FINAL','R48_OFFICIAL_PRODUCTS','R48_HISTORY_SEARCH','R48_PRODUCT_IMAGE',
       'R48_HISTORY_SUBSET','R48_ROOT_MOMENTS_FACTOR','AR_TOTAL']
comparison={name:{'before_seconds':phase(old,name),'after_seconds':phase(new,name)} for name in names}
parts=['R51_CLOSURE_PREPARATION','R48_OFFICIAL_PRODUCTS','R51_CLOSURE_POSTSEARCH_AND_PUBLICATION']
assert abs(sum(phase(new,n) for n in parts)+new['phases']['R51_PRODUCT_CLOSURE_FINAL']['exclusive_seconds']-
           phase(new,'R51_PRODUCT_CLOSURE_FINAL'))<1e-6
(evidence/'real_timing_comparison.json').write_text(json.dumps(comparison,indent=2))
lines=[
    '产品闭合优化与验证（2026-10-19）',
    '',
    '运行语义保持：ratio=2；整网联合 NIS 不计算；每历元独立 AR。没有跨历元缓存整数或复用历史 RHS。',
    '未启动一小时实验。只完成一个受控真实历元后主动退出130；并非整小时成功退出0。',
    '首次真实探针因新增重放目录root/rx权限不匹配退出134，未提交产品；失败记录保留且不计为通过。修正探针新目录所有权并预检后，使用同一二进制重新验收。',
    '正常AR结束后的退出收尾中做过一次约5秒栈采样，显示yaml-cpp释放与malloc堆整理；不在AR计时区间内，整进程墙钟不用于严格性能对照。',
    '',
    '旧270.708秒的可证拆分（嵌套分项不可相加）：',
    '闭合总计270.708424 = 正式产品192.162892 + 产品基编译2.769817 + 账本预搜索2.594211 + 外围未细分73.178108 + 其他约0.003396。',
    '正式产品内：条件路由146.715248（其中搜索前整数像60.561487，余86.153762当时未细分）；整数行选择27.311639；根矩4.254545；其他已计时约1.4198；自身未细分12.461685。',
    '独立精确计数器另记7次产品像121.647474秒，包含在上述路由等路径内，不是额外耗时；不能把所有耗时归给首次60.6秒。',
    '旧记录未计时部分无法事后补出精确函数秒数。本次已为HNF、Smith原点、kernel、投影、lift、条件化及闭合外围增加父子计时，保留每次实际矩阵rows/columns/nnz/bits样本。',
    '',
    '代码变化：',
    '1. HNF/Smith/kernel找到首个最小绝对值1后结束主元扫描；只对单位Smith主元省略恒真的整除扫描。非单位主元仍严格保留整除/同余检查。主元、坐标和顺序不变。',
    '2. 仿射域的Smith原点求解同时提供秩和Smith因子，不再额外分解转置系统做零成员查询。单RHS沿原操作带恢复特解，不构造无用完整右变换。',
    '3. 产品条件化把逐特征向量GEMV改成同序批量GEMM；遗留Eigen分解在管理范围内BLAS=1，周围GEMM仍使用多核。',
    '4. 同一不可变FLOAT后验内最多缓存两项数值验证，完整核对约束行/RHS/alpha。证明、风险、整数可行性不因命中而跳过；不按历元号或尺寸盲复用。',
    '',
    '正式编译完成；39/39回归测试通过（新增产品闭合测试，保留确定性周跳/预算取消等测试；未重跑10项长耗时周跳蒙特卡罗）。',
    '220组随机系统×两种产品像路径的新旧输出逐字节一致，覆盖特解、完整kernel、Smith因子、搜索投影/偏移、成员组合；另以独立解析奇偶模型和LDLT条件高斯公式验证整数像与完整协方差。',
    '实际历元稀疏整数输入逐项重放，对照冻结旧后端的全部精确输出；不是调用两遍同一核心函数充作独立实现。',
    '',
    '真实首历元版本对照（非重复成对整历元性能试验）：']
for name,record in comparison.items():lines.append(f"{name}: {record['before_seconds']:.6f} -> {record['after_seconds']:.6f} 秒")
lines+=['','本次闭合直接子阶段（可以相加）：']
for name in parts:lines.append(f"{name}: {phase(new,name):.6f} 秒")
lines.append(f"闭合自身未分配：{new['phases']['R51_PRODUCT_CLOSURE_FINAL']['exclusive_seconds']:.6f} 秒")
trace=(evidence/'one_epoch_probe/Network-case-202419900.TRACE').read_text()
exact=re.search(r'R48_EXACT_TIMING phase=R48_PRODUCT_IMAGE seconds=(\S+) calls=(\d+)',trace)
assert exact
lines += [f"正式产品精确计数器：{exact[2]}次整数像累计{float(exact[1]):.6f}秒，旧记录7次累计121.647474秒；包含在上述直接子阶段内。",
          '剩余热点：闭合后处理/发布，以及正式产品路径自身的未细分开销；本次未宣称整个系统已达到30秒实时要求。',
          '本次232条产品记录及15312项完整联合协方差与修改前首历元逐值一致，发布/分段标志无变化；这是数值等价验证，不是外部精度评价。']
numeric=next(r for r in kernels['runs'] if r['mode']=='marginal')
lines += ['',f"2350维/2233约束条件化内核：{numeric['baseline']['kernel_seconds']:.6f} -> {numeric['current']['kernel_seconds']:.6f} 秒；完整均值/协方差最大差 {numeric['maximum_numeric_difference']:.3g}。",
          f"该独立内核进程峰值RSS：{numeric['baseline']['peak_rss_mib']:.2f} -> {numeric['current']['peak_rss_mib']:.2f} MiB，批量GEMM需要额外工作区；没有宣称所有优化都降峰值内存。",
          '新旧独立头内核测试均使用实际生产宏与链接库、O2构建；完整PEA为原Release优化构建，内核秒数不可当作整历元秒数。']
lines += ['',f"受控800列整数像内核5次中位数：{kernels['frame_median_seconds']['baseline']:.6f} -> {kernels['frame_median_seconds']['current']:.6f} 秒。不是完整AR加速倍数。",
          f"真实首历元峰值RSS：{old['peak_rss_mib']:.2f} -> {new['peak_rss_mib']:.2f} MiB（包含输入加载/FLOAT，不是单独闭合峰值）。峰值略升，不宣称节省整进程内存。",
          '完整4755维FLOAT均值、状态顺序及22610025项协方差逐元素比较相同，详见full_float_compare.log。',
          '一个历元不证明一小时连续性、长期内存稳定性或产品外部精度；不能把高固定率直接当作产品精度。',
          '约束实际改变后的新域仍会重新精确构造；本次移除同一域内冗余转置分解与重复数值验证，并未宣称所有跨版本分解都已增量化。',
          '详见one_epoch_probe/real_epoch_audit.json、real_timing_comparison.json、kernel_comparison.json、replay_comparison.json、product_epoch_comparison.json及TRACE。']
(evidence/'REPORT.txt').write_text('\n'.join(lines)+'\n')
manifest={'scope':'one-epoch closure optimization; not a completed one-hour AR',
          'ratio':2,'joint_nis':False,'epoch_only':True,'sources':sources,
          'binary_before_sha256':digest(evidence/'pea_before'),'binary_after_sha256':digest(evidence/'pea_after'),
          'comparisons':comparison,'products':products,'files':{}}
for path in sorted(evidence.rglob('*')):
    if path.is_file() and path.name not in ('MANIFEST.json','finalize.log'):
        manifest['files'][str(path.relative_to(evidence))]={'bytes':path.stat().st_size,'sha256':digest(path)}
(evidence/'MANIFEST.json').write_text(json.dumps(manifest,indent=2))
print(evidence/'REPORT.txt')
