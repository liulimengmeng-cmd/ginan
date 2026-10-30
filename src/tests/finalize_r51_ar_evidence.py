"""Freeze this optimization's sources, exact before/after diff and evidence.
Usage: python3 finalize_r51_ar_evidence.py ARTIFACT_DIRECTORY
"""
from pathlib import Path
import datetime
import difflib
import hashlib
import json
import re
import shutil
import statistics
import sys

root = Path(__file__).resolve().parents[2]
artifact = Path(sys.argv[1]).resolve()
assert artifact.name == "ar_optimization_20261019" and (artifact / "before").is_dir()
results = json.loads((artifact / "paired_final/results.json").read_text())
condition_profiles = []
for line in (artifact / "condition_adaptive_comparison.log").read_text().splitlines():
    if line.startswith("CONDITION "):
        condition_profiles.append(dict(re.findall(r"(\w+)=(\S+)", line)))
condition_medians = {name: statistics.median(float(p["seconds"]) for p in condition_profiles
    if (p.get("adaptive") == "1") == enabled)
    for name, enabled in (("unmanaged_12", False), ("adaptive_12_1", True))}
ctest = (artifact / "ctest_adaptive_final.log").read_text()
assert "100% tests passed, 0 tests failed out of 38" in ctest
assert "Built target unit_tests" in (artifact / "build_shutdown_final.log").read_text()
assert all(c.get("byte_identical", True) for c in results["comparisons"])

old_files = {
    "CMakeLists.txt": "src/tests/CMakeLists.txt",
    "main.cpp": "src/cpp/pea/main.cpp",
    "ppp_ambres.cpp": "src/cpp/pea/ppp_ambres.cpp",
}
for name in (artifact / "before").glob("*.hpp"):
    old_files[name.name] = "src/cpp/common/" + name.name
new_files = [
    "src/cpp/common/zhangArThreading.hpp",
    "src/cpp/common/zhangPrivateUnitMembership.hpp",
    "src/cpp/common/zhangEpochStopRequest.hpp",
    "src/tests/unit/r51_ar_hotspot_benchmark.cpp",
    "src/tests/unit/r51_ar_optimization_contract.cpp",
    "src/tests/unit/r51_ar_condition_profile.cpp",
    "src/tests/unit/r51_epoch_stop_contract.cpp",
    "src/tests/unit/r51_checkpoint_probe_verify.cpp",
    "src/tests/run_r51_ar_hotspots.py",
    "src/tests/audit_r51_ar_hotspots.py",
    "src/tests/run_r51_ar_one_epoch_probe.py",
    "src/tests/audit_r51_ar_real_epoch.py",
    "src/tests/compare_r51_product_epoch.py",
    "src/tests/finalize_r51_ar_evidence.py",
]
patch = []
for old_name, relative in old_files.items():
    before = (artifact / "before" / old_name).read_text().splitlines(True)
    after = (root / relative).read_text().splitlines(True)
    patch.extend(difflib.unified_diff(before, after, fromfile="before/" + relative,
                                      tofile="after/" + relative))
for relative in new_files:
    patch.extend(difflib.unified_diff([], (root / relative).read_text().splitlines(True),
                                      fromfile="/dev/null", tofile="after/" + relative))
(artifact / "ar_optimization.patch").write_text("".join(patch))


def digest(path):
    value = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


source_hashes = {}
for relative in list(old_files.values()) + new_files:
    target = artifact / "source_current" / relative
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(root / relative, target)
    source_hashes[relative] = digest(target)
shutil.copy2(root / "bin/pea", artifact / "pea_after")
build = Path("/home/rx/GINAN/build-r51-handoff-20260930")
build_hashes = {}
for name in ("flags.make", "link.txt"):
    destination = artifact / "build_inputs" / name
    destination.parent.mkdir(exist_ok=True)
    shutil.copy2(build / "cpp/CMakeFiles/pea.dir" / name, destination)
    build_hashes[name] = digest(destination)
libraries = {str(p): digest(p) for p in (
    Path("/home/rx/GINAN/vcpkg_installed/linux/x64-linux/lib/libopenblas.a"),
    Path("/home/rx/GINAN/vcpkg_installed/linux/x64-linux/lib/liblapack.a"))}
probe_path = artifact / "one_epoch_probe/probe.json"
probe = json.loads(probe_path.read_text()) if probe_path.exists() else None
real_path = artifact / "one_epoch_probe/real_epoch_audit.json"
real = json.loads(real_path.read_text()) if real_path.exists() else None
product_comparison=json.loads((artifact/"product_epoch_comparison.json").read_text())
assert product_comparison["identical_product_keys"] and product_comparison["identical_covariance_keys"]
assert not any(product_comparison["publication_flag_changes"].values())
float_line=next(line for line in (artifact/"real_full_float_compare.log").read_text().splitlines()
                if line.startswith("REAL_FLOAT_FULL_COMPARE "))
full_float=dict(re.findall(r"(\w+)=(\S+)",float_line))
assert float(full_float["covariance_relative_frobenius"])<2e-11
assert full_float["state_order_identical"]=="1"
timing_comparison={}
if real:
    assert digest(artifact/"pea_after")==probe["binary_sha256"]
    old=json.loads((artifact.parent/"ar_full_process_breakdown_20261019/timing.json").read_text())
    rows=old["per_epoch"][real["committed_epoch_times"][0]]
    for name in ("AR_TOTAL","R51_WL_BLOCK_SEARCH","R51_L1_BLOCK_SEARCH",
                 "R51_PRODUCT_CLOSURE_FINAL","KF_ROOT_FACTOR","R51_WL_PHYSICAL_COMPATIBILITY",
                 "R51_L1_PHYSICAL_COMPATIBILITY","R48_OFFICIAL_PRODUCTS","R48_HISTORY_SEARCH","R48_PRODUCT_IMAGE"):
        selected=[r for r in rows if r["phase"]==name]
        before=sum(r["inclusive_ms"] for r in selected)/1000 if selected else None
        after=real["phases"][name]["inclusive_seconds"] if name in real["phases"] else None
        timing_comparison[name]=dict(old_seconds=before,new_seconds=after)
(artifact/"real_timing_version_comparison.json").write_text(json.dumps({
    "scope":"historical version comparison, not strict paired timing; inclusive phases overlap",
    "phases":timing_comparison},indent=2))
manifest = dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                source_root=str(root), source_sha256=source_hashes,
                binary_before_sha256=digest(artifact / "pea_before"),
                binary_after_sha256=digest(artifact / "pea_after"),
                source_patch_sha256=digest(artifact / "ar_optimization.patch"),
                build_input_sha256=build_hashes, math_library_sha256=libraries,
                ctest_passed=38, ctest_excluded="10 long cycle-slip Monte Carlo tests; deterministic detector tests included; detector unchanged",
                ratio=2, joint_nis=False, epoch_only=True,
                cpu_logical=12, cpu_physical=6, ar_omp_budget=12, ar_blas_budget=12,
                eigensolver_blas_cap=6, benchmark=results, one_epoch_probe=probe,
                final_v2_kernel_confirmation=json.loads((artifact / "v2_confirmation.json").read_text()),
                numerical_ar_object_sha256=digest(build / "cpp/CMakeFiles/pea.dir/pea/ppp_ambres.cpp.o"),
                numerical_ar_object_unchanged_after_shutdown_fix=(digest(build / "cpp/CMakeFiles/pea.dir/pea/ppp_ambres.cpp.o") ==
                    digest(Path("/home/rx/GINAN/ar_optimization_20261019/optimized_v2/product_compiler.o"))),
                eigen_tridiagonal_blas_cap=1,
                condition_regression_guard={"medians_seconds": condition_medians,
                    "records": condition_profiles,
                    "note": "Avoids this iteration's all-12 regression; not a 20x gain over original BLAS=4."},
                real_epoch_audit=real,
                full_real_float_covariance_comparison=full_float,
                full_product_joint_covariance_comparison=product_comparison,
                real_timing_version_comparison=timing_comparison,
                superseded_evidence={
                    "paired/": "Initial parent retained large covariance arrays; RSS contaminated. Only paired_final/ is used.",
                    "baseline_closure*.log except baseline_closure_structural.log":
                        "Early startup-tree fixtures produced an empty target lattice; replaced with explicit distinct controlled trees."},
                integer_audit_unchanged_from_before=(digest(root / "src/cpp/common/zhangIntegerAudit.hpp") ==
                    digest(artifact / "before/zhangIntegerAudit.hpp")))
(artifact / "MANIFEST.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False))
lines = [
    "R51 AR 三项优化与验证（2026-10-19）",
    "",
    "产品闭合：具有独占 ±1 列的循环基直接恢复并逐项核验整数见证；其他情况回退原精确后端。命名关系批量恢复，不再逐目标重复分解。",
    "FLOAT 根：当前 AR 调用内保留不可变 FLOAT 根；严格核对 P 全部元素、状态顺序及历元。大矩阵使用 LAPACK DSYEVD，仍保留原 PSD 门禁。",
    "物理兼容性：通过变量—约束邻接图只核验受新增行影响的连通部分；保留整除、同余、零行矛盾和提交票据门禁。",
    "线程 v2：AR 独立目标和大型 GEMM 默认 12 线程；LAPACK 根分解上限 6；已管理的 WL/L1/KF 主要 Eigen 分解单独限制 BLAS=1。并非所有遗留 Eigen 调用均已改造，未宣称整个AR全核满载。退出后恢复设置，精确预算取消保持有效。",
    "",
    "生产内核新旧对照：每项 5 次串行交错运行，以下为中位秒数。不是完整产品闭合或整历元的加速结论。",
]
for mode in ("closure", "root", "compat"):
    value = results["summary"][mode]
    lines.append(f"{mode}: {value['before']['median_seconds']:.6f} -> {value['after']['median_seconds']:.6f} s; "
                 f"{value['speedup']:.2f}x; peak RSS "
                 f"{value['before']['maximum_peak_rss_mib']:.2f} -> "
                 f"{value['after']['maximum_peak_rss_mib']:.2f} MiB")
lines.extend([
    "",
    "闭合受控图：32 站×20 星、589 个循环坐标、产品秩 19；整数基、见证、搜索坐标/顺序与偏移逐字节一致。",
    "FLOAT 根受控矩阵：1800×1800；比较全部 3,240,000 个协方差元素，不比较非唯一的特征向量基。",
    "兼容性受控矩阵：768 个商变量、32 个独立块后追加矛盾；这类输入加速不能外推到所有变量高度耦合的情况。",
    "回归：38/38 通过，包括真实生产对象链接的 FLOAT/AR、图运行时、周跳 QC、周跳确定性单测、整数/连续性门禁和真实 SIGINT 单测。未重跑 10 项长耗时周跳蒙特卡罗测试。",
    "缓存代价：额外 FLOAT P 和根均为当前调用内有界矩阵，约 2*n*n*8 字节；不按历史长度增长。",
    "ratio=2、整网联合 NIS 关闭、每历元独立 AR 不变。未提交/推送 Git，未启动一小时解算。",
    "首轮 paired/ 的内存记录受计时父进程大数组影响，未采用；正式结果只来自 paired_final/。早期空目标格夹具日志也保留，但不计为通过的性能测试。",
    "v1 单历元探针因全局 BLAS=12 导致小块条件化退化而主动停止（退出137、无提交）；该失败未当作测试通过。v2 使用按数值内核分配线程的策略重新验收。",
    "v2 完成首个 AR 历元（516.077秒）并提交产品，但旧 SIGINT 写 max_epochs 导致 E29 配置溯源失败（退出1）。已改为独立原子停止标志，保留全部溯源检查；主动停止明确返回130并标记一小时未完成，不把未到的检查点伪装为完成。",
    f"1175维/3次32行条件化：全局12线程退化版中位 {condition_medians['unmanaged_12']:.6f} 秒，v2为 {condition_medians['adaptive_12_1']:.6f} 秒。原4线程单次对照0.858385秒，因此这是消除本轮引入的退化，不能算成相对原版20倍新增收益。",
])
if probe:
    lines.extend(["", "单真实历元验收记录：" + str(probe_path),
                  "提交 bundle 数（不等同于历元数）：" + str(len(probe["committed_epochs"])),
                  "退出码：" + str(probe["returncode"]) + "（主动停止的预期值130，不是完整一小时成功退出0）"])
if real:
    lines.extend([
        "实际提交历元：" + ", ".join(real["committed_epoch_times"]),
        f"AR_TOTAL：{real['phases']['AR_TOTAL']['inclusive_seconds']:.6f} 秒",
        f"整进程峰值 RSS（含输入加载与 FLOAT）：{real['peak_rss_mib']:.2f} MiB",
        "最多双频 AR 可用卫星：" + str(max(len(b['dual_frequency_ar_satellites']) for b in real['bundles'])),
        "bundle SHA256、字节数及 AR 发布主数值有限性核验通过。单历元不证明一小时连续性、长期内存稳定性或产品外部精度。",
        "已生成当历元新检查点，并通过生产读取器的包校验、完整KF和全部区段解码；未进行恢复后下一历元重放。原始01:59:30检查点未修改。",
        "真实历元计时不是新旧二进制严格成对重放，不用来宣称整小时加速倍数。",
    ])
    lines.extend(["", "同一首历元与旧冻结记录的版本对照（分项有包含关系，不相加）："])
    for name,values in timing_comparison.items():
        lines.append(f"{name}: {values['old_seconds']:.6f} -> {values['new_seconds']:.6f} 秒")
    lines.extend([
        "HISTORY 路径在本次 epoch-only 模式下接收同历元已接受整数，不是跨历元整数复用；日志 prior_integer_reuse=0 已核对。",
        "局部兼容性只拆分整数方程约束图，不把浮点相关性当作零；完整交叉协方差继续传播。",
        "仍未解决：产品闭合占本次 AR 超过一半。官方产品子路径比旧首历元记录更慢，不能称全部子步骤都加速；精确整数像构造和未细分的闭合/总控开销是后续重点。",
        "完整真实FLOAT比较：4755维、22610025个协方差元素，状态顺序一致，协方差相对Frobenius误差="+full_float["covariance_relative_frobenius"]+"；绝对差为混合原生参数单位，不可标成统一m²。",
        "产品比较：232行记录及15312项完整产品联合协方差；发布/分段标志无变化，改正量最大差="+
            str(product_comparison["maximum_absolute_differences"]["correction_m"])+" m；联合协方差最大差="+
            str(product_comparison["maximum_covariance_difference_m2"])+" m²。",
    ])
(artifact / "FINAL_REPORT.txt").write_text("\n".join(lines) + "\n")
print(json.dumps({"report": str(artifact / "FINAL_REPORT.txt"),
                  "binary_sha256": manifest["binary_after_sha256"]}, indent=2))
