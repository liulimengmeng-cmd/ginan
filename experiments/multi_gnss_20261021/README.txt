GPS + Galileo + BDS-3 输入审计及 A-G 阶段驱动
2026-10-21

实际信号：GPS L1C/L2W；Galileo L1C/L5Q；BDS-3 L1P/L5P。
约束不变：ratio=2；所有整数/产品 NIS 数值计算关闭；每历元独立 AR。
本目录脚本不修改生产源码，不重写原观测/轨钟/旧输出，不下载数据。

1. audit_inputs.py
读取原 180 站服务端清单、全部 RINEX 头；仅对候选站读取 00:00-00:05 小窗口。
生成 input_readiness.json、bds3_whitelist.json、INPUT_READINESS_REPORT.txt。
头有某信号不表示整日有有效观测；首5分钟有非零观测也不是残差/QC/精度验收。
BDS-3 白名单由本地 SINEX IDENTIFIER + 当日有效 PRN 指派，再交 WUM SP3/CLK，
记录来源文件/原始行哈希/有效期，不采用 PRN 数值阈值猜测代际。

2. run_stages.py
默认只 prepare，必须显式 --execute 才运行，并且必须给 --build-approval JSON。
每个 --run-root 必须是新目录，不允许覆盖任何已有实验。
JSON 格式的 overlay.yaml 是标准合法 YAML，使用 yaml-cpp 原生解析器。

A：原180 GPS检查点回归。六个 pinned YAML 不变。在私有 mount namespace 内将
原输出/输出checkpoint目录绑定到新目录，baseline/candidate 顺序运行，到第一个
原子COMMIT后停止新启动的进程组。比较完整 products/covariance/certificate。
需要 root 仅用于 unshare/mount；不是更改宿主全局挂载，也不删除旧输出。
SIGINT 退出是受控单历元，不得声称完整小时自然结束。

B_E/B_C：Galileo-only / BDS3-only 冷启动 FLOAT，00:00-02:05。
C_GE/C_GC：G+E / G+C 冷启动 FLOAT，00:00-02:05。
D：G+E+C 冷启动 FLOAT，00:00-02:05。
E：G+E+C 从00:00 FLOAT预热，02:00开始短 AR，到02:05。
F：独立静态用户验证的明确阻断清单；不是已实现或已运行的用户实验。
   必须有真正多系统服务产品、独立参考坐标、同信号/同基准配置、完整跨系统协方差。
G：通过短 AR 后，冷启动预热再持续一小时 AR；E11故障隔离/换树/资源长期趋势
   必须另有受控试验证据，不能仅据这个运行的正常返回认定通过。

B-G 一律不导入原01:59:30 GPS检查点，不使用 ZHANG_E29_RESTORE_BUNDLE 或 FLOAT_REUSE。
原检查点明确标为 E29_GPS_L1C_L2W_ZHANG_FULL_RANK，不含E/C状态、图和协方差。

执行批准文件必需字段：
{
  "approved_for_execution": true,
  "build_passed": true,
  "mandatory_tests_passed": true,
  "binary_sha256": "本次实际候选程序 SHA256",
  "multi_system_clock_model_review_passed": true,
  "precise_clock_datum_review_passed": true,
  "short_ar_acceptance_passed": false
}
最后一个字段在G阶段必须为true。缺失/false/二进制hash变化一律不执行。
禁止为凑通过率将尚未检查的字段写true；配置生成/编译/测试通过各不等同科学验收。

本机编译环境只读核验：
WSL Ubuntu-22.04-G，Ubuntu22.04.5，12逻辑CPU、约19GiB RAM、16GiB swap。
build: /home/rx/GINAN/build-r51-handoff-20260930
CMAKE_HOME_DIRECTORY: 当前repo/src；Release/O3/OpenMP/OpenBLAS。
默认产物: repo/bin/pea（不是 build/pea 或 build/bin/pea）。
旧冻结verified_v2: /home/rx/GINAN/bridge_continuity_20261020/verified_v2/pea。
构建: cmake --build /home/rx/GINAN/build-r51-handoff-20260930 --target pea -j2
测试: ctest --test-dir /home/rx/GINAN/build-r51-handoff-20260930 --output-on-failure
编译并行2是按约19GiB内存防峰值，不表示求解阶段只使用2核。

示例（prepare，不运行）：
python3 run_stages.py --stage E --audit /path/input_readiness.json \
  --run-root /home/rx/GINAN/local_timing/multi_gnss_20261021_E_01 \
  --binary /mnt/c/Users/rx/Documents/GINAN/ginan-r51-ratio-only-nis-skip-20260923/bin/pea

后续在 root 明确完成编译/阻断测试/模型核对后，再用新目录和 --execute --build-approval。
当前目录存在不意味着任何真实多系统试验已完成。
