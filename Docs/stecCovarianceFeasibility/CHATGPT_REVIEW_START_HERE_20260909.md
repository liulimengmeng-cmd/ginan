# GINAN PPP-AR 当前问题与证据：供 ChatGPT 独立分析

本文件汇总截至 2026-09-09 的工作。请先读此文件，再读报告和原始矩阵，不要把历次讨论中的猜测当成结论。当前未建立可靠固定率，也没有证明电离层解错误。

## 研究对象

- 数据日期：2024-07-17；站点包括 DYNG、GRAZ、MARS、MATE、NICO、SOFI。
- GINAN 未组合解算，估计 STEC；候选按站构造，使用星间单差及双频整数关系，再映射回原始状态。
- CODE 产品作为外部输入；伪距保持 CODE 码 OSB。用户允许未来用 CODE DCB 替代对应码改正，但本轮没有切换，也没有重复叠加。
- 原始 CODE GPS L1C/L2W 相位 OSB sigma 字段为零。人工 3 mm / 10 mm 是敏感性实验，不是 CODE 发布的精度。
- 非零相位 OSB 先验实现为日持久的卫星/信号偏差状态，由站间共享，过程噪声为零。
- 本轮配置 once_per_epoch=true、fix_and_hold=false，整数反馈作用于当历元输出副本，不作为后续主滤波状态继续传播。
- 基线代码提交为 `91c1ceae9ecf1192419d1db83841d1751ef7cb37`；本分支增加快照导出及显式启用的实验联合拒绝检查。

## 已经观察到什么

### A. 人工 10 mm 实验存在跨站候选冲突

20:51:00 GPST，冻结同一个先验、完整 P 和同一组候选：

| 候选 | 行数 | NIS/行 |
|---|---:|---:|
| DYNG | 8 | 2.8239 |
| GRAZ | 6 | 5.0694 |
| NICO | 6 | 3.3913 |
| DYNG+NICO | 14 | 2.6995 |
| DYNG+GRAZ | 14 | 664.7152 |
| GRAZ+NICO | 12 | 719.3033 |
| 全部 | 20 | 592.7271 |

仅为诊断而将跨站协方差块置零后，NIS/行降为 3.6678。没有在估计器中删除相关性。

只将零基索引 12 的 GRAZ L2 单差整数右端从 -241 改成 -240：GRAZ 单站 NIS 从 30.4166 升到 81.9548，网络 NIS 从 11854.5421 降到 106.7847。这证明原逐站拼接向量不是该网络目标的最优整数；没有证明替代值为真值或全局最优，也未将其应用于解算。

### B. 原始 CODE 对照没有同样的跨站暴增，但名义检验几乎全部拒绝

原始 CODE 对照 20:51:00 的 60 行候选 NIS 为 612.9819；删除跨站块前后几乎相同。不能将人工 10 mm 的机制直接泛化到原始产品。

新增实验检查在更新及后验降权之前计算 S=HPHᵀ+R、NIS=vᵀS⁻¹v。复用现有滤波器 sigma=4 的双尾换算，按行数计算名义卡方门限。它没有进行搜索后分布/实际误拒率标定。

| 2504 历元重放 | 有候选历元 | 通过名义检查 | 拒绝 |
|---|---:|---:|---:|
| 人工 10 mm | 2175 | 1465 | 710 |
| 原始 CODE | 2433 | 1 | 2432 |

**被拒绝的是候选，不是产品。这些计数不是正确/错误固定数。**还不能确定原始产品对照的问题来自整数候选、协方差物理标定、观测模型或不适用的名义门限。实验检查默认关闭，需 `GINAN_AR_JOINT_NIS_GATE=1` 显式启用。

### C. 电离层影响的证据边界

原 10 mm AR 异常事件经后验降权后，GRAZ/G14 STEC 相对更新前变化约 -0.211 TECU；这证明输出改变，不证明精度变差。

新增拒绝检查后，10 mm 实验全部 710 个拒绝历元的 30400 个 STEC 值及方差与独立 FLOAT 完全一致；1465 个通过历元的 64295 个 STEC 值及方差与原 AR 完全一致。拒绝快照的完整 STEC 协方差也与先验一致。没有独立电离层真值。

## 推荐阅读顺序

1. [同历元诊断报告](PPPAR_SAME_EPOCH_JOINT_NIS_DIAGNOSIS_20260909.md)
2. [联合拒绝检查实现与验证](PPPAR_JOINT_GATE_VALIDATION_20260909.md)
3. [人工 sigma 全天实验](PPPAR_CODE_PHASE_OSB_NONZERO_SIGMA_RESULTS_20260909.md)
4. [原始 CODE 产品对照](PPPAR_CODE_RAPID_PRODUCT_CONTROL_2024_RESULTS_20260909.md)
5. [关键整数反例](evidence_20260909/experiment_6/joint_integer_counterexample.json)
6. [原始 CODE 拒绝验证](evidence_20260909/experiment_7/0mm_gate_verification.json) 与 [10 mm 拒绝验证](evidence_20260909/experiment_7/10mm_gate_verification.json)

源码入口：[候选生成和反馈](../../src/cpp/pea/ppp_ambres.cpp)、[联合检查函数](../../src/cpp/pea/pppArJointGate.hpp)、[离线矩阵分析](../../scripts/audit_pppar_same_epoch_joint_nis.py)、[C++ 检查测试](../../scripts/test_pppar_joint_gate.cpp)、[重放验证](../../scripts/verify_pppar_joint_gate_replay.py)。

## 已上传的证据与未上传的内容

[evidence_20260909](evidence_20260909) 包含实验 5–7 的配置、运行/构建收据、审计结果、关键 TRACE 摘录，以及实验 6/7 各十份矩阵快照。快照包含 x、P、H、v、R、整数右端和状态键，可独立计算同历元候选的联合/边缘代价。

完整原始 RINEX、产品文件、全时段 TRACE、二进制及编译对象没有打包进 Git。它们的路径和哈希保留在收据中。复制的历史重放脚本含机器绝对路径和对已有构建对象的依赖，不是开箱即用的跨平台完整实验；离线快照分析不依赖这些原始文件。

本仓库根目录下可运行，例如（输出文件须不存在，依赖 NumPy）：

```sh
python scripts/audit_pppar_same_epoch_joint_nis.py \
  Docs/stecCovarianceFeasibility/evidence_20260909/experiment_6/snapshots_10mm/2024-07-17_20-51-00.snapshot \
  --output /tmp/ginan_joint_nis_review.json
```

历史报告中的部分 UTC 标签有误：TRACE 实际时间是 **GPST**，在该日期比 UTC 快 18 秒。请以实验 6/7 报告为准。旧收据及证据文件保持原样，原路径并不表示那些文件能通过 GitHub 访问；应使用此处复制的对应文件。证据文件清单及 SHA256 见 [MANIFEST.json](evidence_20260909/MANIFEST.json)。

## PRIDE 源码对照：已检查，但未做实测

核对官方 PRIDE 提交 `71fef9a7057156dfeae8fa659d37f972596263f7`。其常规流程是单站无电离层组合最小二乘，残差清理后进行宽巷/窄巷固定。窄巷搜索保留完整站内协方差；传统验证检查 ratio 和固定/浮点方差比型指标，当前模板默认 AI 验证；失败可尝试部分固定。它不是当前跨站共享 OSB 网络的同一模型，门限不能直接移植。

- [PRIDE 固定检验](https://github.com/PrideLab/PRIDE-PPPAR/blob/71fef9a7057156dfeae8fa659d37f972596263f7/src/arsig/fixamb_search.f90)
- [部分固定](https://github.com/PrideLab/PRIDE-PPPAR/blob/71fef9a7057156dfeae8fa659d37f972596263f7/src/arsig/candid_amb.f90)
- [OSB 读取](https://github.com/PrideLab/PRIDE-PPPAR/blob/71fef9a7057156dfeae8fa659d37f972596263f7/src/lib/read_bias.f90)：所查路径读偏差值和变化率，没有读取产品 formal sigma；不是共享待估 OSB 状态模型。

尚未完成：GINAN 单站 GRAZ 对照、PRIDE 同站同日实测、联合全局整数搜索、门限标定、独立 STEC 精度验证。此前建议属于下一步计划，不能写成已完成结果。

## 请 ChatGPT 重点分析

1. 候选搜索时的协方差与提交检验的 HPHᵀ+R 是否采用相同状态、单位、规范和统计尺度？
2. 人工 10 mm 的跨站反例能证明什么，不能证明什么？怎样设计保留相关性的部分/联合固定？
3. 对经过整数搜索的候选使用当前名义卡方门限是否合理？原始 CODE 几乎全拒绝，究竟还需要哪些证据才能定位原因？
4. 原始 CODE 单站候选的绝对代价、ratio、浮点残差和协方差应如何检查？不要仅建议增加 OSB sigma。
5. 如何做可解释的 GRAZ 单站 PRIDE/GINAN 对照，同时控制产品匹配和无电离层/未组合模型差异？

请将结论分为源码直接证据、数值直接证据和待验证假设，并给出最小可区分实验。不要把提交率当可靠固定率，不要把 STEC 变化当精度误差。
