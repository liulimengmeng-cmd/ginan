# R51 ratio-only 实验模式

本次按用户明确要求，将整数候选的统计接受方式设为 ratio >= 3。通过环境变量 `ZHANG_R51_RATIO_ONLY=1` 显式开启；正式实验同时开启 `ZHANG_R51_ENABLE=1`，YAML 使用 LAMBDA、solution_ratio_threshold=3。默认不开启时保留原有判据。

## 实施范围

- 普通 LAMBDA：不再以 bootstrap 成功率决定可搜索部分秩；局部 NIS 数值和卡方阈值均不计算。继续使用完整、不同的最优二候选距离和 ratio 判据；候选距离是 ratio 所必需的量，不能把它重新标作已评估的 NIS。
- 网络 R51 WL/L1 分块：Perr/风险预算耗尽不再作为整数接受拒绝条件；根空间联合 NIS 的矩阵构造、特征分解、数值和阈值均跳过。
- 历史网络约束、产品 pending、gauge/ledger：整数候选／产品层 NIS 数值与阈值不再计算，字段为 NaN／NA，不得解读为零或通过。产品继续发布，但不带 NIS 或校准失败概率认证；必须具备上游 ratio 决策来源、精确整数一致性、数值有效性及写出授权。没有上游 ratio 证明的独立标量舍入候选不得凭关闭统计门槛取得产品授权。
- 产品搜索、桥接、历史选择、联合约束提交和最终产品写出：同步解除非 ratio 统计拒绝，避免搜索通过后被后续统计门限重新拒绝。
- 保留协方差数值合法性、确定性零空间冲突检查、完整二候选搜索、精确整数相容性、物理弧/坐标运输、证据来源及测量级 QC。测量级 QC 不属于本次整数／产品 NIS 开关；不得误称为全部统计检验关闭。NIS 驱动的历史候选排序改为确定性索引次序，冲突分支在无独立证据时不强制切换。
- ratio-only 接受证明标记 `R51_RATIO_ONLY_NO_FFR_GUARANTEE`，以 1 表示单次无信息失败上界，不伪造低 Perr；运行日志 `risk_bound_valid=0`。此模式的固定不能声称满足原先校准的失败概率门限。
- 新检查点将 ratio-only、R51 开关、AR 开始时刻、融合权重及数值线程配置纳入身份，恢复不得静默切换。

## 可审计日志

`R51_VALIDATION_POLICY ratio_only=1 ratio_threshold=3 bootstrap_gate=0 nis_gate=0 perr_gate=0 risk_budget_gate=0 risk_bound_valid=0`

`ZHANG_RATIO_GATE` 记录实际 ratio、阈值、完整搜索及接受/拒绝。`ZHANG_R51_ROOT_NIS_POLICY evaluation=SKIPPED_RATIO_ONLY` 和 `ZHANG_R51_BLOCK_ACCEPT root_nis=NA root_threshold=NA root_nis_evaluation=SKIPPED_RATIO_ONLY` 记录根空间 NIS 未执行。其他 NIS 字段的 NaN 也只表示未评估，不能用于产品质量声明。

注意旧 `ZHANG_LAYERED_AR_RESULT ... fixed=0` 是 R51 已跳过的旧搜索分支，R51 实际固定应读取 `ZHANG_R51_BLOCK_STAGE` 和最终提交/历史整数记录，不应据旧日志误报为全网零固定。

## 验证范围

新增生产对象回归：高 ratio 但 bootstrap 失败、高 ratio 但旧模式 NIS 超限、ratio 在 3 上下边界以及等距候选；ratio-only 模式不执行 NIS 数值计算。保留旧模式回归，以及非法协方差、矛盾整数、缺失证明来源拒绝测试。固定产品必须另行检验可用性与独立相位一致性，不能仅凭 ratio=3 断言整数正确。

测试与真实短段结果见本次独立证据目录 `C:/Users/rx/Documents/GINAN/r51_ratio_only_20260917`。不能以编译通过代替测试完成。

## 正式实验

按用户要求复用 2024-07-17 01:59:30 GPST 的纯 FLOAT 检查点，从 02:00:00 开启 ratio-only AR，计算至 02:59:30（120 个新历元）。前 240 个 FLOAT 历元沿用原始结果。

显式 `ZHANG_R51_REUSE_FLOAT_20240717=1` 仅授权一个经过审计的旧 bundle SHA256、旧二进制 SHA256、完整旧配置 SHA256、输入清单 SHA256 组合。当前配置由当前 argv、原文件及有效配置重新计算，仅去掉新版新增的环境身份字段；其余必须逐字相同。线程、融合权重、AR 开始时间也必须符合原运行。02:00:00 及之后的检查点由于 bundle 哈希不同被拒绝。旧读取、ABI、模块、游标、KF 位级往返和两阶段提交校验照常执行；后续捕获记录真实的新二进制与模式身份。原 bundle 不修改，不伪造新二进制来源。

旧输出和旧检查点目录完整归档后在原配置路径生成新结果，避免改变任何配置内容。证据目录保存迁移来源、原 FLOAT 日志及源码差异；本入口是此次指定检查点的受限兼容入口，不是通用跨版本恢复功能。

## 2026-09-23 NIS 数值关闭验证状态

隔离分支的 `r51_ratio_policy_retains_structural_rejections` 已通过 41 项断言；全量 `zhang_full_rank_tests` 为 448/449 测试通过。唯一失败是未修改的 E29 检查点测试严格期待 `CHECKPOINT_PROVENANCE_MISMATCH`，实际返回带 `CONFIG_SHA256` 细节的同类错误。尚未用本次新二进制完成真实历元的固定产品发布、独立相位一致性检查或端到端性能对照；因此不能声称新产品质量或速度已被实测证明，也不能直接替换正在运行的服务端任务。
