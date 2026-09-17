# R51 ratio-only 实验模式

本次按用户明确要求，将整数候选的统计接受方式设为 ratio >= 3。通过环境变量 `ZHANG_R51_RATIO_ONLY=1` 显式开启；正式实验同时开启 `ZHANG_R51_ENABLE=1`，YAML 使用 LAMBDA、solution_ratio_threshold=3。默认不开启时保留原有判据。

## 实施范围

- 普通 LAMBDA/LAMBDA_ALT：不再以 bootstrap 成功率决定可搜索部分秩，不再由局部 NIS 拒绝候选；继续使用现有 LAMBDA 搜索及完整、不同的最优二候选 ratio 判据。
- 网络 R51 WL/L1 分块：Perr/风险预算耗尽、根空间 NIS 超限不再作为整数接受拒绝条件。
- 历史网络约束、产品 pending、gauge/ledger：统计 NIS/Perr 只保留诊断，不因这些门限重新筛掉历史整数。
- 产品搜索、桥接、历史选择、联合约束提交和最终产品写出：同步解除非 ratio 统计拒绝，避免搜索通过后被后续统计门限重新拒绝。
- 保留协方差数值合法性、完整二候选搜索、精确整数相容性、物理弧/坐标运输、证据来源及测量级 QC。网络分块、候选顺序、搜索次数上限不变。
- ratio-only 接受证明标记 `R51_RATIO_ONLY_NO_FFR_GUARANTEE`，以 1 表示单次无信息失败上界，不伪造低 Perr；运行日志 `risk_bound_valid=0`。此模式的固定不能声称满足原先校准的失败概率门限。
- 新检查点将 ratio-only、R51 开关、AR 开始时刻、融合权重及数值线程配置纳入身份，恢复不得静默切换。

## 可审计日志

`R51_VALIDATION_POLICY ratio_only=1 ratio_threshold=3 bootstrap_gate=0 nis_gate=0 perr_gate=0 risk_budget_gate=0 risk_bound_valid=0`

`ZHANG_RATIO_GATE` 记录实际 ratio、阈值、完整搜索及接受/拒绝。`ZHANG_R51_BLOCK_ACCEPT` 记录 root_nis_veto=0、bootstrap_veto=0。

注意旧 `ZHANG_LAYERED_AR_RESULT ... fixed=0` 是 R51 已跳过的旧搜索分支，R51 实际固定应读取 `ZHANG_R51_BLOCK_STAGE` 和最终提交/历史整数记录，不应据旧日志误报为全网零固定。

## 验证范围

新增生产对象回归：高 ratio 但 bootstrap 失败、高 ratio 但 NIS 超限、ratio 在 3 上下边界以及等距候选；新模式应仅让前两类越过统计拒绝。保留旧模式回归，以及非法协方差、矛盾整数、缺失证明来源拒绝测试。

测试与真实短段结果见本次独立证据目录 `C:/Users/rx/Documents/GINAN/r51_ratio_only_20260917`。不能以编译通过代替测试完成。

## 正式实验

新二进制和模式不能直接使用旧二进制绑定的检查点。本次重新从 2024-07-17 00:00:00 GPST 启动：00:00:00–01:59:30 为 FLOAT；02:00:00–02:59:30 为 ratio-only AR，共360历元，180站、30秒间隔。使用新冻结目录与新输出目录，保留之前所有运行结果。
