# 2024 PPP-AR 反馈信息增益与重启对照实验

## 结论

本轮完成了 guarded selected-subset PPP-AR 反馈路径的代码诊断、全天运行和同绝对时窗冷重启对照。工程路径通过结构与提交一致性检查，但科学收益没有通过：反馈对耦合 STEC 协方差迹的一步中位降幅仅为 0.87%–0.94%，远低于本轮采用的 20% 可行性筛选值；收敛后坐标相对 FLOAT 的变化只有亚毫米的千分量级，且连续解与冷重启方向不一致。

因此，本结果只能证明整数伪观测确实进入滤波器，不能证明整数正确，也不能证明 PPP-AR 已改善定位或 STEC。当前不应通过降低 LAMBDA 成功率或 ratio 门槛来提高表面固定率，也不应进入基于“fixed STEC”的空间建模阶段。

## 代码与版本

- 工作树：`D:\tec\ginan-main-full-rank-ar`
- 分支：`codex/full-rank-pppar-datum`，由 `origin/main` 隔离工作树继续开发；未修改原 Zhang 脏工作树。
- 诊断实现提交：`f82bf67` (`audit: trace PPP-AR feedback information gain`)
- PEA SHA-256：`E2F234BD2C3712CA09762DF3C3DD90DCCC8CA994F24A6EA1421641731593DEBD`
- 新诊断记录：每次反馈前后的联合 NIS、状态更新、协方差迹变化，以及 COORDINATE、STEC、AMBIGUITY、TROPOSPHERE、CLOCK 等状态块的全块和实际耦合子集信息增益。
- 审计脚本支持 `--epoch-start/--epoch-end`，用于按绝对时段比较连续解和重启解；联合 NIS 目前仅诊断，没有虚构或启用未经标定的拒绝阈值。

测试与构建：Python 审计单元测试 11 项通过；`gnss_ambiguity_diagnostics_tests` 通过；PEA 构建通过。GitVersion 无法解析 Windows worktree 的 `.git` 指向，因此实验绑定源码提交和二进制 SHA-256，而不是使用其 `untagged--git-error` 版本串。

## 使用的数据

- 日期：2024-07-17（GPS week 2323，day 199），30 s 采样。
- 测站：DYNG、GRAZ、MARS、MATE、NICO、SOFI，共 6 站。
- 观测：对应 6 个 `2024199` 全天 RINEX；冷重启使用从 12:00 切出的 `data_restart_2024199_1200`，绝对评估窗为 12:00–24:00。
- 系统与频点：仅 GPS，公共 L1C/L2W 双频图；G01 因所选 WUM SP3、CLK 和 OSB 不完整而排除。缺少第二频率的观测不会进入双频候选。
- 轨道：`WUM0MGXFIN_20241980000_01D_05M_ORB.SP3`、`...20241990000...`、`...20242000000...`。
- 钟差：`WUM0MGXFIN_20241990000_01D_30S_CLK.CLK`。
- 偏差：`WUM0MGXFIN_20241990000_01D_01D_OSB.BIA`。这是外部 WUM Bias-SINEX/OSB，不是 Zhang 内部 `HOU_OSB_LIKE` 产品。
- ERP：`WUM0MGXFIN_20241990000_01D_01D_ERP.ERP`。
- 天线：`igs20_2303.atx`；坐标外检：`IGS0OPSSNX_20241960000_07D_07D_CRD.SNX`。
- 其他模型：BRDC day 199、GPT2、DE436、IGRF14、海潮/大气潮/极潮加载表。

## 实验设计

| 解算 | 绝对时段 | 历元数 | 用途 |
|---|---:|---:|---|
| Exp2g 连续 PPP-AR | 00:00–24:00 | 2880 | 全天反馈与信息增益 |
| Exp2g 同窗切片 | 12:00–24:00 | 1440 | 与冷重启严格同窗比较 |
| Exp2h 冷重启 PPP-AR | 12:00–24:00 | 1440 | 清除前 12 h 状态历史 |
| FLOAT 配对对照 | 相同绝对时段 | 同上 | 持久无反馈轨迹与坐标对照 |

AR 配置保持 `success_rate_threshold=0.9999`、`solution_ratio_threshold=3`，没有为了提高发生率放宽门槛。所谓“反馈历元”是至少一个通过现有门控的整数子集已调用滤波更新并返回的历元；它不是所有可见模糊度均正确固定的历元。

## 固定/反馈发生率与秩覆盖

| 指标 | 全天连续 | 连续解 12:00–24:00 | 12:00 冷重启 |
|---|---:|---:|---:|
| 有反馈历元 | 1483/2880 (51.49%) | 987/1440 (68.54%) | 492/1440 (34.17%) |
| 有反馈站历元 | 3259/17280 (18.86%) | 2145/8640 (24.83%) | 1027/8640 (11.89%) |
| 提交秩/合格秩历元和 | 10.49% | 14.05% | 7.05% |
| 提交整数伪观测行 | 22530 | — | 7642 |

因此，“目前固定率”不能只报 51.49% 或 68.54%。较不误导的主指标是站历元反馈发生率和提交秩覆盖：严格同窗连续解分别为 24.83% 和 14.05%，冷重启分别为 11.89% 和 7.05%。这仍不是经真值验证的 full ambiguity fix rate。

同窗中，连续解各站首次反馈在窗起点后 77.5–522.0 min；冷重启为 448.5–643.0 min。连续解的能力约有一半依赖 12:00 前的滤波历史。共享候选的规范化整数右端 5506/5506 一致，说明重启差异不是共享行整数值相互矛盾，而是候选可用性和状态历史不同。

## 反馈信息增益

| 同窗中位相对协方差迹降幅 | 连续解 | 冷重启 |
|---|---:|---:|
| STEC 全状态块 | 0.0741% | 0.1159% |
| STEC 实际耦合子集 | 0.8686% | 0.9393% |
| 坐标全状态块 | 1.2196% | 2.4801% |
| 坐标实际耦合子集 | 3.7840% | 8.0704% |

解析一步 shadow 与实际滤波调用的最大状态更新误差约 `1.0e-9`，协方差误差约 `1e-16` 量级，说明信息增益计算与实际更新一致。联合 NIS/行中位数为连续解 1.128、冷重启 1.436；由于没有预注册阈值，本轮不把它当拒绝门。

耦合 STEC 一步增益不到 1%，没有通过 20% 形式可行性筛选。坐标耦合子集的形式增益较大，但外部精度没有同步改善，说明形式协方差收缩不能代替外部准确度验证。

## 收敛后坐标精度

以下“收敛后”严格指预设 2 h burn-in，不是由数据自动检测出的收敛时刻。公平比较窗为 14:00–24:00，共 7200 个站历元。SINEX 坐标直接使用，未推断站速传播。

| 解算 | 对 SINEX 3D RMS | 相对配对 FLOAT |
|---|---:|---:|
| 连续 PPP-AR | 7.210 cm | 改善 0.090 mm (0.125%) |
| 连续 FLOAT | 7.219 cm | 基线 |
| 冷重启 PPP-AR | 9.195 cm | 恶化 0.072 mm (0.078%) |
| 冷重启 FLOAT | 9.188 cm | 基线 |

连续 PPP-AR 与冷重启 PPP-AR 在相同绝对时刻的 3D 差异 RMS 为 8.682 cm；对应 FLOAT–FLOAT 差异为 8.676 cm。两者几乎相同，所以重启敏感性主要来自滤波初始化/历史，而不是整数反馈。反馈对外部 3D RMS 的效应小于 0.1 mm，方向还随是否冷启动改变，不构成可重复的精度收益。

## 失败位置与下一步

第一失败位置主要集中在宽巷阶段：大量站历元的宽巷只形成部分 accepted subset，未覆盖目标秩；冷重启另有 1205 个站历元落在宽巷成功率不足，而同窗连续解只有 24 个。SECOND_D2 的成功率/ratio 拒绝是次要部分。

据此，下一轮不应继续降低整数检验阈值。优先代码方向应是把持久 shadow 分支和状态生命周期拆开，做同一时刻的状态块重置消融（只重置 ambiguity、STEC、clock、ZTD、coordinate），并在新的 2024 日期上预注册检验：耦合 STEC 形式增益至少 20%，外部/交叉站 dSTEC 指标至少改善 10%。任一未通过，就停止把该反馈用于 STEC 空间研究。

## 证据文件

- `experiment_2g_2880e_feedback_information_audit.json`
- `experiment_2g_1200_2400_feedback_information_audit.json`
- `experiment_2h_1440e_restart_1200_feedback_information_audit.json`
- `experiment_2g_2880e_dual_frequency_feedback_audit.json`
- `experiment_2h_1440e_restart_1200_dual_frequency_feedback_audit.json`
- `experiment_2g_2880e_feedback_information_coordinate_integrity_summary.json`
- `experiment_2h_1440e_restart_1200_feedback_information_coordinate_integrity_summary.json`
- `experiment_2gh_same_window_coordinate_consistency_summary.json`
- `experiment_2gh_float_same_window_coordinate_consistency_summary.json`
- `experiment_2gh_restart_candidate_consistency_audit.json`

运行输出：

- `/home/rx/GINAN/inputData/outputs/exp2g_feedback_info_2880e_20260901_f82bf67_bine2f234bd`
- `/home/rx/GINAN/inputData/outputs/exp2h_restart_1200_feedback_info_1440e_20260901_f82bf67_bine2f234bd`
