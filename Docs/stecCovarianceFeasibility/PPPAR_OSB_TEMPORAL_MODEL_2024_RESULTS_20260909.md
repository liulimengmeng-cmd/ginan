# 2024 PPP-AR OSB 时间相关模型修复与重算结果

## 结论

旧实现把 WUM 日解 Bias-SINEX 的形式标准差在每个 30 s 历元中重复当作相互独立的观测噪声。对于同一个 86400 s 产品区间，这会让产品不确定度随历元数错误地按近似 `1/sqrt(N)` 缩小，从而高估整数信息量。

本次修复把卫星/信号相位 OSB 建模为产品区间内的持久 Kalman 状态：产品均值只用于状态初始化，产品方差只用于一次先验初始化，过程噪声为零，后续历元不再重复加入同一份相位 OSB 方差。为了避免共享卫星相位偏差状态在按接收机分块时丢失跨接收机协方差，修复实验关闭了 receiver chunking。

修复后的严格门限结果不是“PPP-AR 已成功”：全天 2880 个历元的反馈历元率为 **0/2880 = 0%**，六站历元反馈率为 **0/17280 = 0%**。所有 17247 个具有完整双频图的站历元都首先失败于宽巷成功率门限，另有 33 个站历元缺少完整双频图。AR 与 FLOAT 坐标逐历元完全一致，说明没有隐蔽整数反馈。

修复后的 14:00–24:00 三维坐标 RMS（相对 IGS 周解 SINEX）为 **0.0252266 m**；这是当前模型下的 FLOAT 精度，不是固定解精度。结果表明：以前出现的高固定表象主要依赖不正确的 OSB 时间白化，当前 WUM CLK+OSB 与六站数据在严格门限下不足以获得可信 PPP-AR。

## Pro 判断与源码判断的对应关系

Pro 对话指出连续参数基准在条件意义下可以闭合，接收机单差变换 `D` 的线性代数形式本身成立；首要风险是把日解 OSB 标准差当成 30 s 独立白噪声。源码检查确认该风险存在于外部 code/phase bias 的观测噪声构造中。本次没有修改 `D`、成功率门限 `0.9999` 或 ratio 门限 `3`，而是修复产品不确定度进入滤波器的时间语义。

## 代码修改

- `f387f18`: 为 code/phase bias 增加 `use_formal_sigma_as_observation_noise`，默认值保持 `true`，以便显式区分“逐历元噪声”和“条件于产品”的诊断模型；增加 `PPP_EXTERNAL_BIAS_APPLICATION` 遥测与审计脚本。
- `b13df92`: 为 phase bias 增加 `use_formal_sigma_as_state_prior`。启用时，Bias-SINEX 均值和方差只初始化卫星/信号相位偏差状态一次；相位产品方差不再逐历元进入 `R`。
- `24e7ed3`: 增加原逐历元模型但关闭 receiver chunking 的控制组，以隔离分块影响。

默认行为没有被静默改变。只有实验配置显式启用新语义。核心修改位于 `src/cpp/common/acsConfig.hpp`、`src/cpp/common/acsConfig.cpp` 和 `src/cpp/pea/ppp_obs.cpp`。

## 对照设计与结果

所有实验使用完全相同的 2024-07-17 数据、站点、信号、轨道、钟差、OSB、门限和整数变换。M1-no-chunk 与 M2 都关闭 receiver chunking，因此二者差异直接隔离相位 OSB 的时间不确定度模型。

| 模型 | 相位 OSB 方差语义 | receiver chunking | 反馈历元率 | 六站历元反馈率 | 候选秩覆盖 | NIS/row 中位数 | shadow 最大状态误差 | 14:00–24:00 AR 3D RMS | FLOAT 3D RMS | AR-FLOAT 3D RMS |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 旧 M1/Exp2g | 每 30 s 独立加入 | 开 | 1483/2880 = 51.49% | 3259/17280 = 18.86% | 10.49% | 1.206 | 1.04e-9 m | 0.07210 m | 0.07219 m | 0.00406 m |
| M1-no-chunk/Exp3e | 每 30 s 独立加入 | 关 | 2190/2880 = 76.04% | 6090/17280 = 35.24% | 19.72% | 2.383 | 2.377 m | 0.03947 m | 0.03985 m | 0.00387 m |
| M0/Exp3a | 条件于 code+phase OSB，均不加产品方差 | 原配置 | 2808/2880 = 97.50% | 16737/17280 = 96.86% | 83.84% | 8.049 | 1.916 m | 0.02489 m | 0.02522 m | 0.00204 m |
| M2/Exp3c | phase OSB 为持久状态先验；code OSB 仍逐历元 | 关 | 0/2880 = 0% | 0/17280 = 0% | 0% | 无反馈 | 无反馈 | 0.02523 m | 0.02523 m | 0 m |

这里的“反馈历元率”只表示至少一个站在该历元出现 `FILTER_CALL_RETURNED_SUBMITTED_UNVERIFIED`；“六站历元反馈率”以 2880×6 为分母。它们不是经过独立整数真值验证的 ambiguity success rate。M0 和 M1-no-chunk 的高数值不能当作可信固定率：M0 的 NIS 和 postfit 迭代警告明显恶化；M1-no-chunk 的 shadow 状态误差达到 2.377 m，说明反馈后的实际解与一次线性预测不一致。

M2 的精度补充：02:00–24:00 三维 RMS 为 0.0271107 m，14:00–24:00 为 0.0252266 m。两窗中 AR-FLOAT RMS 均为 0，因为没有整数伪观测进入解。

## 外部偏差应用审计

M2 AR 与 FLOAT 各检测到 496919 条外部偏差应用记录，唯一键数均为 496919，缺失为 0，重复为 0。相位 L1C/L2W 共 248429 条记录使用 `PERSISTENT_STATE`，可见的 62 个卫星/信号相位偏差状态均恰好初始化一次。code L1C/L2W 共 248490 条记录仍使用 `PER_EPOCH`，这是当前 M2 的明确边界。

M2 未出现 receiver chunking 重叠警告，也未出现 `/AR` 最大 postfit iteration 警告。M0 出现 216 次 `/AR` iteration-limit 警告；M1-no-chunk 与 M2 均为 0。

## 使用的数据

- 日期：2024-07-17，DOY 199，GPS week 2323，30 s，全天 2880 历元。
- 系统与信号：GPS only，L1C/L2W；G01 因 WUM SP3/CLK/OSB 完整性不足而排除。
- 站点：MARS、GRAZ、MATE、DYNG、SOFI、NICO。
- 观测：`MARS00FRA_R_20241990000_01D_30S_MO.rnx`、`GRAZ00AUT_R_20241990000_01D_30S_MO.rnx`、`MATE00ITA_R_20241990000_01D_30S_MO.rnx`、`DYNG00GRC_R_20241990000_01D_30S_MO.rnx`、`SOFI00BGR_R_20241990000_01D_30S_MO.rnx`、`NICO00CYP_R_20241990000_01D_30S_MO.rnx`。
- 精密轨道：WUM `WUM0MGXFIN_20241980000_01D_05M_ORB.SP3`、`WUM0MGXFIN_20241990000_01D_05M_ORB.SP3`、`WUM0MGXFIN_20242000000_01D_05M_ORB.SP3`。
- 精密钟差：WUM `WUM0MGXFIN_20241990000_01D_30S_CLK.CLK`。
- 外部偏差：WUM `WUM0MGXFIN_20241990000_01D_01D_OSB.BIA`。
- ERP：WUM `WUM0MGXFIN_20241990000_01D_01D_ERP.ERP`。
- 广播导航：`BRDC00IGS_R_20241990000_01D_MN.rnx`。
- 天线、参考坐标与模型表：`igs20_2303.atx`、`IGS0OPSSNX_20241960000_07D_07D_CRD.SNX`、IGRF14、DE436、GPT2.5、OLOAD/ALOAD、ocean-pole loading、卫星 metadata 和 yaw bias-rate 表。

完整文件大小和 SHA-256 见 `EXPERIMENT_1_2024_DATA_MANIFEST.md`，本次运行回执 `experiment_3_osb_temporal_model_run_receipt.json` 也冻结了关键输入的 SHA-256。

Bias-SINEX 头部声明 `PARAMETER_SPACING 86400`、`BIAS_MODE ABSOLUTE`、`TIME_SYSTEM G`、GPS clock reference `C1W C2W`、APC `IGS20_2303.ATX`。CLK 为 30 s WUM/PANDA 产品，并声明钟差与相位及 P1/P2 码一致、已应用 `CC2NONCC V6.5 P1C1.DCB`。因此不能在 CLK+absolute OSB 之外再叠加 CAS DCB；那会形成重复改正或混合基准。若要测试 CAS DCB，只能做替换式控制实验。

## 构建、测试与输出

- 源码基线：官方 `origin/main` 的 `7baa32a`，独立分支 `codex/full-rank-pppar-datum`。
- 当前结果提交前 HEAD：`24e7ed3`；核心 M2 修复提交：`b13df92`。
- PEA SHA-256：`3AB177D7EE7C06936D02384368FE0F090CCEE0A16D4E6B6049FADDE8B486B23C`。
- C++ `gnss_ambiguity_diagnostics_tests`：通过。
- Python 外部偏差、反馈信息和坐标审计测试：13 项通过。
- M2 AR 输出：`/home/rx/GINAN/inputData/outputs/exp3c_persistent_phase_osb_pppar_2880e_20260909_b13df92_bin3ab177d7`。
- M2 FLOAT 输出：`/home/rx/GINAN/inputData/outputs/exp3d_persistent_phase_osb_float_2880e_20260909_b13df92_bin3ab177d7`。
- M1-no-chunk AR 输出：`/home/rx/GINAN/inputData/outputs/exp3e_per_epoch_osb_no_chunk_pppar_2880e_20260909_24e7ed3_bin3ab177d7`。
- M1-no-chunk FLOAT 输出：`/home/rx/GINAN/inputData/outputs/exp3f_per_epoch_osb_no_chunk_float_2880e_20260909_24e7ed3_bin3ab177d7`。

## 尚未闭合的边界

当前 M2 只验证了单个 86400 s Bias-SINEX 区间。多日连续处理时，产品区间边界上的状态重置/跳变尚未实现。M2 也没有产品提供的跨卫星、跨信号公共模协方差；code OSB 方差仍按逐历元处理。因此 M2 是纠正最严重时间白化问题的最小实现，不是完整的外部产品协方差模型。

下一步不应降低成功率或 ratio 门限来“制造固定”。应优先获得带有明确 CLK/OSB 协方差与基准说明的同源产品，或在有证据支持时实现公共低秩模式加残差状态（M3），再做多日、重启和独立保留站验证。
