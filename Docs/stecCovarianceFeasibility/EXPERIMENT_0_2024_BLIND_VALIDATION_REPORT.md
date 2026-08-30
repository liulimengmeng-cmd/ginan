# Experiment 0（2024）FLOAT pivot 正式盲测结果（非 fixed-STEC）

## 1. 结论

本次正式盲测不支持“在空间模型训练中使用完整卫星差分 STEC 协方差（full），优于忽略模型网络训练协方差（none）或仅保留其对角项（diagonal）”这一机制主张。

- quiet 日，full 校准后总体 95% 区间覆盖率为 94.048%，但相对 full 未校准区间的平均宽度膨胀 94.612%，超过预注册上限 40%。full 相对 none 的覆盖误差仅改善 0.475 个百分点，区间得分反而恶化 4.873%，H1 类比判据失败。
- storm 日，quiet 冻结的 full 校准转移后覆盖率降至 67.681%，最差留出站覆盖率为 63.988%。full 相对 none 的覆盖误差恶化 26.412 个百分点，区间得分恶化 149.033%，H1 类比判据失败。
- quiet 日 full 相对 diagonal 的平均对数得分改善为 -0.1157 nat/预测，95% 整块 bootstrap 区间为 [-0.1800, -0.0578]，方向明确偏向 diagonal，而不是 full。storm 日点估计也偏向 diagonal，但区间跨 0，不能作方向性结论。
- quiet 到 storm 的 full 校准覆盖损失点估计为 26.367 个百分点，超过预注册的 5 个百分点阈值（按点估计判定），因此必须进行 storm 条件化；该转移量没有单独的 bootstrap CI，不能作显著性表述，也不能声称 quiet 标定可普遍转移。

这不是原始 PPP-AR 固定 STEC 假设 H1–H5 的检验。两日 AR 候选历元均未达到 G2（70%）门槛，正式输出阶段为 `FILTER_POSTERIOR_NO_EPOCH_AR`。因此本报告只评价 FLOAT STEC 协方差传播和区间校准流水线，H5 错误固定风险没有被测试。

## 2. 分支与软件来源

- 工作树：`D:\tec\ginan-main-stec-cov`
- 工作分支：`codex/main-stec-covariance-feasibility`
- 分支基点：官方 `origin/main`，`7baa32a70f819c475b03d0c6834d9a6c90a67c1b`
- 盲测冻结时 HEAD：`77029b9edd82767a2e5e09f207c76803f07a0418`
- PEA 二进制内嵌提交：`1296966`
- PEA SHA-256：`6131cac3b37e160bd351b4b23b5830b0dca5de42a88c71cdd6d8c2034480d9c4`
- 冻结与评估均记录 `git.exe` 后端、正确分支、正确 HEAD 和 tracked-clean 状态。

操作记录（不属于正式 receipt 科学证据链）：旧的 Zhang 工作树及其未提交改动没有被修改。本分支从官方 main 建立，不是从 Zhang 分支派生。

## 3. 使用的 2024 数据

完整逐文件清单、相对路径、字节数和 SHA-256 见 `EXPERIMENT_0_2024_DATA_MANIFEST.md` 及机器清单：

`/home/rx/GINAN/inputData/experiment0_2024/input_manifest.json`

机器清单 SHA-256 为 `7ba3f1e1661cc2bf6592b4cdb8e10e310ea4bd7eb25120cd33e4a8ecc318f45a`，共 53 个唯一文件、1,590,946,109 字节：26 个 RINEX 观测文件、16 个逐日产品、9 个共享辅助模型和 2 个空间天气选日文件。

### 3.1 日期、系统、信号和采样

- quiet：2024-05-08（DOY 129），max Kp = 2，min Dst = +3 nT。
- storm：2024-05-11（DOY 132），max Kp = 9，min Dst = -406 nT。
- 每日 GPST 00:00:00–23:59:30，30 s 采样，理论 2880 历元。
- 仅处理 GPS；信号为 `C1C/L1C/C2W/L2W`。

### 3.2 观测站

- 空间模型/校准站（10）：ARMC、BATH、HOB2、MCHL、MOBS、STR1、SYDN、TID1、WGGA、YARR。
- 完全留出站（3）：STR2、BALL、PARK。
- 每站分别使用 quiet 和 storm 两天 RINEX，共 13 × 2 = 26 个文件，来源为 Geoscience Australia GNSS Data Centre。
- quiet 模型站数据用于空间模型和 κ 冻结；storm 不重新估计 κ。三个留出站不进入空间模型拟合、κ 选择或方法选择。其 RINEX 必须被 manifest、QC 和留出 PEA 读取，但 freeze 前不读取由留出 PEA 生成的 STEC/SD 输出；这些输出只进入唯一一次正式统计评估。

### 3.3 精密产品与逐日产品

每个实验日使用：

- WUM `WUM0MGXFIN` 最终轨道：目标日前一日、目标日和后一日共 3 个 5 min SP3；两天共 6 个文件。
- WUM `WUM0MGXFIN` 30 s 卫星钟：每日 1 个 CLK；两天共 2 个文件。
- WUM `WUM0MGXFIN` ERP：每日 1 个；两天共 2 个文件。
- WUM `WUM0MGXFIN` Bias-SINEX/OSB：每日 1 个 BIA；两天共 2 个文件。
- IGS 广播导航：每日 1 个 BRDC RINEX，由 Geoscience Australia API 获取；两天共 2 个文件。
- Geoscience Australia 站坐标 APR/SINEX：每日 1 个 AUT SINEX；两天共 2 个文件。

处理采用 PRECISE-only 策略；未在缺口处切换广播轨道/钟、其他分析中心产品或插值产品。WUM BIA 的卫星码与相位 OSB 均作为输入使用；它不是本项目内部 Zhang `HOU_OSB_LIKE` 产品。

### 3.4 共享辅助模型

共 9 个共享文件：

- `igs20_2303.atx`（IGS20_2303 天线/APC）；
- `igs_satellite_metadata.snx`；
- `sat_yaw_bias_rate.snx`；
- `igrf14coeffs.txt`；
- `DE436.1950.2050`；
- `gpt_25.grd`；
- `OLOAD_GO.BLQ`；
- `ALOAD_GO.BLQ`；
- `opoleloadcoefcmcor.txt`。

它们分别支持天线改正、卫星元数据/姿态、二三阶电离层地磁模型、行星星历、GPT2 对流层以及海潮/大气潮/海极潮加载。

### 3.5 空间天气选日文件

- GFZ Kp JSON：`gfz_kp_2024-05-08_2024-05-12.json`。
- Kyoto provisional Dst：`kyoto_dst_provisional_2024-05.txt`。

Kp/Dst 只用于预先选择 quiet/storm 日期，不进入 PEA 状态估计、空间模型或评分特征。

### 3.6 已知缺口与明确未使用的数据

- quiet：G32 CLK 在 23:56:00–23:59:30 缺 8 个历元。
- storm：G02–G32 CLK 在 03:35:00–03:39:30 共同缺 10 个历元。
- 因上述全日完整性缺口，正式 product audit 顶层为 `validation_passed=false`，其中 `metadata_and_eligibility_passed=true`、`full_day_completeness_passed=false`。实验按预注册规则排除缺口并保留 2880 历元分母，不把产品描述为“全日完整”。
- storm covariance 审计预登记 7 个 `NO_STATES` 历元：GPS week 2313、TOW 531420–531600，30 s 间隔。
- 最终空间评分另按预注册优先级 `product gap -> invalid model -> invalid heldout` 排除不可评分历元：quiet 排除 8 个产品缺口历元和 77 个无可用 heldout 预测行历元，保留 2795 个历元；storm 排除 10 个产品缺口、1 个 invalid-model 和 56 个 invalid-heldout 历元，保留 2813 个历元。该分类与 raw covariance 审计的 7 个 `NO_STATES` 不是同一计数口径，重叠历元只记入优先级最高的一类。
- Galileo、GLONASS、BDS、QZSS、SBAS、SLR、外部 GIM 和 OBX 姿态产品未使用。
- 本次 FLOAT 模型站与留出站运行均为 `ambiguity_resolution: false`、`mode: OFF`；没有 AR 固定解或 AR pseudo-observation。这一描述不适用于先前用于判定 G2 的 AR 候选 replication。
- quiet 排除 G01、G10；storm 排除 G01。
- 预下载但未选择的 CODE rapid 产品不在 53 文件正式清单中，也未进入配置或评分。

## 4. 评分对象与协方差臂

本实验不评分绝对 STEC 或 VTEC。每个历元的响应量是留出站某颗卫星相对同站参考星的 FLOAT slant STEC 差分（`SD_STATE`，单位 TECU）。空间平面使用截距、东西向/1000 km 和南北向/1000 km 三个基函数；目标设计行是目标星与参考星 IPP 映射设计行之差。因此，报告中的梯度只是这一指定空间模型下的代理量。

三个协方差臂只改变模型网络空间拟合时对训练 `Q_SD` 的使用：

- none：OLS，忽略模型网络训练 `Q_SD`；
- diagonal：GLS，只使用 `diag(Q_SD)`；
- full：GLS，使用完整模型网络 `Q_SD`。

`none` 不表示最终预测区间没有不确定度。三臂在验证时都加入留出 FLOAT `SD_STATE` 的边际方差。quiet 模型站 LOSO 校准中，full 臂还保留同一网络解内已知的预测—验证交叉协方差，且该交叉项不随 κ 缩放；正式 model/heldout 来自独立 PEA 解，二者之间没有可用交叉协方差，按 0 处理。

## 5. 盲测与可复现性

### 5.1 Analysis-blind amendment 历史

第一版 freeze `spatial_model_freeze_e0c614d.json`（SHA-256 `7f9f7b8eb68c0139bf9bf33e1d36a5d63f1db5a286023f32134fc6a65e6cda01`）创建时尚未读取留出输出。随后第一版 quiet 留出 PEA 完成。正式评分前的源码审计发现：full 臂的模型站 LOSO 校准遗漏已知 test–train 交叉协方差、三臂可能使用不同校准行，以及代码/二进制/执行/替换策略绑定不足。

在任何留出残差、覆盖率、NLL、区间得分或 `SD_STATE` 估值被打印或计算前，自动验证器已经打开第一版 quiet covariance 矩阵，检查序列化、对称性、PSD 和 `D C D^T` 重构，并只暴露 2880/2880 有效计数；另有 header/META/DATUM 结构检查。第一版 quiet 输出保留在 `/home/rx/GINAN/inputData/outputs/exp0_2024_quiet_heldout_float_full_2880e_1296966`，但永久排除于评分，也没有运行第一版 storm 留出解。

随后依据 `EXPERIMENT_0_2024_BLIND_AMENDMENT.md` 修复并提交分析代码，以新名称冻结 V2、创建全新的 quiet/storm 输出根并执行唯一一次评分。因此本实验是 **analysis-blind replacement**：方法修复不依赖留出估值或评分，但不能宣称“任何程序从未打开第一版 quiet 留出文件”。

### 5.2 V2 冻结与正式运行

在读取任何 replacement `blindv2` 留出 STEC/SD 输出前，模型、κ、阈值、输入、代码和执行计划已冻结：

- freeze：`/home/rx/GINAN/inputData/experiment0_2024/spatial_model_freeze_blindv2_20260830a.json`
- freeze SHA-256：`dfb15d010fe0281c023280b85937a49676a32289df06b1f17a2f99e4225a2378`
- freeze receipt：`/home/rx/GINAN/inputData/experiment0_2024/spatial_model_freeze_blindv2_20260830a.receipt.json`
- freeze receipt SHA-256：`a3949e7c81a002f297680db00f84f0f769848d5942c7df45cac50ac219a56831`
- freeze 标记：`heldout_output_read=false`

quiet 模型站 LOSO 校准使用相同的 182,510 条严格共同记录和 24 个一小时块。完整网格为 2880 个历元，排除 8 个注册产品缺口后保留 2872 个校准历元；其中 1100 个历元至少丢弃了一条不满足几何条件的候选行，但不是 1100 个整历元失效：

| 协方差臂 | κ | 模型站校准覆盖率 | 模型站平均宽度（TECU） |
|---|---:|---:|---:|
| none | 5.3351 | 93.638% | 18.901 |
| diagonal | 15.5530 | 93.176% | 17.512 |
| full | 17.2714 | 92.708% | 18.772 |

正式运行回执：

- quiet：`/home/rx/GINAN/inputData/outputs/exp0_2024_quiet_heldout_float_blindv2_20260830a_1296966/run.receipt.json`，SHA-256 `e7f91bf748c19025c820470434666e014a42096f1d6c1893c1117af7b87db810`，2880/2880 covariance 历元有效。
- storm：`/home/rx/GINAN/inputData/outputs/exp0_2024_storm_heldout_float_blindv2_20260830a_1296966/run.receipt.json`，SHA-256 `d1e4962ff4fb07aa06e35cd2d0f5a44e17cb22e972dd4409664d7f00de821d9c`，2873/2880 有效，7 个异常键与预登记完全一致。

唯一一次评估使用连续一小时历元块、有放回重采样、1000 次 bootstrap、seed 20240508。评估前永久 claim 为 `/home/rx/GINAN/inputData/experiment0_2024/spatial_validation_report_blindv2_20260830a.claim.json`，SHA-256 `1e910b971a69bdf962695de98b6e1356f2429c1db1582fd5924b2c74e700f9e3`。评估记录 `inputs_unchanged=true`。

正式输出：

- report：`/home/rx/GINAN/inputData/experiment0_2024/spatial_validation_report_blindv2_20260830a.json`
- report SHA-256：`88539c82840a311df38b22b94e755ee259b5faa09a8d84e206615794f11111ef`
- evaluation receipt：`/home/rx/GINAN/inputData/experiment0_2024/spatial_validation_report_blindv2_20260830a.receipt.json`
- evaluation receipt SHA-256：`15fa03b577c13b13ff3598db1e73397fbf0ec098f58699e0d09e79801cd32a91`

本 Markdown 是正式评估完成后依据上述 machine report 和 receipt 编写的人类可读摘要，不是盲测输入，也没有被既有 evaluation receipt 反向绑定；机器 JSON 与 V2 receipt 是数值和执行完整性的权威记录。

三种协方差臂在同一天使用完全相同的评分行：quiet 为 53,648 行、2795 个可用历元；storm 为 50,283 行、2813 个可用历元。

## 6. 正式结果

以下均为三个留出站合并后的校准结果。RMSE 使用同一日、同一协方差臂的点预测；覆盖率、宽度、NLL 和区间得分评价相应不确定度。

| 日 | 臂 | 覆盖率 | 最差站覆盖率 | 平均宽度（TECU） | 宽度膨胀 | 区间得分（TECU） | NLL（nat/预测） | RMSE（TECU） |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| quiet | none | 96.427% | 95.883% | 17.825 | 80.987% | 21.550 | 2.426 | 3.804 |
| quiet | diagonal | 94.932% | 94.458% | 16.368 | 78.458% | 20.691 | 2.384 | 3.758 |
| quiet | full | 94.048% | 93.174% | 17.860 | 94.612% | 22.600 | 2.500 | 3.826 |
| storm | none | 95.907% | 93.919% | 96.659 | 112.330% | 117.877 | 3.975 | 20.459 |
| storm | diagonal | 67.269% | 66.486% | 45.910 | 24.483% | 281.285 | 19.795 | 23.680 |
| storm | full | 67.681% | 63.988% | 47.201 | 27.985% | 293.552 | 20.707 | 24.650 |

注意：storm none 的覆盖率接近 95%，但其平均宽度达到 96.659 TECU，且相对未校准 none 膨胀 112.330%。这不是高效或高精度结果，不能只依据覆盖率判为成功。

### 6.1 H1 诊断类比：full 对 none

预注册条件为：覆盖误差减少至少 5 个百分点，或区间得分降低至少 10%；同时平均宽度增加不超过 25%。

覆盖误差定义为 `|coverage_percent - 95|`；“覆盖误差减少”是基线臂覆盖误差减去 full 覆盖误差，不是两臂原始覆盖率之差。下表区间均为一小时连续块、1000 次重采样得到的 95% percentile block-bootstrap 区间。

| 日 | 覆盖误差减少（点；95% percentile block-bootstrap 区间） | 区间得分降低（%；95% percentile block-bootstrap 区间） | 宽度增加（%；95% percentile block-bootstrap 区间） | 判定 |
|---|---:|---:|---:|---|
| quiet | 0.475 [-3.272, 2.200] | -4.873 [-10.774, -0.139] | 0.196 [-5.812, 6.027] | 失败 |
| storm | -26.412 [-35.685, -17.126] | -149.033 [-203.040, -85.305] | -51.168 [-57.921, -43.128] | 失败 |

负的“区间得分降低”表示 full 的区间得分更差。quiet 与 storm 的相应区间均完全低于 0。

### 6.2 H2 诊断类比：full 对 diagonal

| 日 | 覆盖误差减少（点；95% percentile block-bootstrap 区间） | diagonal NLL − full NLL（nat/预测；95% percentile block-bootstrap 区间） | 解释 |
|---|---:|---:|---|
| quiet | -0.884 [-2.441, 0.798] | -0.1157 [-0.1800, -0.0578] | 覆盖无实质改善；对数得分明确偏向 diagonal |
| storm | 0.412 [-1.999, 2.738] | -0.9122 [-5.0501, 2.3248] | 覆盖改善未达 3 点；对数得分方向不确定 |

H2 没有预注册的综合 pass boolean。现有证据至少不支持 full 优于 diagonal。

### 6.3 H3 诊断类比：full 校准覆盖

预注册条件为总体覆盖 93%–97%、最差站覆盖至少 90%、相对未校准宽度膨胀不超过 40%。

| 日 | 总体覆盖 | 最差站覆盖 | 宽度膨胀 | 判定 |
|---|---:|---:|---:|---|
| quiet | 94.048%（通过） | 93.174%（通过） | 94.612%（失败） | 失败 |
| storm | 67.681%（失败） | 63.988%（失败） | 27.985%（通过） | 失败 |

### 6.4 H4 诊断类比：quiet 到 storm 的转移

full 校准覆盖损失为 26.367 个百分点，超过 5 点阈值；`requires_storm_conditioning=true`。quiet κ 和空间校准不能作为 storm 条件下的通用校准。

## 7. 科学解释边界

结果直接支持的结论是：在本次澳大利亚 13 站、GPS L1/L2、两个 2024 实验日、FLOAT 卫星差分 STEC 留出站设计下，当前“以完整模型网络 `Q_SD` 进行空间 GLS 拟合”的方法与 quiet 标定没有改善留出预测区间；storm 转移尤其失败。

结果不能证明：

- 非对角 STEC 协方差在物理上没有信息；
- PPP-AR 固定 STEC 下 full 协方差一定无效；
- WUM Bias-SINEX 产品或 GINAN AR 整体无效；
- 留出 FLOAT STEC 是独立电离层真值；
- 两天结果可推广到全年、全球网络或其他星座。

留出 FLOAT 解本身带噪，平面东西/南北系数只是模型依赖的梯度代理。当前失败可能来自协方差传播、空间模型、quiet κ 的非平稳转移、storm 几何/噪声变化或这些因素的组合；本实验不能在它们之间唯一归因。

## 8. 下一步

1. 不以当前结果主张 full 优于 none/diagonal，也不进入论文主结果包装。
2. 先对 storm 进行只使用模型站的预先可见诊断，分解 full/diagonal/none 的预测方差、交叉协方差项、κ 非平稳性和几何/卫星构成变化；不得反向调整本次已完成的盲测。
3. 若继续 PPP-AR 主假设，必须先解决 G2 固定历元不足，并建立错误固定、float 对照、重启/重新初始化和坐标完整性控制，再以新的预注册留出数据检验原始 H1–H5。
4. 使用更多 2024 quiet/storm 日期和独立区域复验；当前两天只能作为可行性筛查。
