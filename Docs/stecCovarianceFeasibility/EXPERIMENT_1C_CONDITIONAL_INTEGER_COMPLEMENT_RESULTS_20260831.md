# Experiment 1C：条件整数补空间诊断结果（2024-07-17）

## 结论

本轮修复了“部分 LAMBDA 只保留已接受尾部行，无法继续检查剩余整数子空间”的诊断缺口，并确认历史 `至少 3 维` 门槛会遮蔽可通过 ratio 检验的 1 维余量。它没有完成全秩 PPP-AR。

480 历元前向解中，第一阶段有 199 个 `SUBMITTED_UNVERIFIED` 历元、4193 条精确独立提交行；全部为 GPS L1C/L2W 宽巷型，单频/窄巷行为 0。条件补空间探针在 116 个历元接受 383 条附加独立行，但这 383 条也全部是宽巷型。组合整数秩最大 29，而逐历元目标秩远高于此；组合秩缺口中位数仍为 56，0 个历元达到全秩。

因此当前代码只能说明：

- 条件均值、Schur 补协方差和未使用的幺模行可以被一致地构造；
- 允许 1--2 维诊断搜索后，确实能发现被原硬门槛遮蔽的宽巷候选；
- 通用 LAMBDA 在当前浮点协方差中仍只识别宽巷子空间，未建立完整 per-signal 整数基准；
- 第二阶段必须继续保持 `PROBE_ONLY_NOT_SUBMITTED`，不能反馈到滤波器。

## 代码与来源

- 上游：`https://github.com/GeoscienceAustralia/ginan.git`
- `origin/main` 基线：`7baa32a70f819c475b03d0c6834d9a6c90a67c1b`
- 实验分支：`codex/full-rank-pppar-datum`
- 本次源码提交：`f9826eddbc8e3c99bd71b2d2bf121d32235208f2`
- PEA SHA-256：`3c815eb2a05ce0a71c646ef5022b8357c9079828a9cf11aac03aa5130c3effe9`

这不是上游 GINAN 原封不动的二进制。该分支从官方 `origin/main` 分叉，并包含 STEC 协方差、接收机单差整数坐标、2024 实验驱动、模糊度反馈修复和本次条件补空间诊断等修改。WSL 内的 CMake 版本脚本不能解析 Windows worktree 的 gitdir，因此程序启动信息显示 `untagged--git-error`；实验收据以源码提交号和二进制 SHA-256 为准。

本轮关键修改：

1. 保留部分 LAMBDA 选择前的完整幺模变换；
2. 使用已接受整数创新计算剩余整数坐标的条件均值和 Schur 补协方差；
3. 把第二阶段候选映射回原始模糊度状态，并用精确有理数消元检查它与第一阶段的独立性；
4. 正式第一阶段继续保持最少 3 维；仅诊断探针允许最少 1 维，以检查历史硬门槛；
5. 每条第二阶段结果强制标记 `PROBE_ONLY_NOT_SUBMITTED`。

## 运行配置

- 日期：2024-07-17（DOY 199）
- 窗口：00:00:00--03:59:30
- 间隔：30 s
- 历元：480
- 测站：MARS、GRAZ、MATE、DYNG、SOFI、NICO
- 整数信号：GPS L1C、L2W
- 配置：`Docs/stecCovarianceFeasibility/experiment_1c_full_rank_diagnostic.yaml`
- 输出：`/home/rx/GINAN/inputData/outputs/exp1c_lowdim_probe_480e_20260831_f9826ed_bin3c815eb2`

完整命令、输出校验和和声明边界见 `experiment_1c_480e_run_receipt.json`。

## 使用的数据

输入清单共有 23 个文件，清单 SHA-256 为 `28a3d3531869625ef5b858b338721866e12a6b6ca5f37a1deee87e1d946dad88`。每个文件的大小与 SHA-256 见 `experiment_1_2024_input_manifest.json`。

### 观测数据

| 测站 | 2024 DOY 199、30 s RINEX 3 文件 |
|---|---|
| MARS | `data/MARS00FRA_R_20241990000_01D_30S_MO.rnx` |
| GRAZ | `data/GRAZ00AUT_R_20241990000_01D_30S_MO.rnx` |
| MATE | `data/MATE00ITA_R_20241990000_01D_30S_MO.rnx` |
| DYNG | `data/DYNG00GRC_R_20241990000_01D_30S_MO.rnx` |
| SOFI | `data/SOFI00BGR_R_20241990000_01D_30S_MO.rnx` |
| NICO | `data/NICO00CYP_R_20241990000_01D_30S_MO.rnx` |

### 精密产品与参考文件

| 用途 | 文件 |
|---|---|
| WUM 30 s 精密钟 | `products/WUM0MGXFIN_20241990000_01D_30S_CLK.CLK` |
| WUM OSB Bias-SINEX | `products/WUM0MGXFIN_20241990000_01D_01D_OSB.BIA` |
| WUM 精密轨道，前一日 | `products/WUM0MGXFIN_20241980000_01D_05M_ORB.SP3` |
| WUM 精密轨道，实验日 | `products/WUM0MGXFIN_20241990000_01D_05M_ORB.SP3` |
| WUM 精密轨道，后一日 | `products/WUM0MGXFIN_20242000000_01D_05M_ORB.SP3` |
| 地球自转参数 | `products/WUM0MGXFIN_20241990000_01D_01D_ERP.ERP` |
| 广播导航 | `products/BRDC00IGS_R_20241990000_01D_MN.rnx` |
| IGS 周解坐标 SINEX | `products/IGS0OPSSNX_20241960000_07D_07D_CRD.SNX` |
| 天线相位中心 | `products/igs20_2303.atx` |
| 卫星元数据 | `products/tables/igs_satellite_metadata_2203_plus.snx` |
| GPS 姿态参数 | `products/tables/sat_yaw_bias_rate.snx` |
| 行星历表 | `products/tables/DE436.1950.2050` |
| GPT2 网格 | `products/tables/gpt_25.grd` |
| 海潮负荷 | `products/tables/OLOAD_GO.BLQ` |
| 大气负荷 | `products/tables/ALOAD_GO.BLQ` |
| 海洋极潮负荷 | `products/tables/opoleloadcoefcmcor.txt` |
| IGRF-14 地磁系数 | `products/tables/igrf14coeffs.txt` |

未使用 Zhang 内部 `HOU_OSB_LIKE`、GIM 或 ATT/OBX 产品。G01 因实验日所选 WUM 精密轨道、钟和两种目标 OSB 信号记录不完整而被排除。

## 第一阶段：实际提交结果

| 指标 | 结果 |
|---|---:|
| AR 尝试历元 | 480 |
| `SUBMITTED_UNVERIFIED` 历元 | 199 |
| 提交行 | 4193 |
| 精确秩和 | 4193 |
| 重复/线性相关行 | 0 |
| 宽巷型行 | 4193 |
| 单频/窄巷行 | 0 |
| 全秩历元 | 0 |
| 秩缺口 min / median / max | 38 / 60 / 77 |

与既有正式 primary trace 比较：199 个共同提交历元、997 个历元-接收机组、11223 个隐含宽巷星对差全部一致，整数差不一致数为 0。这证明新增诊断没有改变第一阶段结果；它不证明这些宽巷整数独立正确。

## 第二阶段：条件补空间探针

| 指标 | 结果 |
|---|---:|
| 有条件补空间的历元 | 199 |
| ratio 接受历元 | 116 |
| ratio 拒绝历元 | 9 |
| 成功率拒绝历元 | 74 |
| 探针接受行 | 383 |
| 宽巷型行 | 383 |
| 单频/窄巷行 | 0 |
| 其他结构行 | 0 |
| 组合秩 min / median / max | 6 / 25 / 29 |
| 组合秩缺口 min / median / max | 38 / 56 / 76 |
| 接受 ratio min / median / max | 10.726 / 39.451 / 248224.217 |

全部 383 条行在原始模糊度坐标中均与对应第一阶段行精确线性独立，但它们没有扩展到窄巷子空间。短窗中原先被 `至少 3 维` 拒绝的 1 维候选在允许诊断搜索后通过 ratio，说明硬门槛会减少可见宽巷行；即使去掉这一遮蔽，完整整数秩仍未形成。

## 坐标与协方差检查

- 前向原始 STEC 协方差：480/480 有效；
- 前向星间差 STEC 协方差：480/480 有效；
- 平滑原始 STEC 协方差：480/480 有效；
- 平滑星间差 STEC 协方差：480/480 有效。

与独立 FLOAT 控制的 2880 个共同测站-历元坐标全部匹配覆盖。AR-FLOAT 三维差最大 0.011676 m、RMS 0.002745 m，与既有 primary 结果相同，因为第二阶段没有反馈。没有预注册坐标阈值，这只是完整性筛查，不是科学通过。

运行仍出现 BLOCK IIIA 姿态模型回退、偏差外推、若干测站/卫星二频码观测缺失并回退单频，以及 G01 缺少精密星历等警告。这些警告限制模型解释，但不能单独解释为什么所有可接受整数方向都是宽巷。

## 为什么加入 WUM 产品后仍没有形成完整 PPP-AR

WUM 精密钟和 OSB 已进入本次处理，但产品只能消除或约束其定义范围内的卫星端偏差；它不会自动增加浮点解中窄巷方向的信息量，也不会替代接收机端整数基准构造。当前接收机-信号星间单差变换理论上给出了完整的 per-signal 整数坐标，但 LAMBDA 的高置信方向持续落在 L1C-L2W 宽巷子空间。条件化掉一部分宽巷后，下一批高置信方向仍是其他宽巷图边；当宽巷子空间接近饱和后，剩余 per-signal/窄巷方向没有同时满足成功率与ratio条件。

因此此前“加入产品效果不好”的直接原因不是产品完全没有被读入，而是把“产品已应用”误当成了“完整整数可观测性和正确整数基准已建立”。现有证据支持产品对宽巷一致性的作用，不支持窄巷或全秩固定。

## 下一步

不应直接把条件补空间探针改成反馈。下一阶段应显式实现和验证顺序整数结构：

1. 先构造并验证每接收机的完整宽巷图基，而不是让部分 LAMBDA 任意选宽巷行；
2. 在固定宽巷条件下，建立严格整数、幺模且单位一致的窄巷/per-signal 余量变换；
3. 对每条窄巷候选记录产品覆盖、弧段、周跳、条件方差、成功率、ratio 和原始状态映射；
4. 只有在独立重启、错误固定对照、FLOAT 对照、坐标完整性和独立 2024 日期均通过后，才考虑反馈与 PPP-AR 接受声明。

当前状态：`DIAGNOSTIC_ONLY_FULL_INTEGER_RANK_NOT_ACHIEVED`。
