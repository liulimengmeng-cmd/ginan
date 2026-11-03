# 服务端 UNKNOWN 修复候选与用户端源码（2026-11-03）

本目录是当前候选的有效源码入口。分支原有顶层 `src/` 和旧快照保留；当前候选使用此目录的 `server/src/`。

- `server/src/`：全球 134 站 GPS/Galileo/BDS 服务端，950 个源码、构建和测试文件。
- `server/test_support/`：GEC 原生测试的生成声明与来源，3 个文件。
- `user/src/`：US20 最终用户端有效源码，757 个文件，本次逐文件 SHA 核对未变。
- `SOURCE_MANIFEST.json`：1711 条源码及测试支持文件的字节数、SHA256；原有附加字段完整保留。
- `VALIDATION_AR1H_PROFILE_20261103.json`、`AR1H_PROFILE_20261103.md`：当前有限启动 profile 的修改和验证边界。
- `VALIDATION_DATUM_UNKNOWN_20261103.json`、`DATUM_UNKNOWN_REPAIR_20261103.md`：提交 `27bcca`、旧 `00a` producer 的修复和历史验证记录，原字节保留。
- `BASELINE_27BCCA_SOURCE_MANIFEST.json`、`BASELINE_27BCCA_PROVENANCE.json`：该提交的旧清单和来源元数据原字节副本。

## 当前有限启动 profile

当前候选 binary SHA256：`8900c00a35db8d75cb90e90f02a584998764f06a0388ce8f5ab2ab5a78ffe4bc`。

相对 `27bcca5520e1e63574193c3219201f3754730e0f` 仅修改三个源码文件：有限本地启动 profile header、其 policy test 和九节恢复 auditor。数值算法源码未变。原两个 profile 保留；新增精确的 GPST 2024-07-17 02:00:00–02:59:30 direct AR 和 e241/e242 probe，使用独立本地候选及输出根。profile、producer、时间范围、AR 起点和 capture 列表都受原权限门约束。

生产构建通过：5 个单元重编、153 个对象按依赖和身份核验复用、158 个对象重新链接。51 项 native 本次全部重新编译和链接，其中 50 项本次实际执行、1 项仅构建九节 auditor；零旧执行采纳。此结果只证明构建和受控测试通过。

元数据记录时，四门真实整进程恢复验证正在进行。原 e240 九节真实往返、正式一小时 AR 完成情况和科学精度尚未验收。native auditor 构建及 fixture 不能当作真实恢复或正确固定的证明。整网联合 NIS 关闭、ratio=2、每历元独立 AR、WSL20GiB 保持原要求。

## 历史修复与边界

`27bcca` 的旧 binary SHA256 为 `00a16f65eee965038b672b5ad773468b5798d86b1690cecf1771902b6474f657`。该阶段实现最终 QC 后已认证替代路径闭环、增量证据事务、完整 FLOAT 目标运输与同历元收据组合、E11 已发布端点预测、PPP_AR schema7 及历史证据定位索引。其旧构建和 native 采纳链见未改动的 `VALIDATION_DATUM_UNKNOWN_20261103.json`；不能把旧验证记录当成本次 8900 producer 的验证。

完整均值、全部交叉 Q、精确整数和发布门保留。历史索引绑定完整父检查点，不以零交叉协方差恢复在线目标，也不授予当前整数权限。E11 旧检查点无已接受速度历史的冷启动仍可能拒绝；尚不能声称全部 UNKNOWN 消失。真实网络耗时、内存及外部科学精度仍需实际验证。

## 用户端与依赖

用户端保持最终实时版本 binary 身份 `bd977ff2e7e07ef21cb746662fce1f76eafe04073d8cd35be55454ec287146f0` 对应源码，本次未重编或重跑用户实验。两端不是同一次集成构建。

`validation/build_and_test_v4.py` 是旧 27bcca/00a 验收流水线的原字节副本，含冻结依赖和本地证据路径，作为历史复现记录保留。本次构建的独立计划、结果及资格证明由新验证文件中的 SHA256 绑定。原 GINAN 的编译器、BLAS/Eigen/Boost/MongoDB、SOFA/IERS 等外部依赖仍需准备；本次未声称从零冷构建。

没有上传观测、轨钟偏差产品、输出、检查点实体、二进制、对象、大型运行日志、凭据或原始对话。没有访问服务器或替换原生产程序。更早中断的基线状态保留在 Git 历史及 `PROVENANCE.json` 的 `baseline_server` 中。
