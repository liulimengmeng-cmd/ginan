# 本地一小时 AR 有限启动 profile

此次修改仅为已完成 UNKNOWN 修复的 `27bcca` 源码新增精确本地启动能力及验证覆盖；没有修改数值算法。

| 项目 | 精确范围 |
| --- | --- |
| 候选程序 | `/home/rx/GINAN/datum_unknown_ar1h_20261103_01/candidate/pea` |
| 正式输出根 | `/home/rx/GINAN/datum_unknown_ar1h_20261103_01/run` |
| 正式范围 | GPST 2024-07-17 02:00:00–02:59:30；全局 max_epochs=360；ARstart=1405216800 |
| 正式 capture | e300 02:29:30、e360 02:59:30 |
| probe 输出根 | `/home/rx/GINAN/datum_unknown_ar1h_20261103_01/validation_e241/run` |
| probe capture | e241 02:00:00、e242 02:00:30；相同 ARstart |
| 原科学设置 | ratio=2；全部 joint NIS 关闭；每历元独立 AR；WSL20GiB |

旧两个有限 profile、原 e240 源检查点和固定输入身份不变。选中 profile 的 producer 路径绑定权限目标，不增加任意路径、任意时长或软链接绕过。policy native 覆盖新正例及错路径、时长、ARstart、capture 等拒绝情形。

九节 auditor 只为精确新 probe 增加路径授权。PPP_AR 源格式可为 6 或 7，但 before/after 必须同 schema、payload 完全相同；不把 schema6 与 schema7 的内容直接判等。其余九节、33 个 core 字段保持原严格比较及已有明确边界。

本次 8900 producer 完成 158 对象链接（5 fresh、153 verified reuse），51 native 通过（50 本次实际执行、1 本次仅构建、0 adoption）。真实四门恢复验证在元数据记录时进行中，正式小时及科学验收未完成。完整证据身份和限制见 `VALIDATION_AR1H_PROFILE_20261103.json`。历史 27bcca/00a 记录及原 manifest/provenance 字节副本单独保留。
