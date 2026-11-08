# 单向发射因果对照

用户批准的独立验证方法，默认关闭。配置须同时显式指定：

```yaml
mesh:
  particle_coupling: conservative-source
  self_field_feedback: false
  one_way_source_diagnostic: true
energy_ledger:
  enabled: true
  discrete_operator_audit: true
```

保留原 EC 的外场子步公式、真实粒子轨迹及其 F3 电流沉积、参考初场和供场、Cowan、CPML、MPI 迁移和负载均衡。只有粒子所受自洽网格 E/B 被设为零。确定的规定轨迹只需一次 trial；仍通过原事务恢复、真实子步 replay 和一次 CPML 提交流程写入科学输出。原实际 grid work 为零。未启用该模式的原 EC 固定点及所有场、粒子公式保持。

该受迫电流试验不要求粒子支付场的能量，故原粒子加场账本的 `closure_valid` 必须为零。Maxwell 离散账本继续开启；不能将规定源造成的物理能量不闭合称为推进故障，也不能用补充控制器功掩盖生产方法的残差。

## 补充诊断及口径

每 rank 在 reference directory 写 `one-way-source-rank-NNNNN.csv` 和 `.complete`，逐实际完成的 field step 写一行：

- `passive_grid_work_J`：实际由本 rank 沉积的 active segment 上，用真正 commit 后的 F3 匹配 E 时间均值计算的被动 `q E·dx`。
- `controller_work_J`：上项的相反数，代表抵消这一主动电流段的自场受力所需机械控制器功。
- `absolute_passive_grid_work_J`：逐段被动总功的绝对值之和。
- `passive_{x,y,z}_J` 和 `active_segments`：三个分量及拥有的实际沉积段数。

该 CSV 的范围是 **active deposited segments**；frozen、CPML carrier 和纵向 handoff remainder 仍按旧约束账单独记录，未在这里重复加总。因此它本身不是全粒子全场闭合账。控制器值由独立的真实场取样得到，并非原残差倒推。

若已有 `APRL_PARTICLE_FORCE_PROBE`/`_IDS` 选择，原 force v1 二进制格式完全保持，grid E/B 仍表示实际施加的零场。额外输出 `one-way-virtual-force-rank-NNNNN.csv` 和 `.complete`，逐选中实际 replay 的 push 保存前后状态和不反馈的真实 F3 E/B、三个虚拟功分量。此 trace 使用完整 pusher chord 的 rank-open-clip 口径，不能冒充上述 active deposited segment 口径。它适合研究相位相关的被动 Ex 功，而不将实际 grid work 为零误认为场本身消失。

两类补充 CSV 均为 boost 系 SI 量，普通粒子 proper velocity 为 `u=gamma*v/c`。三个 `virtual_B` 是与原 EC 对应的半时间层磁场；场只被读取。Completion marker 为 `1 final_step`。若未启用选中力 probe，不生成 virtual-force 文件。

## 模块验证

`checks.cpp` 链接当前 `direct_eb_core`，接受一个独立结果目录，可用 1/2/4 MPI rank 运行。它检查：

1. 两组不同的非零自场均不能改变规定外场轨迹，逐位对上独立外场 EC 子步重算；有/无静磁两种情形。
2. 实际 grid work 为零，而真实 F3 电流仍会驱动 Maxwell 电源项。
3. 真正 replay/commit 后的被动功非零，独立 native `J·(Eold+Enew)/2` 与 CSV 相符；三分量、控制器及完成步一致。
4. 正常反馈模式仍对强自场产生非零响应，作为负控。

这是方法隔离及诊断一致性验证，不是 FEL 辐射物理验收。7 周期结果必须先与同输入、同长度的原 PIC 配对，不可直接与 37 周期 Genesis 能量相比；两个输出能量之差也不等于纯反馈辐射能，因为复场存在相干交叉项。
