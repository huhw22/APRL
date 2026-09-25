# 实验室系探测器 HDF5 输出

探测器文件独立于每 rank 的轨迹文件。它们使用普通串行 HDF5，并且只由 MPI rank 0 打开；没有 rank 共享 parallel-HDF5 文件句柄。

全部数值采用 SI，事件在写入前都变换到实验室系。按配置正常停止时 `complete=1`；交互模式协调中断或耗尽最大 `mesh.duration` 时为 `0`。读取器只能使用相应 `committed_*` 标量给出的已提交前缀。

## 场探测面

每个场面生成 `<name>.h5`，包含 `/field_plane`：

- `time_s`：float64，形状 `[samples]`，实际实验室采样时刻；
- `fields`：复合记录，形状 `[samples,ny,nx]`，成员为 `electric_V_per_m[3]` 和 `magnetic_T[3]`；
- `committed_samples`：uint64，可读样本前缀；
- `complete`：uint8 完成标志。

组属性包括 `plane_z_m`、`x_first_m`、`y_first_m`、`dx_m`、`dy_m`、`nx`、`ny`、`reference_entrance_z_m`、`reference_distance_m`、`reference_rho_guard_m` 和 `reference_gamma_guard`。x、y 为单元中心坐标。E/B 从 Yee 位置在空间上配准，不分配三维共点副本。B 在已存 leapfrog 半时刻采样，`magnetic_time_stagger` 属性明确记录这一点。

场面保存 Maxwell 网格场，包括自洽场和未来的 Maxwell 注入种子场；规定磁元件场不会复制到这一面向辐射的输出，后处理若需要可从输入卡重新获得。

文件始终保存原始总场。`particle_background_status` 指示是否生成弹道参考 companion，`particle_background_validation` 指示是否启用双面比较。`external_background_subtracted` 当前为 `none`；探测器 API 已接受可空的实验室系本底采样器，未来可扣除解析注入激光而不修改传播网格。

`particle_current_policy` 记录普通物理粒子或实验性 C2 五次退休电流。退休文件还含 `particle_retirement_entrance_z_m` 与 `particle_retirement_exit_z_m`。这些都是审计元数据；写出时不会暗中扣除基线。

若开启与退休层绑定的频率保护，`retirement_frequency_protection` 表示初始化检查已通过，并保存：

- `retirement_protected_minimum_photon_energy_eV`；
- `retirement_protected_minimum_cycles`；
- `retirement_required_length_m`；
- `retirement_observed_cycles_at_minimum_energy`；
- `retirement_causal_guard_distance_m`；
- `retirement_causal_guard_entrance_z_m`。

这些属性使用户无需重新打开粒子输入，就能审计手动长度和实际最大粒子 gamma 检查。

文件不重复保存 Poynting 通量；后处理可由 `S=E cross B/mu0` 计算。`field_reconstruction` 同时计算原始场与扣除粒子本底后的前向功率和时间积分能量，见[粒子本底场重构](FIELD_RECONSTRUCTION.md)。`field_plane_analysis` 只使用电场构造前向真空角谱、积分能谱、Stokes 量及窗口系综互谱密度，可以直接读取本组或 `/reconstructed_field`，见[场面频谱与相干性分析](FIELD_PLANE_ANALYSIS.md)。

## 弹道粒子本底参考

对名为 `<name>` 的场面，显式启用 `particle_background` 会生成 `<name>-ballistic-reference.h5`，其中 `/ballistic_reference` 含：

- `records`：粒子向下游穿过自动推导左边界时的一条复合记录；
- `committed_records`：uint64 可读记录前缀；
- `complete`：uint8 完成标志。

记录成员为 `particle_id`、`source_id`、`time_s`、`position_m[3]`、`proper_velocity[3]`、`charge_C`、`mass_kg`、`weight` 和 `crossing_direction`。载荷格式与粒子面相同，但用途不同：它初始化一条虚拟匀速世界线，不冻结或替代真实粒子。

属性含参考入口 `plane_z_m`、`field_detector_z_m`、保护几何、`role`，以及明确的 `affects_particle_push=false`、`affects_maxwell_current=false` 契约。只有 rank 0 打开文件。参考事件复用探测器有界 MPI 批次，每个粒子只写一次，而不是按轨迹节拍写入。

### 可选双面验证

当 `particle_background.validation.enabled: true` 时，同组增加：

- `validation_records`：每个按 `particle_id` 匹配的参考入口/场面记录；
- `committed_validation_records`：uint64 可读验证前缀。

每条验证记录包含真实场面时刻、位置和固有速度，直线预测时刻和位置，以及横向位置误差、带符号到达时间误差、固有速度相对变化和方向误差。入口状态仍保存在 `records`，可按 ID 连接而不重复。

正常或协调中断关闭时，属性报告匹配数、未匹配入口/出口、重复入口、无效预测以及 RMS/最大误差。部分运行可能仅因粒子尚未到达场面而有未匹配入口，读取器仍必须服从两个 committed 前缀。该路径只用于小束团，并由 `maximum_particles` 限制。

## 粒子探测面

每个粒子面生成 `<name>.h5` 和 `/particle_plane`：

- `records`：可扩展复合数据集；
- `committed_records`：uint64 可读前缀；
- `complete`：uint8 完成标志。

每条记录包含 `particle_id`、`source_id`，`time_s`，`position_m[3]`（z 精确等于面位置），`proper_velocity[3]`，`charge_C`、`mass_kg`、`weight`，以及 `crossing_direction`（`+1` 向下游，`-1` 向上游）。过面事件在完整粒子步的实验室时空中线性定位；若粒子反向并再次过面，可以合法出现多次。

## I/O 与轨迹隔离

固定实验室面在相对论变换系 MPI z 薄层中移动时，场归属会改变；当前拥有者只向 rank 0 发送一个瞬时 x-y 样本。粒子事件只在确实过面的步聚合，并以有界点对点批次传输。rank 0 串行写探测器，其他 rank 从不打开文件。

探测器写入器不会调用、刷新或修改轨迹写入器；两者的数据集、缓冲、完成标志和文件名完全独立。正常结束时先关闭轨迹，再关闭探测器，不存在两条输出路径抢写。

### 高 gamma 能损分析

粒子面和轨迹输出的逆变换使用光前动量分量，避免粒子 gamma 与 boost gamma 接近时两个大数相减。但远小于束流总能量的能损仍不能通过两个分别舍入的总能量相减获得。

对于两个粒子面，应按 `particle_id` 连接选定的下游过面记录，由三个固有速度分量计算各自 gamma，先逐粒子形成

```text
delta_E_i = mass_kg_i * c^2 * (gamma_in_i - gamma_out_i),
```

再用补偿或扩展精度求和。`charge_C` 与 `mass_kg` 已包含宏粒子缩放，不能再次乘诊断用相对 `weight`。方向、重复过面、粒子丢失和不完整前缀必须显式处理。

粒子能损与一个场面收集的辐射不是普适硬等式：有限孔径、其他边界通量、缺失配对粒子、规定外场和不完整时间窗都会改变账本。因此能量闭合是问题相关的验证报告，而不是模拟启动条件。
