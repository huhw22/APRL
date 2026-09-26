# 场探测面的弹道参考区

## 目的与范围

实验室系场面保存传播后的总 Maxwell 场。在相对论电子束附近，总场同时包含目标辐射和束团的束缚带电粒子场。场探测器因此支持一个紧凑参考，供后处理重构粒子匀速延续的场，并评估或扣除其贡献。

该机制只用于诊断。真实粒子仍由自洽推进器处理，照常沉积电流、跨 MPI rank 迁移并响应场；没有粒子被冻结或删除，Maxwell 内核也不改变。

## 仅左侧几何

对实验室位置 `z_D` 的场面，启动时推导

```text
rho_guard = sqrt(size_x^2 + size_y^2)
gamma_guard = maximum initial laboratory particle gamma
L_reference = gamma_guard * rho_guard
z_F = z_D - L_reference
```

使用完整横向网格对角线是有意的保守选择。粒子向下游穿过 `z_F` 时提供一个实验室系状态，定义虚拟匀速世界线。区间 `[z_F,z_D]` 称为场面的参考区，但它不是新的模拟域，也不分配体场数组。

只需左侧区间。参考用于预测场面采样事件处的带电粒子贡献；该事件之后的粒子运动不能改变已经采样的事件，因此对称右侧缓冲只会增加成本，不增加该探测器信息。这仍是需要针对束团、探测孔径和漂移长度验证的近似，不是严格的辐射/束缚场数学投影算符。

`z_F` 与 `z_D` 都是固定实验室平面，在相对论变换计算系中则是运动的时空曲面。探测器把每条已完成粒子线段的端点变换到实验室系，并在其中插值过面事件，避免错误地把两个面当成计算系内同时静止的固定 z 平面。

## 放置规则

- 磁作用区之间不得重叠；
- 磁作用区不得与 `[z_F,z_D]` 重叠；
- 不同场面的参考区可以互相重叠；
- 粒子探测面不参加排斥检查；
- `after-last-element` 对场面使用 `z_D`，不使用 `z_F`；
- `z_D` 下游可以存在元件，只禁止与左参考区重叠。这与绑定真实粒子退休、必须成为最后束线元件的 `retirement_frequency_protection` 场面不同。

重叠会在初始化阶段被拒绝。错误信息根据最后一个磁作用区出口、参考长度和一个相对论变换系纵向单元余量，给出建议的最小场面位置。

## YAML 配置

```yaml
detectors:
  enabled: true
  directory: output/detectors
  field_planes:
    - name: radiation
      z: 0.4
      rhythm: 0.0002
      buffer_samples: 2
      compression: 0
      particle_background:
        enabled: true
        buffer_records: 16384
        compression: 0
        validation:
          enabled: false
          maximum_particles: 100000
```

参考文件名为 `<name>-ballistic-reference.h5`。`particle_background` 默认关闭，必须显式开启。省略或关闭时，不创建参考文件、过面事件、缓冲、通信或参考区；总场文件 `<name>.h5` 不受影响。

## 小规模双面验证

可选验证在 `z_F` 和 `z_D` 各采样同一真实粒子一次，由 rank 0 按 ID 配对，并把第二个状态与入口直线传播比较。入口固有速度为 `u` 时：

```text
gamma = sqrt(1 + |u|^2)
v = c u / gamma
Delta t = (z_D - z_F) / v_z
r_pred = r_F + v Delta t
```

文件为每个匹配粒子保存一条验证记录，并用属性汇总横向位置误差、带符号和绝对到达时间误差、固有速度相对变化及速度方向变化。

这不是采样轨迹；通信与存储规模是每粒子两个过面事件，而不是粒子数乘时间步。但 rank-0 配对表仍随宏粒子数增长，所以验证默认关闭，启动时强制 `maximum_particles`（指宏粒子数量，不是代表电子权重）。关闭后不建立配对表、验证数据集，也不产生第二次场面粒子事件。

比较只衡量所选参考距离内无外力延续的准确度，不能单独证明后续场扣除完全精确。验证顺序为：先检查双面轨道误差，再仅对很小束团开启完整轨迹，并比较轨迹远场与预期解析谱。

## 本底扣除边界

主场文件保存原始总 Maxwell 场。可空的实验室系本底采样接口为解析种子激光扣除保留；第一版未接入该采样器，并记录 `external_background_subtracted=none`。规定波荡器场不会写入辐射场面，因此无需扣除。

把本参考与原始 E/B 面组合成清洁场的独立 C++ 工具见[粒子本底场重构](FIELD_RECONSTRUCTION.md)。它也可以读取与场面同位置的粒子探测面；通常优先使用弹道参考，因为无需启用第二个探测器。
