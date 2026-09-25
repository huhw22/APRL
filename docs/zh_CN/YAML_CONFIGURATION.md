# YAML 输入卡完整规范

模拟器接受一个普通 YAML 文档。物理解析前，每个 mapping 都按本文档模式检查；未知键是硬错误，并报告 YAML 行号、映射路径和合法键。即使某个可选块已关闭，其中的拼写错误仍会被捕获。五个后处理器也采用相同严格规则。必需量缺失或无效时同样给出所在行；程序避免依赖不透明的物理默认值。

## 单位

```yaml
units:
  length: micrometer
  time: picosecond
```

长度接受 `m`、`mm`、`um`/`micrometer`、`nm`；时间接受 `s`、`ms`、`us`、`ns`、`ps`、`fs`、`as`。这些单位作用于 YAML 中数值型长度和时间，核心只转换一次，内部使用 SI。粒子 HDF5 位置始终是米，不继承 YAML 单位。

## 输出安全与运行身份

```yaml
output:
  overwrite: false
  manifest: output/run-manifest.yaml
```

该顶层块控制本次模拟的全部输出。默认存在任一目标文件即报错；有意整体替换时必须设 `overwrite: true`。rank-0 清单保存原始输入卡、唯一运行 ID、源码/构建信息、MPI 环境和粒子输入身份；相同身份嵌入全部 HDF5，见[运行身份与覆盖保护](RUN_PROVENANCE.md)。

## 网格与 boost

```yaml
mesh:
  field_solver: cowan-z
  cells: [40, 40, 60]
  cell_size: [1.0, 1.0, 0.2]
  center: [0.0, 0.0, 0.0]
  duration: 0.01
  boost_gamma: 2.0
  particle_steps_per_undulator_period: 32
  maximum_particle_substeps: 4096
```

`field_solver` 必填：`cowan-z` 为 z 优先可控色散模板，`yee` 为原始回归路径；两者共享交错 E/B 网格，都不分配 A/phi。Cowan 在 Faraday 计算中直接应用横向平滑，不保存三份平滑体场。

`cowan-z` 要求 `dx>=dz`、`dy>=dz`，在分配场和读粒子前检查，并由 `dz` 给出最小横向 `cell_size`。时间步为 `dt=dz/c`，使已分辨真空模沿 z 精确无色散。启动日志报告纵横比平方、全部 Cowan 系数以及每波长 16 单元时横向轴相/群速度比；后者只是诊断，不表示任意斜模无色散。

当前 TF/SF 种子边界修正仍是 Yee 形式，因此 `cowan-z` 与非空 `incident_waves` 会在读粒子前被拒绝。通用 Cowan TF/SF 完成前，种子注入只能用于 Yee 回归；CPML 仍拒绝非空入射波。

`cells`、`cell_size`、`center` 定义计算系盒子。`cells` 是权威整数数目，每轴至少 3 格；z 还要求每 MPI rank 至少 2 格。总尺寸只由

```text
extent_i = cells_i * cell_size_i
```

得到。旧 `lengths`、`resolution` 被直接拒绝。MPI z 薄层由整数商和余数构造，rank 计数最多相差 1，不用浮点除法确定边界。启动打印网格数、尺寸、推导范围和每 rank z 计数范围。

`duration` 是计算系最大时间保护，不是请求传播距离。变换粒子后打印 `[duration-estimate]`：`reference-center-z` 使用与停止谓词相同的惯性参考世界线精确估计；`after-last-element` 把每个初始粒子弹道外推到最后作用区出口并取最慢者，因此磁力可使其失准。若按步取整的估计超过配置值，只警告并给出 YAML 时间单位建议，不自动加长。

`particle_steps_per_undulator_period` 是最短波荡周期的最少 Boris 采样。读入并 boost 后，程序由任一粒子单个 Maxwell 步的最大实验室 z 前进量选择最小整数子步数。`maximum_particle_substeps` 是正成本上限，默认 4096；超出时硬错误并给出所需上限和场步/z 单元细化建议。

子步会沿轨道重采样解析实验室器件；网格 E/B 在外层场步内保持一次采样，电荷守恒沉积和探测输出使用场步起点到最终位置的弦。因此它提高规定器件轨道积分，不提高 Maxwell Nyquist、不分辨步内辐射电流，也不允许更粗辐射网格。示例中的 1 只用于冒烟，见[粒子子步](PARTICLE_SUBCYCLING.md)。

`boost_gamma` 为整次运行选择一个恒定惯性系，可逐卡修改，但粒子进出元件时不会改变。随时间变化的 boost 是非惯性坐标，会破坏共享 Maxwell 网格及实验室元件/探测面世界线。平面波荡器中静磁偏转保持总实验室 gamma，却降低平均纵向速度；运行 `undulator_resonance INPUT.yaml [LAB_GAMMA]` 可得 `gamma_lab/sqrt(1+K^2/2)` 推荐，比较当前 boost 并估算穿过特征波荡器的盒内 z 漂移。

## 可选辐射分辨率预检查

```yaml
radiation_resolution:
  enabled: true
  maximum_photon_energy_eV: 124.0
  warning_grid_points_per_wavelength: 8.0
  warning_maxwell_samples_per_cycle: 8.0
  warning_detector_samples_per_cycle: 8.0
```

该块声明本次运行希望解释的最高实验室光子能量，不假设 FEL 共振、束团谐波或磁元件类型。程序按傍轴前向 `+z` 模把频率/波长变换到计算系，并报告：

- 每个计算系波长的纵向网格数；
- 每个计算系光周期的 Maxwell/电流沉积样本数；
- 每个场面经完整 Maxwell 步量化后的保守实验室采样间隔和每目标周期样本数。

只有不可避免的 Nyquist 失败为硬错误：网格和 Maxwell 更新都必须严格超过每波长/周期 2 样本；场面也必须严格超过实验室目标周期 2 样本。错误给出限制 `dz`、`dt` 或 `rhythm`。三个警告阈值默认 8，低于只警告，不改变算法。Cowan-z 对已分辨轴向真空模相速精确，正式余量仍应由细化扫描确定。

省略或关闭时不从波荡器/束团猜测频率，不施加频段限制，也不增加高阶采样、场状态、通信或时间循环工作。预检查不要求傍轴载波满足 `dx/dy<wavelength`；横向包络、角接受、孔径、宏粒子噪声、CPML 反射、时间窗和斜向色散仍需问题相关收敛。

## 场边界

边界策略必填。无种子正式测试使用紧凑非分裂 CFS-CPML：

```yaml
boundary:
  type: cpml
  cells: [8, 8, 8]
  polynomial_order: 3.0
  target_reflection: 1.0e-8
  kappa_max: 8.0
  alpha_fraction: 0.0
```

`cells` 是 x/y/z 每面的厚度；0 会关闭该轴两面且不分配历史，非零至少 2。相对两层之间必须留下非 PML 内部区。当前 z 分解还要求完整 z PML 位于端点 rank，非法 rank 数在读粒子前失败。

`target_reflection` 用于计算最大电导率，必须在 `(0,1)`，不保证每个离散模都达到相同反射率。`kappa_max>=1`、正 `polynomial_order`、非负 `alpha_fraction` 控制剖面。示例 `alpha_fraction=0` 在当前近轴前向脉冲测试中反射更低；非零频移保留为低频/倏逝研究的实验选项，必须针对目标谱重新验证。日志报告实际总辅助内存和每 rank 最大值。

无吸收回归使用：

```yaml
boundary:
  type: pec
```

PEC 不分配边界历史。数值构造见[Cowan-z 与 CPML](MAXWELL_COWAN_CPML.md)。

## 初始粒子自场

```yaml
initial_self_field:
  enabled: true
  model: relativistic-poisson
  relative_tolerance: 1.0e-10
  maximum_iterations: 10000
```

块可省略，默认如上。启用后以生产 CIC 形状沉积 boost 后束团，求和 MPI 共享顶点平面，并在第一个物理步前求离散相对论 Poisson。默认 `relativistic-poisson` 由代表质量加权平均轴向速度构造 Vay 刚性束椭圆算符，同时初始化 Yee E 与相应 B。`electrostatic-poisson` 仅保留作 B=0 控制回归。所得 Yee 边 E 用相同电荷沉积检查；不收敛或后残差超容差是硬错误。`maximum_iterations` 必须为正，容差严格位于 `(0,1)`。

CPML 下静态求解零势面是**CPML 内表面**，CPML 内 E/B 从零开始，避免非辐射 Coulomb 尾给零记忆变量充电后回流。若任一 CIC 电荷权重触及该面则停止；应把粒子至少移动到 CPML 内一整个单元，或增大物理 padding，不能静默丢电荷。

标量势和四个分布式 CG 工作薄层只在初始化存在并计入峰值内存，推进前释放，所以传播仍为纯 E/B。静态边界是有限域近似，正式工作需扫描束团到 CPML 距离。刚性模型假定共同轴向速度；速度展宽、包络失配和宏粒子噪声仍需收敛。关闭只用于控制回归并警告，见[初始自场](INITIAL_SELF_FIELD.md)。

## 束团参考与输入

所有粒子位置相对一个实验室束团参考：

```yaml
beam:
  reference:
    input_plane_z: -2.0
    initial_center_z: -1.8
```

粒子输入面、磁元件和探测面共享同一绝对实验室 z；不强制任何元件位于零。HDF5-v3/v4 中 `input_plane_z` 是 Elegant 观察面坐标，`initial_center_z` 是直线向前投影后、共同实验室时刻束团参考坐标。v3 正纵向偏移向下游；v4 保存固定面带符号到达时间，并按每粒子自身速度重构。v1/v2 和生成高斯已是共同时间快照，忽略 `input_plane_z`。

v3/v4 要求每个重构粒子位于 `input_plane_z` 或其下游；随后要求束头加一个等效实验室 boost-z 单元仍位于首个磁元件**作用区**入口前。这给出 `initial_center_z` 合法区间；若为空，错误同时报告两端并建议移动 Elegant 面、调整边缘/作用区、缩短束团或细化 `dz`。

Lorentz 变换把计算系时间零锚定到下游束头事件，其他粒子相对同步，整个变换束团再平移到数值盒中心。相应实验室事件之间可能跨米级，只作为相对同时性诊断，不解释为物理入口漂移，也不把 `initial_center_z` 人为移到更上游。

正式 HDF5 输入：

```yaml
  input:
    type: hdf5
    file: ../examples/particles/example_particles.h5
    electrons: 4.0
    position_offset: [0.0, 0.0, 0.0]
```

相对路径以 YAML 所在目录为准。`electrons` 始终是全部记录代表的物理电子总数。v2/v3/v4 的正 `macro_weight` 会全局归一到 `electrons`，每个宏粒子电荷和质量同比缩放；v1 隐含等权。可选 `position_offset` 用 YAML 长度单位，x/y 平移输入面坐标；v4 的 z 在到达时间同步后平移相对快照。二进制规范见[粒子 HDF5 输入](PARTICLE_INPUT_HDF5.md)。

仅用于小型集成测试的确定性高斯：

```yaml
  input:
    type: generated-gaussian
    electrons: 8.0
    macroparticles: 8
    gamma: 4.0
    direction: [0.0, 0.0, 1.0]
    center: [0.0, 0.0, 0.025]
    sigma_position: [0.02, 0.02, 0.01]
    sigma_proper_velocity: [0.0, 0.0, 0.0]
    random_seed: 17
```

`proper_velocity` 指无量纲 `gamma*v/c`。生成器基于计数器，因此 MPI rank 数变化不改变全局粒子。它只提供无相关高斯，不是束流制备模型。

## 源与磁元件

入射波通过闭合 TF/SF 面注入并由 Maxwell 推进。轮廓有 `plane`、`truncated-plane`、`gaussian`、`super-gaussian` 及相应 `standing-*`。振幅用 `peak_electric_field_V_per_m` 或 `normalized_amplitude`：

```yaml
sources:
  incident_waves:
    - profile: gaussian
      position: [0.0, 0.0, -4.0]
      direction: [0.0, 0.0, 1.0]
      polarization: [1.0, 0.0, 0.0]
      normalized_amplitude: 0.001
      wavelength: 0.8
      radius: [5.0, 5.0]
      order: [0, 0]
      envelope:
        type: gaussian
        center_time: 0.0
        duration: 0.02
        carrier_phase_rad: 0.0
```

包络类型为 `neumann`、`gaussian`、`secant`、`flat-top`、`inverse-gaussian`；适用时可设 `rising_cycles` 和双值 `inverse_gaussian_sigma`。当前 Cowan+CPML 主路线不支持通用注入，这一节只定义已保留接口和 Yee 回归。

磁元件是规定实验室场，只在粒子事件处变换，不复制到 Maxwell 网格。它们可完全省略；无磁元件时保留用户束团中心和 Lorentz 同步，只跳过首磁作用区余量检查，支持纯探测器和无元件自由传播。

```yaml
  magnetic_elements:
    - type: planar-undulator
      characteristic: true
      strength_parameter: 0.1
      period: 10.0
      periods: 1
      entrance_z: 0.0
      polarization_angle_rad: 0.0
      gaussian_fringe: true
      fringe_relative_cutoff: 1.0e-9
```

`uniform-dipole` 改用 `length`、`field_T`。首入口取所有磁铁最小物理 `entrance_z`，不依赖 YAML 顺序，也不强制为零。

平面波荡器保留散度相容解析端场思想，同时把无限支撑改为有限：Gaussian 乘五次紧支撑 taper，横向分量由纵向包络推导，使入口和出口均保持 `div(B)=0`。`fringe_relative_cutoff` 是放置紧支撑边缘的原 Gaussian 值，必须在 `(0,1)`；`gaussian_fringe: false` 时作用区与物理区一致。

`characteristic: true` 是 `undulator_resonance` 使用的可选元数据。多个平面波荡器时必须恰好标一个，选择共振预检查和恒定 boost 推荐所用元件。

## 运行策略

运行策略全局生效，与轨迹/探测器选择独立。可随时停止的小服务器模式：

```yaml
runtime:
  mode: interactive
  stop_check_interval_steps: 16
```

`interactive` 安装最小 SIGINT/SIGTERM 处理器，每隔配置个完整场步由各 rank 协调一次停止标志。停止不会打断 Maxwell 或粒子更新，而是在完整步后退出、提交并关闭全部输出，标记为不完整但可读；轨迹每 `flush_every_samples` 使前缀耐久。即使关闭轨迹也可用。

超算批处理：

```yaml
runtime:
  mode: throughput
  resource_monitor:
    enabled: true
    progress_interval_steps: 1000
    calibration_steps: 1
    memory_safety_factor: 1.25
    time_safety_factor: 1.25
```

`throughput` 为默认，不装信号处理器、不做停止 MPI 轮询、关闭周期轨迹耐久刷新；正常物理停止仍关闭输出，从而避免测试路径同步和文件系统成本。`local-test`、`hpc` 分别是 `interactive`、`throughput` 别名。

`resource_monitor` 独立于停止政策。启用时模型估算场、CPML、粒子、探测器、轨迹和 MPI halo 内存；无入射波时还基准测试 `calibration_steps` 个不改变状态的零场 Maxwell 更新及有界粒子推进/沉积样本。安全系数至少为 1。`progress_interval_steps: 0` 关闭周期报告但保留启动和最终测量；其他值使 rank 0 定期写立即刷新的 `[resource]` 单行，适合 `sbatch` 捕获。时间估算不含文件系统竞争、迁移负载变化和提前停止，正式申请应在目标机器校准。

兼容旧卡时，若无 `runtime.mode` 可用 `trajectory.mode` 作全局别名；新卡应使用前者，同时给出不同值是错误。

## 运行时能量账本

```yaml
energy_ledger:
  enabled: true
  directory: output/example
  filename: energy-ledger.h5
  sample_interval_steps: 10
  buffer_records: 64
  compression: 0
  warning_relative_tolerance: 0.01
```

它独立于探测器/轨迹，记录计算系粒子动能、物理区 E/B、六个 CPML 内表面通量、移除粒子动能和规定器件做功。完整体场/粒子统计只按 `sample_interval_steps` 计算，边界功率和规定做功每场步累积。只有 rank 0 写小 HDF5；省略或关闭完全移除工作与存储。

同一记录含实验室 gamma 均值/rms、线性纵向 chirp 和去 chirp rms；它们位于等计算系时间片，不能与计算系能量项相加。严格加速器入口/出口或切片诊断用固定实验室粒子面。方程和 schema 见[能量账本](ENERGY_LEDGER.md)。

## 停止策略

除最大保护 `mesh.duration` 外必须选择一个物理停止；在物理条件前耗尽 duration 是错误，并使轨迹输出不完整。

所有有效粒子越过最后元件作用区出口：

```yaml
stop:
  mode: after-last-element
```

元件可为磁器件、场面或粒子面，故纯漂移中探测器可成为唯一元件。若完全没有已启用元件，该模式在读卡时即拒绝，因为没有有限目标；应使用 `reference-center-z`，避免作业只运行到 duration 才晚失败。

CPML 内表面是物理粒子边界。轨迹和探测器参与在第一次 CPML 入口终止，余下场步由不输出弹道载流子完成并按 CPML 电导率衰减；外表面以虚拟向外电流导出残余电荷。真实粒子立即从 all-particles 谓词移除，载流子不延迟停止。若无有效粒子，模式停止并明确报告，而不声称束团越过下游边界。若需粒子离开后的场振铃，使用 `reference-center-z`。

该粒子边界没有 YAML 开关；无 CPML 面使用直接电荷守恒外边界。正式运行不分配面大小粒子通量数组，也无边界 I/O，结束报告 CPML 入口、直接逃逸和残余载流子，详见[粒子开放边界](PARTICLE_OPEN_BOUNDARY.md)。

参考中心到固定实验室 z：

```yaml
stop:
  mode: reference-center-z
  z: 2.5
```

参考中心从 `initial_center_z` 出发，沿计算系原点世界线运动。目标必须严格位于初始中心及所有当前元件作用区下游，允许器件后观察或完全无元件传播。

停止统一使用束线元件范围。磁器件和两类实验室探测面都参加首末边界；场面停止边界是采样面，不是左诊断区入口。粒子面不受元件重叠规则约束，也可位于所有磁器件之前。

## 实验性粒子退休

退休层是场面弹道本底重构的替代方案：

```yaml
particle_retirement:
  enabled: true
  entrance_z: 5.0
  length: 2.0
```

一个场面可绑定手动频谱保护：

```yaml
      retirement_frequency_protection:
        enabled: true
        minimum_photon_energy_eV: 50.0
        cycles: 2
```

位置使用 YAML 长度单位且固定在实验室系。粒子到 `entrance_z` 后停止真实推进、普通轨迹节拍和探测器过面；紧凑弹道载流子沿瞬时速度继续，并以

```text
w(s) = 1 - 10 s^3 + 15 s^4 - 6 s^5,   s in [0,1]
```

衰减电流，在 `entrance_z+length` 移除。只影响粒子电流，已有 E/B 和求解器不变，无右侧诊断区。

该方法明确为实验性。taper 降低立即删除的高频冲击，但域内净电荷移除不严格连续；必须使用匹配零辐射基线，并在 E/B 振幅层面相减，见[退休功率比较](FIELD_POWER_COMPARISON.md)。

保护不会选择或修改 length。读粒子后以最大实验室 gamma 检查

```text
beta_max = sqrt(1 - 1/gamma_max^2)
delay = length/c * (1/beta_max - 1)
observed_cycles = minimum_photon_energy_eV * delay / h
observed_cycles >= cycles.
```

光子能量是必须保护频段的下边缘，不一定是波荡器中心。场面还保留

```text
extent_x = cells_x * cell_size_x
extent_y = cells_y * cell_size_y
rho_guard = sqrt(extent_x^2 + extent_y^2)
causal_distance = gamma_max * rho_guard
causal_z = detector_z - causal_distance.
```

退休出口加一个等效实验室 z 单元必须位于 `causal_z` 左侧；无效长度/位置在场分配前退出并建议长度或场面 z。受保护场面必须是最后束线元件，任何磁器件、场面或粒子面都不能在其下游。场面区可互相重叠，粒子面仍不参加区域排斥；保护区与磁作用区按弹道参考相同规则排斥。

只有一个场面可绑定全局退休，且其 `particle_background.enabled` 必须为 false。弹道参考是独立替代，不改变真实粒子或 Maxwell，因此对其采样面下游元件无此限制。

初始化还强制：入口至少比最后磁作用区（含边缘支撑）靠后一个等效实验室 z 单元；所有场面至少比 taper 出口靠后一个单元；频率保护时长度覆盖请求周期；出口在因果区上游且受保护场面为最后元件；至少存在一个场面；停止必须为 `reference-center-z` 且位于最后场面下游。粒子面可位于退休入口作小规模诊断；其他场面可重叠但都要位于退休出口后。关闭退休不分配载流子向量，也无活动载流子工作。

## 实验室探测面

程序只接受固定实验室场面和粒子过面；旧程序的计算系追踪面与移动诊断面被明确废弃。

```yaml
detectors:
  enabled: true
  directory: output/example/detectors
  field_planes:
    - name: exit-fields
      z: 2.5
      rhythm: 0.002
      buffer_samples: 2
      compression: 0
      particle_background:
        enabled: true
        buffer_records: 16384
        compression: 0
        validation:
          enabled: false
          maximum_particles: 100000
  particle_planes:
    - name: exit-particles
      z: 2.5
      buffer_records: 4096
      compression: 0
```

`particle_background` 记录虚拟直线参考供离线解析扣除；`retirement_frequency_protection` 审计手动全局退休层，两者不能在同一场面或同一退休运行同时开启。`particle_background` 默认关闭，必须显式开启；省略时不创建 companion、过面事件、缓冲、通信或参考区。

粒子面是实验室 z 处零长度元件。场面的物理和停止位置也是 z；仅当显式开启粒子本底时，它拥有紧邻左侧的诊断作用区。任一种面都可成为 `after-last-element` 最后元件；`reference-center-z` 必须位于采样面下游。名字全局唯一，只能含字母、数字、`.`、`-`、`_`。

读完实验室粒子后，场面推导

```text
extent_x = cells_x * cell_size_x
extent_y = cells_y * cell_size_y
rho_guard = sqrt(extent_x^2 + extent_y^2)
gamma_guard = maximum initial laboratory particle gamma
reference_length = gamma_guard * rho_guard
reference_z = detector_z - reference_length.
```

使用完整横向宽度对角线是保守规则。真实粒子向下游过 `reference_z` 时记录一次实验室位置、固有速度、电荷、质量和权重，定义虚拟直线；不改变推进、电流、MPI 或 Maxwell。主场文件仍为原始总场，companion 用于分析机卷积/扣除。

`particle_background.validation` 只用于测试：粒子在场面再采样一次，rank 0 按 ID 配对并报告横向位置、到达时间、固有速度和方向误差。`maximum_particles` 是强制宏粒子数保护；超过即退出，且验证要求本底已开启。关闭验证时无配对表、验证数据集和额外场面粒子事件。

独立排斥规则：磁作用区互不重叠；磁作用区不与场面左参考区重叠；场面参考区之间可重叠；粒子面不参加。非法放置在场分配前退出并建议最小下游 z，详见[场面弹道参考区](FIELD_DETECTOR_REFERENCE.md)。

场面在所有 x-y 单元中心存实验室 E/B。`rhythm` 是实验室最小间隔，实际时刻因完整 Maxwell 步量化而另存。`buffer_samples` 控制 rank-0 HDF5 批次，通常应小，因为每样本是完整 x-y 面。Yee B 仍错开半步，元数据明确记录，不用第二份全域场隐藏。

粒子面只在粒子线段实际穿过固定 z 时保存记录，含方向；不改变本底参考、轨迹节拍或轨迹记录。`buffer_records` 只控制独立粒子面写入器。

省略 `detectors`、列表为空或 `enabled: false` 时，不构造管理器、不分配缓冲、不打开文件、不进入探测 MPI。启用时只有当前拥有场面的 rank 采样，粒子事件有界批量发送，rank 0 独占追加 HDF5，详见[探测器输出](DETECTOR_OUTPUT_HDF5.md)。

## 轨迹输出

```yaml
trajectory:
  enabled: false
  directory: output/example
  basename: particles
  rhythm: 0.002
  buffer_records: 64
  flush_every_samples: 1
  compression: 0
```

轨迹用于少粒子验证而非正式辐射主路径。关闭时不打开文件、不保留记录缓冲、不周期采样，也跳过终端记录转换。开启时每 MPI rank 一个 HDF5；`rhythm` 使用 YAML 时间单位。全局交互模式每 `flush_every_samples` 提交可读前缀，throughput 只在批次满或正常关闭时刷新。`compression` 范围 0–9，0 最省 CPU。关闭时可省略 `rhythm`。格式 v2 额外记录节拍间精确 CPML 入口或直接域逃逸；数值载流子从不写入。见[轨迹 HDF5](TRAJECTORY_OUTPUT_HDF5.md)。
