# 分层测试与数值分析流程

测试体系分为软件回归与物理数值分析。软件回归用于确认程序接口和数值路径保持一致；物理数值分析用于评估特定束团、频段和网格下的收敛性。物理分析不采用跨问题通用的通过阈值。

## 一次运行全部必需测试

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_COMPILER=mpic++ -DBUILD_TESTING=ON
cmake --build build -j
cmake --build build --target verify_required
```

也可以直接使用 CTest：

```bash
ctest --test-dir build --output-on-failure -L required
```

所有输入卡只使用很小的网格和少量宏粒子。它们是运行前回归门禁，不是实际 EUV/FEL 数值实验，也不替代生产参数下的网格、宏粒子、时间窗和吸收层收敛研究。测试输出只写入 `build/verification/work/`，不会污染源代码目录。

## 四个层级

### 1. 运行测试（必需）

`required.runtime.smoke` 从 YAML 装载开始，完成一次无外场、无磁元件的自由传播，验证主程序能够初始化、推进并由实验室系参考中心停止条件正常结束。

### 2. 功能测试（必需）

| 预定功能 | 自动测试 | 检查内容 |
|---|---|---|
| Lorentz 变换和 Elegant 重构 | `required.lorentz-elegant` | 时空、动量和 E/B 往返；HDF5-v3 纵向偏移重构；HDF5-v4 定平面到达时刻重构；MPI 权重质心 |
| 单粒子 Boris 相位收敛 | `required.boris-phase` | 均匀磁场回旋相位的二阶收敛和 gamma 守恒 |
| 离散电荷连续性 | `required.charge-continuity` | 跨多个网格的轨迹分段及离散连续性残差 |
| Cowan 色散 | `required.cowan-dispersion` | z 轴真空波的无色散传播、平滑系数归一化、解析相速和群速 |
| CPML 反射 | `required.cpml-reflection` | 同一脉冲的 PEC 对照与 CPML 反射幅度抑制 |
| MPI 单核/多核一致性 | `required.mpi.*` | 1 rank 和 2 ranks 的粒子探测面记录按源粒子 ID 逐项比较 |
| 探测器节拍和停止条件 | `required.detector-stop.*`、`required.stop.*` | 场采样间隔、粒子过面、完成标记、末元件停止，以及空束线错误配置拒绝 |
| HDF5 各版本 | `required.hdf5-versions` | 粒子输入 v1-v4、权重兼容、坐标语义及 MPI hyperslab 完整覆盖 |
| 能量账本和 K=0 基线 | `required.energy-ledger-k0.*` | K=0 账本结构、提交前缀、完成标记、有限值、时序及零解析外场功 |

另有 `required.yaml.unknown-key*`，保证主程序以及五个后处理器的输入卡拼写错误不会被静默忽略。`required.output-safety.*` 验证默认拒绝覆盖，而不是只检查配置解析。

完整配置注册 37 项必需测试。`required.postprocess.*` 使用微型模拟输出，分别验证原始场探测面直接进入频谱/相干性分析，以及自由区直线参考经过近场重构后进入同一分析工具。两条路线均检查实际读取的 HDF5 类型，不以程序退出码代替数据验证。测试同时覆盖少量粒子轨迹辐射、信号/零辐射基线功率比较和实验室系能量闭合报告，并检查 HDF5 完成标记、关键数值有限性、报告字段和后处理覆盖保护。主程序与后处理之间的格式接口因此属于发布门禁。探测器测试还核对场与粒子文件的 `run_id` 一致性，并确认运行清单保存原输入卡。

这里的 `K=0` 必需测试只验证账本的数据通路和一个宽松的基线守恒量级。即使关闭初始自场，运动电荷也会在推进后沉积电流并建立网格场，所以它不是“永远没有场”的测试。它不声称有限电荷束团的近场、动能和辐射能量已经在任意生产参数下收敛。

### 3. 物理自洽分析（可选、非门禁）

配置时显式开启：

```bash
cmake -S . -B build -DAPRL_ENABLE_PHYSICS_ANALYSIS=ON
cmake --build build --target verify_physics
```

该层运行一对小型 `K=0`/非零 `K` 示例，生成完整运行时能量账本，并调用 `energy_ledger_report` 输出粒子动能、内部场能、六面净通量和解析外场做功的变化。CTest 只要求模拟和报告工具成功完成，不对闭合误差设置通过/失败阈值；结果用于观察数值尺度、符号和收敛趋势。

对正式问题，应在固定物理束团下分别扫描网格、CPML 厚度、横向孔径、时间窗、宏粒子数和粒子子步进，再比较：

```text
Δ(有效粒子动能 + 已移除粒子动能 + 计算域内场能)
+ 六面向外场能通量
- 解析外场做功
= 数值余额
```

单个前向场探测面只覆盖部分通量，不能单独当作全局能量守恒的右端。

### 4. 实验或学术比较（未纳入自动测试）

`verification/academic/` 保留为空，且不注册为 CTest。新增论文或实验对照时，测试定义必须声明参考数据、单位、观测量、误差定义和适用参数。外部结果在完成适用性审查前不得纳入核心程序硬门禁。

## 选择性运行

```bash
ctest --test-dir build -N
ctest --test-dir build --output-on-failure -R 'required\.cpml-reflection'
ctest --test-dir build --output-on-failure -R 'required\.mpi\.'
ctest --test-dir build --output-on-failure -L smoke
ctest --test-dir build --output-on-failure -L functional
```

必需测试失败的构建不得用于生产运行。可选物理分析中的数值偏差仅作为诊断材料，其物理可接受性必须结合具体问题的参数收敛结果判定。
