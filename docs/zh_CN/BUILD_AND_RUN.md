# 构建、粒子转换与运行

## 依赖

- CMake 3.16 或更新；
- 支持 C++11 的编译器；
- MPI；
- HDF5 C 1.10 或更新；多 rank 粒子输入必须使用 parallel HDF5；
- yaml-cpp 0.6 或更新；
- 构建原生 Elegant SDDS 转换器时需要官方 SDDS C 库/工具包；
- 场重构和场面频谱/相干性分析需要 FFTW3 及其线程库。

在新机器上或更换编译器/MPI module 后，应先运行编译型依赖体检。完整安装、离线集群和卸载流程见[安装、卸载与依赖维护](../INSTALLATION.md)。

## 构建

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_COMPILER=mpic++
cmake --build build -j
```

主构建默认优先并强制 parallel HDF5。严格单 rank 的本地版本可使用 `-DAPRL_REQUIRE_PARALLEL_HDF5=OFF`，但所得程序不能用多 MPI rank 读取粒子 HDF5。

构建得到：

- `build/aprl`；
- `build/particle_text_to_hdf5`；
- 找到 SDDS 时的 `build/elegant_sdds_to_hdf5`；
- `build/undulator_resonance`；
- `build/energy_ledger_report`；
- `build/lab_frame_energy_estimate`；
- 默认还会统一构建五个后处理器。

正式运行前执行轻量发布门禁：

```bash
cmake --build build --target verify_required
```

它覆盖启动、Lorentz/Elegant、Boris、电荷连续性、Cowan、CPML、MPI、探测器/停止、HDF5 兼容、K=0 账本以及端到端后处理。非门禁的物理能量示例单独保存。测试矩阵和解释见[分层测试](../TESTING.md)。

SDDS 转换器有意链接官方实现，而不是不完整的自写二进制解析器。若 SDDS/Elegant 使用推荐的同级目录布局，CMake 会自动发现；否则：

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_COMPILER=mpic++ \
  -DSDDS_ROOT=/path/to/built/SDDS
```

## 随主程序构建的诊断工具

总结完整或协调中断的运行时能量账本：

```bash
./build/energy_ledger_report output/run-name/energy-ledger.h5
```

可选第二个参数选择从 0 开始的已提交记录，同时报告初始化至该样本的变化和最终摘要。工具只读取 committed 前缀，打印初末账本、场能比例、平均 gamma 与能散分解，见[运行时能量账本](ENERGY_LEDGER.md)。

不加载模拟粒子或场，生成实验室系解析尺度估计：

```bash
./build/lab_frame_energy_estimate config/lab_frame_energy_estimate.yaml
```

示例给出归一化 50 A 束缚场与波荡器辐射尺度，见[实验室系能量诊断](LAB_FRAME_ENERGY_DIAGNOSTICS.md)。

选择退休层保护频段前，检查特征平面波荡器共振：

```bash
./build/undulator_resonance config/generated_gaussian.yaml
./build/undulator_resonance config/example.yaml 1000
```

当 HDF5 束流卡没有 `beam.input.gamma` 时，第二个参数提供实验室 gamma。只有一个平面波荡器时自动选择；存在多个时必须恰好把一个设为 `characteristic: true`。工具报告轴上冷束基频的一维增益中心估计，并给出

```text
boost_gamma_recommended = gamma_lab/sqrt(1+K^2/2),
```

与 `mesh.boost_gamma` 比较，再把速度失配换算为穿过特征波荡器核心时的计算系 z 漂移（米和网格数）。推荐值只选择一个恒定惯性系，运行中绝不改变 boost。自由漂移、边缘场、能散、发射度和集体演化仍需预留盒子余量。

如果场面启用 `retirement_frequency_protection`，同一工具会读取最低光子能量和周期数，报告手动退休长度要求及预检查通过/失败。工具从不修改输入卡。

## 独立后处理器

轨迹辐射工具可独立构建：

```bash
cmake -S postprocess/trajectory_radiation -B build-radiation \
  -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=mpic++
cmake --build build-radiation -j
```

它读取每 rank 轨迹，重组迁移粒子历史，由 rank 0 写一个 HDF5，见[轨迹远场工具](TRAJECTORY_RADIATION.md)。

原始场/粒子参考重构：

```bash
cmake -S postprocess/field_reconstruction -B build-field-reconstruction \
  -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=mpic++
cmake --build build-field-reconstruction -j
./build-field-reconstruction/field_reconstruction \
  postprocess/field_reconstruction/example.yaml
```

它在一个分析节点上使用多线程 FFTW，不修改模拟输出。输入假设、输出和收敛要求见[粒子本底场重构](FIELD_RECONSTRUCTION.md)。

匹配零辐射基线功率比较：

```bash
cmake -S postprocess/field_power_compare -B build-field-power-compare \
  -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=mpic++
cmake --build build-field-power-compare -j
./build-field-power-compare/field_power_compare \
  postprocess/field_power_compare/example.yaml
```

它从两个匹配场面流式读取小型横向批次，进行场振幅相减与频段滤波，并可与单次轨迹远场比较峰值功率和能量，见[粒子退休与功率比较](FIELD_POWER_COMPARISON.md)。

完整场面频谱/相干性分析：

```bash
cmake -S postprocess/field_plane_analysis -B build-field-plane-analysis \
  -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=mpic++
cmake --build build-field-plane-analysis -j
./build-field-plane-analysis/field_plane_analysis \
  postprocess/field_plane_analysis/example.yaml
```

它接受原始探测场、扣除粒子本底的重构场，或匹配信号/基线对；创建输出前检查工作内存和输出大小，见[场面频谱与相干性分析](FIELD_PLANE_ANALYSIS.md)。

粒子/辐射能量闭合报告：

```bash
cmake -S postprocess/energy_closure -B build-energy-closure \
  -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=mpic++
cmake --build build-energy-closure -j
./build-energy-closure/energy_closure \
  postprocess/energy_closure/example.yaml
```

它按 ID 连接两个粒子面，稳定计算高 gamma 动能差，可逐粒子扣除匹配 `K=0`，并与场面或轨迹远场能量比较，见[粒子/场能量闭合](ENERGY_CLOSURE.md)。

## 转换粒子输入

对 Elegant 原生输出，优先使用直接转换器：

```bash
./build/elegant_sdds_to_hdf5 \
  --input elegant-watch.sdds \
  --page 1 \
  --output particles.h5
```

它通过官方库读取一页中的 `x,xp,y,yp,t,p`，存在时保留 `particleID`，写 HDF5 v4 固定面过面事件，不产生中间文本。只有真正的逐行正相对权重列才能通过 `--weight-column NAME` 选择；Elegant 页参数 `Charge` 不是这种列。精确映射和多页规则见[粒子 HDF5 输入规范](PARTICLE_INPUT_HDF5.md)。

对旧六列文本导出：

```bash
./build/particle_text_to_hdf5 \
  --input particles.txt \
  --output particles.h5 \
  --length-unit micrometer
```

文本转换器扫描输入两次并分块写入，内存不随完整粒子数增长。HDF5-v3 文本列为 `x_plane y_plane zeta ux uy uz [macro_weight]`：x/y 位于固定 Elegant 面，`zeta` 是重构共同时间束团的带符号纵向偏移；可选第七列为正相对权重，默认 1。接受空行和 `#` 注释。`beam.reference.input_plane_z` 与 `initial_center_z` 必须位于同一实验室坐标系。正式运行只集体读取转换后的 HDF5，每 rank 一个连续区间。

重建项目自带示例：

```bash
./build/particle_text_to_hdf5 \
  --input examples/particles/example_particles.txt \
  --output examples/particles/example_particles.h5 \
  --length-unit micrometer
```

## 运行

```bash
./build/aprl config/example.yaml
mpirun -n 4 ./build/aprl config/example.yaml
```

`config/generated_gaussian.yaml` 是无需粒子文件的测试输入。小服务器调试选择全局 `runtime.mode: interactive`；不间断批处理选择 `throughput`。这是全局策略，所以即使关闭轨迹、只开启探测面，交互模式信号停止仍然有效。

示例卡开启 rank-0 资源报告。时间循环前打印每 rank 峰值内存估计、聚合模型内存、未压缩输出上限、校准秒/步和墙钟时间上限；长运行按间隔打印紧凑进度，结束时打印实测当前/峰值驻留内存及循环/总墙钟时间。每条记录是带 `[resource]` 前缀并立即刷新的单行，普通批处理重定向即可：

```bash
sbatch --output=run-%j.log run.sh
# run.sh 最终例如执行：
srun ./build/aprl config/production.yaml
```

启动计时只是本地短微基准，不保证调度时间；并行文件系统竞争、后续粒子迁移和提前物理停止均不在预测内。申请大型资源前，应在目标机器和 MPI 分解上校准内存/时间安全系数。

紧凑 I/O 示例开启一个实验室场面和一个粒子面，但因小盒子没有足够自由漂移，显式关闭弹道参考。文件位于 `output/example/detectors/`，只有 rank 0 写入。无需探测器的运行可删除 `detectors` 或设 `detectors.enabled: false`，此时不构造缓冲和通信。格式见[探测器 HDF5 输出](DETECTOR_OUTPUT_HDF5.md)。

只有需要直线粒子本底重构的正式场面才应显式开启 `particle_background`。启动时由横向网格对角线和最大输入粒子实验室 gamma 计算左参考边界，若与磁作用区重叠则拒绝并建议修正场面位置，见[弹道参考区](FIELD_DETECTOR_REFERENCE.md)。小规模验证可再启用 `particle_background.validation.enabled`；它只记录比较真实轨道和虚拟直线所需的两次过面，不要求完整轨迹，并由强制粒子数上限防止误用于生产束团。

若变换后的初始束团不适合纵向盒子，求解器会明确退出。HDF5-v3/v4 还检查每条记录是否从配置 Elegant 面向前投影。有磁元件时，重构束团还必须位于第一个磁作用区前；纯探测器或无元件卡不会凭空引入该约束。头部锚定的 Lorentz 同步跨度单独报告，不被当作真实入口漂移距离。

示例使用 `mesh.field_solver: cowan-z`，要求 `dx>=dz`、`dy>=dz`；不合法时在分配场和读取粒子前退出，并给出横向网格建议。日志记录 Cowan 系数与横向真空色散诊断。紧凑 CFS-CPML 会报告实际辅助内存。当前 CPML 与非空 TF/SF `incident_waves` 不兼容，Cowan TF/SF 需等待注入面推广，详见[Cowan-z 与 CPML](MAXWELL_COWAN_CPML.md)。真实粒子轨迹终止于 CPML 内表面，随后仅以紧凑、不输出的载流子继续并匹配衰减电流；外表面执行电荷守恒残余清理，见[粒子开放边界](PARTICLE_OPEN_BOUNDARY.md)。轨迹 v2 结构及精确 CPML 入口事件见[轨迹 HDF5 输出](TRAJECTORY_OUTPUT_HDF5.md)。
