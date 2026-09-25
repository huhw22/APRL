# 粒子本底场重构

## 目的

`field_reconstruction` 是离线 C++/FFTW 工具，读取：

1. 一个原始实验室系场面 HDF5；
2. 该场面的 `*-ballistic-reference.h5`，或位于完全相同 z 的实验室系粒子面；
3. 一份 YAML 后处理卡。

工具写出新的 HDF5，其中包含均匀时间采样、已扣除粒子本底的 E/B 面，以及原始/清洁 Poynting 诊断；输入文件不会改变。它独立于模拟器，设计为在低成本分析服务器上运行。

## 模型与规模

工具不会对每个粒子、像素和时刻逐项求和。每次向下游过面只用 cloud-in-cell 权重沉积一次到 `(t,y,x)` 网格；三个实源网格经 FFTW 变换后，通过解析 Maxwell 传递函数得到五个本底场分量。因此主复杂度为

```text
O(number of crossing records) + O(N log N),  N = Nt * Ny * Nx,
```

而不是 `O(particles*samples*pixels)`。FFTW 线程数由 YAML 配置。内存正比于补零后的场立方体，而非粒子数；启动日志与输出属性公开全部补零和模型参数。

纵向速度为 `v_z` 的粒子穿过观察面时，匀速解的傅里叶空间标量分母为

```text
D = kx^2 + ky^2 + omega^2/(gamma^2 v_z^2).
```

对窄能散、主要纵向传播的束团，实现使用电荷加权平均 `1/(gamma^2 v_z^2)` 构造 `D`，并分别沉积与 `q`、`q/v_z`、`q/(gamma^2 v_z^2)` 成正比的源作为分子。所得场是 Maxwell 方程的匀速解，不是把静电场粘贴到平面上。

它仍是模型，不是 PIC 历史精确回放。单次过面记录无法恢复先前加速场、网格色散、CPML 历史或初始场误差；这些会有意留在清洁场中，作为辐射或数值内容。以下硬保护使近似边界显式化：

- 最大相对 gamma 展宽；
- 最大横向 beta 均方根；
- 落在补零 FFT 立方体之外的最大电荷比例；
- 用户控制的横向平滑和 t/x/y 补零。

任一保护超限即退出。必须扫描补零和平滑，直到目标辐射量收敛。

## 输入卡

```yaml
input:
  field_file: ../../output/run/detectors/radiation.h5
  particle_file: ../../output/run/detectors/radiation-ballistic-reference.h5
  require_complete: true
  particle_read_chunk: 65536

model:
  fft_threads: 8
  padding_factor: [2, 2, 2]   # [laboratory time, y, x]
  transverse_smoothing_m: 0.0
  maximum_relative_gamma_spread: 0.05
  maximum_rms_transverse_beta: 0.05
  maximum_outside_charge_fraction: 0.001

output:
  file: reconstructed-fields.h5
  compression: 0
  overwrite: false
```

路径相对输入卡解析。弹道参考记录在沉积前从左入口直线传播至场面；`/particle_plane` 输入必须与场面同位置。两类输入都只使用每个粒子 ID 的第一次合法下游过面。模拟器粒子 ID 是稠密、从 1 开始的，因此重复检测用精确位图；`10^7` 粒子约需 1.25 MB，而哈希表通常需数百 MB。重复和无效记录都会计数。

场面时间由 Maxwell 步量化，通常不严格均匀。工具先把原始 E/B 线性重采样到覆盖相同 committed 时间范围的均匀网格，再执行 FFT，并把新时间轴和重采样契约写入输出。

## 输出

`/reconstructed_field` 包含：

- `time_s`，`[samples]`，均匀实验室时间；
- `electric_V_per_m`，`[samples,ny,nx,3]`；
- `magnetic_T`，`[samples,ny,nx,3]`；
- `raw_signed_power_W`、`cleaned_signed_power_W`；
- `raw_forward_power_W`、`cleaned_forward_power_W`；
- `complete`，仅在全部场分量和诊断写完后置 1。

带符号功率是 `(E cross B)_z/mu0` 在横向的积分；前向功率使用 `max((E cross B)_z/mu0,0)`，二者都在存储的单元中心求积。`raw_signed_energy_J` 和 `cleaned_signed_energy_J` 是带符号平面通量的梯形时间积分，用于实验室系入口/出口控制体账本。

`raw_forward_energy_J`、`cleaned_forward_energy_J` 和 `relative_forward_energy_change` 给出前向时间积分及相对变化，直接度量所选孔径和时间窗内的强度污染；它本身不是相对解析波荡器谱的误差条。

原始探测器在 leapfrog 半步时刻保存 B；当前工具保留该错位，而解析本底按标记样本时刻计算。该残差应随主时间步缩小，必须纳入收敛测试。

## 实际使用

正式计算优先使用弹道参考文件，因为每个粒子只需一个紧凑记录。双面验证只对小束团开启，以验证直线延续。之后可用同位置粒子面替换重构输入；两种清洁输出应在轨道验证误差内一致。

大文件先用补零 `[1,1,1]` 估计内存和墙钟时间，再提高因子。补零降低周期卷绕，却按三个因子的乘积放大 FFT 内存。压缩减少磁盘但增加分析 CPU 时间。
