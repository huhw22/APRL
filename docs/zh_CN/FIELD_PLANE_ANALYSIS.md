# 场探测面频谱与相干性分析

## 目的与范围

`field_plane_analysis` 是面向正式大宏粒子计算的固定实验室场面后处理器。它在分析节点上独立使用多线程 FFTW 与串行 HDF5，不需要粒子轨迹。输入可以是：

- 模拟器原始 `/field_plane`；
- 已扣除粒子本底的 `/reconstructed_field`；
- 上述任一种信号场加几何和时间完全相同的零辐射基线，在电场振幅层面相减。

输出正频率前向角谱、积分能谱、偏振/Stokes、相干/涨落分解、选定角度互谱密度（CSD）和选定双频 CSD。轨迹辐射只保留为小粒子验证路线。

## 辐射模型

对每个保留的时间 FFT 频点，工具在探测孔径上对实验室电场做二维傅里叶变换：

```text
kx = 2 pi qx/(Nx dx),  ky = 2 pi qy/(Ny dy),  k = 2 pi f/c,
nx = kx/k, ny = ky/k, nz = sqrt(1-nx^2-ny^2).
```

只有 `nx^2+ny^2<1` 是前向传播模。复电场投影到与 `n=(nx,ny,nz)` 正交的局部水平/垂直基。连续变换

```text
A(kx,ky,f) = integral E(x,y,t) exp[-i(kx x + ky y + 2 pi f t)] dx dy dt
```

对应普通正频率角谱能量密度

```text
d^2 W/(dE_gamma dOmega)
  = 2 epsilon0 c k^2 nz^2 |A_transverse|^2 / ((2 pi)^2 h_eV).
```

Nyquist 自共轭频点的单边因子为 1 而非 2。一个 `nz` 来自法向 Poynting 投影，另一个来自 `dkx dky=k^2 nz dOmega`。离散立体角权重为 `dkx dky/(k^2 nz)`，密度乘权重求和得到通过记录平面的已接受前向模能量。

这是下游真空中的前向辐射角谱模型，不是通用近场分解。工具故意不使用 Yee B，因为探测 B 错开半步，且重构文件的解析粒子 B 按标记时刻计算。纵向/倏逝内容和后向波被排除。场面应放在真空区，并重复孔径、横向分辨率、探测距离和时间步收敛。

## 输入卡

```yaml
input:
  field_file: ../../output/run/detectors/radiation.h5
  # 退休路线可选：
  zero_radiation_baseline: ../../output/baseline/detectors/radiation.h5
  require_complete: true

analysis:
  photon_energy_band_eV: [50.0, 100.0]
  transverse_window: none
  angular_zero_padding: [2, 2]  # [y,x]
  time_windows:
    enabled: true
    interval_s: [5.0e-15, 25.0e-15]
    duration_s: 4.0e-15
    step_s: 2.0e-15

calculation:
  frequency_block: 8
  spatial_batch_points: 32
  fft_threads: 8
  maximum_working_mib: 4096
  maximum_output_gib: 64.0

coherence:
  spatial_photon_energy_eV: [75.0]
  spatial_reference_angles_rad: [[0.0, 0.0]]
  temporal_photon_energy_eV: [65.0, 70.0, 75.0, 80.0, 85.0]
  temporal_reference_angles_rad: [[0.0, 0.0]]

output:
  file: field-plane-analysis.h5
  compression: 0
  overwrite: false
```

路径相对分析卡解析。原始场和基线必须有相同 committed 时间与几何。量化的探测器时间会线性重采样到均匀网格，并报告以一个时间步为单位的最大修正。目标频段超过所得 Nyquist 光子能量是硬错误。

`angular_zero_padding` 只插值 FFT k 网格，不增加孔径或物理角分辨率。`transverse_window: hann` 可压制硬孔径边缘，但改变已接受场且不做能量修正；应与 `none` 对比，不能只选择视觉更平滑者。

`frequency_block` 用内存换取较少的重复读场和时间 FFT；`spatial_batch_points` 限制每次 HDF5 读取。工具在创建输出前估算峰值工作内存和 HDF5 大小，两个用户上限都是硬保护。

## 时间窗系综

`time_windows.enabled: false` 时，完整记录就是一个矩形样本，其 CSD 必然秩一，不能说明部分相干。

开启后，`interval_s` 内每个完整 Hann 窗作为一个系综样本。持续时间和步长会量化到探测器样本。窗振幅乘 `sqrt(N/sum(w^2))`，使平稳信号的积分谱能量无偏；窄频段仍可能漏掉 Hann 旁瓣。重叠窗互相相关，不等价于同数量独立 shots。

只有检查 `sample_energy_spectrum_J_per_eV` 与 `sample_band_energy_J` 后才能选择区间。启动、饱和或脉冲衰减落在区间内会表现为相干性降低；这可能是物理问题，但不是平稳时间平均。

## 输出

根组为 `/field_plane_analysis`；全部产品刷新后 `complete` 才置 1。主坐标包括 `photon_energy_eV`、`frequency_Hz`、`omega_rad_per_s`、`wavelength_m`、`kx_rad_per_m`、`ky_rad_per_m`。每个频率的角度定义保存在 `angular_coordinates` 属性；`solid_angle_weight_sr=0` 表示倏逝或非前向频点。

| 数据集 | 形状 | 含义 |
|---|---|---|
| `mean_angular_electric_field_spectral` | `[f,ky,kx,2]` complex | 系综平均水平/垂直孔径谱，`V s m` |
| `mean_angular_spectral_energy_density_J_per_eV_sr` | `[f,ky,kx]` | 平均角强度 |
| `coherent_angular_spectral_energy_density_J_per_eV_sr` | `[f,ky,kx]` | 系综平均复场的强度 |
| `fluctuation_angular_spectral_energy_density_J_per_eV_sr` | `[f,ky,kx]` | 平均减相干强度 |
| `mean_stokes_spectral_energy_density_J_per_eV_sr` | `[f,ky,kx,4]` | 约定下的 I,Q,U,V |
| `solid_angle_weight_sr` | `[f,ky,kx]` | 离散 `dOmega` 权重 |
| `mean_energy_spectrum_J_per_eV` | `[f]` | 角积分能谱 |
| `sample_energy_spectrum_J_per_eV` | `[sample,f]` | 平稳性/shot 诊断 |
| `coherent_fraction_spectrum` | `[f]` | 角积分的平均场强度/平均强度 |
| `global_degree_of_transverse_coherence` | `[f]` | 角度与偏振空间的基无关 `Tr(W^2)/Tr(W)^2` |
| `ensemble_gram_matrix_J_per_eV` | `[f,sample,sample]` complex | 与完整角 CSD 有相同非零本征值的加权 Gram 矩阵 |
| `sample_band_energy_J` | `[sample]` | 保留频段的矩形频点求和 |

`/spatial_coherence` 保存选定能量的 `angular_cross_spectral_density_J_per_eV_sr[energy,reference,ky,kx,2,2]` 和基不变的 `spectral_degree_of_coherence_squared`。每个频率分别把请求参考角映射到最近传播 FFT 方向，并保存实际角度与索引。

若 `a_s` 是样本 s 经立体角和偏振加权的角场，则

```text
G_st = <a_s,a_t>/N_samples.
```

`G` 与完整 CSD 有相同非零本征值，所以 `Tr(G^2)/Tr(G)^2` 是全局横向相干度，对较小样本矩阵对角化即可得到相干模权重。它不同于只度量系综平均场、且对共同载波相位抖动敏感的 `coherent_fraction_spectrum`。

`/temporal_coherence` 保存 `two_frequency_cross_spectral_density_J_per_eV_sr[reference,f1,f2,2,2]` 及归一化相干度。这是一阶时间相干信息；在致密、均匀频率选择上作相应双频逆变换可得到互相干函数。

选定参考 CSD 以有界成本给出相干算符切片。Gram 本征值给出全局相干模权重，但若要重构角向本征函数，还需要逐窗口角场；第一版为控制存储而不重复保存这些场。
