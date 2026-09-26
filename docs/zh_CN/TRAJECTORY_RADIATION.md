# 粒子轨迹到远场辐射工具

## 范围

`postprocess/trajectory_radiation` 是独立 C++/MPI 程序，不链接模拟内核，也不读取 Maxwell 网格。输入是一组或多组实验室系轨迹 HDF5；输出是单个 HDF5，包含复偏振远场谱、频谱/角能量密度、Stokes 量、角积分能谱，以及可选 shot 系综或时间窗的一阶相干矩阵。

这种分离是有意的：昂贵模拟只保存一次可复用粒子历史，之后可在小型分析服务器上改变频率范围、角网格、观察距离和系综统计，而无需重跑粒子动力学。

工具只计算辐射场，不含 Coulomb 或速度近场；匀速粒子输出严格为零。

## 辐射模型

对观察方向 `n` 定义

```text
F(beta) = n cross (n cross beta) / (1 - n dot beta).
```

采样路径在实验室时空中按分段直线处理。每个内部轨迹记录处，无量纲弦速度由 `beta_before` 变为 `beta_after`，带电荷权重的谱振幅相干累积为

```text
A(omega,n) = sum_particles q sum_internal_knots
             [F(beta_after) - F(beta_before)]
             exp{i omega [t - n dot r/c]}.
```

这是分段直线路径的远区 Liénard–Wiechert 加速度场端点表示。人为首末启动/停止端点被抑制，即假定记录范围外路径继续惯性运动，所以每个粒子至少需要三条记录。CPML 和退休终止事件可作为合法末点，但非物理载流子不参与。抑制末端点也使退休记录等价于真实轨迹的惯性延续。

程序用 `long double` 计算推迟时间、相位约化、本地相干累加和 MPI 归约，输出为 float64。轨迹采样仍必须分辨速度变化和规定磁场；正式结果必须同时做频率网格与轨迹节拍收敛。

距离 R 处保存的复电场谱为

```text
E_tilde(omega,n;R) = exp(i omega R/c) A / (4 pi epsilon_0 c R),
E_tilde = integral E(t) exp(+i omega t) dt.
```

正频率谱角能量为

```text
d^2 W / (d omega d Omega)
  = |A|^2 / (16 pi^3 epsilon_0 c).
```

实现依据 A. G. R. Thomas, *Phys. Rev. ST Accel. Beams* **13**, 020702 (2010)，<https://doi.org/10.1103/PhysRevSTAB.13.020702>。第一版内核仅实现分段直线端点形式，不提供高阶路径插值；现有输出契约与高阶插值兼容。

## 构建与运行

```bash
cmake -S postprocess/trajectory_radiation \
      -B build-radiation \
      -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_CXX_COMPILER=mpic++
cmake --build build-radiation -j

./build-radiation/trajectory_radiation radiation.yaml
mpirun -n 16 ./build-radiation/trajectory_radiation radiation.yaml
```

每个输入文件只由一个分析 rank 打开。记录按 `particle_id` 重分配，因此跨模拟 rank 迁移的粒子会先重组完整历史。只有分析 rank 0 写输出。频率和 theta-y 分块限制临时场内存。

可选时间平均一次只计算一个窗口，并只累积请求的相关量，不保存每个窗口的完整场。启动日志报告累加器大小，超过 `maximum_accumulator_mib` 会在轨迹分析前拒绝。一个分析 rank 所拥有的轨迹在当前 shot 内保持内存；分块限制辐射数组而非轨迹集。单 rank 不受 MPI 消息长度限制；多 rank 时单个点对点重分配桶必须适合 MPI 有符号 `int` 字节数（约 2 GiB），达到前应增加分析 rank。

## YAML 输入

完整例卡：`postprocess/trajectory_radiation/example.yaml`。

```yaml
input:
  require_complete: true
  read_chunk_records: 65536
  shots:
    - name: shot-000
      files:
        - particles-rank-00000.h5
        - particles-rank-00001.h5

observation:
  axis: [0.0, 0.0, 1.0]
  horizontal: [1.0, 0.0, 0.0]
  distance_m: 100.0
  theta_x_rad: {min: -0.001, max: 0.001, count: 101}
  theta_y_rad: {min: -0.001, max: 0.001, count: 101}

spectrum:
  photon_energy_eV:
    min: 100.0
    max: 1000.0
    count: 181
    spacing: linear

calculation:
  frequency_block: 16
  theta_y_block: 4
  minimum_records_per_particle: 3

output:
  file: radiation-output.h5
  compression: 0
  overwrite: false
```

路径相对辐射卡解析。一个 shot 是一个统计独立实现，可含任意数量每-rank 轨迹。接受轨迹格式 1、2。`require_complete: false` 可读取有意中断但已提交的前缀。参与场平均或双频 CSD 的 shots 必须共享有意义的实验室时间原点；有意随机到达时间是物理量，而意外逐文件时间重置会抹掉相位敏感系综量。

中心 `axis` 和投影后的 `horizontal` 定义正交基，网格方向为

```text
n = normalize(axis + tan(theta_x) horizontal + tan(theta_y) vertical).
```

光子能量支持 `linear` 或 `log`，角轴为线性。

## 准平稳时间平均

稳定或缓慢演化的脉冲区段可用短时谱作为平稳/准平稳过程样本：

```yaml
time_average:
  enabled: true
  interval_s: [5.0e-15, 25.0e-15]
  window_duration_s: 4.0e-15
  window_step_s: 2.0e-15
  photon_energy_eV: [450.0, 475.0, 500.0, 525.0, 550.0]
  reference_angles_rad:
    - [0.0, 0.0]
  maximum_accumulator_mib: 1024
```

每个 shot/window 对是一个样本。固定 Hann 窗作用于约化观测者时间

```text
u = t_lab - n dot r_lab/c.
```

省略共同传播延迟 `R/c`，所以 `interval_s` 接近辐射发射时刻，不含任意探测距离。每个观察方向独立判断窗口归属。日志报告每个 shot 中心轴内部折点的实际 u 范围，是选择区间的第一实际边界；角边缘范围可略有不同。

区间应位于稳定平台。包含启动、饱和瞬态或脉冲衰减会把确定性包络变化与相干损失一起测量。重叠窗口让采样更平滑但互相相关，不等于同数量独立实现。

`/time_average` 包含：

| 数据集 | 形状 | 含义 |
|---|---|---|
| `spatial_cross_spectral_density` | `[ref,f,y,x,2,2]` complex | 时间窗 CSD `mean(conj(E_ref,a) E_target,b)` |
| `spectral_degree_of_coherence_squared` | `[ref,f,y,x]` | 基不变电磁 `mu_EM^2` |
| `mean_window_spectral_energy_density` | `[f,y,x]` | 每 Hann 窗算术平均能谱 |
| `reference_mean_window_spectral_energy_density` | `[ref,f]` | 参考角对应平均 |
| `reference_window_spectral_energy_density` | `[shot,window,ref,f]` | 紧凑逐窗平稳性诊断 |
| `two_frequency_cross_spectral_density` | `[ref,f1,f2,2,2]` complex | 双频时间窗相关 |
| `solid_angle_quadrature_weight` | `[y,x]` | 切线角网格 `dOmega` 权重 |

保留完整偏振矩阵，后续可重新归一或投影偏振而无需重跑轨迹。逐窗参考强度能显示启动、衰减和平台漂移，又不保存每窗全角场。双频 CSD 含一阶时间相干信息；相应逆傅里叶变换得到互相干函数，建议使用致密均匀频率网格。

设 `reference_angles_rad: all` 会把每个角网格点都作为参考，生成相干模分解所需的完整角 CSD。存储随角点数平方增长，通常应缩小角网格或只选代表参考。累加器内存上限防止意外多 GiB 分配。

偏振相干模分解时，把角度和偏振展平为一个索引，对 `sqrt(dOmega_i) W_ij sqrt(dOmega_j)` 做对称化后对角化；保存的立体角权重已含切线坐标 Jacobian 和梯形边缘因子。

时间窗平均描述选定脉冲部分的未分辨变化，与 shot-to-shot 统计及完整脉冲确定性谱分开存储。

## HDF5 输出

`/far_field` 的 `format_version=1`，含标量 `complete`。绝对谱能量归一和当前 40/80/160 步单电子峰值强度收敛见[数值验证](VALIDATION.md)；只看峰位置一致不能验证强度。

主坐标为 `photon_energy_eV`、`omega_rad_per_s`、`frequency_Hz`、`wavelength_m`、`theta_x_rad`、`theta_y_rad`；`observation_basis` 保存中心、水平和垂直基向量。

全部逐-shot 场、系综摘要、谱和相干产品写入并刷新后 `complete` 才置 1。中断文件可含有用诊断块，但不能视为完整辐射结果或重启文件。

| 数据集 | 形状 | 含义 |
|---|---|---|
| `electric_field_spectral` | `[shot,f,y,x,2]` complex | 逐-shot 水平/垂直 `E_tilde`，V s/m |
| `spectral_energy_density` | `[shot,f,y,x]` | `d²W/(dω dΩ)`，J s/sr |
| `stokes_spectral_energy_density` | `[shot,f,y,x,4]` | I,Q,U,V，同能量单位 |
| `mean_electric_field_spectral` | `[f,y,x,2]` complex | 系综平均场 |
| `mean_spectral_energy_density` | `[f,y,x]` | 系综平均强度 |
| `coherent_spectral_energy_density` | `[f,y,x]` | 系综平均场的强度 |
| `fluctuation_spectral_energy_density` | `[f,y,x]` | 平均强度减相干强度 |
| `energy_spectrum` | `[shot,f]` | 角积分 `dW/dω`，J s |
| `mean_energy_spectrum` | `[f]` | 系综平均 `dW/dω` |
| `mean_energy_per_log_frequency` | `[f]` | `dW/dln(ω)`，J |
| `band_energy_J` | `[shot]` | 配置频率和角网格内积分，J |

属性 `mean_band_energy_J`、`band_min_photon_energy_eV`、`band_max_photon_energy_eV` 提供匹配场功率比较所需标量频段积分。

复数类型含 float64 `real`、`imag`。偏振索引 0 是投影水平基，1 补成右手横向基；Stokes 约定写入属性。角积分使用两个切线角坐标的精确 Jacobian 和梯形权重；任一角轴只有一点时不可积分。

## 时间与空间相干性

一阶相干性是系综量：

```text
W_ab(p,q) = mean_shots[conj(E_a(p)) E_b(q)].
```

单个确定性 shot 的外积必然秩一，不能独自证明部分相干，所以始终保留逐-shot 复场。多个统计独立 shots 可生成选定 CSD：

```yaml
coherence:
  spatial_photon_energy_eV: [500.0]
  spatial_reference_angles_rad:
    - [0.0, 0.0]
  temporal_photon_energy_eV: [450.0, 475.0, 500.0, 525.0, 550.0]
  temporal_reference_angles_rad:
    - [0.0, 0.0]
```

`/coherence/spatial_cross_spectral_density` 保存选定能量处参考角和全部角点之间的 2×2 偏振 CSD；`/coherence/temporal_cross_spectral_density` 保存各参考角上选定频率对之间的 2×2 CSD。请求坐标映射到最近实际网格点，并同时保存所选坐标。

CSD、互强度及其与谱密度的关系可参考 “Coherence properties of the high-energy fourth-generation X-ray synchrotron sources,” *J. Synchrotron Radiat.* **26** (2019)，<https://doi.org/10.1107/S1600577519013079>。

## 宏粒子与收敛注意事项

- `charge_C` 已包含宏粒子代表电荷并相干求和，不能再次乘诊断 `weight`；
- 只有输入实现和宏粒子模型确实代表物理 shot noise 时，系综涨落才有物理意义；对粗糙平滑分布重采样主要测到数值宏粒子噪声；
- 首末端点抑制假定惯性延续，应在辐射元件前后记录足够漂移，使保留内部速度变化覆盖完整加速区；
- 只有两条记录的粒子远场为零是有意行为：单弦没有已分辨内部加速度；
- 必须依次细化轨迹节拍、角间距、角孔径和光子能量间距；FEL 相位尤其敏感于轨迹节拍和推迟时间精度；
- 这是远区工具；近场、有限孔径 Fresnel 传播、光学元件和探测器响应需另行处理。
