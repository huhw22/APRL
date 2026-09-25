# 运行时粒子/场能量账本

可选运行时账本是双粒子面诊断的全局补充。它在相同的**相对论变换系时间**切片上审计自洽模拟，不尝试从一个下游实验室平面反推域内近场储能。

## 守恒方程

无入射波边界源、无粒子退休时，账本计算

```text
K_accounted = K_active + K_removed

R = (K_accounted - K_initial)
  + (U_field - U_field_initial)
  + E_field_out
  - W_prescribed.
```

`K` 是计算系宏粒子动能，使用稳定恒等式

```text
gamma - 1 = |u|^2 / (gamma + 1)
```

累加。`U_field` 是六个 CPML 入口内物理区的单元中心 E/B 能量，不含 CPML 单元。`E_field_out` 是穿过六个内表面的向外 Poynting 通量时间积分；即使完整账本记录采样较稀，也在每个 Maxwell 步用梯形规则累积。PEC 方向相应通量设为零。

物理粒子进入 CPML 或离开物理区时，其插值过面动能从 `K_active` 转入累计 `K_removed`。这是物质能通量，不同于电磁 Poynting 通量，避免粒子逃逸表现为无解释能量损失。

实验室磁元件不存于 Maxwell 网格。静态理想磁铁在实验室不做功，但 Lorentz 变换后在计算系具有电场分量，所以推进器在每个 Boris 子步沿路径用梯形线积分累积 `W_prescribed`。省略该项会让波荡器计算在计算系看似不守恒。

`R` 保存为 `balance_residual_J`。算法不会强制其为零：电流沉积/场采样、Boris 推进、共点能量采样和半步 Yee B 并非严格代数能量守恒离散。残差必须随网格、宏粒子数、初始场 padding 与器件采样细化而收敛。

## 配置与成本

```yaml
energy_ledger:
  enabled: true
  directory: output/run-name
  filename: energy-ledger.h5
  sample_interval_steps: 10
  buffer_records: 64
  compression: 0
  warning_relative_tolerance: 0.01
```

块缺失或 `enabled: false` 时，不存在边界功率循环、粒子做功诊断、MPI 归约、缓冲或文件。开启后：

- 每个场步在本地累积六面功率和规定器件做功；
- 完整体场能量、粒子矩和 MPI 归约只在每 `sample_interval_steps` 执行；
- 只有 rank 0 打开文件；
- `buffer_records` 控制小型追加缓冲；
- 交互中断提交可读文件，`complete=0`；正常配置停止为 `complete=1`。

`warning_relative_tolerance` 默认 0.01。结束时若 `abs(relative_balance_initial)` 超限只发警告，不丢弃昂贵运行，并建议扫描网格/时间步、宏粒子和初始场 padding。`closure_valid=1` 只表示方程所需源类别均被记录，不表示数值残差通过阈值。TF/SF 入射波由于边界做功未入账会使其为 0；实验性退休开始后也为 0，因为退休电流不严格连续。

## HDF5 契约与报告工具

根组为 `/energy_ledger`。读取器只使用 `records[0:committed_records]`；`complete` 标记正常物理完成。每条记录包含：

- 步数和计算系时间；
- 活动/已移除宏粒子数及活动代表电子数；
- 活动与累计移除粒子的动能和总能量；
- 物理区场能；
- 累计向外场能总量及 `x-,x+,y-,y+,z-,z+` 六面分量；
- 累计规定源做功和两种残差归一化；
- 场能相对动能加场能、以及相对含静能粒子总能量的比例；
- 加权实验室与计算系 gamma 矩；
- 投影实验室相对 rms 能散、线性纵向 chirp 和扣除最佳加权线性 gamma-z 趋势后的 rms 展宽。

```bash
./build/energy_ledger_report output/run-name/energy-ledger.h5
./build/energy_ledger_report output/run-name/energy-ledger.h5 170
```

第二种形式额外打印从 0 开始的某个 committed 样本及其相对初始化的变化，可用于检查邻近固定实验室入口面的账本样本，不需要 Python。

## 能散增长的含义

能散是方差，不是守恒方程中的独立能量库。粒子可通过自洽场交换能量，使 rms 展宽增长而均值和总能量几乎不变；相干辐射也可降低均值，同时相关能量调制提高展宽。

投影 `sigma_gamma_lab` 至少混合：

1. 相干纵向相关/chirp，由加权最小二乘斜率 `linear_chirp_gamma_per_m` 表示；
2. 扣除线性趋势后的 `uncorrelated_sigma_gamma_lab`。

该分解适合调试，但位于等计算系时间超曲面。相对同时性意味着它不是完整等实验室时间束流诊断。加速器质量的入口/出口相空间和切片能散仍需固定实验室粒子面。尤其不能把报告的实验室粒子能量变化与计算系场能直接相加，它们位于不同超曲面。

## Gamma-1000 控制结果

当前 52×52×400、128 宏粒子、三周期控制先使用一个电子总电荷，再使用 `10^6` 电子，均启用满足 Gauss 定律的初始自场；无粒子离开物理区。下表为最终停止时计算系账本项：

| 电荷/K | delta K (J) | delta U field (J) | field out (J) | prescribed work (J) | residual (J) | residual/exchange |
|---|---:|---:|---:|---:|---:|---:|
| 1e / 0 | 1.32e-32 | 4.05e-25 | -2.08e-25 | 0 | 1.97e-25 | 0.321 |
| 1e / 0.5 | 5.51e-24 | 4.89e-24 | -7.31e-26 | 9.96e-24 | 3.70e-25 | 0.0181 |
| 1e6 e / 0 | 1.32e-14 | 3.92e-13 | -2.08e-13 | 0 | 1.97e-13 | 0.321 |
| 1e6 e / 0.5 | 1.14e-14 | 4.87e-12 | -7.32e-14 | 4.44e-12 | 3.62e-13 | 0.0385 |

近零 `K=0` 的交换归一化百分比条件很差，尽管其绝对残差更小。受驱动大电荷算例闭合到交换项和的约 3.9%，尚未达到正式容差；它支持账本和符号约定，但仍要求网格/宏粒子/padding 收敛。

`10^6` 电子时，初始场能占含粒子静能总量的 `5.53e-5`；最终 `K=0` 为 `6.01e-5`，`K=0.5` 为 `1.15e-4`。相反，场能除以**计算系动能加场能**接近 1，因为所选 boost 几乎消除了平均束流运动；两个分母回答不同问题，并不矛盾。

受驱动大电荷组末态投影 `sigma_gamma=4.376e-2`，线性去趋势后为 `1.097e-2`，约 93.7% 投影方差由线性纵向相关解释；匹配 `K=0` 仅产生 `sigma_gamma=1.43e-4`。所以该控制中大部分能散是相干纵向相空间结构，不是新能量库，也不能自动解释为不可逆切片能散增长。
