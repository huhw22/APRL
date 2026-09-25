# 粒子退休与匹配功率比较

## 目的与限制

粒子退休是一条诊断路线，用于在不修改 E/B 传播内核的情况下获得更干净的下游场面。真实粒子在全部磁作用区之后的固定实验室入口终止，一个小型弹道载流子继续沉积剩余直线电流，以 C2 五次函数平滑衰减，并在探测器前消失。

这并不是对净电荷的严格电荷守恒移除。平滑性只控制数值冲击的频谱，不能证明冲击可忽略。只有当零辐射功率相对目标辐射足够小，并随参数收敛而下降时，才能接受该路线。

第一条基线是第二次模拟：网格、CPML、束团、随机种子、退休入口/长度、探测器、停止和采样时刻全部相同，只把辐射磁强度设为零。比较工具在 E/B 振幅层面扣除基线，因而同时测到直束近场以及退休、初始化和边界伪影。它是筛选基线；若磁铁改变了退休入口的粒子分布，就不是严格反事实。未来若从入口状态回放，才能更精确隔离退休伪影。

## 构建与运行

```bash
cmake -S postprocess/field_power_compare -B build-field-power-compare
cmake --build build-field-power-compare -j
./build-field-power-compare/field_power_compare \
  postprocess/field_power_compare/example.yaml
```

输入卡指定信号场、匹配零辐射场、可选 `field_reconstruction` 清洁结果、可选轨迹远场、光子能量频段、小型横向读批次和输出文件。所有映射拒绝未知键；除非显式 `output.overwrite: true`，否则拒绝已有输出。

两个场文件必须具有完全相同的平面几何和实验室样本时刻。实际时刻不要求严格均匀：宽带功率保留原时刻，有限记录 DFT 前把各场分量线性重采样到共同均匀网格。内存只包含小横向批次的时间序列和单点四个滤波分量，不复制整个探测面；模拟过程中不增加 MPI 或 I/O。

若提供 `analytic_electron_reconstruction`，其几何、样本数和实验室时间范围必须相同；均匀时间轴可以不同于原始探测器的不均匀内部样本。工具把两条路线都重采样到同一频段分析轴，并计算退休估计与解析直线电子扣除之间的场级残差。

## 定义

信号 `s` 和零辐射基线 `0` 必须先在场振幅层面相减：

```text
E_delta = E_s - E_0
B_delta = B_s - B_0
P_delta(t) = integral max(((E_delta cross B_delta)_z)/mu0, 0) dx dy.
```

工具也保存带符号 Poynting 功率。绝不使用 `P_signal-P_baseline`，因为那会丢失干涉项，也不是差场的功率。

每个横向点的 DFT 先清零频段外频点再逆变换 E/B。输出包含信号、基线和差场的宽带/限带带符号与前向功率曲线、梯形积分能量和峰值前向功率。

限带功率有两种定义。旧的瞬时曲线直接使用实滤波场。与轨迹严格比较时，先形成 Hilbert 正交分量，再计算

```text
P_cycle(t) = integral 0.5 Re((E + i H[E]) cross
                             conjugate(B + i H[B]))_z / mu0 dx dy.
```

该解析信号周期平均功率与 `trajectory_band_power_W` 约定相同。实场瞬时峰值可以接近周期平均的两倍而两者都正确。旧 `field_to_trajectory_*` 属性仅为兼容；验证应使用 `cycle_averaged_field_to_trajectory_*`。

提供 `trajectory_far_field` 时，它必须只有一个 shot、线性频率间隔、相同光子能量边界，并且两个角轴都至少两点。工具从 `electric_field_spectral` 重建周期平均功率，在配置立体角内积分，把周期时间原点移到峰值，并把时间积分归一化到 `band_energy_J`，然后报告场/轨迹峰值与能量比例。

解析重构保持为独立路线。工具报告其宽带、限带和周期平均功率，以及相对退休差场和轨迹的比例，并构造

```text
E_residual = E_retirement_difference - E_analytic_reconstruction
B_residual = B_retirement_difference - B_analytic_reconstruction.
```

这能区分清洁场本身的分歧，避免用两个分别平方的功率差误导判断。轨迹角网格必须代表与有限场面相同的辐射接受度；程序不会猜测这一关系，因为近场平面通常不能映射为单源点和简单角矩形。

## HDF5 输出

`/power_comparison` 包含：

- `time_s`、`band_time_s`；
- `{signal,baseline,difference}_{signed,forward}_power_W`；
- `{signal,baseline,difference}_band_{signed,forward}_power_W`；
- `{signal,baseline,difference}_band_cycle_averaged_{signed,forward}_power_W`；
- 对应前向能量和峰值属性；
- `baseline_to_difference_band_energy_ratio`、`baseline_to_difference_band_peak_power_ratio`；
- 可选轨迹相对时间、功率、能量/峰值、旧瞬时比例及严格周期平均比例；
- 有解析重构时的原始、限带和周期平均功率，限带 `method_residual` 曲线/指标，退休/解析比例和解析/轨迹比例；
- 探测面积、z、能量上下限、Nyquist 能量、定义与完成标志。

## 接受流程

不存在普适百分比。每个物理问题都应：

1. 确认探测器 Nyquist 高于目标上限，完整脉冲落在时间窗内；
2. 要求零辐射基线峰值和能量相对差场足够小；
3. 在小粒子测试中，用真正匹配的角接受度比较差场与轨迹的周期平均峰值和能量；固定总电荷提高宏粒子数，因为稀疏相干束团可能让积分能量变化到量级 1，而峰值仍看似合理；
4. 对同一辐射问题关闭退休并执行 `field_reconstruction`，检查退休/解析比例和场级残差，单独扫描解析横向平滑；
5. 至少比较两个退休长度和探测器间距；目标功率/能量应稳定，基线伪影应下降或保持次要；
6. 扫描纵向分辨率/节拍、横向分辨率、孔径、CPML 厚度及频率/角网格。

观察角 `theta` 对应横向波长约 `lambda/sin(theta)`；即使 z 分辨网格在轴附近一致，dx/dy 过大仍会漏掉宽角功率，不能把这种失败归因于退休。
