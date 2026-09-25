# 粒子/场能量闭合诊断

`postprocess/energy_closure` 比较宏粒子在两个固定实验室粒子面之间损失的能量与 `field_plane_analysis` 积分的前向辐射能量。它是只读独立 C++ 工具；不运行时不会增加模拟通信、内存或输出。

报告严格采用实验室观察者，不把粒子面与相对论变换系运行账本混合。格式版本 3 报告两个平面的实验室平均 gamma 和 rms 能散、零磁场控制交换、信号减控制后与器件相关的粒子变化、未由所收集前向频段解释的部分，以及可选双面纵向 Poynting 账本。解释和 50 A 快速尺度工具见[实验室系能量诊断](LAB_FRAME_ENERGY_DIAGNOSTICS.md)。

## 推荐实验

入口粒子面放在第一个磁作用边界前，出口面放在最后边缘场后；场面必须保存完整脉冲。物理算例和匹配零辐射算例使用相同网格、束团、自场、边界、探测节拍和随机种子。单纯波荡器测试可设 `strength_parameter: 0` 作为基线。

用各自弹道参考清理每个场面，再通过 `zero_radiation_baseline` 把清洁信号与清洁基线送入 `field_plane_analysis`。最后向本工具提供两对粒子面和该分析文件，完整示例为 `postprocess/energy_closure/example.yaml`。所有映射拒绝未知键，已有报告默认拒绝覆盖。

更强的能量检查还需在入口和出口粒子面同位置重构原始总场，并提供四个可选 `*_field_reconstruction` 路径。报告计算

```text
Delta W_z = [(W_exit - W_entry)_signal
             - (W_exit - W_entry)_K=0],
W_plane = integral dt integral_A (E cross B)_z / mu0 dA.
```

`W_plane` 带符号，向负 z 的能量为负；使用截断前向功率会破坏控制体平衡。原始总场是首要守恒量。清洁值只作模型诊断，因为在二次 Poynting 运算前扣除估算粒子本底并不是精确能量恒等式。

少粒子验证中，`field_analysis` 也可指向完整 `/far_field` 轨迹辐射输出。报告会记录标量能量来自场面还是独立轨迹路线；后者是交叉检查，不替代正式场探测器。

```bash
cmake -S postprocess/energy_closure -B build-energy-closure
cmake --build build-energy-closure -j
./build-energy-closure/energy_closure closure.yaml
```

## 稳定的高 gamma 计算

工具按 `particle_id` 连接下游过面；重复下游记录被拒绝。`require_all_particles: true` 时缺少任一入口或出口也是硬错误。宏粒子 `mass_kg` 已含代表粒子缩放，不能再次乘诊断 `weight`。

逐粒子避免两个已舍入 gamma 直接相减：

```text
gamma_in - gamma_out
  = (|u_in|^2 - |u_out|^2) / (gamma_in + gamma_out),
gamma = sqrt(1 + |u|^2).
```

随后以扩展精度计算 `mass_kg*c^2*delta_gamma` 并补偿求和。存在匹配基线时，信号减基线也在逐粒子层面完成。这比 gamma 约 1000 时两个束流总能量相减安全得多。

## 闭合含义

静态磁元件且无注入激光时，规定磁场在实验室系不做功。扣除匹配数值基线后，粒子动能损失应与发射电磁能量处于可比尺度。单场面只测其孔径和分辨频段内的前向部分，通常不应大于粒子损失。

匹配 `K=0` 粒子变化作为空间电荷加数值控制报告，不能直接命名为近场储能；精确近场/辐射拆分需要闭合实验室曲面或等实验室时刻三维场。因此 `device_associated_particle_loss-collected_forward_radiation` 只报告为未解析部分，而不强行归因。

比较不是普适精确恒等式。残差还含横向/后向边界辐射、频段外能量、有限时间/孔径损失、初末束缚自场差、探测器前 CPML 吸收及离散误差。正式结论前必须扫描粒子数、网格、孔径、时间窗、探测距离、子步、CPML 余量和场重构补零。

密集束团中，匹配 `K=0` 不保证信号和基线在出口具有相同束缚/自场能量。如果基线动能变化已可比或大于预期辐射损失，应直接报告，不能把差值自动解释为辐射。

运行时 `energy_ledger` 在等计算系时间片上审计

```text
particle kinetic change
+ domain electromagnetic-energy change
+ flux through every physical/CPML boundary
- prescribed-source work = 0.
```

见[运行时能量账本](ENERGY_LEDGER.md)。实验室粒子面与场面是等位置比较；相对同时性禁止把这些实验室粒子项直接加到计算系储能上。本工具不会从一个下游平面虚构缺失的全局项。

## 探测时间窗收敛示例

两个场面都必须包含完整脉冲。固定起止位置时增加纵向网格数不会延长实验室时间窗，只会改变模拟盒；应把束团起点前移、停止位置后移，并检查功率曲线首尾样本。

一个 `gamma=1000`、128 宏粒子、`10^6` 电子、三周期 `K=0.5` 测试配合 `K=0` 得到：

| 实验室量 | 结果 |
|---|---:|
| 基线修正粒子损失 | `5.919605 nJ` |
| 匹配原始带符号出口减入口通量 | `5.874187 nJ` |
| 原始纵向闭合 | `99.2328%` |
| 粒子减原始通量残差 | `0.045418 nJ` |
| 振幅扣基线后的 0–300 eV 前向能量 | `5.430281 nJ` |
| 前向频段占粒子损失 | `91.7338%` |
| 原始带符号通量减前向频段 | `0.443906 nJ` |

最短时间窗只闭合 54.57%；前移束团起点后为 62.45%，停止延至 0.35 m 后为 83.38%。仅把 `Nz` 加倍而不延长运行仍为 62.45%，说明早期缺口主要来自时间窗截断。剩余 0.77% 原始残差可包含横向通量、两面间储能、尾部截断和离散误差；相对所选前向频段的 8.27% 还含非传播/束缚场、被排除频率或角度及场分解交叉项。
