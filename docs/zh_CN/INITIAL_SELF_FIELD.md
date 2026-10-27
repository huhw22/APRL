# 初始粒子自场

## 目的

第一次物理推进前，Maxwell 状态必须满足已加载粒子电荷对应的离散 Gauss 约束。非零束团配零初始 E 会产生并非粒子辐射的数值启动脉冲。默认初始化因此一次性构造刚性束流 E/B，随后释放全部势和求解工作数组；推进状态仍为直接 SI E/B，A/phi 从不进入时间推进。

## 默认离散构造

1. 每个计算系粒子用与电荷守恒轨迹沉积器相同的一阶 CIC 形状，把电荷沉积到 Yee 顶点。
2. MPI z 接口重复顶点的贡献相加；每个全局 z 平面只有一个归约所有者，点积不会重复计算共享平面。
3. 计算模拟系中以代表质量为权的平均轴向速度 `beta_bar`，再由无显式矩阵的分布式共轭梯度求解内部顶点标量势：

   ```text
   -[Dx^2 + Dy^2 + (1-beta_bar^2) Dz^2] phi = rho / epsilon0
   ```

   这是 Vay 构造中的轴向刚性束相对论 Poisson 方程；当 boost 跟随束团平均速度时连续退化为普通 Poisson。
4. 在 E/B 各自 Yee 位置构造

   ```text
   E = -grad(phi) + beta_bar (beta_bar dot grad(phi))
   B = beta_bar cross E / c
   ```

   E 存在时间零；B 用刚性平移 `z -> z+beta_bar*c*dt/2` 计算到 leapfrog 时刻 `-dt/2`。第二次 CIC 沉积和真正跨 rank 散度检查验证 `div(E)=rho/epsilon0`，通过后才能继续。
5. CPML 运行以 CPML 内表面作为静态零势边界；CPML 场单元与辅助变量从零开始，避免静态 Coulomb 尾被装入 CPML 而记忆变量仍为零造成启动能流。PEC 回归则以网格外表面为零势边界。
6. 求解器要求完整 CIC 电荷严格位于静态边界内部，并报告束团每个 RMS 所含网格数和质心到边界的 RMS 余量；分辨不足或余量小于 4 RMS 会警告。

标量势只是满足约束的临时数值变量。CG 期间同时存在四个顶点薄层（`phi`、残差、方向、算符像）；电荷格只在沉积和后检查阶段分配。启动报告迭代数、Gauss 残差、电荷、平均速度、电/磁能、RMS 分辨率、padding 和每 rank 临时内存。

## Cowan 准稳态参考场

可选的 `model: cowan-relativistic-poisson` 使用 `mesh.field_solver: cowan-z` 对应的空间差分算符，配合其他 Maxwell 求解器会明确报配置错误。该选项不改变默认模型、粒子位置和动量、held 释放条件、网格尺寸或传播设置。

以 `D*` 表示 Cowan 修正梯度，`D` 表示普通 Yee 梯度，`v = beta_bar*c` 为计算系中按代表质量加权的平均轴向速度，`a = v*v*epsilon0*mu0`。分布式 CG 求解并构造

```text
-div[(D* - a z_hat Dz) phi] = rho / epsilon0
E = -D* phi + a z_hat Dz phi
B = v epsilon0 mu0 z_hat cross E
```

Gauss 后检查仍使用相同 CIC 电荷和普通 Yee 电场散度；B 仍通过刚性平移插值到 `-dt/2`。修正梯度采用生产 Cowan 模板，见[Cowan-z 与 CPML](MAXWELL_COWAN_CPML.md)。

该模型给出与 Cowan 空间算符相容的准稳态参考场，并非精确的离散运动平衡。有限时间步平移、速度展宽、有限边界以及后续正常空间电荷膨胀仍会导致粒子和场交换能量。应在匹配的小工况中对照原 held 初始化的能量交换、束团膨胀和能散；差异应足够小以满足目标观测量，而不要求逐位吻合。初始化器不新增经验性的 held 吻合阈值，CG 收敛与 Gauss 后检查仍是必要数值约束。Cowan 的磁场约束使用修正散度，不能仅以普通 Yee `div(B)` 是否为零作为对应验收条件。

CG 前要求 `min(S_x) > 64*epsilon`、`min(S_y) > 64*epsilon` 和 `min(S_z)-a > 64*epsilon`；`epsilon` 为浮点机器精度，`S_i` 为方向 `i` 的横向平滑符号。各符号是两个横向余弦变量的双线性式，因此在 `[-1,1]^2` 四个角点计算其最小值。这是与网格格数无关的保守正定充分条件，并非必要条件：某些完整椭圆算子谱实际为正的有限网格也可能被拒绝。报错会给出方向及下界，不代表物理束流无效；现有网格形状和 boost 选择还须满足该 CG 模型的适用条件。

## 输入

```yaml
initial_self_field:
  enabled: true
  model: relativistic-poisson
  relative_tolerance: 1.0e-10
  maximum_mean_transverse_beta: 1.0e-8
  maximum_iterations: 10000
```

该块可省略，默认即为以上值。旧对照模型使用 `model: electrostatic-poisson`，Cowan 参考场使用 `model: cowan-relativistic-poisson`。启用求解是硬预检查：不收敛或后检查失败会给出修正建议并终止。只有受控零场/回归比较才应关闭，并会打印醒目警告。

`maximum_mean_transverse_beta` 控制两种刚性束模型在 boost 后计算系内按代表质量加权的平均横向 `|β|` 上限，默认仍为 `1e-8`。放宽它不会旋转刚性束初场，也不会改变电荷沉积；比较不同工况时应记录并说明所用上限。

## 物理范围

该步骤给出具有共同轴向速度的刚性分布、满足 Gauss 约束的场；它不是任意速度展宽的精确推迟 Liénard–Wiechert 解，也不能把任意采样相空间变为 Vlasov 平衡。有限、无约束电子束在自由空间没有静态平衡；必须通过匹配 `K=0` 漂移和宏粒子/网格/padding 收敛，把真实空间电荷膨胀与初始化误差分开。

静态面零势是有限域近似，不是无限空间边界。定量运行必须扫描束团到 CPML 入口的横向和纵向距离。小 Gauss 残差只证明离散约束，不证明盒子或采样分布物理充分。

使用 `mesh.cells: [..., auto]` 时，设置 `mesh.longitudinal.self_field_padding_rms: 0.0` 表示只在观测有界物理区上求解。实际采样到的全部 CIC 电荷仍被保留，完整的迭代收敛和后检查也不会取消；这不是把初场改成零。所得初态在框内满足离散 Gauss 约束，但框外 Coulomb 尾被有限边界截断。manifest 会同时记录初场启用状态与请求的 RMS padding，避免把框内截断和扩展边界工况混淆。

## 入口漂移验证

有限长度 `K=0` 测试使用 `gamma=boost_gamma=1174`、高斯峰值电流 50 A、`sigma_x=sigma_y=60 um`、4096 宏粒子代表 `10^6` 电子。重构实验室束团参考保持在 `-50 mm`；首个磁铁物理入口 `34 mm`，紧凑边缘作用区始于 `10.0238 mm`。头部锚定同步把共同时间实验室束头置于 `-49.9987 mm`，到作用边界余量 `60.0225 mm`，而一个单元实验室余量为 `11.74 mm`。变换后束团只占计算系 z 向 `3.2666 mm`，适合配置盒子。

由于相对同时性，相应实验室事件跨越 `3.8349` 光米；这只作为诊断，不是真实入口漂移，也不再强迫束团参考放到 `-2 m`。对 HDF5-v3/v4 Elegant 输入，独立物理要求是每条记录从 `input_plane_z` 向前推进到重构快照；存在磁铁时，快照还必须完整位于第一个磁作用区前。纯探测器和无元件漂移没有人为磁入口约束。

更新后的 `-50 mm` CPML 运行初始化场能 `0.202850 pJ`，正常到达 `200 mm` 参考停止；末态场能 `0.203168 pJ`，累计 CPML 接口通量 `-3.92e-6 pJ`。亚百分比变化是收敛目标，不是物理入口漂移要求。

把静态边界移到 CPML 入口之前，同一测试场能变化约 17.5%，并有 `-0.035 pJ` 从吸收层回流，说明大启动效应来自 CPML/静态尾初始化不匹配，而非普通自由传播。剩余约 0.1% 与粗网格账本残差相当，仍只是收敛目标。

## 参考文献

- J.-L. Vay, “Simulation of beams or plasmas crossing at relativistic velocity”, *Physics of Plasmas* 15, 056701 (2008), <https://doi.org/10.1063/1.2837054>。
- [WarpX relativistic electrostatic solver documentation](https://warpx.readthedocs.io/en/26.06/theory/models_algorithms/electrostatic_pic.html)。
