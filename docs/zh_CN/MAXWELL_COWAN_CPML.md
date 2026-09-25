# Cowan-z Maxwell 内核与 CPML 边界

## 可控色散推进

主求解器保留普通 Yee 格式的交错 E/B 存储，只修改 Faraday 更新：

```text
Dt B = -curl'(E)
Dt E =  c^2 curl(B) - J/epsilon0
```

对方向 `i` 的导数，`curl'` 使用 `D_i S_i`；`S_i` 是位于另外两个横向方向上的九点平滑模板。平滑量就地计算，不会保留三份全域场副本。

令 `lambda` 不大于任一网格间距，并定义

```text
r_i = lambda^2 / delta_i^2
D   = r_x r_y + r_y r_z + r_z r_x .
```

实现的 Cowan 系数为

```text
beta_i     = r_i/8 * (1 - r_x r_y r_z / D)
delta_ij   = r_i r_j * (1/16 - r_i r_j/(8 D))
alpha_i    = 1 - 2 beta_j - 2 beta_k - 4 delta_jk .
```

最后一个关系保证每个平滑模板系数和为 1，因此常量场既不衰减也不放大。使用这些系数时，真空色散关系可以分解为有界乘积。解析稳定条件为：每个方向都满足 `lambda <= delta_i`，并且 `c*dt <= lambda`。

对 z 优先求解器取 `lambda=dz`。所以 `dx>=dz`、`dy>=dz` 与 Cowan 系数公式及 `c*dt<=dz` 必须同时满足；只满足前两个网格不等式并不能构成一般稳定性证明。程序取 `c*dt=dz`，使已分辨真空波沿 z 精确无色散传播。这位于解析稳定边界，因此长期、粒子耦合和边界耦合回归仍是必需的。

立方网格下所有 `r_i=1`，代入得到

```text
beta  = (1/8)(1-1/3) = 1/12
delta = 1/16 - 1/24  = 1/48
alpha = 1 - 4 beta - 4 delta = 7/12 .
```

这三个数只适用于立方网格。各向异性正式网格会得到不同系数，并打印在启动日志中。

参考：B. M. Cowan 等，[Generalized algorithm for control of numerical dispersion in explicit time-domain electromagnetic simulations](https://doi.org/10.1103/PhysRevSTAB.16.041303)，2013。

## CFS-CPML 耦合

边界采用非分裂复频移卷积 PML。方向 `i` 上的拉伸导数为

```text
D_i -> D_i/kappa_i + psi_i
psi_i(n) = b_i psi_i(n-1) + c_i D_i(n) .
```

对 Cowan Faraday 更新，`D_i(n)` 指完整的修正差分 `D_i S_i(E)`；Ampere 更新使用普通交错 `D_i(B)`。因此 PML 不会在边界单元中悄悄退回 Yee 磁场模板。

电导率、`kappa` 和可选频移按多项式渐变。`target_reflection` 决定连续剖面的最大电导率，它只是设计参数，不能保证所有波长、入射角和层数下都得到同一离散反射率。

项目提供的 FEL 示例使用零频移 `alpha_fraction=0`。在当前近轴前向脉冲回归中，它比最初测试的非零频移产生更小反射。非零 alpha 仍可用于有针对性的低频或倏逝场研究，但不能假定总是更优。

每个场分量有两个卷积历史，共十二类；但每段历史只在相应导数坐标处于 PML 薄层时分配。关闭的方向和非 PML 内部点不占历史存储。最外表面由 PEC 终止。

当前限制：

- CPML 开启时拒绝 TF/SF 注入种子波；未来需要把闭合注入面放在 CPML 内部并补齐修正项；
- 每个 z 向 CPML 薄层必须完整位于端点 MPI rank 内，程序会在读取粒子前检查，并建议修改 rank 数、z 网格数或厚度；
- CPML 用于吸收场；粒子采用[考虑 CPML 的开放边界](PARTICLE_OPEN_BOUNDARY.md)，而不是新的物理材料模型。

参考：J. A. Roden 与 S. D. Gedney，[Convolution PML: An efficient FDTD implementation of the CFS-PML for arbitrary media](https://doi.org/10.1002/1098-2760(20001205)27:5%3C334::AID-MOP14%3E3.0.CO;2-A)，2000。
