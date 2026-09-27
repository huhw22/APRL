# 独立 SASE 初始束团生成器

## 适用范围

`sase_bunch_generator` 为 SASE 启动研究生成 APRL HDF5-v2 共同实验室时刻粒子
快照。它是独立预处理工具，不属于模拟内核；运行时内核不会生成、重新归一化或修改
已加载的 shot noise。

纵向密度为包含整数个共振波长的平顶分布，横向位置为高斯分布。可选归一化发射
度、Twiss alpha 和非相关 gamma 展宽共同给出高斯横向及能量相空间。

本工具不模拟注入器、压缩器、尾场、相干微聚束不稳定性、阴极不均匀或量子光子
发射。若上游相关性重要，应使用 Elegant 粒子输入。

## 为什么不能直接随机撒宏粒子

`N_e` 个独立物理电子在一个谐波上的散粒噪声满足

```text
sqrt(<|b|^2>) = 1/sqrt(N_e).
```

若把远少于物理电子数的 `N_macro` 个等权宏粒子完全随机地放在纵向相位上，数值
微聚束约为 `1/sqrt(N_macro)`，远大于真实 shot noise，进而错误地缩短 SASE 启动、
指数增益和饱和距离。

生成器采用 Fawley 的多 beamlet、多谐波 quiet-start 方法。每个共振波长单元含
`B` 个独立 beamlet；同一 beamlet 的成员具有完全相同的
`(x,y,ux,uy,gamma)`，仅纵向相位不同。要求正确加载至第 `H` 阶谐波时，每个
beamlet 使用 `2(H+1)` 个成员：

```text
u_quiet(j,k) = u0(k) + j*pi/(H+1),
delta_u(j,k) = sum[m=1..H] (
  a(m,k)*cos(m*u_quiet) + b(m,k)*sin(m*u_quiet)),
a(m,k), b(m,k) ~ Normal(0, 2/(m^2*N_beamlet)).
```

quiet-start 消除阶数取 `H+1`，使所有实际加载的谐波都处于 Fawley 方差公式的
内部阶数，同时压制会反向耦合到第 `H` 阶的相邻非物理谐波。不同 beamlet、不同
波长单元使用互相独立的随机系数。

正确的局部统计量为

```text
<N_cell*|b_m,cell|^2> = 1,  m=1..H.
```

生成文件会把每个目标谐波的实测值写入 HDF5 元数据。
只有 `1..H` 阶被有意赋予物理 shot-noise 统计；更高阶只是 quiet-start 结构，
不能当作物理预测。`H` 至少应覆盖 APRL 正式计算能够分辨或计划诊断的最高谐波，
并同时对 `H` 和 `B` 做收敛扫描。

## 横向与能量相空间

每个 beamlet 用独立标准正态变量生成横向高斯相空间。给定 rms 尺寸
`sigma_x`、归一化 rms 发射度 `epsilon_nx` 和 Twiss 参数 `alpha_x`：

```text
x  = sigma_x*r1,
ux = epsilon_nx/sigma_x*(-alpha_x*r1 + r2).
```

y 方向相同，并满足 `det(Cov[x,ux]) = epsilon_nx^2`。零发射度对应零横向固有
速度。gamma 按指定 rms 展宽生成，随后精确重构正向固有速度

```text
uz = sqrt(gamma^2 - 1 - ux^2 - uy^2).
```

## 构建与运行

该目标默认构建，也可在不改变 APRL 内核的情况下关闭：

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_COMPILER=mpic++ \
  -DAPRL_BUILD_SASE_BUNCH_GENERATOR=ON
cmake --build build --target sase_bunch_generator -j
```

示例：

```bash
./build/sase_bunch_generator \
  --output particles-sase.h5 \
  --electrons 6.24e8 \
  --length-m 30e-6 \
  --resonant-wavelength-m 100e-9 \
  --gamma 1000 \
  --sigma-gamma 0.2 \
  --sigma-x-m 30e-6 \
  --sigma-y-m 30e-6 \
  --emittance-n-x-m 1.0e-6 \
  --emittance-n-y-m 1.0e-6 \
  --harmonics 3 \
  --beamlets-per-wavelength 64 \
  --seed 12345
```

实际束长取最接近的整数个共振波长，并在终端摘要和 HDF5 元数据中同时报告。宏
粒子数不是独立输入，而由下式确定：

```text
N_macro = N_wavelength_cells * B * 2*(H+1).
```

生成 x 射线输入前必须先检查该数量；短波长、长束团可能产生很大的粒子文件。
工具以有界 HDF5 分块流式写出，默认拒绝已有文件，只有显式 `--overwrite` 才覆盖。

生成文件使用单位相对宏权重。APRL 输入卡必须重复完全相同的物理电子总数：

```yaml
beam:
  input:
    type: hdf5
    file: particles-sase.h5
    electrons: 6.24e8
    position_offset: [0.0, 0.0, -0.05]
```

HDF5-v2 是共同实验室时刻相对快照，因此该文件不使用
`beam.reference.input_plane_z`；`position_offset` 的 z 分量负责把束团放入 APRL
实验室束线坐标。

## 诊断与验证

`/particles` 组记录：

- 请求和实际束长；
- 代表电子数、平均电流、每波长单元和每 beamlet 的物理电子数；
- 波长单元、beamlet、成员和宏粒子总数；
- 相空间参数与随机种子；
- 每阶谐波的全局 `|b_m|`、理论全局 rms，以及平均
  `N_cell*|b_m,cell|^2`；
- 全部记录和诊断刷新后才置为 1 的 `complete`。

内置统计与谐波消除自检：

```bash
./build/sase_bunch_generator --self-test
```

正式 SASE 仍需使用多个独立 seed，并扫描每波长 beamlet 数。一个生成文件只代表
一个物理 shot，不是系综平均。

## 参考文献

1. C. Penman and B. W. J. McNeil, “Simulation of input electron noise in the
   free-electron laser,” *Optics Communications* 90, 82-84 (1992),
   DOI `10.1016/0030-4018(92)90333-M`。
2. W. M. Fawley, “Algorithm for loading shot noise microbunching in
   multidimensional, free-electron laser simulation codes,” *Phys. Rev. ST
   Accel. Beams* 5, 070701 (2002), DOI
   `10.1103/PhysRevSTAB.5.070701`。
