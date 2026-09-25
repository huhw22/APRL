# 粒子 HDF5 输入规范

## 用途与坐标归属

这是正式计算使用的粒子输入契约。格式紧凑，可由 MPI rank 直接切片读取，并且不依赖 YAML 的显示单位。版本 4 描述粒子穿过一个固定 Elegant 观察面的事件，并保留其到达时间坐标。YAML 输入卡把该平面和重构后的束团参考位置放在同一个用户定义实验室坐标系中；不要求任何磁元件或探测器定义全局 `z=0`。

## 必需布局

```text
/particles                         group
  @format_version = 4              signed integer attribute
  /records                         1-D compound dataset, N > 0
```

`/particles/records` 的每条记录包含：

| 字段 | HDF5 类型 | 含义 |
|---|---:|---|
| `plane_position_m[2]` | 小端 float64 | 固定输入面上的 x、y，单位 m |
| `arrival_time_offset_s` | 小端 float64 | 过面时刻减去文件加权参考到达时刻，单位 s |
| `proper_velocity[3]` | 小端 float64 | 无量纲 `gamma*v/c` |
| `source_id` | 小端 uint64 | 稳定的源记录标识 |
| `macro_weight` | 小端 float64 | 正的相对宏粒子权重 |

所有浮点量必须有限，`macro_weight` 必须为正。`source_id` 不要求连续，但在重新生成输入时应保持稳定，以便后续轨迹分析连接记录。文本转换器使用从 1 开始的输入行号；SDDS 转换器在存在 `particleID` 时原样保留，否则使用所选页内行号。

转换器由 Elegant 标准 `(xp,yp,p)` 列精确计算固有速度，不使用傍轴近似：

```text
uz = p/sqrt(1+xp^2+yp^2)
ux = xp*uz
uy = yp*uz.
```

其中 `p=beta*gamma` 为动量模。启动时，模拟器用每条记录自己的速度和过面时间构造共同实验室时刻快照。共同时间的选择使宏权重加权纵向质心位于 `beam.reference.initial_center_z`，并保留到达时间与纵向速度的相关性：

```text
T = [(initial_center_z-input_plane_z)/c + <beta_z*tau>_w] / <beta_z>_w
r_i(T) = r_i(tau_i) + beta_i*c*(T-tau_i),
```

其中 `tau=arrival_time_offset_s`。每个粒子都必须满足正 `uz` 和 `T-tau_i>=0`；否则程序在场分配前退出，并报告最小合法中心位置。该弹道同步是模拟器唯一执行的 Elegant 到快照输运，它不会重新追踪上游 Elegant 晶格。

版本 3 仍可读取旧式 `position_m=(x_plane,y_plane,zeta)` 表示，使用参考速度近似 `zeta=-v_reference*(t-t_ref)`，并以 `ux/uz`、`uy/uz` 推进横向坐标。版本 1、2 是旧的共同时间相对快照；版本 1 没有 `macro_weight`，按单位权重解释。版本 1、2 不使用 `input_plane_z`。新的 Elegant 转换应使用版本 4。

SDDS 转换文件还保存所选页、协议版本、源数据的二进制/ASCII 模式、绝对参考到达时刻，以及 ID 和权重选择说明。读取器以格式版本和复合字段的名称/类型为准，描述属性仅供检查。

## MPI 读取行为

对 N 条记录和 P 个 rank，程序给每个 rank 分配一个连续且尽量等长的 hyperslab，并通过 collective MPI-IO 读取；rank 0 不会读取并重新分发整个文件。因此多 rank 运行必须使用 parallel HDF5。串行 HDF5 只允许单 rank，其他情况会明确失败。

YAML 的 `beam.input.electrons` 表示完整数据集代表的物理电子总数。版本 2–4 中记录 i 代表

```text
electrons_i = beam.input.electrons * macro_weight_i / sum(macro_weight).
```

全局权重和通过 collective 计算。宏粒子电荷与质量同比缩放，保持物理荷质比不变。启动时还以扩展精度独立求和所代表电子数；若归一化或荷质比不再保持则退出。输出 `charge_C` 是权威的代表电荷，输出 `weight` 只保留输入相对权重。

版本 4 中，横向 `position_offset` 作用在 Elegant 输入面；其 z 分量在同步后加到重构相对快照上，避免在文件读取器中混用秒和米。

## 直接转换 Elegant SDDS

`elegant_sdds_to_hdf5` C++ 工具链接官方 SDDS C 库，直接读取 Elegant 原生文件，不生成或解析中间文本。必需列为 `x[m], xp, y[m], yp, t[s], p[m$be$nc]`。标准可选 `particleID` 会被精确保留，包括 SDDS5 的 `long64`/`ulong64`；缺失时使用所选页从 1 开始的行号。Elegant 通常使用等电荷宏粒子，所以默认权重为 1；若源文件确实含逐行正相对权重，可显式选择其列。

SDDS 文件可包含多页，例如不同 pass、step 或 bunch。转换器有意只写一页，默认第一页，`--page` 用于显式选择其他页，绝不静默拼接。`SDDS5` 指协议版本 5，而不是数据模式；`&data mode=binary` 与 `mode=ascii` 都合法。官方库负责协议、行/列优先、字节序、压缩和页面布局。

```bash
./build/elegant_sdds_to_hdf5 \
  --input elegant-watch.sdds \
  --page 1 \
  --output particles.h5

# 仅当源文件确实具有不等逐行权重时：
./build/elegant_sdds_to_hdf5 \
  --input elegant-watch.sdds \
  --output particles.h5 \
  --weight-column macroparticleWeight
```

Elegant 页级 `Charge` 参数不是逐粒子权重。总物理粒子数仍由 YAML 的 `beam.input.electrons` 权威定义；只有真实逐行相对权重列才能传给 `--weight-column`。

## 文本转换输入

旧工具 `particle_text_to_hdf5` 仍接受：

```text
x_plane  y_plane  zeta  ux  uy  uz  [macro_weight]
```

位置列使用命令行 `--length-unit`；固有速度始终为 `gamma*v/c`。可选第七列为正相对权重，省略时为 1。工具为了向后兼容写格式版本 3。它验证每个数值、拒绝多余列、忽略空行和行尾 `#` 注释，并分块有界写入。命令见[构建、转换与运行](BUILD_AND_RUN.md)，最小示例见 `examples/particles/example_particles.txt`。
