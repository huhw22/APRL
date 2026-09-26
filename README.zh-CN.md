# APRL——洛伦兹参考系中的加速器粒子与辐射

[English](README.md)

APRL 的英文全称为 **Accelerator Particles and Radiation in Lorentz
Frames**。APRL 是采用 C++/MPI 实现的加速器粒子与辐射模拟程序，主可执行
文件为 `aprl`。电磁推进直接使用国际单位制 E/B；标量势和矢量势不属于时间推进
状态。

## 适用范围

第一版提供以下计算路径：

- 相对论粒子与 Maxwell 场自洽推进；
- 沿指定 z 轴的 Cowan/CKC 可控色散传播；
- 紧凑 CFS-CPML 场边界和满足电流连续性的粒子开放边界；
- 由 Lorentz 变换接入的实验室系解析磁元件；
- HDF5 粒子输入、确定性高斯测试输入和基于官方 SDDS 库的 Elegant SDDS5
  转换；
- 固定实验室系场探测面与粒子探测面；
- 可在独立分析服务器运行的场面频谱、偏振、相干性和能量后处理；
- 用于少量粒子验证的轨迹辐射工具，以及可选的全局能量账本。

第一版未实现 Cowan/CPML 上的通用 TF/SF 激光或种子场注入。依赖该边界源的
seeded-FEL 与 laser-modulation 计算不属于已验证范围。

## 数值模型

- E/B 在交错 Yee 网格上以国际单位制推进。
- 默认 `cowan-z` 格式消除轴向真空传播的相位误差；标准 Yee 格式用于回归。
- 粒子使用相对论 Boris 推进；解析磁元件可启用粒子子步进。
- 初始粒子自场由分布式 CIC 相对论 Poisson 求解生成，并通过离散 Gauss 约束
  检查。
- 网格仅由整数单元数和单元尺寸定义；物理尺寸与 MPI 薄层偏移不依赖浮点单元
  计数。
- CPML 内表面同时是物理粒子边界。物理轨迹在该处终止，轻量且不输出的载流子
  在吸收层内完成连续电流衰减。

## 辐射计算流程

主要辐射输出为固定实验室系场探测面。正式运行保存原始自洽 Maxwell 场；解析
磁元件场不写入探测器文件。

场分析支持两条路径：

1. 直接分析原始场面，其中保留粒子束缚场；
2. 根据探测器左侧弹道参考记录或同位置粒子面，重构并扣除匀速粒子本底后分析。

粒子本底在 E/B 振幅层面扣除，随后计算 Poynting 功率和谱能量。弹道参考区不得
与磁作用区重叠。完整轨迹仅用于少粒子验证，不作为大规模宏粒子计算的主要输出。

独立后处理程序如下：

| 程序 | 功能 |
|---|---|
| `field_reconstruction` | 重构并扣除匀速粒子本底场 |
| `field_plane_analysis` | 计算角谱、能谱、Stokes 量和相干性 |
| `trajectory_radiation` | 根据少量粒子实验室系轨迹计算远场 |
| `field_power_compare` | 对实验性粒子退休路线执行信号/基线振幅比较 |
| `energy_closure` | 比较粒子面动能变化与场面辐射能量 |
| `energy_ledger_report` | 汇总运行时粒子/场能量账本 |

## 输入与输出

- 输入卡采用严格校验的 YAML 映射，未知键直接报错。
- 粒子输入采用版本化 HDF5，保存 SI 坐标、归一化固有速度、稳定粒子标识和正的
  相对宏粒子权重。
- 每次运行生成 YAML 清单，并在全部科学 HDF5 中写入统一的 `run_id`、配置摘要、
  源码版本和清单路径。
- 默认禁止覆盖已有科学输出；整体替换必须显式设置 `output.overwrite: true`。
- 探测器和能量账本仅由 MPI rank 0 写入；可选轨迹输出采用每 rank 一个文件。

## 运行模式

`runtime.mode: interactive` 用于本地测试，支持协调处理 SIGINT/SIGTERM 并关闭为
可读取的不完整输出。`runtime.mode: throughput` 用于批处理集群，关闭信号轮询和
周期耐久刷新。未启用的诊断不分配缓冲、不打开文件，也不进入对应通信路径。

可选资源监测器输出内存与文件规模估计、校准秒/步、运行进度、实测峰值驻留内存
和墙钟时间，格式可直接由批处理日志收集。

## 数值验证要求

必需回归测试覆盖配置读取、Lorentz 变换、Elegant 事件重构、Boris 收敛、离散
电荷连续性、Cowan 色散、CPML 反射、MPI 一致性、探测器与停止条件、HDF5 兼容、
输出保护、两条场分析路径和运行时能量账本数据链路。

回归通过只表示软件接口与数值路径一致，不表示具体物理问题已经收敛。正式结果
必须扫描场网格和时间步、宏粒子数、粒子子步、CPML 厚度与反射、横向孔径、探测
时间窗、初始场边界距离和能量账本残差。实验性粒子退休路线默认关闭，并要求匹配
零辐射基线。

## 构建与测试

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_COMPILER=mpic++ -DBUILD_TESTING=ON
cmake --build build -j
cmake --build build --target verify_required
./build/aprl config/example.yaml
```

多 rank 粒子输入默认要求 parallel HDF5。完整依赖、安装和集群运行流程见下列说明。

## 中文文档

### 运行与文件规范

- [构建、转换与运行](docs/zh_CN/BUILD_AND_RUN.md)
- [安装、卸载与依赖维护](docs/INSTALLATION.md)
- [YAML 输入卡规范](docs/zh_CN/YAML_CONFIGURATION.md)
- [粒子 HDF5 输入规范](docs/zh_CN/PARTICLE_INPUT_HDF5.md)
- [运行身份与覆盖保护](docs/zh_CN/RUN_PROVENANCE.md)
- [实验室系探测器 HDF5 输出](docs/zh_CN/DETECTOR_OUTPUT_HDF5.md)
- [实验室系轨迹 HDF5 输出](docs/zh_CN/TRAJECTORY_OUTPUT_HDF5.md)

### 数值方法

- [Cowan-z Maxwell 内核与 CPML 边界](docs/zh_CN/MAXWELL_COWAN_CPML.md)
- [满足 Gauss 约束的初始粒子自场](docs/zh_CN/INITIAL_SELF_FIELD.md)
- [粒子子步进](docs/zh_CN/PARTICLE_SUBCYCLING.md)
- [考虑 CPML 的粒子开放边界](docs/zh_CN/PARTICLE_OPEN_BOUNDARY.md)
- [场探测面的弹道参考区](docs/zh_CN/FIELD_DETECTOR_REFERENCE.md)

### 辐射与能量分析

- [粒子本底场重构](docs/zh_CN/FIELD_RECONSTRUCTION.md)
- [场面频谱与相干性分析](docs/zh_CN/FIELD_PLANE_ANALYSIS.md)
- [粒子轨迹远场工具](docs/zh_CN/TRAJECTORY_RADIATION.md)
- [粒子/场能量闭合诊断](docs/zh_CN/ENERGY_CLOSURE.md)
- [运行时粒子/场能量账本](docs/zh_CN/ENERGY_LEDGER.md)
- [实验室系能量诊断](docs/zh_CN/LAB_FRAME_ENERGY_DIAGNOSTICS.md)
- [实验性粒子退休与功率比较](docs/zh_CN/FIELD_POWER_COMPARISON.md)

### 测试与发布状态

- [测试框架](docs/TESTING.md)
- [数值验证记录](docs/zh_CN/VALIDATION.md)
- [发布状态与生产计算要求](docs/zh_CN/RELEASE_AUDIT.md)

示例输入位于 [`config/`](config/) 和 [`postprocess/`](postprocess/) 目录。
