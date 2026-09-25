# APRL——洛伦兹参考系中的加速器粒子与辐射

[English](README.md)

**APRL** 的英文全称是 **Accelerator Particles and Radiation in Lorentz
Frames**，中文名称为“洛伦兹参考系中的加速器粒子与辐射”。名称直接概括了
程序的三个核心对象：加速器粒子、辐射与 Lorentz 参考系模拟。主可执行文件名为
`aprl`。

## 基本思想

这是一个面向自由电子激光、预聚束电子辐射及相关相对论束流问题的轻量级 C++/MPI 程序。程序在交错 Yee 网格上直接以国际单位制推进电场与磁场；主 Maxwell 内核采用沿 z 方向优先的 Cowan/CKC 可控色散格式，使配置传播轴上的真空波精确无色散，标准 Yee 格式保留为回归选项。粒子由自洽网格场推进，实验室系解析磁元件和未来的注入场通过明确的 Lorentz 变换进入计算系；A/phi 不参与推进状态。

网格几何仅以整数网格数 `mesh.cells` 与单元尺寸 `mesh.cell_size` 为准，总尺寸通过乘法得到。MPI 沿 z 的不等长分区只使用整数商和余数构造，避免由“总长度除以分辨率”引入的不稳定取整。

主要辐射输出是固定在实验室系的场探测面。场探测面可以显式启用一个仅用于诊断的左侧直线参考区：真实粒子继续接受完整自洽推进，而每个粒子只记录一次过面状态，供后处理重构并扣除匀速带电粒子本底。小规模验证还可以在探测面再次采样真实粒子，与直线外推逐粒子配对比较，而无需开启完整轨迹。粒子探测面与上述机制相互独立。

完整粒子轨迹只用于少粒子验证和调试，因为对约 `10^7` 个宏粒子的正式计算，其存储和读写代价不可接受。大规模模拟只负责轻量推进与必要诊断，频谱、相干性和远场后处理可以在较便宜的分析服务器上独立运行。

所有科学输出默认禁止覆盖。每次运行生成一份清单，并在全部 HDF5 输出中写入统一的 `run_id`、配置摘要和源码版本，详见[运行身份、清单与覆盖保护](docs/zh_CN/RUN_PROVENANCE.md)。

项目提供五个独立 C++ 后处理器：

- `trajectory_radiation`：从少量粒子的实验室系轨迹计算带偏振的远场角谱、积分能谱、Stokes 量及可选互谱密度；
- `field_reconstruction`：结合原始实验室 E/B 探测面与直线参考或同位置粒子面，以三维 FFT 重构匀速带电粒子场并写出扣除后的场；
- `field_plane_analysis`：从原始场或重构场提取前向传播角谱、积分能谱、Stokes 量、相干/涨落分解、空间与双频互谱密度；
- `field_power_compare`：在实验性粒子退休路线中，对信号与零辐射基线先做 E/B 振幅相减，再比较瞬时功率、频段能量与少粒子轨迹远场；
- `energy_closure`：以粒子 ID 配对实验室系粒子探测面，使用逐粒子 gamma 差和匹配 `K=0` 基线分析粒子动能变化与场能量。

另一条明确标记为实验性的路线会在磁系统后停止真实粒子推进，并用短程数值载流子把电流平滑衰减到零。它不修改 Maxwell 内核，默认关闭，而且由于域内净电荷消失并非严格连续，必须配合零辐射基线和收敛研究使用。频率保护模式允许手动指定最低受保护光子能量和过渡周期数，并在启动时检查退休长度及因果保护距离。

## 适用场景

- 预聚束电子在波荡器等磁元件中的辐射；
- 无注入激光时的 FEL、相干辐射和高次谐波数值研究；
- 相对论变换系中的电子运动与自洽 Maxwell–粒子模拟；
- 面向远场、频谱、偏振和空间相干性的实验室系场面输出；
- 直接 SI E/B 格式、Cowan 色散控制和 CPML 边界的数值实验。

当前第一版的目标范围不包括 Cowan+CPML 内核上的通用激光/种子场 TF/SF 注入。该功能预留给后续版本；在其完成前，不能把本版描述为完整的 seeded-FEL 或 laser-modulation 求解器。

## 运行与资源原则

运行可以在两种停止条件中选择：所有仍有效粒子越过最后一个有限作用区，或实验室系参考中心到达指定 z。场面和粒子面都属于束线元件；场面的直线参考区仅位于其左侧，粒子面在几何上始终为零长度。

`runtime.mode: interactive` 适合小服务器调试，能够响应 SIGINT/SIGTERM 并关闭成可读取但 `complete=0` 的输出；`throughput` 面向批处理超算，去掉信号轮询和周期耐久刷新。场面、粒子面和能量账本只由 MPI rank 0 写入，轨迹则每个 rank 一个文件，避免多 rank 抢写同一 HDF5 文件。

粒子输入支持正的相对宏粒子权重，并按输入卡的总电子数归一化；测试用高斯束团由总电子数与宏粒子数生成等权重样本。原生 Elegant SDDS5 转换器依赖官方 SDDS 库，可直接保留固定平面到达时间、斜率、动量、粒子 ID 和可选权重，并写成 HDF5 v4，主程序随后把各粒子事件同步到共同实验室时刻再做 Lorentz 变换。

可选能量账本记录相对论变换系中的活动/逃逸粒子动能、物理区 E/B 储能、六个 CPML 内表面通量、解析磁元件做功以及实验室系能散统计。关闭后不增加循环、集合通信和文件；开启后只产生一个小型 rank-0 HDF5 流。

## 数值状态与限制

- 无种子 Cowan 运行支持紧凑 CFS-CPML；通用 Cowan/CPML 种子场注入尚未实现。
- 可选的分布式 CIC/相对论 Poisson 初始化能够生成满足离散 Gauss 定律的初始 E/B，自此不保留任何 A/phi 状态；仍需对宏粒子数、边界距离和求解容差做收敛研究。
- Boris 粒子子步进只细化解析磁元件的轨道采样，不提高 Maxwell 带宽，也不增加探测器 Nyquist 能量。
- CPML 内表面同时是物理粒子边界。真实轨迹在此终止，随后由不输出的弹道载流子在 CPML 内平滑电流并于外边界完成残余电荷清理。
- 单个前向场面不能独立闭合全局能量账本；正式结果必须同时研究储能、其余边界通量、孔径、网格、CPML、初始场和宏粒子收敛。
- 粒子退休法是实验功能；稳定的第一版主路线是直接分析原始场，或使用直线粒子本底重构后再分析。

## 中文文档导航

- [构建、转换与运行](docs/zh_CN/BUILD_AND_RUN.md)
- [安装、卸载与依赖维护](docs/INSTALLATION.md)（原文已为中文）
- [分层测试和数值分析流程](docs/TESTING.md)（原文已为中文）
- [YAML 输入卡完整规范](docs/zh_CN/YAML_CONFIGURATION.md)
- [运行身份、清单与覆盖保护](docs/zh_CN/RUN_PROVENANCE.md)
- [Cowan-z 内核与 CPML 边界](docs/zh_CN/MAXWELL_COWAN_CPML.md)
- [满足 Gauss 定律的初始粒子自场](docs/zh_CN/INITIAL_SELF_FIELD.md)
- [实验室系能量诊断与 50 A 尺度估算](docs/zh_CN/LAB_FRAME_ENERGY_DIAGNOSTICS.md)
- [粒子子步进及其场步长限制](docs/zh_CN/PARTICLE_SUBCYCLING.md)
- [考虑 CPML 的粒子开放边界](docs/zh_CN/PARTICLE_OPEN_BOUNDARY.md)
- [粒子 HDF5 输入规范](docs/zh_CN/PARTICLE_INPUT_HDF5.md)
- [实验室系轨迹 HDF5 输出](docs/zh_CN/TRAJECTORY_OUTPUT_HDF5.md)
- [粒子轨迹远场工具](docs/zh_CN/TRAJECTORY_RADIATION.md)
- [实验室系探测器 HDF5 输出](docs/zh_CN/DETECTOR_OUTPUT_HDF5.md)
- [场探测面的弹道参考区](docs/zh_CN/FIELD_DETECTOR_REFERENCE.md)
- [粒子本底场重构](docs/zh_CN/FIELD_RECONSTRUCTION.md)
- [场面频谱与相干性分析](docs/zh_CN/FIELD_PLANE_ANALYSIS.md)
- [粒子/场能量闭合诊断](docs/zh_CN/ENERGY_CLOSURE.md)
- [运行时粒子/场能量账本](docs/zh_CN/ENERGY_LEDGER.md)
- [粒子退休与匹配功率比较](docs/zh_CN/FIELD_POWER_COMPARISON.md)
- [数值验证状态](docs/zh_CN/VALIDATION.md)
- [当前发布审计与生产门禁](docs/zh_CN/RELEASE_AUDIT.md)

英文文档仍是同步维护的另一语言版本。项目名、主可执行文件名与安装目录均统一使用 APRL 标识。
