# 中文文档总览

本目录给出 APRL 第一版全部英文说明文件的中文对应版本。APRL 的英文全称为 “Accelerator Particles and Radiation in Lorentz Frames”，中文名称为“洛伦兹参考系中的加速器粒子与辐射”。

安装与测试流程原文已经使用中文，因此直接维护唯一版本，避免双份文档随版本演进发生偏差：

- [安装、卸载与依赖维护](../INSTALLATION.md)
- [分层测试与数值分析流程](../TESTING.md)

其余说明按主题列在下面。英文原文仍保留在上一级 `docs/` 目录中，配置键、命令、文件名和 HDF5 数据集名不翻译，以便输入卡、程序输出和说明严格对应。

## 首次使用

1. [构建、粒子转换与运行](BUILD_AND_RUN.md)
2. [YAML 输入卡完整规范](YAML_CONFIGURATION.md)
3. [粒子 HDF5 输入规范](PARTICLE_INPUT_HDF5.md)
4. [发布状态与生产计算要求](RELEASE_AUDIT.md)
5. [数值验证状态](VALIDATION.md)

## 推进内核与边界

- [Cowan-z Maxwell 内核与 CPML 边界](MAXWELL_COWAN_CPML.md)
- [粒子子步进](PARTICLE_SUBCYCLING.md)
- [初始粒子自场](INITIAL_SELF_FIELD.md)
- [考虑 CPML 的粒子开放边界](PARTICLE_OPEN_BOUNDARY.md)

## 输出、探测器与重构

- [运行身份、清单与覆盖保护](RUN_PROVENANCE.md)
- [实验室系探测器 HDF5 输出](DETECTOR_OUTPUT_HDF5.md)
- [实验室系粒子轨迹 HDF5 输出](TRAJECTORY_OUTPUT_HDF5.md)
- [场探测面的弹道参考区](FIELD_DETECTOR_REFERENCE.md)
- [粒子本底场重构](FIELD_RECONSTRUCTION.md)
- [场探测面频谱与相干性分析](FIELD_PLANE_ANALYSIS.md)
- [粒子轨迹到远场辐射工具](TRAJECTORY_RADIATION.md)
- [粒子退休与匹配功率比较](FIELD_POWER_COMPARISON.md)

## 能量诊断与物理解释

- [粒子/场能量闭合诊断](ENERGY_CLOSURE.md)
- [运行时粒子/场能量账本](ENERGY_LEDGER.md)
- [实验室系能量诊断](LAB_FRAME_ENERGY_DIAGNOSTICS.md)

## 文档维护约定

- 说明文件采用客观技术陈述，不记录对话、开发过程或未确定方案；未实现能力只列为适用范围限制或预留接口。
- 改动输入键、默认值、硬停止或警告条件时，必须同步更新中英文 `YAML_CONFIGURATION.md`。
- 改动 HDF5 布局、单位或属性时，必须同步更新对应输入/输出规范。
- 改动数值方法、适用边界或已知限制时，必须同步更新内核说明、验证状态和发布审计。
- 测试结果只记录已经实际运行并保留复现实例的结论；计划中的测试不得写成已验证能力。
- 修改正式名称或缩写时，必须同步更新源码命名空间、文件格式标识、构建选项、安装路径与双语文档。
