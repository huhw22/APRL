# 新版离散能量账本：成功测试与验收说明

整理日期：2026-10-30。测试程序版本：`c7fcd20`；程序及输入的完整哈希保存在证据目录。

新版账本已完成单节点、双节点、已知场脉冲、三维制造场、故障识别和含粒子自由漂移的验收。诊断开关前后的原始场、粒子结果与既有账本科学字段逐位一致。含粒子小工况的场算子累计核对差约为 **{{FIELD_MAX}} J**，总核对差的绝对值为 **{{TOTAL_ABS}} J**，相当于净粒子能量变化绝对值的 **{{TOTAL_PERCENT}}%**。

这些测试验证的是**离散记账、程序操作和诊断开关的一致性**。时间离散与 Cowan 算子贡献仍明确保留，不把它们解释为辐射或真实外场做功。本文能量结果均在 **boost 系**，单位为 J。

## 1. 文件夹怎么使用

```text
README.md                    本说明，图片使用相对路径
figures/                     本文使用的成功测试结果图及 PDF
tests/01_single_node/        单节点成功验收入口
tests/02_two_nodes/          双节点成功验收入口
tests/03_reader_checks/      汇总读取器检查
evidence/                    完整原始输出、日志、冻结卡、验收 JSON
scripts/                     整理时的测试与汇总脚本快照
programs/                    本次测试程序入口及版本索引
metadata/                    迁移清单、文件身份与图片数据来源
```

测试结果实际存放在本目录的 `evidence/` 下；`tests/` 提供分类入口。工作区原路径保留兼容链接，冻结记录中的绝对路径仍可读取。本文和 `figures/` 可单独复制阅读；执行脚本仍依赖 APRL 工作区及原 MPI/HDF5 运行环境。

## 2. 已通过的测试矩阵

| 测试 | 输入与运行范围 | 通过标准 | 已完成配置 |
|---|---|---|---|
| 已知真空场脉冲 | 预置 E/B，无真实粒子、无沉积电流；28 步 | 与已知两格平移解、原生端点通量及离散账本相符 | 单节点 1／4 MPI；双节点 4 MPI |
| 三维制造场 | 非均匀 E/B、指定 J，无真实粒子；4 步 | Cowan 项被实际激发；逐点增量预测和场账本满足舍入界 | 单节点 1／4 MPI；双节点 4 MPI |
| 场值扰动负控 | 对一个 E 节点施加已知改动；1 步 | 改动必须被核对差和逐点预测检出 | 单节点 1／4 MPI；双节点 4 MPI |
| 含粒子自由漂移 | 31,675 个高斯型功能束团粒子，无外磁元件；257 步 | 完整结束、全粒子退出覆盖、诊断开关不改变科学输出、数值项目可独立核对 | 单节点、双节点各一对开／关，均为 16 MPI × 8 线程 |
| 既有软件回归 | 原 required 测试集 | 77 项全部通过 | 使用本次新程序重新执行 |
| 汇总读取器 | 已完成结果及隔离的派生 fixture | 采样步长正确；缺失、不完整、损坏或错误参考系必须拒绝 | 5 项全部通过 |

专项开关对照共 **4 次整程序运行**。制造场、脉冲、场值扰动共 **9 组模块测试**；77 项回归单列计数，不与前面的专项运行混算。

| 作业 | 节点 | 主要内容 | 结束状态 |
|---|---|---|---|
| 19027 | node3 | 77 项回归、6 组模块、单节点整程序开／关对照 | COMPLETED，0:0 |
| 19028 | node3、node5 | 3 组跨节点模块、双节点整程序开／关对照 | COMPLETED，0:0 |

原始入口：[单节点记录](tests/01_single_node/COMPLETION.json)、[双节点记录](tests/02_two_nodes/COMPLETION.json)。

## 3. 已知真空场脉冲

这一测试**先放入一个已知场脉冲，再让它传播**。没有粒子，也没有粒子电流。横向网格为 $160\times160$，纵向为 64 格，保持生产设计的 `dx=dy=5 µm`、`dz≈1.254326 µm` 和 `c·dt=dz`。核对孔径位于中央 $16\times16$ 单元。

场脉冲向 +z 传播，窗口每步向 −z 移一格，因此相对于窗口，已知解每步平移两格。场值与这一解逐步比较；原生面交换还与端点公式独立比较。三种布局的最大场平移相对偏差均约为 $2.2\times10^{-13}$。

![已知真空脉冲：场平移、面通量及场算子核对](figures/01_vacuum_pulse.png)

图中三项分别表示场值平移相对偏差、原生面交换核对偏差、场算子归一核对差。后两项的归一尺度以测试程序中的能量／交换尺度为准，不是粒子能量变化。[PDF](figures/01_vacuum_pulse.pdf)

详细记录：[单节点 1 rank](tests/01_single_node/execution/pulse-1rank/summary.json)、[单节点 4 rank](tests/01_single_node/execution/pulse-4rank/summary.json)、[双节点 4 rank](tests/02_two_nodes/execution/pulse/summary.json)。

## 4. 三维非均匀场与指定电流

该测试构造非均匀 E/B/J，使三维 curl、Cowan 修改项和 MPI 共享面都参与计算。它使用真实场推进程序，但 J 由测试函数指定，**不包含真实粒子**。纵格为 65，4 rank 时分片长度不完全相同，也检查了不等长分片。

通过条件包括：独立重算的 E/B 增量与实际增量相符；场算子账目在基于字段尺度的浮点舍入界内；Cowan 项非零，避免只验到退化的一维情形。

![三维制造场：增量预测与场算子核对](figures/02_manufactured_fields.png)

归一化场算子核对差约在 $10^{-18}$ 量级，电场增量预测的相对偏差约在 $10^{-16}$ 量级。[PDF](figures/02_manufactured_fields.pdf)

详细记录：[单节点](tests/01_single_node/execution/manufactured-4rank/summary.json)、[双节点](tests/02_two_nodes/execution/manufactured/summary.json)。

## 5. 故意改动场值，验证诊断能识别问题

测试在一个 E 节点上人为加入已知扰动，然后检查新账本能否识别。这里的“通过”指**故障被成功检出**，而不是接受被改动后的字段。

![场值扰动负控：正常状态与已知扰动的识别结果](figures/03_perturbation_detection.png)

人为扰动后，场算子归一核对差约为 $10^{-9}$，电场增量预测相对偏差约为 $10^{-4}$，与正常算术偏差明确分离。单 rank、4 rank 以及跨节点测试均成功识别。[PDF](figures/03_perturbation_detection.pdf)

详细记录：[单节点负控](tests/01_single_node/execution/injected-4rank/summary.json)、[双节点负控](tests/02_two_nodes/execution/injected/summary.json)。

## 6. 含粒子的无外磁场自由漂移

### 输入和设置

这一测试包含 **31,675 个真实推进的宏粒子**，不是无粒子场测试。输入为固定种子 `281028` 的高斯型功能束团：带截断、位置—动量相关，并设置两个尾部见证粒子，用来覆盖注入和入窗路径。横向位置—动量相关系数约为 −0.7334、−0.6546。

它是按既有束流矩参数构造的合成输入。此处的验收目标是自由漂移程序、粒子—场交换和账本的一致性；输入定义不应简化成理想冷高斯束团的物理基准。

| 设置 | 本次值 |
|---|---|
| 外部磁元件／入射波 | 均无 |
| 自场反馈 | 开启 |
| 参考场 | 既有快速准稳态参考方法 |
| boost γ | 929.1305492858543 |
| 活动场窗 | 160 × 160 × 289 |
| 活动窗 lower／upper 保护格 | 64／128 |
| 横向 CPML | x、y 各 10 格；z 无 CPML |
| 运行步数 | 257 |
| 整程序布局 | 16 MPI × 8 线程 |

计算卡：[单节点开启审计](tests/01_single_node/on/run.yaml)、[双节点开启审计](tests/02_two_nodes/on/run.yaml)。

### 核对结果

单节点和双节点分别运行严格审计关闭／开启两种配置。除了诊断开关和输出位置，物理设置保持相同。

| 输出 | 验收结果 |
|---|---|
| 原始场面 E/B | 开关前后逐位一致 |
| 三个粒子面 | 逐 ID 数据逐位一致，粒子覆盖完整 |
| 各 rank 轨迹 | 逐位一致 |
| 既有 v5 账本科学字段 | 逐位一致，计时字段单列 |
| 新场算子账本 | 累计核对差约 $10^{-31}$ J |
| 新总核对差 | 绝对值 {{TOTAL_ABS}} J，约净粒子能量变化绝对值的 {{TOTAL_PERCENT}}% |

![含粒子自由漂移：单节点与双节点的新账本核对结果](figures/04_particle_free_drift.png)

上图只展示本次新账本的场算子核对差；下图展示总核对差与独立粒子—网格交换差。单、双节点结果基本重合。这个小工况中，总剩余主要对应粒子—网格交换核对差。[PDF](figures/04_particle_free_drift.pdf)

证据：[单节点科学输出比较](tests/01_single_node/PHYSICS_BITWISE_REVIEW.json)、[双节点科学输出比较](tests/02_two_nodes/PHYSICS_BITWISE_REVIEW.json)、[v5 科学字段比较](tests/01_single_node/V5_LEDGER_BITWISE_REVIEW.json)。

日志：[单节点](tests/01_single_node/on/aprl.log)、[双节点](tests/02_two_nodes/on/aprl.log)。

## 7. 读取器与结果完整性

| 检查 | 成功结果 |
|---|---|
| 每 7 步采样 | 与完整记录中对应步逐条一致 |
| 缺少一个 rank 的文件 | 拒绝汇总 |
| rank 文件标记未完成 | 拒绝汇总 |
| 数值项被改坏 | 拒绝汇总 |
| 参考系被标为 lab | 拒绝按 boost 账本汇总 |

读取器测试使用单独的派生 fixture，不计作新的 PIC 仿真。原始输入保持不变。[读取器验收记录](tests/03_reader_checks/READER_REVIEW.json)

## 8. 启用成本

严格审计默认关闭。开启后会复制活动窗字段、独立重算算子，并进行额外 B 面交换和逐步写入。

| 配置 | 关闭时 s/步 | 开启时 s/步 | 本次增量 |
|---|---:|---:|---:|
| 单节点，16 MPI × 8 | {{SINGLE_OFF}} | {{SINGLE_ON}} | +{{SINGLE_COST}}% |
| 双节点，16 MPI × 8 | {{MPI_OFF}} | {{MPI_ON}} | +{{MPI_COST}}% |

![严格审计开关的耗时与内存测量](figures/05_audit_cost.png)

耗时取第 40–240 步，每种节点分配各一对，无 profiling，保留完整小工况输出；不外推为生产规模成本。内存为日志中采样到的同时总 RSS，不把各 rank 不同时刻的峰值相加。双节点申请独占两节点，整程序实际仍使用 128 个线程。[PDF](figures/05_audit_cost.pdf)

不开启时不分配新增场数组、不写新增 CSV、不发新增 B 面通信；本轮没有单独证明关闭后的严格零耗时。

## 9. 程序、开关与汇总入口

程序入口：[本次已验程序](programs/aprl)。完整身份见 [程序索引](programs/INDEX.json)。

在原卡的 `energy_ledger` 下加入：

```yaml
energy_ledger:
  enabled: true
  discrete_operator_audit: true
```

数据写入原 `energy_ledger.directory` 中的 `discrete-field-rank-*.csv`，并有完整标记。汇总入口为 [report.py](scripts/discrete_energy_ledger/report.py)。脚本快照仅供查阅和复现；实际运行当时的冻结脚本也完整保存在各测试目录的 `scripts/` 中。

现有启用条件为：移动窗、账本开启、未滤波电流、原横向 CPML 物理孔径；checkpoint/restart、入射波、横向窗外 ballistic 组合会明确拒绝。

原始验收总记录：[FINAL.json](evidence/discrete_energy_ledger_20261029/FINAL.json)。整理校验：[VALIDATION.json](metadata/VALIDATION.json)。
