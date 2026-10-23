# 通用粒子计算性能标定

本基准使用 APRL 内置 `generated-gaussian`，不读取外部粒子文件，也不复算大型 ELEGANT/SASE 生产工况。每张卡定义一个独立的、自洽高斯电子束漂移计算，随机种子固定为 20261023；同一张卡在基线、候选和所有 MPI/OpenMP 布局间逐字节一致。

| 卡 | 宏粒子数 | 网格单元 | 目的 |
|---|---:|---|---|
| balanced | 262144 | 48×48×1024 | 首次计时校准及综合计算 |
| field-heavy | 32768 | 64×64×2048 | 增大场更新工作量 |
| particle-heavy | 1048576 | 32×32×1024 | 增大粒子推进及沉积工作量 |

三个工况都表示 100000 个物理电子，使用 CPML、自洽场反馈、相对容差 1e-9 的初场 Poisson 求解及固定探测器和能量记录。不同卡的粒子权重、网格和工作量不同，不能把它们之间的时间比当成版本加速。运行不得修改已经测量的卡；需要新规模时创建新的独立卡。

覆盖范围包括初始化、普通粒子推进、场插值、电流沉积、Maxwell/CPML 更新、MPI 粒子/场通信和 HDF5 诊断。粒子集中于束团附近，沿 z 切分时存在真实负载不均衡；它不代表均匀分布的理想粒子扩展性，也不涵盖 ELEGANT 文件解析、held 注入、外部磁铁或激光等所有可选路径。此前 required/held 回归另行覆盖这些实现的一致性。本测试不证明物理收敛或生产 SASE/ELEGANT 加速。

## 运行顺序

先接受工作区已有 smoke 和 required 结果，再在计算节点执行。模板部署到 `tests/generic_particles/`，另生成 `prepared-manifest.sha256`，同时沿用 `tests/cluster_acceptance/` 的构建和运行库清单。不得覆盖历史结果目录或修改原验收清单来绕过校验。

`run.sbatch` 申请一个独占 128 CPU 节点、32 GiB 内存，作业限时 30 分钟；每个探针或模拟另限时 600 秒。首次校准用 balanced、`GENERIC_MODE=calibration`、`GENERIC_LAYOUTS=128x1`。若物理终点、输出及同 rank 对照通过且时间适中，再做 `GENERIC_MODE=throughput`（默认 3 次重复），并单独运行 `GENERIC_MODE=profile`（候选计时周期 64，默认 1 次）。不要批量提交尚未校准的大测试。

提交前设置工作区 `cluster.env`，并设置以下环境变量：

- `APRL_SMOKE_RESULT`：已接受的 smoke 绝对结果目录。
- `APRL_REQUIRED_RESULT`：已接受的 required 绝对结果目录。
- `GENERIC_CASE`：本目录某张 YAML 卡部署后的绝对路径。
- `GENERIC_MODE`：`calibration`、`throughput` 或 `profile`。
- `GENERIC_LAYOUTS`：默认 `128x1,64x2,32x4`，各候选布局固定 128 CPU。
- `GENERIC_REPEATS`：可选；正式吞吐至少 3 次。

所有步骤先检查实际 MPI/OpenMP 核绑定。新结果自动写入 `results/cluster_acceptance/generic-particles-<job>-<time>-<unique>/`，失败也保留清单、启动日志和部分报告。`run.py --plan` 仅创建新计划目录，既不生成粒子也不运行程序，可在管理节点审阅。

## 判定与计时

新旧程序路径和哈希来自独立构建记录。候选使用 128×1、64×2、32×4；每种 rank 数另有单线程基线作为正确性控制。固定资源速度比仅以基线 128×1 为参照；64×1、32×1 基线消耗的 CPU 更少，其速度比单独标为固定 rank 加线程。

每次运行必须达到卡上的物理终点，HDF5 输出完整、提交长度匹配、数值有限；粒子编号恰为 1..N 且只记录一次，活动粒子数量不变且无损失。新旧同 rank 的所有数值输出要求逐字节相同，重复运行亦同。跨 rank 差异另用原有 rtol=1e-9、scaled_atol=1e-12 报告，不能放宽或和实现一致性混淆。

同时报告包含初始化的 APRL 总墙钟、包含最终输出关闭的模拟循环墙钟，以及每步时间。粒子/网格吞吐由逻辑步数归一化，不是 FLOPS 或逐核计数。启动器进程墙钟另列。只有关闭 profiling 的重复运行用于吞吐判断；单次校准不作为稳定加速结论。

Profiling 中通信时间包含等待其它 rank；`field_communication` 已含于 `field_total`，粒子通信已含于迁移总时间，不能重复相加。阶段最大值可能来自不同 rank，不能相加重建总墙钟。记录粒子工作量 min/max，以区分 OpenMP 效果与束团沿 z 的负载不均衡。
