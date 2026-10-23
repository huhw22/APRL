# APRL 4 m SASE：1.35 nm、宏粒子权重 1 测试

这个目录是独立的单工况测试，不会修改上一轮四工况结果，也不会重新编译 APRL。它只运行：

| case | `dz_lab` | 波荡器 | 束团 | beamlets/λ | 宏粒子 | 总电子数 | 权重 |
|---|---:|---:|---:|---:|---:|---:|---:|
| `dz1p35-B0768` | 1.35 nm | 170 × 2.34 cm = 3.978 m | 40 µm top-hat SASE | 768 | 18,130,944 | 18,130,944 | 1 |

其余物理与数值设置保持上一轮 `dz1p35-B0192` 不变：相同 seed、横向网格、自动纵向框、64/128 保护格、20-cell z-CPML、相对论 Poisson 初场、能量账本、入口/出口粒子面以及两节点 MPI/watchdog。

## 为什么不是恰好 1800 万电子

40 µm 有效束长对应 2951 个共振波长单元，`harmonics=3` 的 quiet-start 每个 beamlet 有 8 条记录。为了让 256 个 MPI rank 严格均分，beamlets 必须取 32 的倍数。最近的严格单位权重点为：

```text
2951 × 768 × 8 = 18,130,944
```

因此总电荷和平均电流比 18,000,000 电子参考工况高 0.7275%。自动报告会同时给出原始场能比和按 `N_e²` 归一后的场能比，避免把这项小电荷变化误判成数值收敛。

## 运行前检查

服务器上的 APRL 应当已经重新编译，并同时提供：

```text
/public/home/wfollowm/aprl_dev/aprl-runtime-lab/bin/aprl
/public/home/wfollowm/aprl_dev/aprl-runtime-lab/bin/sase_bunch_generator
```

如果实际目录不同，修改 `campaign.env` 中的 `APRL_PROGRAM_DIR`，或提交时临时指定：

```bash
APRL_PROGRAM_DIR=/path/to/bin ./submit.sh --plan
```

首次放到服务器后执行：

```bash
cd ~/APRL_Project/APRL_SASE_4m_weight1_2nodes
chmod +x submit.sh analyze.sh prepare_python_env.sh run_campaign.sbatch scripts/*.sh scripts/*.py
./prepare_python_env.sh
./submit.sh --plan
```

`--plan` 的关键输出应当是：

```text
dz1p35-B0768: ... records=18,130,944 e/macro=1 e/beamlet=8 estimated_auto_nz=45,568
```

若 `records`、`e/macro` 或 `estimated_auto_nz` 不一致，不要提交。

随后先做同样的两节点小预检：

```bash
./submit.sh --preflight
```

确认预检完成后提交正式单工况：

```bash
./submit.sh
```

脚本默认让 Slurm 自行选择两个节点。需要指定节点时：

```bash
NODELIST=node4,node7 ./submit.sh
```

需要排除节点时：

```bash
EXCLUDE_NODES=node8 ./submit.sh
```

## 资源预算

基于上一轮同一 1.35 nm 场网格的 B96/B192 实测线性外推：

- 预计墙钟时间：约 18.9 小时，正常范围约 19–21 小时；默认 walltime 为 48 小时。
- 预计峰值总 RSS：约 82 GiB；保守按低于 100 GiB 估计。
- 输入束团约 1.1 GiB；入口和出口粒子面各约 1.7 GiB。
- 完整持久数据约 4.5–5 GiB；建议至少预留 10 GiB，若可能保留失败归档则预留 15–20 GiB。

场网格与 B192 完全相同，所以计算仍主要由场推进占据。新增风险主要是大粒子文件的共享文件系统 I/O，不是节点内存。

生成器会给出：

```text
warning: only 8 physical electrons per beamlet; Fawley's small-phase approximation is marginal
```

这是这个极限权重测试的预期警告，不表示任务失败。

## 看进度

提交后 `submit.sh` 会打印 job id 和日志路径。通常可查看：

```bash
squeue -j JOBID
tail -f logs/aprl-sase4m-w1-JOBID.out
```

正式 APRL 输出位于：

```text
work/runs/dz1p35-B0768/aprl.log
```

每 50 步的 `[resource]` 行可直接看到当前/峰值内存、实测单步耗时和剩余时间上界。

## 自动分析与主要判据

任务完成后会自动运行分析，也可以手工重跑：

```bash
./analyze.sh
```

结果位于 `work/results/`：

- `report.md`：主要结论及与旧 `dz1p35-B0192` 的定量对照；
- `case-summary.csv`：单工况汇总；
- `reference-comparison.csv`：内场、残差、时长和内存的 new/reference 比值；
- `comparison-diagnostics.png`：能量闭合、内场、平均粒子能量和非相关能散；
- `microbunching.png`：入口/出口的一波长 slice-RMS bunching；
- `resource-timeseries.csv`：运行时间和内存记录；
- `summary.json`：完整机器可读结果。

优先看三件事：

1. `peak abs(R)/E0` 是否低于 0.05；这是独立的能量闭合门槛。
2. 相对旧 B192 的 `final field energy normalized by Ne²` 是否明显小于 1；这直接检验权重 1 是否继续压低错误内场。
3. `elapsed time` 和 `peak total RSS` 的 new/reference 比值。运行时增长若约 10% 左右，就支持“当前工况仍由场网格主导，继续增加宏粒子成本不大”。

旧 B192 参考值已写入 plan，仅用于报告：峰值残差 0.0599181、末态内场能 3.27415 nJ、时长 17.1337 h、峰值总 RSS 78.727 GiB。它不会参与或影响新仿真的推进。

注意：这是单 seed 数值试验。即使内场和能量闭合改善，也不能单独作为 SASE 统计正确性的最终证明；而且本卡没有固定场探测面，不能用它认证远场辐射谱。

## 续跑与重跑

脚本会跳过指纹匹配且三个 HDF5 输出都完整的工况。若上次未完成，它会归档旧目录并从头重跑，不做 checkpoint 续算。

强制重跑：

```bash
FORCE_RERUN=1 ./submit.sh
```

使用另一个工作目录做隔离测试：

```bash
CAMPAIGN_WORK_DIR="$PWD/work_retry" ./submit.sh
```
