# Elegant 粒子输入的 APRL 4 m 仿真

把 `elegant_sdds_to_hdf5` 转换出的 `particles.h5` 放在本目录。该文件必须是 HDF5 v4 的 Elegant 固定观察面事件，观察面位于原 RAD 波荡器**名义入口前 50 mm**。脚本直接读取此文件，不复制数 GiB 粒子数据。

波荡器保持原 `APRL_SASE_4m_weight1_2nodes` 设计：170 个 23.4 mm 周期（3.978 m）、`K=1.092714657005`、名义入口 `z=0`、高斯边缘场、1.35 nm 实验室纵向网格、104×104 横向网格、自动纵向框和相对论 Poisson 初场。默认向 `ssmb` 分区申请 2 节点 × 128 ranks，最多 144 小时。

粒子诊断面位于入口 `z=-0.03 m` 和波荡器内部 `z=1、2、3 m`。原出口 `z=4.03 m` 的普通粒子面改由场探测器自带的**弹道参考粒子面**承担，不再重复写一份同位置的完整粒子文件。场探测器保存未经近场扣除的原始 E/B 随时间变化，位于出口粒子面下游：其位置按输入粒子的最大 Lorentz γ 和横向计算域尺寸计算，约为 `z=4.89 m`，准确数值见每次生成的 `input-plan.json`。APRL 要求弹道参考面至场面的整个区间都在磁体边缘场外，所以两个面不能都放在 `z=4.03 m`。

运行目录的 `detectors/` 中会出现 `entry-particles.h5`、`after-1m-particles.h5`、`after-2m-particles.h5`、`after-3m-particles.h5`、`exit-field-ballistic-reference.h5`（出口粒子过面记录）和 `exit-field.h5`（原始 E/B 场，`/field_plane/time_s` 与 `/field_plane/fields`）。后两者供**仿真结束后**进行粒子近场重构；生成的 `reconstruction.yaml` 是后处理输入卡，作业不会自动运行 `field_reconstruction`。重构输出预设为 `detectors/exit-field-reconstructed.h5`。参考粒子文件使用 `/ballistic_reference` 格式，与普通粒子面的 `/particle_plane` 分组不同；现有分析程序若硬编码后者，需要相应调整。

场每个实验室 Maxwell 步采一次，约 `4.50 as`，对卡中声明的最高目标光子能量约有 8 次采样/周期。原始场覆盖整个 104×104 横向网格，单是场文件每记录 100 fs 就约需 **10.7 GiB**（未压缩）；对目前约一亿粒子的输入，四个普通粒子面、一个参考粒子面和原始场合计可能需要**数十至上百 GiB** 输出空间，实际大小取决于运行时采样数和过面粒子数。请先确认配额和文件系统余量。

Elegant 的输出面是过面事件，而不是同一时刻的粒子快照。APRL 会从该面向前投影到共同实验室时刻；脚本流式检查所有过面时间并把初始束团中心设在输出面下游的**最小合法位置再加 20 µm**。因此初始中心通常仍在名义入口前约 49.9 mm，但不严格等于 50.000 mm。它也在提交前检查束团头部上界是否位于边缘场及一个纵向网格安全距离之前。

若输入每条 `macro_weight` 都是 1，默认 `beam.input.electrons` 等于文件记录数。这对应先前宏粒子权重为一个电子的 Elegant 束团。若文件采用其他相对权重，提交前必须提供真实物理电子总数：`ELECTRONS=... ./submit.sh`。HDF5 的相对权重之和本身不一定等于物理电子数。

在登录节点，进入本目录并激活含 `numpy` 与 `h5py` 的 Python 环境后：

```bash
PYTHON_BIN=python3 ./submit.sh --plan
./submit.sh
```

`--plan` 只检查输入并生成 `runs/<时间>/run.yaml`、`input-plan.json` 与 `reconstruction.yaml`，不申请节点。正式运行输出在另一独立的 `runs/<时间>/`，其中 `aprl.log` 是 APRL 的标准输出，`detectors/` 保存粒子及场探测结果。每次计划或提交都会流式读一遍大文件，但内存占用受 262,144 条记录的块大小限制。提交后不要修改或替换 `particles.h5`；作业启动时会核对文件大小和修改时间。

如 APRL 安装位置不同，设置 `APRL_BIN=/absolute/path/to/aprl`。可沿用旧工况的 `MODULES="..."`、`NODELIST=...` 或 `EXCLUDE_NODES=...`。如果站点的 Slurm MPI 插件不是 `pmix_v3`，可设置 `SRUN_MPI=pmix_v4` 或 `SRUN_MPI=pmi2`。物理电子总数、输入 HDF5 格式与 WATCH 面位置应在正式提交前核对。

现在可显式选择每个 MPI 进程的计算线程数。默认仍为 **2 节点 × 128 MPI 进程/节点 × 1 线程**，原输入文件、网格、粒子数、诊断面和采样频率均保持不变。`THREADS_PER_RANK` 同时设置 Slurm 的 `--cpus-per-task` 与 `OMP_NUM_THREADS`；`NODES`、`TASKS_PER_NODE` 仅改变计算资源分配。每节点请求的 CPU 数为 `TASKS_PER_NODE × THREADS_PER_RANK`，必须适合实际节点资源。例如，在每节点 128 个可用 CPU 的分配上比较 2 个线程/进程：

```bash
NODES=2 TASKS_PER_NODE=64 THREADS_PER_RANK=2 ./submit.sh --plan
NODES=2 TASKS_PER_NODE=64 THREADS_PER_RANK=2 APRL_PROFILE_INTERVAL=100 ./submit.sh
```

使用支持计算线程的新 APRL 二进制。启动日志应显示 `Compute threads per MPI rank: 2`；仅设置环境变量不能使旧二进制产生并行计算。默认 `OMP_PLACES=cores`、`OMP_PROC_BIND=close`、`CPU_BIND=cores`，关闭 OpenMP 动态线程数，BLAS 等库保持单线程。需要按站点拓扑调节时可显式覆盖前述绑定选项；不要让 MPI 进程数乘线程数超过分配 CPU 数。

`execution-plan.txt` 记录请求布局；`execution-environment.json`、`cpu-topology.txt` 和 `aprl.log` 中各 rank 的 `[launch]` 行记录线程、网络选择和实际允许的 CPU。`APRL_PROFILE_INTERVAL=N` 每 N 步输出阶段耗时（0 或未设置关闭）。原有 `OMPI_MCA_btl=self,vader,tcp` 默认及用户覆盖机制保持原样，脚本不自动切换 UCX、网卡或网络协议。不同并行布局是否更快应看服务器实测；先用 APRL 源码内 `verification/performance` 的小型 SASE 加载测试验证新旧二进制的一致性。

作业结束后还会核对物理终止标记、诊断 HDF5 的完成/提交标记和非空记录数，并写出 `completion-check.json`。这一检查只读元数据，不扫描大型场数组；完整数值新旧对照由小型性能测试承担。
