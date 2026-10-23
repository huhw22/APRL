# APRL 计算性能优化：集群会话交接

## 用户目标及不可变条件

用户主要使用 ELEGANT 粒子输入的 APRL 脚本，主要性能测试参考 SASE 工况。目标是改善 MPI/OpenMP 在集群上的计算效率；不允许为了加速修改物理模型或计算精度。

新旧程序的粒子数、宏粒子权重、物理电子数、网格范围和间距、求解器、电流沉积、粒子子步、场反馈、边界、终止条件、探测器及输出频率、初场精度必须一致。不要把降低粒子数、缩短物理过程或减少诊断称为实现优化。不同规模的独立测试可以先规划，但每组新旧比较使用同一输入。

用户关注 held 冻结启动对负载不均衡的影响。冻结粒子减少单个 rank 的工作并不保证总时间下降；应观察最慢 rank、活跃/冻结粒子分布和阶段计时，而不是只看总粒子数。

## 代码与已完成的工作

- 远程项目根：`/public/home/wfollowm/APRL_Workspace`。先读根目录 `AGENTS.md` 和 `WORKSPACE_LAYOUT.md`。
- Git 远程：`git@git.tsinghua.edu.cn:huhw22/aprl_dev.git`。
- 候选分支：`wfollowm_dev/compute-performance-20261023`。准确交付提交见根目录 `HANDOFF.md` 的交付记录，并用 `git -C source/aprl rev-parse HEAD` 确认。
- 固定对照基线：`7a8ff35f5784df8868f24efe56f24132e4c5e678`，位于 `source/aprl-baseline`。

已实现网格/CPML 更新和独立粒子准备的 OpenMP 并行、held 不变速度缓存、工作缓冲复用、空/禁用功能跳过、减少每步统计集合通信，以及磁场 halo 通信与内部电场更新重叠。电流沉积和粒子事件仍保持原串行顺序；MPI 和 HDF5 调用保留在初始化线程。`APRL_PROFILE_INTERVAL` 控制分项计时，默认关闭。

ELEGANT 启动脚本已支持 `THREADS_PER_RANK`、匹配的 Slurm CPU 申请和线程绑定记录。其物理卡保持原样。4 m SASE 生产模板是原脚本的版本化副本，目前仍固定一线程；不能仅设置 `OMP_NUM_THREADS` 就声称它已经支持混合并行。需要混合布局时，统一修改提交资源与运行资源配置，并保持物理参数不变。

## 已有验证及其边界

本机已完成 required CTest 70 项、关闭 OpenMP 的构建与 10 项运行控制检查；held 全流程回归在 4 ranks、1/2 threads、10410 粒子条件下，13 组输出逐字节一致，包含 282 次跨 rank 解冻事件。同 rank 旧版对照也通过。

小型 SASE-style harness 已完成 24 次关闭 profiling、8 次开启 profiling 的本机运行；所有同 rank 新旧数值比较通过。默认输入 2560 条记录、19 步，使用 HDF5-v2，不启用 held；它用于集成及实现一致性，不是生产 SASE 增益实验。held 由单独回归覆盖。

本机无绑定短测试不是集群加速结论。旧基线已经存在跨 MPI rank 数的差异；报告中的 `implementation_equivalence_pass` 与 `cross_rank_consistency_pass` 必须区分，不能放宽容差掩盖差异。更多细节和历史本机数字见 `source/aprl/docs/PERFORMANCE.md`。本机报告的原二进制哈希不能用于新编译的集群程序。

本次交付阶段只推送、同步源码和测试脚本、建立工作区；没有因此提交大作业，也不声称完成了集群性能测量。服务器待重新编译并实测。

## 已核实的服务器环境

- MPI C++：`/usr/bin/mpic++`。
- CMake：`/public/home/wfollowm/software/cmake-3.31.12/bin/cmake`。
- Python：`/public/home/wfollowm/micromamba/envs/physics/bin/python3`；运行前再确认 `numpy`、`h5py` 可导入。
- 依赖前缀：`/public/home/wfollowm/aprl_dev/aprl-runtime-lab`。
- 环境脚本：该前缀下的 `share/aprl/runtime-env.sh`。
- Slurm 分区沿用 `ssmb`；`srun --mpi=list` 已列出 `pmix_v3`、`pmix_v4`、`pmi2` 等。默认先沿用原脚本的 `pmix_v3`；插件存在仍不等于多节点运行已验证。

不能混用不同 MPI 实现的编译器、运行器和并行 HDF5。保留 `Release`、`APRL_REQUIRE_PARALLEL_HDF5=ON`；新版开启 `APRL_ENABLE_OPENMP=ON`。核对 CMake 检测结果和运行日志实际线程数，不能仅凭环境变量判断 OpenMP 已启用。不要擅自改 UCX、网卡或 MPI 网络默认值。

## 新会话先做什么

1. 核对工作区、Git HEAD/状态和交付记录，确认旧目录和既有结果完整。
2. 读取 `cluster.env`，确认编译工具及 Python 依赖，核对 `sinfo`、CPU 拓扑及实际可用资源。
3. 用同一服务器工具链分别构建基线和候选，保留版本、构建参数、依赖和二进制哈希。
4. 先做小型同输入对照及线程绑定检查，再做完整 required 回归；通过后设计固定 CPU 预算的集群性能对比。
5. 完成小测试后再使用 SASE/ELEGANT 真实输入扩大规模。不要直接沿用生产脚本的 2 节点、256 ranks、48/144 小时默认值启动大作业。

只准备小测试的两个构建，不提交：

```bash
cd /public/home/wfollowm/APRL_Workspace
source ./cluster.env
bash ./tests/performance/prepare-and-submit.sh --prepare-only
```

小规模集成检查的提交入口（去掉 `--prepare-only` 才会提交）：

```bash
cd /public/home/wfollowm/APRL_Workspace
source ./cluster.env
PROFILE=16 REPEATS=1 bash ./tests/performance/prepare-and-submit.sh
```

该入口固定申请 1 节点、2 MPI slots、每 slot 2 CPU、30 分钟，测试 1/2 ranks × 1/2 threads。输出写入每次新建的 `results/performance/` 子目录；Slurm 日志位于 `tests/performance/`。再次运行使用新的结果目录，不复用已生成计划或已有结果。

吞吐测试使用同一物理输入，设置 `PROFILE=0 REPEATS=3`。需要不同资源布局时直接使用源码中的 `verification/performance/run.sbatch`，显式设置 `APRL_SOURCE_DIR`、`APRL_BASELINE`、`APRL_CANDIDATE`、`SASE_GENERATOR`、`PYTHON_BIN`、`PERFORMANCE_OUTPUT_DIR`，并按最大布局申请资源。测试矩阵的 ranks 和 threads 是笛卡尔积；不能把 `1,2,4 × 1,2,4` 当成只比较三组等 CPU 布局。

固定 CPU 预算可分别比较 4×1、2×2、1×4，或在实际允许 128 CPU/节点时比较 128×1、64×2、32×4。同 rank 增加线程会使用更多 CPU，必须与固定资源收益分开报告。多节点布局同时检查每 rank 的 CPU 绑定和节点分布。

快速入口使用 `BUILD_TESTING=OFF`，只构建 APRL 和生成器。完整回归另用 `build/verification`，开启 `BUILD_TESTING=ON` 并完成全量构建，在计算节点串行调度测试：

```bash
ctest --test-dir /public/home/wfollowm/APRL_Workspace/build/verification \
  -L required --output-on-failure -j 1
```

该命令只能在已准备好的 Slurm 分配内运行；required 测试含 4 ranks × 2 threads 的 held 回归，至少覆盖 8 CPU。MPI 测试启动器需符合本站运行环境；不能拿小性能入口的 2-rank 配额直接运行全套测试。

## 工况与结果判读

ELEGANT 工作副本在 `tests/elegant_4m`，`particles.h5` 经工作区 `inputs/` 链接到原 HDF5-v4 固定观察面数据。`cluster.env` 指定新 `APRL_BIN`。先核对输入语义及生成计划，不以 HDF5-v2 小 SASE 输入替代。生产诊断可能产生大量输出，扩大规模前检查实际空间与预计输出。

SASE 工作副本在 `tests/sase_4m`，环境覆盖 `APRL_PROGRAM`、`APRL_PROGRAM_DIR`、`SASE_GENERATOR` 和 `CAMPAIGN_WORK_DIR`，避免调用旧安装。不要照搬旧工作目录的已完成标记来跳过新程序验证。

验收需同时检查物理终点、完整且提交的 HDF5、有限数值、粒子守恒、终点粒子恰好记录一次，以及同 rank 新旧结果的一致性。单看程序退出码或文件存在不够。

`field_communication` 包含在 `field_total` 内，`particle_communication` 包含在 `particle_migration_total` 内；不要把父子耗时相加。MPI 计时包含等待慢 rank 的时间，不等于纯网络开销；不同阶段的最大值也可能来自不同 rank，不能相减重建精确墙钟时间。用 profiling 定位原因，用关闭 profiling 的重复运行给出吞吐结论。
