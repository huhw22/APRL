# APRL 集群验收：用户手动提交

本目录只提供分阶段作业脚本，不会自动提交。当前工作区为
`/public/home/wfollowm/APRL_Workspace`，实际使用副本位于
`tests/cluster_acceptance/`。每个作业在 `results/cluster_acceptance/` 下新建独立目录。
作业输出会打印完整目录；成功时才创建 `PASS`，不能仅凭 Slurm 退出码或文件存在验收。

## 当前准备状态

- 基线计算程序来自 `7a8ff35f5784df8868f24efe56f24132e4c5e678`。
- 候选计算程序来自 `ac58e67e33d6b4ecfd12c3d00f4a9947d0e8500d`。
- 三套构建分别为 `build/baseline`、`build/candidate`、`build/verification`。
  候选与回归构建的 APRL 可执行文件 SHA256 相同。
- GCC 9.4.0、Open MPI 4.0.3、并行 HDF5 1.14.6、CMake 3.31.12；Release，候选启用 OpenMP。
  Python 为 `cluster.env` 中的 physics 环境；其 h5py 自带 HDF5 用于独立结果分析进程。
- held 比较器增加 HDF5 新版 `H5Ovisit` 参数兼容，未更改科学数据比较规则、容差或计算实现。
- 已完成编译、脚本语法和计划检查；required 清单为 70 项。尚未执行任何集群科学测试。
- 新增模板及比较器修复独立提交，因此源码当前 HEAD 会前进；待测计算程序仍为上述原始交付版本。
  `build-record.json` 分别记录编译核心版本、当前源码版本、绝对程序路径、依赖和哈希。

## 第一阶段：现在只提交这一项

```bash
cd /public/home/wfollowm/APRL_Workspace
source ./cluster.env
sbatch --chdir="$APRL_WORKSPACE/tests/cluster_acceptance" \
  "$APRL_WORKSPACE/tests/cluster_acceptance/smoke.sbatch"
```

资源为单节点 2 tasks × 2 CPU、8 GiB 内存，30 分钟。逐 rank、逐线程检查实际 CPU/物理核绑定，
随后运行原有小型 SASE-style harness：1/2 ranks × 1/2 threads，1 次重复，profile=16。
它使用原有 2560 条记录、19 步的独立集成工况；没有改变真实 ELEGANT/SASE 的物理卡。
检查物理终点、HDF5 完整提交、有限值、粒子守恒和同 rank 逐位一致；已有跨 rank 差异另列。

运行结束后将作业号或结果目录交给 Codex，审阅 `driver.log`、绑定记录及
`harness/report.json`，再继续下一阶段。不要提前串联提交后续作业。

## 第二阶段：审阅 smoke 后手动提交

```bash
cd /public/home/wfollowm/APRL_Workspace
source ./cluster.env
export APRL_SMOKE_RESULT=/绝对路径/第一阶段实际结果目录
sbatch --chdir="$APRL_WORKSPACE/tests/cluster_acceptance" \
  "$APRL_WORKSPACE/tests/cluster_acceptance/required.sbatch"
```

单节点 4 tasks × 2 CPU、16 GiB 内存，45 分钟。先检查 4×2 绑定，再在分配的 CPU 上串行调度全部
70 项 required（包含候选 held 4 ranks × 1/2 threads），追加基线 4×1 的 held 对照。
三次 held 使用同一原始卡和 10410 粒子输入；检查两个粒子探测面、完整输出、
40 µm 物理终点、粒子唯一记录与守恒，不缩短过程、不降低精度。

CTest 的输出路径在配置时写入构建目录，故脚本通过 `.cluster-run-started` 实施单次使用。
退出时复制 `verification/work` 和 `Testing` 到新结果目录，失败也保留记录。
若要重试，由下一轮工作准备新的回归构建和对应清单；不得删除标记后覆盖旧输出，
也不得移动运行中的目录。

## 第三阶段：前两阶段审阅通过后再测试固定 CPU 布局

同样先设置 `APRL_SMOKE_RESULT`，再设置：

```bash
export APRL_REQUIRED_RESULT=/绝对路径/第二阶段实际结果目录
```

下面是三种独立手动提交方式，每次 4 CPU、8 GiB 内存，profile=0，每组重复 3 次。建议一次提交一组，
审阅后再提交下一组；不要用 ranks 与 threads 的笛卡尔积替代它们。

```bash
sbatch --ntasks=4 --cpus-per-task=1 --chdir="$APRL_WORKSPACE/tests/cluster_acceptance" \
  "$APRL_WORKSPACE/tests/cluster_acceptance/fixed-4cpu.sbatch"
sbatch --ntasks=2 --cpus-per-task=2 --chdir="$APRL_WORKSPACE/tests/cluster_acceptance" \
  "$APRL_WORKSPACE/tests/cluster_acceptance/fixed-4cpu.sbatch"
sbatch --ntasks=1 --cpus-per-task=4 --chdir="$APRL_WORKSPACE/tests/cluster_acceptance" \
  "$APRL_WORKSPACE/tests/cluster_acceptance/fixed-4cpu.sbatch"
```

脚本核验前两阶段 PASS 和同一 manifest。每组的一线程基线只作为同 rank 一致性控制；
固定 4 CPU 性能参考只能使用 baseline 4×1。报告合并前核对三组 input SHA256 相同，
并去除输入文件绝对路径差异后比较物理卡。跨作业节点可能不同，应记录节点并考虑噪声，
需要严格同节点时在用户提交参数中明确同一节点。
该短工况用于验证测量流程，不能据此宣称生产 SASE/ELEGANT 加速。

## 构建和记录

两个性能程序使用既有入口的 `--prepare-only` 编译；该参数不可省略。
完整验证另以 `BUILD_TESTING=ON`、`APRL_BUILD_POSTPROCESSORS=ON` 全量构建，
MPI 启动器指定 `tests/cluster_acceptance/mpi-launch.sh`，进程数参数为 `-n`，pre/postflags 为空。
启动包装必须用于所有 MPI 测试，含不转发 MPIEXEC_PREFLAGS 的旧回归。
绑定探针用同一个 `/usr/bin/mpic++ -std=c++11 -O2 -fopenmp` 编译到
`build/cluster-tools/binding_probe`。具体本次参数、日志和编译标志保存在根目录 `logs/` 与三个构建中。

完成模板同步、编译和提交后，使用 physics Python 执行：

```bash
"$PYTHON_BIN" "$APRL_WORKSPACE/tests/cluster_acceptance/freeze-build.py" "$APRL_WORKSPACE"
```

此命令仅创建 `build-record.json` 和 `manifest.sha256`，拒绝覆盖已有清单。
运行前各脚本验证程序、共享库、卡、检查器、脚本和配置哈希。任何变更都需要重新审阅和新一轮记录。
不要为了通过哈希检查而删除历史清单或伪造 PASS。

## 后续真实输入测试

2026-10-23 只读核对的 ELEGANT 输入为 41,732,775 条 HDF5-v4 固定观察面记录，
2,671,808,608 字节；链接目标仍为旧实验目录中的原文件。它是已有的 40 µm 裁剪数据，
本次未裁剪。现有工况 README 的“约一亿”描述不适用于当前文件。
小工况通过后仍需准备真实输入流式验收及输出空间评估，再设计较大固定 CPU 预算测试。
不直接使用生产脚本默认的 2 节点/256 ranks/48 或 144 小时资源；
SASE 生产脚本目前仍硬编码一线程。MPI 网络、UCX 和系统服务配置保持原状。
