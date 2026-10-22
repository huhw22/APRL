# 独立 10 µm、50 A 自由漂移启动实验

这里的输入由 APRL 自身生成，**不读取 Elegant 导出**。实验用于比较完整三维冻结漂移与“横向位置保持至启动面”两种启动处理；不放波荡器。生成器在固定实验室输入面产生 HDF5 v4 粒子事件：纵向均匀分布，横向和动量按照旧 `beam.ele` 的 Twiss 参数独立抽样。默认全长 10 µm、电流 50 A、中心动量 600 MeV/c、`emit_x=emit_y=10⁻⁹ m`、`beta_x=4.527431766 m`、`beta_y=4.77124831 m`、`sigma_dp=1.3×10⁻⁴`、随机种子 20261019。约 1.668 pC 由 10,409,720 个电子表示。

## 本地小规模测试

从仓库根目录执行。构建需要与运行时相容的 MPI、HDF5 和 yaml-cpp；若仅单进程且只有串行 HDF5，可在 CMake 配置时显式使用 `-DAPRL_REQUIRE_PARALLEL_HDF5=OFF`。

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j4 --target aprl uniform_lab_plane_bunch_generator energy_ledger_report
mkdir -p work/uniform10um
./build/uniform_lab_plane_bunch_generator --plan --macro-weight-electrons 10
./build/uniform_lab_plane_bunch_generator \
  --output work/uniform10um/particles-10um-50A-w10.h5 \
  --macro-weight-electrons 10 --seed 20261019
./build/uniform_lab_plane_bunch_generator \
  --verify work/uniform10um/particles-10um-50A-w10.h5
```

权重 10 的文件约 1,040,972 行，但 `beam.input.electrons` 仍为 10,409,720；每个宏粒子的物理电荷是 10 个电子。若改动生成器的长度、电流或权重，要从生成器输出/文件属性读取新的 `represented_electrons`，并作为下面的 `--electrons`，不能机械沿用本例数值。

两张历史本机 γ=2 卡 [`pilot-g2-held-w10.yaml`](pilot-g2-held-w10.yaml) 和 [`pilot-g2-ballistic-w10.yaml`](pilot-g2-ballistic-w10.yaml) 曾用于修复前的较长漂移试运行，其固定输出目录已经占用；**不要直接重新运行**。它们说明了带初场、自场反作用、每步能量账、单个场面且不写大粒子探测面的设置。修复后的默认测试用 `g2-short`：束团中心从实验室 z=10 µm 到 z=30 µm，即目标约 20 µm 自由漂移；下面命令会生成全新的卡和输出目录。旧 `g2` 配置仍可作为到 z=100 µm 的长漂移参考。

```bash
python3 experiments/uniform10um/prepare.py \
  --input work/uniform10um/particles-10um-50A-w10.h5 \
  --output-dir work/uniform10um/g2-suite-w10 \
  --profile g2-short --electrons 10409720
mkdir -p work/uniform10um/g2-suite-w10/runs/g2-short-held-field
mpiexec -n 4 ./build/aprl \
  work/uniform10um/g2-suite-w10/cards/g2-short-held-field.yaml \
  > work/uniform10um/g2-suite-w10/runs/g2-short-held-field/run.log 2>&1
```

`prepare.py` 默认写六组卡：`held`/`ballistic` ×（初场开且反作用开、初场关且反作用开、初场开但反作用关）。先运行上面两个“初场开且反作用开”的配对组；需要分离零场启动和轨道反作用时再逐组运行控制组。不要同时盲目提交六个大任务。若只需部分组，使用 `--cases held-field,ballistic-field`。γ=100 权重 10 组可用 `--profile g100` 生成，但在其 `dz=2 µm` 下每场步约对应实验室 397 µm 粒子行程，**不能用于判定 20 µm 启动瞬态**，也不是 γ≈929 的收敛证明。改变输入文件后生成新卡，不要改旧卡对应的运行记录。只有明确要保存每个过面粒子时才加 `--particle-plane`：权重 1 时单面可超过 1 GB/组。

启动原点曾被发现重复加上 `initial_center_z`；修复并重编译之前的数值**一律作废**。修复后权重 1 的 `g2-short` 主配对已在新目录重跑：两组均正常结束，10,409,720 个宏粒子全部过启动面。横向保持/完整三维的场能初值分别为 `1.061901394/1.061898464×10⁻⁹ J`，峰值为 `1.061991436/1.061988506×10⁻⁹ J`；最大单步正增量分别为 `9.004209/9.004167×10⁻¹⁴ J`。在相同计算系时间比较，最大场能差为较大峰值的 `2.759×10⁻⁶`，没有检测到病态增长。权重 10 的重复试验给出相近的差异比例 `2.735×10⁻⁶`。这些只是 **γ=2 的粗步长 pilot**，不是 γ≈929 的验证：`dz=1 µm` 时每个场步对应实验室粒子约 3.732 µm 行程，20 µm 仅约 5.36 次场更新，尚需较细 `dz` 收敛扫描。两组结束时间不同，不能直接把末态场能作相同时间差。冻结/边界模式的 `closure_valid=0`，权重 1 两组能量交换归一化残差警告约 6.50%/6.67%，初始 Poisson 边界约 3.92 RMS；后续仍需初场边界、网格、CPML、孔径和权重扫描。

权重 10 的四组控制卡也已正常跑完，全部粒子过启动面。关闭初场但保留粒子反馈时，两模式从零场增长到约 `3.70/3.94×10⁻¹¹ J`，最大单步正增量约 `2.37/2.41×10⁻¹² J`，明显大于有初场组的约 `9.0×10⁻¹⁴ J`；零场启动故意违反初始 Gauss 约束，**不能把这段增长直接解释成自激不稳定**。保留初场而关掉粒子反作用时，账本报告精度下的场能历程与有反作用组几乎相同，说明这段 γ=2 短漂移内反馈影响尚未分辨出来，而不是证明反馈永远可忽略。零场组原始交换残差警告约 18.5%/19.8%，`closure_valid=0`；控制组也不是严格的能量守恒检验。完整原始输出在本机 `work/uniform10um/corrected-w10/`，不随源码包传输。

还做了一组权重 10、γ=2 的细网格对照：`g2-short-fine` 将 `dz` 从 1 µm 降至 0.5 µm、`nz` 从 256 加至 512、z-CPML 从 12 加至 24 格，保持盒长与 CPML 物理厚度不变。两模式均正常结束且全部粒子激活；横向保持/完整三维的初始场能为 `1.067558505/1.067554745×10⁻⁹ J`，最大单步正增量为 `2.228685/2.228673×10⁻¹⁴ J`，约为粗网格的四分之一。相同时间模式差的最大值占场能峰值 `3.52×10⁻⁶`，粗网格为 `2.73×10⁻⁶`；初始场能绝对值也改变约 0.5%。**这还不是收敛证据**，只是表明粗步长的单步尖峰对步长敏感。若要重现细网格卡，可运行：

```bash
python3 experiments/uniform10um/prepare.py \
  --input work/uniform10um/particles-10um-50A-w10.h5 \
  --output-dir work/uniform10um/g2-fine-repeat \
  --profile g2-short-fine --electrons 10409720 \
  --cases held-field,ballistic-field
```

本机参考细网格结果在 `work/uniform10um/g2-fine-w10/analysis.json`；该目录使用实验时手工改出的 `g2-short` 卡名，新 `g2-short-fine` profile 复现的是其物理配置，并在输出名中显式区分网格。

用已构建的 `energy_ledger_report` 读取完整场能历程，无需 Python `h5py`：

```bash
python3 experiments/uniform10um/analyze.py \
  work/uniform10um/g2-suite-w10/runs/g2-short-held-field \
  work/uniform10um/g2-suite-w10/runs/g2-short-ballistic-field \
  --json-out work/uniform10um/pilot-summary.json \
  --csv-out work/uniform10um/pilot-energy.csv
```

分析器逐条读取已提交的账本记录，给出初值、峰值、末值、最大单步正增量，并在相同计算系时间插值比较两组。CSV 可自行绘图；它不需要载入完整粒子 HDF5。`boundary_source_abs_work_J` 代表边界源活动量级，**不是**外部近场缺失比例；`closure_valid=0` 的原始残差不是严格整束能量守恒证明。

## γ≈929、1 m 无波荡器漂移候选卡

`g929-1m` 是针对启动方法的**首轮大束测试**，不是已经完成的物理收敛结论。束团初始中心在实验室 z=10 µm，末端场面和可选粒子面在 z=1,000,010 µm，故中心到探测面恰好 1 m。卡使用 γ=929.130549、`dz′=5 µm`、自动 `observation-bounded` 纵向盒、初场 4σ 留白和两端各 384 个物理保护格；权重 10 的轻量几何预检得到 6,144 个纵向网格单元，正式运行仍须看启动日志的网格和内存估计。`mesh.duration=50 ps` 只是硬上限。停止条件是 `after-last-element`，会等待所有**仍有效**粒子通过末端面。预计整束通过约需 2,040 个场步（约 34 ps 计算系），实际以启动时 `[duration-estimate]` 和停止日志为准。场面以 1 fs 实验室时间间隔记录原始 E/B；它不是已经减去粒子近场的辐射场。

```bash
python3 experiments/uniform10um/prepare.py \
  --input work/uniform10um/particles-10um-50A-w10.h5 \
  --output-dir work/uniform10um/g929-1m-w10 \
  --profile g929-1m --electrons 10409720 \
  --cases held-field,ballistic-field --particle-plane
```

`prepare.py` 只写卡，不运行 APRL。全粒子 1 m 截面必须用 `--particle-plane`：权重 1 的单组记录约 1.2 GB，可按粒子 ID 与输入惯性轨道配对。整束末态能量账可监测自场量级及突增，**不能当作全体粒子恰好走了 1 m 的截面**：尾粒子到面时，头粒子已继续漂移。正常退出还要核对 `analysis.json` 的 `all_input_particles_past_last_element=true`、`activated_macroparticles` 等于输入行数、`still_frozen_macroparticles=0`、无逃逸粒子，且粒子面 HDF5 完整并有预期记录数；单凭 APRL 进程退出码不够。若一些粒子横向逃逸或纵向自动盒报错，仍可保留结果诊断，但不可称为完整 1 m 截面。

仿真结束后可选用 `analyze_particle_plane.py` 将过面 x/y 和方向斜率与各自输入粒子的无场直线外推配对；它按 `source_id` 分块读取，给出均值、RMS 与抽样绝对误差分位数，不把约千万粒子同时载入内存。此后处理需要 Python `h5py`、NumPy 和至少约 350 MB 临时磁盘，**不参与 Slurm 仿真**。例如：

```bash
python3 experiments/uniform10um/analyze_particle_plane.py \
  --input work/uniform10um/particles-10um-50A-w10.h5 \
  --plane work/uniform10um/g929-1m-w10/runs/g929-1m-held-field/detectors/after-1m.h5 \
  --plane work/uniform10um/g929-1m-w10/runs/g929-1m-ballistic-field/detectors/after-1m.h5 \
  --json-out work/uniform10um/g929-1m-w10/particle-plane-comparison.json
```

在 γ≈929、`dz′=5 µm` 下，单粒子每场步约走 5.71 mm 实验室距离，1 m 约 175 步；整束因固定实验室输入事件的相对论同步跨度，需约 2,040 步。`dz′/γ≈5.38 nm` 是前向光波的实验室等效采样刻度，不是粒子每步行程，且不满足旧 4 m 卡约 1.35 nm 的光波采样要求。这组卡不用于辐射频谱结论。自动盒保证设计的粒子世界管和前向因果路径，不证明末态全空间场或近场边界已收敛；初始 Poisson、横向开边界注入、CPML 和宏粒子权重仍须扫描。轻量几何预检发现 ballistic 模式初始约 5.96% 粒子在整个横向网格外，held 模式约 0.03%；两组初始自场实际纳入的源粒子数因此不同，场差不能只解释为启动面后的轨道动力学差。关闭初场的控制组会跳过 4σ 初场留白，自动网格因而可能改变，不可在未对齐网格的情况下将差异全归于初场。现有冻结/开边界能量账 `closure_valid=0`，其残差只作数量级诊断。

## held 启动的 1 m 网格/边界敏感性测试

这组测试**只比较 held-field 与 held-field**：相同的现有 `particles.h5`、初始自场和粒子反馈，观察面仍在实验室中心后 1 m。它检验保持横向位置这一启动近似对出口总场及粒子末态的数值敏感性，不以 ballistic 启动为“物理真值”。默认新跑三个变体（每个都记录完整的 1 m 粒子面）：

| 变体 | 相对参考卡的唯一目的 | 横向物理半宽 |
|---|---|---:|
| `wide` | 104² → 128²，CPML 仍 10 格，检查横向口径 | ±270 µm |
| `thick-cpml` | 124²、CPML 20 格，物理区宽度不变，只加厚吸收层 | ±210 µm |
| `fine-z` | `dz′` 5 → 2.5 µm，纵向 CPML/guard/alignment 格数同步加倍 | ±210 µm |

还可用 `--variants base,wide,thick-cpml,long-guard,fine-z` 加跑两组可选控制：`base` 与参考卡相同，用于测量重新编译或 MPI 布局变化造成的重复运行差；`long-guard` 将纵向两端物理保护区从各 384 格加到各 768 格，用于检查纵向边界对末场的影响。它们默认不跑，以免重复花费。`fine-z` 和加宽网格的内存消耗会明显增加，也可先用 `--variants wide` 单独提交试跑，其余变体各用新的 `SUITE` 目录提交。

下例中的 `OLD` 指已有**正常完成的 held-field** 作业根目录，其中有 `particles.h5` 和 `runs/g929-1m-held-field/`；`SUITE` 必须是不存在的新目录。先用 `--plan` 核查绝对路径和资源，再去掉它提交。`NODES=2`、每节点 64 rank 仅示范你之前使用的资源量；新变体可能增加网格内存，先结合启动日志的 `[resource]` 估算确认节点容量。脚本只读取现有 H5，不再生成或复制粒子，不依赖 Slurm spool 工作目录。

```bash
SRC=/public/home/wfollowm/aprl_dev
OLD=/absolute/path/to/existing/held-job
SUITE=/absolute/path/to/new/held-convergence-suite
NODES=2 TASKS_PER_NODE=64 MPI_RANKS=128 APRL_SRC_DIR="$SRC" \
  "$SRC/experiments/uniform10um/submit_held_convergence.sh" \
  --reference "$OLD/runs/g929-1m-held-field" \
  --input "$OLD/particles.h5" --output-dir "$SUITE" \
  --electrons 10409720 --plan
# 确认路径和资源后，原命令去掉最后的 --plan 才会提交。
```

提交后看 `$SUITE/status.tsv`、`job-status.txt` 和每组 `runs/g929-1m-held-field/run.log`；`complete` 必须同时有正常停止、完整能量账和两个探测器文件存在。文件内部的完整标记会在后处理时再核对。**作业不自动跑 Python 后处理**，以免计算节点缺 `h5py` 时把成功的仿真误报失败。全部完成后，在有 NumPy/h5py 的环境中比较同一固定中心 ±160 µm 区域、共同实验室时间窗的 E/B 和 Poynting：

```bash
python3 "$SRC/experiments/uniform10um/analyze_held_convergence.py" \
  --reference "held=$OLD/runs/g929-1m-held-field" \
  --variant "wide=$SUITE/wide/runs/g929-1m-held-field" \
  --variant "thick=$SUITE/thick-cpml/runs/g929-1m-held-field" \
  --variant "fine=$SUITE/fine-z/runs/g929-1m-held-field" \
  --roi-half-width-um 160 \
  --json-out "$SUITE/held-field-comparison.json" \
  --csv-out "$SUITE/held-field-traces.csv"
```

报告中的 `relative_E/B_spacetime_L2` 和相同时间窗的出口场能流用于判断收敛；`tail_unresolved=true` 表示最后样本尚未衰减（默认阈值为峰值 5%），`after-last-element` 的时间窗可能截掉后续场，累计通量不能据此定论。场面是**未减去粒子近场的总 Maxwell 场**，Poynting 不是纯辐射功率。`U_field` 属于 boosted 系整个物理盒，横向盒大小改变后，绝对值不能直接拿来当共同空间区域的收敛判据；`closure_valid=0` 也仍不是严格能量闭合。`fine-z` 的 E/B 时间交错误差不同，这个测试不能代替目标光子能段的辐射频谱收敛。逐粒子末态差另用 `compare_particle_planes.py` 按 `source_id` 比较；两组相对于无场直线轨迹的 ~10 nm 偏差，不等同于两组之间的差值：

```bash
python3 "$SRC/experiments/uniform10um/compare_particle_planes.py" \
  --reference "$OLD/runs/g929-1m-held-field/detectors/after-1m.h5" \
  --variant "$SUITE/wide/runs/g929-1m-held-field/detectors/after-1m.h5" \
  --scratch-dir "$SUITE" \
  --json-out "$SUITE/particle-wide-vs-held.json"
```

把 `wide` 换成 `thick-cpml`、`fine-z` 再各跑一次。`dxy_um` 是两组逐粒子横向位移差的模；`dr_um` 是两组到轴线距离的差，不是同一量。这个配对器只比较两面都出现的首个下游过面事件，两个文件都缺失的粒子无法由它发现；还应使用上面的 `analyze_particle_plane.py` 核对每组相对于原始 `particles.h5` 的完整性。临时磁盘约需 34 字节/最大 `source_id`（一千万粒子约 350 MB），不会把全部记录载入 RAM。

## γ≈929 几何预检与生产限制

固定实验室面输入在高 boost 的 `t_box=0` 伸展约 11.4 mm，而从第一到最后一个粒子过启动面的计算系时间差约 31 ps。因此不能把一个只覆盖微米自由漂移的盒子或几十飞秒的运行当作完整整束启动。曾经拟议的 5.12 mm 纵向网格不足以容纳初始分布。即便增大到预检的 `3072×5 µm = 15.36 mm`，也只验证了初始几何，没有证明整段演化期间粒子、近场和观测面都远离纵向 CPML。

更根本的是时间分辨率：γ≈929、此束能量下，`dz=5 µm` 的 Cowan 步长为 16.6782 fs（计算系），粒子每步约前进 5.7136 mm（实验室系）；目标 20 µm 的自由漂移只占约 0.00350 个 Maxwell 步。这与旧 4 m 设计中“1.35 nm 实验室等效纵向分辨率”不矛盾：旧卡使用约 1.254 µm 的 boost 系网格，1.35 nm 是前向光波的实验室等效采样尺度；本预检 5 µm 网格按同一口径约为 5.38 nm。粒子沿其世界线跨一个场步的实验室位移则是另一个量，`Δz_lab = γ_boost (β_boost + β′_z) dz_box`。即使粒子过探测面的记录可在一步两端之间插值，场源/PIC 仍按 Maxwell 步更新，不能据此称 20 µm 内的自场变化已被解析。粒子子步或扩大 `nz` 都不能代替场的时间分辨率。加上 `ballistic_outside_transverse` 与 `particle_retirement` 在现有配置中互斥，现阶段**没有经验证的 γ929 整束短漂移带场服务器卡**，也不会在本包里给出貌似可运行的生产卡。高 boost 后续必须另行设计可解析短漂移的坐标/时间或诊断方案、完整整束停止条件、纵向因果盒和收敛扫描，不能直接把以下几何卡当成仿真卡。

以下只做准确的初始横向孔径计数，**初场关闭、无反作用、无探测器、无能量账，不能检验场瞬态**。运行器在 APRL 完成几何及初始化日志后发送 SIGINT；最多一个 Maxwell 步，正常受控退出应是状态码 0，日志出现 `Direct E/B run stopped after a complete field step`。不要直接运行这些 16 fs 卡，否则会因未达到物理停止位置而触发预期的 `mesh.duration` 耗尽错误。固定例卡仅适用于本机路径；服务器用生成器重写绝对路径。

```bash
python3 experiments/uniform10um/geometry_preflight.py \
  --input work/uniform10um/particles-10um-50A-w10.h5 \
  --output-dir work/uniform10um/geometry-g929-w10 \
  --electrons 10409720 --run
```

修复启动原点重复加 `initial_center_z` 之后，用完整权重 1 输入（10,409,720 宏粒子）作过一次几何预检：横向保持模式 `t_box=0` 物理区外 43,893 个（0.421654%）、全网格外 2,855 个（0.0274263%）；完整三维模式分别为 1,322,041 个（12.7001%）、623,284 个（5.98752%）。原先修复前的高比例已作废。这里的“物理区外”包括横向 CPML 和全网格外；初始 Poisson 可选择一次性考虑仍在 CPML 内的粒子，但全网格外粒子无法进入初场源。这些数值只适用于本次随机输入、±210 µm 物理区半宽和 ±260 µm 全网格半宽，**是初始几何而非带场稳定性验证**。生产配置或源码变化后要重做预检。横向保持本身是非物理近似，启动面的横向电流跃变可能诱发场，必须再看带场对照。

## 服务器移植

若服务器不能直接 `git pull` 当前未提交的开发版本，可只传源码包，不传约 667 MB 的权重 1 HDF5：

```bash
experiments/uniform10um/package_server.sh
# 将输出的 tar.gz 传到服务器，在专用目录解压
tar -xzf aprl-uniform10um-source-YYYYMMDD-HHMMSS.tar.gz
cd aprl-uniform10um
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j4 --target aprl uniform_lab_plane_bunch_generator energy_ledger_report
mkdir -p work/uniform10um
./build/uniform_lab_plane_bunch_generator \
  --output work/uniform10um/particles-10um-50A-w1.h5 \
  --macro-weight-electrons 1 --seed 20261019
./build/uniform_lab_plane_bunch_generator \
  --verify work/uniform10um/particles-10um-50A-w1.h5
```

权重 1 是 10,409,720 行，文件约 667 MB；服务器上重新生成可避免跨机器搬运大输入。源码包无二进制和 HDF5，服务器必须用其自己的 MPI/HDF5 环境编译。生成器、运行卡、分析器分别使用 `--plan`、`--verify`、独立输出目录，便于给每次运行留可核查的记录。结果归档时请填写 [`REPORT_TEMPLATE.md`](REPORT_TEMPLATE.md)。

服务器上权重 1 的 γ=2、约 20 µm 实验室自由漂移配对提交可参考 [`submit_g2_pair.sbatch`](submit_g2_pair.sbatch)：**从仓库根目录执行** `sbatch experiments/uniform10um/submit_g2_pair.sbatch`。脚本进入 `SLURM_SUBMIT_DIR`（避免 Slurm spool 目录定位错误），默认申请一节点 8 rank，用 `mpiexec -n 8` 顺序跑两个案例，不假设 `srun --mpi` 与已编译 MPI 相容；可在提交环境设 `APRL_MPI_RANKS` 为不超过分配数的整数。先加载与构建 APRL 时完全一致的 MPI/并行 HDF5 模块，再提交。Slurm 的 `slurm-<jobid>.out/.err` 是调度/脚本日志，各组 APRL 日志在各自 `runs/g2-short-*/run.log`，场面和能量账也分别在相应目录。每次以作业 ID 建新目录；若目录已存在就退出，绝不覆盖旧结果。这是模板，不会自行提交任务，也不代表 γ=2 可外推到 γ≈929。
