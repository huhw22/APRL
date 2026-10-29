# 球形束团空间电荷：首轮全程序物理对照

本轮使用已验收 `build/discrete-energy-ledger-20261029/aprl`，不改计算核心、
快参考供场或默认入口。沿用生产 boost γ=929.1305492858543、原 dz/dt、
dx=dy=5 µm、横向 CPML 和 lower/upper 64/128 格底座。

束团平均 lab γ约1174，保证当前移动窗要求的 held 纵向速度为非负。
三个独立小工况：冷均匀球、相干向内速度斜率的均匀球、冷球形高斯。
各用32768个带权安静采样点表示1e8电子；初始每轴rms约44.72µm。
参考以球对称 Maxwell 场的径向相对论动量方程独立积分，监测壳层顺序，
并由世界线与固定lab面的交点生成输入位置、到达时间、动量。
不把不同粒子的释放事件重新定义成同一个冷初态。
首轮观测到2.2t0，覆盖已知会聚转折及开始回弹；未改束团电荷或尺寸。

原能量账本与严格离散审计均开启。无外磁、无吸收粒子层、无快照。
只输出粒子面、稀疏轨迹和账本；本轮重点是粒子动力学与能量交换。
主比较逐ID重建束团静止系事件；能量变化相对于同一ID的入口事件。
固定lab面并非共同静止系时刻，理论按相同事件采样。

```bash
python3 verification/performance/cases/space_charge_physics/prepare.py \
  --workspace /public/home/wfollowm/APRL_Workspace --root /absolute/new/result
sbatch --output=/absolute/new/result/slurm-%j.out \
  /absolute/new/result/scripts/run.sbatch /absolute/new/result
```

单个独占128CPU节点，实际16MPI×8，先验证128不同物理核绑定，再顺序执行。
每工况1200秒硬限，作业65分钟上限、不重排、不指定节点名。
保留600/60秒MPI-free收尾watchdog及现有作业内self,vader,tcp。
生命周期、完整标记、入口映射、有限值和完整粒子覆盖首先核对。
物理偏差直接报告；运行成功不自动宣称物理验收，也不自动放宽容差。

所有脚本冻结复制到结果目录。`model.py`、`prepare.py`、`analyze.py`
各负责理论、生成卡与输入、闭合后对照分析。`launcher.py`原样复用既有
限时启动和绑定检查函数。首轮无网格/粒子收敛扫描，后续按结果决定。

当前模板的账本图与`energy_budget.csv`默认以`D=max|ΔK_particle|`为分母显示
带符号百分数，原J数据另存`energy_budget_absolute.csv`，分母说明在
`normalization.json`。该显示修改不改变账本计算或物理判据；已启动作业的
冻结脚本保持，使用独立只读后处理补绘。
