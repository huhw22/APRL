目的与适用范围
  从已验 v5 的无磁自由漂移开始，观察实际推进、沉积和边界操作，不改
  Cowan、Boris、沉积、CPML、参考生成、物理卡或默认程序。仅在独立分支
  和独立 binary 增加默认关闭的只读诊断。

入口
  APRL_GRID_WORK_AUDIT=/已有空目录
    每 rank CSV：物理孔径内原生 Yee 边/分别插值到单元中心的 J·E。
    同时保留 Eold、(Eold+Enew)/2、Enew，Enew 取 MPI 完成后的值，
    左端 Ex/Ey k=0 使用覆盖前缓存。不会通过残差反推出网格功。
  APRL_FIELD_STAGE_AUDIT=/已有空目录
    要求同时启用上项。逐操作记录场能、六面功率、粒子/载流子数、
    全网格 max|J|、右端横向 E 面在移窗前的变化。只读 API 接收 const fields。
  限 moving window、原能量账本启用、无 current filter、无 checkpoint。
  未设环境变量时没有上述场遍历、文件或新增 MPI；本轮没有为关闭开销
  作新的吞吐估计，不能声称严格零成本。开启后增加观察开销，只用于诊断。

源文件
  src/grid_work_audit.h / field_stage_audit.h：所有观察和 CSV。
  moving_window_fields.{h,cpp}：可选 const-field 回调。
  simulation.cpp：按原阶段调用观察器。所有场/粒子赋值和 MPI 顺序保留。

当前程序一步实际执行
  1. 入窗调度；清空本步 J；原 Cowan Faraday 更新 B。
  2. 横向 CPML 修正 B；安装只读粒子采样 halo。
  3. 原 Boris 推物理粒子；冻结粒子/弹道载流子也按原路径沉积 J；迁移。
  4. 合并 MPI 共享 J；原 Ampere 更新 E（本卡无 filter）。
  5. 横向 CPML 修正 E；rank0 覆盖左端 Ex/Ey k=0 为参考 n+1 值。
  6. MPI 完成上下 slab 的 E 接口更新；完成本步源功和六面能流诊断。
  7. 移窗 -1 个纵格：右端裁掉，左端插入参考 E/B，CPML histories 随移；
     清 J，迁移物理粒子和载流子，超过右端的载流子裁掉。
  8. 输出探测面、轨迹及账本。
  z 无 CPML。最右侧原生 Ex/Ey k=nz 不在 Ampere 内区更新循环，物理孔径
  内横向 CPML 也不改它；随后这一面在移窗中被丢弃。它是当前裁剪规则，
  不应仅因这一规则就推断可信内部区被污染。

量的定义（全部 boost 系 J，无 FFT/Hann/扣底）
  S = W_Boris_grid - W_frozen - W_carrier
  G = S - W_mesh
  F = delta_U + Qsix - W_remap - W_inlet - W_boundary + W_mesh
  R_field = G + F
  单元中心/原生网格两个 W_mesh 分别保留，不把更小结果替换原 R。
  Emean 的依据是电场平方能量差的代数恒等式，非按结果挑选时间窗：
    (Enew²-Eold²)/2 = (Enew+Eold)(Enew-Eold)/2。
  这不证明现用单元中心 U、混合 E/B 时间层和 cell-centre E×B 代理
  构成 Cowan+裁剪的严格离散能量恒等式。

验证
  prepare_run.py：复制已闭合 free-drift-r1-new-on 物理卡，只改输出路径；
    冻结 binary/input/卡/脚本/绑定程序。--field-stages 开第二观察器。
  run.sbatch/run.py：1 exclusive128CPU，16MPI×8，无节点 pin，20min 硬限；
    binding → 77required → 观察器解析检验 → 20 无磁耦合探针 → 257步重放；
    检查 H5 complete/finite、逐 ID 覆盖、全部原始场/粒子/轨迹逐位一致，
    v5 所有科学 ledger 字段逐位一致（计时字段单列）。
  deposition_work_probe.cpp：用真实输入第0条动量、synthetic static E，
    独立比较实际 Boris 功、q E_gather·dx、沉积 J·E。无 Maxwell 推进；
    确认显式轨道/沉积不是普遍严格能量共轭，但不外推生产误差百分比。
  grid_work_observer_checks.cpp：均匀场解析功和覆盖前端点缓存检验。
  analyze_work.py：累计与每步绝对差，保留两种空间口径、三个时间配对。
  analyze_stages.py：独立遍历前后真实 U，对原 inlet/remap/Q 逐步核验；
    真正无源要求物理粒子=0、载流子=0、全网格 max|J|=0。
  vacuum_translation_probe.cpp / vacuum.sbatch：原 Cowan+xy CPML+移窗代码，
    无粒子无 J。160²×64，生产 dx/dy/dz/dt，中心 16² 孔径，28步。
    已知右行脉冲在此移动窗每步平移2格，直接核对每步 E/B。
  analyze_vacuum.py：用闭式平移剖面（不迭代 Maxwell）独立计算同一诊断。
    不同宽度改变物理脉冲，不能称固定物理网格收敛。

2026-10-29 已完成结果
  results/energy_residual_origin_20261029/
  19022/node5：77required、20探针及第一只读重放通过。
    61→180（全释放、无退出）：native Emean 净差0.11425%，
    每步绝对差累加0.25657%，最大单步0.62471%；old E 比较差18.804%。
  19023/node5：77required和分阶段重放通过，物理输出和原 v5 账本逐位。
    第239步物理粒子全部退出；245→257 才是严格零粒子/载流子/J。
    此无源段 R 仍变化 -2.150087e-18J；入口覆盖/移窗变化与账本吻合，
    物理孔径 CPML 的直接 U 改变量约舍入级。
  19024：原核真空脉冲正确平移，最大相对偏差2.16214e-13；
    8格脉冲现诊断仍报约 -5.52183% 初始场能差额，闭式模型复现。
    16/32格脉冲分别 -1.54737%/-0.397867%；不是生产误差估计。

结论边界
  证明现诊断会在正确无源传播/裁剪时产生差额，不能把原32%直接认定
  粒子反馈错误、边界制造真实能量或 Cowan 不稳定。原实束工况全程
  -3.02833e-17J 仍未全部定量解释，也未做全面3D物理验收。
  下一步应推导与同一离散时空采样/Cowan/裁剪匹配的能量和能流账，
  保留原诊断对照；不应先按此残差更换推进或边界方法。
  本轮不改默认、不新跑生产、不推送远程；全部原结果保留。
