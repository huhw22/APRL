# Puffin / APRL 同样本源、取场及传播定位

结果目录：`/public/home/wfollowm/APRL_Workspace/results/puffin_aprl_source_propagation_20261113_r1/`。
固定 w64-a、物理种子2026111201、1954743粒子；无新种子或PIC。
19768/19769/19771为闭合数据只读分析；唯一物理控制19770仍是同一粒子、同一Puffin程序的7周期，现已完成。

## 结论

原来31.47pJ对24.22pJ的差异不能继续当作一项单纯的物理增益差。
本轮确认了Puffin磁场类型传递错误、输出取场口径不同，以及APRL当前横向网格下可量化的离轴传播相位误差。
入口和第7周期的粒子横向复源在共同表示下非常接近；这两个截面尚不等于完整沿程4D沉积电流，因此不宣称所有源误差已排除。

| 量 | Puffin | APRL/PIC | 含义 |
| --- | ---: | ---: | --- |
| 旧单段配置，节点频带能量 | 31.47pJ | 24.22pJ | 旧Puffin磁场未匹配，撤销基准资格 |
| 显式匹配磁场后，节点频带能量 | 29.29pJ | 24.22pJ | 原生节点与插值探测像素仍不同 |
| 相同物理像素、各自原生取场 | 22.56pJ | 24.22pJ | Puffin低6.85%；不改原场，不拟合幅度 |
| 7周期平均滑移 | 10.8108fs | 10.7896fs | 明显不同的旧值为11.5044fs |

均为全记录Rect±5%、±560µm，实际相同物理频点。
不同取场行是**不同离散观测量**，不是辐射能量的补偿阶梯；不能把22.56pJ宣布为Puffin的唯一正确能量。
同点取场后的载频复场空间相关度约0.919，说明总量接近仍未等于复场全部一致。

![分环节诊断](/public/home/wfollowm/APRL_Workspace/results/puffin_aprl_source_propagation_20261113_r1/stage_diagnosis.png)

![共同CIC横向复电流源与同像素频带场强](/public/home/wfollowm/APRL_Workspace/results/puffin_aprl_source_propagation_20261113_r1/transverse_source_and_field.png)

上排是载频横向电流源幅值，下排是±5%频带能流密度积分；两排物理量不同。两程序在每排共用色标。

## 1. 修正真实运行路径，而非只检验磁场函数

Puffin `acc_lattice.f90:setupMods` 的无lattice分支分配 `zundtype_arr` 后没有赋值。
随后 `initUndulator` 把 `zUndType_G` 设为该未初始化数组成员，`PopulateUndulatorFromGlobals` 再把它复制到实际RHS使用的段对象。
因此主输入的 `zundType='aprl-ramp77'` 并不保证进入该自定义函数。
此外，`qUndEnds=false` 会在模块入口/出口调用 `matchIn/matchOut`，适用于无端部周期轨道的人工匹配，与已包含平滑入口的外磁场冲突。

上轮的独立磁场函数探针，以及**在段入口处理之前**的原生导入逐位核对，都不能验证这条完整路径。此前将它称作“匹配外磁场”的判断过早；本轮已在历史报告入口标记撤销。

本轮不改数值求解器，只改控制卡：

```
lattFile='matched.latt'
qUndEnds=.true.
```

```
UN 'aprl-ramp77' 7 1.0 0.0 80 0.0 1.0 0.0 0.0
```

`qUndEnds=true` 在此用于关闭人工入口/出口匹配；自定义外磁场函数直接实现APRL原77周期、一周期平滑端部，不采用Puffin默认两周期端部。
7周期是观察截面，不把77周期磁场的物理出口提前。

准入补上两项端到端核对：

- 用同一已编译Puffin库实际执行 `setupMods → initUndulator → PopulateUndulatorFromGlobals → getBFields`，段类型必须为 `aprl-ramp77`；在1/4周期场值为0.14644660940672616 B0，不能落入默认sin模型的1 B0。
- 实际运行保存第1步粒子快照：350µm后平均Δux为−0.000791786，与平滑入口公式−0.000791674相符；所有初始7列仍逐位保持。

闭合证据：`physical-route-control-r1/{route.log,PHYSICAL_ROUTE_VALIDATED.json,COMPLETE.json}`。
使用显式物理步索引选取首次快照，避开Puffin模块末尾的 `matchOut/correctTrans` 人工重心/动量调整。
`puffin_dump()` 选择同一真实步的最早快照，不能把末尾重心修正后的dump当作该截面的真实粒子状态。
已有生产模板也改为显式lattice，防止后续复用错误路径；原冻结卡、原始数据、旧程序未改。

## 2. 横向复电流源

对实验室固定面上的真实到达事件计算

\[
\widetilde{J_x}(x,y,\omega;z)
=\sum_p q_p\frac{u_{xp}}{u_{zp}}
 e^{-i\omega(t_p-z/c)}\delta^2(r_\perp-r_{\perp p}).
\]

这是电流过面时间傅里叶系数，单位C/m²。先让两套粒子通过**同一CIC横向表示**，将物理源变化与各自原生沉积基函数分开。
保留原宏电荷；共同参考时钟来自PIC输入重建中已记录的 `reference_time_s`，没有按输出拟合时间或相位。

修正路径后，第7周期：

- 全频带CIC横向Jx平方范数比Puffin/PIC为约1.0102，逆uz密度代理约1.0006。
- 载频Jx复空间相关度0.9999895；这检查相位图案，而不只检查全束BF幅度。
- 两套束流横向rms约121.01/121.02µm和94.79/94.80µm；源图基本重合。

旧Puffin横向源幅值也接近，但有约−2.74rad共同载频相位差。由真实平均到达时间算出的相位差−2.74048rad与之相符，指向滑移/入口模型，不能用BF幅度相同排除这个问题。
修正后相位差约−0.0812rad，其中已知代码常数 `c_code=299792448` 与物理c的传播时钟差对应−0.08361rad；不是新的共同自场调制证据。这个常数差没有用于改变任何能量或拟合输出。
零发射度轴上外磁解析轨道给10.77615fs滑移，作为独立量级参照，不能充当完整热束的精确预测。

数据：`physical-route-control-r1/analysis/{source_summary.csv,direct_source_modes.csv,source-*.npz}` 与 `phase-analysis/{arrival_moments.csv,source_phase.csv,TIMING.json}`。
只测两个真实粒子截面，未保存/重建每个Maxwell步的完整4D源，因此不能把这些源范数直接乘到输出功率上。

## 3. 原生沉积基函数不同

实际源码：APRL `longitudinal_shape.h::visit` 中Jx沿x为边上的常数基，沿y为线性基，沿z为线性基加三点核；Puffin `getInterps_3D/getSource_3D_kernel` 使用三线性节点源。
在同一已闭合PIC粒子事件上，取无穷小轨迹段极限，单独计算横向基函数的响应：

| 横向谱范围 | APRL NGP-x/CIC-y源功率 / 同粒子CIC-x/CIC-y源功率 |
| --- | ---: |
| k⊥≤10000 m⁻¹ | 0.9706 |
| k⊥≤20000 m⁻¹ | 1.0300 |
| k⊥≤30000 m⁻¹ | 1.0627 |
| k⊥≤40000 m⁻¹ | 1.1017 |

三点核对轴上20格辐射模的独立源功率响应为 `cos^4(π/20)=0.951655`。
这是确定存在的离散源差别，不是错误修复系数。
表中未包含实际有限轨迹段、时间平均、纵向网格相位、别名与反馈；不据它完成全输出能量的百分比归因。
见 `kernel-analysis/{RESULT.json,transverse_source_basis.csv}`。

## 4. 场取样与像素位置

Puffin存的是同一节点格点上的两个横向实场分量（实部Ex、虚部−Ey）；APRL探测面是在Yee分量上插值到单元中心。
本工况两个横向格点相差半格。在同一PIC像素处使用Puffin自己的线性取场，两个横向方向的权重均为1/2：

\[
E_{\rm Puffin}(i+\tfrac12,j+\tfrac12)
=\tfrac14(E_{i,j}+E_{i+1,j}+E_{i,j+1}+E_{i+1,j+1}).
\]

实际坐标核对半格权重误差小于6.2e-15，直接四节点求值而非任意低通。
对应傅里叶响应为 `cos(kx dx/2) cos(ky dy/2)`。
APRL Ex/Jx本身沿x在边中心，探测面只需沿y平均；Ey相应沿x平均。因此不能把节点强度积分与已经插值的像素强度积分视为同一离散量。

修正磁场后，各观察算子产生：

| 同一闭合Puffin场的求值方式 | ±5%能量 / pJ |
| --- | ---: |
| 原生节点求强度后积分 | 29.29490 |
| 仅用傅里叶重建到同一坐标 | 29.05407 |
| 仅施加APRL横向探测平均响应 | 25.59455 |
| Puffin原生双线性取场到PIC像素 | 22.55998 |

对照PIC原探测面为24.21877pJ。最后一行比PIC低6.85%，并不意味着它是一个新的“正确Puffin总能量”。
它定位了**比较结果对离散取场算子的敏感性**；原生场仍原样保留。
这四行也不是独立修正项，不能逐行相减解释物理能量去了哪里。
100µm Gaussian横向场模在仅匹配APRL横向取场响应后，两程序能量8.8163/8.8146pJ接近；更宽横向相干模式仍明显不同，不宣布复场已一致。

## 5. 场传播：已确认离轴Cowan相位差

在固定**实验室频率**下，用实际Cowan九点平滑符号Cj求离散根。
令 `q=kz' dz`、`φ=ω' dt`、`c dt=dz`，正确的两系关系是

\[
\phi=\frac{\omega_{lab}\Delta t}{\Gamma_b}-\beta_bq,
\qquad k_{z,lab}=\beta_b\frac{\omega_{lab}}c+\frac{q}{\Gamma_b\Delta z}.
\]

\[
\sin^2(\phi/2)=\frac{\Delta z^2}{4}
\sum_j\left[\frac{2\sin(k'_j\Delta j/2)}{\Delta j}\right]^2 C_j(k').
\]

没有把实验室小角度直接当作boost坐标中的传播角度。
Puffin对应傍轴 `kz=k−k⊥²/(2k)`。
同一真空模式传播0.196m，相对精确真空色散的相位误差为：

| k⊥沿y / m⁻¹ | APRL Cowan / rad | Puffin傍轴 / rad |
| --- | ---: | ---: |
| 0 | 0 | 0 |
| 10000 | 0.01576 | 1.17e-7 |
| 20000 | 0.17464 | 1.87e-6 |
| 30000 | 0.79831 | 9.49e-6 |
| 40000 | 2.37807 | 3.00e-5 |

为排除仅是推导/读码错误，直接编译现有 `eb_field.cpp/boostframe.cpp`，让真实Cowan/Yee更新体推进六个精确交错平面波。
实际E更新误差最大2.13e-14、cB误差最大3.35e-15，验证了上述离散符号。没有修改APRL生产程序或计算方法。
探针使用的C++与头文件也逐项对应原生产程序构建时的冻结源，见 `VERIFICATION.json`。
证据 `COWAN_ACTUAL_KERNEL_PROBE.{csv,json}`，完整角度/对角方向见 `analysis/vacuum_dispersion.csv`。

正确配置的Puffin整个保留横向场在该频带内，约38.6%的傅里叶能量位于k⊥>20000 m⁻¹，约26.5%在>30000 m⁻¹。
这是整张存储面上的模分解，含离散源/有限域效应，不是所选孔径的远场角分布。
它说明离轴相位差涉及实际存在的模式；不将该比例直接当辐射误差占比。
轴上无色散不保证当前40µm横向网格的三维复场已经足够准确。

## 6. 目前能说到哪一步

已经确认：旧Puffin参考配置不匹配；共同源的两个截面幅值/空间相位接近；原生基函数和输出取场算子不同；APRL离轴真空传播有可重复的相位误差。
不能把剩下6.85%或某一条宽模式差全部算到传播，也不能将源核、插值与色散的几个数字直接相乘“闭合”原30%。
尚未用同一完整4D逐步源交叉驱动两套传播器，也没有网格/步长收敛验收；近场和辐射未被严格分离。
后续若继续应优先量化实际模式上的传播/取场收敛，而不是追加种子稳定性检查。
原31.47pJ图表保留为历史响应并加失效说明；修正后29.29pJ与同点22.56pJ分别保留各自定义，不覆盖成单一“修正功率”。

两张最终PNG为 `stage_diagnosis.png`、`transverse_source_and_field.png`，已目视；本轮不导出PDF。
全程未读取运行中H5、未改变原冻结任务/输入/default/cluster.env、未分析19743/19744、未推送。

## 7. 复核与重现入口

`control.py/control.sbatch` 是唯一修正磁场路由的同样本物理控制；`analyze.py`、`phase.py`、`kernel.py` 只分析已闭合父数据，需受限Slurm allocation。
`report.py RESULTS` 从小型缓存重绘两张图；`verify.py RESULTS` 独立复核冻结身份、缓存积分及实际核探针结果，不提交或推进仿真。
初始只读 `PLAN.json` 的 `new_physical_simulations=0` 只描述19768原计划；后续物理控制另存子目录PLAN，最终 `RESULT.json/JOB_LINKS.json` 明确记录唯一控制19770，不覆盖原计划。
所有原生数据保持；不声称重新逐字节扫描大型原始场H5。
