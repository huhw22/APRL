# 纵向动量负偏差定位（2026-10-30）

目标是定位球形束团相对独立参考的平均纵向动量偏置，暂不开展elegant对照。
不改变Cowan、Boris、沉积、CPML算法、快参考生成或held释放事件。

先复用三组闭合v5账本：用参考世界线与程序共同boost时间超平面求交，
比较全部已释放、尚无移除的同一粒子集合的总动量与动能。相对boost为
Gamma_rel=cosh(acosh(Gamma_beam)-acosh(Gamma_program))，不用固定lab面平均
代替同时刻总量。参考交点满足 t*+beta_rel z*/c=指定常数，逐壳层/极角求解。

诊断程序只在环境变量APRL_PARTICLE_MOMENTUM_PROBE和名义DT显式设置时输出。
每20个场步记录实际采样E/B、源ID、时间/本次dt、位置及推进前后u；每个
MPI/OpenMP线程单独160字节定长记录、直接write，无MPI/HDF5或缓存未落盘。
原状态产生表达式保持。独立Python Boris回放必须复现实际后状态。
按负电子荷拆分纵向电冲量q dt Ez/(m c)和剩余磁旋转冲量。
参考场在实际静止系时间、理论粒子位置取值，并非对扰动PIC分布重解Poisson；
这一近似必须与实测位置误差共同解释，不能直接称所有场差都为程序bug。

小对照均使用原uniform输入/权重、16MPI×8、相同dt/dx/dz，保留入口粒子面，
省略尚未到达的下游诊断面，缩短物理stop到500步并重定向输出。通过reference-center-z正常闭合，
不是把duration耗尽或未完成输出当成功。保留原reference完整网格和源几何。

- base：160²，lower guard64，活动Nz1332，与原结果前500步科学账本逐位核对。
- inlet_guard_256：只将lower guard64→256、活动Nz增加192格，右侧物理位置保持，
  检查左端供场/裁剪影响的因果延迟。原源参考域不改，缩短时段后查询仍有覆盖。
- transverse_240：横向160²→240²、dx/dy与CPML厚度保持，移动横向边界位置。
  这也会改变快参考的有限横向求解域，因此须检查早期差异，不把它单独当纯CPML开关。

预检为77required、100步探针off/on科学输出逐位，并确认确实写出探针。
前两次卡保留了未到达的下游粒子面，stop被全部element约束初始化拒绝，均未推进；
40步本身还位于入口上游。77required已独立通过，最终卡省略下游诊断面，
使用合法的100步停点，并静态核对stop约束；仅复用同binary哈希的回归记录。
每作业1独占128CPU节点、不指定节点、25分钟Slurm硬限，单case1200秒，
原600/60秒MPI-free watchdog和既有self,vader,tcp保持。不改default。
独立目录、冻结脚本/输入/binary哈希；计算失败或数据不闭合则拒绝结论。

```bash
python3 verification/performance/cases/longitudinal_momentum_origin/prepare.py \
  --workspace /public/home/wfollowm/APRL_Workspace --build /absolute/new/build \
  --root /absolute/new/result
sbatch /absolute/new/result/scripts/run.sbatch /absolute/new/result preflight
# 预检明确通过后，分别提交base/inlet_guard_256/transverse_240。
```

后续已实施同参考物理区域的厚CPML控制（prepare_matched_domain.py）：
横向192²/CPML26格，对照160²/10格，二者物理区均±350um，标量内部节点均139。
原rho/phi逐数核查相对L2为2.03e−14/3.45e−16。500步lab能量误差−98.634keV，
原−98.630keV；移远接口到±550um则−18.651keV。单步原场核探针在原接口首个
electric curl生成+114.903kV/m，电场CPML修正不改变该接口增量；移远接口峰46.517kV/m。
这定位了有限Dirichlet参考场截断与随后开放边界推进之间的初态匹配问题。

最终力后处理必须扣除geometry.box_translation_m再从box坐标变换。早期冻结
r3分析漏减该空间原点，仅影响参考场取样时间；最终review/补齐并核对直接lab
变换，还同时在实际PIC坐标评估连续参考密度。共同时间总动量、实际采样E/B、
Boris重放和首步场核证据不受该后处理修订影响。未改变原核心物理方法。
