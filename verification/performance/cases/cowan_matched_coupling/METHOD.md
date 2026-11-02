# 为什么修改粒子取场

以下均在 boost 系。令原 Yee 差分为 $D_i$，Cowan 差分为 $D_i^*=S_iD_i$，$S_i$ 为现有 Cowan 系数确定的对称九点算子。$C,C^*$ 分别表示两套 curl。在周期网格上它们及 $S_i$ 对易；有限物理边界须另记交换，不能直接援用周期恒等式。

定义

$$M_E=\operatorname{diag}(S_yS_z,S_zS_x,S_xS_y),\qquad M_B=\operatorname{diag}(S_x,S_y,S_z).$$

这些算子由当前分辨率和 Cowan 系数唯一确定，没有按结果拟合的常数。逐分量可验证

$$M_E C=(C^*)^\dagger M_B.$$

例如 $(x,y)$ 项两侧都是 $-i\kappa_z s_ys_z$。这里 $s_i$ 是 $S_i$ 的 Fourier 符号。由此，原推进

$$B^{n+1/2}=B^{n-1/2}-hC^*E^n,\quad E^{n+1}=E^n+c^2hCB^{n+1/2}-hJ^{n+1/2}/\epsilon_0$$

满足配套离散能量

$$H_C^n=\frac{\epsilon_0}{2}\langle E^n,M_EE^n\rangle+\frac{1}{2\mu_0}\langle B^{n-1/2},M_BB^{n+1/2}\rangle.$$

电项作平方差，磁项展开相邻两步双线性差，并用上述伴随关系消去 curl 项，得到

$$H_C^{n+1}-H_C^n=-h\langle J^{n+1/2},M_E(E^{n+1}+E^n)/2\rangle.$$

原电流由整条轨道弦沉积。对同一弦取沉积的伴随取场算子 $G_J$，令

$$E_{p,\mathrm{self}}=G_JM_E(E^{n+1}+E^n)/2.$$

原平均速度粒子推进满足 $\Delta K_p=q_p\Delta\boldsymbol{x}_p\cdot E_p$；磁力作功为零。因此闭合周期盒中，粒子电功与上述场能损失相同，只有原非线性迭代容差、舍入及其他独立外功。代码把 $S_jS_k$ 合并为最多 45 个局部权重，避免嵌套读取 81 项；不解新增全空间方程。

## 纵向磁冲量

设 $A_z$ 是相邻 z 层的平均。与 Cowan 配套的参考场动量为

$$P_{z,C}=\epsilon_0\langle A_zE_x,S_y\bar B_y\rangle-\epsilon_0\langle A_zE_y,S_x\bar B_x\rangle,\quad \bar B=(B^{n-1/2}+B^{n+1/2})/2.$$

固定一次 Ampere 源更新时，其电流动量交换为 $h\langle J_x,A_zS_yB_y\rangle-h\langle J_y,A_zS_xB_x\rangle$。用相同轨道弦沉积的伴随分别取 $B_{y,\mathrm{eff}},B_{x,\mathrm{eff}}$，即可与粒子磁冲量

$$\Delta p_{z,p}^{B}=q_p(\Delta x_p B_{y,\mathrm{eff}}-\Delta y_p B_{x,\mathrm{eff}})$$

配对。整个 $\boldsymbol v\times\boldsymbol B$ 保持；$B_z$ 仍用原取样。此阶段目标是 boost 变换需要的纵向动量，并未证明另外两个动量分量严格守恒。外加波荡器场的取样保持原方法。

## 为什么仍须做物理测试

**电功配对和磁冲量配对，不自动保证电力冲量配对。** 场动量关系的电项涉及

$$F_{E,z}^{n}=\langle\rho^n,A_z(M_EE^n)_z\rangle,\qquad\rho^n=\epsilon_0\nabla\cdot E^n.$$

但粒子电力冲量是 $\sum_p q_ph_p(E_{p,\mathrm{self}})_z$。整弦电流的伴随取样与节点电荷的伴随取样不同；有限网格和粒子采样下，两者不必相等。`periodic.cpp` 独立测量这项差，不把它吸收到场端通量或外功，也不通过修改粒子末态强行归零。Cowan 的磁约束是 $\nabla^*\cdot B$，测试分解保留对应项。

第一轮粗加载反例正显示：配套 boost 能量接近闭合时，纵向电力冲量仍能显著漂移。更密加载的轴向例又会因对称性使 $E_z$ 极小，因此增设有 $E_z$ 的斜向波。粒子采样加密后误差变小不能代替固定物理问题的空间/时间收敛；也不能说明实际生产束团已满足辐射精度要求。

lab 组合差额以

$$R_{\rm lab}=\Gamma(R_{K'}+\beta cR_{P'_z})$$

报告，并与实际粒子能量变化比较，而不除以巨大的初始总动能。$H_C$ 和 $P_{z,C}$ 是配套离散泛函；原生 E/B 的物理能量及原输出保留。配套泛函闭合只是一项必要检查，不是原生场、Lorentz 协变性或辐射能量已通过验收的证明。

## 保留原未加权诊断时的功差

原读取器曾隐含假设粒子网格功等于 raw $J\cdot\bar E$。本候选不再满足这个假设。因此另列

$$N_{\rm gather}=W_{\rm particle,grid}-W_{\rm boundary}-W_{\rm frozen}-W_{\rm carrier}-W_{J,\rm applied,mid}.$$

各项来自实际粒子电功、已有虚拟源功和真实场更新采用的电流中点功，均独立测量，不由总差额倒推。原能量恒等式的数值通道须加上此项；它是取场与场能定义的差，不能称为真实外功或辐射。未扣除它的原始差额继续显示。修复读取器后保持原 `1e-6 eV/e` 检查门限，并保留省略此项会失败的负控值。这项读取器修复仍不等于完整移动边界的 $H_C$ 账本已完成。

作为方法背景，[Chen 等人的原始论文](https://arxiv.org/abs/1101.3701)也明确区分能量/电荷守恒与动量并非精确守恒。本文 Cowan 权重由本程序算子独立推导；不把该论文当作本候选已被验证的证据。
