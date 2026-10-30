# 物理对照的公式、推导与取样方法

日期：2026-10-30。本文是[阶段性结论记录](README.md)的理论说明，解释图中参照值怎样得到。以下推导不调用APRL的场推进、沉积或粒子推进，已有仿真结果不因补写本文而改动。

需要区分三个层次：Maxwell方程与球对称给出的场关系、用于理解尺度的非相对论解析解、实际生成对照数据的相对论球壳世界线。最后一层再与固定lab探测面求交，才能生成表中的理论均值和能散。

| 对照 | 理论对象 | 怎样得到图中的参照 |
|---|---|---|
| 真空脉冲 | Maxwell一维行波解 | 在原生场样本的位置与时间层计算平移后的场形 |
| 冷均匀球 | 均匀初密度、径向冷初态 | 相对论球壳轨迹；自相似解析式用于解释膨胀与能量尺度 |
| 向内聚焦均匀球 | 同一球形密度、已知线性内向初速度 | 改变初始径向动量，仍解相同球壳方程 |
| 冷球形高斯 | 高斯初密度、各半径包围电荷不同 | 各球壳独立积分，监测交叉，不强制后续保持高斯 |
| 固定距离电荷扫描 | 同一初始球形几何，不同总电荷 | 完整球壳解逐面求交；另用弱电荷展开解释尺度 |
| 有限移动窗收支 | Poynting定理及移动控制体积 | 物理通量与粒子做功，加上实际离散时间／算子／覆盖通道分别核对 |

<a id="notation"></a>
## 1. 记号、参考系与物理假设

电子电荷为 $-e$，$e>0$；$m=m_e$；$k_e=1/(4\pi\epsilon_0)$。$N$是物理电子数，$Q=Ne>0$表示总电荷的绝对值，实际总电荷为 $-Q$。$F_0(r_0)$是初始半径 $r_0$ 内的电子数比例，取值0至1。

带撇号的 $t',\boldsymbol r',\boldsymbol p'$ 在束团初始平均静止惯性系S′中定义。它在膨胀后仍是同一个惯性系，不是逐时刻重新定义的局部静止系。束团中心相对lab以 $V=\beta c$ 沿+z运动，$\Gamma=(1-\beta^2)^{-1/2}$，本测试约1174。程序boost系使用另一组 $\Gamma_b\simeq929$，不能与S′混同。

空间电荷参照假设孤立真空、单一电子种类、连续平均场、球对称密度与纯径向流，不额外加入入射自由波。它不描述宏粒子颗粒噪声、碰撞或单电子辐射反作用。场关系无需内部低速假设；下文自相似闭式解另行要求 $|v'|\ll c$。

<a id="vacuum"></a>
## 2. 真空脉冲为什么应平移两格

对仅依赖z、t的横向场 $E_x,B_y$，无源Maxwell方程为

$$\partial_t E_x=-c^2\partial_zB_y,\qquad
\partial_tB_y=-\partial_zE_x.$$

将第一式对t求导，再用第二式消去B，得到

$$\partial_t^2 E_x=c^2\partial_z^2E_x.$$

因此任意形状 $f$ 都能构成沿+z的行波

$$E_x(z,t)=f(z-ct),\qquad B_y(z,t)=f(z-ct)/c.$$

窗口坐标记为 $\xi=z-v_wt$，则场形变成 $f[\xi-(c-v_w)t]$。本测试 $v_w=-c$ 且 $c\Delta t=\Delta z$，故每步场形相对窗口沿+z平移两格。参照代码分别在E和B的原生交错样本上评价此场形，保留B的时间／位置错位，不能先把两者当成同点同刻的数组。

行波有 $u=\epsilon_0E_x^2$、$S_z=E_xB_y/\mu_0=cu$；这给出传播方向和能流符号的独立检查。端点离散面交换与实际交错样本配对核对，其定义仍须与连续 $S_z$ 区分。[Maxwell方程及行波的教材来源：费曼第18章§18–4](https://www.feynmanlectures.caltech.edu/II_18.html#Ch18-S4)。

<a id="sphere-maxwell"></a>
## 3. 球对称自场：位移电流为何不能省略

球对称下 $\boldsymbol E'=E_r'(r,t')\hat{\boldsymbol r}$。以**带符号**的包围电荷

$$Q_{\rm enc}(r,t')=4\pi\int_0^r\rho(s,t')s^2\,ds$$

代入Gauss定律和电荷连续性，得到

$$E_r'=\frac{Q_{\rm enc}}{4\pi\epsilon_0r^2},\qquad
\partial_{t'}Q_{\rm enc}=-4\pi r^2j_r'.$$

所以在固定r处

$$\epsilon_0\partial_{t'}E_r'=-j_r',\qquad
\nabla\times\boldsymbol B'=\mu_0\boldsymbol j'
+\frac1{c^2}\partial_{t'}\boldsymbol E'=0.$$

电流与位移电流恰好抵消。再结合球对称、无磁单极及所选无外加自由波的边界条件，得到 $\boldsymbol B'=0$。这正是[费曼第18章§18–2，式18.5–18.7](https://www.feynmanlectures.caltech.edu/II_18.html#Ch18-S2)讨论的情形。

因此球对称连续模型的径向场是完整Maxwell解，并未通过“删掉位移电流”把传播速度设为无穷。它与因果性相容；球壳内部的电荷变化必须由穿过该位置的径向电流造成。这个特殊结论不能推广为任意非球形束团的瞬时Poisson场都精确。

在该理想模型S′中 $\boldsymbol S'=\boldsymbol E'\times\boldsymbol B'/\mu_0=0$，仍可通过 $\boldsymbol j'\cdot\boldsymbol E'$ 交换场能和粒子动能。变换到lab后一般有磁场和Poynting流，它们本身不等于新增的辐射。

<a id="shell-equations"></a>
## 4. 实际对照采用的相对论球壳方程

壳层交叉前，同一球壳内的粒子数不变，包围电荷为 $-QF_0(r_0)$。电子受到的径向力向外：

$$\frac{dp_r'}{dt'}=\frac{k_eNe^2F_0(r_0)}{r^2},\qquad
\frac{dr}{dt'}=\frac{p_r'}{m\sqrt{1+(p_r'/mc)^2}}.$$

第二式明确保留相对论速度与动量关系。即使初始动量为负、球壳先向内运动，同一组方程仍适用。以 $\gamma_r=\sqrt{1+(p_r'/mc)^2}$ 表示内部运动的洛伦兹因子，有

$$\frac{d(mc^2\gamma_r)}{dt'}=v_r'\frac{dp_r'}{dt'},\qquad
\frac{d}{dt'}\left(mc^2\gamma_r+
\frac{k_eNe^2F_0(r_0)}r\right)=0.$$

于是任一壳层的动能变化满足

$$\boxed{\Delta\kappa'(r_0,t')=
k_eNe^2F_0(r_0)\left(\frac1{r_0}-\frac1{r(t')}\right).}$$

这是壳层轨迹的第一积分，不要求内部非相对论。它不是说所有外层电荷的标势都不变；它只利用外层球壳对内部的合力为零。

为对应现有参考求解器，定义

$$\eta=\frac{k_eNe^2}{mc^2\sigma},\quad
t_0=\sqrt{\frac{m\sigma^3}{k_eNe^2}},\quad
X=r/\sigma,\quad \tau=t'/t_0,\quad
P=p_r'/(mc\sqrt\eta).$$

则实际积分的无量纲方程为

$$\frac{dX}{d\tau}=\frac{P}{\sqrt{1+\eta P^2}},\qquad
\frac{dP}{d\tau}=\frac{F_0}{X^2}.$$

冷球取 $P(0)=0$；聚焦球先指定 $v_r'(0)=h_0r_0$，再用完整的 $p=mv/\sqrt{1-v^2/c^2}$ 设置初动量，本组 $h_0t_0=-0.2$。对照用64个径向积分节点、16个极角节点及32个方位角，共32768个带权样本，同时检查径向排序和壳层第一积分。[已冻结的实际参考求解器](physical_benchmarks/space_charge_scaling_20261030/scripts/model.py)。

<a id="uniform"></a>
## 5. 冷均匀球：由力推到半径、动能和能散

初始均匀球半径 $R_0$，有 $F_0=(r_0/R_0)^3$。在**内部运动非相对论**极限，令 $r=a(t')r_0$，由上一节力方程除以 $mr_0$ 得

$$\ddot a=\frac{\omega^2}{a^2},\qquad
\omega^2=\frac{k_eNe^2}{mR_0^3},\qquad
a(0)=1,\quad\dot a(0)=0.$$

右侧与壳层标签 $r_0$ 无关，因而所有壳层共用同一个a，这就是自相似的原因。乘以 $\dot a$ 后积分：

$$\frac12\dot a^2=\omega^2(1-1/a).$$

再分离变量，得到可解的标量时间关系

$$t'\sqrt{2\omega^2}
=\int_1^a\sqrt{\frac{s}{s-1}}\,ds
=\sqrt{a(a-1)}+\operatorname{asinh}\sqrt{a-1}.$$

径向概率测度为 $dF_0=3r_0^2dr_0/R_0^3$，所以

$$\langle r_0^2\rangle=3R_0^2/5,\quad
\langle r_0^4\rangle=3R_0^4/7,\quad
\sigma_x(t')=aR_0/\sqrt5.$$

将 $F_0$ 与 $r=ar_0$ 代入壳层能量第一积分：

$$\Delta\kappa'=\kappa_{\max}\frac{r_0^2}{R_0^2},\qquad
\kappa_{\max}=\frac{k_eNe^2}{R_0}(1-1/a).$$

由上述二阶、四阶矩直接得到

$$\langle\Delta\kappa'\rangle=\frac35\kappa_{\max},\quad
\operatorname{std}(\Delta\kappa')=\sqrt{\frac{12}{175}}\,\kappa_{\max},\quad
\frac{\operatorname{std}}{\operatorname{mean}}=\sqrt{\frac4{21}}\simeq0.43644.$$

这些闭式系数限于非相对论冷均匀球、S′同时刻及相对共同冷初态的动能增量。固定lab面既有不同过面时刻，又有纵向动量变换，不能直接套用这个比值。实际表中仍使用第4节相对论方程，而非强制精确自相似。[均匀球解析路线的文献对照：Grech等，§II B 2与式27](https://arxiv.org/pdf/1105.0409)。

<a id="focusing"></a>
## 6. 向内聚焦均匀球：为什么可以先损失动能再恢复

在同一非相对论尺度方程中改初值为 $\dot a(0)=h_0<0$，积分常数随之改变：

$$\frac12\dot a^2=\frac12h_0^2+\omega^2(1-1/a).$$

均匀球的静电场能可从球内外的 $E^2$ 积分得到。球内 $|E|=k_eQr/R^3$，球外 $|E|=k_eQ/r^2$，故

$$U(R)=\frac{\epsilon_0}2\int_0^\infty E^2\,4\pi r^2dr
=\frac{k_eQ^2}{2}\left(\frac1{5R}+\frac1R\right)
=\frac{3k_eQ^2}{5R}.$$

同时粒子总径向动能为 $K'=3NmR_0^2\dot a^2/10$。结合上式，有

$$U(t')=U_0/a,\qquad
K'(t')-K'(0)=U_0(1-1/a),\qquad U_0=\frac{3k_eQ^2}{5R_0}.$$

当 $a<1$ 时，内向粒子克服排斥力，动能降低、场能升高，符合 $\Delta K'=-\Delta U$。令转折处 $\dot a=0$，得到

$$a_{\min}=\frac1{1+h_0^2/(2\omega^2)}.$$

这解释了聚焦、转折与回弹的物理方向。它只适用于相干线性聚焦的非相对论自相似模型；实际相对论壳层按第4节各自演化，不要求所有壳层严格同时转折。更一般的内部温度也会破坏共同转折。

<a id="gaussian"></a>
## 7. 冷球形高斯：包围电荷、轨迹与密度

初始每轴RMS为 $\sigma$，取

$$n_0(r)=\frac{N}{(2\pi)^{3/2}\sigma^3}e^{-r^2/(2\sigma^2)}.$$

将其在球内积分，令 $x=r_0/\sigma$，由分部积分得到

$$F_0(r_0)=\sqrt{\frac2\pi}\int_0^x y^2e^{-y^2/2}dy
=\operatorname{erf}(x/\sqrt2)-\sqrt{\frac2\pi}xe^{-x^2/2}.$$

因此初场 $E_r'(r,0)=-k_eQF_0(r)/r^2$；每个壳层把自己的 $F_0(r_0)$ 代入第4节方程。非相对论冷极限下，每个壳层也可令 $r=a(r_0,t')r_0$，得到

$$A(r_0)=\frac{k_eNe^2F_0(r_0)}{mr_0^3},\quad
\ddot a=A(r_0)/a^2,\quad
t'\sqrt{2A}=\sqrt{a(a-1)}+\operatorname{asinh}\sqrt{a-1}.$$

这里A依赖 $r_0$，所以不同壳层没有共同的a，不能假定密度以后仍为高斯。由粒子数守恒 $n(r,t')r^2dr=n_0(r_0)r_0^2dr_0$，得到

$$n(r,t')=n_0(r_0)\frac{r_0^2}{r^2}
\left(\frac{\partial r}{\partial r_0}\right)^{-1},\quad
\sigma_x^2(t')=\frac13\int r(r_0,t')^2\,dF_0.$$

一旦 $\partial r/\partial r_0=0$，单流密度可出现尖峰；交叉后包围粒子数不再固定，必须停止使用这一固定F₀参考。这个约束不是球形Maxwell场关系失效，而是原来的壳层标签不再对应恒定包围电荷。[非均匀球壳与高斯交叉的文献对照：Zerbe等，§IV–V及附录B3](https://arxiv.org/pdf/1712.08629)。

全空间自能也可独立积分。逐层组装电荷的相互作用能为

$$U=k_eQ^2\int_0^1\frac{F_0}{r(F_0,t')}\,dF_0.$$

初始高斯代入 $F_0(x)$、$dF_0=\sqrt{2/\pi}\,x^2e^{-x^2/2}dx$ 后，积分值为

$$U_0=\frac{k_eQ^2}{\sigma}\int_0^\infty\frac{F_0(x)}x\,dF_0(x)
=\frac{k_eQ^2}{2\sqrt\pi\sigma}.$$

中间的无量纲积分用一次分部积分即可核对：$dF_0=\sqrt{2/\pi}\,x^2e^{-x^2/2}dx$，故

$$\int_0^\infty\frac{F_0}x\,dF_0
=\sqrt{\frac2\pi}\int_0^\infty F_0xe^{-x^2/2}dx
=\frac2\pi\int_0^\infty x^2e^{-x^2}dx
=\frac1{2\sqrt\pi}.$$

也可将 $E=-k_eQF_0/r^2$ 代入 $\epsilon_0\int E^2dV/2$ 独立核对。相同初始RMS的均匀球需取 $R_0=\sqrt5\sigma$，不能把两者的几何半径直接当成同一个量。

<a id="lab-planes"></a>
## 8. 从世界线到固定lab面：理论值怎样与输出对应

先在S′生成共同初态，再对每个粒子沿其固定方向积分径向世界线。设该粒子极角余弦为 $\mu_i$，$z_i'(t')=\mu_i r_i(t')$。采用与当前参考求解器等价的平移原点，令中心通过lab入口面0 m对应S′时刻 $T_0$：

$$z_{\rm lab}=\Gamma[z'+\beta c(t'-T_0)],\qquad
t_{\rm lab}=\Gamma[(t'-T_0)+\beta z'/c].$$

在每个固定lab面 $z=z_j$，求解标量方程

$$\boxed{\mu_i r_i(t'_{ij})+\beta c(t'_{ij}-T_0)=z_j/\Gamma.}$$

得到事件 $t'_{ij}$ 后，分别写入该事件的横向位置、lab时间和lab动量：

$$x_{ij}=r_i(t'_{ij})\hat r_{i,x},\quad y_{ij}=r_i(t'_{ij})\hat r_{i,y},$$

$$\gamma_{ij}^{\rm lab}=\Gamma(\gamma_{ij}'+\beta u_{z,ij}'),\quad
u_{z,ij}^{\rm lab}=\Gamma(u_{z,ij}'+\beta\gamma_{ij}'),\quad
\boldsymbol u_{\perp,ij}^{\rm lab}=\boldsymbol u_{\perp,ij}'.$$

这里 $\boldsymbol u=\boldsymbol p/(mc)$。同一面上不同粒子的 $t'_{ij}$ 一般不同，因此不是S′同时刻快照。输入粒子已经经过理论从共同初态到入口事件的预演化，不能在粒子各自进入时把动量再归零。

逐ID、同权重计算

$$\Delta K_{ij}^{\rm lab}=mc^2(\gamma_{ij}^{\rm lab}-\gamma_{i0}^{\rm lab})
=\Gamma mc^2(\Delta\gamma_{ij}'+\beta\Delta u_{z,ij}'),$$

$$\overline{\Delta K}_j=\sum_iw_i\Delta K_{ij},\qquad
\sigma_{\Delta K,j}^2=\sum_iw_i(\Delta K_{ij}-\overline{\Delta K}_j)^2,\qquad\sum_iw_i=1.$$

采用总体加权标准差，不用样本ddof修正。能散变化还包括协方差：

$$\operatorname{Var}(\Delta K_{\rm lab})=\Gamma^2m^2c^4\left[
\operatorname{Var}(\Delta\gamma')+\beta^2\operatorname{Var}(\Delta u_z')
+2\beta\operatorname{Cov}(\Delta\gamma',\Delta u_z')\right].$$

这就是为什么不能只将静止系均值或能散乘上Γ，也不能由球形同时刻的总纵动量为零，推断固定lab面平均的 $\Delta u_z'$ 必为零。

本记录保留两套参照：**固定面主口径**对理论自己的过面事件求值；**实际事件诊断口径**则在PIC输出事件对应的S′时刻评价理想世界线。两者时刻通常不相同，先前−620 eV与−1048 eV的区别就在这里，不能混成一次误差变化。

场也应在同一事件变换。S′中B′=0时，$E_\parallel=E_\parallel'$、$\boldsymbol E_\perp=\Gamma\boldsymbol E_\perp'$、$\boldsymbol B=\Gamma\boldsymbol V\times\boldsymbol E'/c^2$。lab磁场与能流可以只是束缚场随束团移动的表现。[场变换的教材对照：费曼第26章](https://www.feynmanlectures.caltech.edu/II_26.html)。

<a id="weak-charge"></a>
## 9. 为什么冷球的弱电荷平均信号近似为Q平方

这一节自行展开固定lab面事件，解释尺度；表中的主理论值始终来自完整相对论壳层解。所有粒子正向穿面，比较同一粒子集合。固定 $R_0,\Gamma,T_0,z_j$，令 $\lambda=\omega^2\propto Q$，在 $\lambda t'^2\ll1$、内部速度远小于c的极限，从 $\ddot a=\lambda/a^2$ 依次解得

$$a(t')=1+\frac\lambda2t'^2-\frac{\lambda^2}{12}t'^4+O(\lambda^3),\quad
v_r'=r_0\left(\lambda t'-\frac{\lambda^2}3t'^3\right)+O(\lambda^3).$$

一阶冲量为O(Q)，但一阶平均纵向项因 $\langle z_0\rangle=0$ 抵消。要算留下的二阶量，固定面事件的变化也必须保留。定义

$$T_j=T_0+\frac{z_j}{\Gamma\beta c},\qquad
s_j=T_j-\frac{z_0}{\beta c}.$$

将a展开代入过面方程，得到

$$t'_{ij}=s_j-\lambda\frac{z_0s_j^2}{2\beta c}+O(\lambda^2).$$

再将此时间代入 $K_{\rm lab}-\text{常数}=\Gamma m(v'^2/2+\beta cv_z')+O(\lambda^3)$，记 $\Delta s^n=s_1^n-s_0^n$，有

$$\frac{\Delta K_{\rm lab}}{\Gamma m}
=\lambda\beta cz_0\Delta s
+\lambda^2\left[\frac{r_0^2}2\Delta s^2-\frac{z_0^2}2\Delta s^2
-\frac{\beta cz_0}3\Delta s^3\right]+O(\lambda^3).$$

均匀球有 $\langle r_0^2\rangle=3R_0^2/5$、$\langle z_0^2\rangle=R_0^2/5$，所有关于 $z_0$ 的奇次矩为零。令 $A=T_1^2-T_0^2$，展开 $s_j=T_j-z_0/(\beta c)$ 后得到

$$\langle z_0\Delta s\rangle=0,\quad
\langle r_0^2\Delta s^2\rangle=\frac35R_0^2A,\quad
\langle z_0^2\Delta s^2\rangle=\frac15R_0^2A,\quad
\langle z_0\Delta s^3\rangle=-\frac{3R_0^2A}{5\beta c}.$$

所以入口0 m至出口 $z_1$ 的二阶平均值可分成三部分：

| 来源 | 对 $\Gamma m\lambda^2R_0^2(T_1^2-T_0^2)$ 的系数 |
|---|---:|
| 内部动能 $mv'^2/2$ | $3/10$ |
| 过面时间改变带来的纵向动量项 | $-1/10$ |
| 膨胀后场强减弱造成的二阶纵向动量项 | $1/5$ |
| 合计 | $2/5$ |

于是，在上述固定初态和固定面条件下，首项为

$$\boxed{\langle\Delta K_{\rm lab}\rangle
=\frac25\Gamma m\omega^4R_0^2(T_1^2-T_0^2)+O(Q^3)\ \propto Q^2.}$$

此前辅助脚本曾写作 $(3+\beta^2)/10$ 的系数；重新完整保留固定面事件修正后应为 $4/10$，旧式不应当作任意β的精确首项。本次 $\Gamma\simeq1174$ 时两者相对差约 $1.81\times10^{-7}$，对最弱档的辅助估计仅约 $5.27\times10^{-8}$ eV；不影响0.1806 eV偏差或Q平方结论。原[辅助估计记录](physical_benchmarks/space_charge_scaling_20261030/WEAK_LIMIT_CHECK.json)保持历史原样，本文另列修正，所有完整相对论基准及PIC表值均不变。另用独立径向相对论ODE在β=0.6、0.9和当前束流速度验证弱场极限系数，并核对均匀球能量矩与高斯自能的两种积分：[理论公式核对记录](metadata/THEORY_CHECKS_20261030.json)。

另一方面，如果平均场／冲量偏差含未抵消的一阶项 $\delta\langle\Delta u_z'\rangle\propto Q$，则 $\delta\langle\Delta K_{\rm lab}\rangle\simeq\Gamma\beta mc^2\delta\langle\Delta u_z'\rangle$ 可以近似正比Q，弱信号相对偏差便可能随1/Q增大。这是允许的尺度分析；实测指数约0.941并未唯一定位一阶项来自哪一操作。

<a id="energy-budget"></a>
## 10. 能量收支公式，以及它能检验什么

从Ampère–Maxwell方程点乘E、Faraday方程点乘 $\boldsymbol B/\mu_0$ 后相加，利用

$$\nabla\cdot(\boldsymbol E\times\boldsymbol B)
=\boldsymbol B\cdot\nabla\times\boldsymbol E
-\boldsymbol E\cdot\nabla\times\boldsymbol B,$$

得到Poynting定理

$$\partial_tu+\nabla\cdot\boldsymbol S=-\boldsymbol J\cdot\boldsymbol E,\quad
u=\frac{\epsilon_0E^2}2+\frac{B^2}{2\mu_0},\quad
\boldsymbol S=\frac{\boldsymbol E\times\boldsymbol B}{\mu_0}.$$

磁力不做功，因为 $\boldsymbol v\cdot(\boldsymbol v\times\boldsymbol B)=0$；粒子自场做功与右侧的场能损失符号相反。孤立球在S′全空间有 $d(K'+U)/dt'=0$，上述球形自能与能量第一积分互相对应。[局部能量守恒及推导的教材来源：费曼第27章§27–2至§27–3](https://www.feynmanlectures.caltech.edu/II_27.html#Ch27-S3)。

对移动体积 $V(t)$，边界速度为 $\boldsymbol v_w$，Reynolds公式给出

$$\frac{d}{dt}\int_{V(t)}u\,dV
=-\oint_{\partial V(t)}(\boldsymbol S-u\boldsymbol v_w)\cdot\boldsymbol n\,dA
-\int_{V(t)}\boldsymbol J\cdot\boldsymbol E\,dV.$$

因此“穿过移动边界”的场能交换不是单独的 $\boldsymbol S\cdot\boldsymbol n$，还含边界扫过的能量存量。粒子若也按体积统计，须同样计入相对移动边界的粒子动能流入流出；若用稳定cohort，则将保留和已移除粒子一起记账。入口重写与离散裁剪要按实际操作的前后场能差记录，不能漏记或重复记为额外做功。

连续定理还不能直接替代交错显式程序的离散恒等式。例如任一原生场分量有

$$\frac{\epsilon_0}2\bigl(\|E+\Delta E\|_h^2-\|E\|_h^2\bigr)
=\epsilon_0\langle E,\Delta E\rangle_h+\frac{\epsilon_0}2\|\Delta E\|_h^2.$$

其中内积含相应网格权重，磁场有对应式；时间层选择还会改变J应与哪个E配对。把实际curl更新、Cowan修正、CPML、入口覆盖和移窗各阶段代入，当前场账本使用的带符号核对式为

$$R_F=\Delta U_h+Q_{\rm curl}+W_{J,\rm mid}
-T_E-T_B-C_{\rm pair}-C_{\rm Cowan}
-\Delta U_{\rm CPML}-\Delta U_{\rm inlet}-\Delta U_{\rm remap}.$$

这里 $Q_{\rm curl}$ 是向外的离散curl面交换，$T_E,T_B,C_{\rm pair},C_{\rm Cowan}$ 是实际离散时间／算子通道。它们的定义来自实际更新，不是新引入的物理能量来源，也不能加总成“辐射”。场账本闭合支持场算子收支一致；仍要单独检查粒子获得的功与网格电流失去的功是否对应。

最后，全空间球形自能不等于有限矩形窗内的场能。对全部电荷包在半径L内的球形分布，L外场能为 $k_eQ^2/(2L)$；均匀球在自身半径R外仍有总自能的5/6。有限窗的存量、通量、源项必须在同一个体积和同一参考系内核对。这里的离散账本与独立球壳动力学是互补的检验。

<a id="scope"></a>
## 11. 怎样读本记录中的“理论符合”

真空平移检验已知行波；球壳对照检验连续平均场空间电荷动力学；固定面统计检验同一物理探测事件口径；能量账本检验实际离散收支。球对称下位移电流抵消是完整场关系，自相似与弱电荷展开有额外近似条件，不能互相替代。

因此，当前球形测试保留主要膨胀与聚焦行为，是明确的阶段性证据；仍存在的均值、能散及粒子—网格差额也须保留。非球形长束团、有磁辐射、近场分离和更弱辐射能损需要各自的适用性及收敛检查，不能由这里的球形特殊解直接推出精度保证。
