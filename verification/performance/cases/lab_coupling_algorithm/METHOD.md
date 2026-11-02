# 低成本能量—动量配对：方法与边界

本轮是在独立测试副本里实际改变粒子自磁场取样，生产程序保持。
不是把数值余差加回粒子，也不是重新定义账本后声称推进已修复。

## 1. 当前半隐式方法已经解决什么

[Chen等的半隐式电磁PIC论文](https://arxiv.org/html/1903.01565v2)把显式leapfrog场与中点粒子—电流耦合结合，用一致取场/沉积保证电功交换。论文的场守恒推导使用Yee算子；不能直接当作当前Cowan的证明。

本程序已采用局部平均速度方程

\[
\Delta u=\frac{qh}{mc}(E+\bar v\times B),\qquad
\Delta x=h\bar v,\qquad
\bar v=c\frac{u_1+u_0}{\gamma_1+\gamma_0}.
\]

因此每个子步有\(mc^2\Delta\gamma=qE\cdot\Delta x\)。同一Maxwell步内自E/B样本固定，外磁场照原规则子循环；自电功与自磁力冲量分别对整步位移望远镜求和。

电场已有与整步直线电流沉积伴随的轨迹平均，磁场原来仍是轨迹中点CIC。电功配对不自动保证磁力动量配对。

## 2. 自磁场取样的可实现修正

以原生交错格点的纵向动量函数为参照，记A_z为相邻z层平均，则

\[
cP_z=\epsilon_0c\langle A_zE_x,B_y\rangle_h
      -\epsilon_0c\langle A_zE_y,B_x\rangle_h.
\]

Ampere电流增量\(\delta E=-\Delta tJ/\epsilon_0\)带来的场动量变化为

\[
\delta(cP_z)_J=-c\Delta t\Delta V\sum_g
       [J_x(A_z^*B_y)-J_y(A_z^*B_x)]_g.
\]

原电流沉积沿整步直线弦\(x(s)=x_0+s\Delta x\)，使用沿电流方向零阶、横向一阶的权重。把这些相同权重记为D_x、D_y，令

\[
B_y^{(p)}=\int_0^1 D_x^*[A_z^*B_y](x(s))\,ds,
\qquad
B_x^{(p)}=\int_0^1 D_y^*[A_z^*B_x](x(s))\,ds.
\]

这使真实粒子自磁力冲量满足

\[
cI_{z,B}=qc(\Delta x B_y^{(p)}-\Delta yB_x^{(p)})
        =-\delta(cP_z)_J.
\]

仍使用完整三维\(\bar v\times B^{(p)}\)，所以磁力无功，不需要投影粒子动量或事后能量修正。B_z保留原插值；本轮只针对boost方向的动量配对，没有声称三方向动量都严格配对。

轨迹在格内的权重积分是两个线性因子的积分：中点乘积加\(\pm\Delta\eta\Delta\zeta/12\)。跨格时按原直线弦分段。并未将电流沉积改为粒子子步折线。

该恒等式对完整物理域内的活动弦成立。本轮真实全程序分析只用全束释放、尚无退出的区间。碰到真实外边界裁剪仍回退原处理；没有宣称完整注入/载流子/退出已验收。

### MPI与成本

跨rank弦的伴随取样需要两层Bx/By。第一轮全程序验证识别了缺少第二层，修正版把原磁halo的Bx/By消息加长为两层，Bz与电场halo保持原范围。消息次数和tag保持，未改MPI网络或停止方式。

新增计算随活动粒子数和已有耦合迭代数线性增长，没有新全局场迭代。标量取场函数在制造场微测试中约为原函数2.76倍；这不是整步倍数。真实短卡迭代次数保持，整步成本另用同节点ABBA测试，不能外推全量工况。

## 3. Cowan能量项为什么不能只靠多迭代粒子

本机Faraday使用修改后的curl C_*，Ampere使用标准curl C。
[WarpX官方NSFDTD说明](https://warpx.readthedocs.io/en/latest/theory/models_algorithms/explicit_em_pic.html)给出相同结构，并采用修改后的磁Gauss约束\(\nabla^*\cdot B=0\)。标准\(\nabla\cdot B\)不为零不能直接等同真实磁荷。

下列是本轮对本机算子的独立推导与数值核验，不声称论文已为本程序证明。

在均匀周期域的Fourier表示中，设标准导数符号为k_i，Cowan平滑符号为s_i，\(A=\operatorname{diag}(s_x,s_y,s_z)\)，则

\[
C=i[k]_{\times},\quad C_*=i[Ak]_{\times}
       =\det(A)A^{-1}CA^{-1}.
\]

因此局部空间权重

\[
M_E=\operatorname{diag}(s_ys_z,s_zs_x,s_xs_y),\qquad M_B=A
\]

满足\(M_EC=C_*^\dagger M_B\)。这给出配套的场能量形式；leapfrog的磁能还需采用前后半时刻乘积。

但场电流功是\(J\cdot M_EE\)，不是现有的\(J\cdot E\)。所以单改账本不够：若保留原Cowan并使用这个配对，还必须对应改变粒子取场、重新验证空间电荷与辐射反馈。周期域的配套纵向动量含\(E_xs_yB_y-E_ys_xB_x\)，也不是现有无权动量。

512组任意Fourier方向的伴随关系误差约2.06e−16；满足各自Gauss约束时，配套真空纵向动量关系误差约1.79e−16。这只是算子关系，不是束团物理精度。

## 4. 可继续试验的局部配对显式场模板

另一个实验方向是保留标准Ampere和现有时间步，让Faraday成为\(-\nabla\times SE\)，其中S是同一标量、自伴随、局部平滑算子。

\[
B_+=B_- -\Delta t\nabla_h\times SE_n,\qquad
E_{n+1}=E_n+c^2\Delta t\nabla_h\times B_+-\Delta tJ/\epsilon_0.
\]

定义

\[
H_S^n=\frac{\epsilon_0}{2}\langle E_n,SE_n\rangle_h+
       \frac1{2\mu_0}\langle B_-,B_+\rangle_h.
\]

在周期均匀域中，S与curl可交换，因而

\[
\Delta H_S+\Delta t\langle J,S(E_n+E_{n+1})/2\rangle_h=0.
\]

粒子应使用与原电流伴随的\(S(E_n+E_{n+1})/2\)取场，才能将这个场功与粒子功对上。这里需要真实改变场推进与取场；不能仅用新能量函数隐藏原数值项。

为满足生产纵向细网格和\(c\Delta t=\Delta z\)，本轮推导的候选为

\[
S(k)=1-q_z(r_xq_x+r_yq_y),\quad
q_i=\sin^2(k_i\Delta i/2),\quad r_i=(\Delta z/\Delta i)^2.
\]

物理空间可写成13个格点的局部模板：

\[
S=I-\frac{\Delta z^4}{16}
       (L_xL_z+L_yL_z),
\]

L_i为标准二阶离散Laplacian。沿z轴\(q_x=q_y=0\)，S=1，保持Courant一步的轴向无色散。该项对长波是四阶小量。

设\(t=r_xq_x+r_yq_y\)、\(z=q_z\)，真空CFL符号为\((z+t)(1-tz)\)。当\(r_x+r_y\le\sqrt2-1\)时，它随z单调增加，最大值为\(1-t^2\le1\)。生产比例和为0.125867，满足这个充分条件，且S正。

这里只证明传播符号界；Nyquist/Courant极限仍是边缘情形。半时刻交叉能量的闭合不能替代高频模稳定性和物理精度检查。

独立局部有限差分原型（没有FFT推进或全空间迭代）已做160步周期真空与32步给定电流检查：能量关系2.03e−15、纵向动量/初始能量6.87e−16、给定电流功关系1.97e−15。尚无PIC粒子、CPML、注入或快参考，不能接管生产程序。

## 5. 推荐验证顺序

先保留本轮磁取场候选作为对照，继续在独立小原型比较原Cowan兼容权重与局部配对场模板。需要同时看lab收支、物理场能近似、相空间、已有冷球/高斯理论对照、辐射色散和耗时。然后才接入原CPML、快参考、移动窗与多节点，并复验出入口过程。

只把机器精度的恒等式检查列为程序/数学验收；真实辐射与空间电荷对照仍单独判断。本轮不启动新的4.5 m大仿真。
