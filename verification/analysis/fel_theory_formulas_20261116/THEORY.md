# 谢明公式与空间电荷抑制聚束的理论推导

讨论两件事：谢明公式如何估计 FEL 的功率增益长度，以及纵向空间电荷如何在冷束近似下产生聚束恢复力。谢明的三维修正是数值拟合公式；空间电荷部分则从屏蔽势方程推到密度调制的演化方程。

## 一 谢明公式

### 1 共振条件与耦合参数

考虑平面波荡器，周期为 \(\lambda_u\)，峰值偏转参数为 \(K\)，平均电子能量为 \(\gamma m_ec^2\)。定义

\[
k_u=\frac{2\pi}{\lambda_u},\qquad
A_u=1+\frac{K^2}{2}.
\]

超相对论电子的周期平均轴向速度近似为

\[
\overline\beta_z\simeq1-\frac{A_u}{2\gamma^2}.
\]

令辐射相位与波荡器相位之和为

\[
\psi=(k_r+k_u)z-\omega_rt,\qquad \omega_r=ck_r.
\]

共振时 \(d\psi/dz=0\)，于是

\[
k_u\simeq\frac{k_rA_u}{2\gamma^2},\qquad
\boxed{\lambda_r=\frac{\lambda_uA_u}{2\gamma^2}}.
\]

平面波荡器的基频耦合采用

\[
a_w=\frac K{\sqrt2},\qquad
\xi=\frac{a_w^2}{2(1+a_w^2)}
=\frac{K^2}{4+2K^2},
\]

\[
JJ=J_0(\xi)-J_1(\xi),\qquad A_w=a_wJJ.
\]

其中 \(J_0,J_1\) 是第一类 Bessel 函数；\(a_w\) 是 RMS 约定，不能直接与峰值 \(K\) 混用。

### 2 一维 Pierce 参数与功率增益长度

设电流为 \(I\)，圆高斯束在单个横向方向的 RMS 尺寸为 \(\sigma\)。Alfvén 电流为

\[
I_A=\frac{4\pi\epsilon_0m_ec^3}{e}\simeq17045\ \mathrm A,
\]

其中 \(e>0\) 是元电荷大小。采用谢明原文的归一化，

\[
\boxed{
\rho=\left[
\frac I{I_A}
\left(\frac{\lambda_uA_w}{2\pi\sigma}\right)^2
\frac1{(2\gamma)^3}
\right]^{1/3}}
\]

也可等价写成

\[
\rho^3=\frac{I}{I_A}
\frac{(KJJ)^2}{16\gamma^3(k_u\sigma)^2}.
\]

冷束、共振、一维线性化的 Maxwell–Lorentz 方程经过场幅和能量调制归一化后，可写为

\[
\frac{d\mathcal A}{dz}=g_Fb,\qquad
\frac{db}{dz}=-ig_Fp,\qquad
\frac{dp}{dz}=-g_F\mathcal A,
\qquad g_F=2\rho k_u.
\]

这里 \(\mathcal A\) 是归一化复场幅，\(b=\langle e^{-i\psi}\rangle\) 是聚束因子，\(p\) 是按这组方程归一化的能量谐波。与第二部分的定义对应时，常系数下 \(p=(\kappa/g_F)g\)。消去 \(b,p\)，得到

\[
\mathcal A'''=ig_F^3\mathcal A.
\]

取 \(\mathcal A\propto e^{\mu z}\)，增长支满足

\[
\mu^3=ig_F^3,\qquad
\operatorname{Re}\mu=\frac{\sqrt3}{2}g_F.
\]

因此功率的增长率是场幅实增长率的两倍：

\[
P\propto|\mathcal A|^2
\propto e^{\sqrt3g_Fz}=e^{z/L_0},
\]

\[
\boxed{L_0=\frac1{\sqrt3g_F}
=\frac{\lambda_u}{4\pi\sqrt3\rho}}.
\]

这里 \(L_0\) 是功率增益长度，场幅的 e 倍增长长度为 \(2L_0\)。

### 3 衍射 发射度和能散的三维修正

设几何 RMS 发射度为 \(\epsilon\)，归一化发射度为 \(\epsilon_n\)，匹配 beta 函数为 \(\beta\)。超相对论近似下

\[
\epsilon\simeq\frac{\epsilon_n}{\gamma},\qquad
\sigma^2=\beta\epsilon.
\]

定义衍射长度和三个无量纲参数：

\[
L_R=\frac{4\pi\sigma^2}{\lambda_r},\qquad
\eta_d=\frac{L_0}{L_R},
\]

\[
\eta_e=\frac{L_0}{\beta}\frac{4\pi\epsilon}{\lambda_r},\qquad
\eta_g=\frac{4\pi L_0}{\lambda_u}\frac{\sigma_E}{E}.
\]

\(\sigma_E/E\) 是非相关 RMS 相对能散。这些参数分别比较一维增益长度与衍射、横向相位混合和能量失谐的尺度。谢明对三维线性增长率的数值解作拟合，得到

\[
\begin{aligned}
\Lambda={}&0.45\eta_d^{0.57}
+0.55\eta_e^{1.6}
+3\eta_g^2\\
&+0.35\eta_e^{2.9}\eta_g^{2.4}
+51\eta_d^{0.95}\eta_g^3\\
&+5.4\eta_d^{0.7}\eta_e^{1.9}
+1140\eta_d^{2.2}\eta_e^{2.9}\eta_g^{3.2},
\end{aligned}
\]

\[
\boxed{L_g=L_0(1+\Lambda)}.
\]

七项中的系数和指数来自拟合，不是由上一节的一维三阶方程逐项解析推出的；一维方程给出基础尺度，三维修正描述该尺度的偏移。

常用饱和功率近似为

\[
\boxed{P_{\mathrm{sat}}\simeq
1.6\rho\left(\frac{L_0}{L_g}\right)^2P_{\mathrm{beam}}},\qquad
P_{\mathrm{beam}}\simeq I\frac{\gamma m_ec^2}{e}.
\]

若以有效种子 \(P_{\mathrm{seed}}\) 描述指数增长，则粗略有

\[
s_{\mathrm{sat}}\simeq
L_g\ln\frac{P_{\mathrm{sat}}}{P_{\mathrm{seed}}}.
\]

\(s\) 是累计增益磁长。实验室饱和位置还要加上已经经过的段间空隙；随机 SASE 的有效种子与启动过程也影响饱和位置，不能仅由 \(\rho\) 唯一确定。

### 4 使用条件与本例代入

上述三维拟合采用纵向均匀、横向和能量分布为高斯的束流，圆束、相同横向发射度、恒定匹配 beta，并对线性增长的失谐作优化。有限束长、滑移、包络失配和完整空间电荷反馈不由这一拟合自动给出。

本例采用

\[
I=53.27545\ \mathrm A,\quad
\gamma=298.26088,\quad K=2.060075,
\]

\[
\lambda_u=0.028\ \mathrm m,\quad
\lambda_r=491.317\ \mathrm{nm},\quad
\frac{\sigma_E}{E}=2\times10^{-4}.
\]

若明确以 x 方向参数构造匹配圆束近似，

\[
\sigma=132\ \mu\mathrm m,\qquad
\epsilon_n=2.5\ \mu\mathrm m,
\]

则得到

\[
\rho\simeq0.00284435,\qquad
L_0\simeq0.452278\ \mathrm m,\qquad
\boxed{L_g\simeq0.670464\ \mathrm m}.
\]

实际束流为椭圆束且包络变化。这是明确约定下的圆束估算；若改用几何平均束斑和发射度，\(L_g\simeq0.612546\ \mathrm m\)，不能把两种约定当作实际三维解的误差上下界。

## 二 纵向空间电荷如何作用于 BF

### 1 从相位滑移得到密度与能量谐波的关系

在第一节的共振相位中加入小能量偏移 \(\delta\gamma\)。将轴向速度展开到一阶，有

\[
\frac{d\psi}{dz}\simeq\kappa\delta\gamma,\qquad
\kappa=\frac{k_rA_u}{\gamma^3}
\simeq\frac{2k_u}{\gamma}.
\]

定义

\[
b=\left\langle e^{-i\psi}\right\rangle,\qquad
g=\left\langle\delta\gamma e^{-i\psi}\right\rangle.
\]

对 \(b\) 求导，得到

\[
\boxed{\frac{db}{dz}=-i\kappa g}.
\]

这说明能量调制通过轴向速度差转化成密度调制。

### 2 从屏蔽势得到纵向场

采用横向归一化椭圆高斯分布

\[
f(x,y)=\frac{1}{2\pi\sigma_x\sigma_y}
\exp\left(-\frac{x^2}{2\sigma_x^2}
-\frac{y^2}{2\sigma_y^2}\right),\qquad
\int f\,d^2r=1.
\]

小密度扰动写成

\[
n=n_0 f\left(1+be^{i\psi}+b^*e^{-i\psi}\right).
\]

忽略慢变纵向包络，取 \(v_z\simeq c\)，对应的电荷复谐波近似为

\[
\rho_k(x,y)\simeq-\frac Ic\,b f(x,y).
\]

以周期平均准静态模型估计纵向非辐射场：

\[
\gamma_z=\frac\gamma{\sqrt{A_u}},\qquad
a=\frac{k_r}{\gamma_z},
\]

\[
(\nabla_\perp^2-a^2)\Phi_k=-\frac{\rho_k}{\epsilon_0},\qquad
E_{z,k}=-\frac{ik_r}{\gamma_z^2}\Phi_k.
\]

在近共振估算中用 \(k_r\) 近似密度相位的纵向波数 \(k_r+k_u\)。对横向作 Fourier 变换，

\[
\widetilde f(q_x,q_y)
=\exp\left[-\frac{\sigma_x^2q_x^2+\sigma_y^2q_y^2}{2}\right],
\]

\[
\widetilde\Phi_k
=-\frac{Ib}{c\epsilon_0}
\frac{\widetilde f}{q_\perp^2+a^2}.
\]

以粒子横向分布作平均，得到

\[
\overline E_{z,k}=\int fE_{z,k}\,d^2r
=i\frac{Ik_r}{c\epsilon_0\gamma_z^2}\,b\,Q,
\]

\[
Q=\int\frac{d^2q}{(2\pi)^2}
\frac{|\widetilde f|^2}{q_\perp^2+a^2}.
\]

用积分恒等式

\[
\frac1{q_\perp^2+a^2}
=\int_0^\infty e^{-t(q_\perp^2+a^2)}\,dt
\]

完成两个高斯积分，再令 \(u=a^2t\)，可得

\[
Q=\frac{\mathcal R}{4\pi a^2\sigma_x\sigma_y},
\]

\[
\boxed{
\mathcal R=\sigma_x\sigma_y\int_0^\infty
\frac{e^{-u}\,du}
{\sqrt{(\sigma_x^2+u/a^2)(\sigma_y^2+u/a^2)}}}.
\]

在宽束极限 \(a\sigma_x,a\sigma_y\gg1\) 下，\(\mathcal R\to1\)。圆束时，令 \(X=a^2\sigma^2\)，则

\[
\mathcal R=Xe^X E_1(X),\qquad
E_1(X)=\int_X^\infty\frac{e^{-t}}t\,dt.
\]

这里使用束流加权平均场，不是轴上场；对宽高斯束，两种选法给出的恢复力系数存在约两倍差别。

### 3 从纵向电场得到聚束恢复力

电子的能量变化率近似为

\[
\frac{d\gamma}{dz}\simeq-\frac{eE_z}{m_ec^2}.
\]

对复谐波投影，并在线性冷束近似下忽略二阶能量矩项，得到

\[
\frac{dg}{dz}=-iSb+F_{\mathrm{FEL}},
\]

\[
S=\frac{4\pi Ik_r}{I_A\gamma_z^2}Q
=\boxed{\frac{I\mathcal R}{I_Ak_r\sigma_x\sigma_y}}.
\]

其中 \(F_{\mathrm{FEL}}\) 表示辐射场与波荡器共同产生的能量谐波驱动。先将该驱动置零，固定束斑与系数，则

\[
b''=-i\kappa g'
=-\kappa Sb,
\]

\[
\boxed{b''+k_{\mathrm{SC}}^2b=0},\qquad
\boxed{
k_{\mathrm{SC}}^2=\kappa S
=\frac{IA_u\mathcal R}
{I_A\gamma^3\sigma_x\sigma_y}}.
\]

一般初值的解为

\[
b(z)=b_0\cos(k_{\mathrm{SC}}z)
-i\frac{\kappa g_0}{k_{\mathrm{SC}}}
\sin(k_{\mathrm{SC}}z).
\]

当初始只有密度调制而没有能量谐波，即 \(g_0=0\) 时，

\[
\boxed{\frac{b(L)}{b_0}=\cos(k_{\mathrm{SC}}L)},
\]

\[
\frac{b(L)}{b_0}\simeq
1-\frac12\frac{IA_u\mathcal R}
{I_A\gamma^3\sigma_x\sigma_y}L^2
\quad(k_{\mathrm{SC}}L\ll1).
\]

恢复力符号为正，说明纯密度调制会先被削弱。这是可逆的相空间振荡；若初始能量谐波 \(g_0\) 的相位不同，也可能先增强聚束。

### 4 变化束斑的推广

令 \(h=ig\)，则无辐射驱动时可直接积分

\[
\boxed{
\frac{db}{dz}=-\kappa(z)h,\qquad
\frac{dh}{dz}=S(z)b}.
\]

若消成二阶方程，必须保留变系数项：

\[
b''-\frac{\kappa'}\kappa b'+\kappa Sb=0.
\]

入口 Twiss 协方差为

\[
\Sigma_i(0)=\epsilon_i
\begin{pmatrix}
\beta_i&-\alpha_i\\
-\alpha_i&(1+\alpha_i^2)/\beta_i
\end{pmatrix},\qquad
\epsilon_i\simeq\frac{\epsilon_{ni}}\gamma.
\]

在线性光学近似下

\[
\Sigma_i(z)=M_i(z)\Sigma_i(0)M_i(z)^{\mathsf T},\qquad
\sigma_i^2(z)=[\Sigma_i(z)]_{11}.
\]

首段估算对 x 取漂移输运，对 y 加入平面波荡器的周期平均自然聚焦，

\[
k_{\beta y}^2(z)\simeq\frac{K(z)^2k_u^2}{2\gamma^2},
\]

再将 \(\sigma_x(z),\sigma_y(z)\) 代回 \(S(z)\)。这仍假设高斯形状和给定线性包络，未自洽求解所有横向空间电荷与非线性形变。

### 5 三点核的影响

三点核系数为 \((1/4,1/2,1/4)\)，其 Fourier 符号为

\[
F_3(\theta)=\frac14e^{-i\theta}+\frac12+\frac14e^{i\theta}
=\cos^2\frac\theta2.
\]

若源沉积与匹配取场各作用一次，则在单个未混叠空间模式的简化估计中

\[
\boxed{S_{F3}\simeq F_3(\theta)^2S
=\cos^4\frac\theta2\,S}.
\]

设沿 z 的 boost 速度为 \(\beta_Bc\)，Lorentz 因子为 \(\Gamma_B\)。纵向波数按

\[
k'=\Gamma_B\left(k-\beta_B\frac\omega c\right)
\]

变换。对前向光波 \(\omega=ck\)，有

\[
k'_{\mathrm{rad}}=\Gamma_Bk(1-\beta_B).
\]

对随束密度模式 \(\omega\simeq\overline\beta_zck\)，有

\[
k'_{\mathrm{SC}}\simeq
\Gamma_Bk(1-\beta_B\overline\beta_z).
\]

在本例接近平均束速的 boost 中，后者波数约为前者的两倍。因此前向光波约 20 格每周期时，空间电荷密度模式约为 10 格每周期，

\[
\theta_{\mathrm{SC}}\simeq\frac{2\pi}{10},\qquad
F_3(\theta_{\mathrm{SC}})^2\simeq0.81813.
\]

这个因子只表示明确的三点核源与取场响应，没有包含完整离散 Maxwell 算子、别名与所有横向响应，也不能直接移用于五次补偿方法。

### 6 参数代入与约百分之五的尺度

取

\[
I=53.27545\ \mathrm A,\quad
\gamma=298.26088,\quad K=2.060075,
\]

\[
\sigma_x=132\ \mu\mathrm m,\quad
\sigma_y=75\ \mu\mathrm m,\quad
\lambda_r=491.317\ \mathrm{nm},\quad L=2.156\ \mathrm m.
\]

周期平均核给出

\[
\gamma_z\simeq168.8041,\qquad
\mathcal R\simeq0.980397,
\]

\[
\kappa\simeq1.504719\ \mathrm{m^{-1}},\qquad
k_{\mathrm{SC}}\simeq0.190839\ \mathrm{m^{-1}}.
\]

固定束斑的纯密度模式为

\[
\frac{b(L)}{b_0}
=\cos(0.411448)\simeq0.916543,
\]

即幅度降低约 \(8.35\%\)。代入变化束斑，再加入三点核，既有模型积分分别给出

\[
\left.\frac{b(L)}{b_0}\right|_{\mathrm{变化束斑}}
\simeq0.936890,
\]

\[
\boxed{
\left.\frac{b(L)}{b_0}\right|_{\mathrm{变化束斑+F3}}
\simeq0.948270},
\]

对应约 \(6.31\%\) 和 \(5.17\%\) 的幅度抑制。

### 7 与辐射能量的关系

即使取最简单的同相、恒耦合辐射积分，也应沿整个磁长积分密度调制：

\[
\mathcal A(L)\propto\int_0^L b(z)\,dz.
\]

在常系数纯密度模式中

\[
\boxed{\frac{W_{\mathrm{SC}}}{W_0}
=\left[\frac{\sin(k_{\mathrm{SC}}L)}{k_{\mathrm{SC}}L}\right]^2},
\]

不是出口 \(|b(L)/b_0|^2\)。本例前者约为 \(0.94483\)。

若把空间电荷恢复项放回冷束单模 FEL 方程，在同一归一化中可写为

\[
\mathcal A'=g_Fb,\qquad
b'=-ig_Fp,\qquad
p'=-g_F\mathcal A-i\frac{k_{\mathrm{SC}}^2}{g_F}b,
\]

从而

\[
\boxed{\mathcal A'''+k_{\mathrm{SC}}^2\mathcal A'
=ig_F^3\mathcal A},\qquad
\boxed{\mu^3+k_{\mathrm{SC}}^2\mu=ig_F^3}.
\]

这给出空间电荷与 FEL 反馈共同存在时的冷束单模色散关系；真实有限脉冲的能量还依赖沿程复电流、衍射和多模干涉。

### 8 近似条件与结论含义

空间电荷推导采用小调制、冷束、近共振周期平均、准静态屏蔽势、给定横向高斯分布和线性包络。通常远长于共振波长推导的波荡器阻抗，不能直接当作本例共振波长附近的精确阻抗；这里的横向核是量级模型，没有包含完整波荡器 Green 张量。

纯密度单模结果表明，纵向空间电荷可以自然产生约 \(5\%\) 至 \(8\%\) 的聚束抑制尺度。实测局部 RMS BF 是随机多模统计量，不等于这个相干模式传递；公式支持空间电荷作为解释候选，不能单凭数值接近完成唯一归因。

## 参考文献

- [Ming Xie PAC 1995 原文](https://proceedings.jacow.org/p95/ARTICLES/TPG/TPG10.PDF)：Pierce 参数约定与三维拟合公式。
- [Huang 与 Kim 2007 FEL 理论综述](https://journals.aps.org/prab/pdf/10.1103/PhysRevSTAB.10.034801)：线性 FEL 增长和增益修正。
- [Wang Litvinenko 与 Webb 2012](https://journals.aps.org/prab/pdf/10.1103/PhysRevSTAB.15.120701)：周期平均空间电荷增长尺度。
- [Geloni 等 DESY 07 087](https://arxiv.org/abs/0706.2280)：波荡器纵向阻抗及长波近似的适用条件。

空间电荷的有限椭圆高斯横向投影及本例参数代入为上述近似下的模型推导，不是文献完整三维场解的直接转录。
