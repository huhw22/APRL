# 参考场的横向等效表面电流（实验选项）

`mesh.moving_window.reference_surface_current: true` 默认关闭，与
`reference_through_cpml` 互斥。此次用户允许调整参考场在横向边界的处理；
不调整内部Cowan/Boris/粒子电流沉积和原CPML系数/更新。

原参考势在横向CPML内缘有限Dirichlet边界外置零，界面的法向E、切向B
具有跳跃，而参考Ez为零。原Ampere更新将B的跳跃变为Ez启动脉冲。
在四个界面的原生Jz位置增加预设电流

    Js,z = [(By(i,j,k)-By(i-1,j,k))/dx
             -(Bx(i,j,k)-Bx(i,j-1,k))/dy] / mu0

其中B只读取已准备的参考场。沿z按既有translationBeta平移到下一半步，
再进行线性插值。不会读取实际总B来构造电流，不钳制实际E/B。
因此额外传播扰动仍受原线性Maxwell/CPML算子作用；这不意味着有限
参考场本身已被证明是无限空间的精确物理解，需要与独立空间电荷理论比较。

新增边界数值源的做功 `-dt*dV*sum(w*Js*(Ez_old+Ez_new)/2)` 并入
既有 `boundary_source_work_J`，并另写逐rank CSV。原生场能权重在面为1/2、
角为1/4。账本闭合仅验证记账，不能单独当成动力学正确性证据。

准备阶段流式写表面电流，约O(Nz*(Nx+Ny))数据；运行时两个32层缓存，
不保留整个参考场。每步仅访问活动窗的横向表面。关闭选项没有新文件、
额外通信或电流循环。当前测试口径为Cowan-z、无电流滤波、仅横向CPML。
无电流滤波是本原型的使用条件；有滤波时还需处理源的滤波和做功口径，
不能把本轮结果外推过去。

测试顺序：77required及关闭100步科学输出逐位；uniform500步开启冲量探针
和严格账本；原核首步脉冲因果检查；改善后才运行三类完整短束团及双节点。
原数据、A方案失败结果、默认程序及cluster.env保持。

原理背景：表面场跳跃可用等效电/磁表面电流表示，见
https://www.ee.cityu.edu.hk/~amhwong/papers/%28AMHW%20Ang%202019%29%20Equivalene%20Priciple%20Based%20Active%20Metasurfaces%20EMTS.pdf
式(1)。当前参考势在面上为零，切向参考E为零，因而本原型只补电流项。
上述连续理论不自动证明本离散实现严格等价，首步与独立束团理论对照才是
本轮验证依据。源是Maxwell中的J项、并不硬指定总E/B的区别也见
https://meep.readthedocs.io/en/stable/Python_User_Interface/#source 。
