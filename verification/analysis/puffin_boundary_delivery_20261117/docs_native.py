"""PPT documentation for native absorption versus the large-domain reference."""
import json
from pathlib import Path
from common import table

ASSETS=[('01_propagation_reference','大域参考中主束离开共同视野'),
 ('02_propagation_native','原生吸收削弱主束并留下弱返回'),
 ('03_absorption_amount','终点全域积分场强剩余0.2961%'),
 ('04_central_full_history','原生吸收与参考的中央场强完整变化'),
 ('05_native_return_zoom','后期中央返回峰0.4569%'),
 ('06_interior_field_feedback','固定入射尺度的中央复场差6.76%'),
 ('07_reverse_x_modes','原生吸收例负kx模式份额'),
 ('08_native_return_xy','返回峰处的原生横向场'),
 ('09_reference_xy','同位置同色标的大域参考'),
 ('10_central_error_xy','中央误差复场的平方')]

def main(root):
    root=Path(root);s=json.loads((root/'SUMMARY.json').read_text());n=s['small-native']
    (root/'README.md').write_text(f'''# Puffin 原生吸收与大域参考对照

2026-11-17 · 实际原生传播 · PPT 交付

## 主要结果

同一倾斜光束在原生吸收小域和大域参考中传播。原生吸收明显削弱到达横向边缘的主束，但仍有弱返回进入中央观察区。

在 z=2.1 m 时，原生小域的全域积分场强只剩入射值的 **{100*n['final_total_norm2_over_initial']:.4f}%**，净移除 **{100*n['net_norm2_removed_fraction']:.4f}%**。在较早的 z=1.505 m，中央返回峰值为入射全域量的 **{100*n['late_max_central_norm2_over_initial']:.4f}%**，对应中央复场范数差约 **{100*n['max_central_field_L2_change_over_initial']:.2f}%**。两个位置的统计不能相减。

## 1 先看主束如何被吸收

[中文线性视频](figures/zh-CN/native_propagation_linear.mp4) · [English](figures/en/native_propagation_linear.mp4)

左栏是大域参考的共同视野，右栏是原生吸收。两栏使用同一输入、相同传播距离及固定色标。开始时光束一致，随后向右到达边缘；右栏的主束明显变暗。左栏光束离开视野后，仍存在于更大的计算域中，所以其全域保留量接近100%。

中央方框为 x/y=±520 µm 观察区，原生小域中的竖虚线标记 x方向边缘处理区。每帧增加的是真实传播距离，不是人为移动输入位置。

![大域参考传播](figures/zh-CN/01_propagation_reference.png)

![原生吸收传播](figures/zh-CN/02_propagation_native.png)

两张图横轴为 x、纵轴为传播距离 z，每行均为一次原生推进的 y=0 强度剖面。大域参考亮带离开后，中央应接近空；原生小域中则有较弱的返回带经过中央。

## 2 再看弱返回如何影响内部

[中文弱场视频](figures/zh-CN/native_propagation_log.mp4) · [English](figures/en/native_propagation_log.mp4)

线性色标适合看主束吸收，较弱残留接近黑色；对数版使用另一套固定色标，使弱返回可见。两栏始终共用尺度，没有逐帧或按各栏峰值重新缩放。

![中央场强完整变化](figures/zh-CN/04_central_full_history.png)

原光束离开后，大域中央场降到接近零。原生小域的后期弱返回在下图放大显示。

![中央弱返回](figures/zh-CN/05_native_return_zoom.png)

在61个保存的快照中，后期中央峰值出现在 z=1.505 m，约为入射全域积分场强的0.4569%；同位置大域参考约为入射量的1.16e-14。此时原生小域全域保留约0.619%。终点 z=2.1 m 的全域剩余0.296%对应更晚的位置，不与这个中央峰值矛盾。

## 3 全域吸收量与内部复场差分别统计

![全域吸收量](figures/zh-CN/03_absorption_amount.png)

全域曲线说明吸收确实移除了绝大多数场；但判断内部是否接近正确参考，还要比较中央区域的复场。

![内部复场差](figures/zh-CN/06_interior_field_feedback.png)

逐点减去同一 z 的大域复场后，中央场差的最大 L2 范数为入射全域场 L2 的6.76%。后期参考中央已接近零，因而0.457%的场强平方量可以对应约6.76%的场范数；二者描述不同量。

分母始终是固定的入射全域场范数，不是后期近零的局部参考。6.76%不是每个点的相对误差，也不是 FEL 增益误差。快照每3.5 cm保存一次，这些峰值是保存样本中的最大值。

![原生返回横向分布](figures/zh-CN/08_native_return_xy.png)

![同位置的大域参考](figures/zh-CN/09_reference_xy.png)

两张图在 z=1.505 m 使用同一色标和相同中央方框。[中央误差场图](figures/zh-CN/10_central_error_xy.png)显示两复场之差的平方，不是反射能量，也不是两幅强度图直接相减。

## 4 返回场与反向模式分开表述

![横向反向模式份额](figures/zh-CN/07_reverse_x_modes.png)

原生吸收例的负 kx 模式平方范数份额最大约5.11e-8，分母为入射全域平方范数。它很小，说明不能把全部中央场差直接当作大量传播方向反转的反射。Puffin不包含完整后向 z 波，这个横向模式指标也不是沿 −z 的物理反射率。

## 测试条件与适用范围

本包采用19834已闭合的原生吸收、大域及更大域参考。小域81×81、横向间隔40 µm、每侧16节点边缘处理；参考网格401×241及481×321。光束强度标准差200 µm，起始 x=−300 µm，倾角约2.346 mrad；传播300步、2.1 m，无粒子及 FEL 放大。

使用原生单波长周期纵向分支隔离横向边界，因此没有测试纵向末端掩膜。此前纵向表达式问题留在[源码背景](background/README.md)，没有参与本次中央返回结果。

积分场强是同频前向模型中的场强平方积分，可作为功率/能量流代理；不是完整 Maxwell 能量账本。当前受控测试不能直接换算为六段/十一段 FEL 的功率、增益或粒子反馈误差。

19834正常COMPLETED 0:0。两个大域参考在中央的差约6.94e-13，大域与独立高斯参考的差约6.94e-13，均按固定入射尺度统计。

[PPT 制作指令](PPT_HANDOFF.md) · [素材索引](PPT_ASSET_MAP.csv) · [结果摘要](SUMMARY.json) · [数据说明](data/README.md) · [验证记录](VALIDATION.json) · [目视记录](VISUAL_REVIEW.json)
''')
    (root/'PPT_HANDOFF.md').write_text('''# PPT 制作指令

本包围绕原生吸收与大域参考制作4页，主图和视频都使用双栏对照。以实际传播和内部弱残留为主，不需要在主讲中推导公式。

## 第1页 主束到达边缘后被吸收

播放 figures/zh-CN/native_propagation_linear.mp4：左栏为大域共同视野，右栏为原生吸收。两栏同一输入、同一z和固定色标。说明光束向右传播，原生边界使其明显变暗。参考视野变暗只表示光束离开这一块区域，大域全域中的光仍然存在。

静态备选为01_propagation_reference和02_propagation_native。各图横轴x、纵轴传播距离z；保持相同色标。

## 第2页 吸收很强 仍需检查内部弱返回

主图03_absorption_amount：终点z=2.1m只剩0.2961%入射积分场强，净移除99.7039%。再播放native_propagation_log.mp4，用固定对数色标看线性色标下很弱的返回，保留所有色标和范围说明。

## 第3页 中央残留与参考的差别

用04_central_full_history和05_native_return_zoom。中央±520µm在z=1.505m的返回峰值为入射全域量0.4569%；同位置大域参考接近零。原生小域在这一位置的全域保留约0.619%。不要把较晚终点的0.296%与中央峰值相减。

并列08_native_return_xy与09_reference_xy，同一z、同一色标和相同观察方框。10_central_error_xy是复场差的平方，不称作反射能量。

## 第4页 内部场差及结论范围

用06_interior_field_feedback：中央复场L2差最高约6.76%，分母是固定入射全域L2，不是每点相对误差或FEL增益误差。07_reverse_x_modes可作补充：负kx模式最大5.11e-8，不能把所有场差统称物理反射。

保留范围：有意让倾斜光束到达横向边界、无粒子、没有FEL放大；纵向使用单波长周期控制，因此此前纵向掩膜问题没有参与本次残留。不能把本测试直接作为生产FEL误差修正。

## 文件与画面

中文素材在figures/zh-CN，英文同源在figures/en；PNG/PDF同源，PDF适合插图。四段视频为中英文各线性/对数一段，双栏1280×760、61个实际快照、5fps、12.2秒。每帧同一尺度，没有移动输入或幅度拟合。

读取根目录README、SUMMARY与data中的保留统计即可；按PPT_ASSET_MAP选材。背景文件是源码补充，不替代实际传播。此次交付为素材和说明，不是PPT本体。
''')
    (root/'SUMMARY_EN.md').write_text('''# Native Puffin absorption versus a wide-domain reference

The two-column movies compare the same incident field at the same propagation distance with one fixed scale. Watch the linear version for removal of the main beam, then the log version for weak returning fields.

At z=2.1m the native case retains0.2961% of the incident full-domain integrated intensity, a net removal of99.7039%. At the earlier z=1.505m the central ±520µm aperture has a late returning peak of0.4569% of the incident squared norm. Its complex-field difference reaches6.76% of the initial full-domain field L2 norm. These are different measurements and propagation positions.

The native negative-kx modal squared norm peaks at5.11e-8 of input. This is not a backwards-z reflectivity. The test is particle-free, deliberately intercepts a transverse boundary and does not measure production FEL gain error or longitudinal end-mask effects.
''')
    (root/'data/README.md').write_text('''# 数据口径

comparison.csv保留原生吸收、大域和更大域参考各61个快照，共183行；history.csv保留对应逐步统计，共903行。相同case、相同位置的数值逐字段保持。

visual_arrays.npz只包含坐标、reference/native共同复场快照和两者的y=0逐步剖面。数组是原数据的直接选择，未改变数值、相位或幅度。RESULT.json与SUMMARY.json保留相关结果及参考验证；元数据中的delivery_projection标明内容选择。

中央积分使用±520µm节点梯形面积权重。复场差按初始全域场L2定标；所有比值的入射分母固定。负kx统计来自完整81²原生小域，不对裁剪的参考做方向统计。
''')
    rows=[]
    for lang in ['zh-CN','en']:
        for stem,meaning in ASSETS:
            rows.append(dict(language=lang,asset=f'figures/{lang}/{stem}.png',vector=f'figures/{lang}/{stem}.pdf',meaning=meaning,source='19834 retained native/reference fields'))
        for mode in ['linear','log']:
            rows.append(dict(language=lang,asset=f'figures/{lang}/native_propagation_{mode}.mp4',vector='',meaning='主束吸收' if mode=='linear' else '弱返回和内部影响',source='61 retained actual snapshots; fixed common scale'))
    table(root/'PPT_ASSET_MAP.csv',rows)
