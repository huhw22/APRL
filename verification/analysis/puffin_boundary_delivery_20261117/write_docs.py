"""Explain the measured propagation process without requiring equations."""
import json
import sys
from pathlib import Path
from common import table

def main(root):
    if (Path(root)/'SELECTION.json').is_file():
        from docs_native import main as write_selected
        return write_selected(root)
    root=Path(root);s=json.loads((root/'SUMMARY.json').read_text());n=s['small-native'];off=s['small-no-absorption']
    text=f'''# Puffin 边界吸收与计算域内返回场

2026-11-17 · 实际原生传播结果 · PPT 主交付修订

## 主要结果

这次直接推进同一束倾斜高斯光束，让它到达横向边界。**关闭吸收时，光束从右侧离开后又从左侧绕回；原生吸收削弱了主束，但仍有弱场返回中央区域。** 因而“总量几乎守恒”和“边界已经吸收绝大多数场”，都不足以单独说明内部光场正确。

原生吸收在 z=2.1 m 时净移除了 **{100*n['net_norm2_removed_fraction']:.4f}%** 的全域积分场强。但在较早的 **z={n['central_return_peak_z_m']:.3f} m**，中央区域的返回场强积分仍达入射全域值的 **{100*n['late_max_central_norm2_over_initial']:.4f}%**，对应中央复场范数差约 **{100*n['max_central_field_L2_change_over_initial']:.2f}%**。此时大域参考的中央场已经接近零。

这些是有意让光束撞到横向边界的无粒子测试结果；不是六段或十一段 FEL 的误差估计。

## 1 先看线性视频 主束如何被吸收

[播放中文线性视频](figures/zh-CN/native_propagation_linear.mp4) · [English](figures/en/native_propagation_linear.mp4)

三栏按同一传播距离同步播放：左栏为大域参考的同一块视野，中栏为原生吸收边界，右栏为关闭吸收的控制。先看到三束相同的光向右移动；接近右边缘后，中栏主束明显变暗。右栏则在右侧离开以后，从左侧重新出现。

每帧推进的是真实传播距离 z，输入中心没有被人为重新放置。全部帧和三栏共用同一强度色标。左栏只显示与小域相同的区域；光束离开该视野后，大域全域中的光仍然存在，所以它的“全域保留”数字仍接近100%。

视频中的中央方框为 x/y=±520 µm 观察区；小域中的虚线表示边缘处理区位置。下方分别显示全域保留量与中央积分量，二者都除以同一个入射全域值。

## 2 再看对数视频 弱返回场如何影响内部

[播放中文弱场视频](figures/zh-CN/native_propagation_log.mp4) · [English](figures/en/native_propagation_log.mp4)

线性色标适合看主束被吸收，但较弱的返回会接近黑色。对数版使用另一个固定色标，专门使弱返回可见；它仍然没有逐帧放大或把各栏峰值拉到一样高。

在后期，左栏参考视野应接近空。右栏的强返回穿过中央方框，中栏则有弱得多的返回。先前看不清的“对计算域内部的影响”，就是这些本不应存在于该区域的场，以及它们相对大域参考的复场差。

## 3 空间传播图把离开与返回连起来

下面每一行都是一次实际原生推进后的 y=0 强度剖面。横轴是 x，纵轴向上增加传播距离 z。三张图使用同一固定对数色标。

![大域中光束离开视野](figures/zh-CN/01_propagation_reference.png)

大域图中，亮带向右上方延伸并离开视野，随后中央区域变暗。

![关闭吸收后从另一边绕回](figures/zh-CN/03_propagation_no_absorption.png)

关闭吸收时，亮带从右侧离开，又从左侧进入，之后继续向右上方传播。这是有限 FFT 数组的周期绕回：它没有因到达数组边缘而消失，也没有主要改成向左传播。

![原生吸收后的弱返回](figures/zh-CN/02_propagation_native.png)

原生吸收显著压低返回场，但对数图仍能看到一条较弱的返回带穿过中央区域。光束具有有限宽度，边缘和中央区域会先接触到它的尾部；图中的几何中心位置不能当作所有影响首次发生的精确时刻。

独立的周期像高斯解析叠加也复现了关闭吸收时的传播。与原生场的最大范数差约4.37e-6，且几乎全部在初始有限高斯尾部截取时已存在，沿程变化仅约2.2e-15。这支持将强返回识别为周期绕回，而不是新产生的能量。

## 4 吸收掉多少 与中央还剩多少分开看

![全域积分场强变化](figures/zh-CN/04_absorption_amount.png)

大域参考和关闭吸收例的全域积分场强都几乎保持不变。原生吸收例在 z=2.1 m 时只剩入射量的 **{100*n['final_total_norm2_over_initial']:.4f}%**。这是边界确实在移除场的直接证据。

![中央区域完整变化](figures/zh-CN/05_central_full_history.png)

中央区域开始都能看到原来的入射光束。它向右离开后，大域参考降到接近零；关闭吸收例随后重新升高，后期峰值为入射全域值的 **{100*off['late_max_central_norm2_over_initial']:.3f}%**。全域总量守恒并没有阻止错误的场重新进入观测区。

![原生吸收后的中央弱返回放大](figures/zh-CN/06_native_return_zoom.png)

这张图只放大 z≥1.3 m 的后期，显示原生吸收例的弱返回。它在 z={n['central_return_peak_z_m']:.3f} m 达到 **{100*n['late_max_central_norm2_over_initial']:.4f}%**，而同一时刻大域参考仅约入射量的 **{n['reference_central_at_return_peak']:.2e}**。

“终点全域剩余0.296%”和“中央返回峰值0.457%”对应不同传播位置，不能把二者相减或要求后者小于前者。在中央返回峰值位置，原生小域的全域保留量约为0.619%。

## 5 对内部复场的影响可以达到几个百分点

![原生吸收例的内部复场差](figures/zh-CN/07_interior_field_feedback.png)

这张曲线先逐点减去相同 z 的大域复场，再统计中央区域的差，分母固定为入射全域场的 L2 范数。原生吸收例的最大值为 **{100*n['max_central_field_L2_change_over_initial']:.2f}%**；关闭吸收例最大约 **{100*off['max_central_field_L2_change_over_initial']:.2f}%**。

原生吸收例在保存的 z=1.05 m 样本中，中央场差首次达到这一入射尺度的1%以上；前一个 z=1.015 m 样本尚未达到。场快照每3.5 cm保存一次，不能将该记录位置当作连续过程的精确阈值交点。

中央返回峰值时，大域参考已接近零，因而0.457%的返回场强平方量可以对应约6.76%的场范数。两个数字描述不同量，不矛盾。这里的6.76%也不是每个网格点的相对误差，更不是 FEL 增益误差。

![原生弱返回场在横向上的分布](figures/zh-CN/09_native_return_xy.png)

![同一位置的参考场](figures/zh-CN/10_reference_xy.png)

两张图在同一 z 和同一固定色标下比较，中央方框相同。另附[中央误差场图](figures/zh-CN/11_central_error_xy.png)，展示两复场之差的平方；它不是两幅强度图直接相减，也不是反射能量。

## 6 绕回与反射如何区分

![横向反向模式](figures/zh-CN/08_reverse_x_modes.png)

还可以检查传播方向。关闭吸收例的负 kx 模式份额始终约4e-12，与初始值相同；但它的中央返回非常强。原生吸收例的负 kx 份额最大也只有约 **{n['negative_kx_max_fraction']:.3g}**，相对于入射全域平方范数。

因此，本测试看到的主要返回是周期绕回，不是大量向 −x 传播的反射。这个横向模式权重也不是沿 −z 的物理反射率：Puffin 的模型本来不包含完整的后向 z 波。PPT 应使用“周期绕回与弱返回场”，不要把全部内部场差统称为反射能量。

## 7 测试条件与结论范围

同一无粒子倾斜高斯光束，强度标准差200 µm，起始中心 x=−300 µm，倾角约2.346 mrad；横向间隔40 µm，传播300步、共2.1 m。小域81×81节点，每侧16节点处理区，名义每步7 mm；另有401×241大域和481×321更大域参考。

为隔离横向边界，本控制使用原生单波长周期纵向分支、每波长20个间隔，关闭纵向末端掩膜。没有粒子、FEL放大或新的边界修复。先前发现的纵向掩膜表达式问题仍在[源码背景附件](background/README.md)，但本测试没有测它的影响。

本文“积分场强”指横向场强平方的积分，在这个同频前向模型中作为功率/能量流代理。所有比例固定使用共同的入射全域值；没有按每幅图或每一时刻重新归一化。这里没有完整 Maxwell 能量账本，也没有量测一个独立的有限脉冲后向反射能量。

该受控结果说明：当光束到达边界时，吸收可以移除绝大多数场，同时仍留下可影响内部复场的弱返回。它不能直接换算为现有六段/十一段 FEL 的功率或增益误差，也没有包含粒子对返回场的响应。

## 来源和交付使用

19834正常 COMPLETED 0:0，总用时4分44秒。原生平面波预检、四个完整传播例均闭合；大域与更大域在中央的场差约6.94e-13，大域与独立傍轴高斯参考的差约6.94e-13，均使用固定入射尺度。关闭吸收例的全域范数沿程漂移约8.76e-14。

- [PPT 制作指令](PPT_HANDOFF.md)：按“离开—吸收—绕回—内部影响”的顺序排版。
- [PPT 素材索引](PPT_ASSET_MAP.csv)：图像、视频、数据和每张图能支持的结论。
- [精确结果摘要](SUMMARY.json) · [原统计](data/RESULT.json) · [数据说明](data/README.md)。
- [验证记录](VALIDATION.json)与[目视记录](VISUAL_REVIEW.json)保留检查范围。

旧版规定位置扫描视频未重复收入新版主素材，旧 ZIP 保持。源码背景附件只是补充；新的PPT应以本页的实际传播证据为主。
'''
    (root/'README.md').write_text(text)
    handoff='''# 转交给 PPT 制作任务的指令

请用本包替换旧 Puffin 边界说明的主素材，以实际传播过程说明吸收、周期绕回和内部光场变化。不要把旧版规定位置扫描视频当作传播视频。建议5页，主讲不需要数学推导；公式和纵向掩膜问题只放背景附件。

## 第1页 一束光到达边界以后会发生什么

先播放 `figures/zh-CN/native_propagation_linear.mp4`。三栏为同一输入、同一z、同一色标：大域参考的同一视野、81²原生吸收、81²关闭吸收。指出主束向右离开、原生边界使它变暗、关闭吸收后从左侧重新出现。不要把参考视野变暗说成大域能量消失。所有帧来自实际原生推进。

静态备选：`01_propagation_reference`、`02_propagation_native`、`03_propagation_no_absorption`。单独使用各PNG/PDF，统一色标。横轴x、纵轴传播距离z；图中的返回带是关键。

## 第2页 吸收确实有效 但全域总量不足以判断内部场

主图 `04_absorption_amount.png`：原生吸收在z=2.1m后只剩入射积分场强的0.2961%，净移除99.7039%。关闭吸收几乎保持100%，却出现绕回。不要把这张图叫作完整Maxwell能量账本。

可并列 `05_central_full_history.png`：没有吸收时，中央原光束离开后又升到入射全域量的73.218%。这是错误位置的场返回，不是能量凭空增加。

## 第3页 弱残留对内部区域仍有影响

播放 `native_propagation_log.mp4`。这个版本使用固定对数色标，看弱返回；没有按每帧或每栏重新缩放。主图用 `06_native_return_zoom.png`：后期中央返回峰值0.4569%，z=1.505m，大域参考在这里几乎为零。

用 `09_native_return_xy.png` 和 `10_reference_xy.png` 显示同一平面的差别，保留共同色标和中央方框。`11_central_error_xy.png` 是复场差的平方，不要称作两幅强度直接相减或反射能量。

## 第4页 小平方量仍可产生几个百分点的场差

主图 `07_interior_field_feedback.png`：原生吸收例中央复场L2差最大6.76%，分母为初始全域场L2。不能写成各点误差6.76%、FEL增益差6.76%或独立反射幅度系数。

0.457%是z=1.505m的中央积分量；0.296%是z=2.1m的全域剩余量，时刻不同，不能互相相减。在中央峰值时，全域剩余约0.619%。

## 第5页 主因是周期绕回 不应把全部返回叫作物理反射

用 `08_reverse_x_modes.png`：关闭吸收时反向横向模式几乎没有增长，却有很强绕回；原生吸收例的负kx份额最大约5.11e-8。配合空间图说明主要返回仍保持向+x传播。这个指标不是沿−z物理反射率。

页脚保留范围：本测试故意让一束约2.346mrad倾斜光到达横向边界，无粒子，纵向使用单波长周期控制。它证明本测试的边界残留与内部场影响，尚未测出六段/十一段FEL误差，也没测试纵向掩膜问题。不能把99.7%净移除改写为“完全无反馈”。

## 文件和画面规范

中文素材在 `figures/zh-CN`，英文同源素材在 `figures/en`；PNG与PDF同源，PDF适合矢量插图。线性/对数两个视频都为61个实际快照、12.2秒、1680×760。保持色标、单位和限定文字；如需缩写，将完整指标定义放在讲稿中。

主文 README.md 已按观看过程解释；精确值在 SUMMARY.json，CSV/NPZ在 data。制作时以新版实际传播证据为主，不把 background 中的单次算子例子拿来替代。不要新增仿真、填造反射率或声称已修复生产程序。
'''
    (root/'PPT_HANDOFF.md').write_text(handoff)
    (root/'SUMMARY_EN.md').write_text('''# Puffin absorption, wrap-around and interior-field contamination

The source-free native propagation test completed successfully as job19834. With absorption disabled, the tilted beam exits the right edge and wraps back through the left edge. Its late central-aperture integrated intensity reaches73.218% of the common incident full-domain value, although the whole-domain norm is almost conserved.

Native absorption removes99.7039% of the integrated intensity by z=2.1m. A weak returning field nevertheless reaches the central ±520µm aperture: the late peak is0.4569% of the incident full-domain squared norm at z=1.505m. The corresponding central complex-field difference reaches6.76% of the initial full-domain field L2 norm. These statistics refer to different quantities; the terminal remaining global value0.2961% and the earlier central peak must not be subtracted.

Reverse-x mode weight is very small, peaking at5.11e-8 of input in the native case. Spatial trajectories and independent periodic-image propagation support wrap-around as the dominant returning mechanism, rather than a large reversal of transverse direction. This is not a backwards-z reflectivity measurement.

Watch the fixed linear-scale movie first to see absorption, then the fixed log-scale movie to see weak returns. All61 snapshots are from actual native propagation; input positions are not manually scanned, and there is no per-frame amplitude rescaling. The wide reference is displayed through the same spatial crop as the small cases.

This intentionally boundary-intercepting test contains no particles or FEL gain. It isolates transverse boundaries using a native one-wavelength periodic longitudinal mode, and does not test the longitudinal end-mask expression. Its results cannot be assigned as six- or eleven-segment production FEL errors. Use PPT_HANDOFF.md and the matching English figures/videos to prepare slides.
''')
    (root/'data/README.md').write_text('''# 数值和数组口径

- comparison.csv：四例61个共同快照的统计，使用共同初始全域积分量作分母。中央积分为±520µm节点上的梯形面积权重。
- history.csv：四例各301步的全域及中央积分。
- RESULT.json：19834原始只读分析结果，字节保持。
- visual_arrays.npz：真实传播的81²共同复场快照、每步y=0强度剖面和原坐标；reference/native/no_absorption分别对应大域、原生吸收、关闭吸收。
- periodic_image_check.csv：独立周期像高斯解析叠加与无吸收原生场的差；这是解释验证，不是另一组新传播仿真。

“范数平方”“积分场强”与前向模型同频功率代理成正比；不是全频实验室系能量账本。roi_field_L2_difference的分母是初始全域场L2，不是后期近零的局部参考场。负kx份额仅在完整81²小域傅里叶基底中计算，未把裁剪的大域场做方向统计。

所有比值分母固定，无时刻拟合、相位对齐或幅度匹配。原始逐步二进制和可执行程序不放入PPT交付包。
''')
    rows=[]
    desc=[('01_propagation_reference','大域中光束向右离开同一视野'),('02_propagation_native','原生吸收削弱主束，仍有弱返回'),
       ('03_propagation_no_absorption','周期绕回从左侧进入，仍向右传播'),('04_absorption_amount','终点净移除99.7039%积分场强'),
       ('05_central_full_history','无吸收中央后期返回73.218%'),('06_native_return_zoom','原生吸收的弱中央返回0.4569%'),
       ('07_interior_field_feedback','固定入射尺度的中央复场差6.76%'),('08_reverse_x_modes','负kx模式很小，区别方向反转与绕回'),
       ('09_native_return_xy','峰值时原生返回横向分布'),('10_reference_xy','同z、同色标的大域参考'),('11_central_error_xy','中央复场差的平方，不是反射能量')]
    for lang in ['zh-CN','en']:
        for stem,meaning in desc:rows.append(dict(language=lang,asset=f'figures/{lang}/{stem}.png',vector=f'figures/{lang}/{stem}.pdf',meaning=meaning,source='19834 native closed fields'))
        for mode in ['linear','log']:
            rows.append(dict(language=lang,asset=f'figures/{lang}/native_propagation_{mode}.mp4',vector='',meaning='主束吸收' if mode=='linear' else '弱返回和内部影响',source='61 native snapshots; same fixed scale'))
    table(root/'PPT_ASSET_MAP.csv',rows)

if __name__=='__main__':main(sys.argv[1])
