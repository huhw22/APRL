"""Create the portable explanation and PPT handoff, without changing old deliveries."""
import csv
import json
import re
import shutil
import sys
import zipfile
from pathlib import Path
from prepare import sha,save

def main(root):
    root=Path(root).resolve();r=json.loads((root/'RESULT.json').read_text());t=r['tail_mask']
    text=f'''# Puffin 边界处理与场和能量差异说明

2026-11-17 · 面向 PPT 制作的独立交付

## 主要结论

Puffin 采用前向、傍轴的辐射模型，场衍射通过傅里叶变换推进，边缘另用掩膜和谱滤波处理。它没有使用 APRL 的 CPML。Puffin 边界影响应理解为有限域中的场形、相位、模式和能量可能受到改变；不能把所有场差都称为物理反射，更不能移用 APRL 的反射率。[1][2]

本次确认了一个具体实现问题：纵向掩膜的总长度表达式是 `nx - 1*dx`，而坐标定义给出的几何长度是 `(nx-1)*dx`。直接调用现有生产库得到的曲线确实发生错位。当前十一段参数下，掩膜第一节点权重为 **{t['native_first']:.9f}**，而仅修正长度括号的几何参考在该点接近零；两曲线最大绝对差为 **{t['max_native_vs_parentheses_difference']:.6f}**。[3]

这一差异是掩膜权重的差，不是反射率。它对已完成六段或十一段 FEL 结果的影响尚未量测，既有结果也没有因此自动作废。当前仍没有 Puffin 专用的入射/返回通量反射率实测。

## 1 Puffin 的边界在做什么

Puffin 官方模型忽略沿 −z 传播的后向波，采用傍轴近似。它可以出现横向模式改变、有限域回卷或边缘处理污染；这些现象与完整 Maxwell 模型中沿 −z 返回的反射波须分开描述。[1]

当前十一段采用非周期时间域设置 `meshType=0`。这与“单波长周期 FEL 模式”不同。但衍射算子使用有限数组的离散傅里叶变换，其数值表示具有周期性，需由足够的缓冲范围和吸收处理降低边缘影响。不能仅凭用了 FFT 就断言实际中心场已被回卷污染。[2]

当前横向有 81 × 81 节点，间隔 40 µm，节点范围为 ±1600 µm；对应 FFT 周期为 NΔx = 3240 µm。每侧 16 个边界节点，掩膜在约 |x|、|y| > 1000 µm 的区域起作用。已有 FEL 功率比较孔径是 ±520 µm，位于掩膜直接作用区域以内。

![生产库横向掩膜](figures/zh-CN/01_transverse_mask.png)

图中 M 是进入边缘滤波分支的场权重，中央约为零、边缘趋近一。它不是每步损失的能量比例，也不是反射率。

## 2 实际操作并非简单把边缘场乘一个衰减系数

在一次衍射传播之后，将原场记为 A，实际边界操作可写成

\\[
M=M_\\perp+M_{{z2}}(1-M_\\perp),\\qquad
A_{{out}}=(1-M)A+\\mathcal F^{{-1}}\\!\\left[D\\,\\mathcal F(MA)\\right].
\\]

实现中的谱系数为

\\[
D(\\bar k_x,\\bar k_y,\\bar k_{{z2}})=
\\exp\\!\\left[-h\\beta\\frac{{|\\bar k_x|+|\\bar k_y|}}{{\\sqrt{{|2\\bar k_{{z2}}|}}}}\\right],
\\qquad \\bar k_{{z2}}\\ne0.
\\]

横线表示 Puffin 的无量纲坐标和波数；h为无量纲衍射步长。这里的β是输入项 `sBeta`，是吸收参数，与束流速度或Lorentz变换中的β无关。波数为零的分支在这一步不乘该指数。代码中的 `ffact` 是 FFT 正反变换的尺度因子，不是可调物理吸收率；相关注释不能替代实际表达式。[2]

因此，边缘场先被取出，再在傅里叶空间按波数和频率滤波，最后与保留分支相加。改变边缘模式可能同时改变场形和相位。即使某个内部节点 M=0，滤波分支的逆变换仍是空间非局域操作，不能据此保证内部节点完全不受影响。

同样，D² 是某个边缘分支傅里叶分量的平方范数系数，不能直接当作整束能量保留率。有限束乘上 M 后会产生额外谱分量，再与另一支相干叠加。

以当前 β=1、每四分之一波荡器周期的一次名义衍射步（约7 mm）及参考辐射频率代入，得到：

| 物理 kx（1/m），ky=0 | 场系数 D | 该分量平方范数系数 D² |
|---:|---:|---:|
'''
    for m in r['modal_factors']:
        text+=f"| {m['kx_per_m']} | {m['D_amplitude']:.9f} | {m['D2_norm2']:.9f} |\n"
    text+='''
![边缘分支谱系数](figures/zh-CN/03_modal_factors.png)

这是一组代入实际参数的算子系数，不表示 FEL 主模式必然有这些角度或损失。低频过滤 `qFilter` 是另一个操作，不应与这里的边缘吸收合成一个反射指标。

## 3 场差和能量差应如何分别定义

与一个相同输入、相同传播模型的大域参考场比较，令 δA=A−Aref，可以分别定义

\[
\\epsilon_A=\\frac{\\|\\delta A\\|_2}{\\|A_{ref}\\|_2},\\qquad
\\epsilon_U=\\frac{\\|A\\|_2^2-\\|A_{ref}\\|_2^2}{\\|A_{ref}\\|_2^2}.
\]

两者满足

\[
\\epsilon_U=\\frac{2\\operatorname{Re}\\langle A_{ref},\\delta A\\rangle+\\|\\delta A\\|_2^2}{\\|A_{ref}\\|_2^2}.
\]

其中存在干涉项，所以一般没有“能量误差等于场误差平方”的关系。即使能量积分一致，复场也可能存在显著相位或空间分布差异。

例如 A=e^{iφ}Aref 时，能量完全相同，而场差为 2|sin(φ/2)|。φ=60° 时，复场相对差可达100%。这是严格数学示例，不是 Puffin 已测得60°相位误差。

![相位示例](figures/zh-CN/07_phase_example.png)

只有在同阻抗、匹配模式，并且入射与返回信号已分离的条件下，才可用独立振幅反射系数的模平方计算对应能量反射率。对一般两程序场差、总场峰值变化或边缘算子误差，都不能直接这样转换。

## 4 图像和视频展示的是一次边界响应

为显示边缘处理如何同时改变场和平方范数，采用独立规定的单频复场：横向强度高斯标准差80 µm，kx=20000 1/m；将其中心位置从0扫到1400 µm。每个输入只经过一次上式的横向边界操作，不含纵向掩膜、衍射传播、粒子源或FEL放大。每个比例均以该输入自己的范数为分母。

| 输入中心 x0（µm） | 复场 L2 变化 | 平方范数净变化 |
|---:|---:|---:|
'''
    for q in r['prescribed_cases']:
        text+=f"| {q['x0_um']:.0f} | {q['complex_field_L2_difference_percent']:.6g}% | {q['energy_proxy_change_percent']:+.6f}% |\n"
    text+='''
![不同规定输入位置的一次响应](figures/zh-CN/06_field_vs_energy.png)

[中文视频](figures/zh-CN/prescribed_boundary_response.mp4) · [English video](figures/en/prescribed_boundary_response.mp4)

视频两栏分别显示规定输入复场与一次操作后的实部，灰色区域表示横向掩膜区域。所有帧使用固定场轴范围；标题中的 x0 是人为规定的输入位置，不是时间或传播距离。每帧独立计算，不能将视频解读为波包在 Puffin 中撞上边界后的真实传播。

此处平方范数在固定单频、同一网格下作为场强平方代理。它不是完整实验室系能量账本，也没有测出返回方向的独立能量通量。边缘处范数降低符合吸收操作的目的，本身不构成数值错误的证明。

## 5 一个已确认的纵向掩膜实现问题

纵向坐标由 xj=jΔx 构造，末节点位置应为 (N−1)Δx；源码中却写成

```text
Lf = (nx-1*dx)
x0 = Lf - Lb/2
```

按 Fortran 运算优先级，第一行是 nx−dx。当前 N=96000、Δx约0.003141592741，得到的 Lf约95999.996858，而几何长度约301.589762。这个问题同时可在官方公开源码中看到；本交付的数值以本地冻结生产库的函数输出为准。[3]

![纵向长度表达式问题](figures/zh-CN/04_longitudinal_length_issue.png)

末端掩膜包含37节点、36间隔，物理宽度约0.884394 µm。当前原生曲线在起点从内部零权重跳到约0.665177。改变总节点数20个，即仅增减一个参考波长的网格范围，也会改变同样末端距离处的曲线。

![改变一个波长的网格长度](figures/zh-CN/05_longitudinal_grid_sensitivity.png)

蓝色几何参考仅把表达式换为 (nx−1)×dx，其余表达式保持。这不是已经验收的完整边界修复。尤其该参考曲线本身也不是 CPML，不能仅修正一处括号就宣布无反射。

掩膜权重差不能直接换算为场误差或能量反射率。还需知道实际场是否到达纵向末端、其频谱和横向波数，以及后续传播如何把改变传回比较区域。当前交付未据此解释六段或十一段的增益、饱和位置和能量差，也没有改动其结果或生产程序。

## 6 关于反射，目前可以和不可以得出的结论

可以确认：Puffin 的边界算法与 CPML不同；其边缘处理依赖波数和频率；场差和能量差应分别测量；纵向掩膜存在上述长度表达式问题。

尚不能确认：生产工况究竟有多少反射、回卷或中心场污染；该实现问题贡献了多少 FEL 增益或功率差；单次边界操作的变化能否直接代表多步自洽计算的误差。

若以后量测，应以同一规定输入、相同网格间隔和传播模型，比较现有小域、充分大域参考及仅修改边界的版本。在同一内部区域分别统计复场差、场强分布差、能流/能量差；若称横向反射，还须分离向内返回的波矢或通量。这个后续控制尚未执行。

前一条 APRL 测试中的0.000461是CPML反射电场峰值比；其能量代理约2.66×10⁻⁷。两个数字都不能作为 Puffin 的反射率放进PPT。

## 来源与数据

1. [Puffin 官方模型说明](https://ukfels.github.io/puffinDocs/)：前向、傍轴以及所忽略的物理项。
2. [官方 diffraction.f90](https://github.com/UKFELs/Puffin/blob/master/puffin/lib/diffraction/diffraction.f90)；本地实际版本的[相关行及SHA](evidence/diffraction.f90.txt)。主公式逐行对应现有冻结实现。
3. [官方 masks.f90](https://github.com/UKFELs/Puffin/blob/master/puffin/lib/diffraction/masks.f90)；本地实际版本的[相关行及SHA](evidence/masks.f90.txt)。掩膜CSV由现有编译库直接输出。

计算结果见 [RESULT.json](RESULT.json)。原生横向掩膜、纵向掩膜、系数表及规定输入扫描在 [数据索引](data/README.md)；[独立验证](VALIDATION.json)保留数值检查与适用范围。PPT的页面安排、素材路径和可直接转交的指令见 [PPT_HANDOFF.md](PPT_HANDOFF.md)。
'''
    (root/'README.md').write_text(text)
    handoff='''# 给 PPT 制作任务的指令

请只根据本包制作 Puffin 边界处理的说明页，不重新仿真，不把 APRL 的 CPML 反射率标成 Puffin 反射率。建议5页，中文讲述，图片亦有英文同源版本。所有公式和数字的主依据为 README.md 与 RESULT.json。

## 第1页 Puffin 的边界处理方式

讲述：Puffin采用前向、傍轴辐射模型；FFT衍射后，通过掩膜取出边缘场再进行谱滤波。它与CPML的处理不同。用 `figures/zh-CN/01_transverse_mask.png` 和 `02_mask_profile.png`，选择其中一张作为主图，另一张作补充。M只表示分支权重，不是反射率。

## 第2页 边缘操作如何改变复场

主公式：Aout=(1−M)A+F⁻¹[D F(MA)]。用 `03_modal_factors.png`。举例kx=20000 1/m时，一次约7mm步的边缘支路D=0.959379、D²=0.920409。必须保留“边缘支路系数”这个限定，不得写成整束损失7.96%。

可播放 `figures/zh-CN/prescribed_boundary_response.mp4`，时长9.4秒；或并排使用 `08_prescribed_input.png` 与 `09_single_boundary_output.png`，两张图使用同一绝对色标。视频的输入位置扫描不是时间演化，必须保持画面中的原限定文字。不要把它改名为“Puffin反射实测”。

## 第3页 场差大不代表能量损失同样大

用 `06_field_vs_energy.png`，x0=1280µm处一次操作复场L2差约2.5913%，平方范数净变化约−4.4818%；两者定义不同，不能直接平方换算。可用 `07_phase_example.png` 补充严格数学例子：仅相位差60°时场差100%、能量相同。明确它不是Puffin实测相位误差。

主公式：ΔU=2 Re〈Aref,δA〉+∥δA∥²。解释干涉项，区分独立反射波与一般场误差。

## 第4页 一个具体的纵向掩膜问题

用 `04_longitudinal_length_issue.png`。将原表达式 nx−1×dx 与几何参考 (nx−1)×dx 写清。当前生产库直接输出的起始权重约0.665177，几何参考接近零。最大权重差约0.812905，不是81.29%反射率。不要把几何参考标为已经验收的修复。

若需要显示对计算域长度的敏感性，使用 `05_longitudinal_grid_sensitivity.png`。三条曲线只是同一函数在总节点数增减20时的数值，不是三次FEL运行。

## 第5页 结论的范围

已确认的是算法形式与纵向掩膜表达式问题；尚未测得生产工况的反射率，也未证实它造成六段/十一段的功率差。因此用“需要单独量测边界误差”作为结论，不写“Puffin有X%反射”“六段差异已归因边界”或“十一段结果作废”。

完整Puffin反射率量测需要同输入的大域参考，以及向内返回模式/通量的分离。现有PPT只说明已得到的证据，不补造这类实验结果。

## 包内文件的使用

- `figures/zh-CN` 与 `figures/en`：中英文同源PNG/PDF与视频；PDF适合矢量插图。
- `data`：原生掩膜、谱系数和规定输入的一次响应表。
- `evidence`：实际源码相关行、生产参数和来源SHA。
- `VALIDATION.json`、`VISUAL_REVIEW.json`、`SHA256SUMS`：交付检查记录。

只交付PPT制作素材与说明；本包不是PPT本体，也不是新FEL仿真的完整科学数据包。
'''
    (root/'PPT_HANDOFF.md').write_text(handoff)
    (root/'data/README.md').write_text('''# 数据说明

- `native_transverse_mask.csv`：现有编译库 getMask 输出，81×81节点。
- `native_tail_masks.csv`：现有编译库 getZ2Mask 输出，三种N各37节点。
- `tail_mask_comparison.csv`：原生值及只更改几何长度括号的参考曲线。
- `modal_factors.csv`：边缘支路D与D²，不是整束能量或反射率。
- `prescribed_position_scan.csv`：141个独立规定输入的一次横向边界响应，x0不是时间或传播距离。

`operator_illustration.npz`在本地分析目录保留完整复场数组；为减少PPT交付体积，不放进ZIP。复场公式与所有用于图表的数值记录均在包内。没有生产科学H5或可执行程序。
''')
    shutil.copyfile(Path('/public/home/wfollowm/APRL_Workspace/results/puffin_eleven_four_slip_20261113_r2/source/LICENSE.txt'),root/'evidence/PUFFIN_LICENSE.txt')
    (root/'SUMMARY_EN.md').write_text('''# Puffin boundary treatment and interpretation of field errors

Puffin uses a forward, paraxial radiation model with FFT diffraction and a masked spectral boundary operator. The boundary weight M selects a branch; it is not a reflection coefficient. The implemented operation is Aout=(1−M)A+F⁻¹[D F(MA)]. Field changes, net squared-norm changes, and separately measured reflected energy are distinct quantities.

The existing production library confirms a longitudinal length-expression issue: `nx-1*dx` differs from the coordinate-based length `(nx-1)*dx`. At Nz2=96000 the first boundary-mask value is about0.665177, rather than near zero for the parentheses-only geometric reference. This reference is not an accepted fix. The production FEL impact has not been measured.

Figures and video show source formulas, native mask-function values and a single boundary operation on independent prescribed inputs. The video is a position scan, not a Puffin propagation run. No Puffin reflectivity measurement or explanation of the six/eleven-segment power difference is claimed. APRL CPML reflection measurements must not be attributed to Puffin.
''')
    # No source, executable, production H5 or large internal operator arrays in the handoff ZIP.
    include=[]
    for p in root.rglob('*'):
        if not p.is_file():continue
        rel=p.relative_to(root)
        if rel.parts[0]=='verification' or p.suffix=='.npz':continue
        if p.name in ['SHA256SUMS','PACKAGE.json']:continue
        include.append(p)
    if not (root/'VALIDATION.json').exists() or not (root/'VISUAL_REVIEW.json').exists():
        print('Documents ready; packaging awaits calculation and visual review records.');return
    manifest={str(p.relative_to(root)):sha(p) for p in sorted(include)}
    (root/'SHA256SUMS').write_text(''.join(v+'  '+k+'\n' for k,v in manifest.items()))
    target=root.with_suffix('.zip')
    with zipfile.ZipFile(target,'x',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for p in include+[root/'SHA256SUMS']:z.write(p,root.name+'/'+str(p.relative_to(root)))
    with zipfile.ZipFile(target) as z:
        assert z.testzip() is None
        for k,v in manifest.items():
            import hashlib
            assert hashlib.sha256(z.read(root.name+'/'+k)).hexdigest()==v
    info=dict(complete=True,zip=str(target),bytes=target.stat().st_size,sha256=sha(target),
         members=len(include)+1,member_sha256_verified=True,crc_verified=True,
         production_impact_or_reflectivity_measured=False)
    save(root/'PACKAGE.json',info)
    save(root.parent/'LATEST_PUFFIN_BOUNDARY_TRANSFER.json',info)
    print(json.dumps(info,ensure_ascii=False,indent=2))

if __name__=='__main__': main(sys.argv[1])
