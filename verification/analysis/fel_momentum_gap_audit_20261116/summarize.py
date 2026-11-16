"""Summarize the existing-data audit; preserve earlier reports and raw data."""
import hashlib
import json
import re
import subprocess
import sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager, ticker

LD = np.longdouble


def read(p):
    return json.loads(Path(p).read_text())


def save(p, value):
    Path(p).write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + '\n')


def sha(p):
    h = hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda: f.read(4194304), b''):
            h.update(b)
    return h.hexdigest()


def main(out, font):
    out, font = Path(out).resolve(), Path(font).resolve()
    campaign = out.parent
    previous = campaign / 'momentum-ledger-20261116-r1'
    result = read(out / 'RESULT.json')
    old = read(previous / 'RESULT.json')
    independent = read(out / 'INDEPENDENT_PORT_CAP_CHECK.json')
    side = np.genfromtxt(out / 'INDEPENDENT_SIDES.csv', delimiter=',', names=True, dtype=LD)
    actual = old['fixed_lab_worldtube']['final_components_J']
    side_checks = []
    for short, key in zip(['xminus_J', 'xplus_J', 'yminus_J', 'yplus_J'],
                         ['transverse_xminus_signed_J', 'transverse_xplus_signed_J',
                          'transverse_yminus_signed_J', 'transverse_yplus_signed_J']):
        error = LD(side[short]) - LD(actual[key])
        assert abs(error) < LD('1e-18')
        side_checks.append(dict(channel=short, independent_J=float(side[short]),
                                previous_J=actual[key], difference_J=float(error)))
    save(out / 'INDEPENDENT_SIDE_CHECK.json', dict(passed=True,
        algorithm='Piecewise affine vertical bounds with exact two-point Gaussian integration; no polygon moments',
        interior_cells=int(side['interior_cells']), boundary_cells=int(side['boundary_cells']),
        maximum_absolute_difference_J=max(abs(v['difference_J']) for v in side_checks), channels=side_checks))

    a = np.load(previous / 'boost/cache.npz')['records']
    G, beta = LD(a['boost_gamma'][0]), LD(a['boost_beta'][0])
    with np.load(previous / 'whole/energy-momentum/whole_cache.npz') as z:
        cum = lambda key: np.r_[LD(0), np.cumsum(z[key], dtype=LD)]
        ge, gb = G * beta * cum('fixed_end_GaussE_J'), G * beta * cum('fixed_end_GaussB_J')
        ce, cp = G * cum('fixed_end_C_H_J'), G * beta * cum('fixed_end_C_P_J')
    pairing = np.genfromtxt(out / 'energy_momentum_pairing.csv', delimiter=',', names=True, dtype=LD)
    longitudinal_with_charge = pairing['longitudinal_pairing_after_static_external_cancel_J'] - ge
    detail = dict(
        longitudinal_pairing_with_native_charge_force_J=float(longitudinal_with_charge[-1]),
        native_charge_electric_impulse_Gamma_beta_J=float(ge[-1]),
        native_magnetic_divergence_Gamma_beta_J=float(gb[-1]),
        Cowan_energy_Gamma_J=float(ce[-1]), Cowan_longitudinal_momentum_Gamma_beta_J=float(cp[-1]),
        caution='GaussE is the measured native charge-electric force; GaussB and Cowan are explicit native operator terms. These use full boost-time slices, not the fixed laboratory worldtube.')
    save(out / 'NATIVE_FORCE_AND_OPERATOR_DETAIL.json', detail)

    # Audit the frozen source actually associated with the reused binary,
    # rather than silently treating the current development tree as its code.
    workspace = campaign.parents[1]
    frozen_root = workspace / 'results/sparc_six_full_ledger_20261112_r3'
    frozen = read(frozen_root / 'FROZEN.json')
    names = ['source_coupling.cpp', 'source_coupling_math.h', 'longitudinal_shape.h',
             'simulation_source_coupling.cpp', 'whole_cowan_budget.h', 'yee_sbp_geometry.h',
             'paired_momentum_geometry.h', 'lab_worldtube_audit.h', 'boostframe.cpp',
             'eb_sources.cpp', 'particle_energy_audit.h']
    source_hashes = {}
    for name in names:
        p = frozen_root / 'source/src' / name
        assert sha(p) == frozen[str(p)], name
        source_hashes[str(p)] = sha(p)
    coupling = (frozen_root / 'source/src/source_coupling.cpp').read_text()
    assert 'if(!linearNotch&&!linearBinomial3)return openOrbitMagnetic' in coupling
    convergence = (frozen_root / 'source/src/simulation_source_coupling.cpp').read_text()
    assert 'force<1e-13&&work<1e-8' in convergence
    save(out / 'SOURCE_REVIEW.json', dict(frozen_source_SHA256_verified=source_hashes,
        actual_shape='quintic-compensated', particle_coupling='conservative-source',
        E_gather_and_J_deposition_use_matching_shaped_orbit_weights=True,
        quintic_B_gather_uses_native_trilinear_midpoint=True,
        convergence_checks_force_fixed_point_and_energy_work=True,
        strict_total_longitudinal_momentum_constraint_not_implemented=True,
        exact_total_momentum_conservation_not_established_by_method_name=True,
        qualification='The previous momentum-matching label denotes the compensated longitudinal-shape option, not a proven simultaneous global energy/momentum invariant.'))

    font_manager.fontManager.addfont(str(font))
    family = font_manager.FontProperties(fname=str(font)).get_name()
    plt.rcParams.update({'font.family':[family,'DejaVu Sans'], 'font.size':12,
        'axes.titlesize':16,'axes.spines.top':False,'axes.spines.right':False,
        'text.parse_math':False,'axes.unicode_minus':False,'pdf.fonttype':42})
    x = a['time_box_s'] * LD('1e12')
    plot_records = []
    for language in ['zh-CN', 'en']:
        zh = language == 'zh-CN'
        directory = out / 'figures' / language
        directory.mkdir(parents=True, exist_ok=True)
        for name, curves, title, unit in [
            ('energy_pairing', [(pairing['energy_pairing_after_static_external_cancel_J']*1e9,
              '电功配对差' if zh else 'Electric-work pairing gap')],
             '五次补偿方法：扣除静态外磁场后的能量配对' if zh else 'Compensated quintic: energy pairing after static-magnet cancellation',
             '能量（nJ）' if zh else 'Energy (nJ)'),
            ('longitudinal_pairing', [
              (pairing['longitudinal_pairing_after_static_external_cancel_J']*1e6,
               '粒子与电流纵向配对差' if zh else 'Particle and current longitudinal pairing gap'),
              (longitudinal_with_charge*1e6,
               '再计入原生电荷电力后的配对差' if zh else 'Gap including native charge-electric force')],
             '五次补偿方法：纵向冲量配对差' if zh else 'Compensated quintic: longitudinal impulse pairing gap',
             '投影后的能量量纲（µJ）' if zh else 'Projected energy units (µJ)')]:
            fig, ax = plt.subplots(figsize=(10,6), layout='constrained')
            for i, (y, label) in enumerate(curves):
                line, = ax.plot(np.asarray(x,float),np.asarray(y,float),label=label,
                                color=['#2679ae','#d37b28'][i],linestyle=['-','--'][i])
                assert np.array_equal(line.get_ydata(),np.asarray(y,float))
            ax.axhline(0,color='#888888',lw=.7);ax.grid(alpha=.2)
            ax.set(title=title,xlabel='计算系模拟时间（ps）' if zh else 'Computational-frame time (ps)',ylabel=unit)
            ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y,pos:f'{y:g}'))
            ax.legend(frameon=False)
            for suffix in ['png','pdf']:
                fig.savefig(directory/(name+'.'+suffix),dpi=180)
            plt.close(fig)
            plot_records.append(dict(language=language,name=name,axes=1))
    save(out / 'PLOTS.json', dict(figures=plot_records, no_subscripts_or_superscripts=True))

    text = '''# 六段五次补偿方法：实验室系差额定位

本轮仅检查19796_1已闭合数据与其冻结代码；没有推进粒子或Maxwell场，没有修改原程序、卡或科学输出。只读后处理19809和19811均正常结束。

## 当前判断

**已排除后处理漏项、加和符号及端口/侧面/时间端面积分算法造成这5.675265µJ差额。当前代码没有严格的粒子加场总纵向动量守恒保证，现有全局配对账本也有非零差额。证据更指向场—粒子耦合及离散场定义，而不是简单统计错误；尚不能把37.56%全部归为局部推进器错误。**

此前“动量匹配”这个简称过强。真实开关是`longitudinal_shape: quintic-compensated`，两种方法共用conservative-source平均速度推进器。这一分支改了纵向电荷/电流形函数与伴随电场取场，并不等于总纵向动量严格守恒。旧图/原记录保留，这里纠正其物理含义。

## 为什么理论推论成立，当前实现却不能直接套用

同一组离散定义、同一时空控制域和一致源项下，若总能量与纵向动量都严格闭合，则变换后的残差满足：

    R实验室 = Gamma [ R能量 + beta R纵向cP ]

因此两项同时为零，实验室系能量也闭合。这里仅验证了计算系的能量配对；没有满足另一项的严格全局配对条件。此外，原生Cowan离散场能与连续场u、S的重构定义不同，不能把两者混成同一个守恒对象。

实际冻结代码显示：

- `source_coupling.cpp`第29–36行：电场使用与电流沉积相同的高阶轨迹权重；这保证电功的伴随配对。
- 同文件第38–39行：五次补偿分支的磁场仍进入`openOrbitMagnetic`，是原生三线性中点取样，没有改成高阶电流沉积在场动量公式下的伴随磁场取样。
- `source_coupling_math.h`第239行起的`momentum()`求解局部平均速度Lorentz方程；函数名中的momentum不是“粒子加场总动量守恒约束”。
- `simulation_source_coupling.cpp`第42–45行：外层检查场/轨道固定点与电功缺陷，未检查粒子纵向冲量是否与离散场反冲量相等。
- `paired_momentum_geometry.h`中的`gaussInteriorElectricImpulse`和`momentum(dj,B)`分别给出原生电荷电力与电流磁力配对；它们使用特定Yee平均。高阶Ez能量取样、原生B中点取样并不自动满足这套冲量伴随关系。

这些代码事实说明不能仅凭“高阶核+能量守恒”推出完整三维/纵向动量守恒。也不能仅将B核换高阶就声称问题已解决：电荷电力、J叉B、Cowan和边界配对都需共同验证。

## 现有数据把交换差缩小到纵向冲量

上一轮在计算系等时面上的实验室系组合差为+3.907792766µJ，其中未扣外场的能量部分约+370.380021424µJ，纵向动量部分约−366.472228657µJ。不能把这两个大数直接当成自场耦合误差。

本卡只有静态实验室磁场、没有入射波或外电场。每个外场样本来自同一Lorentz变换，平均速度同时用于力与位移，所以在离散方程的代数层面有：

    外场电功prime + beta 外场纵向cP冲量prime = 0

据此消去外磁场的对应功/冲量，得到：

| 量 | 末值 µJ |
|---|---:|
| 已记录外场boost功乘Gamma | +370.380203429 |
| 扣外场后的电功配对差 | −0.000182006 |
| 扣外场后的粒子/电流纵向配对差 | +3.907974772 |
| 两者加和 | +3.907792766 |
| 再计入原生Gauss电荷电力后的纵向配对差 | +3.934780242 |

电功差还独立用`Gamma × (particle_grid_work − 累计WJ)`重算，得到−0.000182006µJ，与上表一致。粒子能量的Lorentz关系误差约1.1e-18J，也不能解释微焦耳级差额。

**上述外场冲量通过静态磁场关系作代数消元，不是新测得的逐粒子外场冲量。** 原计算没有保存完整子步电力/磁力轨迹，因此不能再把3.934780µJ精确拆成Ez取场、v叉B取场和外磁推进各占多少。当前已定位到纵向场—粒子配对层面，尚未唯一定位到某一条局部推进语句。

这个全局差额还可能含注入/载流子的未单列纵向约束冲量，不能全部称为活跃粒子的推进误差；37.56%的固定lab管已经排除源支持所在区域，两个范围仍须区分。

![电功配对](figures/zh-CN/energy_pairing.png)

![纵向冲量配对](figures/zh-CN/longitudinal_pairing.png)

实际不存在“全束已注入且无人退出”的共同窗口：首次退出在第10225步，全束注入完成到第16938步。因此不伪造一个无边界对照。已记录差额在首次退出前累积0.977216µJ，注入/退出重叠期累积2.203782µJ，全注入后仍累积0.726794µJ。它不是只在最后一个退出事件突然出现；边界的独立力分量仍不能从现有记录全部分离。

## 独立检查统计：原37.56%仍成立

原实验室管为z=0.03→14.94625m、横向±840µm、完整公共迟滞时间门、全频总场。粒子净输运−15.108648µJ，场端口+14.009310µJ，四侧+4.744129µJ，场时间端面差+2.030474µJ，加和+5.675265µJ，即37.5630%。完整粒子覆盖与空粒子端盖的检查保持。

本轮用不同算法复算全部已记录场表面：

- 轴向端口和时间端面：直接作逆Lorentz直线、按时间/空间格点切段，用Simpson公式积分每段二次函数，没有复用原多边形矩算法。21项最大差8.21e-21J。
- 四侧：独立用仿射上下界与二点Gauss积分，处理243638721个内部格和77704个边界格；不使用原多边形裁切/矩公式。四项均与原数值一致，具体差值见JSON。
- 已保存的早期与晚期无粒子区域单独作全表面收支：末面差分别+0.000208390µJ与−0.000695894µJ。无粒子判据来自完整过面事件及亚光速轨道上tau单调增加。这是一个不同域的场统计检查，不是完整差额的分解，也不是所有有粒子区域的误差上界。

独立核对证明后处理积分几何没有制造5.675265µJ。它只验证已保存的电流型观测量积分，不证明从原生Yee/Cowan场重构连续u、S没有离散差。

## 仍未分清的场定义差异

原生场动量账本保留了Cowan及Gauss项。本轮末值：Gamma×Cowan能量项+11.983387µJ，Gamma×beta×Cowan动量项−1.880926µJ，Gamma×beta×GaussE为−0.026805µJ，GaussB为−0.733823µJ。这些项本来就不是零，不能把其全局算子恒等式闭合当成连续Poynting守恒已经成立。

worldtube观测器先将错位E/B平均到共同位置与时刻，再形成平方和乘积；原生离散能量/动量使用其原生位置上的乘积及配对时间层。“先平均后平方”和“原生平方后平均”之间有方差/相关项。代码中的等时B平均、光锥端盖u−S纵/c和Lorentz变换符号检查没有发现错误，但现有二维面与六条纵向观测量不保留重建任意固定lab管原生Cowan修正所需的完整局部三维信息。

因此不能把5.675265−3.907793=1.767472µJ称作统计误差：两者的时空域、孔径、场定义与源项通道不同，尚无同域的严格桥接。37.56%保留为连续场观测账本的实测残差；当前证据尚不支持将其全部定量归入局部粒子推进。

## 下一步定位边界

现有数据检查已完成。若随后修复，应先按实际提交子步记录粒子自场Ez冲量、v叉B冲量和静态外磁场冲量，与相同控制域的原生场反冲量分别比较；同时验证连续场重构到原生Cowan账本的同域修正。不能把缺少的力项用残差补回，也不能用账本补偿强行抹零。本轮没有提交这些新物理测试。

## 证据

- [纵向/电功配对表](energy_momentum_pairing.csv)
- [现有账本分段结果](PAIRING.json)
- [原生力与算子项](NATIVE_FORCE_AND_OPERATOR_DETAIL.json)
- [独立端口与端盖检查](INDEPENDENT_PORT_CAP_CHECK.json)
- [独立四侧检查](INDEPENDENT_SIDE_CHECK.json)
- [无粒子场区检查](PARTICLE_FREE_FIELD_CHECK.json)
- [实际过面的时间范围](PARTICLE_TIME_BOUNDS.json)
- [冻结源码核对](SOURCE_REVIEW.json)
- [原完整账本报告](../../../fel_six_comparison_20261116_r3/README.md)
'''
    # The prior report sits directly under results, two directories up.
    text = text.replace('../../../fel_six_comparison_', '../../fel_six_comparison_')
    (out / 'README.md').write_text(text)
    for target in re.findall(r'!?\[[^\]]*\]\(([^)]+)\)', text):
        assert (out / target).exists(), target
    for path, before in read(out / 'PROTECTED_OUTPUTS.json').items():
        s = Path(path).stat()
        assert before == dict(size=s.st_size,mtime_ns=s.st_mtime_ns,inode=s.st_ino)
    for path, expected in read(out / 'FROZEN_ANALYSIS.json').items():
        assert sha(path) == expected
    save(out / 'VALIDATION.json', dict(passed=True, no_new_physics=True,
        original_output_identities_unchanged=len(read(out/'PROTECTED_OUTPUTS.json')),
        frozen_runtime_source_files_verified=len(source_hashes),
        independent_ports_and_caps_passed=True, independent_four_sides_passed=True,
        native_energy_and_longitudinal_momentum_pairing_separated=True,
        exact_full_worldtube_residual_causal_decomposition_available=False,
        actual_particle_force_substeps_missing=True, PNG_visual_review_pending=True,
        PDF_separately_rendered=False))
    print('AUDIT_REPORT_READY',out,flush=True)


if __name__ == '__main__':
    main(*sys.argv[1:])
