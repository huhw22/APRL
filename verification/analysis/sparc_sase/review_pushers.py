"""Interpret finite-pulse gain and independently check the two completed pushers."""
import csv
import hashlib
import json
import sys
from pathlib import Path

import numpy as np


def read(p):
    return json.loads(Path(p).read_text())


def write(p, value):
    Path(p).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')


def theory(top, review):
    reference=top/'explicit_review_20261104/radiation_phase_r1'
    gain=read(reference/'GAIN_AUDIT.json')['cases']['after']
    rows=[]
    for j in range(5):
        values=dict(interval=f'{j+1}-{j+2}',full_bunch_band_Lg_m=gain['adjacent']['rectangular_band_J'][j]['power_gain_length_magnet_m'])
        for name,key in [('peak_slice','fixed_initial_peak_current_slice'),('central95','fixed_central_current_ge_0.95_of_input_peak')]:
            ratio=gain['local_power_exact_segment_exit_fits'][key]['adjacent_ratio'][j]
            values[name+'_Lg_m']=float((77*.028)/np.log(ratio))
        rows.append(values)
    write(review/'THEORY_GAIN.json',dict(rows=rows,theory_approximate_Lg_m=.67,
      metric='Adjacent exact segment exits; cumulative magnetic length; power or band energy e-folding',
      finite_pulse_theory_acceptance=False))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size':12,'axes.spines.top':False,'axes.spines.right':False,'lines.linewidth':1.8})
    fig,ax=plt.subplots(figsize=(9.0,4.5),layout='constrained')
    for key,label,color,marker in [('full_bunch_band_Lg_m','Full-pulse band energy','#24648c','o'),
      ('central95_Lg_m','Fixed central region: input I >= 95% peak','#8e4b75','^'),
      ('peak_slice_Lg_m','Fixed input peak-current slice','#d47928','s')]:
        ax.plot(np.arange(5),[r[key] for r in rows],marker+'-',color=color,label=label)
    ax.axhline(.67,color='#666666',ls='--',label='Literature local estimate: 0.67 m')
    ax.set_xticks(np.arange(5),[r['interval'] for r in rows])
    ax.set_xlabel('Adjacent undulator segments');ax.set_ylabel('Effective power gain length (m)')
    ax.set_title('GENESIS: gain inferred from exact segment exits');ax.grid(alpha=.18);ax.legend()
    fig.savefig(review/'theory_gain_diagnostic.png',dpi=210);plt.close(fig)
    (review/'THEORY_GAP.md').write_text(r'''# GENESIS与理论的差距

目前没有完成整短束的理论验收。修正相位后，末第4–6磁段：整束矩形频带能量有效增益长度0.827 m、固定入口中央高流强区域（I≥95%峰值）平均功率0.782 m；比文献0.67 m分别长约23.5%、16.7%。单个入口峰流强切片0.668 m仅占约0.108%电荷，不能用它替代前两种验收口径。

文献0.67 m是匹配束流的近似局部基模功率增益尺度。Huang–Kim式79–81对失谐做优化，给出Ming Xie近似及其适用范围；它与对应数值本征模解通常仍有10%以内拟合差，因此不能把此估计当作整短束绝对输出功率的5%标准。

$$P(s)\propto e^{s/L_g},\qquad L_{g,0}=\frac{\lambda_u}{4\pi\sqrt3\rho},\qquad L_g^{j,j+1}=\frac{\Delta s}{\ln(W_{j+1}/W_j)}.$$

这里s只计累计有磁长度。当前已经排除间隙，不能再次用间隙坐标换算消掉余差。图中每个比值均由确切段出口计算；频带矩形能量与固定区域的全包络功率是不同指标，分别标明。

![各段有效增益长度](theory_gain_diagnostic.png)

整束相邻段的有效增益长度约1.424、1.015、0.901、0.850、0.806 m，仍逐段缩短。固定峰值单片第2–5段拟合约1.095 m，末第4–6段才接近0.668 m。因此目前观测到明显的启动及模式选择过渡，不能假定六段中已形成稳定的整脉冲渐近增长平台。

中央95%区域入口流强只有约50.627–53.275 A。仅用rho随I的三分之一次方变化，冷一维Lg最多增加约1.7%，不足以解释16.7%差距。中央区域长约152 μm，全程磁段名义滑移约227 μm，光场跨越中央区域，固定切片并非始终处在独立且恒定的局部环境。有限束长及滑移影响确实需要考虑，但目前尚未定量证明它们解释全部余差。

本次保留峰值流强，将原文12.67 ps、450 pC束团缩短为2.271 ps、80.675 pC；因此也不是原实验整脉冲功率曲线的直接复现。当前可以用两种PIC推进在同一输入下的差异判断算法敏感性，仍不能凭GENESIS绝对值判定哪一种已满足5%物理精度。

另需分清功率量级：此前约17–23 MW估计是局部饱和峰功率。六磁段尚未饱和，当前GENESIS出口原生瞬时峰值约2.75 MW；同一频带的Hann记录平均值约0.458 MW。瞬时峰值、加窗记录平均值及饱和理论值是三种不同量，不能直接按其差额判程序错误。这里只澄清口径，没有把当前绝对功率视为理论验收通过。

来源：[Campbell等，FEL2019 TUP050](https://proceedings.jacow.org/fel2019/papers/tup050.pdf)；[Huang与Kim，2007，式6及79–81](https://journals.aps.org/prab/pdf/10.1103/PhysRevSTAB.10.034801)。逐项旧新拟合见../explicit_review_20261104/radiation_phase_r1/GAIN_AUDIT.json。上述原因为有限束影响的证据与待检验解释，没有给余差做残差反推或功率拟合。
''')


def review_result(top, review):
    out=review/'radiation_comparison_r1'
    assert read(out/'COMPLETE.json')['complete']
    result=read(out/'COMPARISON.json');old=Path(result['reference_directory'])
    explicit=read(old/'pic.json');energy=read(out/'pic.json');gen=read(old/'genesis.json')
    checks=[]
    band=read(old/'RADIATION_RESULT.json')['band_eV']
    for a,b,g in zip(explicit,energy,gen):
        name=b['plane'];T=b['fft_samples']*b['dt_s']
        with np.load(old/(name+'-spectra.npz')) as ac,np.load(out/(name+'-spectra.npz')) as bc,np.load(old/(name+'-genesis.npz')) as gc:
            for window in ['rect','hann']:
                values=[];peaks=[]
                for cache,key in [(ac,window+'_native_flux_J_per_eV'),(bc,window+'_native_flux_J_per_eV'),(gc,window+'_J_per_eV')]:
                    x=cache['energy_eV'];m=(x>=band[0])&(x<=band[1]);y=cache[key]
                    values.append(float(y[m].sum()*(x[1]-x[0])))
                    peaks.append(float(x[np.flatnonzero(m)[np.argmax(y[m])]]))
                row=next(r for r in result['rows'] if r['plane']==name and r['window']==window)
                for v,key in zip(values,['explicit_W','energy_coupled_W','GENESIS_W']):
                    assert abs(v/T/row[key]-1)<1e-11
                assert abs(values[1]/values[0]/row['energy_over_explicit']-1)<1e-11
                checks.append(dict(plane=name,window=window,explicit_W=values[0]/T,
                  energy_W=values[1]/T,GENESIS_W=values[2]/T,
                  energy_vs_explicit_pct=100*(values[1]/values[0]-1),
                  energy_over_GENESIS=values[1]/values[2],native_spectral_peaks_eV=peaks))
    for manifest in [review/'FROZEN.json',top/'compact-production-r1/FROZEN.json']:
        assert all(hashlib.sha256(Path(p).read_bytes()).hexdigest()==h for p,h in read(manifest).items())
    write(out/'VALIDATION.json',dict(arithmetic_validation_passed=True,physical_acceptance=False,
      checks=checks,images_visually_inspected=False,frozen_inputs_unchanged=True,
      note='Spectral L1 includes amplitude difference; peak columns in comparison.csv use E-only proxy, native peaks independently reported here'))
    far=[r for r in checks if r['plane']=='far'];a,b=far
    maxdelta=max(abs(r['energy_vs_explicit_pct']) for r in checks)
    text=f'''# 能量守恒耦合版与显式版：SPARC短束辐射对照

本工况中，两种推进的共振频带读数差异显著：最终出口Hann记录平均功率由显式版{b['explicit_W']/1e6:.3f} MW变为能量耦合版{b['energy_W']/1e6:.3f} MW，相差{b['energy_W']/b['explicit_W']:.3f}倍。不能沿用“只改善能量记账、辐射基本不变”的假设；这也不直接证明哪一版更接近真实物理。

19294_1能量耦合版正常完成，23984步、八个场面全部闭合，与19294_0显式版的各面时间向量和横向网格完全相同。原相位修正GENESIS参考直接复用，本轮不重跑物理仿真。GENESIS与理论的余差见[单独说明](../THEORY_GAP.md)，仍不把它当作5%绝对真值。

|15.246 m，同频带与孔径|显式版|能量耦合版|GENESIS|能量版相对显式版|
|---|---:|---:|---:|---:|
|矩形窗记录平均功率|{a['explicit_W']/1000:.3f} kW|{a['energy_W']/1000:.3f} kW|{a['GENESIS_W']/1000:.3f} kW|{a['energy_vs_explicit_pct']:+.3f}%|
|能量归一Hann记录平均功率|{b['explicit_W']/1000:.3f} kW|{b['energy_W']/1000:.3f} kW|{b['GENESIS_W']/1000:.3f} kW|{b['energy_vs_explicit_pct']:+.3f}%|

两种推进的八面、两种窗口中，最大频带功率相对变化为{maxdelta:.3f}%。这是同一噪声实现下的算法敏感性；是否满足5%需要同时查看功率、谱形和位置依赖，不能只选最终一个面。完整数值见comparison.csv。谱差L1含幅度改变，不是各谱独立归一化后的形状误差。

从磁段出口到其后0.15/0.30 m，Hann下能量版/显式版比值约6.518/6.499/6.467，差异在三个面上保持，不能仅用矩形窗记录边缘解释。最终Hann原生能流谱峰从约2.50196 eV移到2.52195 eV；名义共振为2.52351 eV，GENESIS为2.51520 eV。峰位更接近名义共振仍不能代替功率和增长机制的验收。

![四宫格能谱](spectra_comparison.png)

![功率随位置变化](radiation_power.png)

频带{band[0]:.6f}–{band[1]:.6f} eV；孔径2.4×2.4 mm²；同一3.930766 ps记录。平均功率是频带能量/记录时长，不是瞬时峰值。PIC是原始总E/B，未扣束缚场；GENESIS是辐射包络。同一Hann窗不等于完成近场分离。CSV峰位继承横向E²代理定义，VALIDATION.json另外从实际E×B谱独立取峰。

两卡除输出目录外，唯一物理方法差为boris→conservative-source，纵向形函数均linear。两组严格账本均开启；称“能量耦合版”并不等于全lab能量与动量的物理验收已通过。

显式在node1耗7:25:18，能量版在node3耗12:41:07；全程平均步时1.113→1.903 s，增加约71%。这是各一次、不同节点、包含诊断与I/O的整场成本，不能当作纯迭代算法成本。各rank峰值RSS之和约11056→19331 MiB，不代表同一瞬间的整作业RSS。

数据与审查：COMPARISON.json、comparison.csv、VALIDATION.json、COMPLETE.json及每面小NPZ缓存。分析任务19301只读已关闭原始数据；原卡、输入、程序及运行中的高阶版本保持。
'''
    (out/'README.md').write_text(text)
    print(json.dumps(far,indent=2))


if __name__=='__main__':
    top=Path(sys.argv[2]).resolve();review=top/'energy_review_20261104'
    if sys.argv[1]=='theory':theory(top,review)
    else:review_result(top,review)
