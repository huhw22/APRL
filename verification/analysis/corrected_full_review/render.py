#!/usr/bin/env python3
"""Compose the closed-run review from small saved statistics and FFT arrays."""
import argparse,csv,hashlib,json,re
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

H=4.135667696e-15
COLORS=['#657b9e','#cc6727','#327e74']
LABELS=['18938 | fixed window','18993 | previous moving window','19098 | corrected moving window']

def save(p,d):p.write_text(json.dumps(d,ensure_ascii=False,indent=2,allow_nan=False)+'\n')

def metrics(directory):
    s=json.loads((directory/'summary.json').read_text());a=np.load(directory/'spectra.npz')
    energy=a['photon_energy_eV'];de=H/(s['samples']*s['sample_interval_s']);duration=s['samples']*s['sample_interval_s'];out={}
    for window in ['rect','hann']:
        flux=np.asarray(a[window+'_native_flux_J_per_eV'],dtype=np.longdouble);bands={}
        for tag,lo,hi in [('all',0,np.inf),('below_1_eV',0,1),('90_to_93_eV',90,93),('91p3_to_91p7_eV',91.3,91.7)]:
            mask=(energy>=lo)&(energy<hi);value=float(np.sum(flux[mask],dtype=np.longdouble)*de)
            bands[tag]={'energy_J':value,'record_mean_power_W':value/duration,'bins':int(mask.sum())}
        out[window]={'bands':bands,'negative_bins':int(np.count_nonzero(flux<0)),
          'below_1_eV_fraction':bands['below_1_eV']['energy_J']/bands['all']['energy_J']}
        # A local continuum estimate, not a near-field subtraction or a line fit.
        sides=((energy>=90.5)&(energy<91.3))|((energy>=91.7)&(energy<92.5))
        target=(energy>=91.3)&(energy<91.7)
        design=np.c_[np.ones(int(sides.sum())),energy[sides]-91.5]
        coeff=np.linalg.lstsq(design,np.asarray(flux[sides],dtype=float),rcond=None)[0]
        baseline=coeff[0]+coeff[1]*(energy[target]-91.5)
        base_energy=float(np.sum(baseline,dtype=np.longdouble)*de)
        out[window]['local_continuum']={'sidebands_eV':[[90.5,91.3],[91.7,92.5]],
          'target_band_eV':[91.3,91.7],'baseline_energy_J':base_energy,
          'signed_excess_J':bands['91p3_to_91p7_eV']['energy_J']-base_energy,
          'scope':'Linear fit of neighbouring spectral density only. No source arrays changed; not a physical near-field/radiation separation.'}
    return {'summary':s,'arrays':{k:a[k] for k in a.files},'metrics':out}

def plot_spectra(cases,figs):
    fig,axes=plt.subplots(1,3,figsize=(15,4.8),layout='constrained')
    for ax,(lo,hi,title) in zip(axes,[(.03,1000.,'Total-field spectrum'),(50.,115.,'50–115 eV'),(88.,95.,'91.5 eV neighbourhood')]):
        mins=[];maxs=[]
        for case,color,label in zip(cases,COLORS,LABELS):
            energy=case['arrays']['photon_energy_eV'];mask=(energy>=lo)&(energy<=hi)
            for window,style,lw in [('rect','-',1.1),('hann','--',.8)]:
                y=case['arrays'][window+'_native_flux_J_per_eV'][mask]*1e9
                mins.append(np.min(y));maxs.append(np.max(y))
                ax.plot(energy[mask],y,style,color=color,lw=lw,label=label+(' | rect' if window=='rect' else ' | Hann'))
        scale=max(abs(min(mins)),abs(max(maxs)))
        ax.set_yscale('log' if min(mins)>0 else 'symlog',**({} if min(mins)>0 else {'linthresh':max(scale*1e-7,1e-20)}))
        if lo<1:ax.set_xscale('log')
        ax.set(xlim=(lo,hi),xlabel='Photon-energy coordinate [eV]',ylabel='Native signed flux [nJ/eV]',title=title)
        ax.axvline(91.5,color='#777777',lw=.6,ls=':');ax.grid(alpha=.2)
    handles,labels=axes[0].get_legend_handles_labels()
    fig.legend(handles,labels,loc='outside lower center',ncol=3,fontsize=8)
    fig.suptitle('Free drift at 0.5 m | full-plane raw lab E/B | identical physical aperture\nSolid: rectangular; dashed: energy-normalized Hann; DC included in totals',fontsize=11)
    fig.savefig(figs/'spectra_comparison.png',dpi=180);fig.savefig(figs/'spectra_comparison.pdf');plt.close(fig)

def plot_particles(root,figs,p,energy):
    first,last=p['planes'];a=np.load(root/'science/particle_distribution.npz');edges=a['edges_loss_eV']/1000
    fig,axes=plt.subplots(2,2,figsize=(12,8),layout='constrained')
    ax=axes[0,0];ax.stairs(a['count']/last['records']*100,edges,fill=True,color=COLORS[-1],alpha=.75)
    ax.axvline(last['mean_loss_eV']/1000,color='#222222',ls='--',lw=.9)
    ax.set(xlabel='Per-particle lab energy loss [keV]',ylabel='Particles per displayed bin [%]',
      title=f'Matched IDs: mean {last["mean_loss_eV"]/1000:.4f} keV; std {last["std_loss_eV"]/1000:.3f} keV')
    ax.text(.03,.95,'Central 99.8% displayed; positive = loss',transform=ax.transAxes,va='top',fontsize=8)
    ax=axes[0,1];values=[first['projected_std_keV'],last['projected_std_keV']]
    ax.bar(['Entry −0.05 m','Exit 0.5 m'],values,color=[COLORS[0],COLORS[-1]],width=.55)
    for i,v in enumerate(values):ax.annotate(f'{v:.3f}',(i,v),ha='center',xytext=(0,4),textcoords='offset points',fontsize=9)
    ax.set(ylabel='Projected energy std [keV]',title='Total beam energy spread');ax.set_ylim(0,max(values)*1.15)
    ax=axes[1,0];labels=list(first['bunching']);x=np.arange(len(labels))
    ax.bar(x-.18,[first['bunching'][k]*100 for k in labels],width=.36,color=COLORS[0],label='Entry')
    ax.bar(x+.18,[last['bunching'][k]*100 for k in labels],width=.36,color=COLORS[-1],label='Exit')
    ax.set(xticks=x,xticklabels=[k.replace('_',' ') for k in labels],ylabel='Bunching factor [%]',title='Two fixed wavelengths; arrival-time Fourier sums');ax.legend(frameon=False)
    ax=axes[1,1];factor=1e9
    vals=[energy['lab_change_from_boost_kinetic_J']*factor,energy['lab_change_from_boost_momentum_J']*factor,
      float(energy['exit_minus_entry_decimal_J']['particle_kinetic_lab_J'])*factor]
    ax.bar(['Γ ΔK_boost','Γβ Δ(cPz_boost)','Lab total'],vals,color=[COLORS[0],COLORS[1],COLORS[-1]],width=.6)
    for i,v in enumerate(vals):ax.annotate(f'{v:+.4f}',(i,v),ha='center',va='top',xytext=(0,-4),textcoords='offset points',fontsize=9)
    ax.axhline(0,color='#555555',lw=.7);ax.set(ylabel='Energy change of full bunch [nJ]',title='Lab energy change: kinetic and momentum terms');ax.set_ylim(min(vals)*1.25,0)
    for ax in axes.flat:ax.grid(axis='y',alpha=.15)
    fig.suptitle('19098 | all 41,732,775 particles | 0.55 m from entry to exit',fontsize=12)
    fig.savefig(figs/'particles_and_energy.png',dpi=180);fig.savefig(figs/'particles_and_energy.pdf');plt.close(fig)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,required=True);ap.add_argument('--workspace',type=Path,required=True)
    ap.add_argument('--reuse-figures',action='store_true');a=ap.parse_args()
    root=a.root.resolve();science=root/'science';assert json.loads((science/'COMPLETE.json').read_text())['passed']
    figs=root/'figures'
    if not a.reuse_figures:figs.mkdir(exist_ok=False)
    else:
        assert all((figs/name).is_file() for name in ['particles_and_energy.png','particles_and_energy.pdf','spectra_comparison.png','spectra_comparison.pdf'])
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False})
    bas=a.workspace/'results/moving_free_drift_spectra_20261029'
    dirs=[bas/'old_free_drift',bas/'new_free_drift',science/'spectrum'];cases=[metrics(d) for d in dirs]
    for case in cases:
        assert case['summary']['physical_shape_yx']==[140,140] and case['summary']['samples']==30019
        assert np.isclose(case['summary']['physical_area_m2'],.49e-6,rtol=1e-13)
    comparison={'cases':{str(j):case['metrics'] for j,case in zip([18938,18993,19098],cases)},
       'identical_frequency_grids':all(np.array_equal(cases[0]['arrays']['frequency_Hz'],case['arrays']['frequency_Hz']) for case in cases[1:]),
       'same_full_record_timestamps':all(all(case['summary'][k]==cases[0]['summary'][k] for k in ('time_start_s','time_end_s','samples')) for case in cases[1:]),
       'source_cache_sha256':{str(d/name):hashlib.sha256((d/name).read_bytes()).hexdigest() for d in dirs for name in ('summary.json','spectra.npz')}}
    save(root/'spectral_comparison.json',comparison)
    p=json.loads((science/'particle_statistics.json').read_text());energy=json.loads((science/'energy_planes.json').read_text());ledger=json.loads((science/'ledger_summary.json').read_text())
    if not a.reuse_figures:plot_spectra(cases,figs);plot_particles(root,figs,p,energy)
    first,last=p['planes'];normal=ledger['normalization'];old=cases[1]['metrics'];new=cases[2]['metrics'];fixed=cases[0]['metrics']
    rows=[]
    for j,case in zip([18938,18993,19098],cases):
        for w in ['rect','hann']:
            for tag,values in case['metrics'][w]['bands'].items():rows.append({'simulation_job':j,'window':w,'band':tag,**values})
    with (root/'spectral_bands.csv').open('w') as f:
        writer=csv.DictWriter(f,fieldnames=rows[0].keys());writer.writeheader();writer.writerows(rows)
    duration=cases[2]['summary']['samples']*cases[2]['summary']['sample_interval_s']
    E=lambda case,w,band:case[w]['bands'][band]['energy_J']
    power_error=None
    times=np.load(science/'energy_plane_series.npz')
    dt=cases[2]['summary']['sample_interval_s'];powers=times['exit_physical_power_W']
    from_energy_plane=float(np.sum(powers,dtype=np.longdouble)*dt)
    from_fft=E(new,'rect','all');power_error=(from_energy_plane-from_fft)/max(abs(from_fft),1e-300)
    interpretation={
      'mean_loss_eV':last['mean_loss_eV'],'std_loss_eV':last['std_loss_eV'],
      'projected_spread_relative_change':last['projected_std_keV']/first['projected_std_keV']-1,
      'full_field_energy_ratio_vs_previous_moving':from_fft/E(old,'rect','all'),
      'full_field_energy_ratio_vs_fixed':from_fft/E(fixed,'rect','all'),
      'raw_vs_scalar_loss_difference_eV':last['mean_loss_eV']-energy['scalar_endpoint_mean_loss_eV'],
      'energy_plane_rect_sum_J':from_energy_plane,'fft_rect_time_sum_J':from_fft,'relative_difference':power_error}
    assert abs(power_error)<2e-10
    save(root/'NUMERICAL_CROSSCHECKS.json',interpretation)
    W=lambda case,w,band:case[w]['bands'][band]['record_mean_power_W']
    lines=['# 修正后全量自由漂移0.5 m：本轮结果','',
      '本轮使用当前两项入窗／参考面修正、原精细网格及Cowan/Boris/沉积/横向CPML方法。入口−0.05 m，出口与场面0.5 m，实际两面距离0.55 m；全量41732775个单位权重电子，总电荷约6.6863 pC。',
      '',f'**完整粒子面的平均lab能损为{last["mean_loss_eV"]:.3f} eV/电子，逐粒子能损标准差{last["std_loss_eV"]/1000:.4f} keV；boost严格总账本差额约为净粒子变化的{abs(normal["final_total_residual_pct_fixed_scale"]):.5f}%。**原生总场谱仍需按总场解释，不能把这些粒子能损全部归为辐射。',
      '', '## 粒子：相同ID在两面的变化','',
      '| 量 | 入口−0.05 m | 出口0.5 m |','|---|---:|---:|',
      f'| 完整记录数 | {first["records"]} | {last["records"]} |',
      f'| 平均γ | {first["mean_gamma"]:.9f} | {last["mean_gamma"]:.9f} |',
      f'| 总体能散std(K)，keV | {first["projected_std_keV"]:.6f} | {last["projected_std_keV"]:.6f} |',
      f'| 相对总体能散 | {first["relative_projected_spread"]*100:.7f}% | {last["relative_projected_spread"]*100:.7f}% |',
      f'| 横向RMS x／y，µm | {first["rms_x_um"]:.4f}／{first["rms_y_um"]:.4f} | {last["rms_x_um"]:.4f}／{last["rms_y_um"]:.4f} |',
      f'| 到达时间等效RMS长度，µm | {first["rms_arrival_length_um"]:.8f} | {last["rms_arrival_length_um"]:.8f} |',
      '', '能损定义为Δγ=γ_entry−γ_exit，正值表示损失。按真实粒子ID配对，得到：','',
      f'- avg(Δγ) = **{last["mean_loss_gamma"]:.12g}**，std(Δγ) = **{last["std_loss_gamma"]:.12g}**；两者无量纲。',
      f'- std／avg = **{last["std_over_mean_loss"]:.6f}**。整束净lab能量损失为 **{last["total_loss_lab_J"]*1e9:.9f} nJ**。',
      f'- 高精度标量端点输运给出{energy["scalar_endpoint_mean_loss_eV"]:.9f} eV/电子，与逐粒子插值动量口径相差{interpretation["raw_vs_scalar_loss_difference_eV"]:.3g} eV；两套诊断相符。',
      f'- 总体能散变化约{interpretation["projected_spread_relative_change"]*100:+.6f}%；它与逐粒子能量变化的5.35 keV展宽是不同统计量。',
      '', '总体能散满足Var(K_out)=Var(K_in)+Var(loss)−2Cov(K_in,loss)。因此逐粒子有正负能量交换，并不要求原本约1.31 MeV的总体能散显著增加。横向仍在会聚，不能把这束含位置—动量相关的生产输入当成冷球膨胀。','',
      '| 固定波长 | 入口BF | 出口BF |','|---|---:|---:|']
    for k in first['bunching']:lines.append(f'| {k.replace("_"," ")} | {first["bunching"][k]:.8g} | {last["bunching"][k]:.8g} |')
    lines += ['', 'BF由固定lab面的到达时间计算；这里仅报告两个指定波长，不代表整个波段的最大微聚束。',
      '', '![全粒子能损、总体能散、BF及能量变换](figures/particles_and_energy.png)',
      '', '[逐粒子统计](science/particle_statistics.json) · [能损分布](science/particle_distribution.npz)',
      '', '## 两套能量账本分别看','',
      '### boost移动窗内部收支','',
      f'全程固定分母D=max|同cohort粒子ΔK_boost| = **{normal["fixed_denominator_J"]*1e12:.9f} pJ**，本次恰在终态达到峰值。',
      f'终态严格总核对差 **{ledger["summary"]["final"]["discrete_total_residual_J"]:.9g} J**，为D的 **{normal["final_total_residual_pct_fixed_scale"]:+.6f}%**；最大累计场算子核对差约 **{ledger["summary"]["max_abs_cumulative_field_operator_residual_J"]:.3g} J**，总余项主要在粒子—网格交换。',
      '', '![boost系相对能量收支](science/ledger/energy_budget.png)',
      '', '这说明本次全量漂移的离散场账本与粒子功能够核对到上述水平。它不是连续空间电荷解的全部误差，也不是lab辐射能损的误差上界；严格账本包含时间／Cowan数值通道。日志旧残差字段不采用这些完整离散通道，不能与这里混用。',
      '', '[账本汇总](science/ledger_summary.json) · [相对逐步表](science/ledger/energy_budget.csv) · [原J表](science/ledger/energy_budget_absolute.csv)',
      '', '### lab出入口输运','',
      '| 量 | 入口 | 出口或出口−入口 |','|---|---:|---:|',
      f'| 全束动能 | {float(energy["planes"][0]["totals_decimal_J"]["particle_kinetic_lab_J"])*1e3:.12f} mJ | {float(energy["planes"][1]["totals_decimal_J"]["particle_kinetic_lab_J"])*1e3:.12f} mJ |',
      f'| 粒子动能变化 | — | {float(energy["exit_minus_entry_decimal_J"]["particle_kinetic_lab_J"])*1e9:+.9f} nJ |',
      f'| 物理孔径有符号总场过面能量 | {energy["planes"][0]["physical_field_energy_lab_J"]*1e9:.6f} nJ | {energy["planes"][1]["physical_field_energy_lab_J"]*1e9:.6f} nJ |',
      f'| 两面的场过面能量之差 | — | {energy["physical_field_exit_minus_entry_J"]*1e9:+.6f} nJ |',
      '',f'两面粒子boost动能差为{float(energy["exit_minus_entry_decimal_J"]["particle_kinetic_box_J"])*1e12:+.9f} pJ，cPz差为{float(energy["exit_minus_entry_decimal_J"]["particle_cpz_box_J"])*1e12:+.9f} pJ；按ΔK_lab=Γ_b[ΔK_boost+β_bΔ(cPz_boost)]分别贡献{energy["lab_change_from_boost_kinetic_J"]*1e9:+.9f}、{energy["lab_change_from_boost_momentum_J"]*1e9:+.9f} nJ，合成正确的lab粒子能量变化。',
      '', '这两面的场能量积分是穿过平面的有符号总场能流，包含束缚场；不是平面上某时刻的场能存量。两面记录覆盖不同lab时段，未封闭对应世界管的侧流和体积存量，不能把上表粒子差加场差称为lab守恒残差。没有把无磁控制直接扣成纯辐射。',
      '', '[两面完整输运汇总及高低位能量和](science/energy_planes.json)',
      '', '## 场频谱：与两份旧漂移结果比较','',
      '三组均为0.5 m原生lab总E/B，140×140物理孔径、面积0.49 mm²，完整30019样本，记录约135.18 fs。保留DC、原Yee时间层，矩形与能量归一Hann定义沿用原分析；不扣底、不移除谱线、不调整幅度。',
      '', '| 工况 | 矩形窗全频总场能量，nJ | <1 eV占比 | 91.3–91.7 eV矩形，fJ | 同频带Hann，fJ |',
      '|---|---:|---:|---:|---:|']
    for j,case in zip([18938,18993,19098],cases):
        m=case['metrics'];lines.append(f'| {j} | {E(m,"rect","all")*1e9:.6f} | {m["rect"]["below_1_eV_fraction"]*100:.5f}% | {E(m,"rect","91p3_to_91p7_eV")*1e15:.6f} | {E(m,"hann","91p3_to_91p7_eV")*1e15:.6f} |')
    lines += ['',f'当前全频矩形总场能流为前一移动窗18993的 **{interpretation["full_field_energy_ratio_vs_previous_moving"]:.6f}倍**，相对旧固定窗18938为 **{interpretation["full_field_energy_ratio_vs_fixed"]:.6f}倍**。这些比较用于判断总场背景变化，不是纯辐射产额比较。',
      '',f'91.3–91.7 eV频带的记录平均能流功率：矩形 **{W(new,"rect","91p3_to_91p7_eV"):.9g} W**，Hann **{W(new,"hann","91p3_to_91p7_eV"):.9g} W**；均包含频带本底，不能仅凭频带积分称其为辐射峰功率。',
      '',f'当前该频带Hann／矩形能量比为 **{E(new,"hann","91p3_to_91p7_eV")/E(new,"rect","91p3_to_91p7_eV"):.6f}**，前一移动窗仅为 **{E(old,"hann","91p3_to_91p7_eV")/E(old,"rect","91p3_to_91p7_eV"):.6f}**，此前矩形记录边缘造成的强窗函数敏感性明显减轻。不过，当前Hann频带本底相对18993约为 **{E(new,"hann","91p3_to_91p7_eV")/E(old,"hann","91p3_to_91p7_eV"):.6f}倍**，不能说所有高频背景都降低了。',
      '',f'用邻近90.5–91.3／91.7–92.5 eV两段作局部线性连续底估计，当前矩形目标带本底为{new["rect"]["local_continuum"]["baseline_energy_J"]*1e15:.6f} fJ、带符号超额为{new["rect"]["local_continuum"]["signed_excess_J"]*1e15:+.6f} fJ；Hann超额为{new["hann"]["local_continuum"]["signed_excess_J"]*1e15:+.6f} fJ。这与图中平滑谱底一致，没有明显孤立正窄峰证据；它只是局部谱底估算，没有完成物理近场扣除。',
      '', '![三份自由漂移场谱](figures/spectra_comparison.png)',
      '', '[频带CSV](spectral_bands.csv) · [频谱数值与对比口径](spectral_comparison.json) · [新谱原始数组](science/spectrum/spectra.npz)',
      '', '## 本轮结论与边界','',
      f'1. 当前全量自由漂移给出的平均能损为 **0.5623 keV/电子**，本组没有出现数百keV的全束平均能损；逐粒子交换展宽仍有 **5.351 keV**。这是生产输入自由漂移的实测基准，不能与不同球形工况的理论增能直接相减。',
      f'2. 严格boost收支的最终余差约为净粒子变化的 **{abs(normal["final_total_residual_pct_fixed_scale"]):.4f}%**，且场算子差处于舍入量级。该结果支持本组离散记账一致，不能据此保证所有粒子动力学或微小辐射信号都同样精确。',
      f'3. 总场全频能流相对前一移动窗下降约 **{(1-interpretation["full_field_energy_ratio_vs_previous_moving"])*100:.3f}%**，回到与旧固定窗相近的量级（仍高约{(interpretation["full_field_energy_ratio_vs_fixed"]-1)*100:.2f}%）。91.5 eV邻域当前主要表现为连续谱底，窗函数敏感性减轻，但Hann高频本底没有全面下降。',
      '4. 旧18993没有保存与本次同口径的出入口粒子面，无法给逐ID能损做直接前后比较。当前场谱与旧谱的变化可以对照；有磁19099尚未闭合，不用其运行中输出给出辐射增量或修正是否全面成功的结论。',
      '', '## 运行与数据完整性说明','',
      '19098模拟步COMPLETED 0:0，模拟推进约14336秒，最终30018步；调度总耗时3:59:25。作业最终FAILED 1:0来自收尾检查器错误假设source_id必须是1..N。真实source_id沿用裁切前SDDS编号，最大104097203；独立全量核验已确认两面各41732775个particle_id无重复、无缺失，source_id逐行对应原输入，所有科学数据完整关闭且有限。原FAILED状态和原completion-check.json保留。',
      '', '日志末行的平均γ和能散只覆盖尚在活动集合中的71476个粒子；本文使用全部41732775个过面事件，二者不可直接混比。本次开启严格逐步诊断并增加粒子面，因此约0.472秒/步不是关闭审计后的纯推进速度。',
      '', '[独立完成核验](science/COMPLETE.json) · [数值交叉核对](NUMERICAL_CROSSCHECKS.json) · [FFT独立核对](science/spectrum_independent_review.json)']
    lines += ['', '原始来源：[19098仿真日志](../r2/free_drift_0p5m/aprl.log) · [冻结计算计划与输入身份](../r2/free_drift_0p5m/run-plan.json) · [原收尾失败记录](../r2/free_drift_0p5m/completion-check.json)。']
    (root/'README.md').write_text('\n'.join(lines)+'\n')
    save(root/'FINAL.json',{'status':'REVIEWED_CLOSED_FULL_DRIFT','simulation_job':19098,'simulation_complete':True,'legacy_batch_status':'FAILED 1:0 after completed solver',
      'particle_statistics':p,'energy_planes':energy,'ledger_normalization':normal,'spectra':comparison,
      'numerical_crosschecks':interpretation,'physics_accuracy_or_pure_radiation_accepted':False})
    print(json.dumps(interpretation,ensure_ascii=False,indent=2))

if __name__=='__main__':main()
