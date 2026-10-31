#!/usr/bin/env python3
"""Present 19099's closed field/ledger and explicitly conditional particle data."""
import argparse,csv,hashlib,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.constants import e
import render as shared

JOBS=[18932,18992,19099,19098]
COLORS=['#657b9e','#cc6727','#8546a0','#327e74']
LABELS=['18932 | fixed | undulator 4.5 m','18992 | previous moving | undulator 4.5 m',
        '19099 | corrected | undulator 4.5 m','19098 | corrected | drift 0.5 m']


def plots(cases,root,p,energy,ledger):
    figs=root/'figures';figs.mkdir(exist_ok=False)
    fig,axes=plt.subplots(2,2,figsize=(14,9),layout='constrained')
    for ax,(lo,hi,title) in zip(axes.flat,[(.03,1000,'Full positive-energy coordinate'),(.03,10,'Low-frequency total field'),
                                        (50,115,'50–115 eV'),(90.5,92.5,'91.5 eV neighbourhood')]):
        values=[]
        for case,color,label in zip(cases,COLORS,LABELS):
            energy_axis=case['arrays']['photon_energy_eV'];mask=(energy_axis>=lo)&(energy_axis<=hi)
            for win,style,lw in [('rect','-',1),('hann','--',.85)]:
                y=case['arrays'][win+'_native_flux_J_per_eV'][mask]*1e9;values.extend([float(y.min()),float(y.max())])
                ax.plot(energy_axis[mask],y,style,color=color,lw=lw,label=label+(' | rect' if win=='rect' else ' | Hann'))
        ax.set_yscale('log' if min(values)>0 else 'symlog',**({} if min(values)>0 else {'linthresh':max(max(abs(v) for v in values)*1e-8,1e-20)}))
        if lo<1:ax.set_xscale('log')
        ax.set(xlim=(lo,hi),xlabel='Photon-energy coordinate [eV]',ylabel='Native signed flux [nJ/eV]',title=title)
        if hi>91:ax.axvline(91.5,color='#777777',lw=.6,ls=':')
        ax.grid(alpha=.2)
    handles,labels=axes[0,0].get_legend_handles_labels()
    fig.legend(handles,labels,loc='outside lower center',ncol=2,fontsize=8)
    fig.suptitle('Full-plane raw lab E/B | 140 × 140 aperture, 0.49 mm²\nSolid: rectangular; dashed: energy-normalized Hann; native time records; DC retained in totals',fontsize=11)
    fig.savefig(figs/'spectra_comparison.png',dpi=180);fig.savefig(figs/'spectra_comparison.pdf');plt.close(fig)
    shared.plot_particles(root,figs,p,energy,title='19099 | 41,732,699 matched arriving IDs | 4.55 m from entry to exit',
                          exit_color=COLORS[2],cohort_label='matched arriving cohort')
    # The agreed peak-change denominator can be much larger than the terminal net
    # change in a magnet. Show that second scale alongside, rather than hide it.
    data=np.genfromtxt(root/'science/ledger/energy_budget_absolute.csv',delimiter=',',names=True)
    terminal=abs(ledger['normalization']['particle_change_final_J']);assert terminal>0
    fig,axes=plt.subplots(2,1,figsize=(11,7),sharex=True,layout='constrained')
    for ax,den,label in [(axes[0],ledger['normalization']['fixed_denominator_J'],'Percent of peak |ΔK_boost|'),
                         (axes[1],terminal,'Percent of terminal |ΔK_boost|')]:
        for name,title,color,style in [('total_residual_J','Total exchange/accounting gap','#c34f37','-'),
                                     ('particle_grid_exchange_difference_J','Particle-grid exchange difference','#8546a0','--'),
                                     ('field_operator_residual_J','Field operator residual','#327e74',':')]:
            ax.plot(data['step'],100*data[name]/den,label=title,color=color,ls=style,lw=1)
        ax.set(ylabel=label);ax.axhline(0,color='#777777',lw=.6);ax.grid(alpha=.2);ax.legend(fontsize=8,frameon=False)
    axes[0].set_title('19099 | boost energy diagnostic: identical residual on two fixed scales')
    axes[1].set_xlabel('Maxwell step')
    fig.savefig(figs/'ledger_exchange_scales.png',dpi=180);fig.savefig(figs/'ledger_exchange_scales.pdf');plt.close(fig)


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,required=True);ap.add_argument('--workspace',type=Path,required=True)
    ap.add_argument('--reuse-figures',action='store_true');a=ap.parse_args();root=a.root.resolve();s=root/'science'
    done=json.loads((s/'COMPLETE.json').read_text());assert done['passed'] and not done['full_input_exit_coverage_passed']
    p=json.loads((s/'particle_statistics.json').read_text());energy=json.loads((s/'energy_planes.json').read_text());ledger=json.loads((s/'ledger_summary.json').read_text())
    coverage=json.loads((s/'COVERAGE.json').read_text());first,last=p['planes'];normal=ledger['normalization']
    base=a.workspace/'results/moving_four_case_spectra_20261029'
    dirs=[base/'old_undulator',base/'new_undulator',s/'spectrum',
          a.workspace/'results/corrected_full_comparison_20261030/free_drift_review_19098_r1/science/spectrum']
    cases=[shared.metrics(d) for d in dirs]
    for item in cases:
        assert item['summary']['physical_shape_yx']==[140,140]
        assert np.isclose(item['summary']['physical_area_m2'],.49e-6,rtol=1e-13)
    assert all(c['summary']['samples']==33487 for c in cases[:3]) and cases[3]['summary']['samples']==30019
    E=lambda case,win,band:case['metrics'][win]['bands'][band]['energy_J']
    W=lambda case,win,band:case['metrics'][win]['bands'][band]['record_mean_power_W']
    current=cases[2];previous=cases[1];fixed=cases[0];drift=cases[3]
    target='91p3_to_91p7_eV';duration=current['summary']['samples']*current['summary']['sample_interval_s']
    raw_vs_scalar=last['mean_loss_eV']-energy['scalar_endpoint_mean_loss_eV'];assert abs(raw_vs_scalar)<.01
    plane_series=np.load(s/'energy_plane_series.npz');plane_rect=float(np.sum(plane_series['exit_physical_power_W'],dtype=np.longdouble)*current['summary']['sample_interval_s'])
    flux_error=(plane_rect-E(current,'rect','all'))/abs(E(current,'rect','all'));assert abs(flux_error)<2e-10
    # Use the saved plane series rectangle consistently; its separate trapezoid
    # total differs slightly at a nonzero final endpoint.
    peak={}
    for win in ['rect','hann']:
        xx=current['arrays']['photon_energy_eV'];yy=current['arrays'][win+'_native_flux_J_per_eV'];ix=np.flatnonzero((xx>=90)&(xx<93));k=ix[np.argmax(yy[ix])]
        peak[win]={'photon_energy_eV':float(xx[k]),'spectral_density_J_per_eV':float(yy[k])}
    scalar_loss=-float(energy['exit_minus_entry_decimal_J']['particle_kinetic_lab_J'])
    correction=float(p['missing_entry_totals_J']['particle_kinetic_lab_J'])
    final_gap=ledger['summary']['final']['discrete_total_residual_J'];grid=float(ledger['final_native_ledger_record']['particle_grid_work_J'])
    direct_coupling=ledger['summary']['final']['particle_grid_coupling_gap_J']
    trace=np.genfromtxt(s/'ledger/energy_budget_absolute.csv',delimiter=',',names=True)
    gap_threshold_steps={str(frac):int(trace['step'][np.flatnonzero(abs(trace['total_residual_J'])>=frac*abs(final_gap))[0]]) for frac in [.1,.5,.9]}
    checks={'simulation_job':19099,'mean_loss_eV':last['mean_loss_eV'],'std_loss_eV':last['std_loss_eV'],
      'raw_vs_scalar_loss_difference_eV':raw_vs_scalar,'rectangle_energy_plane_J':plane_rect,
      'rectangle_FFT_J':E(current,'rect','all'),'plane_fft_relative_difference':flux_error,
      'projected_spread_relative_change':last['projected_std_keV']/first['projected_std_keV']-1,
      'missing_entry_kinetic_correction_J':correction,'naive_cohort_mismatch_extra_loss_eV':correction/(last['records']*e),
      'spectral_peaks_90_to_93_eV':peak,'target_band_to_particle_loss_rect':E(current,'rect',target)/scalar_loss,
      'target_band_to_particle_loss_hann':E(current,'hann',target)/scalar_loss,
      'field_total_ratio_vs_previous_moving':E(current,'rect','all')/E(previous,'rect','all'),
      'field_total_ratio_vs_fixed':E(current,'rect','all')/E(fixed,'rect','all'),
      'strict_gap_over_abs_terminal_boost_particle_change':final_gap/abs(normal['particle_change_final_J']),
      'strict_gap_over_abs_net_grid_work':final_gap/abs(grid),
      'strict_gap_equals_particle_grid_gap_difference_J':final_gap-direct_coupling,
      'first_steps_reaching_fraction_of_final_absolute_gap':gap_threshold_steps,
      'target_band_hann_energy_per_matched_electron_eV':E(current,'hann',target)/(last['records']*e),
      'field_operator_residual_must_not_be_used_as_particle_physics_error_bound':True}
    assert abs(checks['strict_gap_equals_particle_grid_gap_difference_J'])<1e-18
    shared.save(root/'NUMERICAL_CROSSCHECKS.json',checks)
    comparison={'cases':{str(j):v['metrics'] for j,v in zip(JOBS,cases)},'case_labels':dict(zip(map(str,JOBS),LABELS)),
      'undulator_identical_frequency_grids':all(np.array_equal(cases[0]['arrays']['frequency_Hz'],v['arrays']['frequency_Hz']) for v in cases[1:3]),
      'undulator_identical_sample_metadata':all(all(v['summary'][k]==cases[0]['summary'][k] for k in ('time_start_s','time_end_s','samples')) for v in cases[1:3]),
      'free_drift_has_its_own_shorter_record_and_frequency_grid':True,
      'source_cache_sha256':{str(d/name):hashlib.sha256((d/name).read_bytes()).hexdigest() for d in dirs for name in ('summary.json','spectra.npz')}}
    shared.save(root/'spectral_comparison.json',comparison)
    rows=[]
    for job,case in zip(JOBS,cases):
        for win in ['rect','hann']:
            for band,values in case['metrics'][win]['bands'].items():rows.append({'simulation_job':job,'window':win,'band':band,**values})
    with (root/'spectral_bands.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False})
    if not a.reuse_figures:plots(cases,root,p,energy,ledger)
    drift_p=json.loads((a.workspace/'results/corrected_full_comparison_20261030/free_drift_review_19098_r1/science/particle_statistics.json').read_text())['planes'][-1]
    text=['# 19099：修正后有磁4.5 m全量结果','',
      f'**到达出口的{last["records"]:,}个粒子，平均lab能损{last["mean_loss_eV"]/1000:.6f} keV，逐粒子能损标准差{last["std_loss_eV"]/1000:.6f} keV。**有磁组的几百keV平均变化仍然存在；当前自由漂移19098的平均能损为{drift_p["mean_loss_eV"]/1000:.6f} keV。两组传播范围不同，这不是控制变量严格相减，也不能把差值直接叫辐射损失。','',
      f'**本次有磁组不能宣布能量交换验收通过。**严格boost总差额为{final_gap*1e9:+.9f} nJ，几乎全在粒子取场做功与沉积电流的场功之间；其绝对值相当于最终净boost粒子变化的{abs(checks["strict_gap_over_abs_terminal_boost_particle_change"])*100:.3f}%。场算子自身的离散恒等式仍闭合。','',
      '## 运行和完整性','',
      '19099使用node6、node8两节点，32 MPI×8、160×160×4099活动网格，当前charge-matched-inlet程序，两项reference_surface_current／inlet_charge_match修正均开启。入口−0.05 m，粒子出口和场面4.5 m，两面距离4.55 m；原物理输入、分辨率与推进方法保持。',
      '模拟步COMPLETED 0:0，33486步，推进约12小时22分；调度总耗时12:23:44，2026-10-31 06:30:24结束。纯progress平均约1.3244秒/步，包含本次启用的诊断及输出，不当作关闭严格审计后的速度。','',
      f'原作业最终FAILED 1:0：本次编号规则已通过，收尾失败来自出口比入口少{coverage["missing_exit_count"]}条记录。入口{coverage["original_input"]:,}、出口{coverage["particle_exit_records"]:,}；两个文件事务均完整关闭，场的33487帧全部完整、有限。原失败记录保留，未放宽成全输入覆盖通过。',
      f'日志另有{coverage["transverse_exit_total"]}次横向CPML退出（x−124、x+25、y−11、y+13，纵向0）。这与部分粒子过面前退出、部分过面后退出相容；日志没有逐ID边界事件，不能证明这76个ID与退出计数一一对应。','',
      f'本文粒子统计只用出入口共同ID，入口也同步去掉未到出口的76个ID。它们在入口的动能共{correction*1e9:.9f} nJ；若直接拿全部入口减幸存者出口，会额外虚增约{checks["naive_cohort_mismatch_extra_loss_eV"]:.6f} eV/电子的表观能损。未知出口状态没有补零或外推。','',
      '[缺失出口ID及入口状态](science/missing_exit_particles.csv) · [覆盖说明](science/COVERAGE.json) · [原收尾记录](../r3_2nodes/undulator_4p5m/completion-check.json)','',
      '## 共同ID粒子统计','',
      '| 量 | 入口−0.05 m，共同ID | 出口4.5 m |','|---|---:|---:|',
      f'| 记录数 | {first["records"]} | {last["records"]} |',
      f'| 平均γ | {first["mean_gamma"]:.9f} | {last["mean_gamma"]:.9f} |',
      f'| 总体能散std(K)，keV | {first["projected_std_keV"]:.6f} | {last["projected_std_keV"]:.6f} |',
      f'| 相对总体能散 | {first["relative_projected_spread"]*100:.7f}% | {last["relative_projected_spread"]*100:.7f}% |',
      f'| RMS x／y，µm | {first["rms_x_um"]:.4f}／{first["rms_y_um"]:.4f} | {last["rms_x_um"]:.4f}／{last["rms_y_um"]:.4f} |',
      f'| 到达时间等效RMS长度，µm | {first["rms_arrival_length_um"]:.7f} | {last["rms_arrival_length_um"]:.7f} |','',
      '逐粒子定义Δγ=γ_entry−γ_exit，正值表示损失，std为总体标准差ddof=0；本输入各宏粒子恰代表一个电子。','',
      f'- avg(Δγ) = **{last["mean_loss_gamma"]:.12g}**；std(Δγ) = **{last["std_loss_gamma"]:.12g}**，均无量纲。',
      f'- std／avg = **{last["std_over_mean_loss"]:.9f}**；共同粒子总lab能损 **{last["total_loss_lab_J"]*1e6:.9f} µJ**。',
      f'- 高低位标量端点输运给出{energy["scalar_endpoint_mean_loss_eV"]:.9f} eV/电子，与逐ID动量口径相差{raw_vs_scalar:.3g} eV。',
      f'- 总体能散变化{checks["projected_spread_relative_change"]*100:+.6f}%；逐粒子能损展宽与总体能散是不同统计量。','',
      '| 固定波长 | 入口BF | 出口BF |','|---|---:|---:|']
    for k in first['bunching']:text.append(f'| {k.replace("_"," ")} | {first["bunching"][k]:.9g} | {last["bunching"][k]:.9g} |')
    text+=['','BF按固定lab面的到达时间求复相位和；两个固定波长不代表切片BF或最强聚束频率。',
      '', '![共同ID的能损、能散、BF和能量变换](figures/particles_and_energy.png)','',
      '[全量共同ID统计及原全入口统计](science/particle_statistics.json) · [能损分布](science/particle_distribution.npz)','',
      '## boost账本：关键问题在粒子—网格能量交换','',
      f'固定分母仍遵守既定定义D=max|稳定同cohort粒子ΔK_boost|，本组为 **{normal["fixed_denominator_J"]*1e9:.9f} nJ**；波荡器外场使中途动能变化很大，终态净变化却只有 **{normal["particle_change_final_J"]*1e9:+.9f} nJ**。因此只说余差是峰值D的{normal["final_total_residual_pct_fixed_scale"]:+.6f}%会掩盖它相对终态小能量交换很大的事实。','',
      '| 终态量，boost系 | 数值 |','|---|---:|',
      f'| 严格总差额 | {final_gap*1e9:+.9f} nJ |',
      f'| 粒子—网格交换差 | {direct_coupling*1e9:+.9f} nJ |',
      f'| 场算子累计最大核对差 | {ledger["summary"]["max_abs_cumulative_field_operator_residual_J"]:.6g} J |',
      f'| 粒子Boris功记账差 | {float(ledger["final_native_ledger_record"]["particle_work_balance_residual_J"]):.6g} J |',
      f'| 粒子受到的网格场净功 | {grid*1e9:+.9f} nJ |',
      f'| 与最终净boost粒子变化相比 | {checks["strict_gap_over_abs_terminal_boost_particle_change"]*100:+.6f}% |',
      f'| 与网格场净功相比 | {checks["strict_gap_over_abs_net_grid_work"]*100:+.6f}% |','',
      '交换差直接由粒子Boris电踢中记录的网格功，扣除约束／供入源功，再减实际场推进采用的累计J·(E_old+E_new)/2得到。它不是把总残差重新命名；场算子、时间配对、Cowan、边界重映射、供场和CPML通道已独立量测，分解恒等式通过。',
      '本组差额主要落在粒子采样场与沉积电流的能量交换通道；场边界覆盖遗漏不再是这份严格差额的主要来源。这里定位到通道，尚未定位到具体的时间配对、空间插值、粒子子步或其他实现因素；没有修改推进来强行闭合。','',
      f'累计差额首次达到终态绝对值10%、50%、90%，分别在第{gap_threshold_steps["0.1"]}、{gap_threshold_steps["0.5"]}、{gap_threshold_steps["0.9"]}步。曲线在较长演化区间持续积累，不表现为仅一个启动尖峰。','',
      '![严格差额相对峰值和相对终态净变化](figures/ledger_exchange_scales.png)','',
      '![按统一峰值分母绘制的完整收支](science/ledger/energy_budget.png)','',
      '全局boost账本跟踪全部原始粒子，包括移除者；固定lab面统计只跟踪到达出口的共同ID，两个时空终点不同，不能要求净变化逐位相同。以上差额也不能乘一个Γ就当作lab辐射误差上界，还缺动量误差的对应分量。',
      '', '[严格汇总](science/ledger_summary.json) · [相对逐步数据](science/ledger/energy_budget.csv) · [J单位原表](science/ledger/energy_budget_absolute.csv)','',
      '## lab两面输运','',
      f'共同ID的boost动能差为{float(energy["exit_minus_entry_decimal_J"]["particle_kinetic_box_J"])*1e12:+.6f} pJ，cPz差为{float(energy["exit_minus_entry_decimal_J"]["particle_cpz_box_J"])*1e12:+.6f} pJ。按ΔK_lab=Γ_b[ΔK_boost+β_bΔ(cPz_boost)]，两项分别贡献{energy["lab_change_from_boost_kinetic_J"]*1e9:+.6f}、{energy["lab_change_from_boost_momentum_J"]*1e9:+.6f} nJ，合成{float(energy["exit_minus_entry_decimal_J"]["particle_kinetic_lab_J"])*1e9:+.6f} nJ；动量项占主导。',
      f'入口／出口物理孔径总场有符号过面能量为 **{energy["planes"][0]["physical_field_energy_lab_J"]*1e9:.6f}／{energy["planes"][1]["physical_field_energy_lab_J"]*1e9:.6f} nJ**，差{energy["physical_field_exit_minus_entry_J"]*1e9:+.6f} nJ。这是总场能流积分，含束缚场；两面时段未封闭侧流和体积存量，不能把它与粒子损失之差称为lab守恒残差。','',
      f'两面在线能流采用梯形积分，FFT全频总量采用全部帧的矩形ΣPΔt；出口两者相差{(plane_rect-energy["planes"][1]["physical_field_energy_lab_J"])*1e9:.9f} nJ，来自非零端点。改用同一矩形口径时，独立能量面的时域积分与FFT总量相对差为{flux_error:.3g}。','',
      '[端点高低位能量及独立洛伦兹复核](science/energy_planes.json)','',
      '## 总场谱：旧有磁两组、新有磁和新自由漂移','',
      '四组均取140×140物理孔径、0.49 mm²；原始lab E/B保留Yee时间层、DC和带符号能流。矩形及能量归一Hann定义与前次一致：Hann除以sqrt(mean(w²))，不逐曲线匹配幅度、不扣除自由漂移、不改变频率格。三组有磁33487点，漂移30019点，分别使用完整记录。','',
      '| 作业／工况 | 矩形全频总场能量，nJ | <1 eV占比 | 91.3–91.7 eV矩形，nJ | 同带Hann，nJ |','|---|---:|---:|---:|---:|']
    for job,case in zip(JOBS,cases):text.append(f'| {job}／'+('自由漂移0.5 m' if job==19098 else '有磁4.5 m')+f' | {E(case,"rect","all")*1e9:.6f} | {case["metrics"]["rect"]["below_1_eV_fraction"]*100:.5f}% | {E(case,"rect",target)*1e9:.9f} | {E(case,"hann",target)*1e9:.9f} |')
    text+=['',f'当前有磁90–93 eV范围最大谱峰坐标：矩形 **{peak["rect"]["photon_energy_eV"]:.9f} eV**、Hann **{peak["hann"]["photon_energy_eV"]:.9f} eV**。91.3–91.7 eV完整带积分对应记录平均功率，矩形 **{W(current,"rect",target):.6f} W**，Hann **{W(current,"hann",target):.6f} W**；记录时长约{duration*1e15:.6f} fs，不是波荡器全束流平均功率定义。',
      f'与前一移动窗18992相比，该频带矩形变化 **{(E(current,"rect",target)/E(previous,"rect",target)-1)*100:+.6f}%**，Hann变化 **{(E(current,"hann",target)/E(previous,"hann",target)-1)*100:+.6f}%**。当前Hann目标带是新漂移带积分的 **{E(current,"hann",target)/E(drift,"hann",target):.6g}倍**；比较的是保留本底的带积分。',
      f'用90.5–91.3、91.7–92.5 eV两侧线性拟合局部连续底，当前矩形／Hann目标带底估计为{current["metrics"]["rect"]["local_continuum"]["baseline_energy_J"]*1e15:.6g}／{current["metrics"]["hann"]["local_continuum"]["baseline_energy_J"]*1e15:.6g} fJ；带符号超额为{current["metrics"]["rect"]["local_continuum"]["signed_excess_J"]*1e9:+.9f}／{current["metrics"]["hann"]["local_continuum"]["signed_excess_J"]*1e9:+.9f} nJ。这里只是局部谱底估计，不构成物理近场分离。',
      f'当前全频矩形总场能流为前一移动窗的{checks["field_total_ratio_vs_previous_moving"]:.6f}倍、旧固定窗的{checks["field_total_ratio_vs_fixed"]:.6f}倍。目标Hann频带能量折合约 **{checks["target_band_hann_energy_per_matched_electron_eV"]:.6f} eV/电子**，只相当于共同粒子能损的 **{checks["target_band_to_particle_loss_hann"]*100:.6f}%**；全频总场和粒子能损均不能称为纯辐射。此对照同时含参考／入口修正及MPI分配等配置差异，不唯一归因于某一项改动。','',
      '![四组原生总场频谱](figures/spectra_comparison.png)','',
      '[频带积分CSV](spectral_bands.csv) · [频谱比较及缓存身份](spectral_comparison.json) · [当前完整谱数组](science/spectrum/spectra.npz)','',
      '## 本轮可以确定的结论','',
      '1. 模拟完整结束；场和账本已闭合文件，原始输入到粒子出口的覆盖少76个ID，单列保留，未伪装成全覆盖验收通过。',
      f'2. 到达出口的平均lab能损仍为 **{last["mean_loss_eV"]/1000:.3f} keV**，不能用自由漂移的小损失说明有磁组问题已解决。',
      f'3. 严格场算子账本闭合，但粒子—电流能量交换差达 **{abs(final_gap)*1e9:.6f} nJ（boost）**，本轮有磁能量交换尚未通过。后续应优先定位这个通道，同时保持原核方法，不能仅以峰值分母百分比小来忽略。',
      '4. 场谱显示的带积分、局部谱底、粒子能损和全频总场能量分别报告；未把其中任一量当成纯辐射总能，也未把严格余差当成谱线误差上界。','',
      '## 分析与来源','',
      '本轮仅对已关闭文件做只读后处理，没有新PIC、没有重编程序、没有修改物理卡或原始结果。共同ID正负控、全量ID唯一性／来源映射、端点标量输运、严格场逐rank恒等式、FFT Parseval及独立复核、能流时域—频域矩形积分均检查；这些通过代表分析完整性，不代表物理精度通过。PNG目视检查另记；PDF由同源输出。','',
      '[分析完整性](science/COMPLETE.json) · [共同ID算术检查](SUBSET_CHECKS.json) · [数值交叉核对](NUMERICAL_CROSSCHECKS.json) · [FFT独立复核](science/spectrum_independent_review.json)',
      '', '[原始日志](../r3_2nodes/undulator_4p5m/aprl.log) · [冻结运行计划](../r3_2nodes/undulator_4p5m/run-plan.json) · [当前漂移报告](../free_drift_review_19098_r1/README.md)']
    (root/'README.md').write_text('\n'.join(text)+'\n')
    shared.save(root/'FINAL.json',{'status':'REVIEWED_CLOSED_UNDULATOR_WITH_EXIT_SUBSET_AND_LARGE_EXCHANGE_GAP',
      'simulation_job':19099,'analysis_job':done['analysis_job'],'solver_complete':True,'original_batch_status':'FAILED 1:0',
      'full_input_exit_coverage_passed':False,'particle_statistics':p,'energy_planes':energy,
      'ledger_normalization':normal,'spectral_comparison':comparison,'numerical_crosschecks':checks,
      'physical_energy_exchange_accepted':False,'pure_radiation_separation_accepted':False,'new_pic_simulations':0})
    print(json.dumps(checks,ensure_ascii=False,indent=2))


if __name__=='__main__':main()
