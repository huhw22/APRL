"""Validate and present the completed19783 calculation and repaired analysis."""
import csv
import json
import shutil
import sys
from pathlib import Path
import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from common import read,save,sha,frozen


def main(root):
    root=Path(root);o=root/'review-20261114-r1';a=o/'analysis';r=read(a/'RESULT.json');rows=r['segments'];last=rows[-1]
    assert r['analysis_complete'] and len(rows)==6
    frozen(root/'FROZEN.json');protected=read(o/'PROTECTED_OUTPUTS.json')
    for name,stat in protected.items():
        st=Path(name).stat();assert stat==dict(size=st.st_size,mtime_ns=st.st_mtime_ns,inode=st.st_ino)
    with (a/'band_statistics.csv').open() as f:bands=list(csv.DictReader(f))
    checks=r['field_checks'];validation=dict(passed=True,physical_acceptance=False,
        production_frozen_hashes=len(read(root/'FROZEN.json')),closed_H5_inode_size_mtime_unchanged=len(protected),
        native_import=read(root/'NATIVE_VALIDATED.json'),binding=read(root/'BINDING.json'),
        original_physics_exit_code=read(root/'main/puffin.log.launch.json')['exit_code'],
        complete_particle_records=1734656,constant_physical_charge_C=rows[-1]['charge_C'],
        maximum_Parseval_relative_error=max(float(b['Parseval_relative_error']) for b in bands),
        maximum_native_float32_profile_replay_error=max(c['native_float32_profile_replay_relative_error'] for c in checks.values()),
        maximum_double_vs_native_integrated_energy_error=max(c['native_power_integration_relative_error'] for c in checks.values()),
        original_analysis_failure='A double-precision threshold was applied to the native single-precision trapezoidal accumulator. The physical run exited0.',
        repair='Keep the independent high-precision integral; separately replay every original float32 trapezoidal addition and compare to the stored native profile.',
        full_physical_energy_balance_measured=False,PDF_separately_rendered=False)
    info=read(root/'INPUT_VALIDATED.json');precheck=[]
    for path in sorted((root/'import-check').glob('main_electrons_*.h5')):
        with h5py.File(path) as f:
            b=f['electrons'][:];step=int(f['runInfo'].attrs['iCsteps'])
        mean=np.sum(b[:,5].astype(np.longdouble)*b[:,6].astype(np.longdouble))/np.sum(b[:,6],dtype=np.longdouble)
        delta=(mean-np.longdouble(info['initial_mean_gamma_scaled']))*np.longdouble(info['gamma_ref'])*np.longdouble('9.1093837139e-31')*299792458**2/np.longdouble('1.602176634e-19')
        precheck.append(dict(file=path.name,step=step,mean_energy_change_eV=float(delta)))
    assert max(abs(x['mean_energy_change_eV']) for x in precheck)<1e-5
    validation['field_feedback_off_precheck_gamma']=precheck
    validation['initial_imported_seed_field_energy_J']=r['initial_forward_field']['all_grid_retained_forward_field_energy_J']
    save(o/'VALIDATION.json',validation)
    comparison=r['comparison'];loss=-last['particle_delta_K_lab_J'];field=last['retained_forward_field_gain_all_grid_J'];gap=loss-field
    summary=dict(analysis_complete=True,physical_acceptance=False,physics_job=19783,physics_exit_code=0,
        original_batch_state='FAILED1:0 during automatic postprocessing',readonly_review_job=19789,
        final_mean_loss_eV=-last['mean_delta_K_eV_per_electron'],final_particle_loss_J=loss,
        final_all_grid_retained_field_J=field,final_aperture520_field_J=last['field_energy_520um_J'],
        particle_loss_minus_retained_field_J=gap,unrecorded_channel_fraction_of_loss=gap/loss,
        final_full_cohort_comparison=comparison,no_new_physics_run=True)
    save(o/'RESULT_REVIEW.json',summary)
    z=[x['actual_lab_z_m'] for x in rows]
    with (a/'native_energy_history.csv').open() as f:dense=list(csv.DictReader(f))
    fig,axes=plt.subplots(2,2,figsize=(13,8),layout='constrained')
    axes[0,0].plot(z,[-x['mean_delta_K_eV_per_electron']/1000 for x in rows],'o-',label='Puffin')
    colors=['#d38a35','#6a9270'];labels=['APRL: energy + three-point','APRL: energy + momentum']
    for row,color,label in zip(comparison[:2],colors,labels):
        axes[0,0].scatter([z[-1]],[-row['mean_delta_eV']/1000],marker='s',color=color,label=label,zorder=4)
        axes[0,1].scatter([z[-1]],[row['exit_aperture520_field_J']*1e6],marker='s',color=color,label=label+' (includes near field)',zorder=4)
    axes[0,0].set(xlabel='Lab z (m)',ylabel='Mean particle loss (keV/electron)',title='Full-cohort laboratory particle losses')
    axes[0,1].plot([float(d['lab_z_m']) for d in dense],[float(d['all_grid_native_field_energy_J'])*1e6 for d in dense],label='Puffin: all grid, retained modes')
    for key,label,style in [('field_energy_520um_J','Puffin: aperture ±520µm','o-'),('Rect_resonant_energy_520um_J','Puffin: Rect ±5%, aperture ±520µm','.--')]:
        axes[0,1].plot(z,[x[key]*1e6 for x in rows],style,label=label)
    axes[0,1].set(xlabel='Lab z (m)',ylabel='Exit-plane field energy (µJ)',title='Unwindowed energy and Rect spectral band')
    axes[1,0].plot(z,[-x['particle_delta_K_lab_J']*1e6 for x in rows],'o-',label='Particle energy loss')
    axes[1,0].plot(z,[x['retained_forward_field_gain_all_grid_J']*1e6 for x in rows],'s-',label='Retained forward-field gain, all grid')
    axes[1,0].set(xlabel='Lab z (m)',ylabel='Energy (µJ)',title='Puffin particle loss and currently retained field')
    axes[1,1].plot(z,[-x['particle_plus_retained_field_change_J']*1e6 for x in rows],'o-',color='#a75a57',label='Particle loss - retained field gain')
    axes[1,1].set(xlabel='Lab z (m)',ylabel='Difference (µJ)',title='Cumulative filter/absorber losses not recorded')
    axes[1,1].text(.04,.92,'This is not a complete physical\nworldtube closure residual.',va='top',transform=axes[1,1].transAxes,fontsize=9)
    for ax in axes.flat:ax.grid(alpha=.2);ax.axhline(0,color='grey',lw=.6);ax.legend(fontsize=7)
    fig.suptitle('Six-segment comparison | same historical source events and macro charges\nAPRL retains PIC near fields; Puffin is forward/paraxial with its native filter. No amplitude fitting.')
    fig.savefig(o/'comparison_review.png',dpi=175);fig.savefig(o/'comparison_review.pdf');plt.close(fig)
    lines=['# 六段Puffin对照：19783已完成物理推进，后处理已补齐','',
        '19783物理推进正常退出0，完成42720步，物理阶段约5小时12分，节点node5、128 MPI×1。整批次的FAILED1:0来自原自动后处理检查，原日志和失败标记保留；只读19789已完成全部统计，没有重跑物理、修改原程序或冻结卡。','',
        '同一份1734656条原始事件及非等权宏电荷，实际原生导入七列逐位一致，场反馈关闭首步磁轨迹通过。旧相关phase-clone结构特意保留，只作19743/19744的历史对照。','',
        '| 方法 | 总粒子能损（µJ） | 平均能损（keV/真实电子） | 出口±520µm场能（µJ） |',
        '| --- | ---: | ---: | ---: |']
    for c in comparison:
        lines.append(f"| {c['method']} | {-c['particle_delta_K_lab_J']*1e6:.5f} | {-c['mean_delta_eV']/1000:.5f} | {c['exit_aperture520_field_J']*1e6:.5f} |")
    lines+=['',
        'Puffin同孔径场能比APRL三点核低约7.28%、比动量匹配高约6.73%；粒子能损则比三点核大约81.3%、比动量匹配大约133.9%。场能数值接近不代表场模型或完整谱已经一致。APRL数值沿用原完整cohort末态账本，包含退场和仍在域内的尾粒子；其第六段覆盖/粒子时间盖缺项仍保留。Puffin为完整束团的第六段出口。','',
        '![对照图](comparison_review.png)','',
        '## Puffin分段结果','',
        '| 段数 | 平均粒子能损（keV） | 全网格保留场能（µJ） | ±520µm保留场能（µJ） | 全束BF |',
        '| --- | ---: | ---: | ---: | ---: |']
    for x in rows:
        lines.append(f"| {x['segment']} | {-x['mean_delta_K_eV_per_electron']/1000:.5f} | {x['retained_forward_field_gain_all_grid_J']*1e6:.6f} | {x['field_energy_520um_J']*1e6:.6f} | {x['global_b1']:.6f} |")
    lines+=['',
        f"Puffin末端束流损失{loss*1e6:.5f} µJ，全网格当前保留场能{field*1e6:.5f} µJ，差{gap*1e6:.5f} µJ，占粒子能损的{100*gap/loss:.2f}%。这次没有累计记录横向吸收和滤波去除的能量，不能把这个差值直接认定为完整守恒误差，也不能据此宣称Puffin账本已闭合。",'',
        '前3段粒子累计能损约1.007/2.101/3.178 µJ，而保留场只有0.0144/0.0196/0.0336 µJ；差额早于后期显著FEL增长就已存在。现有记录足以显示这一点，尚不能唯一分解过滤、吸收及粒子—场耦合等通道。关闭场反馈的原1周期预检中，平均gamma从步骤0到80保持在约5×10⁻¹¹ eV的舍入范围内，未见该预检中的外磁场单独造成平均能损；它不能替代完整六段的力分项审计。','',
        '## 后处理检查为什么停止','',
        '实际编译源码av_write.f90的m_trapz及m_trapz2D均声明为默认real函数。m_trapz在每次累加后都会舍入到float32，即使输入和H5场是双精度。旧脚本误用2e−11相对阈值比较独立高精度积分与原生功率，第一段差7.51e−11而拒绝，未损坏物理输出。','',
        '只读修复保留独立高精度场积分，同时对每个纵向点逐次复现原生单精度累加，最大功率轮廓相对差7.24e−16；六段积分能量与原生值的最大相对差2.20e−9。这个细小的输出精度差异不能解释µJ级能量差。没有修改Puffin源码或原冻结分析脚本。','',
        '实场单边Parseval最大相对差2.87e−16；28项原冻结SHA通过，570个闭合H5的inode/size/mtime保持。原始N/电荷全程保持、选各步第一个物理dump，避免末尾correctTrans后的粒子快照。几何输出相对名义面最多晚0.35mm，第六段实际14.94625m。','',
        '## 统计口径和限制','',
        '所有粒子平均能损均用物理电荷加权，采用与APRL一致的电子静质量和c，分母为原代表电子数约6.1543×10⁸；不是每宏粒子的总能变。Puffin内部打印的电子数采用其旧q_e常数，约相差1.03×10⁻⁷相对量；映射保留物理总电荷，不作跨程序幅度拟合。','',
        '能量账本和主图使用完整时域的未加窗场积分。±5%共振带同时保存Rect和原RMS-Hann统计；Hann是窗口加权量，不能作为未加窗总场能的物理子集直接相加。末端±520µm下Rect带能7.00638µJ，RMS-Hann带能15.10255µJ，后者较大来自窗口权重，不能拿来闭合粒子能损。','',
        'Puffin保留前向傍轴横向场，含原低于0.3ωr的滤波，不含PIC完整空间电荷/速度近场；原事件到磁入口的1cm在Puffin中作弹道映射。APRL场值含近场且时间支持受原终止范围限制。尚未做严格相同完整时空支持和同频带PIC谱对照，不宣布哪一个程序物理正确或误差已唯一定位。','',
        '## 文件入口','',
        '- [分段粒子能损与场能](analysis/energy_gain.csv)','- [同孔径末端对照](analysis/final_comparison.csv)',
        '- [每周期原生场能历史](analysis/native_energy_history.csv)','- [Rect/Hann共振带统计](analysis/band_statistics.csv)',
        '- [全部统计及场检查](analysis/RESULT.json)','- [复核记录](VALIDATION.json)','- [总结果](RESULT_REVIEW.json)','',
        '原analysis/、analysis.log、FAILED.json和生产脚本全部保留；有效完整后处理在review-20261114-r1/analysis/。PNG已目视检查，PDF同源导出未另渲染。本轮仅只读后处理，不追加新物理工况或主动监控其他任务。']
    (o/'README.md').write_text('\n'.join(lines)+'\n')
    save(root/'REVIEW_COMPLETE_20261114.json',dict(complete=True,physical_tracking_complete=True,
        physics_job=19783,original_batch_failed_in_analysis=True,readonly_analysis_job=19789,
        result_directory=str(o),original_failed_marker_preserved=True,no_physics_rerun=True,physical_acceptance=False))
    print(json.dumps(summary,indent=2),flush=True)


if __name__=='__main__':main(sys.argv[1])
