"""Overall energy scales and relative residuals of two completed SPARC runs."""
import csv
import json
import sys
from pathlib import Path

import numpy as np

from report_discrete import load, calculate
from relative_ledger import construct, write_tables

LD = np.longdouble


def save(p,v):
    Path(p).write_text(json.dumps(v,indent=2,allow_nan=False)+'\n')


def main(root,out):
    root,out=Path(root),Path(out);out.mkdir(exist_ok=False)
    reports={};curves={}
    for case in ['boris','energy']:
        assert json.loads((root/case/'PIC_COMPLETE.json').read_text())['run']['exit_code']==0
        rows,sums,ranks,steps=load(root/case/'ledger',root/case/'ledger/energy-ledger.h5')
        index,q,summary=calculate(rows,sums,steps)
        absolute,pct,meta=construct(rows,sums,index,q)
        d=out/case;write_tables(d,index,absolute,pct,meta)
        a=lambda n:np.asarray(rows[n],dtype=LD)
        cumulative=lambda n:np.r_[LD(0),np.cumsum(sums[n],dtype=LD)][index]
        residual=q['discrete_total_residual_J']
        K0=LD(rows['particle_kinetic_active_J'][0])+LD(rows['particle_kinetic_removed_J'][0])
        E0=LD(rows['particle_total_active_J'][0])+LD(rows['particle_total_removed_J'][0])
        U0=LD(sums['native_begin_J'][0])
        D=K0+U0
        Klab0=LD(rows['particle_kinetic_lab_active_J'][0])+LD(rows['particle_kinetic_lab_removed_J'][0])
        numerical=(cumulative('time_electric_J')+cumulative('time_magnetic_J')+
          cumulative('current_time_pairing_J')+cumulative('cowan_operator_J'))
        paired_grid=a('particle_grid_work_J')-a('boundary_source_work_J')-a('frozen_constraint_work_J')-a('carrier_constraint_work_J')
        G=LD(rows['boost_gamma'][0]);beta=LD(rows['boost_beta'][0])
        transform_error=a('delta_particle_kinetic_lab_J')-G*(a('delta_particle_kinetic_box_J')+beta*a('delta_particle_cPz_box_J'))
        assert np.max(abs(transform_error))<1e-11*Klab0
        peak=np.max(abs(residual));end=residual[-1]
        final={name:float(v[-1]) for name,v in q.items()}
        report=dict(case=case,steps=steps,mpi_ranks=ranks,physical_samples=len(rows),
          initial_boost_particle_kinetic_J=float(K0),initial_boost_particle_total_with_rest_J=float(E0),
          initial_native_field_J=float(U0),initial_boost_dynamic_energy_J=float(D),
          initial_lab_particle_kinetic_J=float(Klab0),
          final_boost_particle_change_J=float(a('delta_particle_kinetic_box_J')[-1]),
          final_lab_particle_change_J=float(a('delta_particle_kinetic_lab_J')[-1]),
          final_boost_cPz_change_J=float(a('delta_particle_cPz_box_J')[-1]),
          final_residual_J=float(end),maximum_abs_residual_J=float(peak),
          maximum_abs_residual_step=int(index[np.argmax(abs(residual))]),
          final_residual_pct_initial_dynamic=float(100*end/D),
          maximum_abs_residual_pct_initial_dynamic=float(100*peak/D),
          final_residual_pct_initial_including_rest=float(100*end/(E0+U0)),
          final_residual_pct_abs_final_particle_change=float(100*end/abs(a('delta_particle_kinetic_box_J')[-1])),
          final_residual_pct_abs_net_particle_grid_work=float(100*end/abs(paired_grid[-1])),
          final_residual_pct_peak_particle_change=meta['final_total_residual_pct_fixed_scale'],
          maximum_abs_residual_pct_peak_particle_change=meta['maximum_absolute_total_residual_pct_fixed_scale'],
          final_particle_change_boost_pct_initial_kinetic=float(100*a('delta_particle_kinetic_box_J')[-1]/K0),
          final_particle_change_lab_pct_initial_kinetic=float(100*a('delta_particle_kinetic_lab_J')[-1]/Klab0),
          lab_particle_energy_component_J=float(G*a('delta_particle_kinetic_box_J')[-1]),
          lab_particle_momentum_component_J=float(G*beta*a('delta_particle_cPz_box_J')[-1]),
          maximum_particle_lorentz_identity_error_J=float(np.max(abs(transform_error))),
          final_numerical_time_and_Cowan_channel_J=float(numerical[-1]),
          final_numerical_channel_pct_initial_dynamic=float(100*numerical[-1]/D),
          final_native_field_change_J=float(absolute['native_field_change_J'][-1]),
          final_net_particle_grid_work_J=float(paired_grid[-1]),
          maximum_abs_field_operator_residual_J=summary['max_abs_cumulative_field_operator_residual_J'],
          maximum_serialization_roundoff_bound_pct_initial_dynamic=float(100*LD(summary['max_serialization_roundoff_bound_J'])/D),
          strict_accounting_summary=summary,final_accounting_channels=final,
          closed_lab_total_residual_available=False,
          lab_missing='Field longitudinal momentum, Maxwell stress, and complete source/boundary four-momentum exchanges',
          physical_conservation_accepted=False)
        reports[case]=report;save(d/'SUMMARY.json',report)
        values=dict(boost_residual_pct_initial=100*residual/D,
          boost_particle_change_pct_initial=100*a('delta_particle_kinetic_box_J')/K0,
          lab_particle_change_pct_initial=100*a('delta_particle_kinetic_lab_J')/Klab0,
          numerical_channel_pct_initial=100*numerical/D,
          boost_residual_pct_peak_change=100*residual/LD(meta['fixed_denominator_decimal_J']))
        with (d/'overview_relative.csv').open('w') as f:
            writer=csv.writer(f);writer.writerow(['step']+list(values))
            writer.writerows([[int(step)]+[str(v[j]) for v in values.values()] for j,step in enumerate(index)])
        np.savez_compressed(d/'overview.npz',step=index,**{k:np.asarray(v,dtype=float) for k,v in values.items()})
        curves[case]=(index,values)
        print(case,json.dumps(report,indent=2),flush=True)
    save(out/'RESULT.json',dict(cases=reports,physical_conservation_accepted=False,
      primary_denominator='Initial total particle kinetic energy plus native field energy; rest mass excluded',
      lab_total_residual_available=False))
    render(out,reports,curves)
    document(out,reports)
    save(out/'COMPLETE.json',dict(complete=True,lab_total_residual_available=False,
      physical_conservation_accepted=False))


def render(out,reports,curves):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size':12,'axes.spines.top':False,'axes.spines.right':False,'lines.linewidth':1.8})
    fig,axes=plt.subplots(2,1,figsize=(10.5,7.7),sharex=True,layout='constrained')
    for case,label,color,style in [('boris','Explicit','#24648c','-'),('energy','Energy coupled','#8e4b75','--')]:
        index,v=curves[case]
        axes[0].plot(index,v['boost_residual_pct_initial'],style,color=color,label=label)
        axes[1].plot(index,v['lab_particle_change_pct_initial'],style,color=color,label=label)
    axes[0].set_title('(a) Boost-frame strict accounting residual')
    axes[0].set_ylabel('Residual / initial (kinetic + field) (%)')
    axes[1].set_title('(b) Lab particle kinetic-energy change; NOT a conservation residual')
    axes[1].set_ylabel('Particle change / initial kinetic (%)')
    for ax in axes:
        ax.axhline(0,color='#888888',lw=.7);ax.grid(alpha=.2);ax.legend()
    axes[-1].set_xlabel('Maxwell step');axes[-1].set_xlim(0,23984)
    fig.savefig(out/'overall_relative.png',dpi=210,bbox_inches='tight',pad_inches=.12);plt.close(fig)


def document(out,reports):
    a,b=reports['boris'],reports['energy']
    text=f'''# 两版整体账本与相对差额

主分母为同一初始总动能加场能，约{a['initial_boost_dynamic_energy_J']*1e6:.6f} μJ。最终严格离散账本余差：显式版{a['final_residual_J']*1e6:.6f} μJ，占{a['final_residual_pct_initial_dynamic']:.4f}%；能量耦合版账面值约{b['final_residual_J']:.3g} J，占{b['final_residual_pct_initial_dynamic']:.3g}%，已经处于浮点记账尺度。这是计入时间/Cowan数值通道后的未配平差额，不是连续物理能量精度证书。lab完整守恒余差因输出缺项，不能可靠给百分数。
'''
    text+='''

本次优先回答“差额占总体多少”。boost主分母采用初始全部粒子的总动能加初始原生场能，不含静质量能；两版分母相同。另保留原相对账本的D=max|同cohort粒子动能变化|和原始J值，避免只因总量大就忽略实际交换中的相对差。

## boost系：可计算的整体离散账本

采用与当前推进、边界和供场配套的严格离散账本。它包含已量测的时间/Cowan数值通道，不能据其余差很小就认定连续物理能量已收敛。

|量|显式版|能量耦合版|
|---|---:|---:|
'''
    entries=[
      ('初始粒子总动能（μJ）','initial_boost_particle_kinetic_J',1e6),
      ('初始场能（J）','initial_native_field_J',1),
      ('主分母：动能+场能（μJ）','initial_boost_dynamic_energy_J',1e6),
      ('最终粒子动能变化（μJ）','final_boost_particle_change_J',1e6),
      ('最终严格账本余差（J）','final_residual_J',1),
      ('最终余差 / 初始动能+场能（%）','final_residual_pct_initial_dynamic',1),
      ('全程最大绝对余差 / 同分母（%）','maximum_abs_residual_pct_initial_dynamic',1),
      ('最终余差 / 实际净粒子动能变化绝对值（%）','final_residual_pct_abs_final_particle_change',1),
      ('最终余差 / 净粒子—网格功绝对值（%）','final_residual_pct_abs_net_particle_grid_work',1),
      ('时间与Cowan数值通道净值（J）','final_numerical_time_and_Cowan_channel_J',1),
      ('上述数值通道 / 初始动能+场能（%）','final_numerical_channel_pct_initial_dynamic',1)]
    for label,key,scale in entries:text+=f"|{label}|{a[key]*scale:.9g}|{b[key]*scale:.9g}|\n"
    text+='''
若把静质量能也计入“总体”，初始boost粒子总能约48.1203 μJ，得到的相对差会缩小约7倍。因此两种分母明确分开，主表不靠静质量能放小差额。完整含静能比例见各版SUMMARY.json。

时间/Cowan通道是对所选原生能量定义的离散修正，已计入本次账本；它们不是新增真实外功，也不能直接当作连续物理误差估计。其净值在两版仍不小，因此“能量耦合版未配平差额约零”只说明配套离散收支闭合，不能推出全部物理能量与lab四动量已经守恒。

## lab系：已有粒子整体量，缺少闭合场账本

当前生产输出不能重建完整lab总守恒余差：缺少场的纵向动量、Maxwell应力及完整边界/源四动量交换。严格boost能量残差不能简单乘boost gamma当作lab残差；固定lab入口出口两面也未构成完整闭合世界管。

以下是可独立核对的同cohort粒子变化，不是守恒残差：

|量|显式版|能量耦合版|
|---|---:|---:|
'''
    for label,key,scale in [
      ('初始lab粒子总动能（mJ）','initial_lab_particle_kinetic_J',1e3),
      ('最终lab粒子动能变化（μJ）','final_lab_particle_change_J',1e6),
      ('该变化 / 初始lab粒子总动能（%）','final_particle_change_lab_pct_initial_kinetic',1),
      ('lab完整守恒差额 / 总体','closed_lab_total_residual_available',None)]:
        if scale is None:text+=f"|{label}|现有数据不足|现有数据不足|\n"
        else:text+=f"|{label}|{a[key]*scale:.9g}|{b[key]*scale:.9g}|\n"
    text+=r'''
粒子变换关系为

$$\Delta K_{\rm lab}=\gamma_{\rm boost}\left(\Delta K'+\beta_{\rm boost}\Delta cP_z'\right).$$

已有粒子洛伦兹恒等式残差只是变换检查，不是lab能量守恒误差。所有粒子量包含相同初始cohort的保留成员和累计移除事件，属于boost时间切片加移除事件的统计，不冒充同一lab时刻体积分。

![整体相对量](overall_relative.png)

本轮只读已结束的19294_0/19294_1账本；没有修改推进、供场、卡或程序，也不读取运行中的高阶版。每版2400个H5记录、23984步逐rank离散记录均独立检查完整性和算子恒等式。原生J和两类相对量CSV、SUMMARY、RESULT与归一化元数据一起保留。
'''
    (out/'README.md').write_text(text)


if __name__=='__main__':
    main(*sys.argv[1:3])
