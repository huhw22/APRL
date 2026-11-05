"""Lab-plane particle kinetic plus signed field transport, closed cases only."""
import argparse
import csv
import json
from pathlib import Path

import h5py
import numpy as np
from scipy.constants import c, e

from entrance_ballistic import save, table, digest

LD=np.longdouble


def cumulative_at(a,t,aperture):
    """Endpoint estimate consistent with the saved linear-power trapezoid."""
    times=a['time_lab_s']
    tolerance=8*np.finfo(float).eps*max(abs(times[0]),abs(times[-1]),1e-30)
    assert times[0]-tolerance<=t<=times[-1]+tolerance
    t=float(np.clip(t,times[0],times[-1]))
    field='physical_field_energy_lab_J' if aperture=='physical' else 'field_energy_lab_J'
    power='physical_field_power_lab_W' if aperture=='physical' else 'field_power_lab_W'
    j=int(np.searchsorted(times,t,side='right')-1)
    if j==len(a)-1 or t==times[j]:return LD(a[field][j]),'native cumulative endpoint'
    h=LD(t)-LD(times[j]);dt=LD(times[j+1])-LD(times[j])
    if a['field_sample_count'][j+1]-a['field_sample_count'][j]==1:
        value=LD(a[field][j])+LD(a[power][j])*h+LD(.5)*(LD(a[power][j+1])-LD(a[power][j]))*h*h/dt
        return value,'linear-power boundary quadrature estimate'
    fraction=h/dt
    return LD(a[field][j])+fraction*(LD(a[field][j+1])-LD(a[field][j])),'linear cumulative boundary estimate over unsaved sample'


def run(top,out):
    top,out=Path(top).resolve(),Path(out).resolve();out.mkdir(exist_ok=False)
    root=top/'compact-production-r1'
    raw_cache={
      'energy':top/'energy_startup_origin_20261105/r2/analysis/detector_energy_comparison.json',
      'momentum':top/'momentum_startup_assessment_20261105/deferred/analysis/momentum/detector_energy_comparison.json'}
    paths=[]
    for case in ['energy','momentum']:
        assert json.loads((root/case/'PIC_COMPLETE.json').read_text())['run']['exit_code']==0
        paths.append(raw_cache[case])
        for plane in ['entry','und-entry']:
            paths += [root/case/'detectors'/(plane+'-transport-energy.h5'),
                      root/case/'detectors'/(plane+'-particles.h5')]
    before=[digest(p) for p in paths]
    arrays={};metadata={};time_rows=[];planes=[];sums=[];traces={}
    for case in ['energy','momentum']:
        cached={r['plane']:r for r in json.loads(raw_cache[case].read_text())['rows']}
        for plane in ['entry','und-entry']:
            directory=root/case/'detectors'
            with h5py.File(directory/(plane+'-transport-energy.h5'),'r') as f:
                g=f['energy_plane'];assert int(g['complete'][()])==1
                a=g['records'][:];assert len(a)==int(g['committed_records'][()])
                attrs={k:(v.decode() if isinstance(v,bytes) else v.item() if isinstance(v,np.generic) else v)
                       for k,v in g.attrs.items()}
            assert np.all(np.diff(a['time_lab_s'])>0)
            assert np.all(a['field_available']==1)
            assert a['field_gap_count'][-1]==a['outside_field_window_steps'][-1]==0
            with h5py.File(directory/(plane+'-particles.h5'),'r') as f:
                g=f['particle_plane'];assert int(g['complete'][()])==1
                p=g['records'].fields(['time_s','weight','charge_C','crossing_direction'])[:]
                assert len(p)==int(g['committed_records'][()])
            assert np.all(p['crossing_direction']==1)
            assert len(p)==a['particle_downstream_count'][-1]==a['physical_particle_downstream_count'][-1]
            assert a['particle_upstream_count'][-1]==a['physical_particle_upstream_count'][-1]==0
            Ne=np.sum(p['weight'].astype(LD)); q=abs(LD(a['particle_charge_C'][-1]))
            assert abs(q/LD(e)/Ne-1)<1e-12
            tmin,tmax=float(p['time_s'].min()),float(p['time_s'].max())
            assert a['time_lab_s'][0]<=tmin<tmax<=a['time_lab_s'][-1]
            z=float(attrs['z_lab_m']);arrays[case,plane]=a
            metadata[case,plane]=dict(attrs=attrs,Ne=Ne,q=q,tmin=tmin,tmax=tmax,z=z,
                raw_mean=cached[plane]['raw_mean_kinetic_eV'])
            traces[case+'_'+plane+'_tau_ps']=(a['time_lab_s']-z/c)*1e12
            traces[case+'_'+plane+'_field_J']=a['physical_field_energy_lab_J']
            traces[case+'_'+plane+'_power_W']=a['physical_field_power_lab_W']
    for case in ['energy','momentum']:
        # A common retarded-time interval still covers the complete particle cohort.
        bounds=[(float(a['time_lab_s'][0])-metadata[case,plane]['z']/c,
                 float(a['time_lab_s'][-1])-metadata[case,plane]['z']/c)
                for (cc,plane),a in arrays.items() if cc==case]
        tau_lo=max(v[0] for v in bounds);tau_hi=min(v[1] for v in bounds)
        for plane in ['entry','und-entry']:
            a=arrays[case,plane];m=metadata[case,plane];last=a[-1];Ne=m['Ne']
            assert tau_lo+m['z']/c<=m['tmin'] and tau_hi+m['z']/c>=m['tmax']
            for aperture in ['physical','full']:
                prefix='physical_' if aperture=='physical' else ''
                k=LD(last[prefix+'particle_kinetic_lab_J'])+LD(last[prefix+'particle_kinetic_lab_J_low'])
                raw=LD(m['raw_mean'])*Ne*LD(e)
                field=LD(last[prefix+'field_energy_lab_J'])
                planes.append(dict(case=case,plane=plane,z_m=m['z'],aperture=aperture,
                    input_count=int(last[prefix+'particle_downstream_count']),
                    represented_electrons=float(Ne),scalar_particle_K_J=float(k),raw_particle_K_J=float(raw),
                    field_transport_J=float(field),particle_plus_field_scalar_J=float(k+field),
                    particle_plus_field_raw_J=float(raw+field),
                    scalar_particle_eV_per_e=float(k/Ne/LD(e)),field_eV_per_e=float(field/Ne/LD(e)),
                    sum_scalar_eV_per_e=float((k+field)/Ne/LD(e)),
                    raw_minus_scalar_with_field_eV=float(((raw+field)-(k+field))/Ne/LD(e)),
                    raw_minus_scalar_particles_only_eV=float((raw-k)/Ne/LD(e))))
                for label,ta,tb in [('full_record',float(a['time_lab_s'][0]),float(a['time_lab_s'][-1])),
                    ('particle_arrival_support',m['tmin'],m['tmax']),
                    ('common_retarded_window',tau_lo+m['z']/c,tau_hi+m['z']/c)]:
                    fa,qa=cumulative_at(a,ta,aperture);fb,qb=cumulative_at(a,tb,aperture)
                    F=fb-fa
                    sums.append(dict(case=case,plane=plane,aperture=aperture,window=label,
                        t_begin_s=ta,t_end_s=tb,tau_begin_s=ta-m['z']/c,tau_end_s=tb-m['z']/c,
                        scalar_particle_K_J=float(k),raw_particle_K_J=float(raw),field_transport_J=float(F),
                        field_eV_per_e=float(F/Ne/LD(e)),sum_scalar_J=float(k+F),sum_raw_J=float(raw+F),
                        represented_electrons=float(Ne),boundary_method=qa+'; '+qb))
                fp,_=cumulative_at(a,m['tmin'],aperture);fq,_=cumulative_at(a,m['tmax'],aperture)
                time_rows.append(dict(case=case,plane=plane,aperture=aperture,
                    first_field_s=float(last['first_field_time_lab_s']),last_field_s=float(last['last_field_time_lab_s']),
                    particle_first_s=m['tmin'],particle_last_s=m['tmax'],
                    field_before_particles_J=float(fp),field_during_particles_J=float(fq-fp),
                    field_after_particles_J=float(field-fq),
                    after_particle_fraction=float((field-fq)/field),
                    field_samples=int(last['field_sample_count']),saved_records=len(a),
                    field_gaps=int(last['field_gap_count']),outside_steps=int(last['outside_field_window_steps']),
                    peak_positive_power_W=float(np.max(a[prefix+'field_power_lab_W'])),
                    final_power_W=float(a[prefix+'field_power_lab_W'][-1])))
    differences=[]
    for case in ['energy','momentum']:
        for aperture in ['physical','full']:
            for window in ['full_record','particle_arrival_support','common_retarded_window']:
                rr={r['plane']:r for r in sums if r['case']==case and r['aperture']==aperture and r['window']==window}
                x,y=rr['entry'],rr['und-entry'];norm=LD(x['represented_electrons'])*LD(e)
                # Difference of the high+low scalar energies before adding the field.
                k0=metadata[case,'entry'];k1=metadata[case,'und-entry']
                a0=arrays[case,'entry'][-1];a1=arrays[case,'und-entry'][-1]
                prefix='physical_' if aperture=='physical' else ''
                dk=(LD(a1[prefix+'particle_kinetic_lab_J'])-LD(a0[prefix+'particle_kinetic_lab_J'])
                    +LD(a1[prefix+'particle_kinetic_lab_J_low'])-LD(a0[prefix+'particle_kinetic_lab_J_low']))
                df=LD(y['field_transport_J'])-LD(x['field_transport_J'])
                dr=LD(y['raw_particle_K_J'])-LD(x['raw_particle_K_J'])
                differences.append(dict(case=case,aperture=aperture,window=window,
                    delta_scalar_particle_J=float(dk),delta_field_transport_J=float(df),
                    delta_scalar_sum_J=float(dk+df),delta_raw_sum_J=float(dr+df),
                    delta_scalar_particle_eV_per_e=float(dk/norm),
                    delta_field_eV_per_e=float(df/norm),delta_scalar_sum_eV_per_e=float((dk+df)/norm),
                    delta_raw_sum_eV_per_e=float((dr+df)/norm),
                    interpretation='Difference of two end-plane transport sums; not a closed-control-volume residual.'))
    table(out/'plane_sums.csv',planes);table(out/'time_support.csv',time_rows)
    table(out/'window_sums.csv',sums);table(out/'two_plane_differences.csv',differences)
    np.savez_compressed(out/'time_traces.npz',**traces)
    result=dict(cases=['energy','momentum'],new_simulation=False,
        primary_aperture='physical inner rectangle 2.4 mm by 2.4 mm, all input particles covered',
        field_quantity='Signed cumulative lab Poynting transport, total self-consistent E/B; prescribed fields excluded; no near-field subtraction.',
        sum_quantity='Particle kinetic transport plus field energy transport, not instantaneous volume energy.',
        windows='Native full record; exact all-particle min/max support with boundary quadrature estimates; common retarded-time interval covering all particles.',
        time_stagger='Native E/B Yee staggering retained.',
        conservation_limit='Two end planes omit lateral flux, initial/final volume storage and source work; no complete lab conservation claim.',
        rows=planes,differences=differences,time_support=time_rows)
    save(out/'RESULT.json',result)
    after=[digest(p) for p in paths];assert before==after
    for r in planes:
        assert abs(r['raw_minus_scalar_with_field_eV']-r['raw_minus_scalar_particles_only_eV'])<1e-6
    save(out/'VALIDATION.json',dict(complete_sources=True,sources_unchanged=True,
        source_manifest=before,same_field_cancels_from_raw_scalar_difference=True,
        all_particles_in_both_apertures=True,no_field_gaps=True))
    render(out,result,traces);document(out,result)
    save(out/'COMPLETE.json',dict(complete=True,new_simulation=False))
    print(json.dumps([r for r in differences if r['aperture']=='physical'],indent=2),flush=True)


def render(out,result,traces):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(2,2,figsize=(12.6,8.0),layout='constrained')
    for j,(case,title,scale,unit) in enumerate([('energy','Energy coupled',1e6,'microjoule'),
                                              ('momentum','Momentum matched',1e3,'millijoule')]):
        for plane,label,color,line in [('entry','Entrance z=0','#24648c','-'),
                                       ('und-entry','Magnet entrance z=0.01 m','#b18428','--')]:
            key=case+'_'+plane
            axes[0,j].plot(traces[key+'_tau_ps'],traces[key+'_field_J']*scale,line,color=color,label=label)
            support=next(r for r in result['time_support'] if r['case']==case and r['plane']==plane and r['aperture']=='physical')
            last_tau=(support['particle_last_s']-(0 if plane=='entry' else .01)/c)*1e12
            axes[0,j].axvline(last_tau,color=color,ls=':',lw=1)
        axes[0,j].set_title(title+': cumulative signed field transport')
        axes[0,j].set_xlabel('Retarded lab time t - z/c (ps)')
        axes[0,j].set_ylabel('Field energy ('+unit+')');axes[0,j].legend(fontsize=9)
        rr=next(r for r in result['differences'] if r['case']==case and r['aperture']=='physical' and r['window']=='full_record')
        norm=1 if case=='energy' else 1e-6
        values=np.array([rr['delta_scalar_particle_eV_per_e'],rr['delta_field_eV_per_e'],rr['delta_scalar_sum_eV_per_e']])*norm
        axes[1,j].bar(range(3),values,color=['#24648c','#b18428','#444444'],width=.55)
        axes[1,j].axhline(0,color='#777777',lw=.8)
        axes[1,j].set_xticks(range(3),['Particle kinetic','Field transport','Sum'])
        axes[1,j].set_ylabel('0.01 m minus 0 m ('+('eV' if case=='energy' else 'MeV')+' / electron)')
        axes[1,j].set_title(title+': end-plane differences (not closure)')
        if case=='momentum':
            axes[1,j].set_ylim(float(min(values))*1.18,abs(float(min(values)))*.12)
        else:
            axes[1,j].set_ylim(0,float(max(values))*1.18)
        for x,y in enumerate(values):
            axes[1,j].annotate(f'{y:.4g}',(x,y),xytext=(0,5 if y>=0 else -15),
                textcoords='offset points',ha='center',fontsize=10)
    for ax in axes.flat:ax.grid(alpha=.18,axis='y')
    fig.suptitle('Lab particle + field transport: physical aperture, full stored records\n'
                 'Dotted lines mark last particle crossing; upper panels use different energy units',fontsize=12)
    fig.savefig(out/'particle_field_transport.png',dpi=220,bbox_inches='tight');plt.close(fig)


def document(out,result):
    text='# 同面加入场项：粒子动能＋场过面能量\n\n'
    text+='使用lab固定入口0 m和磁入口0.01 m，能量版与动量版均已闭合。场项为累计有符号Poynting能量输运，不是窗口瞬时场存量；J/m截面能量密度积分未当作J相加。主口径为2.4×2.4 mm²物理孔径，全部1,419,264个输入宏粒子均通过该孔径。\n\n'
    text+='![粒子与场输运](particle_field_transport.png)\n\n'
    text+='## 完整记录的两平面输运\n\n'
    text+='|方法|入口场能/J|磁入口场能/J|粒子动能变化/eV每电子|场输运变化/eV每电子|两者之和/eV每电子|\n|---|---:|---:|---:|---:|---:|\n'
    for case in ['energy','momentum']:
        p={r['plane']:r for r in result['rows'] if r['case']==case and r['aperture']=='physical'}
        d=next(r for r in result['differences'] if r['case']==case and r['aperture']=='physical' and r['window']=='full_record')
        text+=f"|{case}|{p['entry']['field_transport_J']:.9g}|{p['und-entry']['field_transport_J']:.9g}|{d['delta_scalar_particle_eV_per_e']:.9g}|{d['delta_field_eV_per_e']:.9g}|{d['delta_scalar_sum_eV_per_e']:.9g}|\n"
    text+='\n动量版入口累计场输运达到约12.36 mJ，已与全束粒子动能约12.25 mJ同量级；能量版入口只有约0.436 μJ。该差异远大于原645 eV的粒子重建差，不能用两版BF相近替代总场核验。场记录排除解析波荡器外场，因此并非把外磁铁场能加进来。\n\n'
    text+='## 粒子通过期间与之后\n\n|方法/面|完整场能/J|粒子首末过面期间/J|最后粒子之后/J|之后占比|\n|---|---:|---:|---:|---:|\n'
    for r in result['time_support']:
        if r['aperture']=='physical':
            total=r['field_before_particles_J']+r['field_during_particles_J']+r['field_after_particles_J']
            text+=f"|{r['case']}/{r['plane']}|{total:.9g}|{r['field_during_particles_J']:.9g}|{r['field_after_particles_J']:.9g}|{100*r['after_particle_fraction']:.4f}%|\n"
    text+='\n所有场面均无field gaps、无窗口缺失；field_sample_count为23985，入口保存23984条记录是诊断输出节拍筛选，不代表积分漏一步。完整记录累计值优先使用在线字段；粒子窗口边缘采用相邻功率的局部积分估计，方法写入window_sums.csv。\n\n'
    text+='动量版巨场并非仅来自粒子离开后的尾段，粒子通过期间仍为mJ量级。数据另给共同retarded-time窗口，确保两面覆盖同一完整粒子cohort，并检查轻微起止时间错位的影响。\n\n'
    text+='## 为什么645 eV差不会被场抵掉\n\n'
    text+=r'在同一面加入相同场项$F$，有$(K_{raw}+F)-(K_{scalar}+F)=K_{raw}-K_{scalar}$。所以原约−645.4 eV/电子差代数上不变；场项改变的是跨面输运比较，不能修复粒子重建口径。'+'\n\n'
    text+='## 守恒范围\n\n'
    text+='两端面输运和之差不等于完整lab守恒残差：还需侧向能流、区域初末场存量以及供场/边界约束源做功。原生Yee E/B时间层保留，未作近场扣除，不能把累计通量都叫辐射。当前结果说明动量版存在非常大的原始场输运，需要进一步定位其场演化、供场和诊断来源，尚未证实是何种物理或数值机制。\n\n'
    text+='数据：plane_sums.csv、two_plane_differences.csv、window_sums.csv、time_support.csv、time_traces.npz、RESULT.json、VALIDATION.json。完整孔径数据作为副口径保留。仅后处理，无新仿真或推进修改。\n'
    (out/'README.md').write_text(text)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('top');p.add_argument('out')
    a=p.parse_args();run(a.top,a.out)
