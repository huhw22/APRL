"""Fixed-lab-plane bunching versus a separate whole-window energy ledger."""
import csv
import json
import sys
from pathlib import Path

import h5py
import numpy as np
import yaml
from scipy.constants import c,h,e

from analyze_explicit import particle_read


def read(p):
    return json.loads(Path(p).read_text())


def save(p,v):
    Path(p).write_text(json.dumps(v,indent=2,allow_nan=False,
      default=lambda x:x.item() if isinstance(x,np.generic) else x.tolist())+'\n')


def quantiles(values,weights):
    order=np.argsort(values);s=np.cumsum(weights[order])
    return np.interp([.1,.5,.9],s/s[-1],values[order])


def run(top,out):
    top,out=Path(top),Path(out);out.mkdir(exist_ok=False)
    root=top/'compact-production-r1';plan=read(root/'plan.json')
    cards={case:yaml.safe_load((root/case/'run.yaml').read_text()) for case in ['boris','energy']}
    planes=cards['boris']['detectors']['particle_planes']
    assert planes==cards['energy']['detectors']['particle_planes']
    for case in cards:assert read(root/case/'PIC_COMPLETE.json')['run']['exit_code']==0
    N=plan['input_particles'];assert N==plan['active_slices']*plan['NPART']
    entry=particle_read(root/'boris/detectors/entry-particles.h5')
    assert np.array_equal(entry['source_id'],np.arange(1,N+1))
    initial_q=abs(entry['charge_C'])
    common=np.ones(N,bool);coverage={}
    for case in cards:
        coverage[case]={}
        for plane in planes:
            with h5py.File(root/case/'detectors'/(plane['name']+'.h5'),'r') as f:
                g=f['particle_plane'];assert int(g['complete'][()])==1
                ids=np.sort(g['records'].fields('source_id')[:])
                assert len(ids)==int(g['committed_records'][()])
            assert ids[0]>=1 and ids[-1]<=N and np.all(np.diff(ids)>0)
            seen=np.zeros(N,bool);seen[ids.astype(int)-1]=True;common&=seen
            coverage[case][plane['name']]=dict(raw_count=len(ids),
                raw_charge_fraction=float(initial_q[seen].sum()/initial_q.sum()))
    ids=np.flatnonzero(common)+1;idx=ids-1
    q=initial_q[idx];Ne=q/e;Q=q.sum();norm=q/Q
    common_report=dict(particles=int(common.sum()),input_particles=N,
      charge_fraction=float(Q/initial_q.sum()),represented_electrons=float(Ne.sum()),
      global_independent_physical_phase_reference=float(1/np.sqrt(Ne.sum())),
      shot_noise_note='Fawley physical-shot-noise input: actual entry bunching is the baseline; independent macro-particle Poisson noise is NOT used')
    save(out/'COHORT.json',dict(common=common_report,coverage=coverage))
    np.savez_compressed(out/'cohort.npz',source_id=ids)
    comparison=read(top/'energy_review_20261104/radiation_comparison_r1/COMPARISON.json')['rows']
    far=next(r for r in comparison if r['plane']=='far' and r['window']=='hann')
    energies=dict(design=h*c/e/plan['physics']['wavelength_m'],
      explicit_peak=far['explicit_peak_eV'],energy_peak=far['energy_peak_eV'])
    wavelengths={name:h*c/e/value for name,value in energies.items()}
    groupings={}
    for width in [5,10,20]:
        group=(idx//(plan['NPART']*width)).astype(int)
        charge=np.bincount(group,weights=q)
        valid=charge>0
        groupings[width]=(group,charge,valid)
    rows=[];maps={case:[] for case in cards}
    for case in cards:
        with h5py.File(root/case/'ledger/energy-ledger.h5','r') as f:
            ledger=f['energy_ledger/records'][:]
            assert int(f['energy_ledger/complete'][()])==1
        with h5py.File(root/case/'detectors/exit-transport-energy.h5','r') as f:
            attrs=f['energy_plane'].attrs
            T0=float(attrs['time_lab_origin_s']);G=float(attrs['boost_gamma']);beta=float(attrs['boost_beta'])
        with np.load(top/'ledger_review_20261105/analysis'/case/'overview.npz') as f:
            assert np.array_equal(f['step'],ledger['step'])
            ledger_R=f['boost_residual_pct_initial'].copy()
        for plane in planes:
            a=particle_read(root/case/'detectors'/(plane['name']+'.h5'))
            ai=a['source_id'].astype(int)-1;keep=common[ai];a=a[keep]
            assert np.array_equal(a['source_id'],ids)
            assert np.array_equal(abs(a['charge_C']),q)
            t=a['time_s'];qt=quantiles(t,q)
            tb=G*((qt.astype(np.longdouble)-np.longdouble(T0))-
                  np.longdouble(beta)*np.longdouble(plane['z'])/np.longdouble(c))
            steps=np.interp(np.asarray(tb,float),ledger['time_box_s'],ledger['step'])
            assert np.all(tb>=ledger['time_box_s'][0]) and np.all(tb<=ledger['time_box_s'][-1])
            residual=np.interp(float(tb[1]),ledger['time_box_s'],ledger_R)
            mapped={name:[] for name in energies}
            for label,photon in energies.items():
                # Face-local common time anchor; global phase changes, |b| does not.
                dt=t.astype(np.longdouble)-np.longdouble(qt[1])
                angle=np.asarray(np.remainder(np.longdouble(2*np.pi*photon*e/h)*dt,
                                              np.longdouble(2*np.pi)),float)
                ph=np.cos(angle)+1j*np.sin(angle)
                global_complex=np.sum(norm*ph)
                result=dict(case=case,plane=plane['name'],z_m=plane['z'],wavelength_label=label,
                  wavelength_nm=wavelengths[label]*1e9,photon_energy_eV=photon,
                  common_count=len(a),common_charge_fraction=common_report['charge_fraction'],
                  raw_count=coverage[case][plane['name']]['raw_count'],
                  global_b=float(abs(global_complex)),phase_anchor_lab_s=float(qt[1]),
                  crossing_step_p10=float(steps[0]),crossing_step_p50=float(steps[1]),crossing_step_p90=float(steps[2]),
                  whole_window_ledger_R_pct_at_crossing_p50=float(residual))
                for width,(group,charge,valid) in groupings.items():
                    real=np.bincount(group,weights=q*ph.real,minlength=len(charge))
                    imag=np.bincount(group,weights=q*ph.imag,minlength=len(charge))
                    local=np.full(len(charge),np.nan+1j*np.nan,dtype=complex)
                    local[valid]=(real[valid]+1j*imag[valid])/charge[valid]
                    rms=float(np.sqrt(np.sum(charge[valid]*abs(local[valid])**2)/Q))
                    mean_abs=float(np.sum(charge[valid]*abs(local[valid]))/Q)
                    assert 0<=rms<=1+1e-12 and abs(global_complex)<=rms+1e-12
                    result['local_rms_'+str(width)+'lambda']=rms
                    result['local_mean_abs_'+str(width)+'lambda']=mean_abs
                    result['global_over_local_rms_'+str(width)+'lambda']=float(abs(global_complex)/rms) if rms else None
                    result['phase_alignment_'+str(width)+'lambda']=float(abs(global_complex)/mean_abs) if mean_abs else None
                    if width==10:
                        mapped[label]=abs(local)
                        result['local_independent_physical_phase_reference']=float(np.sqrt(valid.sum()/Ne.sum()))
                field_name=plane['name'].removesuffix('-particles')
                power=next((r['energy_coupled_W' if case=='energy' else 'explicit_W'] for r in comparison
                    if r['plane']==field_name and r['window']=='hann'),None)
                result['record_mean_Hann_band_power_W']=power
                rows.append(result)
            maps[case].append(mapped)
            print(case,plane['name'],{r['wavelength_label']:(r['global_b'],r['local_rms_10lambda']) for r in rows[-3:]},flush=True)
        np.savez_compressed(out/(case+'-local.npz'),z_m=[p['z'] for p in planes],
          **{label:np.array([m[label] for m in maps[case]]) for label in energies})
    save(out/'RESULT.json',dict(rows=rows,wavelengths_nm={k:v*1e9 for k,v in wavelengths.items()},
      common=common_report,local_slice_definition='Fixed original source-ID groups; 1024 particles per input optical slice; primary width 10 lambda0, controls 5/20',
      continuous_BF_time_history_available=False,physical_acceptance=False,
      scope='Ten fixed-laboratory planes; whole-window ledger time markers are not per-plane conservation residuals'))
    with (out/'bunching_growth.csv').open('w') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    render(out,rows,planes,energies,top)
    document(out,rows,common_report,energies)
    save(out/'COMPLETE.json',dict(complete=True,physical_acceptance=False,only_closed_A_B_data=True))


def render(out,rows,planes,energies,top):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    plt.rcParams.update({'font.size':12,'axes.spines.top':False,'axes.spines.right':False,'lines.linewidth':1.8})
    styles=[('boris','Explicit','#24648c','o-'),('energy','Energy coupled','#8e4b75','^--')]
    fig,axes=plt.subplots(2,2,figsize=(12.4,8.2),layout='constrained')
    configs=[('design','global_b','(a) Global BF: design wavelength'),
             ('design','local_rms_10lambda','(b) Local RMS BF: design wavelength'),
             ('explicit_peak','local_rms_10lambda','(c) Local RMS: fixed explicit outlet peak'),
             ('energy_peak','local_rms_10lambda','(d) Local RMS: fixed energy-coupled outlet peak')]
    for ax,(wave,key,title) in zip(axes.flat,configs):
        for case,label,color,style in styles:
            rr=[r for r in rows if r['case']==case and r['wavelength_label']==wave]
            ax.semilogy([r['z_m'] for r in rr],[max(r[key],1e-14) for r in rr],style,color=color,label=label)
        ax.set_title(title);ax.set_xlabel('Fixed laboratory plane z (m)');ax.set_ylabel('Bunching factor')
        ax.set_ylim(1e-5,1)
        ax.grid(alpha=.2);ax.legend()
    fig.suptitle('SPARC: same-particle bunching growth; fixed material slices',fontsize=15)
    fig.savefig(out/'bunching_growth.png',dpi=210,bbox_inches='tight',pad_inches=.12);plt.close(fig)
    # Keep arrival-event intervals separate from the cumulative whole-window ledger.
    fig,axes=plt.subplots(2,1,figsize=(11,8.0),layout='constrained')
    for case,label,color,style in styles:
        rr=[r for r in rows if r['case']==case and r['wavelength_label']=='design']
        dy=-.13 if case=='boris' else .13
        for j,r in enumerate(rr):
            axes[0].hlines(j+dy,r['crossing_step_p10'],r['crossing_step_p90'],color=color,lw=2)
            axes[0].plot(r['crossing_step_p50'],j+dy,'o' if case=='boris' else '^',color=color,
                         label=label if j==0 else None)
        with np.load(top/'ledger_review_20261105/analysis'/case/'overview.npz') as a:
            axes[1].plot(a['step'],a['boost_residual_pct_initial'],style[-2:] if style.endswith('--') else '-',
                         color=color,label=label)
        for r in rr:
            if r['plane'] in ['p1-particles','p3-particles','p6-particles']:
                x=r['crossing_step_p50'];y=r['whole_window_ledger_R_pct_at_crossing_p50']
                axes[1].plot(x,y,'o' if case=='boris' else '^',color=color)
                label=r['plane'].split('-')[0]+f": b={r['global_b']:.3g}"
                axes[1].annotate(label,(x,y),xytext=(5,-24 if case=='energy' else -19),
                                 textcoords='offset points',color=color,fontsize=9)
    axes[0].set_yticks(np.arange(len(planes)),[p['name'].removesuffix('-particles') for p in planes])
    axes[0].set_title('(a) Same-particle crossings: charge-weighted 10-90% interval and median')
    axes[0].set_xlabel('Boost-frame Maxwell step');axes[0].set_ylabel('Fixed laboratory plane')
    axes[1].set_title('(b) Whole-window boost residual; labels: design BF at median crossing time')
    axes[1].set_xlabel('Boost-frame Maxwell step');axes[1].set_ylabel('Residual / initial kinetic + field (%)')
    for ax in axes:ax.set_xlim(0,23984);ax.grid(alpha=.2);ax.legend()
    axes[0].legend(loc='upper left')
    fig.savefig(out/'bunching_ledger_timing.png',dpi=210,bbox_inches='tight',pad_inches=.12);plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(12.8,5.7),layout='constrained')
    images=[]
    for ax,(case,label,color,style) in zip(axes,styles):
        with np.load(out/(case+'-local.npz')) as a:v=a['design']
        im=ax.pcolormesh(np.arange(v.shape[1]+1),np.arange(v.shape[0]+1),np.ma.masked_invalid(v),
                         norm=LogNorm(vmin=1e-5,vmax=1),cmap='magma',shading='flat')
        images.append(im);ax.set_title(label)
        ax.set_yticks(np.arange(len(planes))+.5,[p['name'].removesuffix('-particles') for p in planes])
        ax.set_xlabel('Fixed initial material-slice group (10 design wavelengths)')
        ax.set_ylabel('Fixed laboratory plane')
    fig.colorbar(images[0],ax=axes,label='Local BF at the design wavelength')
    fig.savefig(out/'local_bunching_maps.png',dpi=210,bbox_inches='tight',pad_inches=.12);plt.close(fig)


def document(out,rows,common,energies):
    far_a=next(r for r in rows if r['case']=='boris' and r['plane']=='far-particles' and r['wavelength_label']=='design')
    far_b=next(r for r in rows if r['case']=='energy' and r['plane']=='far-particles' and r['wavelength_label']=='design')
    a_peak=next(r for r in rows if r['case']=='boris' and r['plane']=='far-particles' and r['wavelength_label']=='explicit_peak')
    text=f'''# 中间粒子面的BF增长与整体账本时间定位

两版局部BF在前两段相近，主要在第4–6段拉开差距；出磁后基本保持。设计波长处far局部RMS：显式版{far_a['local_rms_10lambda']:.6f}、能量版{far_b['local_rms_10lambda']:.6f}。显式版设计波长全局BF虽仅{far_a['global_b']:.6g}，在其实际红移出口峰处，全局BF仍为{a_peak['global_b']:.6f}、局部RMS为{a_peak['local_rms_10lambda']:.6f}，不能称聚束消失。
'''
    text+='''

在线能量账本没有BF字段；本轮从已有10个固定lab粒子面补算。两版、所有面使用同一ID交集和实际电荷权重。BF只用过面时刻和电荷，不用按中间动量插值重算的粒子能量。

'''
    text+=f"共同粒子数{common['particles']:,}，占输入电荷{100*common['charge_fraction']:.6f}%。主局部分组为10个原始光学切片，每切片1024宏粒子，约4.913 μm、16.39 fs；另保留5/20切片敏感性。\n\n"
    text+="三个固定分析波长为："+'、'.join(f"{name}={h*c/e/value*1e9:.6f} nm" for name,value in energies.items())+"。局部材料分组只是名义10个输入光学切片；不要求传播后仍构成严格不重叠的到达时间窗。\n\n"
    text+=r'''
$$b=\frac{|\sum_i|q_i|e^{i\omega t_i}|}{\sum_i|q_i|},
\qquad b_{\rm local,rms}=\sqrt{\frac{\sum_j Q_j|b_j|^2}{\sum_jQ_j}}.$$

固定源ID分组在所有位置保持，不随出口分布重分箱。相位抵消单列为global BF除以电荷加权局部|b|均值：

$$C_\phi=\frac{|\sum_jQ_jb_j|}{\sum_jQ_j|b_j|}.$$

它与global/local RMS不同，后者也受切片幅度不均匀影响，不能当作纯相位指标。全局BF下降不必意味着局部聚束消失。相位先减去面内共同时间锚点，锚点保存于CSV，避免大绝对时刻的载波计算损失精度。

输入为Fawley物理shot-noise加载，不能把独立宏粒子的1/sqrt(Nmacro)当成物理噪声底。实际入口BF为本次实测基线，1/sqrt(Ne)仅作为独立物理相位的参考尺度，不是单种子置信区间。

非设计波长的局部BF还包含有限材料切片窗口的边缘项。因此离设计波长时，必须相对该波长的实测入口基线判断增长，不能把初始局部BF直接称为物理噪声。

![BF增长](bunching_growth.png)

## 设计波长处的中间统计

|固定lab面 / m|显式全局BF|能量版全局BF|显式局部RMS|能量版局部RMS|
|---:|---:|---:|---:|---:|
'''
    for plane in [r['plane'] for r in rows if r['case']=='boris' and r['wavelength_label']=='design']:
        a=next(r for r in rows if r['case']=='boris' and r['wavelength_label']=='design' and r['plane']==plane)
        b=next(r for r in rows if r['case']=='energy' and r['wavelength_label']=='design' and r['plane']==plane)
        text+=f"|{a['z_m']:.3f}|{a['global_b']:.6g}|{b['global_b']:.6g}|{a['local_rms_10lambda']:.6g}|{b['local_rms_10lambda']:.6g}|\n"
    text+=r'''
## 时间与坐标的限制

固定lab面上的整束BF积累了不同时刻的过面事件。boost账本则是每10步、整个活动窗加同cohort移除事件的累计收支。两者不能当作同一瞬间的同一群局部粒子账本。

账本保持原始全部1,419,264粒子的cohort；BF为两版所有面共同可见的1,400,688个ID（99.948538%电荷），没有将原账本重算成这个子集。相近覆盖不代表两个统计量已成为严格同一控制体。

下图给每个面的电荷加权中位过面时刻及10–90%时间区间，再单列全束账本。区间重叠是相对论同时性和长束扫描的实际结果。CSV的whole_window_ledger_R_pct_at_crossing_p50只是中位时刻的整个窗口累计余差，不是“该面粒子自己的守恒误差”，不能据此单独作因果判断。

![过面时间与账本](bunching_ledger_timing.png)

![固定材料切片的局部BF](local_bunching_maps.png)

设计波长与两版已测出口峰值波长分别固定计算，在同一个面板内两版使用相同波长；不会逐位置追峰制造增长曲线。该统计可以描述增长、失谐敏感性与切片相位抵消，单种子、现有分辨率及稀疏探测面仍不能单独证明哪种推进的FEL机制已通过理论验收。

BF增强与已测辐射增强的方向一致，但不能只用某一个末端BF的平方定量解释6.47倍功率差；完整的场增长、频率与相位历史也参与辐射。

完整数据：bunching_growth.csv、RESULT.json、COHORT.json、cohort.npz及每版局部BF矩阵。只读19294_0/19294_1已关闭数据，没有新PIC传播，也没有更改运行中的高阶工况。
'''
    (out/'README.md').write_text(text)


if __name__=='__main__':
    run(*sys.argv[1:3])
