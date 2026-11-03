"""Closed-evidence report: radiation, feedback, phase, cost, and limitations."""
import argparse,csv,json,re,subprocess,sys
from pathlib import Path
import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from prepare import save,require

COLORS={'boris':'.45','energy':'#1764a5','momentum':'#c76b1d','prescribed':'#36804b'}
LABELS={'boris':'A: explicit','energy':'B: energy iteration','momentum':'C: longitudinal pairing',
        'prescribed':'No field feedback'}

def load_json(path):return json.loads(Path(path).read_text())

def main():
    ap=argparse.ArgumentParser();ap.add_argument('root',type=Path);a=ap.parse_args()
    root=a.root.resolve();base=root.parent
    plan=load_json(root/'campaign.json')
    figures=base/'figures';figures.mkdir(exist_ok=True)
    analyzed={}
    for case in plan['cases']:
        d=Path(case['directory']);out=d/'analysis-regularized-r1'
        require((out/'RESULT.json').exists(),'Missing repaired case: '+case['name'])
        run=load_json(d/'PIC_COMPLETE.json');require(run['exit_code']==0,'Unclosed PIC')
        radiation=load_json(out/'radiation.json');particles=load_json(out/'particles.json')
        analyzed[case['name']]=dict(directory=d,output=out,run=run,radiation=radiation,
                                   particles=particles,spectra=np.load(out/'spectra.npz'))
    subprocess.run([sys.executable,str(Path(__file__).parent/'compare_particle_phase.py'),str(root)],check=True)
    phases=load_json(root/'PARTICLE_PHASE_COMPARISON.json')
    reference={}
    for q in plan['genesis_references']:
        d=Path(q['directory'])
        require(load_json(d/'genesis.log.launch.json')['exit_code']==0,'Reference not closed')
        with h5py.File(d/'reference.out.h5','r') as f:
            z=f['Lattice/zplot'][:]+plan['parameters']['entrance_z_m']
            b=f['Beam/bunching'][:];phi=f['Beam/bunchingphase'][:]
            cur=f['Beam/current'][:].reshape(-1,b.shape[-1])[0];weights=cur/cur.sum()
            global_b=np.abs(np.sum(b*np.exp(1j*phi)*weights,axis=1))
            total_power=f['Field/power'][:].sum(axis=1)
            global_energy=f['Beam/Global/energy'][:]
        reference[q['name']]=dict(z_m=z.tolist(),global_bunching=global_b.tolist(),
           summed_slice_power_W=total_power.tolist(),
           pulse_energy_J=(total_power*q['sample']*plan['parameters']['resonant_wavelength_m']/299792458.).tolist(),
           mean_gamma_change=(global_energy-global_energy[0]).tolist())
    save(base/'GENESIS_GLOBAL_COMPARISON.json',reference)
    nominal=analyzed['momentum']
    WC=nominal['radiation'][-1]['results']['Hann']['E_only_band_J']
    bC=nominal['particles'][-1]['bunching']
    TA=analyzed['boris']['run']['seconds'];TC=nominal['run']['seconds']
    rows=[]
    for name,v in analyzed.items():
        rad=v['radiation'][-1];p=v['particles'][-1];seconds=v['run']['seconds']
        relative_phase=next((q for q in phases if q['plane']=='far' and
                           q['left']==name and q['right']=='momentum'),None)
        if name=='momentum':phase=0.
        elif relative_phase:phase=relative_phase['rms_phase_rad']
        else:
            reverse=next(q for q in phases if q['plane']=='far' and
                         q['left']=='momentum' and q['right']==name)
            phase=reverse['rms_phase_rad']
        spec=v['spectra'];reference_spec=nominal['spectra']
        selected=abs(reference_spec['far_Hann_energy_eV']-91.494)<3.
        s=spec['far_Hann_E_only_J_per_eV'][selected]
        sc=reference_spec['far_Hann_E_only_J_per_eV'][selected]
        spectrum_error=float(np.linalg.norm(s-sc)/np.linalg.norm(sc))
        growth_difference=[100*(p['results']['Hann']['E_only_band_J']/q['results']['Hann']['E_only_band_J']-1)
                           for p,q in zip(v['radiation'],nominal['radiation'])]
        rows.append(dict(method=name,PIC_elapsed_s=seconds,
          mean_elapsed_per_field_step_s=seconds/plan['geometry']['field_steps'],
          relative_cost_to_A=seconds/TA,Hann_annulus_band_J=rad['results']['Hann']['E_only_band_J'],
          rectangular_annulus_band_J=rad['results']['rectangular']['E_only_band_J'],
          Hann_relative_difference_from_C_percent=100*(rad['results']['Hann']['E_only_band_J']/WC-1),
          global_bunching=p['bunching'],bunching_relative_difference_from_C_percent=100*(p['bunching']/bC-1),
          phase_RMS_from_C_rad=phase,raw_particle_mean_kinetic_change_eV=p['mean_kinetic_change_eV'],
          exit_spectral_relative_L2_from_C=spectrum_error,
          maximum_abs_plane_band_difference_from_C_percent=max(abs(q) for q in growth_difference),
          spectral_centroid_eV=rad['results']['Hann']['band_centroid_eV']))
    with (base/'FEL_METHOD_COMPARISON.csv').open('w') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    refs=load_json(analyzed['boris']['output']/'GENESIS_APERTURE.json')
    with plt.rc_context({'font.size':9,'axes.spines.top':False,'axes.spines.right':False,
                         'pdf.fonttype':42,'ps.fonttype':42}):
        fig,ax=plt.subplots(2,2,figsize=(7.5,6.4),layout='constrained')
        for name,v in analyzed.items():
            z=[p['z_m'] for p in v['radiation']]
            energy=[p['results']['Hann']['E_only_band_J']*1e9 for p in v['radiation']]
            style=':' if name=='prescribed' else '-'
            ax[0,0].semilogy(z,energy,style+'o',ms=3,label=LABELS[name],color=COLORS[name])
            p=v['particles']
            ax[0,1].plot([q['z_m'] for q in p],[q['bunching'] for q in p],style+'o',
                         ms=3,label=LABELS[name],color=COLORS[name])
            spec=v['spectra'];ev=spec['far_Hann_energy_eV'];s=spec['far_Hann_E_only_J_per_eV']
            window=abs(ev-91.49)<4
            ax[1,0].plot(ev[window],s[window]*1e9,style,label=LABELS[name],color=COLORS[name])
        for tag,style in [('space-refined','--'),('weak',':')]:
            r=[q for q in refs if q['reference']==tag and q['window']=='Hann']
            ax[0,0].semilogy([q['z_m'] for q in r],[q['band_J']*1e9 for q in r],
                             style,color='black',lw=1.1,label='GENESIS '+('strong' if tag=='space-refined' else 'weak x I^2'))
        rr=reference['space-refined']
        ax[0,1].plot(rr['z_m'],rr['global_bunching'],'--',color='black',lw=1.1,label='GENESIS strong')
        ax[0,0].set(xlabel='Lab z (m)',ylabel='Hann band energy (nJ)',
                    title='90-130 um annulus; 3 eV half-band')
        ax[0,0].legend(fontsize=6.8,frameon=False)
        ax[0,1].set(xlabel='Lab z (m)',ylabel='Global bunching',title='Feedback and microbunching')
        ax[0,1].legend(fontsize=7,frameon=False)
        ax[1,0].axvline(91.49,color='black',lw=.8,ls='--')
        ax[1,0].set(xlabel='Photon energy (eV)',ylabel='Hann energy (nJ/eV)',title='Exit spectrum: same aperture')
        selected=rows[:3]
        ax[1,1].bar(range(3),[q['relative_cost_to_A'] for q in selected],
                     color=[COLORS[q['method']] for q in selected])
        ax[1,1].set_xticks(range(3),['A','B','C'])
        ax[1,1].set(ylabel='Total PIC cost / A',title='Same allocation; ledger and output on')
        for i,q in enumerate(selected):
            ax[1,1].text(i,q['relative_cost_to_A']+.03,f"{q['PIC_elapsed_s']/60:.1f} min",ha='center',fontsize=8)
        for ext in ['png','pdf','svg']:fig.savefig(figures/('FEL_three_methods.'+ext),dpi=500)
        plt.close(fig)
        fig,ax=plt.subplots(2,2,figsize=(7.5,6.2),layout='constrained')
        for panel,(folder,band) in enumerate([('synchrotron-r1',(8,12)),('single-electron-r1',(50,100))]):
            for name in ['boris','energy','momentum']:
                d=base/folder/name/'analysis'
                csvdata=np.genfromtxt(d/'aperture_spectra.csv',delimiter=',',names=True)
                x=csvdata['photon_eV'];m=(x>=band[0]*.5)&(x<=band[1]*1.5)
                ax[0,panel].semilogy(x[m],csvdata['PIC_J_per_eV'][m],color=COLORS[name],label=LABELS[name])
                ax[1,panel].semilogy(x[m],csvdata['PIC_J_per_eV'][m]/csvdata['total_LW_J_per_eV'][m],
                                     color=COLORS[name])
            ax[0,panel].semilogy(x[m],csvdata['total_LW_J_per_eV'][m],'--',color='black',label='Complete LW')
            ax[0,panel].set(xlabel='Photon energy (eV)',ylabel='Raw energy (J/eV)',
                             title='Gaussian short bunch' if panel==0 else 'Single electron')
            ax[0,panel].legend(fontsize=7,frameon=False)
            ax[1,panel].axhline(1,color='black',lw=.8,ls='--')
            ax[1,panel].set(xlabel='Photon energy (eV)',ylabel='PIC / complete LW')
            ax[1,panel].axvspan(*band,color='.8',alpha=.2)
        for ext in ['png','pdf','svg']:fig.savefig(figures/('synchrotron_three_methods.'+ext),dpi=500)
        plt.close(fig)
    save(base/'RESULT.json',dict(methods=rows,global_reference=reference,
      source_binary_sha256=plan['binary_sha256'],physical_acceptance=False,
      user_confirmed_engineering_target_percent=5.,
      conclusions=dict(original_explicit_is_not_equivalent_in_this_strong_feedback_case=True,
        energy_vs_longitudinal_Hann_percent=next(q for q in rows if q['method']=='energy')['Hann_relative_difference_from_C_percent'],
        energy_vs_longitudinal_bunching_percent=next(q for q in rows if q['method']=='energy')['bunching_relative_difference_from_C_percent'],
        energy_method_cost_saving_against_C_percent=100*(1-analyzed['energy']['run']['seconds']/TC)),
      caveats=['Finite off-beam radiation proxy, not a pure-radiation full-plane extraction',
               'GENESIS direct model omits transverse space charge and fringe radiation',
               'SR raw-field absolute calibration has not passed; no automatic production change']))
    text=['# 三种推进的辐射与成本对照','',
      '本轮完成了同一程序、同一输入、同一分配的三种推进与关闭反馈对照。**没有自动切换生产设置。**','',
      '## 冷束强反馈 FEL','',
      '选定共振附近 ±3 eV、离轴 90–130 μm 环形孔径，Hann 用相同能量归一；这是选定孔径的辐射代理量，不能称全平面纯辐射总量。','',
      '| 方法 | 完整 PIC 耗时/min | Hann 选频能量/nJ | 全局 bf | 较 C 的能量差 | 谱形 L2 差 | 较 C 的相位 RMS/rad |',
      '|---|---:|---:|---:|---:|---:|---:|']
    for q in rows:
        text.append(f"| {LABELS[q['method']]} | {q['PIC_elapsed_s']/60:.2f} | {q['Hann_annulus_band_J']*1e9:.6g} | {q['global_bunching']:.6f} | {q['Hann_relative_difference_from_C_percent']:+.3f}% | {100*q['exit_spectral_relative_L2_from_C']:.3f}% | {q['phase_RMS_from_C_rad']:.5g} |")
    energy=next(q for q in rows if q['method']=='energy')
    text+=['','![FEL 三种推进](figures/FEL_three_methods.png)','',
      f"能量迭代相对高阶配对节省本轮完整 PIC 时间 {100*(1-analyzed['energy']['run']['seconds']/TC):.2f}%；选定孔径能量差 {energy['Hann_relative_difference_from_C_percent']:+.3f}%，bf 差 {energy['bunching_relative_difference_from_C_percent']:+.3f}%。这些是同一测试的实测数，不外推为所有生产工况的固定比例。",'',
      f"用户确认以约 5% 为工程筛选目标。B 的出口谱形 L2 差为 {100*energy['exit_spectral_relative_L2_from_C']:.2f}%，沿程选频能量最大差为 {energy['maximum_abs_plane_band_difference_from_C_percent']:.2f}%。因此积分能量已经接近，但谱形尚未达到本轮 5% 目标，不能宣布 B 全面通过；bf 和相位差也分别保留。这里仍是以当前 C 为参照的三态比较，不是绝对物理精度认证。",'',
      '原显式版和迭代版在强反馈冷束中不等价。旧全量 19099/19174 目标 Hann 频带仅差约 1.413%，只支持那个工况的目标谱稳定，不能推广到这里。','',
      '### 与独立 FEL 模型的差异','',
      '同一出口、环形孔径、频带和 Hann 权重下，GENESIS 空间加细参考为 14.611 nJ；A/B/C 为 17.226/33.521/32.922 nJ。关闭反馈的 PIC 为 15.088 nJ，弱束 GENESIS 按电流平方放大后为 2.493 nJ。独立参考与 PIC 尚未达到约 5% 的绝对一致性；观察到反馈和谱峰不等于 FEL 物理验收。','',
      'GENESIS 横向网格由 241² 加细到 481² 时，完整脉冲能量改变约 5.2%，这里的环形孔径选频量改变约 7.1%；参考自身尚未收敛到 5% 以下。其时间采样加细影响约 0.26%，不能用这一项代替空间收敛。','',
      '本轮直接 GENESIS 模型缺横向空间电荷。独立冷束 rms 包络估计得到出口 σx/σy≈22.45/20.73 μm；PIC C 实测为 22.72/20.86 μm，GENESIS 为 20.00/18.31 μm。这个没有拟合 PIC 的估计说明模型范围差异有实质影响，但不能据此解释或校正全部辐射幅度差。公式和假设见理论文档。','',
      '## 基本辐射的绝对核对仍未通过','',
      '有限二极铁 raw 总场与完整 LW 有明显偏差；单电子排除了只由集体效应造成的解释。理论归一用 Larmor 角积分和零频端点恒等式独立核对，误差约 10⁻¹⁰。','',
      '| 工况与真实积分频带 | A：raw PIC / 完整 LW | B | C |',
      '|---|---:|---:|---:|',
      '| 高斯短束，8–12 eV | 97.247 | 102.245 | 36.164 |',
      '| 单电子，50–100 eV | 102.287 | 95.936 | 29.528 |','',
      '![同步辐射三种推进](figures/synchrotron_three_methods.png)','',
      '![空间尺度审计](figures/single_electron_spatial_audit.png)','',
      '单电子记录中约 45–50% 的空间谱能量位于横向 0.7 Nyquist 以上；同口径 LW 加速度参照约为 10⁻⁶。相位速度估计提出了数值切伦科夫候选，但小型加细/弱弯转结果并未证明这个单轴解释就是唯一或主要根因。初场/供场、采样和数值传播仍须分开定位。','',
      '固定理论模式投影只用于量测已知辐射形状在记录中的含量，不把参考幅度重新拟合到数据，也不把该投影作为通用近场分离或通过验收。','',
      'B 的额外校准中，单电子 50–100 eV 的 raw PIC/LW 比值由原网格的 95.94，变为横向 x 加细的 490.34、纵向和步长减半的 127.90；弱弯转对照也未恢复一致。它们没有给出单调收敛证据。原冻结小图的频带显示错误另列于勘误，数值积分掩码没有随之改变。','',
      '## 判断范围','',
      '较轻方法是否可用，应同时看与高阶版的差、相位和反馈、独立物理参照以及其模型不确定度。不能只看总束流动能下很小的余差；也不要求为了辐射目标无限追求 lab 粒子能量闭合。','',
      '当前直接 GENESIS 参考不包含横向空间电荷与 fringe 辐射，不能把它与 PIC 的全部差异都记成推进误差。严格账本与输出在所有三方法中保持开启，成本尚不是关闭诊断后的纯推进成本。','',
      '原冻结 FEL 后处理拒绝了全记录中的一个缺采样点。统一的共同迟滞时间区间避开该缺点；修复仅为独立后处理，未改原始数据、计算卡或推进。各面缺采样与时间格点误差仍逐面记录。','',
      '[理论公式与推导](THEORY.md) · [完整数字表](FEL_METHOD_COMPARISON.csv) · [总览与运行记录](README.md) · [显示与后处理勘误](METADATA_ERRATA.md)','']
    (base/'RESULT_REVIEW.md').write_text('\n'.join(text))
    print(base/'RESULT_REVIEW.md')

if __name__=='__main__':main()
