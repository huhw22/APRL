#!/usr/bin/env python3
"""Review compact postprocessing products, without reopening production HDF5 files."""
import argparse
import csv
import json
from pathlib import Path

import numpy as np


def review(job, out):
    out.mkdir(parents=True, exist_ok=False)
    if not (job / 'COMPLETE.json').exists():
        raise ValueError('The full analysis must finish before review')
    particle = json.loads((job / 'particles/particle-summary.json').read_text())
    field = json.loads((job / 'fields/lab_field_summary.json').read_text())
    ledger = json.loads((job / 'fields/boost_energy_summary.json').read_text())
    rows, checks, spectra = [], [], []
    for s in particle['planes']:
        with np.load(job / f"particles/particles-{s['label']}.npz") as a:
            w = a['slice_electrons']
            good = w > 0
            mean = float(np.average(a['slice_mean_energy_MeV'][good], weights=w[good]))
            variance = np.average(a['slice_sigma_energy_keV'][good]**2
                                 + (1000*(a['slice_mean_energy_MeV'][good]-mean))**2,
                                 weights=w[good])
            assert int(a['slice_records'].sum()) == s['records']
            assert np.isclose(mean, s['mean_total_energy_MeV'], rtol=0, atol=1e-9)
            assert np.isclose(np.sqrt(variance), s['sigma_energy_keV'], rtol=1e-10)
            assert np.isclose(w.sum(), s['represented_electrons'], rtol=1e-12)
            assert np.isclose(a['phase_weighted_electrons'].sum(), w.sum(), rtol=1e-12)
            missing = float(s['represented_electrons'] - a['loss_weighted_electrons'].sum())
            # The original frozen run can have a rounding-only endpoint omission.
            # Preserve it explicitly; do not invent which edge the missing count belongs to.
            checks.append(dict(plane=s['label'], slice_and_phase_charge_conserved=True,
                               independent_energy_moments_match=True,
                               loss_histogram_missing_represented_electrons=missing,
                               loss_histogram_missing_fraction=missing/s['represented_electrons']))
            energy, bf = a['spectrum_energy_eV'].copy(), a['spectrum_b'].copy()
            band = np.flatnonzero((energy >= 60) & (energy <= 120))
            peak = band[np.argmax(bf[band])]
            spectra.append((s['z_m'], energy, bf))
        b = s['bunching']
        rows.append(dict(z_m=s['z_m'], mean_total_energy_MeV=s['mean_total_energy_MeV'],
                         projected_rms_MeV=s['sigma_energy_keV']/1000,
                         linear_detrended_rms_MeV=s['linear_detrended_sigma_energy_keV']/1000,
                         recorded_entry_minus_plane_MeV=s['mean_ID_matched_loss_keV']/1000,
                         exact_b1_fixed_wavelength=b['exact_b_h1_h2_h3_fixed_entry_resonance'][0],
                         exact_b2_fixed_wavelength=b['exact_b_h1_h2_h3_fixed_entry_resonance'][1],
                         exact_b3_fixed_wavelength=b['exact_b_h1_h2_h3_fixed_entry_resonance'][2],
                         spectral_peak_60_120_eV=float(energy[peak]),
                         spectral_peak_b_60_120_eV=float(bf[peak]),
                         physical_aperture_fraction=s['physical_xy_aperture_charge_fraction'],
                         internal_plane_interpolation_risk=s['label'] not in ('entry', 'exit')))
    with (out/'plane-comparison.csv').open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)
    result = dict(source_job=str(job.resolve()), full_scan_script_commit='74ef494',
                  fixed_bunching_wavelength_nm=particle['planes'][0]['bunching']['fixed_entry_resonant_wavelength_nm'],
                  planes=rows, checks=checks,
                  field=field, boost_ledger=ledger,
                  diagnostic_audit=str((job.parent.parent/'diagnostic-interpolation-audit.txt').resolve()),
                  scientific_claim_limits=[
                      'Middle-plane energy and spread are recorded values, affected by a documented interpolation-bias mechanism; no numerical correction is possible from these crossing records alone.',
                      'Bunching uses interpolated arrival times; this analysis does not establish time-step convergence.',
                      'Raw field energy is finite-window and finite-aperture; it is not separated complete-pulse radiation.',
                      'The boosted ledger has closure_valid=0 and cannot certify energy closure.',
                      'Tiny omissions in original loss histograms are retained and reported; exact mean loss is unaffected. Future analyzer revision fixes endpoint arithmetic.'
                  ])
    (out/'review-summary.json').write_text(json.dumps(result, indent=2, ensure_ascii=False, allow_nan=False)+'\n')
    plot(out, rows, spectra, result['fixed_bunching_wavelength_nm'],
         particle['planes'][0]['shot_noise_b_rms_reference'])
    plot_field(out, job, field)
    lines = [
        '历史 ELEGANT 4 m 仿真结果：已审阅结论',
        '原始历史结果已完整复制并校验，未修改。全量后处理作业18248正常完成；本审阅只读取小型派生数据。',
        '分析对象为40µm裁剪束段，41,732,775个单位权重电子；五面同ID、source_id、电荷及质量验证通过。',
        '',
        '粒子固定lab面（下列中间面数值含探测插值风险，不是已校正的物理能损/能散曲线）：',
        'z(m)   均能(MeV)   投影RMS(MeV)   入口减该面(MeV)   BF(lambda=13.5532117nm)',
    ]
    for r in rows:
        lines.append(f"{r['z_m']:5.2f}  {r['mean_total_energy_MeV']:12.6f}  {r['projected_rms_MeV']:12.6f}  {r['recorded_entry_minus_plane_MeV']:14.6f}  {r['exact_b1_fixed_wavelength']:.8f}")
    first, last = particle['planes'][0], particle['planes'][-1]
    lines += [
        f"出口记录净能降{last['mean_ID_matched_loss_keV']:.6f}keV/电子，整束合计{last['total_ID_matched_loss_J']*1e6:.6f}µJ。",
        f"出口投影RMS较入口增加{last['sigma_energy_keV']-first['sigma_energy_keV']:.6f}keV。",
        '旧探测器在完整Maxwell步端点插值boost动量，再重算gamma并变到lab；存在确定的向下插值偏差机制，量级与1m/3m异常相符。',
        '该证据尚不能把全部异常定量归因给插值，也不能从保存的过面记录无损修复；不可把中间面的数MeV降能称累计辐射损失。',
        '入口已有显著BF谱峰和周期结构；固定波长BF与移动谱峰回答不同问题，不把全段称为从零噪声产生的微聚束。',
        '1fs切片BF含矩形窗/束流包络和边缘低粒子数效应；去线性趋势RMS也不等于不可逆无关联能散。',
        '',
        '场面：lab z=5.370768546m，记录153.547225fs，160×160格，去CPML物理孔径140×140格。',
        f"原始场峰值 |E|={field['peak_E_magnitude_V_per_m']/1e9:.6f}GV/m，|B|={field['peak_B_magnitude_T']:.6f}T，主要为Ex与By。",
        f"全记录孔径净过面能量={field['recorded_aperture']['net_energy_J']*1e6:.6f}µJ；物理孔径={field['physical_aperture']['net_energy_J']*1e6:.6f}µJ。",
        f"物理孔径峰值净功率={field['physical_aperture']['peak_net_power_W']/1e6:.6f}MW，末态仍为{field['physical_aperture']['final_net_power_W']/1e6:.6f}MW。",
        f"最后一帧最大场强仍为全局峰值的{field['last_frame_peak_E_fraction_of_global_peak']*100:.2f}%。",
        'raw场包含粒子近场；E/B含Yee时间交错。上述能量是已记录时间窗和孔径内的原生Poynting积分，不是完整脉冲纯辐射能，更不能除以束流损失直接称辐射效率。',
        '',
        'boost系账本（共同boost时间片，与固定lab过面统计分开）：',
        f"粒子净动能损失={ledger['particle_kinetic_loss_J']*1e9:.7f}nJ；域内末态场能={ledger['final']['field_energy_interior_J']*1e9:.7f}nJ；域内场能增量={ledger['delta_interior_field_J']*1e9:.7f}nJ。",
        f"累计净流出场能={ledger['final']['outward_field_energy_J']*1e9:.7f}nJ；规定元件做功={ledger['final']['prescribed_work_J']*1e9:.7f}nJ；残差={ledger['final']['balance_residual_J']*1e9:.7f}nJ。",
        '规定磁场在boost系有电场功，不能省略。closure_valid=0；本次只核查账本代数一致，不能宣称严格整体能量守恒。',
        '主粒子面全部或几乎全部在物理横向孔径内；约64.13%框外是更晚共同boost末态的比例，不等于波荡器内部64.13%粒子全程未受场。',
        '',
        '产品说明：同目录含plane-comparison.csv、review-summary.json及审阅版科学图；完整原始派生数据和图在上一层particles/、fields/。',
        '原作业少数loss直方图有极端边界舍入漏计，逐面缺额见review-summary.json；精确均能、能散、同ID平均能损不受影响。新脚本已统一计算顺序并扩展外侧边界供今后使用；本次未伪造补回历史直方图计数。',
    ]
    (out/'analysis-notes.txt').write_text('\n'.join(lines)+'\n')
    return result


def plot(out, rows, spectra, wavelength, noise):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size':10, 'axes.grid':True, 'grid.alpha':.22})
    z=np.array([r['z_m'] for r in rows])
    fig, axes=plt.subplots(1,3,figsize=(15,4.6),constrained_layout=True)
    for ax,key,label in ((axes[0],'recorded_entry_minus_plane_MeV','Recorded entry minus plane energy [MeV]'),
                         (axes[1],'projected_rms_MeV','Recorded projected RMS spread [MeV]')):
        y=np.array([r[key] for r in rows]); ax.plot(z,y,color='0.6',lw=1)
        ax.scatter(z[[0,-1]],y[[0,-1]],color='#2166ac',s=55,label='Entrance / exit')
        ax.scatter(z[1:4],y[1:4],color='#b2182b',marker='x',s=65,label='Internal planes: interpolation risk')
        ax.set_ylabel(label);ax.set_xlabel('Fixed laboratory z [m]')
    axes[0].legend(fontsize=8)
    for h in range(1,4):
        axes[2].semilogy(z,[r[f'exact_b{h}_fixed_wavelength'] for r in rows],'o-',label=f'h={h}')
    axes[2].axhline(noise,color='gray',ls=':',label='Independent-phase RMS reference')
    axes[2].set(xlabel='Fixed laboratory z [m]',ylabel='Bunching factor from recorded arrival times')
    axes[2].legend(fontsize=8)
    fig.suptitle(f'Recorded beam diagnostics | fixed reference wavelength {wavelength:.6f} nm\nInternal energy/spread are not a validated physical loss curve; arrival-time convergence is untested.',fontsize=12)
    for ext in ('png','pdf'):fig.savefig(out/f'particle-recorded-overview.{ext}',dpi=170,bbox_inches='tight')
    plt.close(fig)
    fig,ax=plt.subplots(figsize=(11,4.5),constrained_layout=True)
    for z,e,b in spectra:
        mask=(e>=60)&(e<=120)
        ax.semilogy(e[mask],b[mask],lw=.8,label=f'z={z:g} m')
    ax.axvline(1239.841984332/wavelength,color='gray',ls=':',label='Fixed reference')
    ax.axhline(noise,color='gray',ls='--',alpha=.6)
    ax.set(xlabel='Photon-equivalent frequency hf [eV]',ylabel='Full-bunch |b(f)|',
           title='Recorded arrival-time form factor | initial structure is already present')
    ax.legend(ncol=3,fontsize=9)
    for ext in ('png','pdf'):fig.savefig(out/f'bunching-fundamental-band.{ext}',dpi=170,bbox_inches='tight')
    plt.close(fig)


def plot_field(out, job, field):
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    with np.load(job/'fields/lab_field_maps.npz') as a:
        peak=a['peak_E_magnitude_V_per_m']/1e9
        x,y=a['x_m']*1e6,a['y_m']*1e6
    with np.load(job/'fields/lab_field_timeseries.npz') as a:
        t=a['time_from_first_fs'];power=a['physical_net_power_W']/1e6
    fig,axes=plt.subplots(1,2,figsize=(11.5,4.8),constrained_layout=True)
    dx,dy=field['dx_m']*1e6,field['dy_m']*1e6
    im=axes[0].imshow(peak,origin='lower',extent=[x[0]-dx/2,x[-1]+dx/2,y[0]-dy/2,y[-1]+dy/2],cmap='magma')
    axes[0].add_patch(Rectangle((-350,-350),700,700,fill=False,edgecolor='cyan',lw=.8))
    axes[0].set(xlabel='x [um]',ylabel='y [um]',title='Peak |E| over recorded time; cyan: physical aperture')
    axes[0].grid(False);fig.colorbar(im,ax=axes[0],label='GV/m',shrink=.85)
    axes[1].plot(t,power,lw=.7,color='#2166ac')
    axes[1].set(xlabel='Time from first sample [fs]',ylabel='Net power through physical aperture [MW]',
                title=f"Recorded energy {field['physical_aperture']['net_energy_J']*1e6:.3f} uJ; tail remains nonzero")
    fig.suptitle('Raw laboratory field at z=5.37077 m | particle near field retained; finite recorded window',fontsize=12)
    for ext in ('png','pdf'):fig.savefig(out/f'lab-field-overview.{ext}',dpi=170,bbox_inches='tight')
    plt.close(fig)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--job-dir',type=Path,required=True)
    p.add_argument('--output-dir',type=Path,required=True)
    a=p.parse_args();review(a.job_dir,a.output_dir)
