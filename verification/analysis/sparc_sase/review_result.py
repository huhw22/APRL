"""Recompute delivered radiation comparisons and preserve scope limitations."""
import sys,json,hashlib
from pathlib import Path
import numpy as np
def main(workspace):
    w=Path(workspace);top=w/'results/sparc_sase_comparison_20261103'
    rev=top/'explicit_review_20261104';r=rev/'radiation_aligned_r1';root=top/'compact-production-r1'
    L=lambda p:json.loads(Path(p).read_text())
    H=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
    assert L(r/'COMPLETE.json')['complete']
    pic=L(r/'pic.json');gen=L(r/'genesis.json');result=L(r/'RADIATION_RESULT.json')
    if result.get('cross_code_comparison_withdrawn'):
        raise RuntimeError('Reference comparison withdrawn; use phase_review.py for the corrected reference')
    rows=result['plane_rows'];band=result['band_eV'];checks=[]
    for p,g in zip(pic,gen):
        a=np.load(r/(p['plane']+'-spectra.npz'));b=np.load(r/(p['plane']+'-genesis.npz'))
        T=p['fft_samples']*p['dt_s']
        for label in ['rect','hann']:
            x=a['energy_eV'];m=(x>=band[0])&(x<=band[1])
            v=a[label+'_native_flux_J_per_eV'][m].sum()*(x[1]-x[0])
            assert abs(v/p['results'][label]['native_Poynting_band_J']-1)<1e-12
            x=b['energy_eV'];m=(x>=band[0])&(x<=band[1])
            v2=b[label+'_J_per_eV'][m].sum()*(x[1]-x[0])
            assert abs(v2/g['results'][label]['band_J']-1)<1e-11
            ratio=v/v2
            assert abs(ratio/g['results'][label]['PIC_native_flux_over_GENESIS']-1)<1e-11
            checks.append(dict(plane=p['plane'],window=label,ratio=float(ratio),
              PIC_mean_W=float(v/T),GENESIS_mean_W=float(v2/T)))
    assert all(H(p)==h for p,h in L(root/'FROZEN.json').items())
    assert all(H(p)==h for p,h in L(rev/'ALIGNED_FROZEN.json').items())
    protected=L(w/'results/group_meeting_materials_20261103/REFERENCE_AUDIT.json')['implementation_and_protected_files_sha256']
    assert all(H(p)==h for p,h in protected.items())
    validation=dict(ready_with_caveats=True,physical_acceptance=False,arithmetic_checks=checks,
      max_PIC_parseval_relative=max(p['results'][k]['parseval_flux_relative'] for p in pic for k in ['rect','hann']),
      max_GENESIS_projection_normalization_relative=max(g['normalization_projection_relative_error'] for g in gen),
      minimum_GENESIS_time_coverage=min(g['captured_pulse_energy_fraction'] for g in gen),
      maximum_empirical_missing_sample_band_energy_bound=max(p['results'][k]['empirical_one_missing_sample_energy_bound'] for p in pic for k in ['rect','hann']),
      original_frozen_inputs_unchanged=True,aligned_reference_inputs_unchanged=True,protected_files_unchanged=True,
      images_visually_inspected=['radiation_power.png','explicit_spectra.png'],
      scope='Radiation power only; record average, not instantaneous peak; raw PIC fields retain bound-field contributions')
    (r/'VALIDATION.json').write_text(json.dumps(validation,indent=2)+'\n')
    f=rows[-1];T=pic[-1]['fft_samples']*pic[-1]['dt_s']
    ph=pic[-1]['results']['hann']['native_Poynting_band_J']/T
    gh=gen[-1]['results']['hann']['band_J']/T
    text=f"""# 显式版辐射功率对照（2026-11-04）

当前功率对照未达到约5%的筛选目标。15.246 m面，矩形窗显式/GENESIS为{f['PIC_native_over_GENESIS']:.4f}，同物理Hann窗为{f['PIC_Hann_over_GENESIS']:.4f}。谱峰位置相近：显式{f['PIC_peak_eV']:.6f} eV，GENESIS{f['GENESIS_peak_eV']:.6f} eV；峰位接近不能替代功率验收。

|15.246 m面|显式版|GENESIS|比值|
|---|---:|---:|---:|
|矩形窗平均频带功率|{f['PIC_native_record_mean_W']/1000:.3f} kW|{f['GENESIS_record_mean_W']/1000:.3f} kW|{f['PIC_native_over_GENESIS']:.4f}|
|同物理Hann窗平均频带功率|{ph/1000:.3f} kW|{gh/1000:.3f} kW|{f['PIC_Hann_over_GENESIS']:.4f}|

![功率随位置变化](radiation_power.png)

![出口四宫格能谱](explicit_spectra.png)

频带{band[0]:.6f}–{band[1]:.6f} eV；孔径2.4×2.4 mm²，排除10格横向CPML。平均功率=频带能量/同一{f['record_ps']:.6f} ps记录时长，不是峰值功率。GENESIS按PIC各面的同一物理时间权重处理，没有幅度拟合或电流缩放。

Hann口径下，磁场后0.15/0.30 m的倍率分别为{rows[-2]['PIC_Hann_over_GENESIS']:.4f}/{rows[-1]['PIC_Hann_over_GENESIS']:.4f}，较矩形窗口径稳定。两种窗口涉及实际波形权重及记录边缘，不能单凭加窗差值计算近场能量。PIC仍是完整总E/B，包含束缚场及交叉项；GENESIS为辐射包络。因此目前支持“现有读数尚未与参考对上”，不能将全部差异归因于真实辐射、推进误差或某个唯一机制。B/C尚待同口径比较。

## 辐射参照的坐标核对

GENESIS正滑移使光场向较大切片编号前进，切片按束团尾到头增加；固定lab面的到达时间按头到尾增加。早先参考的活跃切片顺序未对应这项反向关系。本次只反转活跃切片组，各组位置、动量、光学相位和流强数据逐位核对保留，再计算与PIC对应的辐射参考。没有按辐射大小调参，原参考保留。GENESIS物理时刻按切片反向及名义飞行时延转换，复包络频谱符号相应处理。

已取消用户不需要的GENESIS出口粒子任务19298，未闭合结果未使用。后处理19297被独立单位核验拦截：GENESIS内部常数510999.06 eV不能用511000做严格恒等核验。最终19299使用本机常数核对，复用已闭合的PIC频谱，没有重跑PIC。原始场、计算卡及运行中B/C保持。

## 可复核数据

radiation_power.csv、RADIATION_RESULT.json、VALIDATION.json及每面频谱NPZ保存完整读数。逐像素FFT后才对复数交叉谱求和，没有先叠加全平面场。每面最多一个早期缺样使用八邻点插值，并保存留一诊断。全部FFT通过时间域/频域Parseval，GENESIS复场与独立int_xy投影核对单位，交付功率比再从谱缓存积分复核。单种子和粗网格不能作为绝对物理精度证明。
"""
    (r/'README.md').write_text(text)
    (rev/'FINAL.json').write_text(json.dumps(dict(complete=True,radiation_only=True,
      result_directory=str(r),analysis_job=19299,cancelled_unneeded_particle_job=19298,
      physical_acceptance=False),indent=2)+'\n')
    print(json.dumps({k:v for k,v in validation.items() if isinstance(v,(int,float,bool))},indent=2))
if __name__=='__main__':main(sys.argv[1])
