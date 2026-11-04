"""Validate and document the phase-corrected radiation-only comparison."""
import ast
import hashlib
import json
import sys
from pathlib import Path

import numpy as np


def main(workspace):
    w = Path(workspace).resolve()
    top = w / 'results/sparc_sase_comparison_20261103'
    rev = top / 'explicit_review_20261104'
    root = top / 'compact-production-r1'
    r = rev / 'radiation_phase_r1'
    old = rev / 'radiation_aligned_r1'
    repo = w / 'source/aprl'
    read = lambda p: json.loads(Path(p).read_text())
    sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()

    def write(p, obj):
        Path(p).write_text(json.dumps(obj, indent=2, allow_nan=False) + '\n')

    assert read(r / 'COMPLETE.json')['complete']
    assert read(old / 'RADIATION_RESULT.json')['cross_code_comparison_withdrawn']
    phase = read(r / 'PHASE_READER_CONTROL.json')
    assert phase['old_rejected'] and phase['new_nominal_closure_pass']
    pic, gen = read(r / 'pic.json'), read(r / 'genesis.json')
    result = read(r / 'RADIATION_RESULT.json')
    band = result['band_eV']
    checks = []
    assert len(pic) == len(gen) == 8
    for p, g in zip(pic, gen):
        assert p['plane'] == g['plane']
        with np.load(r / (p['plane'] + '-spectra.npz')) as a, np.load(r / (p['plane'] + '-genesis.npz')) as b:
            T = p['fft_samples'] * p['dt_s']
            for label in ['rect', 'hann']:
                x = a['energy_eV']; m = (x >= band[0]) & (x <= band[1])
                energy = a[label + '_native_flux_J_per_eV'][m].sum() * (x[1] - x[0])
                y = b['energy_eV']; m = (y >= band[0]) & (y <= band[1])
                reference = b[label + '_J_per_eV'][m].sum() * (y[1] - y[0])
                assert abs(energy / p['results'][label]['native_Poynting_band_J'] - 1) < 1e-12
                assert abs(reference / g['results'][label]['band_J'] - 1) < 1e-11
                ratio = energy / reference
                assert abs(ratio / g['results'][label]['PIC_native_flux_over_GENESIS'] - 1) < 1e-11
                checks.append(dict(plane=p['plane'], window=label, ratio=float(ratio),
                                   PIC_record_mean_W=float(energy / T), GENESIS_record_mean_W=float(reference / T)))
        assert sha(r / (p['plane'] + '-spectra.npz')) == sha(old / (p['plane'] + '-spectra.npz'))
    with np.load(r / 'p1-genesis.npz') as a, np.load(old / 'p1-genesis.npz') as b:
        assert a.files == b.files
        assert all(np.array_equal(a[k], b[k]) for k in a.files)

    for manifest in [root / 'FROZEN.json', rev / 'ALIGNED_FROZEN.json', rev / 'PHASE_FROZEN.json']:
        assert all(sha(p) == h for p, h in read(manifest).items()), str(manifest)
    protected = read(w / 'results/group_meeting_materials_20261103/REFERENCE_AUDIT.json')['implementation_and_protected_files_sha256']
    assert all(sha(p) == h for p, h in protected.items())
    source_files = [repo / 'verification/analysis/sparc_sase' / name for name in
                    ['aligned_reference.py', 'phase_reference.py', 'phase_review.py', 'review_result.py']]
    source_files += [repo / 'verification/performance/cases/sparc_sase' / name for name in ['common.py', 'run_reference.py']]
    for p in source_files:
        ast.parse(p.read_text(), filename=str(p))
    gain = read(r / 'GAIN_AUDIT.json')
    after = gain['cases']['after']
    late = 'late_segments_4_to_6'
    lg = after['fits']['rectangular_band_J'][late]['magnet']['power_gain_length_m']
    local = after['local_power_dense_fits']
    peak_lg = local['fixed_initial_peak_current_slice'][late]['magnet']['power_gain_length_m']
    central_lg = local['fixed_central_current_ge_0.95_of_input_peak'][late]['magnet']['power_gain_length_m']
    validation = dict(arithmetic_validation_passed=True, physical_acceptance=False,
      reference_phase_error_corrected=True, old_comparison_withdrawn=True,
      phase_chosen_from_nominal_ODE_not_power=True, first_segment_reference_unchanged_bitwise=True,
      original_PIC_spectral_caches_unchanged=True, frozen_inputs_unchanged=True, protected_files_unchanged=True,
      max_PIC_parseval_relative=max(p['results'][k]['parseval_flux_relative'] for p in pic for k in ['rect', 'hann']),
      max_GENESIS_projection_normalization_relative=max(g['normalization_projection_relative_error'] for g in gen),
      minimum_GENESIS_time_coverage=min(g['captured_pulse_energy_fraction'] for g in gen),
      maximum_empirical_missing_sample_band_energy_bound=max(p['results'][k]['empirical_one_missing_sample_energy_bound'] for p in pic for k in ['rect', 'hann']),
      arithmetic_checks=checks, images_visually_inspected=['radiation_power.png', 'explicit_spectra.png'],
      source_sha256={str(p): sha(p) for p in source_files},
      scope='Radiation only; fixed-record mean powers; PIC bound-field contribution retained; single realization, no precision certificate')
    write(r / 'VALIDATION.json', validation)
    far = [x for x in checks if x['plane'] == 'far']
    a, b = far
    text = rf'''# GENESIS参考相位核查与修正（2026-11-04）

旧GENESIS参考卡存在段间相位错误，先前“显式版Hann频带功率高约3.8倍”的比较已撤回。修正后的参考已由19300完成，仅重算GENESIS辐射，原PIC计算、输入与频谱缓存保持。

## 确定的错误及修正

本机GENESIS4.6.4的BeamSolver在漂移段将零波荡器波数替换为 $k_s/(2\gamma_0^2)$，已经包含参考载波补偿。旧参考卡遗漏这项，又按完整物理漂移时延补偿一次，每个间隙造成约−2.662 rad（−152.5°）的名义入口相位误配。问题在参考卡构造，不是已证实的GENESIS推进错误。

对零横向偏移、固定参考能量、关闭辐射反馈的名义轨迹，实际相位方程为

$$
\frac{{d\theta}}{{dz}}=\widetilde k_u-k_s\frac{{1+a_w^2}}{{u_z(\gamma_0+u_z)}},\qquad
u_z=\sqrt{{\gamma_0^2-1-a_w^2}},\qquad
\widetilde k_u=\begin{{cases}}k_u,& k_u\ne0,\\ k_s/(2\gamma_0^2),& k_u=0.\end{{cases}}
$$

GENESIS实际施加的移相为 $-\phi$。因此按每个“磁段+间隙”的积分取 $\phi=\operatorname{{wrap}}\int(d\theta/dz)\,dz$，得到约−0.9577053 rad。只修改五个移相参数，未用输出功率拟合；修正后的名义闭合余差约2.45×10⁻¹⁶ rad。粒子能散、横向轨迹引起的物理相差仍由GENESIS演化，没有强制所有粒子同相。

输入电荷80.675 pC、峰值流强53.275 A、151.9 MeV动能、0.02%相对能散及2.5/2.9 μm归一化发射度的单位均已独立复核。第一个移相器之前的第一磁段辐射缓存逐数组完全一致；其后晶格除相位外的z、dz、aw、ku、qf、slippage逐位保持。旧参考、图像及数值保留并标注撤回。

## 与理论尺度的关系

功率增益长度定义为 $P\propto\exp(s/L_g)$，这里 $s$ 是累计有磁长度；场幅增益长度为其两倍，间隙不算入 $s$。文献的0.67 m是匹配束流的Ming Xie近似尺度，不是当前短束整脉冲功率的严格解析答案。

|末第4–6段的固定统计口径|旧参考|修正后|
|---|---:|---:|
|整束矩形频带能量|1.005 m|{lg:.3f} m|
|固定入口峰流强单切片|1.184 m|{peak_lg:.3f} m|
|固定入口中央区域（I≥95%峰值）平均功率|见GAIN_AUDIT.json|{central_lg:.3f} m|

中央单切片的局部增长已接近0.67 m，但中央有限区域及整短束仍有明显差异。不能用单片结果替代整束验收，也不能把0.83 m与0.67 m的差额直接定性为程序误差。第2–5段预先固定区间、末段区间、精确段出口和密集输出拟合全部保留在GAIN_AUDIT.json；未按接近理论程度筛选区间。原始文献束长12.67 ps、450 pC，本次按用户要求缩短为2.271 ps、80.675 pC并保持峰值流强，滑移及有限脉冲影响不能直接沿用原实验整脉冲曲线。

[Campbell等，FEL2019 TUP050](https://proceedings.jacow.org/fel2019/papers/tup050.pdf)给出SPARC参数与0.67 m估计；本报告采用其局部尺度，不声称复现原实验完整工况。

## 与显式版的更新比较

|15.246 m面，同一频带与孔径|显式PIC|修正GENESIS|PIC/GENESIS|
|---|---:|---:|---:|
|矩形窗记录平均频带功率|{a['PIC_record_mean_W']/1000:.3f} kW|{a['GENESIS_record_mean_W']/1000:.3f} kW|{a['ratio']:.4f}|
|能量归一Hann记录平均频带功率|{b['PIC_record_mean_W']/1000:.3f} kW|{b['GENESIS_record_mean_W']/1000:.3f} kW|{b['ratio']:.4f}|

频带{band[0]:.6f}–{band[1]:.6f} eV，物理孔径2.4×2.4 mm²，记录时长3.930766 ps。功率=频带能量/记录时长，不能与局部瞬时峰功率混用。Hann使用相同物理时刻的权重。没有幅度匹配、近场扣除或粒子能量对比。PIC仍是原始自洽总E/B，GENESIS为辐射包络；单种子、粗网格及重建晶格的物理差异仍存在。

目前不能给出显式推进“辐射精度达5%”的验收结论；参考相位错误已经解决，但完整的参考与PIC辐射验收仍未完成。

![修正后的功率比较](radiation_power.png)

![修正后的四宫格能谱](explicit_spectra.png)

可复核数据：RADIATION_RESULT.json、radiation_power.csv、VALIDATION.json、GAIN_AUDIT.json及各面NPZ。相位推导见../genesis-phase-r1/PHASE_DERIVATION.json，实际晶格闭合见../genesis-phase-r1/PHASE_CLOSURE_VERIFIED.json，旧卡拒绝/新卡接受对照见PHASE_READER_CONTROL.json。
'''
    (r / 'README.md').write_text(text)
    write(rev / 'FINAL.json', dict(complete=True, radiation_only=True, result_directory=str(r),
      analysis_job=19300, cancelled_unneeded_particle_job=19298, physical_acceptance=False,
      reference_phase_error_corrected=True, prior_cross_code_comparison_withdrawn=True,
      prior_result_directory=str(old)))
    jobs = read(top / 'JOB_LINKS.json')
    jobs['reference_phase_correction'].update(state='COMPLETED', exit_code='0:0',
      result_directory=str(r), nominal_phase_closure_verified=True, physical_acceptance=False)
    jobs['current_radiation_result'] = str(r / 'README.md')
    write(top / 'JOB_LINKS.json', jobs)
    audit = read(rev / 'REFERENCE_AUDIT_20261104.json')
    audit.update(corrected_reference_pending=False, corrected_reference_job=19300,
      corrected_gain_length_rectangular_late_m=lg, corrected_fixed_peak_slice_gain_length_m=peak_lg,
      corrected_central_95pct_gain_length_m=central_lg, physical_acceptance=False)
    write(rev / 'REFERENCE_AUDIT_20261104.json', audit)
    print(json.dumps({k: v for k, v in validation.items() if isinstance(v, (int, float, bool))}, indent=2))


if __name__ == '__main__':
    main(sys.argv[1])
