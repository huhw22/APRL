"""Add native, phase-corrected GENESIS bunching to the cached PIC comparison."""
import csv
import hashlib
import json
import sys
from pathlib import Path

import h5py
import numpy as np
from scipy.constants import c

from analyze_explicit import particle_read
from bunching_growth import read, save


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(4 * 1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def combine(current, moment, width):
    starts = np.arange(0, current.shape[1], width)
    charge = np.add.reduceat(current, starts, axis=1)
    complex_sum = np.add.reduceat(current * moment, starts, axis=1)
    local = complex_sum / charge
    total = current.sum(axis=1)
    rms = np.sqrt((charge * abs(local) ** 2).sum(axis=1) / total)
    return local, charge, rms


def entrance_check(top, wavelength, npart, expected_n):
    path = top / 'compact-production-r1/boris/detectors/entry-particles.h5'
    data = particle_read(path)
    assert np.array_equal(data['source_id'], np.arange(1, expected_n + 1))
    charge = abs(data['charge_C'])
    time = data['time_s'].astype(np.longdouble)
    angle = np.asarray(np.remainder(
        (time - time[len(time) // 2]) * np.longdouble(2 * np.pi * c / wavelength),
        np.longdouble(2 * np.pi)), float)
    phase = np.cos(angle) + 1j * np.sin(angle)
    group = (data['source_id'].astype(int) - 1) // (10 * npart)
    weights = np.bincount(group, weights=charge)
    moment = (np.bincount(group, weights=charge * phase.real)
              + 1j * np.bincount(group, weights=charge * phase.imag))
    return dict(global_b=float(abs(moment.sum()) / charge.sum()),
                local_rms_10lambda=float(np.sqrt(np.sum(abs(moment) ** 2 / weights)
                                                       / charge.sum())),
                particles=len(data), full_input=True)


def render(out, pic, genesis, wavelength):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 12, 'axes.spines.top': False,
                         'axes.spines.right': False, 'lines.linewidth': 1.8})
    styles = [('boris', 'Explicit', '#24648c', 'o-'),
              ('energy', 'Energy coupled', '#8e4b75', '^--')]
    fig, axes = plt.subplots(2, 2, figsize=(12.4, 8.3), layout='constrained')
    configs = [('design', 'global_b', f'(a) Global BF: {wavelength * 1e9:.3f} nm'),
               ('design', 'local_rms_10lambda', '(b) Local RMS BF: design wavelength'),
               ('explicit_peak', 'local_rms_10lambda',
                f"(c) PIC only: {pic['wavelengths_nm']['explicit_peak']:.3f} nm"),
               ('energy_peak', 'local_rms_10lambda',
                f"(d) PIC only: {pic['wavelengths_nm']['energy_peak']:.3f} nm")]
    for ax, (wave, key, title) in zip(axes.flat, configs):
        for case, label, color, style in styles:
            rows = [r for r in pic['rows'] if r['case'] == case
                    and r['wavelength_label'] == wave]
            ax.semilogy([r['z_m'] for r in rows], [r[key] for r in rows],
                        style, color=color, label=label)
        if wave == 'design':
            ax.semilogy([r['z_m'] for r in genesis], [r[key] for r in genesis],
                        '.-', color='#39815b', markersize=4, label='GENESIS', zorder=4)
        ax.set_title(title)
        ax.set_xlabel('Laboratory z (m)')
        ax.set_ylabel('Bunching factor')
        ax.set_ylim(1e-5, 1)
        ax.set_xlim(-.7, 16)
        ax.grid(alpha=.2)
        ax.legend()
    fig.suptitle('SPARC: PIC and phase-corrected GENESIS bunching growth', fontsize=15)
    fig.supxlabel('PIC: fixed 99.9485% charge subset; GENESIS: full input. '
                  'Local BF: complex averaging over 10 input optical slices.', fontsize=9)
    fig.savefig(out / 'bunching_growth_with_genesis.png', dpi=210,
                bbox_inches='tight', pad_inches=.12)
    plt.close(fig)


def document(out, top, pic, rows, entry, wavelength):
    tail = rows[-1]
    subset = next(r for r in pic['rows'] if r['case'] == 'boris'
                  and r['plane'] == 'entry-particles' and r['wavelength_label'] == 'design')
    text = f'''# 加入GENESIS的BF沿程增长对比

本图加入相位修正后的GENESIS参考（作业19300），只读取已完成的参考数据；PIC原曲线和原结果保持。

![BF增长对比](bunching_growth_with_genesis.png)

绿色为GENESIS的原生117个z记录，覆盖0–{tail['z_m']:.6f} m。只有设计波长{wavelength * 1e9:.6f} nm的基频复BF可用，因此只加入上排两个面板；下排两个离设计波长的面板仍只有PIC，未用基频输出伪造宽频BF。

GENESIS在前半程的局部BF增长较慢，随后持续增强；最后原生记录z={tail['z_m']:.6f} m，全局BF={tail['global_b']:.8f}，10切片局部RMS={tail['local_rms_10lambda']:.8f}。同一末段附近PIC两版的局部BF分别约0.15、0.35。它们不是完全相同位置的输出，不能据此给出精确出口百分比差。

## 复BF合并的口径

GENESIS的原始局部BF按一个光学切片统计。先把1386个活跃切片反序还原原输入顺序，再合并每10个原始切片（最后一组为6个）。必须先按电流加权复数，再取模，不能直接平均BF幅值：

$$B_j=\\frac{{\\sum_{{s\\in j}}I_s b_s e^{{i\\phi_s}}}}{{\\sum_{{s\\in j}}I_s}},
\\qquad b_{{\\mathrm{{local,rms}}}}=\\sqrt{{\\frac{{\\sum_j I_j|B_j|^2}}{{\\sum_j I_j}}}}.$$

切片间距固定为设计波长，电流权重与电荷权重成正比。输入切片电流在117个诊断位置逐位保持。设计波长处的全局BF为所有切片复数相加后取模。另保存5/20切片合并的敏感性及原始单切片RMS。

## 入口差异来自共同粒子筛选

|入口统计范围|全局BF|10切片局部RMS|
|---|---:|---:|
|GENESIS全部输入|{rows[0]['global_b']:.12g}|{rows[0]['local_rms_10lambda']:.12g}|
|PIC完整入口|{entry['global_b']:.12g}|{entry['local_rms_10lambda']:.12g}|
|前图PIC跨所有面共同ID子集|{subset['global_b']:.12g}|{subset['local_rms_10lambda']:.12g}|

完整PIC入口与GENESIS入口BF的绝对差约10⁻¹²，确认原始种子和复BF口径一致。前图为了沿程追踪同一批粒子，只保留两版、所有面均存在的99.948538%输入电荷。删去少量粒子会破坏初始微弱BF的相消，因此共同子集入口不完全重合；不是原始输入种子不同。本图保留前图PIC数值，没有人为对齐入口幅值。

## 原生诊断位置

以下均为实际GENESIS记录，不是插值出口值。程序未把最后14.89 m当成14.946 m波荡器出口或15.246 m远场面，也没有外推尾部。

|GENESIS原生z / m|全局BF|10切片局部RMS|
|---:|---:|---:|
'''
    for target in [0, 2.166, 4.722, 7.278, 9.834, 12.39, 14.946]:
        selected = [r for r in rows if r['z_m'] <= target + 1e-10][-1]
        text += (f"|{selected['z_m']:.6f}|{selected['global_b']:.8g}|"
                 f"{selected['local_rms_10lambda']:.8g}|\n")
    text += '''
当前图展示增长过程。GENESIS采用光学周期平均的FEL模型，PIC直接推进总E/B；单个种子及现有分辨率的差距仍不能单独证明哪一种推进已通过FEL精度验收。

原[BF与boost账本的时间定位图](../../reviewed/analysis/bunching_ledger_timing.png)保持；GENESIS按lab z输出，没有把它附会成PIC的逐步boost能量账本。

数据：genesis_bunching.csv、genesis_bunching.npz、RESULT.json；原PIC缓存来源及入口核查见VALIDATION.json。无新PIC或GENESIS传播。
'''
    (out / 'README.md').write_text(text)


def run(top, out):
    top, out = Path(top).resolve(), Path(out).resolve()
    out.mkdir(exist_ok=False, parents=True)
    reference = top / 'explicit_review_20261104/genesis-phase-r1'
    assert (reference / 'TRACK_COMPLETE').read_text().strip() == '0'
    assert read(reference / 'PHASE_CLOSURE_VERIFIED.json')['pass_nominal_closure']
    pic_path = top / 'bunching_growth_20261105/reviewed/analysis/RESULT.json'
    pic = read(pic_path)
    plan = read(top / 'compact-production-r1/plan.json')
    source = reference / 'reference.out.h5'
    with h5py.File(source, 'r') as f:
        z = f['Lattice/zplot'][:]
        current = f['Beam/current'][:]
        moment = f['Beam/bunching'][:] * np.exp(1j * f['Beam/bunchingphase'][:])
        wavelength = float(f['Global/lambdaref'][0])
        assert float(f['Global/sample'][0]) == 1
    assert np.isfinite(moment).all() and np.isfinite(current).all()
    assert np.all(current >= 0) and np.all(np.diff(z) > 0)
    assert np.array_equal(current, np.broadcast_to(current[0], current.shape))
    assert wavelength == plan['physics']['wavelength_m']
    active = np.flatnonzero(current[0] > 0)
    assert np.array_equal(active, np.arange(512, 512 + plan['active_slices']))
    current, moment = current[:, active[::-1]], moment[:, active[::-1]]
    total = current.sum(axis=1)
    global_b = abs((current * moment).sum(axis=1)) / total
    rows = [dict(z_m=float(zz), global_b=float(bb)) for zz, bb in zip(z, global_b)]
    arrays = dict(z_m=z, global_b=global_b)
    for width in [1, 5, 10, 20]:
        local, weights, rms = combine(current, moment, width)
        assert np.all(global_b <= rms + 1e-12) and np.all(rms <= 1 + 1e-12)
        arrays[f'local_complex_{width}lambda'] = local
        arrays[f'group_current_{width}lambda_A'] = weights[0]
        arrays[f'local_rms_{width}lambda'] = rms
        for row, value in zip(rows, rms):
            row[f'local_rms_{width}lambda'] = float(value)
    with np.load(top / 'explicit_review_20261104/radiation_phase_r1/genesis-beam.npz') as cached:
        assert np.max(abs(global_b - cached['global_b1'])) < 1e-14
        assert np.max(abs(arrays['local_rms_1lambda'] - cached['local_b1_rms'])) < 1e-14
    entry = entrance_check(top, wavelength, plan['NPART'], plan['input_particles'])
    initial_differences = {key: abs(entry[key] - rows[0][key])
                           for key in ['global_b', 'local_rms_10lambda']}
    assert max(initial_differences.values()) < 1e-10
    np.savez_compressed(out / 'genesis_bunching.npz', **arrays)
    with (out / 'genesis_bunching.csv').open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    save(out / 'RESULT.json', dict(genesis=rows, full_PIC_entrance=entry,
        wavelength_m=wavelength, last_GENESIS_z_m=float(z[-1]),
        PIC_source=str(pic_path), physical_acceptance=False,
        native_GENESIS_records_only=True, off_design_GENESIS_BF_available=False))
    render(out, pic, rows, wavelength)
    document(out, top, pic, rows, entry, wavelength)
    save(out / 'VALIDATION.json', dict(passed=True, raw_genesis_records=len(z),
        input_slice_current_constant=True, active_slices=len(active), reversed_to_original_input_order=True,
        full_input_BF_absolute_differences=initial_differences,
        old_reference_cache_checked=True, PIC_results_reused_unchanged=True,
        no_tail_extrapolation=True, no_off_design_BF_reconstruction=True,
        source_sha256={str(p): sha(p) for p in [source, pic_path, Path(__file__)]}))
    print(json.dumps(dict(output=str(out), last_GENESIS=rows[-1],
                         full_input_BF_absolute_differences=initial_differences), indent=2))


if __name__ == '__main__':
    run(*sys.argv[1:3])
