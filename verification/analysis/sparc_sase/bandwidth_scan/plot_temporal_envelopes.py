#!/usr/bin/env python3
"""Plot cached SPARC power envelopes with explicit time and current-profile axes."""
import argparse
import csv
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def run(root):
    scan = root/'analysis-energy-binomial3/bandwidth-scan-r1'
    result = json.loads((scan/'RESULT.json').read_text())
    output = root/'analysis-energy-binomial3/temporal-envelope-review'
    output.mkdir(exist_ok=False)
    selected = [5, 6, 8, 11]
    colors = ['#24648c', '#8e4b75', '#cf8a24', '#75884b']
    checks = {row['module']: row for row in result['plane_checks']}
    plan = json.loads((root/'plan.json').read_text())
    peak_current = float(plan['physics']['peak_current_A'])
    profiles = []
    maxima = []
    for module in selected:
        with np.load(scan/f'p{module:02d}_envelopes.npz') as data:
            index = np.flatnonzero(np.isclose(data['widths'], .05)).item()
            t = data['time_offset_from_midpoint_s']*1e12
            u = data['u']
            power = data['power_E_W'][index]/1e6
            current = peak_current*np.maximum(1-u*u, 0)
            profiles.append((t, u, power, current))
            peak = int(np.argmax(power))
            maxima.append(dict(module=module, z_m=checks[module]['z_m'],
                full_record_peak_MW=float(power[peak]),
                peak_time_relative_electron_midpoint_ps=float(t[peak]),
                nominal_half_duration_ps=checks[module]['interval']['half_duration_s']*1e12))
    plt.rcParams.update({'font.size': 12, 'axes.spines.top': False,
                         'axes.spines.right': False, 'lines.linewidth': 1.8})
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharex=True, sharey=True,
                             layout='constrained')
    edge = math.sqrt(.2)
    for ax, module, color, (t, u, power, current) in zip(axes.flat, selected, colors, profiles):
        half = checks[module]['interval']['half_duration_s']*1e12
        use = (t >= -2.0) & (t <= 1.6)
        ax.plot(t[use], power[use], color=color, label='±5% optical power envelope')
        ax.axvspan(.1*edge*half, .9*edge*half, color='#aaaaaa', alpha=.18,
                    label='previous central averaging interval')
        ax.axvline(0, color='#555555', linestyle=':', linewidth=1)
        right = ax.twinx()
        right.plot(t[use], current[use], color='#555555', linestyle='--',
                    linewidth=1.3, label='nominal input current profile')
        right.set_ylim(0, 60)
        right.spines['top'].set_visible(False)
        right.spines['right'].set_visible(True)
        right.set_ylabel('nominal current [A]')
        ax.set_ylim(0, 12)
        ax.set_xlim(-2.0, 1.6)
        ax.set_title(f'Module {module}: z = {checks[module]["z_m"]:.3f} m')
        ax.set_xlabel('time relative to local electron-bunch midpoint [ps]')
        ax.set_ylabel('power envelope [MW]')
        ax.grid(alpha=.20)
        if module == 5:
            lines, labels = ax.get_legend_handles_labels()
            lines2, labels2 = right.get_legend_handles_labels()
            ax.legend(lines+lines2, labels+labels2, fontsize=8, loc='upper right')
    fig.savefig(output/'power_vs_time.png', dpi=190)
    fig.savefig(output/'power_vs_time.pdf')
    plt.close(fig)
    with (output/'full_record_peaks.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(maxima[0]))
        writer.writeheader()
        writer.writerows(maxima)
    report = dict(complete=True, source='bandwidth-scan-r1 cached envelopes',
        relative_half_width=.05, time_origin='local electron source-slice midpoint',
        raw_lab_field=True, additional_hann_window=False, full_record_maxima=maxima,
        current_curve='nominal parabolic input profile aligned to the measured local midpoint; not an independent measured output current',
        limitations=['temporal flatness is not required by spatial gain saturation',
            'a later-half time cut does not establish longitudinal steady state',
            'the current trace alone does not quantify the mechanism of the optical decline'])
    (output/'RESULT.json').write_text(json.dumps(report, indent=2)+'\n')
    (output/'README.md').write_text('''# 功率时间包络复查

![完整时间包络](power_vs_time.png)

横轴是各探测面相对于电子束中心切片过面时刻的时间；负时间也属于
有效光学脉冲。实线是±5%带限光学功率包络，灰色区是此前中心区均值的
时间窗口。虚线是按当地电子束中心对齐的名义输入抛物线电流，并非独立
测得的出口电流。包络使用原始总场，没有再乘Hann窗。

本工况将完整抛物线电流压缩到2.776 ps，保持峰值电流；并未从长束团
截取一个平顶中心区。名义11段总滑移时间约1.388 ps，短束团中的电流变化
和有限脉冲作用仍不可忽略。单个SASE脉冲的时间包络不要求平坦，沿z方向
的饱和与时间包络的平台需要分别判断。输入电流、滑移与非线性演化各自
对当前形状的定量贡献尚未分解；没有平台本身不能验证或否定PIC精度。

全记录峰值与此前后半束团峰值不同：例如第11段全记录峰值约11.04 MW，
位于电子束中心前0.531 ps，属于此前后半段筛选排除的时间区域。
''')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('campaign', type=Path)
    run(parser.parse_args().campaign.resolve())
