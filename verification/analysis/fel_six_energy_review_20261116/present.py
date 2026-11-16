"""Standalone bilingual figures and report for the closed six-segment run."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager, ticker


def read(path):
    return json.loads(Path(path).read_text())


def table(path):
    return np.genfromtxt(path, delimiter=',', names=True, dtype=np.longdouble)


def main(review, ledger, out, font, only=None):
    review, ledger, out, font = map(Path, [review, ledger, out, font])
    result = read(review / 'data/RESULT.json')
    energy = read(ledger / 'RESULT.json')
    out.mkdir(exist_ok=only is not None)
    (out / 'data').mkdir(exist_ok=only is not None)
    (out / 'ledger').mkdir(exist_ok=only is not None)
    for name in ['RESULT.json', 'growth.csv', 'initial_bunching.csv', 'comparison_by_segment.csv']:
        shutil.copy2(review / 'data' / name, out / 'data' / name)
    for source, name in [(ledger / 'RESULT.json', 'RESULT.json'),
                         (ledger / 'boost/energy_budget_absolute.csv', 'boost.csv'),
                         (ledger / 'boost/normalization.json', 'boost_normalization.json'),
                         (ledger / 'whole/energy-momentum/whole_lab_slice.csv', 'lab_on_boost_slices.csv'),
                         (ledger / 'lab/energy-momentum/lab_energy_absolute.csv', 'fixed_lab.csv'),
                         (ledger / 'lab/GEOMETRY_CHECK.json', 'GEOMETRY_CHECK.json')]:
        shutil.copy2(source, out / 'ledger' / name)
    font_manager.fontManager.addfont(str(font))
    family = font_manager.FontProperties(fname=str(font)).get_name()
    plt.rcParams.update({'font.family': [family, 'DejaVu Sans'], 'font.size': 12,
        'axes.titlesize': 17, 'axes.labelsize': 13, 'legend.fontsize': 10,
        'axes.spines.top': False, 'axes.spines.right': False, 'axes.unicode_minus': False,
        'text.parse_math': False, 'axes.formatter.use_mathtext': False,
        'figure.facecolor': 'white', 'savefig.facecolor': 'white', 'pdf.fonttype': 42})
    colors = ['#2679ae', '#d37b28', '#57916f', '#a45d88', '#777777', '#8a6a3c', '#359c9d', '#7066a3']
    styles = ['-', '--', '-.', ':']
    names = ['Puffin', 'APRL EC+F3', 'APRL EC+MM']
    zh_names = ['Puffin', '能量守恒 + 三点核', '能量守恒 + 动量匹配']
    plots = []

    def plot(language, name, x, series, title, xlabel, ylabel, log=False, markers=False):
        if only is not None and name != only:
            return
        target = out / 'figures' / language
        target.mkdir(parents=True, exist_ok=True)
        fig, ax = plt.subplots(figsize=(10, 6), layout='constrained')
        arrays = []
        for index, item in enumerate(series):
            label, values = item[:2]
            xx, yy = np.asarray(item[2] if len(item) == 3 else x, float), np.asarray(values, float)
            keep = np.isfinite(xx) & np.isfinite(yy)
            line, = ax.plot(xx[keep], yy[keep], label=label, color=colors[index % len(colors)],
                linestyle=styles[index % len(styles)], linewidth=2 if index < 3 else 1.4,
                marker=['o', 's', '^'][index % 3] if markers else None, markersize=5)
            assert np.array_equal(line.get_ydata(), yy[keep])
            arrays.append(dict(label=label, x=xx[keep].tolist(), y=yy[keep].tolist()))
        ax.set(title=title, xlabel=xlabel, ylabel=ylabel)
        if log:
            ax.set_yscale('log')
        else:
            ax.axhline(0, color='#888888', linewidth=.7)
        ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, pos: f'{y:g}'))
        ax.grid(alpha=.18, which='major')
        ax.legend(frameon=False, ncol=2 if len(series) > 4 else 1)
        fig.canvas.draw()
        hidden = set()
        for axis, limits in [(ax.xaxis, ax.get_xlim()), (ax.yaxis, ax.get_ylim())]:
            low, high = sorted(limits)
            for tick in axis.get_major_ticks() + axis.get_minor_ticks():
                if not low <= tick.get_loc() <= high:
                    hidden.update([id(tick.label1), id(tick.label2)])
        for text in fig.findobj(matplotlib.text.Text):
            if not text.get_visible() or not text.get_text() or id(text) in hidden:
                continue
            assert '$' not in text.get_text() and '_' not in text.get_text()
            b = text.get_window_extent(fig.canvas.get_renderer())
            assert b.x0 >= -2 and b.y0 >= -2 and b.x1 <= fig.bbox.width + 2 and b.y1 <= fig.bbox.height + 2, text.get_text()
        for suffix in ['png', 'pdf']:
            fig.savefig(target / f'{name}.{suffix}', dpi=175)
        plt.close(fig)
        plots.append(dict(language=language, name=name, axes=1, arrays=arrays))

    boost = table(out / 'ledger/boost.csv')
    whole = table(out / 'ledger/lab_on_boost_slices.csv')
    fixed = table(out / 'ledger/fixed_lab.csv')
    dt = read(ledger.parent / 'six-aprl/PLAN.json')['dt_box_s']
    for language in ['zh-CN', 'en']:
        zh = language == 'zh-CN'
        labx = '实验室系位置（m）' if zh else 'Laboratory position (m)'
        boxx = '计算系模拟时间（ps）' if zh else 'Computational-frame time (ps)'
        micro = '带符号能量（µJ）' if zh else 'Signed energy (µJ)'
        for name, key, scale, title, ylabel, use_initial, log in [
            ('power', 'Hann_pm5pct_power_Eproxy_W', 1e-6,
             '后半束团高流区的功率增长' if zh else 'Power in the rear high-current window',
             '中央频带功率（MW）' if zh else 'Windowed resonant power (MW)', False, True),
            ('bunching', 'central_10lambda_rms_b1', 1,
             '后半束团高流区的局部聚束' if zh else 'Local bunching in the rear high-current window',
             '局部 RMS 聚束因子' if zh else 'Local RMS bunching factor', True, False),
            ('particle_energy_loss', 'mean_loss_eV', .001,
             '整束平均粒子能损' if zh else 'Whole-beam mean particle energy loss',
             '平均能损（keV / 电子）' if zh else 'Mean energy loss (keV / electron)', True, False)]:
            series = []
            for i, method in enumerate(names):
                rows = [r for r in result['rows'] if r['method'] == method]
                if use_initial:
                    rows = [r for r in result['initial'] if r['method'] == method] + rows
                x = [r['actual_z_m'] for r in rows]
                series.append(((zh_names if zh else names)[i], [r[key] * scale for r in rows], x))
            plot(language, name, x, series, title, labx, ylabel, log, True)

        def budget(name, x, source, keys, title, xlabel=boxx, unit=micro, scale=1e6):
            series = [(a if zh else b, source[k] * scale * sign) for k, a, b, sign in keys]
            plot(language, name, x, series, title, xlabel, unit, markers=(xlabel == labx))

        t = boost['step'] * dt * 1e12
        budget('boost_balance', t, boost, [
            ('particle_change_J', '粒子能量变化', 'Particle energy change', 1),
            ('field_boundary_source_prediction_J', '场与边界独立预测', 'Independent field and boundary prediction', 1),
            ('total_residual_J', '全部分项加和', 'Sum of all signed terms', 1)],
            '动量匹配：计算系能量收支' if zh else 'Momentum matching: computational-frame energy balance')
        budget('boost_components', t, boost, [
            ('particle_change_J', '粒子能量变化', 'Particle energy change', 1),
            ('native_field_change_J', '原生场能变化', 'Native field energy change', 1),
            ('outward_curl_exchange_J', '离散 curl 边界交换', 'Discrete curl boundary exchange', 1),
            ('prescribed_and_constraint_source_work_J', '规定源与约束功的负值', 'Negative prescribed and constraint work', -1),
            ('inlet_energy_added_J', '入口覆盖能量的负值', 'Negative inlet overwrite energy', -1),
            ('remap_energy_added_J', '移窗能量的负值', 'Negative moving-window remap energy', -1)],
            '动量匹配：计算系场、粒子与边界分项' if zh else 'Momentum matching: field, particle and boundary terms')
        budget('boost_numerical', t, boost, [
            ('time_pairing_channel_J', '时间与电流配对项的负值', 'Negative time and current pairing', -1),
            ('Cowan_channel_J', 'Cowan 项的负值', 'Negative Cowan operator term', -1),
            ('physical_CPML_channel_J', '吸收层项的负值', 'Negative CPML term', -1)],
            '动量匹配：计算系离散数值分项' if zh else 'Momentum matching: computational discrete terms')
        budget('boost_residual', t, boost, [
            ('total_residual_J', '全部分项加和', 'Sum of all signed terms', 1),
            ('particle_grid_exchange_difference_J', '粒子与电流交换差', 'Particle and current exchange gap', 1),
            ('field_operator_residual_J', '纯场算子恒等式余差', 'Field operator identity residual', 1)],
            '动量匹配：计算系账本剩余量' if zh else 'Momentum matching: computational ledger residual',
            unit='带符号能量（aJ）' if zh else 'Signed energy (aJ)', scale=1e18)
        budget('lab_slice_balance', t, whole, [
            ('particle_kinetic_change_lab_J', '实验室系粒子能变', 'Laboratory particle energy change', 1),
            ('field_boundary_source_prediction_lab_J', '场与能量动量通道预测', 'Field energy and momentum prediction', 1),
            ('full_window_Lorentz_accounting_difference_J', '全部分项加和', 'Sum of all signed terms', 1)],
            '动量匹配：计算系等时面上的实验室系能量' if zh else 'Laboratory energy on computational-time slices')
        z = fixed['z_lab_m']
        budget('lab_physical_balance', z, fixed, [
            ('particle_net_transport_J', '粒子净输运（含侧向）', 'Net particle transport including sides', 1),
            ('field_ports_net_J', '场出口减入口', 'Field exit minus entrance', 1),
            ('field_transverse_net_J', '场四侧净能流', 'Net field energy through four sides', 1),
            ('field_time_caps_net_J', '场晚端面减早端面', 'Late minus early field cap', 1),
            ('physical_lab_sum_J', '完整物理收支加和', 'Complete physical balance sum', 1)],
            '动量匹配：固定实验室观测管能量收支' if zh else 'Momentum matching: fixed laboratory worldtube balance', xlabel=labx)
        budget('lab_field_ports', z, fixed, [
            ('field_entrance_signed_J', '入口场能流（带负号）', 'Entrance field flux with outward sign', 1),
            ('field_exit_signed_J', '出口场能流（带正号）', 'Exit field flux with outward sign', 1),
            ('field_ports_net_J', '出口减入口', 'Exit minus entrance', 1)],
            '动量匹配：场的轴向出入能流' if zh else 'Momentum matching: longitudinal field ports', xlabel=labx)
        budget('lab_field_sides', z, fixed, [
            ('transverse_xminus_signed_J', '横向 x 负侧', 'Negative x side', 1),
            ('transverse_xplus_signed_J', '横向 x 正侧', 'Positive x side', 1),
            ('transverse_yminus_signed_J', '横向 y 负侧', 'Negative y side', 1),
            ('transverse_yplus_signed_J', '横向 y 正侧', 'Positive y side', 1),
            ('field_transverse_net_J', '四侧加和', 'Sum of four sides', 1)],
            '动量匹配：四个侧面的场能流' if zh else 'Momentum matching: four transverse field sides', xlabel=labx)
        budget('lab_field_caps', z, fixed, [
            ('time_cap_early_signed_J', '早时间端面（带负号）', 'Early time cap with negative sign', 1),
            ('time_cap_late_signed_J', '晚时间端面（带正号）', 'Late time cap with positive sign', 1),
            ('field_time_caps_net_J', '晚端面减早端面', 'Late minus early time cap', 1)],
            '动量匹配：场的时间端面存量' if zh else 'Momentum matching: field time-cap storage', xlabel=labx)
        budget('lab_particle_transport', z, fixed, [
            ('particle_ports_net_J', '粒子出口减入口', 'Particle exit minus entrance', 1),
            ('particle_lateral_signed_J', '粒子侧向净输运', 'Net lateral particle transport', 1),
            ('particle_net_transport_J', '粒子净输运总和', 'Total net particle transport', 1)],
            '动量匹配：粒子的端口与侧向输运' if zh else 'Momentum matching: particle ports and side transport', xlabel=labx)
    if only is not None:
        old = read(out / 'PLOTS.json')['figures']
        plots = [p for p in old if p['name'] != only] + plots
    (out / 'PLOTS.json').write_text(json.dumps(dict(single_axes_only=True, figures=plots), ensure_ascii=False, indent=2) + '\n')
    (out / 'PROVENANCE.json').write_text(json.dumps(dict(review=str(review.resolve()), ledger=str(ledger.resolve()),
        original_statistical_definitions_preserved=True, no_combined_panel_figures=True,
        no_new_physics=True), indent=2) + '\n')
    print('PRESENTED', out, 'figures', len(plots), flush=True)


if __name__ == '__main__':
    main(*os.sys.argv[1:])
