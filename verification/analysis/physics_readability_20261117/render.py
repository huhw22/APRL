"""Presentation-only revision of closed, delivered CSV/NPZ records.

No H5 reads, particle/field evolution, fitted curves, or residual-based fixes.
Every energy channel is assigned to one of two linear-scale panels. The
original signed sum is also shown on its own scale. Source tables stay intact.
"""
import csv
import hashlib
import json
import sys
import textwrap
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager, ticker
import numpy as np

LD = np.longdouble
COLORS = ['#1769a0', '#cf661f', '#308367', '#9356a0', '#b38917',
          '#618cbb', '#ba6572', '#525b66', '#54a9a8', '#9b793e']
STYLES = ['-', '--', '-.', ':']
METHODS = [('F3', 'linear-binomial3'), ('MM', 'quintic-compensated')]
GROUPS = {'uniform': ('冷均匀球', 'Cold uniform sphere'),
          'converging': ('会聚均匀球', 'Converging uniform sphere'),
          'gaussian': ('冷高斯球', 'Cold Gaussian sphere')}


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda: f.read(4194304), b''):
            h.update(b)
    return h.hexdigest()


def save(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2,
                               allow_nan=False) + '\n')


def read_table(path):
    return np.genfromtxt(path, names=True, delimiter=',', dtype=LD)


def configure(font):
    font_manager.fontManager.addfont(str(font))
    family = font_manager.FontProperties(fname=str(font)).get_name()
    plt.rcParams.update({'font.family': [family, 'DejaVu Sans'], 'font.size': 12,
        'axes.titlesize': 15, 'axes.labelsize': 12, 'legend.fontsize': 10,
        'axes.spines.top': False, 'axes.spines.right': False,
        'text.parse_math': False, 'axes.unicode_minus': False,
        'pdf.fonttype': 42, 'figure.facecolor': 'white',
        'savefig.facecolor': 'white', 'axes.formatter.useoffset': False})


def unit(values):
    peak = max((float(np.max(np.abs(v))) for v in values), default=0.)
    for name, amount in [('J', 1), ('mJ', 1e-3), ('µJ', 1e-6), ('nJ', 1e-9),
                         ('pJ', 1e-12), ('fJ', 1e-15), ('aJ', 1e-18),
                         ('zJ', 1e-21), ('yJ', 1e-24)]:
        if peak >= .7 * amount:
            return name, 1 / amount
    return '1e-27 J', 1e27


def array_info(x):
    a = np.asarray(x, dtype='<f8')
    return dict(shape=list(a.shape), sha256=hashlib.sha256(a.tobytes()).hexdigest(),
                minimum=float(a.min()), maximum=float(a.max()))


class Renderer:
    def __init__(self, root, font):
        self.root = Path(root).resolve()
        self.out = self.root / '07_readable_figures'
        self.out.mkdir(exist_ok=True)
        configure(font)
        self.sources, self.plots, self.books, self.theory = {}, [], [], []

    def table(self, rel):
        p = self.root / rel
        self.sources[str(rel)] = sha(p)
        return read_table(p)

    def finish(self, fig, rel, records, title):
        dest = self.out / rel
        dest.parent.mkdir(parents=True, exist_ok=True)
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        for ax in fig.axes:
            objs = [ax.title, ax.xaxis.label, ax.yaxis.label]
            if ax.get_legend() is not None:
                objs += list(ax.get_legend().texts)
            for obj in objs:
                assert '$' not in obj.get_text() and '_' not in obj.get_text()
                box = obj.get_window_extent(renderer)
                assert box.x0 >= -2 and box.y0 >= -2
                assert box.x1 <= fig.bbox.width + 2 and box.y1 <= fig.bbox.height + 2, obj.get_text()
        for ext in ['png', 'pdf']:
            fig.savefig(dest.with_suffix('.' + ext), dpi=160)
        self.plots.append(dict(path=str(rel.with_suffix('.png')), title=title,
                               data_axes=len(fig.axes), series=records))
        plt.close(fig)

    def curve(self, ax, x, y, label, index, records, key, is_sum=False):
        x, y = np.asarray(x, float), np.asarray(y, float)
        color = '#151b23' if is_sum else COLORS[index % len(COLORS)]
        line, = ax.plot(x, y, label=label, color=color,
                       ls='--' if is_sum else STYLES[index % 4],
                       lw=2.4 if is_sum else 2.,
                       marker=('D' if is_sum else ['o', 's', '^', 'v'][index % 4]) if len(x) < 24 else None,
                       ms=6, mfc='white', mew=1.3)
        assert np.array_equal(line.get_ydata(), y)
        records.append(dict(key=key, label=label, x=array_info(x), y=array_info(y)))

    @staticmethod
    def style(ax, xlabel, ylabel, title):
        ax.set(xlabel=xlabel, ylabel=ylabel, title=title)
        ax.axhline(0, color='#9ba2ac', lw=.85, zorder=0)
        ax.grid(alpha=.18)
        for axis in [ax.xaxis, ax.yaxis]:
            axis.set_major_formatter(ticker.FuncFormatter(lambda v, _: f'{v:.4g}'))
        ax.margins(x=.03, y=.18)

    def ledger(self, category, name, source, x, terms, total, title, frame,
               sum_key, note='', formula=''):
        # terms: (source/formula key, Chinese label, English label, signed joules).
        ys = [np.asarray(t[3], LD) for t in terms]
        total = np.asarray(total, LD)
        assert all(np.isfinite(y).all() and len(y) == len(x) for y in ys + [total])
        reconstructed = sum(ys)
        scale = max(float(max(np.max(abs(y)) for y in ys)), 1e-100)
        error = float(np.max(abs(reconstructed - total)))
        assert error <= max(scale * 3e-12, 1e-30), (name, error, scale)
        amplitudes = [float(np.max(abs(y))) for y in ys]
        major = [i for i, a in enumerate(amplitudes) if a >= .18 * max(amplitudes)]
        minor = [i for i in range(len(terms)) if i not in major]
        # Small sum always has an explicit scale even if all physical terms are large.
        book = dict(category=category, name=name, source=source, sum_key=sum_key,
                    formula=formula, note=note, reconstruction_max_error_J=error,
                    maximum_term_J=scale, final_sum_J=float(total[-1]),
                    sum_range_J=float(np.ptp(total)), major_keys=[terms[i][0] for i in major],
                    minor_keys=[terms[i][0] for i in minor],
                    terms=[dict(key=t[0], cn=t[1], en=t[2], maximum_absolute_J=amplitudes[i])
                           for i, t in enumerate(terms)])
        for lang in ['zh-CN', 'en']:
            zh = lang == 'zh-CN'
            rel = Path(category) / lang / name
            xlabel = ('计算系模拟时间（ps）' if zh else 'Computational-frame time (ps)') if frame == 'boost' else (
                      '实验室系记录面位置（m）' if zh else 'Laboratory detector position (m)')
            # Fixed layout leaves a dedicated bottom band for each legend.
            fig, axes = plt.subplots(1, 2, figsize=(15.2, 7.3))
            fig.subplots_adjust(left=.075, right=.985, bottom=.34, top=.80, wspace=.29)
            fig.suptitle(title[0 if zh else 1], y=.96, fontsize=17)
            fig.text(.5, .89, '每条分量只分配到一侧；两图合起来构成完整收支' if zh else
                     'Every component appears in one panel; both panels form the complete budget', ha='center', fontsize=11)
            records = []
            for ax, indices, small in [(axes[0], major, False), (axes[1], minor, True)]:
                name_unit, factor = unit([ys[i] for i in indices] + ([total] if small else []))
                for i in indices:
                    self.curve(ax, x, ys[i] * factor, terms[i][1 if zh else 2], i,
                               records, terms[i][0])
                if small:
                    self.curve(ax, x, total * factor, '完整加和' if zh else 'Complete sum', 0,
                               records, sum_key, True)
                self.style(ax, xlabel, ('带符号能量' if zh else 'Signed energy') + f' ({name_unit})',
                           ('小分量与加和' if zh else 'Small components and sum') if small else
                           ('主要分量' if zh else 'Major components'))
                handles, labels = ax.get_legend_handles_labels()
                if not zh:
                    labels = [textwrap.fill(s, width=28) for s in labels]
                ax.legend(handles, labels, loc='upper center', bbox_to_anchor=(.5, -.24), ncol=2,
                          frameon=False, columnspacing=1.1, handlelength=2.8, fontsize=9.5)
            fig.text(.5, .025, '线性坐标；两侧单位以轴标签为准，未逐线归一化；加和另附独立图' if zh else
                     'Linear scales; read each axis unit. No per-curve normalization. A separate sum plot is included.',
                     ha='center', fontsize=10)
            self.finish(fig, rel.with_name(name + '_pair'), records, title[0 if zh else 1])
            fig, ax = plt.subplots(figsize=(10, 5.7))
            fig.subplots_adjust(left=.13, right=.97, bottom=.20, top=.81)
            u, fac = unit([total]); records = []
            self.curve(ax, x, total * fac, '完整加和' if zh else 'Complete sum', 0, records, sum_key, True)
            self.style(ax, xlabel, ('带符号能量' if zh else 'Signed energy') + f' ({u})',
                       title[0 if zh else 1] + ('：加和' if zh else ': sum'))
            ax.legend(frameon=False)
            fig.text(.5, .035, '连续物理收支的预期为零；离散完整账本包含单列的数值通道，均保留实测非零值' if zh else
                     'Zero is the continuous-balance target. Discrete budgets include the listed numerical channels.',
                     ha='center', fontsize=9.5)
            self.finish(fig, rel.with_name(name + '_sum'), records, title[0 if zh else 1])
        self.books.append(book)
        print('LEDGER', category, name, flush=True)


def old_boost(r, category, name, rel, title):
    a = r.table(rel)
    keys = [('particle_delta_K_J', '粒子动能变化', 'Particle kinetic-energy change'),
            ('field_delta_U_native_J', '内场能量变化', 'Native field-energy change'),
            ('front_outward_minus_injection_J', '前端净传出（减注入）', 'Front outward minus injection'),
            ('rear_outward_plus_crop_J', '后端净传出（加裁出）', 'Rear outward plus cropped energy'),
            ('side_outward_J', '四侧净能流', 'Net transverse flux'),
            ('negative_constraint_and_external_work_J', '外源及约束功的负值', 'Negative external and constraint work'),
            ('explicit_discrete_correction_J', '独立离散数值项', 'Independent discrete channels')]
    r.ledger(category, name, [str(rel)], a['time_box_s'] * 1e12,
             [(*t, a[t[0]]) for t in keys], a['discrete_accounting_sum_J'], title,
             'boost', 'discrete_accounting_sum_J',
             formula='particle + field + front + rear + sides - source work + numerical channels')


def lab(r, category, name, rel, title, modern=False):
    a = r.table(rel)
    keys = [('particle_net_transport_J' if modern else 'particle_net_K_J', '粒子净能量变化', 'Net particle energy'),
            ('field_ports_net_J' if modern else 'field_end_ports_net_J', '场出口减入口', 'Field exit minus entrance'),
            ('field_transverse_net_J' if modern else 'transverse_total_J', '场四侧净能流', 'Transverse field flux'),
            ('field_time_caps_net_J' if modern else 'time_caps_net_J', '晚减早时间端面', 'Late minus early time cap')]
    if modern:
        keys.append(('particle_time_caps_J', '粒子时间端面存量', 'Particle time-cap storage'))
    terms = [(*t, a[t[0]]) for t in keys]
    total_key = 'physical_lab_sum_J' if modern else 'complete_sum_J'
    r.ledger(category, name, [str(rel)], a['z_lab_m'], terms, a[total_key], title,
             'lab', total_key, formula='net particle + particle caps + field ports + field sides + field time caps')
    # The raw entrance/exit and early/late terms remain explicitly visible.
    detailed = [('particle_net_transport_J' if modern else 'particle_net_K_J', '粒子净能量变化', 'Net particle energy'),
                ('field_entrance_signed_J', '场入口（已含负号）', 'Field entrance (signed)'),
                ('field_exit_signed_J', '场出口（已含正号）', 'Field exit (signed)'),
                ('field_transverse_net_J' if modern else 'transverse_total_J', '场四侧净能流', 'Transverse field flux'),
                ('time_cap_early_signed_J', '早时间端面（已含负号）', 'Early time cap (signed)'),
                ('time_cap_late_signed_J', '晚时间端面（已含正号）', 'Late time cap (signed)')]
    if modern:
        detailed.append(('particle_time_caps_J', '粒子时间端面存量', 'Particle time-cap storage'))
    r.ledger(category, name + '_surfaces', [str(rel)], a['z_lab_m'],
             [(*t, a[t[0]]) for t in detailed], a[total_key],
             (title[0] + '：入口出口分别显示', title[1] + ': signed surfaces'),
             'lab', total_key, formula='net particle + particle caps + signed entrance + signed exit + sides + early + late')


def discrete(r, category, name, rel, x, title):
    a = r.table(rel)
    # All entries below enter the sum with their displayed signs.
    spec = [('particle_change_J', 1, '粒子动能变化', 'Particle kinetic-energy change'),
            ('native_field_change_J', 1, '内场能量变化', 'Native field-energy change'),
            ('outward_curl_exchange_J', 1, '原生旋度边界交换', 'Native curl boundary exchange'),
            ('inlet_energy_added_J', -1, '入口覆盖注入的负值', 'Negative inlet overwrite'),
            ('remap_energy_added_J', -1, '移窗加入能量的负值', 'Negative remap energy added'),
            ('prescribed_and_constraint_source_work_J', -1, '外源及约束功的负值', 'Negative source and constraint work'),
            ('time_pairing_channel_J', -1, '时间配对项的负值', 'Negative time-pairing channel'),
            ('Cowan_channel_J', -1, 'Cowan 通道的负值', 'Negative Cowan channel'),
            ('physical_CPML_channel_J', -1, 'CPML 直接场变化的负值', 'Negative direct CPML field change')]
    terms = [(('+' if sign == 1 else '-') + k, cn, en, sign * a[k]) for k, sign, cn, en in spec]
    r.ledger(category, name, [str(rel)], x, terms, a['total_residual_J'], title,
             'boost', 'total_residual_J', formula='K + U + curl - inlet - remap - source - time - Cowan - CPML')


def spheres(r, theory_only=False):
    for kind, names in GROUPS.items():
        base = Path('01_reviewed_physics/02_space_charge') / kind
        data = []
        for method, folder in METHODS:
            title = (names[0] + ' ' + method, names[1] + ' ' + method)
            if not theory_only:
                old_boost(r, '02_space_charge/' + kind, method + '_boost', base / folder / 'boost_energy_budget.csv',
                          (title[0] + '：计算系完整账本', title[1] + ': computational budget'))
                lab(r, '02_space_charge/' + kind, method + '_lab', base / folder / 'lab_energy_budget.csv',
                    (title[0] + '：实验室系账本', title[1] + ': laboratory budget'))
            data.append(r.table(base / folder / 'theory_comparison.csv'))
        a, b = data
        specs = [('size', 'PIC_sigma_x_m', 'reference_sigma_x_m', '横向 RMS 尺寸', 'Transverse RMS size', 1e6, 'µm', 1e9, 'nm'),
                 ('mean_rest_gain', 'PIC_mean_rest_gain_eV', 'reference_mean_rest_gain_eV', '静止系平均能量变化', 'Mean rest-frame energy change', 1, 'eV', 1, 'eV'),
                 ('std_rest_gain', 'PIC_std_rest_gain_eV', 'reference_std_rest_gain_eV', '静止系能量变化标准差', 'Rest-frame energy-change standard deviation', 1, 'eV', 1, 'eV'),
                 ('mean_lab_gain', 'PIC_mean_lab_gain_eV', 'reference_mean_lab_gain_eV', '实验室系平均能量变化', 'Mean laboratory energy change', .001, 'keV', 1, 'eV')]
        refdiff = dict(case=kind, model='one relativistic spherical shell model, sampled at actual per-ID crossing times',
                       maximum_reference_size_difference_nm=float(max(abs(a['reference_sigma_x_m']-b['reference_sigma_x_m']))*1e9),
                       maximum_reference_mean_gain_difference_eV=float(max(abs(a['reference_mean_rest_gain_eV']-b['reference_mean_rest_gain_eV']))))
        r.theory.append(refdiff)
        for name, pk, tk, cn, en, factor, u, df, du in specs:
            for lang in ['zh-CN', 'en']:
                zh = lang == 'zh-CN'
                fig, axes = plt.subplots(2, 1, figsize=(10.8, 9.3), sharex=True,
                                         gridspec_kw={'height_ratios': [1.3, 1]})
                fig.subplots_adjust(left=.14, right=.97, bottom=.19, top=.89, hspace=.34)
                title = (names[0] + '：' + cn) if zh else (names[1] + ': ' + en)
                fig.suptitle(title, y=.975, fontsize=17)
                records = []
                # One physical model. F3 event samples form the visible black line;
                # the thin envelope includes both actual-event sample sets.
                x = a['lab_z_m']
                axes[0].fill_between(np.asarray(x, float), np.asarray(np.minimum(a[tk], b[tk])*factor, float),
                                     np.asarray(np.maximum(a[tk], b[tk])*factor, float), color='#b5bbc2', alpha=.5)
                self_label = '同一球壳理论' if zh else 'One spherical-shell theory'
                axes[0].plot(np.asarray(x, float), np.asarray(a[tk]*factor, float), color='#252a31', lw=2.2,
                             label=self_label, zorder=1)
                records.append(dict(key=tk, label=self_label, x=array_info(x), y=array_info(a[tk]*factor)))
                for i, (d, method) in enumerate(zip(data, ['F3', 'MM'])):
                    label = ('三点核' if i == 0 else '动量匹配') if zh else ('EC + three-point' if i == 0 else 'EC + compensated quintic')
                    # Hollow, different markers preserve visibility when all lines overlap.
                    for ax, values, key in [(axes[0], d[pk]*factor, pk),
                                            (axes[1], (d[pk]-d[tk])*df, pk+' minus '+tk)]:
                        xx, yy = np.asarray(d['lab_z_m'], float), np.asarray(values, float)
                        ax.plot(xx, yy, color=COLORS[i], marker='o' if i == 0 else 's',
                                ms=9 if i == 0 else 6, mfc='none', mew=1.6,
                                ls='--' if i == 0 else '-.', lw=1.7, label=label, zorder=3+i)
                        records.append(dict(key=key, label=label, x=array_info(xx), y=array_info(yy)))
                self_x = '实验室系记录面位置（m）' if zh else 'Laboratory detector position (m)'
                r.style(axes[0], '', cn + f' ({u})' if zh else en + f' ({u})',
                        '重合表示数值接近' if zh else 'Overlapping curves indicate close values')
                # Size curves should not include a forced zero if it obscures differences.
                if name == 'size':
                    vals = np.r_[a[pk], b[pk], a[tk], b[tk]]*factor
                    lo, hi = float(vals.min()), float(vals.max());pad=max((hi-lo)*.15, .1)
                    axes[0].set_ylim(lo-pad, hi+pad)
                r.style(axes[1], self_x, ('仿真减理论' if zh else 'Simulation minus theory') + f' ({du})',
                        '分别减去同一模型在各自过面事件的取样值' if zh else 'Each method minus the same model at its own crossing events')
                axes[0].legend(frameon=False, loc='best')
                axes[1].legend(frameon=False, loc='best')
                fig.text(.5, .055, '黑线为三点核事件取样，灰带覆盖两组事件取样差；没有为推进方法另设理论或拟合' if zh else
                         'Black: F3-event samples; grey: the two event-sampling values. No method-specific model or fit.',
                         ha='center', fontsize=10)
                r.finish(fig, Path('02_space_charge')/kind/lang/('theory_'+name), records, title)


def main(root, font):
    r = Renderer(root, font)
    base = Path('01_reviewed_physics/01_vacuum_and_short_drift')
    for method, folder in METHODS:
        old_boost(r, '01_short_drift', method+'_boost', base/f'boost_{folder}.csv',
                  (method+' 短漂移：计算系完整账本', method+' short drift: computational budget'))
        lab(r, '01_short_drift', method+'_lab', base/f'lab_complete_{folder}.csv',
            (method+' 短漂移：实验室系完整账本', method+' short drift: laboratory budget'))
    spheres(r)
    base = Path('02_FEL_six/01_energy_curves/data')
    for method, _ in METHODS:
        t = r.table(base/f'{method}_time_curves.csv')['time_box_s']*1e12
        discrete(r, '03_FEL_six', method+'_boost_discrete', base/f'{method}_boost_components.csv', t,
                 (method+' 六段 FEL：计算系完整离散账本', method+' six-segment FEL: discrete computational budget'))
        lab(r, '03_FEL_six', method+'_lab', base/f'{method}_lab_position_curves.csv',
            (method+' 六段 FEL：实验室系完整账本', method+' six-segment FEL: laboratory budget'), True)
    base = Path('06_synchrotron_current')
    result = json.loads((r.root/base/'RESULT.json').read_text())
    for case in result['cases']:
        tag = case['tag'];p = base/'cases'/tag
        a = r.table(p/'boost/physical_components.csv')
        kind = '弯铁' if case['kind']=='bend' else '直线'
        cn = f"{case['method']} {kind} 样本{case['seed']-20261110}"
        en = f"{case['method']} {case['kind']} sample {case['seed']-20261110}"
        spec = [('particle_change_J','粒子动能变化','Particle kinetic-energy change'),
                ('field_change_J','内场能量变化','Field-energy change'),
                ('front_outward_minus_injection_J','前端净传出（减注入）','Front outward minus injection'),
                ('back_outward_plus_cropped_J','后端净传出（加裁出）','Back outward plus cropped energy'),
                ('side_outward_J','四侧净能流','Transverse outward flux'),
                ('moving_interface_J','移窗接口项','Moving-interface term'),
                ('external_and_constraint_work_signed_J','外源及约束功的负值','Negative external and constraint work')]
        r.ledger('04_synchrotron', tag+'_boost_physical', [str(p/'boost/physical_components.csv')],
                 a['time_box_s']*1e12, [(*s,a[s[0]]) for s in spec], a['physical_component_sum_J'],
                 (cn+'：计算系物理分项',en+': computational physical components'), 'boost', 'physical_component_sum_J')
        terms = [(*s,a[s[0]]) for s in spec]
        terms[0] = ('particle_change_J + external_prescribed_work_signed_J', '粒子变化减已记外磁功',
                    'Particle change minus recorded magnetic work', a['particle_change_J']+a['external_prescribed_work_signed_J'])
        terms[-1] = ('external_and_constraint_work_signed_J - external_prescribed_work_signed_J', '其余边界源及约束功的负值',
                     'Negative remaining source and constraint work', a['external_and_constraint_work_signed_J']-a['external_prescribed_work_signed_J'])
        r.ledger('04_synchrotron', tag+'_boost_internal', [str(p/'boost/physical_components.csv')],
                 a['time_box_s']*1e12, terms, a['physical_component_sum_J'],
                 (cn+'：扣除已记外磁功后的交换',en+': exchange after recorded magnetic work'),
                 'boost','physical_component_sum_J',note='Known recorded external magnetic work regrouped algebraically; unchanged sum, not self-field subtraction.')
        discrete(r, '04_synchrotron', tag+'_boost_discrete', p/'boost/complete_components.csv', a['time_box_s']*1e12,
                 (cn+'：计算系完整离散账本',en+': discrete computational budget'))
        lab(r, '04_synchrotron', tag+'_lab', p/'lab/complete_components.csv',
            (cn+'：实验室系完整账本',en+': laboratory budget'), True)
    # Make traceable, compact plot records rather than duplicating all source series.
    save(r.out/'PLOTS.json', dict(no_physics_rerun=True, no_individual_curve_normalization=True,
                                 all_original_points_used=True, figures=r.plots))
    save(r.out/'LEDGERS.json', r.books)
    save(r.out/'THEORY_SAMPLING.json', r.theory)
    save(r.out/'SOURCE_TABLE_SHA256.json', r.sources)
    for rel, expected in r.sources.items():
        assert sha(r.root/rel)==expected,rel
    save(r.out/'RENDER_COMPLETE.json', dict(complete=True, ledgers=len(r.books),
         figures=len(r.plots), source_tables=len(r.sources), source_tables_unchanged=True))


if __name__ == '__main__':
    main(*sys.argv[1:])
