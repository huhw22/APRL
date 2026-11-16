"""Package the reviewed six-segment comparison; never open scientific H5 files.

Prepare writes Chinese presentation figures from the reviewed small JSON/CSV
tables and copies English originals. Finalize follows visual inspection, checks
all local links and source identities again, and builds a checksummed ZIP.
"""
import argparse
import csv
import hashlib
import json
import math
import re
import shutil
import subprocess
import zipfile
from datetime import datetime
from pathlib import Path
from zoneinfo import ZoneInfo


METHODS = ['Puffin', 'APRL EC+F3', 'APRL EC+MM']
NAMES = ['Puffin', '能量守恒 + 三点核', '能量守恒 + 动量匹配']
COLORS = ['#2679ae', '#d37b28', '#57916f']
SHAPES = ['o', 's', '^']
STYLES = ['-', '--', '-.']


def read(path):
    return json.loads(Path(path).read_text())


def save(path, value):
    Path(path).write_text(json.dumps(value, ensure_ascii=False, indent=2) + '\n')


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def identity(path):
    s = Path(path).stat()
    return dict(size=s.st_size, mtime_ns=s.st_mtime_ns, inode=s.st_ino)


def verify(review):
    """Verify prior evidence without redoing physical postprocessing."""
    result = read(review / 'data/RESULT.json')
    validation = read(review / 'VALIDATION.json')
    assert result['complete'] and result['independent_source']
    assert validation['passed'] and result['methods'] == METHODS
    assert len(result['rows']) == 18 and len(result['initial']) == 3
    assert result['source_records'] == 1733612
    for name, expected in validation['hashes'].items():
        assert sha(review / 'data' / name) == expected, name
    for path, expected in read(review / 'data/PROTECTED_OUTPUTS.json').items():
        assert identity(path) == expected, path
    frozen = {}
    for name in ['six/FROZEN_PHYSICS.json', 'six-aprl/FROZEN.json']:
        for path, expected in read(review.parent / name).items():
            assert path not in frozen or frozen[path] == expected
            frozen[path] = expected
    for path, expected in frozen.items():
        assert sha(path) == expected, path
    with (review / 'data/growth.csv').open() as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == 18
    for row, saved in zip(result['rows'], rows):
        assert row['method'] == saved['method']
        assert row['segment'] == int(saved['segment'])
        for key in ['Hann_pm5pct_power_Eproxy_W', 'central_10lambda_rms_b1', 'mean_loss_eV']:
            assert row[key] == float(saved[key])
        assert row['full_cohort_available'] and row['records'] == result['source_records']
        assert row['missing_samples_reconstructed'] == 0
        assert math.isclose(row['power_denominator_s'], result['fixed_window_duration_s'], rel_tol=1e-15)
    return result, dict(
        reviewed_data_hashes_verified=len(validation['hashes']),
        closed_output_identities_verified=result['protected_closed_outputs'],
        frozen_asset_hashes_verified=len(frozen), scientific_H5_opened=False,
        reviewed_rows_and_CSV_exact=True)


def figures(result, target, font):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib import font_manager, ticker

    font_manager.fontManager.addfont(str(font))
    family = font_manager.FontProperties(fname=str(font)).get_name()
    plt.rcParams.update({
        'font.family': [family, 'DejaVu Sans'], 'font.size': 12,
        'axes.titlesize': 16, 'axes.labelsize': 13, 'legend.fontsize': 11,
        'axes.spines.top': False, 'axes.spines.right': False,
        'axes.unicode_minus': False, 'axes.formatter.use_mathtext': False,
        'text.parse_math': False, 'pdf.fonttype': 42, 'ps.fonttype': 42,
        'figure.facecolor': 'white', 'savefig.facecolor': 'white',
    })
    rendered = []
    plotted = []

    def finish(fig, name):
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        labels = []
        hidden_ticks = set()
        for ax in fig.axes:
            for axis, limits in [(ax.xaxis, ax.get_xlim()), (ax.yaxis, ax.get_ylim())]:
                low, high = sorted(limits)
                for tick in axis.get_major_ticks() + axis.get_minor_ticks():
                    if not low <= tick.get_loc() <= high:
                        hidden_ticks.update([id(tick.label1), id(tick.label2)])
        for label in fig.findobj(matplotlib.text.Text):
            if not label.get_visible() or not label.get_text():
                continue
            if id(label) in hidden_ticks:
                continue
            assert not label.get_usetex()
            assert '$' not in label.get_text() and '_' not in label.get_text()
            extent = label.get_window_extent(renderer)
            assert extent.x0 >= -2 and extent.y0 >= -2, label.get_text()
            assert extent.x1 <= fig.bbox.width + 2 and extent.y1 <= fig.bbox.height + 2, label.get_text()
            labels.append(label.get_text())
        for suffix in ['png', 'pdf']:
            path = target / f'{name}.{suffix}'
            fig.savefig(path, dpi=190)
            rendered.append(str(path.name))
        plt.close(fig)
        return labels

    def draw(ax, method, ykey, scale, initial=False):
        index = METHODS.index(method)
        rows = [r for r in result['rows'] if r['method'] == method]
        if initial:
            rows = [r for r in result['initial'] if r['method'] == method] + rows
        x = [r['actual_z_m'] for r in rows]
        y = [r[ykey] * scale for r in rows]
        line, = ax.plot(x, y, color=COLORS[index], marker=SHAPES[index],
                        linestyle=STYLES[index], linewidth=2, markersize=6, label=NAMES[index])
        assert list(line.get_xdata()) == x and list(line.get_ydata()) == y
        plotted.append(dict(method=method, metric=ykey, scale=scale, x=x, y=y))

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.4), layout='constrained')
    for method in METHODS:
        draw(axes[0], method, 'Hann_pm5pct_power_Eproxy_W', 1e-6)
        draw(axes[1], method, 'central_10lambda_rms_b1', 1, initial=True)
    axes[0].set(yscale='log', xlabel='实验室系位置（m）', ylabel='中央频带功率（MW）',
                title='六段出口的功率增长')
    axes[0].yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: f'{x:g}'))
    axes[1].set(xlabel='实验室系位置（m）', ylabel='局部 RMS 聚束因子',
                title='相同中央窗口内的聚束增长', ylim=(0, .7))
    for ax in axes:
        ax.set_xticks([0, 3, 6, 9, 12, 15])
        ax.grid(alpha=.2, which='major')
        ax.legend(loc='upper left', frameon=False)
    fig.suptitle('独立随机束流六段对照', fontsize=20)
    labels1 = finish(fig, 'power_bf_growth')

    fig, ax = plt.subplots(figsize=(9.5, 5.5), layout='constrained')
    for method in METHODS:
        draw(ax, method, 'mean_loss_eV', .001, initial=True)
    ax.axhline(0, color='#777777', linewidth=.8)
    ax.set(xlabel='实验室系位置（m）', ylabel='整束平均能损（keV / 电子）',
           title='六段漂移与波荡过程中的粒子能损')
    ax.set_xticks([0, 3, 6, 9, 12, 15])
    ax.grid(alpha=.2)
    ax.legend(loc='upper left', frameon=False)
    labels2 = finish(fig, 'particle_energy_loss')
    return dict(files=rendered, plotted_arrays=plotted, labels=labels1 + labels2,
                no_subscripts_or_superscripts=True, no_corner_badges=True,
                plot_data_exactly_from_reviewed_rows=True, PDF_separately_rendered=False)


def documents(result, out):
    last = {method: next(r for r in result['rows'] if r['method'] == method and r['segment'] == 6)
            for method in METHODS}
    p, f, m = [last[name] for name in METHODS]
    f_rel = 100 * (f['Hann_pm5pct_power_Eproxy_W'] / p['Hann_pm5pct_power_Eproxy_W'] - 1)
    m_rel = 100 * (m['Hann_pm5pct_power_Eproxy_W'] / p['Hann_pm5pct_power_Eproxy_W'] - 1)
    mf_rel = 100 * (m['Hann_pm5pct_power_Eproxy_W'] / f['Hann_pm5pct_power_Eproxy_W'] - 1)
    summary = [
        '# 六段 Puffin、三点核与动量匹配结果', '',
        '整理日期：2026-11-16。共同独立随机输入：1733612 个等权宏粒子，权重约 355，总电荷 98.6031 pC。', '',
        '## 主要结果', '',
        f'第六段中央频带功率：Puffin 8.1455 MW、三点核 9.0903 MW、动量匹配 8.5825 MW。相对 Puffin，三点核高 {f_rel:.2f}%，动量匹配高 {m_rel:.2f}%；动量匹配比三点核低 {abs(mf_rel):.2f}%。这些是相同窗口中的功率差异，不作为统一理论误差。', '',
        '| 方法 | 第六段中央功率 MW | 局部 RMS 聚束因子 | 整束平均能损 keV/电子 | 整束总粒子能损 µJ |',
        '|---|---:|---:|---:|---:|',
    ]
    for method, name in zip(METHODS, NAMES):
        r = last[method]
        summary.append(f"| {name} | {r['Hann_pm5pct_power_Eproxy_W']/1e6:.6f} | {r['central_10lambda_rms_b1']:.6f} | {r['mean_loss_eV']/1000:.6f} | {r['full_cohort_particle_loss_J']*1e6:.6f} |")
    summary += [
        '', '![功率与局部聚束因子](figures/zh-CN/power_bf_growth.png)', '',
        'Puffin 与三点核的局部聚束因子在第 5 段出口最高，第 6 段回落；动量匹配到第 6 段仍在上升。六个出口采样不能确定连续的饱和位置。', '',
        '![实验室系平均粒子能损](figures/zh-CN/particle_energy_loss.png)', '',
        '能损为正表示粒子失能，为负表示粒子增能。两种 PIC 方法在前四段平均能损为负，后两段转正，图中保留负值。其来源尚未用本组分项能量账本分解。', '',
        '中央频带功率只覆盖固定中央时间窗、±5% 频带和 ±520 µm 孔径；粒子能损覆盖整束。因此不能将两者直接相减作为能量守恒残差。PIC 原始全量账本已记录并闭合，但本包尚未新增其完整分项分析。', '',
        '## 任务对应关系', '',
        '| 方法 | 显示任务号 | Slurm 实体号 | 物理完成步数 | 运行时间 |',
        '|---|---|---|---:|---|',
        '| Puffin | 19795 | 19795 | 42720 | 5:44:20 |',
        '| 能量守恒 + 三点核 | 19796_0 | 19798 | 28199 | 18:44:24 |',
        '| 能量守恒 + 动量匹配 | 19796_1 | 19796 | 28199 | 28:50:22 |', '',
        '两组 PIC 物理进程均返回 0，到达配置停止面，探测器与账本完整关闭。Slurm 批次 FAILED 1:0 来自收尾脚本将实际 28199 步与预期 28198 步作严格相等检查；保留该原始状态，不将它改写成批次成功。', '',
        '## 图与记录', '',
        '- [中文功率与聚束图 PDF](figures/zh-CN/power_bf_growth.pdf)',
        '- [中文粒子能损图 PDF](figures/zh-CN/particle_energy_loss.pdf)',
        '- [英文功率与聚束原图](reviewed/data/power_bf_growth.png)',
        '- [英文粒子能损原图](reviewed/data/particle_energy_loss.png)',
        '- [逐段数值表](reviewed/data/comparison_by_segment.csv)',
        '- [完整增长统计](reviewed/data/growth.csv)',
        '- [入口聚束统计](reviewed/data/initial_bunching.csv)',
        '- [统计口径与公式说明](METHODS.md)',
        '- [原三方法分析报告](reviewed/README.md)',
        '- [原独立复核记录](reviewed/VALIDATION.json)',
        '- [原物理输出闭合复核](reviewed/data/CLOSURE_REVIEW.json)', '',
        '本包保留英文原图，中文图由已复核统计表重新绘制，不含上下标和角标。仅含图、结果与校验记录，不含源代码、大型场文件或粒子文件。',
    ]
    (out / 'README.md').write_text('\n'.join(summary) + '\n')
    methods = f'''# 统计口径与公式说明

## 1. 同一物理输入与比较范围

三组使用同一独立随机六维事件，宏粒子数 {result['source_records']}，实际等权 {result['source_weight']:.12f}，输入 SHA256 为 `{result['input_sha256']}`。六段磁结构、磁端、间隙、四极铁和输出位置一致。三点核和动量匹配是两种替代方法，不能把后者写作三点核再叠加动量匹配，也不由方法名称推断三维动量或实验室系能量账本已闭合。

Puffin 的入口由共同源弹道传播 1 cm 得到；APRL 同时计算这段的自场。Puffin 使用前向傍轴场和原生过滤；APRL 还包含纵向场和速度近场。两者模型差异没有由幅度拟合消除。本结果是该带权独立样本的对照，不能当作每电子一个粒子的 SASE 精度验收。

## 2. 固定中央时间窗口

迟滞时间定义为 tau = t − z/c。在名义抛物线后半束团中选取 I/I峰 ≥ 0.8，再去掉所选区间两端各 10%。每个面按粒子电荷加权中位迟滞时间平移，窗口相对中位值为 [{result['fixed_window_relative_to_charge_median_s'][0]*1e12:.12f}, {result['fixed_window_relative_to_charge_median_s'][1]*1e12:.12f}] ps，固定宽度 T = {result['fixed_window_duration_s']*1e12:.12f} ps。不按辐射亮度选窗，不使用固定粒子 ID 群代替同一个 Eulerian 时间窗口。

## 3. 中央频带功率

共振波长 lambda = {result['wavelength_m']*1e9:.9f} nm，f0 = c/lambda。对窗口中的横向电场乘连续 RMS 归一化 Hann 权重：

    w(tau) = sqrt(8/3) sin²[pi (tau − tau低) / T]
    F(f,x,y) = ∫ w(tau) E横(tau,x,y) exp(−2 pi i f tau) d tau

sqrt(8/3) 来自标准 Hann 权重的平均平方为 3/8，使其在连续固定窗口中具有单位均方权重；它不保证任意脉冲加窗前后能量相同。FFT 使用相应时间步长与 Parseval 归一化，对实场的正负频率成对计入。只保留 |f| 在 [0.95 f0, 1.05 f0] 中的频率，并在横向矩形 |x|、|y| ≤ 520 µm 内积分：

    U频带 = epsilon0 c ∫孔径 dA ∫选定正负频带 |F|² df
    P中央 = U频带 / T

能量除以同一固定 T，不除以各程序整段记录时长。APRL 保留已保存探测器的插值；Puffin 保留原生横向节点，不追加四节点平均。该定义从前向真空波的 S纵 = epsilon0 c |E横|² 得到，作为共同电场功率代理；对于含速度近场的 PIC，它并不等于已分离的纯辐射，也不是原生全频有符号 Poynting 通量。频带、孔径与窗口是比较定义，不是完整守恒控制体。

## 4. 聚束因子

同一中央窗口划分为长度 10 lambda 的局部小组。对每组电荷取正的权重 q：

    b组 = Σ q exp(−2 pi i c tau / lambda) / Σ q
    BF局部RMS = sqrt[Σ Q组 |b组|² / Σ Q组]

复相位相干叠加给出组内聚束，按组电荷平均平方再开根得到局部 RMS，避免不同小组相位相互抵消。全束 global BF 直接对整束复相位求和，含义不同，不能和这里的局部 RMS 或 Puffin 原生 1 lambda 切片幅度混称。入口值与六个出口值均保留在 CSV 中。功率只有六个出口样本，图中连线不增加采样信息。

## 5. 实验室系粒子能量与正负号

所有平均值按粒子电荷加权。共同源的平均 gamma 为 {result['source_gamma_mean']}；各出口包含完整宏粒子数和电荷，因此可以比较同一整束总能量。以共同源为基准：

    平均能损 = [共同源平均gamma − 出口平均gamma] m c²
    粒子平均能量变化 = −平均能损
    整束总能损[J] = (总电荷/e) × 平均能损[eV] × e[J/eV]

正能损表示粒子失能，负能损表示粒子增能；不能取绝对值后相加。使用同一物理输入作基准，不分别把三种方法的入口重置成零。PIC 磁入口相对共同源约有几 eV 的变化也保留。

整束平均能损不等于中央窗窄带辐射能量。要检验完整能量收支，还须按同一控制体及时间范围加入轴向端口、四侧能流、时间端面存量、粒子输运和相应离散算子项。本组已记录这些原始诊断，但当前三方法对照尚未完成其分项后处理，不能把未分析项当作零。

## 6. 数据与核验范围

本包直接复用 comparison-20261116-r2 的 18 行统计及 3 行入口统计，没有重跑物理或重新读取科学 H5。原数据覆盖完整粒子电荷，中央场窗没有补造缺失采样；原 Parseval 最大相对误差 {result['max_Parseval_relative']:.6g}。英文图和原复核 JSON/CSV 按 SHA256 逐字节保留。输出文件身份和生产冻结依赖在整理前后复核。软件闭合与数值一致性不代替网格收敛或物理模型验收。
'''
    (out / 'METHODS.md').write_text(methods)


def prepare(args):
    result, evidence = verify(args.review)
    args.out.mkdir(parents=True, exist_ok=False)
    copied = args.out / 'reviewed'
    copied.mkdir()
    for name in ['README.md', 'VALIDATION.json', 'SOURCE_REVIEW.json', 'JOB_STATUS.txt']:
        shutil.copy2(args.review / name, copied / name)
    shutil.copytree(args.review / 'data', copied / 'data')
    target = args.out / 'figures/zh-CN'
    target.mkdir(parents=True)
    rendering = figures(result, target, args.font)
    documents(result, args.out)
    save(args.out / 'RENDER.json', rendering)
    save(args.out / 'PROVENANCE.json', dict(
        source_review=str(args.review), date_Asia_Shanghai=datetime.now(ZoneInfo('Asia/Shanghai')).isoformat(),
        source_input_sha256=result['input_sha256'], no_new_physics=True,
        evidence=evidence, reviewed_files_sha256={str(p.relative_to(args.out)): sha(p)
            for p in copied.rglob('*') if p.is_file()},
        presentation_script_sha256=sha(__file__), original_English_figures_preserved=True))
    print(json.dumps(dict(prepared=str(args.out), evidence=evidence), ensure_ascii=False))


def finalize(args):
    assert args.png_reviewed, 'Inspect both delivered PNG files before finalizing'
    result, evidence = verify(args.review)
    provenance = read(args.out / 'PROVENANCE.json')
    assert sha(__file__) == provenance['presentation_script_sha256']
    for path, expected in provenance['reviewed_files_sha256'].items():
        assert sha(args.out / path) == expected, path
        assert sha(args.review / Path(path).relative_to('reviewed')) == expected, path
    links = []
    for path in args.out.rglob('*.md'):
        for target in re.findall(r'!?\[[^\]]*\]\(([^)]+)\)', path.read_text()):
            assert (path.parent / target.split('#')[0]).exists(), (path, target)
            links.append(target)
    commit = subprocess.check_output(['git', '-C', str(Path(__file__).resolve().parents[4]),
                                      'rev-parse', 'HEAD'], text=True).strip()
    save(args.out / 'SOURCE_COMMIT.json', dict(commit=commit, script_sha256=sha(__file__), pushed=False))
    save(args.out / 'VALIDATION.json', dict(passed=True, **evidence,
        copied_review_files_byte_identical=True, relative_Markdown_links_verified=len(links),
        PNG_visual_review_passed=True, PDF_separately_rendered=False,
        no_source_code_or_raw_scientific_outputs_in_archive=True,
        full_energy_ledger_postprocessing_complete=False))
    files = sorted(p for p in args.out.rglob('*') if p.is_file() and p.name != 'SHA256SUMS')
    forbidden = {'.py', '.h5', '.bin', '.raw', '.npy', '.npz', '.sbatch'}
    assert not any(p.suffix in forbidden for p in files)
    (args.out / 'SHA256SUMS').write_text(''.join(f'{sha(p)}  {p.relative_to(args.out)}\n' for p in files))
    with zipfile.ZipFile(args.archive, 'x', compression=zipfile.ZIP_DEFLATED) as archive:
        for path in sorted(args.out.rglob('*')):
            if path.is_file():
                archive.write(path, str(Path(args.out.name) / path.relative_to(args.out)))
    with zipfile.ZipFile(args.archive) as archive:
        assert archive.testzip() is None
        assert len(archive.infolist()) == len(files) + 1
        for path in files + [args.out / 'SHA256SUMS']:
            name = str(Path(args.out.name) / path.relative_to(args.out))
            assert hashlib.sha256(archive.read(name)).hexdigest() == sha(path)
    package = dict(complete=True, archive=str(args.archive), bytes=args.archive.stat().st_size,
                   sha256=sha(args.archive), directory=str(args.out), archive_CRC_passed=True,
                   all_members_SHA256_verified=True, source_commit=commit)
    save(args.out.parent / (args.out.name + '_PACKAGE.json'), package)
    print(json.dumps(package, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--review', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--font', type=Path)
    parser.add_argument('--archive', type=Path)
    parser.add_argument('--finalize', action='store_true')
    parser.add_argument('--png-reviewed', action='store_true')
    arguments = parser.parse_args()
    arguments.review = arguments.review.resolve()
    arguments.out = arguments.out.resolve()
    if arguments.finalize:
        assert arguments.archive is not None
        finalize(arguments)
    else:
        assert arguments.font is not None
        prepare(arguments)
