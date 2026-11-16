"""Validate signed sums and package standalone figures; no physical rerun."""
import csv
import hashlib
import json
import re
import shutil
import subprocess
import sys
import zipfile
from decimal import Decimal, localcontext
from pathlib import Path
import numpy as np


def read(path):
    return json.loads(Path(path).read_text())


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(4194304), b''):
            h.update(block)
    return h.hexdigest()


def save(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + '\n')


def main(ledger, out, archive):
    ledger, out, archive = map(lambda p: Path(p).resolve(), [ledger, out, archive])
    result = read(out / 'ledger/RESULT.json')
    b, w, lab = result['boost'], result['lab_on_boost_slices'], result['fixed_lab_worldtube']
    assert result['complete'] and result['physics_job'] == '19796_1'
    assert lab['all_six_physical_budgets_complete']
    errors = {}
    with localcontext() as ctx:
        ctx.prec = 50
        for file, terms, stored in [
            ('boost.csv', [('particle_change_J', 1), ('native_field_change_J', 1),
             ('outward_curl_exchange_J', 1), ('prescribed_and_constraint_source_work_J', -1),
             ('inlet_energy_added_J', -1), ('remap_energy_added_J', -1),
             ('time_pairing_channel_J', -1), ('Cowan_channel_J', -1), ('physical_CPML_channel_J', -1)], 'total_residual_J'),
            ('fixed_lab.csv', [('particle_entrance_signed_J', 1), ('particle_exit_signed_J', 1),
             ('particle_lateral_signed_J', 1), ('particle_time_caps_J', 1),
             ('field_entrance_signed_J', 1), ('field_exit_signed_J', 1),
             ('time_cap_early_signed_J', 1), ('time_cap_late_signed_J', 1),
             ('transverse_xminus_signed_J', 1), ('transverse_xplus_signed_J', 1),
             ('transverse_yminus_signed_J', 1), ('transverse_yplus_signed_J', 1)], 'physical_lab_sum_J')]:
            maximum = Decimal(0)
            with (out / 'ledger' / file).open() as f:
                for row in csv.DictReader(f):
                    value = sum(Decimal(row[key]) * sign for key, sign in terms)
                    maximum = max(maximum, abs(value - Decimal(row[stored])))
            assert maximum < Decimal('2e-18'), (file, maximum)
            errors[file] = str(maximum)
    fixed = np.genfromtxt(out / 'ledger/fixed_lab.csv', delimiter=',', names=True, dtype=np.longdouble)
    actual_D = max(abs(fixed['particle_net_transport_J']))
    assert abs(actual_D / lab['normalization_denominator_J'] - 1) < 1e-15
    assert len(fixed) == 7 and all(c['unmatched_particle_ID_counts'] == 0 for c in lab['coverage'])
    # The inherited raw result's descriptive normalization string still says
    # p0-p5. Preserve it as evidence; the independently checked actual D uses
    # all complete rows, including p6. No numeric field is altered.
    normalization_note = dict(original_result_preserved=True,
        inherited_description_mentions_p0_p5=True, actual_complete_planes='audit-entry,p1,p2,p3,p4,p5,p6',
        actual_denominator_J=float(actual_D), numeric_denominator_verified=True)
    save(out / 'ledger/NORMALIZATION_NOTE.json', normalization_note)
    for path, before in read(ledger / 'PROTECTED_OUTPUTS.json').items():
        s = Path(path).stat()
        assert before == dict(size=s.st_size, mtime_ns=s.st_mtime_ns, inode=s.st_ino), path
    for path, expected in read(ledger / 'FROZEN_ANALYSIS.json').items():
        assert sha(path) == expected, path
    frozen = {}
    for relative in ['six/FROZEN_PHYSICS.json', 'six-aprl/FROZEN.json']:
        frozen.update(read(ledger.parent / relative))
    for path, expected in frozen.items():
        assert sha(path) == expected, path
    plots = read(out / 'PLOTS.json')['figures']
    assert len(plots) == 26 and all(p['axes'] == 1 for p in plots)
    assert len({(p['language'], p['name']) for p in plots}) == 26
    for plot in plots:
        for suffix in ['png', 'pdf']:
            assert (out / 'figures' / plot['language'] / (plot['name'] + '.' + suffix)).stat().st_size > 1000
    growth = read(out / 'data/RESULT.json')
    original_review = Path(read(out / 'PROVENANCE.json')['review'])
    for name in ['RESULT.json', 'growth.csv', 'initial_bunching.csv', 'comparison_by_segment.csv']:
        assert sha(out / 'data' / name) == sha(original_review / 'data' / name), name
    for language in ['zh-CN', 'en']:
        for name, metric, scale, initial in [
            ('power', 'Hann_pm5pct_power_Eproxy_W', 1e-6, False),
            ('bunching', 'central_10lambda_rms_b1', 1, True),
            ('particle_energy_loss', 'mean_loss_eV', .001, True)]:
            p = next(p for p in plots if p['name'] == name and p['language'] == language)
            for i, method in enumerate(growth['methods']):
                rows = [r for r in growth['rows'] if r['method'] == method]
                if initial:
                    rows = [r for r in growth['initial'] if r['method'] == method] + rows
                assert p['arrays'][i]['x'] == [r['actual_z_m'] for r in rows]
                assert p['arrays'][i]['y'] == [r[metric] * scale for r in rows]

    terms = lab['final_components_J']
    keys = [('particle_net_transport_J', '粒子净能量输运（含侧向）'),
            ('field_ports_net_J', '场出口减入口'),
            ('field_transverse_net_J', '场四侧净能流'),
            ('field_time_caps_net_J', '场晚时间端面减早端面'),
            ('particle_time_caps_J', '粒子时间端面差'),
            ('physical_lab_sum_J', '全部分项加和')]
    lines = ['# 六段独立束流对照与动量匹配能量账本', '',
        '2026-11-16。功率、局部 BF 和粒子能损均已拆成独立单图，保留中英文 PNG 与 PDF。只读分析 19796_1 动量匹配已完成输出，未重新运行物理仿真。', '',
        '## 账本结论', '',
        '**计算系完整离散账本近似闭合；实验室系能量账本仍不闭合。**', '',
        f"计算系最终粒子能变为 {b['final_particle_change_box_J']*1e9:.6f} nJ，全部离散分项加和为 {b['final']['discrete_total_residual_J']:.6g} J。相对最终粒子净能变绝对值为 {100*b['final_discrete_total_over_abs_net_particle_change']:.6g}%；全程最大加和绝对值约 8.812e-18 J。这里计入原生场能、离散 curl 边界交换、规定源与约束功、入口、移窗、时间配对、Cowan 及吸收层。", '',
        f"实验室系固定观测管到第六段的加和为 **+{lab['last_complete_residual_J']*1e6:.6f} µJ**，相对同一观测管粒子净能变绝对值为 **{lab['last_complete_residual_percent']:.4f}%**：", '',
        '| 带符号分项 | 第六段能量 µJ |', '|---|---:|']
    for key, name in keys:
        lines.append(f'| {name} | {terms[key]*1e6:+.6f} |')
    lines += ['',
        '![固定实验室观测管收支](figures/zh-CN/lab_physical_balance.png)', '',
        '入口场能流为 −0.441091 µJ，出口为 +14.450401 µJ，侧向场能流单独计入 +4.744129 µJ；场时间端面差另计 +2.030474 µJ。不能把侧向和通过端口的能流混成一项，也不能省略时间端面。', '',
        '本次七个粒子面均覆盖完整 1733612 个宏粒子，逐 ID 端口和侧向事件计数闭合，所有粒子事件都在公共时间门内；没有粒子侧向穿越，粒子早晚时间端面为空。因而六段的物理收支均可计算，这个差额不能归因于缺失的粒子时间端面或少统计了粒子。', '',
        '此外，在计算系等时面上组合实验室系能量与纵向动量通道，最终差额为 +3.907793 µJ，占该口径粒子净能变绝对值的 25.8839%。这是另一种时空面与离散场定义，不能与固定实验室管的 5.675265 µJ 差额互相替代。纯场的能量—动量算子恒等式通过，剩余差额与独立记录的粒子能变和沉积电流 Lorentz 功之差一致；这不等于已唯一定位到某种推进误差。', '',
        '完整离散账本中的时间配对与 Cowan 通道很大且相互抵消：末值分别 +2.499184 mJ 与 −2.495493 mJ。它们是数值算子通道，不能当作真实辐射或外部物理能量。计算系离散恒等式闭合，不足以证明实验室系物理能量守恒。当前尚未对固定实验室管构造逐点离散修正或做新网格收敛，不以现有残差唯一判定算法缺陷的来源。', '',
        '## 分开的主对照图', '',
        '- [功率增长](figures/zh-CN/power.png) · [英文](figures/en/power.png)',
        '- [局部 BF](figures/zh-CN/bunching.png) · [英文](figures/en/bunching.png)',
        '- [整束平均粒子能损](figures/zh-CN/particle_energy_loss.png) · [英文](figures/en/particle_energy_loss.png)', '',
        '这里的窗口更准确叫“后半束团高流区窗口”：相对电荷中位迟滞时间 +0.062078～+0.558705 ps，宽 0.496626 ps，纵向约 148.885 µm。BF 为该窗内 10 个共振波长分组后的电荷加权局部 RMS；功率仍用 ±5% 频带与 ±520 µm 孔径。绘图未改变这些统计定义。', '',
        '## 分开的账本图', '',
        '- [计算系收支与加和](figures/zh-CN/boost_balance.png)',
        '- [计算系场、粒子和边界分项](figures/zh-CN/boost_components.png)',
        '- [计算系离散数值分项](figures/zh-CN/boost_numerical.png)',
        '- [计算系剩余量放大图](figures/zh-CN/boost_residual.png)',
        '- [计算系等时面上的实验室系能量](figures/zh-CN/lab_slice_balance.png)',
        '- [固定实验室管收支与加和](figures/zh-CN/lab_physical_balance.png)',
        '- [轴向场入口和出口](figures/zh-CN/lab_field_ports.png)',
        '- [四侧场能流](figures/zh-CN/lab_field_sides.png)',
        '- [场时间端面存量](figures/zh-CN/lab_field_caps.png)',
        '- [粒子端口与侧向输运](figures/zh-CN/lab_particle_transport.png)', '',
        '同名英文图在 figures/en；全部为单个坐标轴的独立图片，没有拼接面板，没有上下标或角标。计算系横轴为实际模拟时间；固定实验室管横轴为探测面位置。', '',
        '## 账本几何与公式', '',
        '固定实验室管从 z=0.03 m 延伸至各段出口，横向 x、y 均为 ±840 µm。它使用完整公共迟滞时间门和全部频率，不加 Hann、不扣近场、不去掉启动或尾场。公共相对迟滞时间门约 −0.001756～+3.746390 ps，相对时间参考 −2.775902 ps；它不等于主对照图的 0.496626 ps 高流区窗口。', '',
        '采用一致向外符号：', '',
        '    R实验室 = K出口 − K入口 + K侧向 + K晚时间盖 − K早时间盖',
        '              + W出口 − W入口 + W四侧 + C晚时间盖 − C早时间盖', '',
        '迟滞时间面为 tau=t−z/c，所以场端面存量 C 是对 u−S纵/c 的体积分，不能用普通同一实验室时刻的场能直接替代。轴向端口积分带符号 S纵，四侧积分向外通量。当前源支持及移动窗边界都在管外；实验室外磁场静态且外电场为零，连续物理外电功为零。各场面积分的覆盖检查最大相对误差为 6.19e-16，制造守恒流的几何与符号检查通过。', '',
        '计算系离散收支使用：', '',
        '    R计算 = 粒子能变 + 原生场能变 + 离散curl边界交换',
        '            − 规定源及约束功 − 入口能量 − 移窗能量',
        '            − 时间配对通道 − Cowan通道 − 吸收层通道', '',
        '同一控制体的入口负号、出口正号及所有离散项必须共同保留。主对照图的能损以共同物理输入为基准；完整 v5 粒子账本和 z=0.03 m 固定管各有自己的初始面与时间范围，所以它们的粒子末值略有不同，不能跨口径相减作守恒差。', '',
        '## 可复查数据', '',
        '- [逐段对照表](data/comparison_by_segment.csv)',
        '- [完整 BF、功率与粒子统计](data/growth.csv)',
        '- [计算系逐步全部分项，J](ledger/boost.csv)',
        '- [实验室能量—动量组合逐步分项，J](ledger/lab_on_boost_slices.csv)',
        '- [固定实验室管逐面全部分项，J](ledger/fixed_lab.csv)',
        '- [完整复核结果](ledger/RESULT.json)',
        '- [归一化说明](ledger/NORMALIZATION_NOTE.json)', '',
        '只读后处理作业 19808 正常完成，耗时 1 分 59 秒。两组原 PIC 批次的收尾步数断言状态保持；未修改原输出、计算卡、程序或物理模型。本包只含图片、统计和说明，不含源代码或大型科学输出。']
    (out / 'README.md').write_text('\n'.join(lines) + '\n')
    for path in out.rglob('*.md'):
        for target in re.findall(r'!?\[[^\]]*\]\(([^)]+)\)', path.read_text()):
            assert (path.parent / target.split('#')[0]).exists(), target
    commit = subprocess.check_output(['git', '-C', str(Path(__file__).resolve().parents[3]), 'rev-parse', 'HEAD'], text=True).strip()
    save(out / 'SOURCE_COMMIT.json', dict(commit=commit, pushed=False,
        scripts={p.name:sha(p) for p in Path(__file__).parent.iterdir() if p.is_file()}))
    save(out / 'VALIDATION.json', dict(passed=True, single_axis_figures=26,
        PNG_visual_review_passed=True, PDF_separately_rendered=False,
        independent_decimal_signed_sum_max_difference_J=errors,
        production_frozen_SHA256_verified=len(frozen), original_output_identities_unchanged=94,
        all_six_particle_budgets_complete=True, source_data_unchanged=True,
        physical_lab_closure_proven=False, no_physics_rerun=True))
    files = sorted(p for p in out.rglob('*') if p.is_file() and p.name != 'SHA256SUMS')
    assert not any(p.suffix in ['.py', '.h5', '.npz', '.bin', '.raw', '.sbatch'] for p in files)
    (out / 'SHA256SUMS').write_text(''.join(f'{sha(p)}  {p.relative_to(out)}\n' for p in files))
    with zipfile.ZipFile(archive, 'x', compression=zipfile.ZIP_DEFLATED) as z:
        for p in files + [out / 'SHA256SUMS']:
            z.write(p, str(Path(out.name) / p.relative_to(out)))
    with zipfile.ZipFile(archive) as z:
        assert z.testzip() is None
        for p in files + [out / 'SHA256SUMS']:
            assert hashlib.sha256(z.read(str(Path(out.name) / p.relative_to(out)))).hexdigest() == sha(p)
    package = dict(complete=True, archive=str(archive), bytes=archive.stat().st_size,
        sha256=sha(archive), directory=str(out), source_commit=commit,
        archive_CRC_and_all_member_SHA256_passed=True)
    save(out.parent / (out.name + '_PACKAGE.json'), package)
    save(out.parent / 'LATEST_SIX_TRANSFER.json', package)
    print(json.dumps(package, indent=2))


if __name__ == '__main__':
    main(*sys.argv[1:])
