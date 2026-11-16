"""Assemble existing result records and explicit, qualified FEL estimates."""
import csv
import hashlib
import json
import math
import shutil
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.font_manager import FontProperties
from matplotlib.ticker import FixedLocator, FormatStrFormatter
from scipy.special import jv


def sha(p):
    h = hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda: f.read(4194304), b''):
            h.update(b)
    return h.hexdigest()


def read(p):
    return json.loads(Path(p).read_text())


def save(p, value):
    Path(p).parent.mkdir(parents=True, exist_ok=True)
    Path(p).write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + '\n')


def rows(p):
    with Path(p).open() as f:
        return list(csv.DictReader(f))


def table(p, values):
    p.parent.mkdir(parents=True, exist_ok=True)
    with p.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(values[0]))
        w.writeheader()
        w.writerows(values)


def xie(sigma, emit_n, label):
    # Round matched-beam surrogate, with the choice recorded in every row.
    I, IA = 53.27545382794002, 17045.09026574083
    gamma, K, lu, lam = 298.2608843802298, 2.060074987322757, .028, 4.913172351123507e-7
    eps = emit_n / gamma
    beta = sigma * sigma / eps
    aw = K / math.sqrt(2)
    xi = aw * aw / (2 * (1 + aw * aw))
    Aw = aw * (jv(0, xi) - jv(1, xi))
    rho = ((I / IA) * (lu * Aw / (2 * math.pi * sigma)) ** 2 / (2 * gamma) ** 3) ** (1 / 3)
    L0 = lu / (4 * math.pi * math.sqrt(3) * rho)
    d = L0 / (4 * math.pi * sigma * sigma / lam)
    e = L0 / beta * 4 * math.pi * eps / lam
    g = 4 * math.pi * L0 / lu * .0002
    terms = [.45*d**.57, .55*e**1.6, 3*g**2, .35*e**2.9*g**2.4,
             51*d**.95*g**3, 5.4*d**.7*e**1.9, 1140*d**2.2*e**2.9*g**3.2]
    return dict(round_surrogate=label, sigma_m=sigma, emit_n_m=emit_n,
                I_A=I, Alfven_current_A=IA, gamma=gamma, K=K, lambda_u_m=lu,
                lambda_r_m=lam, relative_rms_energy_spread=.0002,
                geometric_emit_m=eps, surrogate_matched_beta_m=beta,
                aw=aw, xi=xi, Aw=Aw, rho=rho, L0_power_gain_m=L0,
                eta_diffraction=d, eta_emittance=e, eta_energy_spread=g,
                **{f'Lambda_term_{i+1}': float(t) for i,t in enumerate(terms)},
                Lambda=sum(terms), Lg_power_gain_m=L0*(1+sum(terms)))


def main(workspace, output):
    ws, root = Path(workspace).resolve(), Path(output).resolve()
    root.mkdir(parents=True, exist_ok=False)
    template = Path(__file__).resolve().parent
    campaign = ws/'results/fel_independent_w355_20261114_r1'
    source, identities = {}, {}

    def protect(src):
        src = Path(src)
        s = src.stat()
        identities[str(src)] = dict(inode=s.st_ino, size=s.st_size, mtime_ns=s.st_mtime_ns)
        source[str(src)] = sha(src)

    def copy(src, dest):
        protect(src)
        dest = root/dest
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(src, dest)
        assert sha(dest) == source[str(src)]

    # Latest result files only: no cache, executable, H5, or source scripts.
    energy = campaign/'two-method-energy-curves-20261116-r1'
    for name in ['README.md','RESULT.json','VALIDATION.json','SOURCE_COMMIT.json','SOURCE_DATA_SHA256.json','PLOTS.json']:
        copy(energy/name, Path('01_energy_curves')/name)
    for sub in ['data','figures']:
        for p in sorted((energy/sub).rglob('*')):
            if p.is_file():
                assert p.suffix in ['.csv','.png','.pdf']
                copy(p, Path('01_energy_curves')/p.relative_to(energy))
    audit = campaign/'momentum-gap-audit-20261116-r1'
    for name in ['RESULT.json','VALIDATION.json','SOURCE_REVIEW.json','PAIRING.json',
                 'INDEPENDENT_PORT_CAP_CHECK.json','INDEPENDENT_SIDE_CHECK.json',
                 'PARTICLE_FREE_FIELD_CHECK.json','NATIVE_FORCE_AND_OPERATOR_DETAIL.json']:
        copy(audit/name, Path('01_energy_curves/audit_records')/name)
    copy(audit/'README.md', Path('01_energy_curves/audit_records/ORIGINAL_REPORT.txt'))
    (root/'01_energy_curves/audit_records/README.md').write_text(
        '# 动量账本差额的既有只读审计\n\n'
        '本目录保留原报告逐字文本和检查记录。端口/端盖、四侧积分已用独立算法重算；'
        '当前五次补偿方法不保证严格总动量不变量，原生电荷电力与电流动量配对仍有差异。'
        '未保存实际完整子步力，不能将全部 37.56% 唯一分摊给某条推进语句。\n\n'
        '[原报告文本](ORIGINAL_REPORT.txt)；[独立端面检查](INDEPENDENT_PORT_CAP_CHECK.json)；'
        '[独立侧面检查](INDEPENDENT_SIDE_CHECK.json)；[审计验证](VALIDATION.json)。'
        '原文中的工作区绝对路径仅作历史来源，不是本转交包的依赖。\n')

    growth = campaign/'comparison-20261116-r2'
    for name in ['comparison_by_segment.csv','growth.csv','initial_bunching.csv','RESULT.json']:
        copy(growth/'data'/name, Path('02_growth_comparison/data')/name)
    copy(growth/'VALIDATION.json', Path('02_growth_comparison/VALIDATION.json'))
    for lang in ['zh-CN','en']:
        for name in ['power','bunching','particle_energy_loss']:
            for ext in ['png','pdf']:
                copy(ws/f'results/fel_six_comparison_20261116_r3/figures/{lang}/{name}.{ext}',
                     Path(f'02_growth_comparison/figures/{lang}/{name}.{ext}'))
    (root/'02_growth_comparison/README.md').write_text(
        '# 当前独立 w355 六段三方法\n\n'
        '同一 1733612 条独立随机等权源，实际宏权重 354.999959。共同窗口宽 '
        '0.496626273 ps，相对每个面的电荷中位迟滞时间平移。Hann ±5% 频带，'
        '±520 µm 孔径，横向 E 场功率代理；BF 是同窗 10 波长局部电荷加权 RMS。'
        '粒子平均能损采用整束共同源基准，允许负值。功率和整束能损不是同一统计域。\n\n'
        '- [功率](figures/zh-CN/power.png)\n'
        '- [BF](figures/zh-CN/bunching.png)\n'
        '- [粒子平均能量损失](figures/zh-CN/particle_energy_loss.png)\n'
        '- [出口数据](data/comparison_by_segment.csv)\n'
        '- [入口数据](data/initial_bunching.csv)\n'
        '- [验证记录](VALIDATION.json)\n\n'
        '每张图独立单轴。英文图在 figures/en，同源 PDF 保留。'
        '动量匹配是五次补偿方法的沿用简称，不能据此推断严格总动量守恒。\n')

    theory = root/'03_theory_notes'
    for name in ['MING_XIE.md','BF_SPACE_CHARGE.md']:
        copy(template/name, Path('03_theory_notes')/name)
    oldsc = ws/'results/sparc_startup_space_charge_theory_20261113_r2'
    for name in ['README.md','RESULT.json','PROVENANCE.json','VERIFICATION.json','VISUAL_CHECK.json',
                 'curves.csv','cold_FEL_models.csv']:
        copy(oldsc/name, Path('03_theory_notes/space_charge')/name)
    startup = ws/'results/sparc_random_startup_20261112_r2/reviewed-r1'
    for name in ['REVIEW.json','VERIFICATION.json','ten_cell_summary.csv','pic_fixed_cohorts.csv',
                 'genesis_matched_endpoints.csv','summary.csv','band_energy.csv']:
        copy(startup/name, Path('03_theory_notes/startup_evidence')/name)
    copy(startup/'README.md', Path('03_theory_notes/startup_evidence/ORIGINAL_REPORT.txt'))
    bucket = ws/'results/sparc_single_seed_gap_20261113_r1'
    for name in ['comparison_summary.csv','matched_cohort_bunching.csv','RESULT.json','VERIFICATION.json']:
        copy(bucket/name, Path('03_theory_notes/startup_evidence/bucket_origin')/name)
    copy(bucket/'README.md', Path('03_theory_notes/startup_evidence/bucket_origin/ORIGINAL_REPORT.txt'))
    (theory/'startup_evidence/README.md').write_text(
        '# 旧独立随机首段 PIC / GENESIS 证据\n\n'
        '此处为旧 w512-a、w512-b、w64-a 配对，不是当前 w355/Puffin 六段统计。'
        '主 BF 使用共同固定粒子群、按生成箱分组的 10 波长局部 RMS，'
        '不同于当前随面平移的 Eulerian 高流区窗口。源均为独立事件，'
        '不属于修复前相关复制加载。三个出口 PIC/GENESIS-off BF 比为 '
        '0.91593、0.92942、0.93731；入口在每对中一致。\n\n'
        '[配对 BF 表](ten_cell_summary.csv)；[原始结论文本](ORIGINAL_REPORT.txt)；'
        '[验证记录](VERIFICATION.json)；[后续切片原点控制原文](bucket_origin/ORIGINAL_REPORT.txt)；'
        '[该控制的 BF 数据](bucket_origin/matched_cohort_bunching.csv)。'
        '切片控制改变了光学能量而没有消除 BF 差，不能将半片结果择作正确参考。'
        '原文本中外部工作区路径仅保留历史来源，'
        '不意味着相关大输出随包转交。本包没有用旧光学能量差或 BF 差证明空间电荷唯一归因。\n')

    oldx = ws/'results/sparc_sase_11segment_2slip_20261111_r1/analysis-energy-binomial3/theory-literature-comparison'
    for name in ['RESULT.json','VALIDATION.json','comparison.csv']:
        copy(oldx/name, Path('03_theory_notes/history')/name)
    copy(oldx/'README.md', Path('03_theory_notes/history/ORIGINAL_REPORT.txt'))
    (theory/'history/README.md').write_text(
        '# 历史谢明比较：相关复制加载，不能作当前独立随机验收\n\n'
        '原分析得到 0.667/0.778 m 斜率和 12.39 m 峰，使用修复前相关复制输入。'
        '后来加载审计发现同横向事件成列重复，不能作为当前独立随机六段或十一段结果。'
        '原文、CSV、JSON 字节保留；JSON 中 current 仅指当时旧工况。'
        '旧论文能量值还涉及论文版本区别，本次不将它作新的定量参照。\n\n'
        '[原文](ORIGINAL_REPORT.txt)；[原数值](RESULT.json)；[原表](comparison.csv)。'
        '当前理论比较请使用 [新谢明说明](../MING_XIE.md)。\n')
    audit_scope = ws/'results/puffin_six_same_ledger_20261113_r1/loading-audit-20261114-r1'
    for name in ['RESULT.json','QUALIFICATION.json']:
        copy(audit_scope/name, Path('03_theory_notes/history/loading_audit')/name)
    sat = campaign/'eleven-single-node-r1/saturation-20261116-r2'
    copy(sat/'SATURATION.json', Path('03_theory_notes/data/current_eleven_saturation.json'))
    copy(sat/'VALIDATION.json', Path('03_theory_notes/data/current_eleven_validation.json'))

    exits = rows(root/'02_growth_comparison/data/comparison_by_segment.csv')
    initial = rows(root/'02_growth_comparison/data/initial_bunching.csv')
    keys = {'Puffin':'Puffin_power_W','F3':'APRL_F3_power_W','MM':'APRL_MM_power_W'}
    gain = []
    for tag,key in keys.items():
        for a,b in zip(exits, exits[1:]):
            lg = 2.156/math.log(float(b[key])/float(a[key]))
            gain.append(dict(method=tag, from_segment=int(a['segment']), to_segment=int(b['segment']),
                             magnetic_length_m=2.156, actual_z_interval_m=float(b['lab_z_m'])-float(a['lab_z_m']),
                             start_power_W=float(a[key]), end_power_W=float(b[key]), effective_power_gain_length_m=lg,
                             difference_from_literature_0p67_percent=100*(lg/.67-1)))
    table(theory/'data/effective_gain_lengths.csv', gain)
    estimates = [xie(132e-6,2.5e-6,'x-direction round surrogate'),
                 xie(math.sqrt(132e-6*75e-6),math.sqrt(2.5e-6*2.9e-6),'geometric-mean round surrogate'),
                 xie(75e-6,2.9e-6,'y-direction round surrogate')]
    table(theory/'data/ming_xie_estimates.csv', estimates)
    bf = []
    for seg in range(7):
        if seg == 0:
            vals = [float(r['central_10lambda_rms_b1']) for r in initial]
            z = .01
        else:
            r = exits[seg-1]
            vals = [float(r[k]) for k in ['Puffin_local_rms_BF','APRL_F3_local_rms_BF','APRL_MM_local_rms_BF']]
            z = float(r['lab_z_m'])
        bf.append(dict(segment=seg, lab_z_m=z, Puffin_BF=vals[0], F3_BF=vals[1], MM_BF=vals[2],
                       F3_difference_percent=100*(vals[1]/vals[0]-1), MM_difference_percent=100*(vals[2]/vals[0]-1)))
    table(theory/'data/bf_relative_to_puffin.csv', bf)
    sc = read(theory/'space_charge/RESULT.json')
    summary = dict(no_new_physics=True, no_scientific_H5_read=True,
        causal_space_charge_attribution_complete=False, current_initial_BF_equal_within_percent=.000052,
        current_first_exit_BF_relative_percent={tag:bf[1][tag+'_difference_percent'] for tag in ['F3','MM']},
        old_F3_only_cold_mode_suppression_percent=100*(1-sc['varying_envelope']['F3_only_density_ratio']),
        circularization_is_illustrative_not_a_bound=True, round_beam_estimates=estimates,
        all_adjacent_gain_lengths=gain,
        not_one_electron_SASE_accuracy_acceptance=True,
        prior_transfer_packages_preserved=True)
    save(root/'RESULT.json',summary)

    plt.rcParams.update({'font.size':13,'axes.titlesize':17,'axes.labelsize':14,
        'axes.spines.top':False,'axes.spines.right':False,'axes.grid':True,
        'grid.alpha':.22,'grid.linewidth':.7,'lines.linewidth':2.3,
        'pdf.fonttype':42,'savefig.dpi':180,'axes.formatter.useoffset':False})
    font = ws/'results/physics_chinese_print_20261112_r1/DroidSansFallback.ttf'
    protect(font)
    plots = []
    curves = rows(theory/'space_charge/curves.csv')
    colors = ['#316AA5','#D78326','#258264']

    def finish(fig, ax, name, lang, series):
        ax.tick_params(labelsize=12)
        if lang == 'zh-CN':
            for obj in [ax.title,ax.xaxis.label,ax.yaxis.label]:
                obj.set_fontproperties(FontProperties(fname=str(font),size=obj.get_fontsize()))
            if ax.get_legend():
                for obj in ax.get_legend().texts:
                    obj.set_fontproperties(FontProperties(fname=str(font),size=11))
        fig.tight_layout(pad=1.8)
        fig.canvas.draw()
        assert len(fig.axes)==1
        visible_ticks = []
        for axis, limits in [(ax.xaxis,ax.get_xlim()),(ax.yaxis,ax.get_ylim())]:
            visible_ticks.extend(t.label1 for t in axis.get_major_ticks()
                                 if min(limits)<=t.get_loc()<=max(limits))
        for obj in [ax.title,ax.xaxis.label,ax.yaxis.label]+visible_ticks:
            if obj.get_visible() and obj.get_text():
                assert '$' not in obj.get_text() and '_' not in obj.get_text(),obj.get_text()
                b = obj.get_window_extent(fig.canvas.get_renderer())
                assert b.x0>=-1 and b.y0>=-1 and b.x1<=fig.bbox.width+1 and b.y1<=fig.bbox.height+1,obj.get_text()
        for ext in ['png','pdf']:
            p = theory/f'figures/{lang}/{name}.{ext}'
            p.parent.mkdir(parents=True,exist_ok=True)
            fig.savefig(p,facecolor='white')
        plots.append(dict(name=name,language=lang,axes=1,series=series))
        plt.close(fig)

    for lang in ['zh-CN','en']:
        zh = lang=='zh-CN'
        fig,ax=plt.subplots(figsize=(10.8,6.5))
        s=[]
        for tag,color,mark,label in zip(keys,colors,['o','s','^'],
                ['Puffin','三点核','动量匹配（五次补偿）'] if zh else ['Puffin','EC + F3','EC + quintic compensation']):
            selected=[r for r in gain if r['method']==tag]
            x=[r['to_segment'] for r in selected];y=[r['effective_power_gain_length_m'] for r in selected]
            ax.plot(x,y,marker=mark,color=color,label=label)
            s.append(dict(method=tag,x=x,y=y))
        ax.axhline(.67,color='#626B75',linestyle='--',label='文献参考 0.67 m' if zh else 'Literature reference: 0.67 m')
        ax.set_yscale('log');ax.set_ylim(.5,9)
        ax.yaxis.set_major_locator(FixedLocator([.5,.67,1,2,4,8]));ax.yaxis.set_major_formatter(FormatStrFormatter('%g'))
        ax.yaxis.set_minor_locator(FixedLocator([]))
        ax.set_xticks([2,3,4,5,6],['1 → 2','2 → 3','3 → 4','4 → 5','5 → 6'])
        ax.set(xlabel='相邻出口区间' if zh else 'Consecutive exit interval',
               ylabel='有效功率增益长度（m，对数轴）' if zh else 'Effective power gain length (m, log scale)',
               title='六段增长斜率与谢明文献尺度' if zh else 'Six-segment growth slopes and Ming Xie reference')
        ax.legend(loc='upper left',frameon=False,fontsize=11)
        finish(fig,ax,'gain_length',lang,s)

        fig,ax=plt.subplots(figsize=(10.8,6.5));s=[]
        names=['固定束斑：冷模式','变化束斑：冷模式','变化束斑与三点核：冷模式'] if zh else [
            'Fixed beam size: cold mode','Evolving beam size: cold mode','Evolving size with F3 only: cold mode']
        for key,color,style,label in zip(['constant_density_transfer','envelope_density_transfer','envelope_F3_only_density_transfer'],
                                      colors,['-','--','-.'],names):
            x=[float(r['z_m']) for r in curves];y=[100*float(r[key]) for r in curves]
            ax.plot(x,y,color=color,linestyle=style,label=label)
            s.append(dict(column=key,x=x,y=y))
        ax.set(xlim=(0,2.156),ylim=(90,100.6),
               xlabel='首段内位置（m）' if zh else 'Position inside the first undulator (m)',
               ylabel='相干密度模式传递（%）' if zh else 'Coherent density-mode transfer (%)',
               title='空间电荷单模估算：未作实测比例拟合' if zh else 'Space-charge single-mode estimate, without fitting measured ratios')
        ax.legend(loc='lower left',frameon=False,fontsize=11)
        finish(fig,ax,'space_charge_transfer',lang,s)

        fig,ax=plt.subplots(figsize=(10.8,6.5));s=[]
        for tag,color,mark,label in zip(['F3','MM'],colors[1:],['s','^'],
                ['三点核','动量匹配（五次补偿）'] if zh else ['EC + F3','EC + quintic compensation']):
            x=[r['segment'] for r in bf];y=[r[tag+'_difference_percent'] for r in bf]
            ax.plot(x,y,marker=mark,color=color,label=label)
            s.append(dict(method=tag,x=x,y=y))
        ax.axhline(0,color='#626B75',linestyle='--',linewidth=1.1)
        ax.set_xticks(range(7))
        ax.set(xlabel='段出口编号（0 为入口）' if zh else 'Undulator exit number (0 = entrance)',
               ylabel='相对 Puffin 的 BF 差异（%）' if zh else 'BF difference relative to Puffin (%)',
               title='入口 BF 一致，首段出口出现差异' if zh else 'Matched entrance BF; differences arise at the first exit')
        ax.legend(loc='lower left',frameon=False,fontsize=11)
        finish(fig,ax,'bf_relative',lang,s)
    save(theory/'PLOTS.json',dict(figures=plots,source_data_only=True,measured_BF_not_overlaid_on_cold_transfer=True))
    save(theory/'REFERENCES.json',dict(primary_sources=[
        dict(title='Ming Xie, PAC 1995',url='https://proceedings.jacow.org/p95/ARTICLES/TPG/TPG10.PDF'),
        dict(title='Huang and Kim, PRSTAB 2007',url='https://journals.aps.org/prab/pdf/10.1103/PhysRevSTAB.10.034801'),
        dict(title='Campbell et al., FEL2019 TUP050',url='https://proceedings.jacow.org/fel2019/papers/tup050.pdf'),
        dict(title='Wang, Litvinenko and Webb, PRSTAB 2012',url='https://journals.aps.org/prab/pdf/10.1103/PhysRevSTAB.15.120701'),
        dict(title='Geloni et al., DESY 07-087',url='https://arxiv.org/abs/0706.2280')],
        full_copyrighted_papers_not_embedded=True))
    (theory/'README.md').write_text(
        '# 理论与证据范围\n\n'
        '- [谢明完整公式、条件和当前对照](MING_XIE.md)\n'
        '- [首段 BF 与空间电荷：已有证据和缺项](BF_SPACE_CHARGE.md)\n'
        '- [完整空间电荷原推导](space_charge/README.md)\n'
        '- [旧独立随机 PIC/GENESIS 配对证据](startup_evidence/README.md)\n'
        '- [修复前相关加载的历史谢明分析](history/README.md)\n'
        '- [主要文献](REFERENCES.json)\n\n'
        '新图为中英文单轴 PNG/PDF，无上下标/角标。旧空间电荷模型、原验证和历史报告保持其日期及适用范围，'
        '不以新打包日期冒充新物理验证。space_charge 中的原记录会引用未随包复制的旧图名；'
        '展示请用本目录 figures 中从相同 CSV 重绘的独立曲线。'
        '谢明表中三种圆束化的差别只反映约定敏感性，不是本椭圆束的误差界。\n')
    (root/'README.md').write_text(
        '# 独立 w355 六段结果：能量账本、增长与理论说明\n\n'
        '2026-11-16 整理。只转交图、统计数据、公式推导和结果复核记录；'
        '没有新的物理仿真，不含代码、H5 或大缓存。此前双方法能量曲线 ZIP 未包含谢明与 BF 空间电荷推导，'
        '本包补齐这些内容，并保留最新三方法增长图与两方法完整能量曲线。\n\n'
        '| 内容 | 本包状态 | 结论边界 |\n'
        '|---|---|---|\n'
        '| 双方法能量曲线及账本 | 已收录中英文图、全部采样 CSV、原复核记录 | 固定 lab 管仍有 26.51% / 37.56% 差额 |\n'
        '| 三方法功率、BF、平均能损 | 已收录独立单图及入口/出口数据 | 共窗光学量与整束账本不是同一统计域 |\n'
        '| 谢明比较 | 已补齐公式、适用条件、三种明确圆束估算和当前出口斜率 | 不是有限脉冲椭圆束的精确验收 |\n'
        '| BF 与空间电荷 | 已补齐完整推导、模型曲线、配对观测与证据说明 | 若干百分点量级相容，尚未完成唯一因果归因 |\n\n'
        '- [两方法完整能量变化曲线](01_energy_curves/README.md)\n'
        '- [动量差额的只读审计](01_energy_curves/audit_records/README.md)\n'
        '- [六段三方法增长、BF 与平均能损](02_growth_comparison/README.md)\n'
        '- [谢明公式与本轮理论对照](03_theory_notes/MING_XIE.md)\n'
        '- [首段 BF 与空间电荷说明](03_theory_notes/BF_SPACE_CHARGE.md)\n'
        '- [理论材料总索引](03_theory_notes/README.md)\n\n'
        '入口 BF 几乎相同；首段出口三点核较 Puffin 低 6.92%，五次补偿低 5.16%。'
        '旧三点核冷模式理论的 5.17% 不能拿来直接验证五次补偿的 5.16%。'
        '所谓动量匹配沿用用户简称，实际为五次补偿，不等于已经证明总动量守恒。\n\n'
        '各图中英 PNG/PDF 均保留，CSV 为可复算来源。原历史文档及 JSON 中的绝对路径/旧图名仅作来源记录；'
        '本包的主 Markdown 导航全部使用包内相对链接。旧相关复制加载结论单独放在 history，'
        '不作为当前独立随机输入的验证。图片与推导材料已汇总；空间电荷的唯一归因没有被标成完成。\n')

    for name in ['fel_six_two_method_energy_curves_20261116_r1.zip','fel_six_puffin_F3_MM_bilingual_20261116_r3.zip']:
        protect(ws/'results'/name)
    save(root/'SOURCE_FILES_SHA256.json',source)
    save(root/'SOURCE_FILE_IDENTITIES.json',identities)
    print(json.dumps(dict(directory=str(root),source_files=len(source),new_standalone_figures=len(plots),
        round_surrogate_gain_lengths_m=[r['Lg_power_gain_m'] for r in estimates],
        first_exit_BF_relative_percent=summary['current_first_exit_BF_relative_percent']),indent=2))


if __name__=='__main__':
    main(*sys.argv[1:])
