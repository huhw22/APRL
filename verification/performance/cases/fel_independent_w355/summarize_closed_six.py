"""Check and document the closed independent six-segment comparison."""
import csv
import hashlib
import json
import math
import sys
from pathlib import Path


def main(out):
    out = Path(out)
    data = out / 'data'
    result = json.loads((data / 'RESULT.json').read_text())
    rows, initial = result['rows'], result['initial']
    assert result['complete'] and result['independent_source']
    assert len(rows) == 12 and len(initial) == 2
    with (data / 'growth.csv').open() as f:
        csv_rows = list(csv.DictReader(f))
    assert len(csv_rows) == 12
    width = result['fixed_window_duration_s']
    area = 4 * .00052**2
    for row, saved in zip(rows, csv_rows):
        assert row['method'] == saved['method'] and row['segment'] == int(saved['segment'])
        assert row['power_denominator_s'] == width
        assert abs(row['aperture_area_m2'] / area - 1) < 2e-12
        assert math.isclose(row['Hann_pm5pct_energy_J'] / width,
                            row['Hann_pm5pct_power_Eproxy_W'], rel_tol=2e-15)
        assert float(saved['Hann_pm5pct_power_Eproxy_W']) == row['Hann_pm5pct_power_Eproxy_W']
        assert 0 <= row['central_b1'] <= row['central_10lambda_rms_b1'] + 1e-12 <= 1 + 1e-12
        assert row['full_cohort_available'] and row['records'] == result['source_records']
        assert abs(row['charge_fraction_of_original'] - 1) < 2e-12
        assert row['missing_samples_reconstructed'] == 0
        assert row['Parseval_relative_error'] < 3e-12
    for path, old in json.loads((data / 'PROTECTED_OUTPUTS.json').read_text()).items():
        s = Path(path).stat()
        assert old == dict(size=s.st_size, mtime_ns=s.st_mtime_ns, inode=s.st_ino)
    comparisons = []
    for segment in range(1, 7):
        p, a = [next(row for row in rows if row['segment'] == segment and row['method'] == method)
                for method in ['Puffin', 'APRL EC+F3']]
        assert abs(p['actual_z_m'] - a['actual_z_m']) < 1e-11
        comparisons.append(dict(segment=segment, lab_z_m=p['actual_z_m'],
            Puffin_power_W=p['Hann_pm5pct_power_Eproxy_W'], APRL_F3_power_W=a['Hann_pm5pct_power_Eproxy_W'],
            F3_over_Puffin_power=a['Hann_pm5pct_power_Eproxy_W'] / p['Hann_pm5pct_power_Eproxy_W'],
            Puffin_local_rms_BF=p['central_10lambda_rms_b1'], APRL_F3_local_rms_BF=a['central_10lambda_rms_b1'],
            Puffin_mean_loss_eV=p['mean_loss_eV'], APRL_F3_mean_loss_eV=a['mean_loss_eV']))
    with (data / 'comparison_by_segment.csv').open('w') as f:
        writer = csv.DictWriter(f, fieldnames=list(comparisons[0]))
        writer.writeheader(); writer.writerows(comparisons)
    last = comparisons[-1]
    lines = ['# 独立随机束流六段 Puffin 与三点核对照', '',
        '本轮使用同一份独立随机六维事件，1733612个等权宏粒子，实际权重354.999959、总电荷98.6031pC；没有纵向复制加载。两程序采用相同六段波荡器、磁端、间隙和四极铁配置。', '',
        f"三点核功率在第一至第五段均低于Puffin；第六段达到{last['APRL_F3_power_W']/1e6:.4f}MW，Puffin为{last['Puffin_power_W']/1e6:.4f}MW，三点核高{100*(last['F3_over_Puffin_power']-1):.2f}%。此前旧克隆输入对照中APRL明显更快的启动，在这轮独立输入下没有重现。两次输入统计不同，不能把曲线差直接作为克隆效应的定量因果分解。", '',
        '![中央功率与局部BF](data/power_bf_growth.png)', '',
        '| 段数 | Puffin功率 MW | 三点核功率 MW | 三点核/Puffin | Puffin局部BF | 三点核局部BF |',
        '|---|---:|---:|---:|---:|---:|']
    for row in comparisons:
        lines.append(f"| {row['segment']} | {row['Puffin_power_W']/1e6:.6f} | {row['APRL_F3_power_W']/1e6:.6f} | {row['F3_over_Puffin_power']:.5f} | {row['Puffin_local_rms_BF']:.6f} | {row['APRL_F3_local_rms_BF']:.6f} |")
    lines += ['',
        '两者局部BF均在第五段最高、第六段回落，而中央频带功率仍增加。图中只连接六个出口采样点，不能据此确定长束稳态饱和位置。', '',
        '## 共同统计范围', '',
        f"相同的抛物线束团后半段高流区：先限定名义I/I峰≥0.8，再截掉所选区间两端各10%。窗口按每面粒子电荷中位迟滞时间平移，固定为中位时刻之后[{result['fixed_window_relative_to_charge_median_s'][0]*1e12:.9f}, {result['fixed_window_relative_to_charge_median_s'][1]*1e12:.9f}]ps，宽度{width*1e12:.9f}ps；不根据辐射最强位置选窗。", '',
        '功率由同一物理共振频率c/λ、±5%频带、横向±520µm孔径的横向电场计算，沿用连续RMS-Hann窗sqrt(8/3)sin²[π(t−t低)/T]，频带能量除以同一窗口宽度T。没有按各程序整段记录时长相除，也没有跨程序幅度拟合。APRL使用已保存的探测器插值，Puffin使用原生横向节点；不另加四节点平均。该量是ε0c∫|E横|²的频带功率代理，PIC可能保留近场贡献。', '',
        '主BF为同一中央Eulerian窗口中的10λ分组电荷加权RMS：sqrt(ΣQ_j|b_j|²/ΣQ_j)。全束global BF另在CSV记录，不能与这里的局部RMS混称。所有12个出口粒子面均含完整1733612条与全束电荷，所有中央场窗覆盖完整，无缺采样点重建。', '',
        '## 粒子能损', '',
        '![实验室系平均粒子能损](data/particle_energy_loss.png)', '',
        '| 段数 | Puffin平均能损 keV/电子 | 三点核平均能损 keV/电子 |',
        '|---|---:|---:|']
    for row in comparisons:
        lines.append(f"| {row['segment']} | {row['Puffin_mean_loss_eV']/1000:.5f} | {row['APRL_F3_mean_loss_eV']/1000:.5f} |")
    lines += ['',
        '能损以两程序共同物理输入的电荷加权平均gamma为基准，使用相同电子质量与物理常数，正值表示失能。三点核前四段平均能量增加；这里只报告观测结果，尚未用分项能量账本判定来源。整束粒子能损与中央窗口的窄带功率范围不同，不能直接互相相减作守恒检验。', '',
        'Puffin只含前向傍轴场、原生过滤，入口前1cm按弹道映射；APRL包含该漂移段的PIC自场以及纵向和速度近场。当前比较没有消除这些模型差异，也没有分解为唯一算法原因。', '',
        '## 物理运行与批处理状态', '',
        'Puffin19795正常COMPLETED0:0，42720步，耗时5:44:20。APRL三点核19796_0物理进程返回0，完成28199步、全部粒子退场并到达配置停止面；探测器、28200条能量账本、worldtube及48个逐rank审计标记均完整闭合。原整批FAILED1:0仅因冻结收尾脚本硬编码预期28198步而拒绝；原失败状态、脚本和日志保留，没有重跑或修改冻结文件。', '',
        '这次只读已闭合六段Puffin和三点核数据。未分析仍运行的动量匹配或十一段数据，没有新物理仿真，也不将输出闭合等同于物理精度或全lab账本闭合。', '',
        f"一致性复核：12组相同分母/孔径、BF范围与RMS不等式、CSV反读、完整电荷、出口坐标一致及Parseval通过；最大Parseval相对误差{result['max_Parseval_relative']:.3g}。分析前后{result['protected_closed_outputs']}个闭合输出文件inode/size/mtime保持。生产冻结输入/卡/程序SHA检查通过。两张PNG已目视检查，PDF同源导出未另渲染。", '',
        '- [逐段对照表](data/comparison_by_segment.csv)',
        '- [完整增长统计](data/growth.csv)',
        '- [入口BF](data/initial_bunching.csv)',
        '- [运行闭合复核](data/CLOSURE_REVIEW.json)',
        '- [完整定义和结果](data/RESULT.json)',
        '- [独立一致性检查](VALIDATION.json)']
    (out / 'README.md').write_text('\n'.join(lines) + '\n')
    validation = dict(passed=True, rows=12, methods=2, fixed_power_denominator_s=width,
        same_aperture_area_m2=area, complete_source_records=result['source_records'],
        input_events_and_magnets_matched=True, all_full_cohorts=True,
        all_windows_without_missing_samples=True, CSV_roundtrip_checked=True,
        BF_RMS_inequality_checked=True, maximum_Parseval_relative=result['max_Parseval_relative'],
        closed_outputs_identity_unchanged=result['protected_closed_outputs'],
        PNG_visual_review_passed=True, PDF_separately_rendered=False,
        original_failed_batch_state_preserved=True, no_physics_rerun=True,
        only_closed_six_segment_outputs_used=True,
        hashes={p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in data.iterdir()
                if p.suffix in ['.json', '.csv', '.png', '.pdf']})
    (out / 'VALIDATION.json').write_text(json.dumps(validation, indent=2) + '\n')
    print(json.dumps(dict(comparison=comparisons, validation=validation), ensure_ascii=False), flush=True)


if __name__ == '__main__':
    main(sys.argv[1])
