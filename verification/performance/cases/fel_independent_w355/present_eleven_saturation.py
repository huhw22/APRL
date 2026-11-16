"""Validate and present the closed eleven-segment saturation observables."""
import csv
import hashlib
import json
import sys
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main(out):
    out = Path(out)
    data = out/'data'
    r = json.loads((data/'RESULT.json').read_text())
    rows, dense = r['segments'], r['native_history']
    assert r['complete'] and len(rows) == 11
    assert [x['segment'] for x in rows] == list(range(1,12))
    assert [x['step'] for x in dense if x['step']%80 == 0] == list(range(0,79201,80))
    with (data/'gain.csv').open() as f: saved = list(csv.DictReader(f))
    assert len(saved) == 11
    for a,b in zip(rows,saved):
        assert a['Hann_pm5pct_power_Eproxy_W'] == float(b['Hann_pm5pct_power_Eproxy_W'])
        assert a['power_denominator_s'] == r['fixed_window_duration_s']
        assert abs(a['aperture_area_m2']/(4*.00052**2)-1) < 2e-12
        assert a['records'] == r['source_records'] and abs(a['charge_fraction_of_original']-1)<2e-12
        assert 0<=a['central_b1']<=a['central_10lambda_rms_b1']+1e-12<=1+1e-12
        assert a['missing_samples_reconstructed']==0 and a['Parseval_relative_error']<3e-12
        assert abs(a['plane_z_roundoff_difference_m'])<1e-9
    protected = json.loads((data/'PROTECTED_OUTPUTS.json').read_text())
    for name, old in protected.items():
        s=Path(name).stat()
        assert old==dict(size=s.st_size,mtime_ns=s.st_mtime_ns,inode=s.st_ino)
    peak = r['native_peak_candidates'][0]
    plan=json.loads((out.parent/'PLAN.json').read_text())
    magnets=[x for x in plan['magnetic_lattice'] if x['type']=='planar-undulator']
    peak_segment=next(i for i,m in enumerate(magnets,1)
                      if m['entrance_z']<=peak['lab_z_m']<=m['entrance_z']+m['period']*m['periods'])
    peak_native_row=next(x for x in dense if x['step']==peak['step'])
    assert peak_native_row['lab_z_m']==peak['lab_z_m']
    assert peak_native_row['central_Hann_native_allgrid_allband_power_W']==peak['power_W']
    band = max(rows,key=lambda x:x['Hann_pm5pct_power_Eproxy_W'])
    crossing = next(x for x in rows if x['Hann_pm5pct_power_Eproxy_W'] >= .9*band['Hann_pm5pct_power_Eproxy_W'])
    previous = rows[crossing['segment']-2]
    summary=dict(complete=True,physical_acceptance=False,definition='First prominent power maximum; exponential-growth rollover and band plateau are reported separately.',
        first_prominent_native_peak_lab_z_m=peak['lab_z_m'],first_prominent_native_peak_power_W=peak['power_W'],
        first_prominent_native_peak_undulator=peak_segment,native_periodic_sampling_interval_m=.028,
        peak_native_quantity='Same temporal Hann gate; full transverse grid and all retained frequencies; not the aperture/band power.',
        highest_sampled_band_exit=band['segment'],highest_sampled_band_exit_z_m=band['actual_z_m'],
        highest_sampled_band_power_W=band['Hann_pm5pct_power_Eproxy_W'],
        band_plateau_exit_segments=[6,7,8],
        first_exit_at_90pct_of_sampled_band_max=crossing['segment'],
        bracket_for_90pct_crossing_m=[previous['actual_z_m'],crossing['actual_z_m']],
        ninety_percent_is_an_optional_marker_not_a_universal_saturation_definition=True,
        end_band_power_W=rows[-1]['Hann_pm5pct_power_Eproxy_W'],
        exact_dense_aperture_band_peak_unavailable=True)
    (out/'SATURATION.json').write_text(json.dumps(summary,indent=2)+'\n')
    fig,axes=plt.subplots(1,2,figsize=(12.5,4.7),layout='constrained')
    axes[0].semilogy([a['actual_z_m'] for a in rows],[a['Hann_pm5pct_power_Eproxy_W']/1e6 for a in rows],'o-',color='#c8792e')
    axes[0].set(xlabel='Lab z (m)',ylabel='Central resonant-band power (MW)',title='Hann ±5%; aperture ±520 µm; exit samples')
    axes[0].annotate(f"U8 exit: {band['Hann_pm5pct_power_Eproxy_W']/1e6:.2f} MW",(band['actual_z_m'],band['Hann_pm5pct_power_Eproxy_W']/1e6),
                     xytext=(-90,-38),textcoords='offset points',arrowprops=dict(arrowstyle='->'),fontsize=9)
    axes[1].plot([a['lab_z_m'] for a in dense],[a['central_Hann_native_allgrid_allband_power_W']/1e6 for a in dense],color='#517f9f')
    axes[1].plot(peak['lab_z_m'],peak['power_W']/1e6,'o',color='#c8792e')
    axes[1].annotate(f"First strong peak: {peak['lab_z_m']:.2f} m (U8)",(peak['lab_z_m'],peak['power_W']/1e6),
                     xytext=(-140,20),textcoords='offset points',arrowprops=dict(arrowstyle='->'),fontsize=9)
    axes[1].set(xlabel='Lab z (m)',ylabel='Native central power (MW)',title='Full grid / all retained frequencies; per-period',xlim=(10,28),ylim=(0,peak['power_W']/1e6*1.22))
    for ax in axes:
        for a in rows:ax.axvline(a['actual_z_m'],color='gray',alpha=.2,lw=.5)
        ax.grid(alpha=.2,which='both')
    fig.suptitle('Puffin | independent w ≈ 355 | eleven segments, four-slip bunch\nCommon central-window selection rule; fixed 0.993 ps gate',fontsize=12)
    fig.savefig(data/'saturation.png',dpi=180);fig.savefig(data/'saturation.pdf');plt.close(fig)
    lines=['# 十一段 Puffin 饱和位置', '',
        '**以首次显著功率峰定位，位于第8段，约z=19–20m。** 每周期原生功率峰在19.358m；按先前Hann±5%、±520µm口径，11个出口中第8段出口20.058m最高，为6.66285MW。', '',
        '增长从第6段末约15m已明显转缓，第7–8段形成带限功率平台。这里区分进入非线性饱和过渡与达到首个显著峰：如果要求给出单个可复核标记，可用原生首峰19.358m并注明统计量；不能把它说成已精确量得同孔径共振带的连续峰位置。', '',
        '![十一段增长与首次显著峰](data/saturation.png)', '',
        '| 段数 | 实际出口z m | Hann共振频带功率 MW | 中央10λ RMS BF |',
        '|---|---:|---:|---:|']
    for a in rows:lines.append(f"| {a['segment']} | {a['actual_z_m']:.5f} | {a['Hann_pm5pct_power_Eproxy_W']/1e6:.6f} | {a['central_10lambda_rms_b1']:.6f} |")
    lines+=['',
        '第6→7段带限功率从5.95854升至6.57585MW，第7→8段仅再增加约1.32%；第9段回落到4.19098MW，第11段回升至4.97301MW。后续表现为饱和后的回落/再增长，不能把最大末端位置当作新的指数增益长度。', '',
        '## 统计范围', '',
        f"与六段对照使用同一选窗规则：抛物线束团后半高流区I/Ipk≥0.8，再截掉所选段两端各10%，每面按电荷中位迟滞时间平移。四倍滑移长束下固定门宽为{r['fixed_window_duration_s']*1e12:.9f}ps，是六段短束窗口的两倍；不将短束曲线拼接成长束结果。", '',
        '左图按实际粒子面中位时间取窗，对原生节点横向实场使用连续RMS-Hann窗、物理中心频率±5%与±520µm孔径，频带能量除以同一T。BF采用相同窗口内10λ电荷加权局部RMS。右图使用原生powerSI每个波荡器周期的记录，同一Hann选窗规则，但包含全部横向网格与全部保留频率；因此10.4056MW原生峰值与左图6.66285MW不是同一个观测量，没有做幅度匹配。', '',
        '原生细曲线按一波长slice电荷插值中位时刻，出口曲线用实际全部粒子中位时刻；两种中位时刻在11个出口最大相差约0.502fs，远小于993.253fs窗宽。原生峰只在每周期记录中定位，间距约2.8cm；精确的连续±5%孔径谱需要额外中间三维场，而本次只有出口场。', '',
        '峰值候选的操作性规则是原生每周期功率的局部最大，prominence至少为整条曲线最大值的10%、候选间隔至少0.5m，以排除磁间隙处的小回落。本工况只筛出19.358m这一突出峰；原始完整曲线同时保存。若采用达到所采样带限峰值90%的另一个标记，第一次出口采样是第7段，交点只能限定在第6和第7段出口之间；这不是通用的物理饱和定义。', '',
        '![局部聚束因子](data/bunching_growth.png)', '',
        'BF第5段最高，其后逐步下降；BF峰与功率峰是不同的观测标记，不能用BF先回落就断言第5段功率已饱和。四倍滑移的抛物线束团也不保证时间功率轮廓存在严格平顶。本次是w≈355独立随机样本，未作为one-electron SASE饱和位置或原论文精度验收。', '',
        '## 完整性', '',
        '19799正常COMPLETED0:0，79200步、19:51:25，沿用同一已验程序和输入，无物理重跑。991个规则周期状态全部存在，另保存特殊出口状态；取每个物理步最早dump，避免最终correctTrans后的重复记录。11组出口均含3467224个宏粒子与完整电荷，窗口无缺帧；固定分母、共同孔径、BF范围、CSV及Parseval检查通过。', '',
        f"最大Parseval相对差{r['max_Parseval_relative']:.3g}，{r['closed_original_outputs']}个已闭合原文件inode/size/mtime保持，原生产冻结清单SHA前后通过。最初后处理在第8段因浮点累加位置差1.294e-11m触发过严的1e-11m检查；保留saturation-20261116-r1失败证据。有效r2仍要求物理步号精确一致，只把坐标舍入检查改为1nm，11面最大实际偏差3.792e-11m；没有改变原数据或统计公式。", '',
        '- [饱和位置与定义](SATURATION.json)', '- [逐段频带功率与BF](data/gain.csv)',
        '- [每周期原生完整历史](data/native_history.csv)', '- [峰值筛选规则和候选](data/NATIVE_PEAK_CANDIDATES.json)',
        '- [完整结果](data/RESULT.json)', '- [独立核验](VALIDATION.json)']
    (out/'README.md').write_text('\n'.join(lines)+'\n')
    validation=dict(passed=True,exit_rows=11,regular_period_states=991,unique_native_states=len(dense),
        all_exit_populations_and_charge_complete=True,common_power_denominator_s=r['fixed_window_duration_s'],
        maximum_Parseval_relative=r['max_Parseval_relative'],closed_output_identities_unchanged=len(protected),
        CSV_roundtrip_checked=True,BF_inequality_checked=True,only_closed_outputs_read=True,
        native_and_band_observables_separate=True,PNG_visual_review_pending=True,PDF_separately_rendered=False,
        hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in data.iterdir() if p.suffix in ['.json','.csv','.png','.pdf']})
    (out/'VALIDATION.json').write_text(json.dumps(validation,indent=2)+'\n')
    print(json.dumps(summary),flush=True)


if __name__=='__main__':main(sys.argv[1])
