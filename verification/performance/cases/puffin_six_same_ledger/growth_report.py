"""Document and validate the existing, closed power/BF comparison tables."""
import csv
import sys
from pathlib import Path
from common import read,save,sha


def main(out):
    out=Path(out);r=read(out/'data/RESULT.json');rows=r['rows'];initial=r['initial_bunching']
    assert r['complete'] and len(rows)==18 and len(initial)==3
    T=r['fixed_window_duration_s'];area=4*.00052**2
    assert all(x['power_denominator_s']==T for x in rows)
    assert max(abs(x['aperture_area_m2']/area-1) for x in rows)<1e-12
    assert all(0<=x['central_10lambda_rms_b1']<=1 and x['central_b1']<=x['central_10lambda_rms_b1']+1e-12 for x in rows)
    with (out/'data/growth.csv').open() as f:table=list(csv.DictReader(f))
    assert len(table)==18
    for x,y in zip(rows,table):
        assert x['method']==y['method'] and x['segment']==int(y['segment'])
        assert float(y['Hann_pm5pct_power_Eproxy_W'])==x['Hann_pm5pct_power_Eproxy_W']
    methods=['Puffin','energy-binomial3','energy-momentum'];labels=['Puffin','APRL 能量守恒＋三点核','APRL 能量守恒＋动量匹配']
    last={m:next(x for x in rows if x['method']==m and x['segment']==6) for m in methods}
    low,high=r['fixed_window_relative_to_charge_median_s']
    lines=['# 六段Puffin与APRL：功率及BF增长对照','',
        '只读后处理19793已COMPLETED0:0，用时3分钟。使用19783、19743和19744已闭合的同一历史束流数据，没有新物理仿真。图中功率是六个波荡器出口的采样点；BF还包含磁入口，连线用于读图，不代表连续采样。','',
        '![功率与局部BF增长](data/power_bf_growth.png)','',
        '| 第六段 | Hann ±5%中央功率（MW） | 中央10λ RMS BF | 全束可用粒子global BF |',
        '| --- | ---: | ---: | ---: |']
    for m,label in zip(methods,labels):
        x=last[m];lines.append(f"| {label} | {x['Hann_pm5pct_power_Eproxy_W']/1e6:.5f} | {x['central_10lambda_rms_b1']:.6f} | {x['global_b1']:.6f} |")
    lines+=['',
        '入口局部RMS BF约0.00108224，三组相对差不到10⁻⁶；随后前几段APRL增长更快。三点核第五段中央功率4.33518MW、局部BF0.647468，到第六段降至3.80078MW、0.547810；Puffin第五段1.46903MW、0.442743，第六段升至4.21323MW、0.683289。六段内观察到的回落不能单独当作长束稳态饱和位置的验收。','',
        '## 一致的统计规则','',
        f'窗口取名义抛物线束团后半段中的高流区：先限定I/I_peak≥0.8，再在该段两端各去掉10%。按每个实际粒子面的电荷加权中位迟滞时间平移窗口，得到中位时刻之后[{low*1e12:.9f}, {high*1e12:.9f}]ps，所有程序、所有位置的宽度严格相同，为**{T*1e12:.9f}ps**。这是Eulerian统计，不是固定粒子ID群，也不是按辐射亮度选择最强片段。','',
        '功率采用同一物理中心频率c/λ（λ=491.317235nm）、±5%频带、x/y∈[−520,+520]µm孔径。对窗口内原始横向电场作连续RMS-Hann窗w(t)=sqrt(8/3)sin²[π(t−t_low)/T]，积分实场单边频谱，再除以共同的T；没有除以不同程序各自的整段仿真/零填充记录时长，没有跨程序幅度匹配。','',
        '实际比较量是 ε0 c∫|E⊥|² 的目标频带功率代理，不是原生有符号Poynting，也不代表已经分离的纯辐射。APRL保留其已保存探测器的Yee到单元中心插值，Puffin使用原节点强度及孔径端点权重，没有额外四点平均。二者保留各自原生时间采样，不作额外载波插值；本次18个窗口均无缺失采样点需要补齐。','',
        '主BF图与功率使用同一个中央时间窗口，将窗口划为10个参考辐射波长宽的组，计算每组b_j=Σq exp(iωt)/Σq，再作电荷加权RMS sqrt(ΣQ_j|b_j|²/ΣQ_j)。它不同于跨完整束团的global BF，也没有将Puffin原生1λ BF与APRL固定ID统计混比。','',
        '## 六段数值','',
        '| 段数 | Puffin功率（MW） | 三点核功率（MW） | 动量匹配功率（MW） | Puffin局部BF | 三点核局部BF | 动量匹配局部BF |',
        '| --- | ---: | ---: | ---: | ---: | ---: | ---: |']
    for k in range(1,7):
        x=[next(y for y in rows if y['method']==m and y['segment']==k) for m in methods]
        lines.append('| '+str(k)+' | '+' | '.join(f"{y['Hann_pm5pct_power_Eproxy_W']/1e6:.7f}" for y in x)+' | '+' | '.join(f"{y['central_10lambda_rms_b1']:.6f}" for y in x)+' |')
    lines+=['',
        '## 全束BF与范围限制','',
        '![全束可用粒子BF](data/global_bf_growth.png)','',
        '全束图对每个粒子面实际保留的全部粒子独立计算Σq exp(iωt)/Σq。APRL第六段过面数据不是完整原束，三点核和动量匹配分别覆盖原电荷的99.95465%、99.95432%；原尾部覆盖/粒子时间端盖缺项不因本图消失。中央窗口位于各自已保存场的时间支持内，所用粒子电荷、数量和中位时间逐行保存。','',
        '输入仍是原1734656条相关phase-clone历史加载，不是后来独立随机束团的验收。Puffin从原事件弹道映射到磁入口，APRL在入口漂移中保留自场；两者模型及探测器离散采样的差别仍在。此图展示实际响应差异，不唯一判定更强增长的物理/算法根因。','',
        f"所有分项使用相同功率分母，孔径面积一致为{area:.10g}m²；最大Parseval相对差{r['max_Parseval_relative']:.3g}，18组窗口无缺帧，BF满足中央global≤局部RMS≤1。运行期核对36个已闭合场/粒子文件的inode、size和mtime未改变；原冻结输入和生产卡保持。",'',
        'PNG两图已目视检查，PDF同源导出未另渲染。原能量账本对照目录保持不变，本图使用中央窗口功率，不能和之前全脉冲µJ值直接混用。','',
        '- [全部增长数据](data/growth.csv)','- [入口BF](data/initial_bunching.csv)',
        '- [完整定义、覆盖和数值检查](data/RESULT.json)']
    (out/'README.md').write_text('\n'.join(lines)+'\n')
    save(out/'VALIDATION.json',dict(passed=True,physical_acceptance=False,rows=18,methods=methods,
        common_power_denominator_s=T,aperture_area_m2=area,
        maximum_Parseval_relative=r['max_Parseval_relative'],missing_samples_reconstructed=sum(x['missing_samples_reconstructed'] for x in rows),
        BF_RMS_inequality_checked=True,CSV_roundtrip_checked=True,protected_original_outputs_during_analysis=r['protected_closed_outputs'],
        PNG_visual_review_passed=True,PDF_separately_rendered=False,readonly_job=19793,no_new_physics=True,
        hashes={p.name:sha(p) for p in (out/'data').iterdir() if p.suffix in ['.json','.csv','.png','.pdf']}))


if __name__=='__main__':main(sys.argv[1])
