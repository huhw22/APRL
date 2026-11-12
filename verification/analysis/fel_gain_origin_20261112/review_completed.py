"""Review closed 19751 results; do not inspect the active ledger campaigns."""
import csv
import hashlib
import json
import math
import os
import sys
from pathlib import Path

import h5py
import numpy as np

from prepare import save, sha


def table(path,rows):
    with Path(path).open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def main(root):
    root=Path(root).resolve();out=root/'analysis/reviewed'
    out.mkdir(exist_ok=False)
    p=json.loads((root/'PLAN.json').read_text())
    data=json.loads((root/'analysis/RESULT.json').read_text())
    gate=json.loads((root/'GENESIS_COMPLETE.json').read_text())
    assert gate['complete'] and data['complete']
    frozen=json.loads((root/'FROZEN.json').read_text())
    assert all(sha(path)==digest for path,digest in frozen.items())
    pic=data['pic_rows'];exits=data['genesis_exits']
    gen={r['module']:r for r in exits if r['case']=='sc-off-fixed'}
    input_n=pic[0]['represented_electrons']
    rows=[]
    for a in pic:
        module=a['module'];b=gen['entry' if module==0 else module]
        rows.append(dict(module=module,z_m=a['z_m'],
            PIC_global_b1=a['global_b1'],GENESIS_global_b1=b['global_b1'],
            PIC_local_rms_b1=a['source_group_local_rms_b1'],GENESIS_local_rms_b1=b['source_group_local_rms_b1'],
            PIC_over_GENESIS_local_b1=a['source_group_local_rms_b1']/b['source_group_local_rms_b1'],
            PIC_common_phase_BF_fraction=abs(a['global_b1']/a['source_group_local_rms_b1'])**2,
            GENESIS_common_phase_BF_fraction=abs(b['global_b1']/b['source_group_local_rms_b1'])**2,
            PIC_missing_input_charge_fraction=1-a['represented_electrons']/input_n))
    table(out/'bunching_comparison.csv',rows)
    adjacent=[]
    for j in range(1,6):
        x=pic[j+1]['source_group_local_rms_b1']/pic[j]['source_group_local_rms_b1']
        y=gen[j+1]['source_group_local_rms_b1']/gen[j]['source_group_local_rms_b1']
        adjacent.append(dict(from_module=j,to_module=j+1,PIC_BF_amplitude_ratio=x,
            GENESIS_BF_amplitude_ratio=y,PIC_log_BF_growth_per_m=math.log(x)/(77*.028),
            GENESIS_log_BF_growth_per_m=math.log(y)/(77*.028),
            definition='BF amplitude growth over cumulative magnetic length; not a radiation-power gain-length measurement.'))
    table(out/'adjacent_bunching.csv',adjacent)
    matched=data['same_band_and_aperture']
    values={(r['case'],r['module']):r['GENESIS_rect_5pct_aperture_J'] for r in matched}
    controls=[]
    for name,a,b in [('ADI repair','longrange-fixed','longrange-installed'),
                     ('Half longitudinal step','sc-off-halfstep','sc-off-fixed'),
                     ('Long-range SC on','longrange-fixed','sc-off-fixed')]:
        controls.append(dict(control=name,numerator=a,denominator=b,
            p1_band_energy_change_percent=100*(values[a,1]/values[b,1]-1),
            p6_band_energy_change_percent=100*(values[a,6]/values[b,6]-1)))
    table(out/'reference_controls.csv',controls)
    checks=[]
    gamma=p['original_plan']['physics']['gamma'];ks=2*np.pi/p['original_plan']['physics']['wavelength_m']
    for case in p['cases']:
        directory=Path(case['directory'])
        complete=json.loads((directory/'COMPLETE.json').read_text())
        assert complete['complete'] and complete['input_import_bitwise']
        assert complete['launch']['exit_code']==0
        assert '[binding] result=PASS failures=0' in (directory/'binding.log').read_text()
        with h5py.File(directory/'reference.out.h5','r') as f:
            z=f['Lattice/z'][:];dz=f['Lattice/dz'][:];aw=f['Lattice/aw'][:]
            ku=f['Lattice/ku'][:];phase=f['Lattice/phaseshift'][:]
            longitudinal=np.sqrt(gamma*gamma-1-aw*aw)
            rate=np.where(ku==0,ks/(2*gamma*gamma),ku)-ks*(1+aw*aw)/(longitudinal*(gamma+longitudinal))
            closure=[]
            for j in range(5):
                z0=.01+j*(77*.028+.4);z1=z0+77*.028+.4
                mask=(z>=z0-1e-8)&(z<z1-1e-8)
                closure.append(float(np.angle(np.exp(1j*(np.sum(rate[mask]*dz[mask])+np.sum(phase[mask]))))))
            assert max(abs(v) for v in closure)<1e-8,(case['name'],closure)
            # Native field-dump markers are at the requested physical exits.
            dumps=sorted(directory.glob('reference.*.fld.h5'),key=lambda f:int(f.name.split('.')[1]))
            positions=[float(z[int(path.name.split('.')[1])]) for path in dumps]
            assert np.allclose(positions,[pic[1]['z_m'],pic[6]['z_m']],rtol=0,atol=2e-10)
            I=f['Beam/current'][0]
            ds=float(f['Global/lambdaref'][0])*float(f['Global/sample'][0])
            power=f['Field/power']
            compare=np.array([power[k].sum()*ds/299792458. for k in [0,len(power)//2,len(power)-1]])
            assert np.allclose(compare,f['Field/Global/energy'][[0,len(power)//2,len(power)-1]],rtol=3e-12,atol=1e-24)
            represented=float(np.sum(I)*ds/(299792458*1.602176634e-19))
            assert abs(represented/input_n-1)<1e-12
        checks.append(dict(case=case['name'],roundtrip_all_active_cells_bitwise=True,
            gap_phase_closure_rad=closure,field_dump_positions_m=positions,
            input_represented_electrons=represented,native_power_energy_check=True))
    result=dict(complete_review=True,job_id=19751,reviewed_cases=checks,
        frozen_files_verified=len(frozen),bunching=rows,adjacent_growth=adjacent,
        same_band_and_aperture=matched,reference_controls=controls,
        reference_import_reinitialization_excluded_for_this_path=True,
        coherent_fraction_definition='|charge-weighted mean complex b|^2 / charge-weighted mean |b|^2 over fixed original 10-cell groups. A BF projection norm fraction, not a radiation-energy fraction.',
        conclusion='Early excess is established by the first two modules; after module 3 GENESIS adjacent local-BF growth is larger. ADI, step size and long-range SC controls do not explain the order-of-magnitude early radiation difference.',
        physical_root_cause_proven=False,new_simulations=False,
        current_ledger_campaigns_inspected=False,
        limitations=['One noise seed; fixed source-ID groups; PIC transverse charge losses explicitly reported.',
            'Current GENESIS short-range space charge is not called. The tested long-range switch cannot rule out 3D space-charge/loading effects.',
            'PIC band energy uses total native E/B flux; GENESIS uses a forward x radiation envelope. Matching aperture and band does not separate near fields.',
            'Original analysis chart mixes continuous full-envelope/full-aperture energy with matched endpoint markers; reviewed chart separates these quantities.'])
    save(out/'REVIEW.json',result)
    render(out,root,p,rows,adjacent,matched)
    write_report(out,result)
    print(json.dumps(dict(review_complete=True,frozen_files_verified=len(frozen),controls=controls,
        first_module_PIC_common_phase_BF_fraction=rows[1]['PIC_common_phase_BF_fraction'],
        first_module_GENESIS_common_phase_BF_fraction=rows[1]['GENESIS_common_phase_BF_fraction']),indent=2))


def render(out,root,p,rows,adjacent,matched):
    os.environ.setdefault('MPLCONFIGDIR',str(out/'mpl-cache'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False,'lines.linewidth':1.8})
    fig,axes=plt.subplots(2,2,figsize=(13.5,9),layout='constrained')
    blue='#337ba6';orange='#b8792e';black='#333333'
    styles={'sc-off-fixed':(blue,'-'),'longrange-fixed':(orange,'--'),
            'longrange-installed':('#a7537d',':'),'sc-off-halfstep':('#657a48','-.')}
    labels={'sc-off-fixed':'GENESIS: SC off, ADI fixed',
            'longrange-fixed':'GENESIS: long-range SC, ADI fixed',
            'longrange-installed':'GENESIS: installed ADI',
            'sc-off-halfstep':'GENESIS: SC off, half step'}
    ax=axes[0,0]
    for case in p['cases']:
        name=case['name']
        with np.load(root/'analysis'/(name+'-growth.npz')) as f:
            color,style=styles[name]
            ax.semilogy(f['z_m'],f['source_group_local_rms_b1'],style,color=color,label=labels[name])
    z=[r['z_m'] for r in rows]
    ax.semilogy(z,[r['PIC_local_rms_b1'] for r in rows],'o-',color=black,label='APRL: energy conserving + F3')
    ax.set_ylabel('Local RMS bunching factor');ax.set_xlabel('Laboratory z (m)')
    ax.set_title('(a) Same input; fixed original 10-cell groups');ax.legend(fontsize=8.5)
    ax=axes[0,1]
    ax.plot(z,[100*r['PIC_common_phase_BF_fraction'] for r in rows],'o-',color=black,label='APRL')
    ax.plot(z,[100*r['GENESIS_common_phase_BF_fraction'] for r in rows],'s--',color=blue,label='GENESIS: SC off, ADI fixed')
    ax.set_ylim(0,100);ax.set_xlabel('Laboratory z (m)')
    ax.set_ylabel('Common-phase fraction of BF norm (%)')
    ax.set_title('(b) Shared phase across original beam groups');ax.legend(fontsize=9)
    ax=axes[1,0]
    for case in p['cases']:
        name=case['name'];rr=[r for r in matched if r['case']==name]
        color,style=styles[name]
        ax.semilogy([0,1],[r['GENESIS_rect_5pct_aperture_J']*1e6 for r in rr],
                    marker='s',ls='None',mfc='none',color=color,markersize=9,label=labels[name])
    ref=[r for r in matched if r['case']=='sc-off-fixed']
    ax.semilogy([0,1],[r['PIC_rect_5pct_aperture_J']*1e6 for r in ref],'o',color=black,markersize=7,label='APRL')
    for j,r in enumerate(ref):
        ax.annotate(f"{r['PIC_over_GENESIS']:.2f} x",(j,r['PIC_rect_5pct_aperture_J']*1e6),
                    xytext=(11,0),textcoords='offset points',va='center',fontsize=10)
    ax.set_xticks([0,1],['First segment exit','Sixth segment exit']);ax.set_xlim(-.35,1.45)
    ax.set_ylabel('Band energy (µJ)');ax.set_title('(c) Rect ±5% band; same ±560 µm aperture')
    ax.legend(fontsize=8.3,loc='center left')
    ax=axes[1,1]
    x=np.arange(len(adjacent))
    ax.plot(x,[r['PIC_BF_amplitude_ratio'] for r in adjacent],'o-',color=black,label='APRL')
    ax.plot(x,[r['GENESIS_BF_amplitude_ratio'] for r in adjacent],'s--',color=blue,label='GENESIS: SC off, ADI fixed')
    ax.axhline(1,color='#777777',ls=':',lw=1)
    ax.set_xticks(x,[f"{r['from_module']} → {r['to_module']}" for r in adjacent])
    ax.set_ylim(0,5.6);ax.set_xlabel('Adjacent segment exits')
    ax.set_ylabel('Local BF amplitude multiplication')
    ax.set_title('(d) Growth per segment; excludes entrance interval');ax.legend(fontsize=9)
    for ax in axes.flat:ax.grid(alpha=.2)
    fig.suptitle('Same-input SPARC comparison: the largest excess forms during start-up',fontsize=15)
    fig.savefig(out/'comparison_reviewed.png',dpi=190);plt.close(fig)


def write_report(out,r):
    first=r['bunching'][1];last=r['bunching'][-1]
    lines=['# 19751：同输入增益对照的完成审查','',
        '**四组Genesis均正常完成；导入后的全部活跃单元六维数组及电流逐位保持。当前差异不来自这次导入重置或重新加噪。**',
        '', '![同输入对照](comparison_reviewed.png)','',
        'APRL取此前已闭合的十一段工况前六段，19743/19744账本任务不在本报告范围。本轮没有新物理仿真。',
        '', '## 直接可比的结果','',
        'BF按固定原始10个光学单元、先合并复数再取模统计；功率表使用完整记录Rect、共振±5%频带、相同±560µm横向孔径。没有幅度匹配或曲线归一化。',
        '', '|位置|APRL局部RMS BF|Genesis SC-off局部RMS BF|APRL/Genesis|',
        '|---|---:|---:|---:|']
    for x in r['bunching']:
        lines.append(f"|{x['module']}|{x['PIC_local_rms_b1']:.8g}|{x['GENESIS_local_rms_b1']:.8g}|{x['PIC_over_GENESIS_local_b1']:.4g}|")
    lines += ['', '|位置|APRL带内能量|Genesis SC-off带内能量|比值|',
              '|---|---:|---:|---:|']
    for x in r['same_band_and_aperture']:
        if x['case']=='sc-off-fixed':
            factor,unit=(1e12,'pJ') if x['module']==1 else (1e6,'µJ')
            lines.append(f"|第{x['module']}段|{x['PIC_rect_5pct_aperture_J']*factor:.6g} {unit}|{x['GENESIS_rect_5pct_aperture_J']*factor:.6g} {unit}|{x['PIC_over_GENESIS']:.5g}|")
    lines += ['', '## 三项Genesis控制的影响','',
        '|改变|第一段同带同孔径能量变化|第六段变化|','|---|---:|---:|']
    for x in r['reference_controls']:
        lines.append(f"|{x['control']}|{x['p1_band_energy_change_percent']:+.5f}%|{x['p6_band_energy_change_percent']:+.5f}%|")
    lines += ['',
        '这三项都不能解释第一段约16倍的带内差异。关闭/开启的是已接线的long-range模型；shortRange()仍未被推进器调用，因此本测试没有排除三维短程自场/加载关联的影响。',
        '', '## 增强发生在哪里','',
        f"入口两者局部BF均约0.00049859。第一段APRL为{first['PIC_local_rms_b1']:.7g}，Genesis为{first['GENESIS_local_rms_b1']:.7g}。",
        '更显著的特征是，APRL出现了跨束团各组共同相位的调制。定义',
        '', '$$C_b=|B_{global}|^2/b_{local,rms}^2.$$','',
        f"第一段APRL的C_b={100*first['PIC_common_phase_BF_fraction']:.3f}%，Genesis为{100*first['GENESIS_common_phase_BF_fraction']:.3f}%；入口两者均约0.748%。",
        '这是固定分组BF平方范数在共同相位分量上的投影比例，**不是辐射能量的相干份额**。',
        '', '第一至第二段的局部BF增幅：APRL 5.057倍、Genesis 2.505倍；第三至第四段开始，Genesis的相邻段BF增幅已经更大。第六段APRL接近饱和/回落而Genesis仍增长。因此目前证据支持“启动时出现了额外共同驱动/相干种子，随后放大”的定位；不能将全部差异描述成全程统一过大的指数增益率。',
        '', '此前16个相位克隆的横向关联、早期纵向场驱动，以及反射加载控制是与此相容的具体机制。当前全输入检查证实108416组均严格保留16重横向与动量状态；但本次没有对当前F3施加加载/自场干预，仍未证明最终力学根因或定量贡献。',
        '', '## 审查与限制','',
        f"核验{r['frozen_files_verified']}项冻结文件、四组32MPI绑定、输入回读、五个间隙的名义相位闭合、第一/六段dump位置及原生功率—脉冲能量关系。两端各8个原生光片的场能占比为0。",
        f"PIC至第六段缺失电荷占输入{100*last['PIC_missing_input_charge_fraction']:.6g}%，已逐面报告；没有通过共同ID筛选人为对齐入口噪声。",
        'APRL的目标带量来自总E/B原生过面能流，Genesis为前向x偏振辐射包络。同孔径同频带不等于完成近场分离；单个噪声种子也不提供统计精度验收。',
        '原analysis/same_input_comparison.png保留为自动初图，其中连续Genesis能量曲线为全孔径、全包络带宽，不能与APRL窄带曲线直接逐点作幅度判断；本审查图将严格匹配的两个出口另列。',
        '', '数据：[审查JSON](REVIEW.json)、[BF表](bunching_comparison.csv)、[相邻段增幅](adjacent_bunching.csv)、[参考控制](reference_controls.csv)。原始数据和初查报告全部保留。']
    (out/'README.md').write_text('\n'.join(lines)+'\n')


if __name__=='__main__':main(sys.argv[1])
