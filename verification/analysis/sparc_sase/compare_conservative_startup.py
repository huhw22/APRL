"""Compare only completed energy and momentum startup analyses, without explicit PIC."""
import argparse
import csv
import json
from pathlib import Path

import numpy as np
import yaml

from entrance_ballistic import save, table, digest
from energy_startup_origin import run


def read_csv(path):
    rows=list(csv.DictReader(Path(path).open()))
    for row in rows:
        for k,v in row.items():
            try: row[k]=float(v)
            except ValueError: pass
    return rows


def card_audit(top):
    root=Path(top)/'compact-production-r1'
    cards={case:yaml.safe_load((root/case/'run.yaml').read_text()) for case in ['energy','momentum']}
    def changed(a,b,path=''):
        if isinstance(a,dict) and isinstance(b,dict):
            result=[]
            for k in sorted(set(a)|set(b)):
                result+=changed(a.get(k),b.get(k),path+'.'+k if path else k)
            return result
        return [] if a==b else [dict(path=path,energy=a,momentum=b)]
    diffs=changed(cards['energy'],cards['momentum'])
    allowed={'output.manifest','mesh.longitudinal_shape','mesh.moving_window.reference_directory',
             'detectors.directory','energy_ledger.directory'}
    assert all(r['path'] in allowed for r in diffs),diffs
    assert cards['energy']['mesh']['particle_coupling']==cards['momentum']['mesh']['particle_coupling']=='conservative-source'
    assert cards['energy']['mesh']['longitudinal_shape']=='linear'
    assert cards['momentum']['mesh']['longitudinal_shape']=='quintic-compensated'
    return dict(differences=diffs,physical_card_differences=['mesh.longitudinal_shape'],
        same_input=True,same_fields_and_detectors=True,
        interpretation='Tests the complete longitudinal-shape coupling change; a difference alone cannot isolate initial reference fields, force gathering, or deposition.')


def compare(top, energy, momentum, out):
    top,energy,momentum,out=map(lambda p:Path(p).resolve(),[top,energy,momentum,out])
    out.mkdir(exist_ok=False)
    dirs={'energy':energy,'momentum':momentum}
    results={case:json.loads((p/'RESULT.json').read_text()) for case,p in dirs.items()}
    for case,p in dirs.items():
        complete=json.loads((p/'COMPLETE.json').read_text())
        assert complete['complete'] and complete['case']==case
        assert results[case]['case']==case and results[case]['no_new_simulation']
    assert results['energy']['all_input_particles']==results['momentum']['all_input_particles']
    assert results['energy']['wavelength_m']==results['momentum']['wavelength_m']
    # Identity checks include the exact frozen particle source hash, not just count.
    manifests={case:json.loads((p/'VALIDATION.json').read_text())['source_manifest'] for case,p in dirs.items()}
    input_sha={case:next(r['sha256'] for r in rows if r['path'].endswith('compact-preflight-r1/particles.h5')) for case,rows in manifests.items()}
    assert len(set(input_sha.values()))==1
    audit=card_audit(top)
    save(out/'CARD_AUDIT.json',audit)
    stages={case:read_csv(p/'stage_bunching.csv') for case,p in dirs.items()}
    harmonic={case:read_csv(p/'harmonics.csv') for case,p in dirs.items()}
    phases={case:read_csv(p/'phase_bins.csv') for case,p in dirs.items()}
    rows=[]
    for case,r in results.items():
        early=json.loads((dirs[case]/'und-entry_stage.json').read_text())
        end=r['p1']
        h=next(x for x in harmonic[case] if x['plane']=='p1' and x['harmonic']==1 and x['quantity']=='phase_residual_within')
        rows.append(dict(case=case,input_particles=r['all_input_particles'],
            entrance_local_b=early['bunching']['actual']['local_rms'],
            p1_global_b=end['bunching']['actual']['global_b'],
            p1_local_b=end['bunching']['actual']['local_rms'],
            p1_phase_alignment=end['bunching']['actual']['phase_alignment'],
            p1_family_only_local_b=end['bunching']['family_common_only']['local_rms'],
            p1_within_only_local_b=end['bunching']['within_family_only']['local_rms'],
            p1_phase_rms_rad=end['residual_phase_rms_rad'],
            p1_within_phase_rms_rad=end['within_phase_rms_rad'],
            p1_phase_first_harmonic_rad=h['global_abs'],
            raw_scalar_gap_und_eV=r['detector_energy_comparison']['rows'][1]['raw_minus_scalar_eV'],
            raw_mean_change_und_eV=r['detector_energy_comparison']['raw_mean_change_eV'],
            scalar_mean_change_und_eV=r['detector_energy_comparison']['scalar_mean_change_eV'],
            endpoint_phase_energy_correlation=r['correlations']['within_phase_vs_energy']))
    table(out/'comparison.csv',rows)
    a,b=rows
    change={key:dict(energy=a[key],momentum=b[key],
        relative_change=float(b[key]/a[key]-1) if a[key] else None)
        for key in ['p1_global_b','p1_local_b','p1_within_phase_rms_rad',
                    'p1_phase_first_harmonic_rad','raw_scalar_gap_und_eV']}
    save(out/'RESULT.json',dict(rows=rows,changes=change,input_sha256=input_sha['energy'],
        case_comparison=['energy','momentum'],explicit_case_read=False,new_simulation=False,
        physical_acceptance=False,
        interpretation='Paired identical-input startup comparison. Endpoint raw energy and time interpolation remain separate limitations; longitudinal-shape sensitivity is not alone a force-source proof.'))
    render(out,rows,stages,phases,harmonic)
    document(out,rows,change)
    save(out/'VALIDATION.json',dict(only_completed_cases=True,case_markers_checked=True,
        input_sha256_equal=True,card_audit_pass=True,originals=[digest(p/'RESULT.json') for p in dirs.values()]))
    save(out/'COMPLETE.json',dict(complete=True,new_simulation=False))


def render(out,rows,stages,phases,harmonic):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from scipy.constants import m_e,c,e
    plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(2,2,figsize=(12.6,8.3),layout='constrained')
    style={'energy':('Energy coupled','#24648c','o-'),'momentum':('Momentum matched','#8e4b75','s--')}
    for case,(label,color,line) in style.items():
        ss=[r for r in stages[case] if r['phase_model']=='actual']
        axes[0,0].semilogy([r['z_m'] for r in ss],[r['local_rms'] for r in ss],line,
            color=color,label=label)
        ss=[r for r in stages[case] if r['plane']=='p1']
        axes[0,1].semilogy(range(4),[r['local_rms'] for r in ss],line,color=color,label=label)
        pp=phases[case]
        axes[1,0].plot([r['phase_rad'] for r in pp],[r['phase_within_mean_rad'] for r in pp],line,
            color=color,ms=3,label=label)
        hh=[r for r in harmonic[case] if r['plane']=='p1' and r['quantity']=='delta_gamma_within']
        axes[1,1].semilogy([r['harmonic'] for r in hh],
            [r['global_abs']*m_e*c*c/e for r in hh],line,color=color,ms=3,label=label)
    axes[0,0].set_xlabel('Fixed laboratory plane z (m)');axes[0,0].set_ylabel('Local RMS BF')
    axes[0,0].set_title('(a) Stored particle planes; no intermediate onset inferred')
    axes[0,1].set_xticks(range(4),['Reference','Family common','Within family','Actual PIC'],rotation=15)
    axes[0,1].set_ylabel('Local RMS BF');axes[0,1].set_title('(b) Phase rearrangement at first segment exit')
    axes[1,0].set_xlabel('Undulator-entrance optical phase (rad)')
    axes[1,0].set_ylabel('Mean within-family phase change (rad)')
    axes[1,0].set_title('(c) Optical-period phase waveform')
    axes[1,1].set_xlabel('Harmonic of entrance optical phase')
    axes[1,1].set_ylabel('Raw-plane energy harmonic (eV)')
    axes[1,1].set_title('(d) Endpoint energy spectrum; diagnostic limitations apply')
    for ax in axes.flat:ax.grid(alpha=.2);ax.legend(fontsize=9)
    fig.suptitle('Same-input startup: energy coupling and momentum matching',fontsize=14)
    fig.savefig(out/'conservative_startup_comparison.png',dpi=220,bbox_inches='tight');plt.close(fig)


def document(out,rows,change):
    a,b=rows
    text='# 同输入首段启动：能量守恒与动量匹配对照\n\n'
    text+='仅在动量版正式计算正常结束、所有读取面闭合后生成；不读取显式版。\n\n'
    text+=f"首段局部RMS BF：能量版{a['p1_local_b']:.9g}，动量版{b['p1_local_b']:.9g}，相对变化{100*change['p1_local_b']['relative_change']:+.3f}%；全局BF相对变化{100*change['p1_global_b']['relative_change']:+.3f}%。\n\n"
    text+='![对照图](conservative_startup_comparison.png)\n\n'
    text+='|口径|能量版|动量版|\n|---|---:|---:|\n'
    for key,label in [('p1_global_b','首段全局BF'),('p1_local_b','首段局部RMS BF'),
        ('p1_phase_alignment','首段相位对齐'),('p1_family_only_local_b','仅组共同迟滞局部BF'),
        ('p1_within_only_local_b','仅组内差别局部BF'),('p1_phase_first_harmonic_rad','组内相位一阶谐波/rad'),
        ('raw_scalar_gap_und_eV','磁入口raw减独立标量能量/eV')]:
        text+=f'|{label}|{a[key]:.9g}|{b[key]:.9g}|\n'
    text+='\n两版粒子源SHA完全相同，物理卡唯一方法差异为纵向形状linear→quintic-compensated；粒子耦合均为conservative-source。初场供入、沉积与反馈的整体离散效果需要一并考虑，不能仅凭变化指派某一项力。\n\n'
    text+='如果动量版BF显著下降，说明早期增强对纵向离散配对敏感；如果BF仍相近，则改善lab能量—动量收支并未消除该启动机制，需要继续查共同场源与入场过程。两种结果均不能仅凭单种子宣布SASE物理验收。\n\n'
    text+='原粒子面整Maxwell段到达时间/三动量插值仍是共同限制，独立标量能量面用于单独核对同面能量重建差；推进收支闭合不自动保证探测插值精度。\n\n'
    text+='每版完整相位分解分别见原能量版报告和../momentum/README.md；数据comparison.csv、RESULT.json、CARD_AUDIT.json及VALIDATION.json。\n'
    (out/'README.md').write_text(text)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('top');p.add_argument('energy');p.add_argument('out')
    a=p.parse_args();root=Path(a.out).resolve();root.mkdir(exist_ok=False)
    # run() checks the closed case marker before opening any momentum HDF5 file.
    run(a.top,root/'momentum',case='momentum')
    compare(a.top,a.energy,root/'momentum',root/'comparison')
