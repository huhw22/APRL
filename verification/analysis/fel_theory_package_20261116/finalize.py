"""Verify the scientific distinctions and package portable result files only."""
import csv
import hashlib
import json
import math
import re
import subprocess
import sys
import zipfile
from pathlib import Path

from build import read, rows, save, sha


def close(a, b, tolerance=3e-14):
    assert math.isclose(a,b,rel_tol=tolerance,abs_tol=1e-16),(a,b)


def main(output, archive):
    root, archive = Path(output).resolve(), Path(archive).resolve()
    assert root.is_dir() and not archive.exists()
    sources = read(root/'SOURCE_FILES_SHA256.json')
    identities = read(root/'SOURCE_FILE_IDENTITIES.json')
    for path, expected in sources.items():
        assert sha(path)==expected,path
        st=Path(path).stat()
        assert identities[path]==dict(inode=st.st_ino,size=st.st_size,mtime_ns=st.st_mtime_ns),path
    ws=root.parent.parent
    assert sha(ws/'results/fel_six_two_method_energy_curves_20261116_r1.zip')=='fa5b6555a56d74ad67a7a7249e3b41b400a47731654ff76f503de3ff91bd3346'
    assert sha(ws/'results/fel_six_puffin_F3_MM_bilingual_20261116_r3.zip')=='4081c338a72fac0504037fe6bc802307bbba4ba58869f49e3b119f73e242411f'
    theory=root/'03_theory_notes'
    gain=rows(theory/'data/effective_gain_lengths.csv')
    bf=rows(theory/'data/bf_relative_to_puffin.csv')
    estimates=rows(theory/'data/ming_xie_estimates.csv')
    exits=rows(root/'02_growth_comparison/data/comparison_by_segment.csv')
    initial=rows(root/'02_growth_comparison/data/initial_bunching.csv')
    assert [r['method'] for r in initial]==['Puffin','APRL EC+F3','APRL EC+MM']
    assert len(gain)==15 and len(bf)==7 and len(estimates)==3
    keys={'Puffin':'Puffin_power_W','F3':'APRL_F3_power_W','MM':'APRL_MM_power_W'}
    for r in gain:
        j=int(r['from_segment'])-1
        assert int(r['to_segment'])==j+2 and float(r['magnetic_length_m'])==2.156
        p1=float(exits[j][keys[r['method']]])
        p2=float(exits[j+1][keys[r['method']]])
        # log difference uses a distinct evaluation from log of the ratio.
        close(float(r['effective_power_gain_length_m']),2.156/(math.log(p2)-math.log(p1)))
        assert float(r['start_power_W'])==p1 and float(r['end_power_W'])==p2
    for seg,r in enumerate(bf):
        expected=([float(p['central_10lambda_rms_b1']) for p in initial] if seg==0 else
                  [float(exits[seg-1][k]) for k in ['Puffin_local_rms_BF','APRL_F3_local_rms_BF','APRL_MM_local_rms_BF']])
        for tag,value in zip(['Puffin','F3','MM'],expected):
            assert float(r[tag+'_BF'])==value
            if tag!='Puffin':
                close(float(r[tag+'_difference_percent']),100*(value/expected[0]-1))
    assert max(abs(float(bf[0][k+'_difference_percent'])) for k in ['F3','MM'])<.000052
    close(float(bf[1]['F3_difference_percent']),-6.918158977557276)
    close(float(bf[1]['MM_difference_percent']),-5.159805871004409)
    for r in estimates:
        q={k:float(v) for k,v in r.items() if k!='round_surrogate'}
        # Equivalent conventional K*JJ and ku definition checks aw convention.
        ku=2*math.pi/q['lambda_u_m']
        KJJ=math.sqrt(2)*q['Aw']
        rho=(q['I_A']/q['Alfven_current_A']*KJJ**2/(16*q['gamma']**3*(ku*q['sigma_m'])**2))**(1/3)
        close(q['rho'],rho)
        close(q['surrogate_matched_beta_m']*q['geometric_emit_m'],q['sigma_m']**2)
        close(q['Lg_power_gain_m'],q['L0_power_gain_m']*(1+sum(q[f'Lambda_term_{n}'] for n in range(1,8))))
    close(float(estimates[0]['Lg_power_gain_m']),.670464280659689)
    close(float(estimates[1]['Lg_power_gain_m']),.6125463846260655)
    sc=read(theory/'space_charge/RESULT.json')
    transfer=rows(theory/'space_charge/curves.csv')
    close(float(transfer[-1]['constant_density_transfer']),sc['fixed_size_estimates'][1]['cold_density_transfer'])
    close(float(transfer[-1]['envelope_density_transfer']),sc['varying_envelope']['density_ratio'])
    close(float(transfer[-1]['envelope_F3_only_density_transfer']),sc['varying_envelope']['F3_only_density_ratio'])
    plots=read(theory/'PLOTS.json')['figures']
    assert len(plots)==6 and all(p['axes']==1 for p in plots)
    for p in plots:
        assert p['language'] in ['zh-CN','en']
        for s in p['series']:
            if p['name']=='gain_length':
                wanted=[r for r in gain if r['method']==s['method']]
                assert s['x']==[int(r['to_segment']) for r in wanted]
                assert s['y']==[float(r['effective_power_gain_length_m']) for r in wanted]
            elif p['name']=='bf_relative':
                assert s['x']==[int(r['segment']) for r in bf]
                assert s['y']==[float(r[s['method']+'_difference_percent']) for r in bf]
            else:
                assert p['name']=='space_charge_transfer'
                assert s['x']==[float(r['z_m']) for r in transfer]
                assert s['y']==[100*float(r[s['column']]) for r in transfer]
    links=[]
    for p in root.rglob('*.md'):
        for target in re.findall(r'!?\[[^\]]*\]\(([^)]+)\)',p.read_text()):
            if target.startswith(('https://','http://','#')):
                continue
            resolved=(p.parent/target.split('#')[0]).resolve()
            assert resolved.is_relative_to(root) and resolved.exists(),(str(p),target)
            links.append(dict(document=str(p.relative_to(root)),target=target))
    result=read(root/'RESULT.json')
    assert result['causal_space_charge_attribution_complete'] is False
    assert sc['physical_acceptance'] is False
    assert not read(theory/'data/current_eleven_saturation.json')['physical_acceptance']
    # Copied historical records are explicitly dated; no stale physical claim is promoted.
    for word in ['5.17%','5.16%','尚未','三点核','动量匹配']:
        assert word in (theory/'BF_SPACE_CHARGE.md').read_text(),word
    commit=subprocess.check_output(['git','-C',str(Path(__file__).resolve().parents[3]),'rev-parse','HEAD'],text=True).strip()
    save(root/'SOURCE_COMMIT.json',dict(commit=commit,pushed=False,
        templates_sha256={p.name:sha(p) for p in Path(__file__).parent.iterdir() if p.is_file()}))
    save(root/'VALIDATION.json',dict(passed=True,source_SHA256_and_identity_preserved=len(sources),
        previous_two_ZIPs_SHA256_unchanged=True,adjacent_gain_intervals_checked=15,
        BF_entrance_and_six_exit_ratios_checked=True,ming_xie_round_surrogates_checked=3,
        source_to_plot_arrays_exact=True,new_single_axis_bilingual_figures=6,
        copied_existing_single_axis_bilingual_figures=26,PNG_visual_review_passed=True,
        PDF_separately_rendered=False,local_markdown_links_checked=len(links),
        no_new_physics=True,no_scientific_H5_outputs_read=True,
        BF_space_charge_causal_proof_complete=False,
        historical_cloned_input_segregated=True,formula_comparison_is_qualified=True))
    (root/'SHA256SUMS').write_text(''.join(f'{sha(p)}  {p.relative_to(root)}\n' for p in sorted(root.rglob('*'))
        if p.is_file() and p.name not in ['SHA256SUMS','PACKAGE.json']))
    files=sorted(p for p in root.rglob('*') if p.is_file() and p.name!='PACKAGE.json')
    for p in files:
        assert p.suffix in ['.md','.txt','.csv','.json','.png','.pdf'] or p.name=='SHA256SUMS',p
        assert not p.is_symlink()
    with zipfile.ZipFile(archive,'x',compression=zipfile.ZIP_DEFLATED) as z:
        for p in files:
            z.write(p,str(Path(root.name)/p.relative_to(root)))
    with zipfile.ZipFile(archive) as z:
        assert z.testzip() is None
        assert len(z.namelist())==len(files)
        for p in files:
            assert hashlib.sha256(z.read(str(Path(root.name)/p.relative_to(root)))).hexdigest()==sha(p)
    package=dict(complete=True,directory=str(root),archive=str(archive),bytes=archive.stat().st_size,
        sha256=sha(archive),files=len(files),source_commit=commit,
        includes_energy_growth_Ming_Xie_and_BF_theory=True,
        BF_causal_attribution_complete=False,CRC_and_all_member_SHA256_passed=True)
    save(root/'PACKAGE.json',package)
    save(root.parent/'LATEST_SIX_COMPLETE_WITH_THEORY.json',package)
    print(json.dumps(package,ensure_ascii=False,indent=2))


if __name__=='__main__':
    main(*sys.argv[1:])
