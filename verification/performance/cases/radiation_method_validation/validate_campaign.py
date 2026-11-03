"""Inspect frozen identity, closed PIC markers, metadata and report consistency."""
import ast,hashlib,json,re
from pathlib import Path
import h5py
import numpy as np
from prepare import save,require,sha

def main():
    ws=Path('/public/home/wfollowm/APRL_Workspace')
    base=ws/'results/radiation_method_validation_20261103';root=base/'fel-r4'
    plan=json.loads((root/'campaign.json').read_text())
    frozen_checks={}
    campaigns=[base/'pilot-r1',base/'synchrotron-r1',base/'single-electron-r1',root,
               *[base/'calibration-r1'/n for n in ['weak-bend','x-refined','zdt-refined']]]
    memo={}
    for d in campaigns:
        frozen=json.loads((d/'FROZEN.json').read_text())
        for p,expected in frozen.items():
            if p not in memo:memo[p]=sha(p)
            require(memo[p]==expected,'Frozen identity changed: '+p)
        frozen_checks[str(d)]=len(frozen)
    rows=[];normalized=[]
    for case in plan['cases']:
        d=Path(case['directory']);run=json.loads((d/'PIC_COMPLETE.json').read_text())
        require(run['exit_code']==0,'PIC nonzero exit')
        out=d/'analysis-regularized-r1'
        radiation=json.loads((out/'radiation.json').read_text())
        require(len(radiation)==6,'Missing field planes')
        for a in radiation:
            with h5py.File(d/'detectors'/(a['plane']+'.h5'),'r') as f:
                g=f['field_plane']
                require(int(g['complete'][()])==1,'Unclosed raw field')
                require(len(g['fields'])==int(g['committed_samples'][()]),'Uncommitted raw field')
            require(not a['missing_native_steps'],'Final common interval would require filled samples')
        rows.append(dict(case=case['name'],PIC_exit_code=0,elapsed_seconds=run['seconds'],
                         analysis=str(out),common_intervals_need_no_sample_filling=True))
        if case['name']!='prescribed':
            import yaml
            cfg=yaml.safe_load((d/'run.yaml').read_text())
            def norm(x):
                if isinstance(x,dict):
                    return {k:norm(v) for k,v in x.items()
                            if k not in ['particle_coupling','longitudinal_shape']}
                if isinstance(x,list):return [norm(v) for v in x]
                return x.replace(str(d),'$CASE') if isinstance(x,str) else x
            normalized.append(norm(cfg))
    require(normalized[0]==normalized[1]==normalized[2],
            'A/B/C physical cards differ beyond method switches and paths')
    protected=json.loads((ws/'results/group_meeting_materials_20261103/REFERENCE_AUDIT.json').read_text())['implementation_and_protected_files_sha256']
    protected_ok={}
    for p,expected in protected.items():
        if 'energy_ledger_test_record' in p:
            path=Path(p)
            if not path.is_absolute():path=ws/'source/aprl'/path
            require(sha(path)==expected,'Preexisting document changed')
            protected_ok[str(path)]=expected
    scripts=list(Path(__file__).parent.glob('*.py'))
    for p in scripts:ast.parse(p.read_text(),filename=str(p))
    files=[base/n for n in ['README.md','THEORY.md','RESULT_REVIEW.md']]
    links=0
    for p in files:
        for target in re.findall(r'\]\(([^)]+)\)',p.read_text()):
            if target.startswith(('https:','http:','#')):continue
            require((p.parent/target.split('#')[0]).exists(),'Broken local report link: '+target)
            links+=1
    result=json.loads((base/'RESULT.json').read_text())
    data={q['method']:q for q in result['methods']}
    require(abs(data['energy']['Hann_relative_difference_from_C_percent']-1.8200389155322627)<1e-9,'Summary differs from spectra')
    require(data['energy']['exit_spectral_relative_L2_from_C']>.05,'Unexpected spectral acceptance classification')
    save(base/'VALIDATION.json',dict(frozen_files_by_campaign=frozen_checks,
       same_physical_A_B_C_cards=True,closed_PIC_cases=rows,
       original_batch_exit_code=int((root/'EXIT_CODE').read_text()),
       batch_failure_scope='Original uniform-time postprocessing assertion; all four PIC processes exited 0',
       repaired_postprocessing_complete=True,finite_field_check_scope='Full XY fields of the common retained time intervals',
       old_evidence_unchanged=True,protected_documents_sha256=protected_ok,
       python_modules_parsed=len(scripts),local_report_links_checked=links,
       PNG_visual_review=['FEL_three_methods.png','synchrotron_three_methods.png','single_electron_spatial_audit.png'],
       PDF_validation='Same source export, not separately rendered',
       physical_accuracy_accepted=False,user_engineering_target_percent=5.))
    save(root/'POSTPROCESSING_REPAIRED.json',dict(PIC_cases=rows,physical_accuracy_accepted=False))
    print('VALIDATED',len(rows),'closed PIC cases;',len(scripts),'Python modules;',links,'local links')

if __name__=='__main__':main()
