"""Freeze independent w355 six-segment paired and eleven-segment Puffin jobs."""
import ast
import copy
import json
import math
import shutil
import sys
from pathlib import Path

W=Path('/public/home/wfollowm/APRL_Workspace')
PUF=W/'results/puffin_eleven_four_slip_20261113_r2'
PIC=W/'results/sparc_six_full_ledger_20261112_r3'


def main(root):
    root=Path(root).resolve();root.mkdir(exist_ok=False)
    template=Path(__file__).resolve().parent
    shutil.copytree(template,root/'scripts',ignore=shutil.ignore_patterns('__pycache__'))
    for name in ['common.py','source_math.py','source_verify.py','puffin_units.py','binding.py']:
        shutil.copy2(PUF/'scripts'/name,root/'scripts'/name)
    sys.path.insert(0,str(root/'scripts'))
    from common import read,save,need,sha,C,E
    from source_math import nominal_electrons,validate_config
    prior=read(PUF/'PLAN.json');old=read(PIC/'PLAN.json')
    need(read(PUF/'FIELD_VALIDATED.json')['passed'] and read(PUF/'SMOKE_VALIDATED.json')['passed'],'Unqualified Puffin binary')
    need(read(PIC/'PREFLIGHT_COMPLETE.json')['complete'],'Unqualified ledger binary')
    need(sha(prior['binary'])==read(PUF/'FROZEN_PHYSICS.json')[prior['binary']],'Puffin binary changed')
    need(sha(old['binary'])==read(PIC/'PREFLIGHT_COMPLETE.json')['binary_sha256'],'APRL binary changed')
    # Reuse the established independent-source generator. Only the filename,
    # field padding and module/output parameters become campaign parameters.
    text=(PUF/'scripts/make_inputs.py').read_text()
    text=text[:text.index('    for directory,evolve in ')]
    text=text.replace("root/'inputs/independent-w64.h5'","root/'inputs/independent-w355.h5'")
    text=text.replace('span+1100+256',"span+plan['slippage_padding_wavelengths']+256")
    text=text.replace('nz=256*next_fast_len(math.ceil(needed/256))',"nz=plan['field_nodes_z2'];need(nz>=needed,'Original field extent does not contain the new source')")
    text+='''    from puffin_deck import deck
    for name in ['import-check','main']:
        deck(root,root/name,plan,scale,nz,preflight=name=='import-check')
    return info
'''
    (root/'scripts/make_inputs.py').write_text(text)
    # Keep the already-qualified six-segment deck writer, including magnetic
    # end ramps, the original continuous external route and redistribution.
    decks=template.parent/'puffin_six_same_ledger/prepare.py'
    original=decks.read_text();tree=ast.parse(original)
    function=next(x for x in tree.body if isinstance(x,ast.FunctionDef) and x.name=='deck')
    (root/'scripts/puffin_deck.py').write_text(ast.get_source_segment(original,function)+'\n')
    for name in ['six','eleven']:
        d=root/name;d.mkdir()
        for sub in ['inputs','import-check','main','preparation']:(d/sub).mkdir()
        (d/'scripts').symlink_to('../scripts',target_is_directory=True)
        plan=copy.deepcopy(prior);segments=6 if name=='six' else 11
        cfg=copy.deepcopy(prior['source_config'])
        duration=old['physics']['duration_s'] if segments==6 else prior['duration_s']
        cfg['profile'].update(full_duration_s=duration,window_start_s=0.,window_stop_s=duration)
        cfg['electrons']=round(nominal_electrons(cfg));cfg=validate_config(cfg)
        outputs=copy.deepcopy(prior['outputs'][:segments])
        periods=math.ceil((outputs[-1]['actual_lab_z_m']-.01)/.028) if segments==6 else prior['propagation_module_periods']
        nz=read(W/'results/puffin_six_same_ledger_20261113_r1/PLAN.json')['field_nodes_z2'] if segments==6 else read(PUF/'INPUT_VALIDATED.json')['field_nodes_z2']
        plan.update(root=str(d),source_config=cfg,target_weight=355,sample_records=round(cfg['electrons']/355),
            segments=segments,outputs=outputs,module_periods=periods,gamma_ref=cfg['phase_space']['gamma_mean'],
            duration_s=duration,source_charge_C=cfg['electrons']*E,source_length_m=C*duration,
            slippage_multiple=2. if segments==6 else 4.,
            field_nodes_z2=nz,slippage_padding_wavelengths=600 if segments==6 else 1100,
            physics_timeout_s=(22 if segments==6 else 70)*3600,
            main_nodes=1 if segments==6 else 2,main_ranks=128 if segments==6 else 256,
            independent_physical_source=True,APRL_input_saved=True,analysis_after_user_wakeup=True,
            original_binary_sha256=sha(prior['binary']),qualified_parent=str(PUF),
            common_drift_m=.01,magnetic_lattice=prior['lattice']['magnets'][:2*segments-1],
            magnetic_entrance='Exact original smooth one-period vector-potential ramp, including envelope derivative; route and real first-step checks required.',
            source_scope='Full original parabola; independent fixed-N six-dimensional source and uniform thinning without replacement; no transverse copies, phase fitting or extra shot noise.')
        save(d/'PLAN.json',plan)
    (root/'six-aprl').mkdir()
    # The production driver is the existing complete-ledger driver. The new
    # hard limit is a resource change; the numerical solver is reused verbatim.
    run=(PIC/'scripts/run_case.py').read_text()
    need(run.count("'170400'")==1,'Unexpected ledger driver timeout')
    (root/'scripts/run_pic.py').write_text(run.replace("'170400'","str(P['physics_timeout_s'])"))
    save(root/'PLAN.json',dict(root=str(root),target_weight=355,seed=prior['source_seed'],
        cancelled_job=19775,cancelled_parent=str(PUF),
        six=str(root/'six'),eleven=str(root/'eleven'),aprl=str(root/'six-aprl'),
        aprl_binary=old['binary'],aprl_binary_sha256=sha(old['binary']),
        puffin_binary=prior['binary'],puffin_binary_sha256=sha(prior['binary']),
        aprl_parent=str(PIC),puffin_parent=str(PUF),
        matched_six_input=True,aprl_methods=['energy-binomial3','energy-momentum'],
        no_core_changes=True,no_rebuild=True,no_default_change=True))
    files=[root/'PLAN.json',root/'six/PLAN.json',root/'eleven/PLAN.json',
        Path(old['binary']),Path(prior['binary']),
        PIC/'PREFLIGHT_COMPLETE.json',PIC/'PLAN.json',PIC/'energy-binomial3/run.yaml',
        PUF/'FIELD_VALIDATED.json',PUF/'SMOKE_VALIDATED.json',PUF/'SOURCE_METHOD_AUDIT.json',
        PUF/'source/puffin/lib/undulator/bfields.f90',
        *[p for p in (root/'scripts').iterdir() if p.is_file()]]
    save(root/'FROZEN_PREPARATION.json',{str(p):sha(p) for p in files})
    save(root/'SOFTWARE_REUSE.json',dict(aprl_validation=str(PIC/'PREFLIGHT_COMPLETE.json'),
        puffin_field_validation=str(PUF/'FIELD_VALIDATED.json'),
        puffin_stock_path_validation=str(PUF/'SMOKE_VALIDATED.json'),
        source_math_sha256=sha(root/'scripts/source_math.py'),
        source_verify_sha256=sha(root/'scripts/source_verify.py'),
        restored_independent_generator=True,no_Cpp_or_Fortran_changes=True))
    print(root,flush=True)


if __name__=='__main__':main(sys.argv[1])
