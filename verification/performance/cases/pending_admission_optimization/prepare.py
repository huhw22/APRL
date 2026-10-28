#!/usr/bin/env python3
"""Freeze pending-admission validation cards and executables; never submit/run MPI."""
import argparse
import copy
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import yaml
from ledger_review import thresholds as ledger_thresholds

HERE = Path(__file__).resolve().parent
BASELINE_SHA = 'c81baefb2953ad04181621e31e0f60644872a562bd5bccbd1f45961e44844f11'
C = 299792458.


def require(test, message):
    if not test: raise RuntimeError(message)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''): h.update(chunk)
    return h.hexdigest()


def save(path, value):
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + '\n')


def redirect(original, directory, input_file=None):
    card = copy.deepcopy(original)
    card['output']['manifest'] = str(directory / 'run-manifest.yaml')
    card['mesh']['moving_window']['reference_directory'] = str(directory / 'reference')
    for key, folder in (('detectors', 'detectors'), ('energy_ledger', 'ledger'), ('trajectory', 'trajectory')):
        if card.get(key, {}).get('enabled'):
            card[key]['directory'] = str(directory / folder)
    if input_file: card['beam']['input']['file'] = str(input_file)
    require(not card.get('checkpoint', {}).get('enabled', False) and not card.get('checkpoint', {}).get('steps'), 'No snapshot in validation cards')
    return card


def bounded_stop(card, geo, steps):
    # Preserve the original physical stop. Put duration one ULP below the
    # smaller of N*dt and N repeated double additions, strictly above stepN-1.
    # This avoids an accidental N+1 step from the while(time < duration) test.
    dt = float(card['mesh']['cell_size'][2]) / C
    accumulated = 0.
    for unused in range(steps): accumulated += dt
    duration = math.nextafter(min(steps * dt, accumulated), 0.)
    require((steps - 1) * dt < duration <= steps * dt, 'Invalid bounded duration')
    card['mesh']['duration'] = duration
    return {'steps': steps, 'dt_box_s': dt, 'nominal_N_dt_s': steps * dt,
            'repeated_addition_N_dt_s': accumulated, 'duration_box_s': duration,
            'physical_stop_unchanged': True, 'expected_exit_code': 1,
            'status': 'TIMING_COMPLETE_PHYSICS_STOP_NOT_REACHED',
            'scope': 'N complete boost steps; known duration exhaustion with closed complete=false diagnostics, not a completed physical simulation.'}


def validate_science_reuse(parent, binary_sha256, small_input_sha256):
    """Verify immutable previous science evidence; never manufacture a new PASS."""
    parent = Path(parent).resolve()
    old_plan = json.loads((parent/'plan.json').read_text())
    passed = json.loads((parent/'SCIENCE_PASS.json').read_text())
    require(passed['pass'] and passed['stage']=='science' and
            passed['binary_sha256']==old_plan['binary_sha256']==binary_sha256,
            'Reused science does not cover both exact binaries')
    require(passed['frozen_manifest_sha256']==sha(parent/'frozen.json'),
            'Reused science frozen identity differs')
    require(old_plan['small_input_sha256']==small_input_sha256 and
            sha(parent/'inputs/held.h5')==small_input_sha256, 'Reused small input differs')
    expected={f'{magnet}-{ranks}x8-{variant}' for magnet in ('undulator','free-drift')
              for ranks in (1,16) for variant in ('baseline','candidate')}
    require(set(passed['cases'])==expected and all(item['pass'] for item in passed['cases'].values()),
            'Reused science did not pass all eight original cases')
    expected_pairs={name for name in expected if name.endswith('-candidate')} | {
        magnet+'-candidate-cross-layout' for magnet in ('undulator','free-drift')}
    require(set(passed['comparisons'])==expected_pairs and all(item['pass'] for item in passed['comparisons'].values()),
            'Reused science comparison coverage incomplete')
    thresholds=json.loads((parent/'THRESHOLDS.json').read_text())
    require(thresholds['field_relative_l2']==2e-12 and thresholds['particle_x_u_t_factor_eps64']==128 and
            thresholds['science_ledger']==ledger_thresholds(31675) and
            thresholds['full_ledger']==ledger_thresholds(41732775), 'Reused numerical gates differ')
    for path, expected_sha in json.loads((parent/'frozen.json').read_text()).items():
        require(sha(path)==expected_sha, 'Prior frozen science artifact changed: '+path)
    evidence=[parent/name for name in ('plan.json','SCIENCE_PASS.json','frozen.json','THRESHOLDS.json')]
    return {'source_root':str(parent),'source_science_pass':str(parent/'SCIENCE_PASS.json'),
        'evidence_sha256':{str(path):sha(path) for path in evidence},
        'binary_sha256':binary_sha256,'small_input_sha256':small_input_sha256,
        'status':'REUSED_PRIOR_SCIENCE_ACCEPTANCE',
        'scope':'Exact prior eight-case science acceptance reused for a fresh same-node full pair; no new local SCIENCE_PASS is fabricated.'}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--workspace', type=Path, required=True)
    p.add_argument('--root', type=Path, required=True)
    p.add_argument('--candidate', type=Path, required=True)
    p.add_argument('--reuse-science-from', type=Path,
                   help='Prepare only fresh full pairs, reusing immutable accepted same-binary science')
    args = p.parse_args(); w, root = args.workspace.resolve(), args.root.resolve()
    require(not root.exists(), 'Use a new campaign directory')
    binaries = {'baseline': w/'build/moving-window-integration-r3-20261028/aprl', 'candidate': args.candidate.resolve()}
    require(sha(binaries['baseline']) == BASELINE_SHA, 'Baseline R3 changed')
    for binary in binaries.values():
        require(binary.is_file() and os.access(binary, os.X_OK), 'Missing executable: '+str(binary))
        require((binary.parent/'runtime-libs').is_dir(), 'Freeze runtime-libs beside each executable')
    root.mkdir(parents=True); scripts = root/'scripts'; scripts.mkdir(); inputs = root/'inputs'; inputs.mkdir()
    for path in HERE.iterdir():
        if path.is_file(): shutil.copy2(path, scripts/path.name)
    review_source = HERE.parent/'moving_window_production/preflight.py'
    shutil.copy2(review_source, scripts/'smoke_review.py')
    shutil.copy2(HERE.parent/'full_scale_followup/run.py', scripts/'duration_review.py')
    source = w/'results/moving_window_integration_20261028/dense-two-node-r3'
    shutil.copy2(source/'inputs/held.h5', inputs/'held.h5')
    reuse = validate_science_reuse(args.reuse_science_from,
        {key:sha(binary) for key,binary in binaries.items()},sha(inputs/'held.h5')) if args.reuse_science_from else None
    small = yaml.safe_load((source/'moving/run.yaml').read_text())
    production = w/'results/moving_window_production_20261028/undulator_4p5m'
    original = yaml.safe_load((production/'run.yaml').read_text())
    geo = json.loads((production/'geometry-plan.json').read_text())
    prior = json.loads((production/'run-plan.json').read_text())
    require(prior['records'] == 41732775 and original['mesh']['cells'] == [160,160,53504], 'Original full input/reference geometry changed')
    require(not original['detectors']['particle_planes'] and not original['trajectory']['enabled'], 'Preserve original full diagnostics')
    cases = {'science': [], 'full': []}; stops = {}
    for ranks in (() if reuse else (1,16)):
        for magnet in ('undulator','free-drift'):
            for variant in ('baseline','candidate'):
                name = f'{magnet}-{ranks}x8-{variant}'; directory=root/'science'/name; directory.mkdir(parents=True)
                card=redirect(small,directory,inputs/'held.h5')
                if magnet=='free-drift': card['sources']['magnetic_elements']=[]
                (directory/'run.yaml').write_text(yaml.safe_dump(card,sort_keys=False))
                cases['science'].append({'name':name,'directory':str(directory),'ranks':ranks,'threads':8,
                    'variant':variant,'magnet':magnet,'profile_interval':0,'steps':257,'timeout_s':240})
    for mode,steps,interval in (('profile',20,10),('throughput',100,0)):
        for variant in ('baseline','candidate'):
            name=f'{mode}-{variant}'; directory=root/'full'/name; directory.mkdir(parents=True)
            card=redirect(original,directory)
            stops[name]=bounded_stop(card,geo,steps)
            require(card['beam']==original['beam'] and card['sources']==original['sources'] and
                    card['mesh']['cells']==original['mesh']['cells'] and card['detectors']['field_planes']==original['detectors']['field_planes'], 'Full physics changed')
            if reuse:
                old_card=yaml.safe_load((Path(reuse['source_root'])/'full'/name/'run.yaml').read_text())
                require(card==redirect(old_card,directory), 'Fresh full card differs from previous pair beyond output/reference paths')
            (directory/'run.yaml').write_text(yaml.safe_dump(card,sort_keys=False))
            cases['full'].append({'name':name,'directory':str(directory),'ranks':16,'threads':8,'variant':variant,
                'magnet':'undulator','profile_interval':interval,'steps':steps,'timeout_s':600})
    binding=w/'build/cluster-tools/binding_probe'
    inp=Path(prior['input_file']); st=inp.stat()
    require((st.st_size,st.st_mtime_ns)==(prior['input_size_bytes'],prior['input_mtime_ns']), 'Input changed since production audit')
    plan={'format':'aprl-pending-admission-validation-v1','workspace':str(w),'binaries':{k:str(v) for k,v in binaries.items()},
        'binary_sha256':{k:sha(v) for k,v in binaries.items()},'binding_probe':str(binding),'cases':cases,
        'source_revision':subprocess.check_output(['git','-C',str(w/'source/aprl'),'rev-parse','HEAD'],text=True).strip(),
        'review_helper_source':str(review_source),'review_helper_sha256':sha(review_source),'small_input_sha256':sha(inputs/'held.h5'),
        'full_source_card':str(production/'run.yaml'),'full_source_card_sha256':sha(production/'run.yaml'),
        'full_input_file':str(inp),'full_input_sha256':prior['input_sha256'],'full_input_size_bytes':st.st_size,
        'full_input_mtime_ns':st.st_mtime_ns,'full_input_records':prior['records'],'full_stops':stops,
        'stage_budget_seconds':1800,'slurm_hard_limit_minutes':40,'full_case_limit_seconds':600,
        'same_mpi_pair':True,'mpi_environment':{'OMPI_MCA_btl':'self,vader,tcp','OMPI_MCA_btl_openib_warn_no_device_params_found':'0'},
        'acceptance':{'field_relative_l2':2e-12,'particle_roundoff_factor_eps64':128,'ledger_contract':'THRESHOLDS.json: source-scale long-double gamma_N and propagated absolute bounds'},
        'limits':['Short validation and early100-step throughput only; not full production radiation/late-stage acceptance.',
                  'Profile20-step pair is diagnostic only; performance uses separate profile-disabled100-step pair.',
                  'No layout scan;1rank is correctness coverage, production comparison stays16x8.',
                  'Every full case reconstructs its own identical reference; initialization is reported separately.',
                  'No production launch or follow-up scheduling is performed by these scripts.']}
    if reuse:
        plan['science_reuse']=reuse
        plan['limits'].append('Prior eight-case same-binary science is reused explicitly; all four full cases run freshly on the new single allocation with no node pin.')
    save(root/'plan.json',plan)
    save(root/'THRESHOLDS.json',{'declared_before_execution':True,
        'field_relative_l2':2e-12,'particle_x_u_t_factor_eps64':128,
        'science_ledger':ledger_thresholds(31675),'full_ledger':ledger_thresholds(41732775)})
    frozen=[root/'plan.json',root/'THRESHOLDS.json',inputs/'held.h5',binding]+list(scripts.iterdir())
    frozen += [Path(item['directory'])/'run.yaml' for stage in cases.values() for item in stage]
    for binary in binaries.values(): frozen += [binary]+[x for x in (binary.parent/'runtime-libs').iterdir() if x.is_file()]
    if reuse: frozen += [Path(path) for path in reuse['evidence_sha256']]
    save(root/'frozen.json',{str(path):sha(path) for path in frozen if path.is_file()})
    print(root/'plan.json')


if __name__=='__main__': main()
