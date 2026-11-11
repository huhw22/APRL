"""Reuse the immutable 11-segment particles; prepare a six-segment ledger pair."""
import copy
import hashlib
import json
import math
import shutil
import sys
from pathlib import Path
import h5py
import numpy as np
import yaml
from scipy.constants import c

W = Path('/public/home/wfollowm/APRL_Workspace')
PREVIOUS = W/'results/sparc_sase_11segment_2slip_20261111_r1'
OBSERVER = W/'results/lab_worldtube_balance_20261112_r1'
TEMPLATE = Path(__file__).resolve().parent


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(4*1024*1024), b''):
            digest.update(block)
    return digest.hexdigest()


def save(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')


def replace_once(text, old, new):
    if text.count(old) != 1:
        raise RuntimeError('Unexpected overlay anchor: '+old)
    return text.replace(old, new, 1)


def main(root):
    root = Path(root).resolve()
    root.mkdir(exist_ok=False)
    scripts = root/'scripts'
    shutil.copytree(TEMPLATE, scripts, ignore=shutil.ignore_patterns('__pycache__'))
    source = root/'source'
    shutil.copytree(OBSERVER/'source', source, ignore=shutil.ignore_patterns('__pycache__','.git'))
    shutil.copy2(scripts/'lab_side_particle_audit.h', source/'src/lab_side_particle_audit.h')
    particle = source/'src/particle_energy_audit.h'
    text = particle.read_text()
    text = replace_once(text, '#include "particle_work_geometry.h"',
                         '#include "particle_work_geometry.h"\n#include "lab_side_particle_audit.h"')
    text = replace_once(text, '  std::vector<EnergyPlaneCrossing> crossings;',
                         '  std::vector<EnergyPlaneCrossing> crossings;\n  std::vector<LabSideParticleEvent> sideCrossings;')
    text = replace_once(text, 'segments.clear(); crossings.clear(); substepCrossings.clear(); observedSubsteps=0;',
                         'segments.clear(); crossings.clear(); sideCrossings.clear(); substepCrossings.clear(); observedSubsteps=0;')
    text = replace_once(text, '    if (recordWork) {\n      Segment s;',
        '    if(recordWork && LabSideParticleAudit::enabled())\n      LabSideParticleAudit::preview(before,after,begin,end,*frame,sideCrossings);\n    if (recordWork) {\n      Segment s;')
    particle.write_text(text)
    simulation = source/'src/simulation.cpp'
    text = simulation.read_text()
    anchor = '       logRoot(communicator_,"[lab-worldtube] observer only; equal-time fields, all laboratory worldtube surfaces; no physical-state writes");}}'
    text = replace_once(text, anchor,
        '       LabSideParticleAudit::open(dir,rank_,(.5*globalGeometry_.nx-config_.boundary.cpml.cells[0]-1)*globalGeometry_.dx,\n'
        '         (.5*globalGeometry_.ny-config_.boundary.cpml.cells[1]-1)*globalGeometry_.dy);\n'
        '       logRoot(communicator_,"[lab-particle-sides] substep crossings, endpoint scalar lab energy, terminal-clipped commit enabled");\n'+anchor)
    text = replace_once(text, '    if(worldtube)worldtube->close(!interrupted_&&configuredStopReached_);',
        '    LabSideParticleAudit::close();\n    if(worldtube)worldtube->close(!interrupted_&&configuredStopReached_);')
    anchor = '            if (detectors_ && hasActiveSegment &&\n                activeStartTime < diagnosticTime)'
    text = replace_once(text, anchor,
        '            if(config_.energyLedger.enabled && hasActiveSegment)\n'
        '              LabSideParticleAudit::commit(result.auditTrace.sideCrossings,\n'
        '                hasTerminalEvent ? diagnosticTime : std::numeric_limits<Double>::infinity());\n'+anchor)
    anchor = '          if (detectors_ && remaining > 0.0)\n            detectors_->captureParticleStep('
    text = replace_once(text, anchor,
        '          if(config_.energyLedger.enabled && remaining > 0.)\n'
        '            LabSideParticleAudit::commit(auditTrace.sideCrossings,\n'
        '              auditTerminal ? auditStop : std::numeric_limits<Double>::infinity());\n'+anchor)
    simulation.write_text(text)
    changed = []
    for path in (OBSERVER/'source').rglob('*'):
        if path.is_file() and '__pycache__' not in path.parts and sha(path) != sha(source/path.relative_to(OBSERVER/'source')):
            changed.append(str(path.relative_to(OBSERVER/'source')))
    assert set(changed) == {'src/simulation.cpp','src/particle_energy_audit.h'}
    (root/'bin').mkdir()
    old = json.loads((PREVIOUS/'plan.json').read_text())
    ready = json.loads((PREVIOUS/'INPUT_READY.json').read_text())
    input_file = Path(ready['input_file'])
    assert sha(input_file) == ready['input_sha256']
    with h5py.File(input_file, 'r') as handle:
        data = handle['particles/records']
        assert len(data) == old['input_particles']
        t = data.fields('arrival_time_offset_s')[:]
        span_input = float(t.max()-t.min())
    base = yaml.safe_load((PREVIOUS/'energy-binomial3/run.yaml').read_text())
    G = float(old['physics']['boost_gamma'])
    beta = math.sqrt(1-1/G**2)
    dz, dt = base['mesh']['cell_size'][2], old['geometry']['dt_s']
    planes = old['geometry']['field_planes'][:6]
    last = planes[-1][1]
    # Same causal coverage formula as the accepted 11-segment campaign, with
    # six undulators and the identical input time support. No beam regeneration.
    duration_lab = span_input+6*77*old['physics']['wavelength_m']/c+(1-beta)*last/c+4e-15
    steps = int(math.ceil(duration_lab*G/dt))
    active = int(math.ceil(last/(G*dz)))+64+128+2
    stop = old['geometry']['stop_reference_z_m']-math.sqrt(G*G-1)*dz*(old['geometry']['steps']-steps)
    magnets = [copy.deepcopy(x) for x in base['sources']['magnetic_elements']
               if x['entrance_z'] < last]
    assert sum(x['type']=='planar-undulator' for x in magnets)==6
    audit_planes = [('audit-entry', .03), *[(n, float(z)) for n,z in planes]]
    cases = []
    for name, shape in [('energy-binomial3','linear-binomial3'),('energy-momentum','quintic-compensated')]:
        directory = root/name
        directory.mkdir()
        for sub in ['detectors','ledger','lab-audit','worldtube']:
            (directory/sub).mkdir()
        card = copy.deepcopy(base)
        card['output']['manifest'] = str(directory/'run-manifest.yaml')
        card['mesh']['duration'] = (steps+.5)*dt
        card['mesh']['longitudinal_shape'] = shape
        card['mesh']['moving_window'].update(cells=active, reference_directory=str(directory/'reference'))
        card['sources']['magnetic_elements'] = magnets
        card['stop'] = dict(mode='reference-center-z', z=stop)
        card['detectors'].update(directory=str(directory/'detectors'),
            field_planes=[p for p in card['detectors']['field_planes'] if p['name'] in dict(planes)],
            particle_planes=[dict(name=n+'-particles',z=z,buffer_records=32768,compression=0)
                for n,z in [('entry',0.),('und-entry',.01),*audit_planes]],
            substep_particle_planes=True,
            energy_planes=[dict(name=n+'-transport', z=z, rhythm=dt/G,
                buffer_samples=64, compression=0) for n,z in [('entry',0.),('und-entry',.01),*audit_planes]])
        card['energy_ledger'].update(directory=str(directory/'ledger'), sample_interval_steps=1,
                                     discrete_operator_audit=True)
        (directory/'run.yaml').write_text(yaml.safe_dump(card, sort_keys=False))
        cases.append(dict(name=name,shape=shape,directory=str(directory),expected_steps=steps))
    # Short full-input off/on runs test diagnostic invariance before production.
    preflight = []
    for name, shape, observer in [('binomial3-off','linear-binomial3',False),
                                   ('binomial3-on','linear-binomial3',True),
                                   ('momentum-on','quintic-compensated',True)]:
        directory = root/'preflight'/name
        directory.mkdir(parents=True)
        for sub in ['detectors','ledger','lab-audit','worldtube']:
            (directory/sub).mkdir()
        card = copy.deepcopy(yaml.safe_load((root/('energy-momentum' if shape=='quintic-compensated' else 'energy-binomial3')/'run.yaml').read_text()))
        card['output']['manifest']=str(directory/'run-manifest.yaml')
        card['mesh']['moving_window']['reference_directory']=str(directory/'reference')
        card['detectors']['directory']=str(directory/'detectors')
        card['energy_ledger']['directory']=str(directory/'ledger')
        # Preserve the full six-segment geometry and its valid downstream stop.
        # The 40-step preflight intentionally ends by the supported duration
        # cutoff, whose closed-but-incomplete output is explicitly checked.
        card['mesh']['duration']=(40-.25)*dt
        card['energy_ledger']['discrete_operator_audit']=observer
        (directory/'run.yaml').write_text(yaml.safe_dump(card,sort_keys=False))
        preflight.append(dict(name=name,directory=str(directory),observer=observer,expected_steps=40))
    plan = dict(complete_preparation=True,root=str(root),previous_campaign=str(PREVIOUS),
        source=str(source),binary=str(root/'bin/aprl'),
        library=str(W/'build/longitudinal_shape_lab_balance_20261110_r4/libdirect_eb_core.a'),
        main_object=str(W/'build/longitudinal_shape_lab_balance_20261110_r4/CMakeFiles/aprl.dir/src/main.cpp.o'),
        input_file=str(input_file),input_sha256=ready['input_sha256'],input_particles=old['input_particles'],
        input_represented_electrons=old['input_electrons'],physics=old['physics'],
        cases=cases,preflight=preflight,steps=steps,active_nz=active,
        reference_nz=base['mesh']['cells'][2],dt_box_s=dt,dz_box_m=dz,
        stop_reference_z_m=stop,segments=6,mpi=16,threads=8,
        worldtube_entrance_z_m=.03,worldtube_exit_planes=planes,
        worldtube_half_width_m=(20-6-1)*base['mesh']['cell_size'][0],
        lab_field_raw_upper_bytes_per_case=steps*active*6*8,
        lab_physical_scope='All optical frequencies, full common retarded-time interval; independently measured end ports, four sides, two null-time caps and side-particle kinetic transport.',
        source_exclusion='Frozen source support <=9*boost_gamma*dz <0.03m; transverse prescribed current sheets at ±560um outside ±520um observer faces; lab devices are magnetostatic.',
        numerical_scope='Native boost discrete ledger plus lab field momentum/stress and all Cowan/CPML/inlet/remap exchanges retained.',
        side_event_scope='Actual committed substep worldlines; scalar endpoint lab energies; outward signs; terminal cutoff; trial trajectories excluded.',
        original_input_unchanged=True,physics_core_unchanged=True,
        observer_source_changes=changed+['src/lab_side_particle_audit.h'])
    assert 9*G*dz < .03
    save(root/'PLAN.json',plan)
    paths=[root/'PLAN.json', input_file, Path(plan['library']),Path(plan['main_object'])]
    paths += [p for p in source.rglob('*') if p.is_file()]
    paths += [root/'scripts'/p.name for p in scripts.iterdir() if p.is_file()]
    paths += [Path(x['directory'])/'run.yaml' for x in cases+preflight]
    save(root/'FROZEN.json',{str(path):sha(path) for path in paths})
    print(json.dumps({k:plan[k] for k in ['root','input_particles','steps','active_nz',
        'stop_reference_z_m','worldtube_half_width_m','lab_field_raw_upper_bytes_per_case']},indent=2))


if __name__=='__main__':main(sys.argv[1])
