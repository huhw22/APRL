#!/usr/bin/env python3
"""Independent lifecycle, unchanged-physics, energy-accounting and timing review."""
import os
for name in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[name] = '1'
import argparse
import csv
import json
from pathlib import Path
import re
import statistics
import h5py
import numpy as np
import yaml
from prepare import canonical, require, save

LD = np.longdouble
EPS = LD(np.finfo(np.float64).eps)


def finite(array):
    a = np.asarray(array)
    if a.dtype.names:
        return all(finite(a[name]) for name in a.dtype.names)
    return bool(np.isfinite(a).all()) if a.dtype.kind in 'fc' else True


def energy_review(path, new, export_csv=True):
    with h5py.File(path) as stream:
        group = stream['energy_ledger']; rows = group['records'][:]
        attrs = {name: value.decode() if isinstance(value, bytes) else value.item() if isinstance(value, np.generic) else value
                 for name, value in group.attrs.items()}
    require('boosted computational frame' in str(attrs.get('frame')), 'Unexpected ledger frame')
    require(len(rows) > 1 and np.all(np.diff(rows['step']) > 0), 'Invalid ledger step sequence')
    def a(name):
        require(name in rows.dtype.names, 'Missing ledger component: ' + name)
        return np.asarray(rows[name], dtype=LD)
    kinetic = a('particle_kinetic_active_J') + a('particle_kinetic_removed_J')
    delta_k = kinetic-kinetic[0]
    field = a('field_energy_interior_J'); delta_u = field-field[0]
    faces = a('outward_field_energy_face_J')
    outward = np.sum(faces, axis=1)
    prescribed = a('prescribed_work_J'); boundary = a('boundary_source_work_J')
    raw = delta_k+delta_u+outward-prescribed
    expected = raw-boundary
    identities = {}
    source_scale = abs(kinetic)+abs(kinetic[0])+abs(field)+abs(field[0])+np.sum(abs(faces),axis=1)+abs(prescribed)+abs(boundary)
    if new:
        inserted = a('moving_inserted_field_J')
        cropped = a('moving_cropped_field_J')
        remap = a('moving_remap_field_change_J')
        correction = a('moving_remap_interface_correction_J')
        overwrite = a('inlet_overwrite_work_J')
        longitudinal = a('longitudinal_boundary_outward_J')
        expected -= remap+overwrite
        source_scale += abs(inserted)+abs(cropped)+abs(remap)+abs(correction)+abs(overwrite)
        identities['remap_equals_insert_minus_crop_plus_interface'] = (remap, inserted-cropped+correction)
        identities['longitudinal_groups_propagation_remap_and_overwrite'] = (longitudinal, np.sum(faces[:,4:],axis=1)-remap-overwrite)
        # Added constrained-current supply fields are handled explicitly once
        # present; they are not inferred from the residual to force a closure.
        boundary_corrected = expected.copy()
        constraint_total = np.zeros_like(expected)
        constraint_absolute_components = np.zeros_like(expected)
        for name in ('frozen_constraint_work_J', 'carrier_constraint_work_J'):
            expected -= a(name); source_scale += abs(a(name)); constraint_total += a(name)
            constraint_absolute_components += abs(a(name))
    identities['outward_equals_sum_of_faces'] = (a('outward_field_energy_J'), outward)
    identities['raw_residual_reconstructed'] = (a('balance_residual_raw_J'), raw)
    identities['corrected_residual_reconstructed'] = (a('balance_residual_J'), expected)
    identity_report = {}
    # This bound covers serialization and cancellation in independently
    # reconstructed sums; it is not a permitted physical conservation error.
    tolerance = 128*EPS*np.maximum(source_scale, LD(np.finfo(np.float64).tiny))
    for name,(reported,reconstructed) in identities.items():
        error = abs(reported-reconstructed)
        passed = bool(np.all(error <= tolerance))
        identity_report[name] = {'pass': passed, 'max_abs_difference_J': float(np.max(error)),
            'max_roundoff_bound_J': float(np.max(tolerance)), 'elements_over_bound': int(np.count_nonzero(error>tolerance))}
    residual = a('balance_residual_J')
    exchange = abs(delta_k)+abs(delta_u)+np.sum(abs(faces),axis=1)+abs(prescribed)+abs(boundary)
    code_exchange = abs(delta_k)+abs(delta_u)+abs(outward)+abs(prescribed)+abs(boundary)
    if new:
        exchange += abs(remap)+abs(overwrite)+constraint_absolute_components
        code_exchange += abs(remap)+abs(overwrite)+constraint_absolute_components
    relative = np.divide(residual, exchange, out=np.zeros_like(residual), where=exchange>0)
    code_relative = np.divide(expected, code_exchange, out=np.zeros_like(expected), where=code_exchange>0)
    reported_relative = a('relative_balance_exchange')
    relative_error = abs(reported_relative-code_relative)
    # Propagate the absolute source-scale summation bound through division.
    ratio_tolerance = np.zeros_like(tolerance)
    nonzero = code_exchange>0
    ratio_tolerance[nonzero] = tolerance[nonzero]/code_exchange[nonzero]*(2+abs(code_relative[nonzero]))+128*EPS
    identity_report['exchange_normalized_residual_reconstructed'] = {
        'pass': bool(np.all(relative_error<=ratio_tolerance)),
        'max_abs_difference': float(np.max(relative_error)), 'max_roundoff_bound': float(np.max(ratio_tolerance)),
        'elements_over_bound': int(np.count_nonzero(relative_error>ratio_tolerance))}
    names = [name for name in rows.dtype.names if rows.dtype[name].shape == () and
             (name.endswith('_J') or name in ('step','time_box_s','closure_valid','removed_macroparticles'))]
    if export_csv:
        with path.with_suffix('.csv').open('w', newline='') as stream:
            writer = csv.writer(stream); writer.writerow(names)
            writer.writerows([[row[name].item() for name in names] for row in rows])
    first_exit = np.flatnonzero(rows['removed_macroparticles']>0)
    windows = {'full': np.ones(len(rows), dtype=bool)}
    if len(first_exit):
        windows['before_first_exit'] = np.arange(len(rows))<first_exit[0]
        windows['from_first_exit'] = np.arange(len(rows))>=first_exit[0]
    phases = {}
    for name, mask in windows.items():
        indices = np.flatnonzero(mask)
        if not len(indices):
            continue
        last = indices[-1]
        phases[name] = {'first_step': int(rows['step'][indices[0]]), 'last_step': int(rows['step'][last]),
            'maximum_abs_residual_J': float(np.max(abs(residual[mask]))), 'last_residual_J': float(residual[last]),
            'maximum_abs_code_exchange_normalized_residual': float(np.max(abs(code_relative[mask]))),
            'maximum_abs_facewise_exchange_normalized_residual': float(np.max(abs(relative[mask])))}
    decomposition = {'raw_final_J': float(raw[-1]), 'raw_max_abs_J': float(np.max(abs(raw))),
        'corrected_final_J': float(expected[-1]), 'corrected_max_abs_J': float(np.max(abs(expected)))}
    if new:
        decomposition.update(boundary_only_final_J=float(boundary_corrected[-1]),
            boundary_only_max_abs_J=float(np.max(abs(boundary_corrected))),
            constraint_work_final_J=float(constraint_total[-1]))
    return {'checks_pass': all(x['pass'] for x in identity_report.values()), 'frame': attrs['frame'],
        'metadata': attrs, 'identities': identity_report, 'closure_valid_values': np.unique(rows['closure_valid']).tolist(),
        'final_components': {name: rows[name][-1].item() for name in names},
        'first_removed_step': int(rows['step'][first_exit[0]]) if len(first_exit) else None,
        'phases': phases, 'final_delta_particle_J': float(delta_k[-1]), 'final_delta_field_J': float(delta_u[-1]),
        'residual_decomposition': decomposition,
        'initial_kinetic_J': float(kinetic[0]), 'initial_field_J': float(field[0]),
        'conservation_acceptance': 'NOT_ASSIGNED: physical residual is reported, not gated by a newly chosen percentage.',
        'roundoff_identity_bound': '128*float64_eps*sum(abs(source energy terms)), including both initial and current kinetic energy.',
        'exchange_normalization': 'Code denominator uses abs(sum of six outward faces), abs(frozen work)+abs(carrier work); facewise denominator replaces abs(sum faces) by sum(abs(faces)).',
        'scope': 'Net signed control-volume energy, boost frame; no FFT/windowing or radiation-only interpretation.'}


def inspect(directory, item, plan):
    log = (directory/'aprl.log').read_text(errors='replace')
    require('Configured stop reached:' in log and 'phase=mpi-finalize-complete' in log, 'Missing normal physical completion')
    require(not re.search(r'Simulation failure|\[shutdown-watchdog\]|DUE TO TIME LIMIT|Segmentation fault|OpenFabrics device',log), 'Unexpected failure/warning')
    moving = re.search(r'Moving window completed: steps=(\d+) particle_exits=(\d+) cropped_carriers=(\d+)',log)
    require(moving and int(moving[1]) == item['steps'] and int(moving[2]) == plan['particles'], 'Incomplete moving lifecycle')
    released = re.search(r'Frozen lab-plane injection: activations=(\d+),.*still frozen=(\d+)',log)
    require(released and int(released[1]) == plan['particles'] and int(released[2]) == 0, 'Incomplete release lifecycle')
    require('[performance]' not in log, 'Profile must be disabled')
    files = sorted(directory.glob('detectors/*.h5'))+sorted(directory.glob('trajectory/*.h5'))+sorted(directory.glob('ledger/*.h5'))
    card = yaml.safe_load((directory/'run.yaml').read_text())
    science = card['trajectory']['enabled']
    expected_files = (20 if science else 1)+(1 if item['ledger'] else 0)
    require(len(files) == expected_files, 'Unexpected output file count')
    terminal_counts = np.zeros(plan['particles'],dtype=np.int64)
    checks = {}
    for path in files:
        with h5py.File(path) as stream:
            require(len(stream) == 1, 'Unexpected HDF5 group count')
            group = stream[next(iter(stream))]
            require(int(group['complete'][()]) == 1, 'Incomplete HDF5: '+str(path))
            for committed, values in (('committed_records','records'),('committed_samples','time_s')):
                if committed in group:
                    require(int(group[committed][()]) == len(group[values]), 'Truncated HDF5')
            for name,dataset in group.items():
                if not isinstance(dataset,h5py.Dataset):
                    continue
                chunks = [dataset[()]] if not dataset.shape else (dataset[i:i+4] if name=='fields' else dataset[i:i+65536]
                    for i in range(0,len(dataset),4 if name=='fields' else 65536))
                for chunk in chunks:
                    require(finite(chunk),'Nonfinite output: '+str(path)+'/'+name)
                    if path.parent.name=='trajectory' and name=='records':
                        final = chunk[chunk['event_type']!=0]
                        require(np.all(final['event_type']==2) and np.all(final['boundary_face']==5), 'Unexpected terminal event')
                        require(np.all((final['source_id']>=1)&(final['source_id']<=plan['particles'])), 'Unexpected terminal ID')
                        np.add.at(terminal_counts,final['source_id'].astype(np.int64)-1,1)
            if path.parent.name=='detectors' and path.stem!='exit-field':
                rows = group['records'][:]
                require(np.array_equal(np.sort(rows['source_id']),np.arange(1,plan['particles']+1)), 'Particle plane coverage differs')
                require(np.all(rows['crossing_direction']==1),'Unexpected backward crossing')
            checks[str(path.relative_to(directory))] = {'complete':True,'bytes':path.stat().st_size}
    if science:
        require(np.all(terminal_counts==1),'Missing/duplicate particle terminal events')
    points = {int(step):float(elapsed) for step,elapsed in re.findall(r'\[resource\] phase=(?:progress|final) step=(\d+) elapsed_s=([\d.eE+\-]+)',log)}
    require(40 in points and 240 in points and item['steps'] in points and points[240]>points[40], 'Missing timing interval')
    result = {'checks_pass':True,'steps':int(moving[1]),'particle_exits':int(moving[2]),'cropped_carriers':int(moving[3]),
        'files':checks,'progress_elapsed_s':points,'steps_40_240_seconds_per_step':(points[240]-points[40])/200,
        'all_steps_seconds_per_step':points[item['steps']]/item['steps'],'profile_enabled':False}
    preparation = json.loads((directory/'reference/preparation.json').read_text())
    result['reference_preparation_seconds'] = preparation['preparationSeconds']
    launched = json.loads((directory/'aprl.log.launch.json').read_text())
    result['launcher_seconds_including_preparation_and_close'] = launched['seconds']
    if item['ledger']:
        result['energy'] = energy_review(directory/'ledger/energy-ledger.h5',item['binary']=='new')
        require(result['energy']['checks_pass'], 'Energy-record algebra does not match independent reconstruction')
        work = result['energy']['final_components']['prescribed_work_J']
        require(work==0. if item['magnet']=='free-drift' else work!=0., 'External-device work coverage differs')
    return result


def compare(first,second):
    require(canonical(yaml.safe_load((first/'run.yaml').read_text())) ==
            canonical(yaml.safe_load((second/'run.yaml').read_text())), 'Comparison cards differ beyond ledger toggle/output paths')
    relative = [p.relative_to(first) for folder in ('detectors','trajectory') for p in sorted((first/folder).glob('*.h5'))]
    require(relative == [p.relative_to(second) for folder in ('detectors','trajectory') for p in sorted((second/folder).glob('*.h5'))], 'Compared outputs differ')
    result = {}
    for rel in relative:
        equal = True; samples = 0
        with h5py.File(first/rel) as fa,h5py.File(second/rel) as fb:
            a,b = fa[next(iter(fa))],fb[next(iter(fb))]
            for name in ('time_s','fields','records'):
                if name not in a:
                    require(name not in b,'Compared dataset absent')
                    continue
                require(name in b and a[name].shape == b[name].shape and a[name].dtype == b[name].dtype,'Compared schema differs')
                # Detector record append order may vary with rank; sort the
                # cohort by immutable source ID. Trajectory rank ownership is fixed.
                if rel.parts[0]=='detectors' and name=='records':
                    x,y = a[name][:],b[name][:]
                    x=x[np.argsort(x['source_id'])];y=y[np.argsort(y['source_id'])]
                    equal &= x.tobytes()==y.tobytes(); samples+=len(x)
                else:
                    for i in range(0,len(a[name]),4 if name=='fields' else 65536):
                        n=4 if name=='fields' else 65536
                        x,y=a[name][i:i+n],b[name][i:i+n]
                        equal &= x.tobytes()==y.tobytes();samples+=len(x)
        result[str(rel)] = {'bitwise_equal':bool(equal),'records_or_samples_compared':samples}
    return {'checks_pass':all(v['bitwise_equal'] for v in result.values()),'files':result}


def stage_review(root,stage,plan=None):
    plan = plan or json.loads((root/'plan.json').read_text())
    items=plan['cases'][stage]
    cases={item['name']:json.loads((Path(item['directory'])/'CASE_REVIEW.json').read_text()) for item in items}
    comparisons={}
    for magnet in ('undulator','free-drift'):
        for repeat in sorted({item['repeat'] for item in items}):
            group=[item for item in items if item['magnet']==magnet and item['repeat']==repeat]
            first=next(x for x in group if x['variant']=='old-on')
            for other in group:
                if other==first:
                    continue
                comparisons[other['name']]=compare(Path(first['directory']),Path(other['directory']))
    result={'format':'aprl-moving-energy-review-v1','stage':stage,'binary_sha256':plan['binary_sha256'],
        'checks_pass':all(x['checks_pass'] for x in cases.values()) and all(x['checks_pass'] for x in comparisons.values()),
        'cases':cases,'physical_output_comparisons':comparisons,'energy_conservation_accepted':False,
        'note':'checks_pass covers completion, exact unchanged physical outputs and accounting identities; it does not declare a small physical energy residual.'}
    if stage=='performance':
        metrics={}
        for magnet in ('undulator','free-drift'):
            samples={v:[cases[x['name']]['steps_40_240_seconds_per_step'] for x in items if x['magnet']==magnet and x['variant']==v]
                     for v in ('old-on','new-on','old-off','new-off')}
            summary={v:{'samples_s_per_step':a,'median_s_per_step':statistics.median(a),'min_s_per_step':min(a),'max_s_per_step':max(a)} for v,a in samples.items()}
            ratios={}
            for name,num,den in (('increment_on_new_vs_old','new-on','old-on'),('disabled_new_vs_old','new-off','old-off'),
                                 ('new_total_diagnostic_on_vs_off','new-on','new-off'),('old_total_diagnostic_on_vs_off','old-on','old-off')):
                paired=[100*(a/b-1) for a,b in zip(samples[num],samples[den])]
                ratios[name]={'paired_percent_changes':paired,'median_paired_percent_change':statistics.median(paired),
                    'ratio_of_medians_percent':100*(statistics.median(samples[num])/statistics.median(samples[den])-1)}
            metrics[magnet]={'configurations':summary,'comparisons':ratios}
        result['performance']=metrics
        result['performance_scope']='Same-node fixed16x8, no profile, complete257-step small case; 40–240 steps exclude initialization and early warmup. Not production overhead.'
    save(root/(stage.upper()+'_REVIEW.json'),result)
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,required=True)
    parser.add_argument('--stage',choices=('science','performance','mpi'),required=True)
    args=parser.parse_args()
    result=stage_review(args.root.resolve(),args.stage)
    require(result['checks_pass'],'Stage review failed')
