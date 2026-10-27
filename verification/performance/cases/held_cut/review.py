#!/usr/bin/env python3
"""Independent offline audit of production-parameter held-inlet/cut probes.

Structural gates are explicit. This script reports scientific differences and
does not invent a tolerance that would turn them into a physical acceptance.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import re

import numpy as np
import yaml


WORK = Path('/public/home/wfollowm/APRL_Workspace')
PRODUCTION = WORK / 'results/elegant_field_only_20261027/production'
C = np.longdouble('299792458')
REST_EV = np.longdouble('9.1093837139e-31') * C*C / np.longdouble('1.602176634e-19')


def require(value, message):
    if not value:
        raise RuntimeError(message)


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for data in iter(lambda: stream.read(1024*1024), b''):
            h.update(data)
    return h.hexdigest()


def save(path, value):
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True,
                                    allow_nan=False) + '\n')


def table(path):
    with Path(path).open(newline='') as stream:
        reader = csv.DictReader(stream)
        require(reader.fieldnames and len(set(reader.fieldnames)) == len(reader.fieldnames),
                'Missing/duplicate CSV columns: '+str(path))
        rows = list(reader)
    require(all(None not in row and all(v is not None and v != '' for v in row.values())
                for row in rows), 'Truncated CSV row: '+str(path))
    result = {name: np.asarray([row[name] for row in rows], dtype=np.longdouble)
              for name in reader.fieldnames}
    require(all(np.all(np.isfinite(a)) for a in result.values()),
            'Non-finite CSV: '+str(path))
    return result


def column(t, *names):
    for name in names:
        if name in t:
            return t[name]
    raise RuntimeError('Missing column '+repr(names)+'; available '+repr(list(t)))


def size(t):
    return len(next(iter(t.values())))


def order_by_id(t):
    ids = column(t, 'id')
    require(len(np.unique(ids)) == len(ids), 'Duplicate particle ID')
    require(np.all(ids == np.floor(ids)), 'Noninteger particle ID')
    order = np.argsort(ids)
    return {k: v[order] for k, v in t.items()}


def rows_for_ids(t, ids):
    source = column(t, 'id')
    lookup = {int(v): i for i, v in enumerate(source)}
    require(all(int(v) in lookup for v in ids), 'Event ID missing from initial cohort')
    indices = np.asarray([lookup[int(v)] for v in ids], dtype=int)
    return {k: v[indices] for k, v in t.items()}


def same_tables(a, b, what):
    require(a.keys() == b.keys(), what+' columns differ')
    require(all(np.array_equal(a[k], b[k]) for k in a), what+' values differ')


def vector(t, names):
    return np.column_stack([column(t, n) for n in names])


def particle_momentum(t):
    # RelativisticParticleSI::properVelocity is FieldVector<Double>, where
    # Double is double. The producer writes 17 significant decimal digits.
    # Recover that original binary64 value before widening for cancellation-
    # resistant energy arithmetic. Do not apply this to long-double delta
    # energy columns, which retain their direct long-double CSV parsing.
    return vector(t, ('ux', 'uy', 'uz')).astype(np.float64).astype(np.longdouble)


def weighted_mean(x, w):
    return np.sum(x*w, dtype=np.longdouble)


def weighted_sigma(x, w):
    # Two passes around the weighted mean preserve a small spread on a large
    # gamma; never subtract two squared large energies.
    d = x-weighted_mean(x, w)
    return np.sqrt(weighted_mean(d*d, w))


def particle_summary(t, initial, boost, beta):
    if size(t) == 0:
        return {'count': 0}
    initial = rows_for_ids(initial, column(t, 'id'))
    charge = np.abs(column(t, 'charge_C'))
    require(np.all(charge > 0), 'Nonpositive particle weight')
    w = charge/np.sum(charge, dtype=np.longdouble)
    u = particle_momentum(t)
    u0 = particle_momentum(initial)
    g = np.sqrt(1+np.sum(u*u, axis=1, dtype=np.longdouble))
    g0 = np.sqrt(1+np.sum(u0*u0, axis=1, dtype=np.longdouble))
    dg = np.sum((u-u0)*(u+u0), axis=1, dtype=np.longdouble)/(g+g0)
    delta = boost*(dg+beta*(u[:, 2]-u0[:, 2]))*REST_EV
    reported = column(t, 'delta_energy_lab_eV')
    lab_gamma = boost*(g+beta*u[:, 2])
    result = {
        'count': size(t), 'absolute_charge_C': float(np.sum(charge)),
        'mean_delta_energy_lab_eV': float(weighted_mean(delta, w)),
        'sigma_delta_energy_lab_eV': float(weighted_sigma(delta, w)),
        'rms_delta_energy_lab_eV': float(np.sqrt(weighted_mean(delta*delta, w))),
        'mean_gamma_lab': float(weighted_mean(lab_gamma, w)),
        'sigma_gamma_lab': float(weighted_sigma(lab_gamma, w)),
        'sigma_energy_lab_eV': float(weighted_sigma(lab_gamma, w)*REST_EV),
        'relative_energy_spread_lab': float(weighted_sigma(lab_gamma, w)/weighted_mean(lab_gamma, w)),
        'reported_delta_energy_max_abs_discrepancy_eV': float(np.max(np.abs(reported-delta))),
        'momentum_changed_from_input': bool(np.any(u != u0)),
        'particle_energy_exchange_nonzero': bool(np.any(delta != 0)),
    }
    for axis in ('x', 'y'):
        coordinate = column(t, axis+'_m')
        before = column(initial, axis+'_m')
        rms, rms0 = weighted_sigma(coordinate, w), weighted_sigma(before, w)
        result.update({axis+'_centroid_m': float(weighted_mean(coordinate, w)),
                       axis+'_rms_m': float(rms), axis+'_initial_rms_m': float(rms0),
                       axis+'_rms_growth_m': float(rms-rms0)})
    return result


def particle_comparison(reference, trial, initial, boost, beta):
    require(np.array_equal(column(reference, 'id'), column(trial, 'id')),
            'Compared particle cohorts differ')
    if size(reference) == 0:
        return {'count': 0}
    charge = np.abs(column(reference, 'charge_C'))
    require(np.array_equal(charge, np.abs(column(trial, 'charge_C'))), 'Particle weights differ')
    w = charge/np.sum(charge, dtype=np.longdouble)
    a = particle_momentum(reference)
    b = particle_momentum(trial)
    ga, gb = (np.sqrt(1+np.sum(v*v, axis=1, dtype=np.longdouble)) for v in (a,b))
    difference = boost*(np.sum((b-a)*(b+a), axis=1, dtype=np.longdouble)/(ga+gb)
                        + beta*(b[:, 2]-a[:, 2]))*REST_EV
    base = particle_summary(reference, initial, boost, beta)
    current = particle_summary(trial, initial, boost, beta)
    return {
        'count': size(reference), 'reference': base, 'trial': current,
        'mean_energy_difference_eV': float(weighted_mean(difference, w)),
        'rms_energy_difference_eV': float(np.sqrt(weighted_mean(difference*difference, w))),
        'max_abs_energy_difference_eV': float(np.max(np.abs(difference))),
        'sigma_energy_difference_eV': current['sigma_energy_lab_eV']-base['sigma_energy_lab_eV'],
        'x_rms_growth_difference_m': current['x_rms_growth_m']-base['x_rms_growth_m'],
        'y_rms_growth_difference_m': current['y_rms_growth_m']-base['y_rms_growth_m'],
        'max_abs_position_difference_m': float(np.max(np.abs(
            vector(trial, ('x_m','y_m','z_box_m'))-vector(reference, ('x_m','y_m','z_box_m'))))),
    }


def field_vector(t, kind):
    aliases = {a: (kind+a, kind+a.lower(), kind+a.lower()+'_V_m',
                    kind+a+'_V_m', kind+a.lower()+'_T', kind+a+'_T',
                    (kind+a).lower()+'_V_m', (kind+a).lower()+'_T')
               for a in 'xyz'}
    return np.column_stack([column(t, *aliases[a]) for a in 'xyz'])


def norm_metrics(reference, trial):
    difference = trial-reference
    numerator = np.sum(difference*difference, dtype=np.longdouble)
    denominator = np.sum(reference*reference, dtype=np.longdouble)
    count = len(reference)
    return {'samples': count,
            'relative_l2': float(np.sqrt(numerator/denominator)) if denominator > 0 else None,
            'vector_rms_difference': float(np.sqrt(numerator/count)) if count else None,
            'reference_vector_rms': float(np.sqrt(denominator/count)) if count else None,
            'max_abs_component_difference': float(np.max(np.abs(difference))) if count else None,
            'all_values_identical': bool(np.array_equal(reference, trial))}


def fields_comparison(base, trial, first_release):
    require(base.keys() == trial.keys() and size(base) == size(trial), 'Field-plane schema/size differ')
    field_names = set()
    for kind in ('E','B'):
        for name in base:
            if name.lower().startswith(tuple((kind+a).lower() for a in 'xyz')):
                field_names.add(name)
    for name in base.keys()-field_names:
        require(np.array_equal(base[name], trial[name]), 'Field-plane event coordinate differs: '+name)
    times = column(base, 'step', 'state_step')
    result = {'all': {}, 'before_first_release': {}, 'after_first_release': {}, 'by_step': []}
    for label, mask in [('all', np.ones(len(times),dtype=bool)),
                        ('before_first_release', times < first_release),
                        ('after_first_release', times >= first_release)]:
        for kind in ('E','B'):
            a,b = field_vector(base,kind),field_vector(trial,kind)
            result[label][kind] = norm_metrics(a[mask],b[mask])
            result[label][kind]['unit'] = 'V/m' if kind == 'E' else 'T'
    for step in np.unique(times):
        mask = times == step
        result['by_step'].append({'step': float(step),
            **{kind:norm_metrics(field_vector(base,kind)[mask],field_vector(trial,kind)[mask])
               for kind in ('E','B')}})
    return result


def value(metrics, *names):
    for name in names:
        if name in metrics:
            return metrics[name]
    raise RuntimeError('Missing metric '+repr(names))


def production_checks(metrics, production, geometry_info, production_log):
    """Use frozen production files rather than a copied expected-value table."""
    mesh, boundary = production['mesh'], production['boundary']
    checks = {
        'boost_gamma': value(metrics,'boost_gamma','gamma_boost') == mesh['boost_gamma'],
        'boost_beta': abs(value(metrics,'boost_beta','beta_boost','window_beta')-geometry_info['beta_boost']) <= 2e-15,
        'cell_size_m': [metrics[k] for k in ('dx_m','dy_m','dz_m')] == mesh['cell_size'],
        'transverse_cells': [metrics['nx'],metrics['ny']] == mesh['cells'][:2],
        'c_dt_equals_dz': abs(metrics['dt_s']*float(C)/metrics['dz_m']-1) <= 2e-15,
        'cpml_cells': list(value(metrics,'cpml_cells')) == boundary['cells'][:2]+[0],
        'cpml_order': metrics['cpml_order'] == boundary['polynomial_order'],
        'cpml_kappa': metrics['cpml_kappa'] == boundary['kappa_max'],
        'cpml_alpha': metrics['cpml_alpha'] == boundary['alpha_fraction'],
        'cpml_reflection': metrics['cpml_reflection'] == boundary['target_reflection'],
        'lower_guard': value(metrics,'lower_guard_cells') == mesh['longitudinal']['lower_guard_cells'],
        'upper_guard': value(metrics,'upper_guard_cells') == mesh['longitudinal']['upper_guard_cells'],
        'initial_self_field_model': value(metrics,'initial_self_field_model') == production['initial_self_field']['model'],
        'initial_tolerance': metrics['initial_tolerance'] == production['initial_self_field']['relative_tolerance'],
        'initial_maximum_iterations': metrics['initial_maximum_iterations'] == production['initial_self_field']['maximum_iterations'],
        'initial_maximum_mean_transverse_beta': metrics['initial_maximum_mean_transverse_beta'] == production['initial_self_field']['maximum_mean_transverse_beta'],
        'activation_plane_lab_m': metrics['activation_plane_lab_m'] == production['beam']['reference']['activation_plane_z'],
        'true_detector_lab_m': metrics['production_detector_lab_m'] == production['detectors']['field_planes'][0]['z'],
        'production_detector_not_claimed': metrics['production_detector_reached'] is False,
        'one_cut_per_step': metrics['shift_count'] == metrics['steps'],
    }
    magnet = production['sources']['magnetic_elements'][0]
    for output,input_key in [('undulator_K','strength_parameter'),('undulator_period_m','period'),
                            ('undulator_periods','periods'),('undulator_entrance_lab_m','entrance_z'),
                            ('undulator_polarization_rad','polarization_angle_rad'),
                            ('undulator_gaussian_fringe','gaussian_fringe'),
                            ('undulator_fringe_cutoff','fringe_relative_cutoff')]:
        checks[output] = metrics[output] == magnet[input_key]
    maximum_advance = metrics['boost_gamma']*metrics['dz_m']*(metrics['window_beta']+metrics['particle_beta_max'])
    expected_substeps = int(np.ceil(mesh['particle_steps_per_undulator_period']*maximum_advance/magnet['period']))
    checks['particle_substeps_derived_from_production'] = metrics['particle_substeps'] == expected_substeps
    match = re.search(r'bunch RMS cells=\(([^)]+)\)',production_log)
    charge = re.search(r'total charge=([-+0-9.eE]+) C',production_log)
    require(match and charge, 'Production log lacks initial RMS or charge')
    rms = [float(x.strip()) for x in match.group(1).split(',')]
    checks['transverse_rms_from_production_log'] = metrics['sigma_position_cells'][:2] == rms[:2]
    checks['original_longitudinal_rms_from_production_log'] = metrics['source_sigma_z_cells'] == rms[2]
    scaled_charge = float(charge.group(1))*metrics['sigma_position_cells'][2]/rms[2]
    checks['line_density_preserving_charge_scale'] = abs(metrics['charge_C']/scaled_charge-1) <= 4e-15
    require(all(checks.values()), 'Production-parameter mismatch: '+repr(checks))
    return checks


def load_world(directory, metrics):
    result = {'name': directory.parent.name+'/'+directory.name, 'path': str(directory),
              'metrics': metrics, 'tables': {}}
    for name in ('initial-particles','release-particles','detector-crossings','final-particles','history','field-plane'):
        path = directory/(name+'.csv')
        require(path.is_file(), 'Missing required scientific output: '+str(path))
        t = table(path)
        if 'id' in t:
            t = order_by_id(t)
        result['tables'][name] = t
    require(size(result['tables']['initial-particles']) == size(result['tables']['final-particles']),
            'Final cohort incomplete')
    return result


def cost(metrics):
    # Return every named elapsed-time measurement without assuming nested
    # timers are additive. Probe documentation defines their scope.
    return {k:v for k,v in metrics.items() if k.endswith(('_wall_s','_seconds','_time_s'))
            and isinstance(v,(int,float))}


def inlet_review(path, metrics, expected_steps, first_release):
    with path.open(newline='') as stream:
        reader=csv.DictReader(stream)
        rows=list(reader)
    require(len(rows) == len(expected_steps)*36, 'Incomplete inlet diagnostic series')
    numeric=[]
    for row in rows:
        require(None not in row and row['scope'] in ('all_xy','physical_xy'), 'Invalid inlet row')
        r={name:(text if name == 'scope' else np.longdouble(text)) for name,text in row.items()}
        require(all(np.isfinite(x) for name,x in r.items() if name != 'scope'), 'Nonfinite inlet diagnostics')
        n,k,a=(int(r[name]) for name in ('step','relative_k','component'))
        require(n in expected_steps and k in (-1,0,1) and 0 <= a < 6, 'Unknown inlet location/time/component')
        require(r['global_k'] == metrics['reference_offset0_cells']-n+k, 'Inlet global native coordinate mismatch')
        target=(np.longdouble(n)-(np.longdouble('.5') if a >= 3 else 0))*np.longdouble(str(metrics['dt_s']))
        require(abs(r['native_time_box_s']-target) <= 8*np.finfo(float).eps*max(abs(target),np.longdouble(str(metrics['dt_s']))),
                'Inlet E/B native time mismatch')
        require(r['candidate_present'] == int(k >= 0), 'Inlet candidate availability mismatch')
        require(r['samples'] > 0 and r['samples'] == int(r['samples']), 'Invalid inlet sample count')
        if k == 0:
            require(r['candidate_reference_difference_L2'] == 0,
                    'New inlet slice does not exactly equal supplied native reference')
        numeric.append(r)
    keys={(int(r['step']),int(r['relative_k']),int(r['component']),r['scope']) for r in numeric}
    require(len(keys) == len(numeric), 'Duplicate inlet diagnostic row')
    aggregate=[]
    for n in expected_steps:
        for k in (-1,0,1):
            for scope in ('all_xy','physical_xy'):
                item={'step':int(n),'relative_k':k,'scope':scope}
                for kind,components in [('E',range(3)),('B',range(3,6))]:
                    selected=[r for r in numeric if r['step'] == n and r['relative_k'] == k and
                              r['scope'] == scope and r['component'] in components]
                    require(len(selected) == 3, 'Incomplete inlet field vector')
                    denominator=sum(r['full_L2']**2 for r in selected)
                    difference=sum(r['reference_full_difference_L2']**2 for r in selected)
                    samples=sum(r['samples'] for r in selected)
                    item[kind]={'reference_full_relative_l2':float(np.sqrt(difference/denominator)) if denominator else None,
                        'reference_full_component_rms_difference':float(np.sqrt(difference/samples)),
                        'reference_full_max_abs_difference':float(max(r['reference_full_difference_max'] for r in selected)),
                        'reference_full_difference_l2':float(np.sqrt(difference)),
                        'full_l2':float(np.sqrt(denominator)),
                        'candidate_full_difference_l2':float(np.sqrt(sum(r['candidate_full_difference_L2']**2 for r in selected))) if k >= 0 else None,
                        'candidate_reference_difference_l2':float(np.sqrt(sum(r['candidate_reference_difference_L2']**2 for r in selected))) if k >= 0 else None}
                aggregate.append(item)
    onset={}; envelope={}
    for kind in ('E','B'):
        for k in (-1,0,1):
            key=kind+'_physical_xy_k'+str(k)
            changed=[a['step'] for a in aggregate if a['relative_k'] == k and a['scope'] == 'physical_xy'
                     and a[kind]['reference_full_difference_l2'] > 0]
            onset[key]=min(changed) if changed else None
            selected=[a for a in aggregate if a['relative_k'] == k and a['scope'] == 'physical_xy']
            envelope[key]={}
            for label,items in [('all',selected),('before_first_release',[a for a in selected if a['step'] < first_release])]:
                valid=[a for a in items if a[kind]['reference_full_relative_l2'] is not None]
                envelope[key][label]={'maximum_relative_l2':max((a[kind]['reference_full_relative_l2'] for a in valid),default=None),
                                     'maximum_component_rms_difference':max((a[kind]['reference_full_component_rms_difference'] for a in valid),default=None)}
    return {'rows':len(rows),'new_slice_exactly_equals_supplied_reference':True,
            'first_sampled_nonzero_reference_difference_step':onset,
            'onset_definition':'Exact stored norm > 0; sampling cadence limits onset resolution and this is not a physical tolerance.',
            'component_norm_scope':'Aggregate of raw native component norms, with stagger-dependent sample counts; not a field energy.',
            'physical_xy_envelope':envelope,
            'by_step_location':aggregate}


def guard_audit(metrics, checks_path):
    checks=json.loads(checks_path.read_text())
    require(checks.get('pass') is True and checks['format'] == 'aprl-held-cut-causality-checks-v1',
            'Causal-envelope experiment failed or has wrong provenance')
    n=checks['steps_per_side']
    require(checks['lower_maximum_native_k'] == list(range(2,n+2)) and
            checks['upper_maximum_native_depth'] == list(range(1,n+1)),
            'Observed causal front differs from documented native bound')
    steps,warm,offset,width,total=(metrics[k] for k in
                                 ('steps','held_prewarm_steps','reference_offset0_cells','cut_nz','full_nz'))
    lower=offset-(2*steps+warm+2)
    upper=total-(offset+width)-(steps+warm+1)
    require(lower > 0 and upper > 0,'Remote endpoint reaches required inlet/retained data')
    require(lower == metrics['reference_causal_lower_margin_cells'] and
            upper == metrics['reference_causal_upper_margin_cells'], 'Producer causal-margin accounting differs')
    return {'lower_margin_native_cells':lower,'upper_margin_native_cells':upper,
            'minimum_requested_lower_global_k':offset-steps-1,
            'maximum_lower_endpoint_influence_native_k':steps+warm+1,
            'independent_check_file':str(checks_path),'check_sha256':digest(checks_path),
            'checked_field_values':checks['field_values_checked'],'checked_history_values':checks['history_values_checked'],
            'lower_dependency_recurrence':'Bt=max(Bt,Et,Ez+1); Bz=max(Bz,Et+1); then Et=max(Et,Bz,Bt+1), Ez=max(Ez,Bt). From initial endpoint index0 this bounds all components by n+1.',
            'upper_dependency_recurrence':'Each magnetic target reads at most the next-z electric plane; electric update adds no further next-z dependency. The upper front retreats at most one native index per complete step.',
            'scope':'Dynamic endpoint influence with fixed prescribed sources; excludes initial elliptic-boundary independence and full particle-feedback envelope.'}


def prefix_review(directory, previous, final_step):
    checks={}
    for role in ('full','cut'):
        for name in ('initial-particles','release-particles','detector-crossings','field-plane','history'):
            current=table(directory/role/(name+'.csv'))
            prior=table(previous/role/(name+'.csv'))
            event_time=column(current,'step')
            if name in ('release-particles','detector-crossings'):
                event_time=event_time+column(current,'fraction')
            current={k:v[event_time <= final_step] for k,v in current.items()}
            if 'id' in current:
                current,prior=order_by_id(current),order_by_id(prior)
            omitted=('continuity_max','charge_balance_max_C') if name == 'history' else ()
            current={k:v for k,v in current.items() if k not in omitted}
            prior={k:v for k,v in prior.items() if k not in omitted}
            same_tables(current,prior,role+'/'+name+' replay prefix')
            checks[role+'/'+name]=True
    return {'previous_case':str(previous),'through_step':final_step,'checks':checks,
            'omitted_history_columns':['continuity_max','charge_balance_max_C'],
            'reason':'The earlier endpoint requested an extra continuity audit at step64; physical history columns remain exactly compared. Final-particle files are at different times and are not replay inputs.'}


def plot(root, worlds, comparisons):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes = plt.subplots(2,2,figsize=(11,8),constrained_layout=True)
    names=list(dict.fromkeys(case['name'].split('/')[0] for case in worlds))
    palette=dict(zip(names,plt.rcParams['axes.prop_cycle'].by_key()['color']))
    for case in worlds:
        t = case['tables']['history']
        x = column(t,'step','state_step')
        name,role=case['name'].split('/')
        for axis, aliases in [(axes[0,0],('mean_delta_energy_lab_eV','mean_delta_lab_eV')),
                              (axes[0,1],('rms_x_m','sigma_x_m'))]:
            try:
                y = column(t,*aliases)
            except RuntimeError:
                continue
            axis.plot(x.astype(float), y.astype(float),label=case['name'],lw=1.3,
                      color=palette[name],linestyle='-' if role == 'full' else '--')
    for comparison in comparisons:
        points = comparison['fields']['by_step']
        for axis,kind in [(axes[1,0],'E'),(axes[1,1],'B')]:
            axis.plot([p['step'] for p in points],
                      [p[kind]['vector_rms_difference'] for p in points],
                      label=comparison['case'],lw=1.3,color=palette[comparison['case'].split('/')[0]])
    axes[0,0].set(ylabel='Mean lab energy change [eV / electron]')
    axes[0,1].set(ylabel='Centered transverse RMS x [m]')
    axes[1,0].set(ylabel='Plane vector RMS E difference [V/m]')
    axes[1,1].set(ylabel='Plane vector RMS B difference [T]')
    for axis in axes.flat:
        axis.set_xlabel('Boost-frame step')
        maximum=max(int(w['metrics']['steps']) for w in worlds)
        ticks=sorted(set(range(0,maximum+1,20))|{maximum})
        axis.set_xticks(ticks)
        axis.grid(alpha=.25)
        handles,_ = axis.get_legend_handles_labels()
        if handles:
            axis.legend(fontsize=8)
    all_events=all(w['metrics']['released'] == w['metrics']['auxiliary_detected'] == w['metrics']['particles'] for w in worlds)
    cohorts=sorted(set(w['metrics']['particles'] for w in worlds))
    description=(str(cohorts[0])+' particles, all released and detected at auxiliary plane') if len(cohorts)==1 and all_events else 'Partial release / auxiliary-plane cohort'
    fig.suptitle('Held inlet / right cut: measured small-case comparison\n'+description+'; not the 4.5 m detector')
    fig.savefig(root/'comparison.png',dpi=180)
    fig.savefig(root/'comparison.pdf')
    plt.close(fig)


def plot_inlets(root, comparisons, first_releases):
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(2,2,figsize=(11,7),layout='constrained')
    colors=plt.rcParams['axes.prop_cycle'].by_key()['color']
    for color,comparison in zip(colors,comparisons):
        if 'inlet' not in comparison:
            continue
        name=comparison['case'].split('/')[0]
        records=comparison['inlet']['by_step_location']
        for column_index,kind in enumerate(('E','B')):
            for k,style in [(-1,'-'),(0,'--')]:
                points=[a for a in records if a['relative_k'] == k and a['scope'] == 'physical_xy']
                axes[0,column_index].plot([p['step'] for p in points],
                    [p[kind]['reference_full_relative_l2'] for p in points],
                    color=color,linestyle=style,label=name+', k='+str(k))
            points=[a for a in records if a['relative_k'] == 1 and a['scope'] == 'physical_xy']
            axes[1,column_index].plot([p['step'] for p in points],
                [p[kind]['candidate_reference_difference_l2']/p[kind]['full_l2'] if p[kind]['full_l2'] else np.nan for p in points],
                color=color,label=name)
            axes[0,column_index].set_ylabel(kind+' reference vs full, relative L2')
            axes[1,column_index].set_ylabel(kind+' candidate vs reference, relative to full L2')
    first=min(first_releases.values())
    maximum=max(p['step'] for c in comparisons for p in c['inlet']['by_step_location'])
    for axis in axes.flat:
        axis.axvline(first,color='.45',linestyle=':',linewidth=1,label='First particle release')
        axis.set_xlabel('Boost-frame step')
        axis.set_xticks(sorted(set(range(0,int(maximum)+1,20))|{int(maximum)}))
        axis.grid(alpha=.25);axis.legend(fontsize=7)
    axes[0,0].set_title('Inlet: reference vs full (k=-1, 0)')
    axes[0,1].set_title('Inlet: reference vs full (k=-1, 0)')
    axes[1,0].set_title('Inside inlet after evolution (k=1)')
    axes[1,1].set_title('Inside inlet after evolution (k=1)')
    fig.suptitle('Native left-inlet diagnostics, physical transverse interior\n'
                 'Read-only full-state comparison; E(n), B(n-1/2); raw norms, not energy')
    fig.savefig(root/'inlet-comparison.png',dpi=180)
    fig.savefig(root/'inlet-comparison.pdf')
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,required=True)
    parser.add_argument('--no-plot',action='store_true')
    parser.add_argument('--previous-campaign',type=Path)
    parser.add_argument('--control-campaign',type=Path)
    parser.add_argument('--causality-checks',type=Path)
    args = parser.parse_args()
    root = args.root.resolve()
    require((root/'COMPLETED.json').is_file(), 'Campaign has not completed')
    plan = json.loads((root/'plan.json').read_text())
    require(plan['cases'], 'No scientific cases')
    production_file = root/'production-run.yaml'
    geometry_file = root/'production-geometry-preflight.json'
    production = yaml.safe_load(production_file.read_text())
    geometry_info = json.loads(geometry_file.read_text())
    log_file = root/'production-aprl.log'
    if not log_file.is_file():
        log_file = PRODUCTION/'aprl.log'
    production_log = log_file.read_text()
    checked, comparisons, worlds, first_releases, reference_moments = [], [], [], {}, {}
    full_replays, groups, prefix_replays, preparation = [], {}, [], []
    files = [Path(__file__).resolve(),root/'plan.json',root/'COMPLETED.json',production_file,geometry_file,log_file]
    files += [root/name for name in ('source-files.json','frozen.json','ANALYSIS_REVISION.json') if (root/name).is_file()]
    previous_cases=[]; control_cases=[]; external_control_replays=[]
    if args.previous_campaign:
        previous_root=args.previous_campaign.resolve()
        require((previous_root/'COMPLETED.json').is_file(), 'Previous campaign incomplete')
        previous_plan=json.loads((previous_root/'plan.json').read_text())
        previous_cases=[(previous_root/s['name'],json.loads((previous_root/s['name']/'metrics.json').read_text()))
                        for s in previous_plan['cases']]
    if args.control_campaign:
        control_root=args.control_campaign.resolve()
        require((control_root/'COMPLETED.json').is_file(), 'Control campaign incomplete')
        control_plan=json.loads((control_root/'plan.json').read_text())
        control_cases=[(control_root/s['name'],json.loads((control_root/s['name']/'metrics.json').read_text()))
                       for s in control_plan['cases']]
    for specification in plan['cases']:
        directory = root/specification['name']
        m = json.loads((directory/'metrics.json').read_text())
        require(m.get('completed') is True, 'Case incomplete: '+str(directory))
        require(m['mode'] in ('live','frozen'), 'Unrecognized held reference mode')
        require(m['boost_gamma'] == production['mesh']['boost_gamma'] and
                m['window_beta'] == geometry_info['beta_boost'], 'Recorded boost parameters differ from production')
        boost = np.longdouble(np.float64(m['boost_gamma']))
        beta = np.longdouble(np.float64(m['window_beta']))
        moment_path = WORK/m['source_moment_sample_path']
        require(digest(moment_path) == m['source_moment_sample_sha256'], 'Source moment sample changed')
        moment = json.loads(moment_path.read_text())
        require(m['source_moment_sample_records'] == moment['sample_records'], 'Source sample size differs')
        require(m['lab_u_mean_parameter'] == moment['lab_u_mean'] and
                m['lab_u_sigma_parameter'] == moment['lab_u_sigma'], 'Source momentum parameters differ from sample')
        require(m['transverse_position_u_correlation_parameter'] ==
                [moment['x_ux_correlation'],moment['y_uy_correlation']],
                'Source transverse correlations differ from sample')
        files.append(moment_path)
        base,case = (load_world(directory/name,m) for name in ('full','cut'))
        worlds.extend((base,case))
        group = tuple(m[name] for name in ('steps','particles','held_prewarm_steps','full_nz','source_moment_sample_sha256'))
        if group in groups:
            previous = groups[group]
            matching = {}
            for csv_file in sorted((directory/'full').glob('*.csv')):
                other = previous/'full'/csv_file.name
                require(other.is_file() and csv_file.read_bytes() == other.read_bytes(),
                        'Full PIC control changed with inlet reference mode: '+str(csv_file))
                matching[csv_file.name] = digest(csv_file)
            full_replays.append({'first_case':previous.name,'second_case':directory.name,
                                 'all_full_control_csv_byte_identical':True,'sha256':matching})
        else:
            groups[group] = directory
        checked.append({'case':directory.name,'production_parameters':production_checks(m,production,geometry_info,production_log),
                        'cost_measurements':cost(m)})
        bt = base['tables']; initial = bt['initial-particles']
        require(size(initial) == m['particles'] > 0, 'Initial cohort count differs')
        release_steps = column(bt['release-particles'],'step')+column(bt['release-particles'],'fraction')
        require(len(release_steps) == m['released'] > 0, 'No released actual PIC particles')
        require(size(bt['detector-crossings']) == m['auxiliary_detected'], 'Auxiliary event count differs')
        first_release = np.min(release_steps)
        first_releases[directory.name] = float(first_release)
        reference_moments[directory.name] = particle_summary(bt['final-particles'],initial,boost,beta)
        require(m['cut_admissions'] > m['initial_cut_active'], 'No source actually entered the left boundary')
        require(m['minimum_cut_pic_clearance_m'] >= 2*m['dz_m'], 'PIC sampled sacrificial boundary region')
        require(m['cut_continuity_max'] < 1e-10 and m['full_continuity_max'] < 1e-10,
                'Existing producer charge-continuity gate failed')
        require(m['cut_discards'] == 0, 'Review of post-discard final particle states is not implemented')
        ct = case['tables']
        same_tables(initial,ct['initial-particles'],'Initial particles')
        same_tables(bt['release-particles'],ct['release-particles'],'Prescribed release events')
        for event in ('detector-crossings','final-particles'):
            require(np.array_equal(column(bt[event],'id'),column(ct[event],'id')),event+' cohorts differ')
        history_steps = column(bt['history'],'step')
        require(np.array_equal(history_steps,column(ct['history'],'step')) and history_steps[-1] == m['steps'],
                'History times or final step differ')
        require(len(np.unique(history_steps)) == len(history_steps), 'Duplicate history step')
        expected_steps = np.asarray(sorted(set(range(0,m['steps']+1,4))|{m['steps']}),dtype=np.longdouble)
        require(np.array_equal(history_steps,expected_steps), 'Incomplete history cadence')
        for role,tables in [('full',bt),('cut',ct)]:
            field=tables['field-plane']
            keys=list(zip(column(field,'step'),column(field,'pixel_i'),column(field,'pixel_j')))
            require(len(set(keys)) == len(keys), 'Duplicate field-plane sample in '+role)
            nx=len(range(m['cpml_cells'][0]+2,m['nx']-m['cpml_cells'][0]-2,4))
            ny=len(range(m['cpml_cells'][1]+2,m['ny']-m['cpml_cells'][1]-2,4))
            require(size(field) == len(expected_steps)*nx*ny, 'Incomplete field-plane image series')
            for step in expected_steps:
                require(np.count_nonzero(column(field,'step') == step) == nx*ny,
                        'Incomplete transverse plane at step '+str(step))
            require(np.all(column(tables['final-particles'],'step') == m['steps']),
                    'Final particle records have inconsistent times')
        overlap = table(directory/'overlap.csv')
        require(np.array_equal(column(overlap,'step'),history_steps), 'Overlap/history sampling differs')
        active_mask = column(bt['final-particles'],'frozen') == 0
        require(np.array_equal(active_mask,column(ct['final-particles'],'frozen') == 0),
                'Full/cut released cohort differs at final time')
        released_full = {name:values[active_mask] for name,values in bt['final-particles'].items()}
        released_cut = {name:values[active_mask] for name,values in ct['final-particles'].items()}
        comparisons.append({'case':case['name'],'reference_case':base['name'],'reference_mode':m['mode'],
            'initial_and_prescribed_release_identical':True,
            'fields':fields_comparison(bt['field-plane'],ct['field-plane'],first_release),
            'particles':{event:particle_comparison(bt[event],ct[event],initial,boost,beta)
                         for event in ('detector-crossings','final-particles')},
            'final_released_cohort_only':particle_comparison(released_full,released_cut,initial,boost,beta),
            'overlap':{name:[float(x) for x in values] for name,values in overlap.items()},
            'cost_measurements':cost(m),
            'candidate_online_including_reference_s':m['cut_preparation_wall_s']+m['cut_advance_wall_s']+
                m['cut_remap_and_feed_wall_s']+m['live_reference_wall_s'],
            'candidate_plus_initial_preparation_s':m['initialization_wall_s']+m['held_prewarm_wall_s']+
                m['cut_preparation_wall_s']+m['cut_advance_wall_s']+m['cut_remap_and_feed_wall_s']+m['live_reference_wall_s'],
            'full_plus_initial_preparation_s':m['initialization_wall_s']+m['held_prewarm_wall_s']+
                m['full_preparation_wall_s']+m['full_advance_wall_s'],
            'cost_sum_scope':'Nonoverlapping measured preparation/advance/feed components; excludes shared diagnostics, uninstrumented construction/output and process overhead. Overall campaign wall includes both full and cut, so it is not a single-method benchmark.'})
        if 'reference_offset0_cells' in m:
            require((directory/'inlet.csv').is_file(), 'Missing native inlet diagnostic')
            comparisons[-1]['inlet']=inlet_review(directory/'inlet.csv',m,expected_steps,first_release)
            checks_path=args.causality_checks.resolve() if args.causality_checks else root.parent/'causality-r1/checks/metrics.json'
            comparisons[-1]['causal_guard_audit']=guard_audit(m,checks_path)
            files.extend((directory/'inlet.csv',checks_path))
            preparation.append({'case':directory.name,'warm_steps':m['held_prewarm_steps'],
                'prepared_native_component_L2':m['prepared_reference_native_component_L2'],
                'prepared_native_component_maximum':m['prepared_reference_native_component_maximum'],
                'prepared_native_history_L2':m['prepared_reference_history_L2'],
                'plane_E_relative_l2':comparisons[-1]['fields']['all']['E']['relative_l2'],
                'plane_B_relative_l2':comparisons[-1]['fields']['all']['B']['relative_l2'],
                'particle_rms_energy_difference_eV':comparisons[-1]['particles']['final-particles']['rms_energy_difference_eV'],
                'inlet_physical_xy_envelope':comparisons[-1]['inlet']['physical_xy_envelope']})
        for previous_directory,previous_metrics in previous_cases:
            if previous_metrics['steps'] < m['steps'] and all(previous_metrics[k] == m[k] for k in
                    ('mode','held_prewarm_steps','particles','full_nz','cut_nz','source_moment_sample_sha256')):
                prefix_replays.append({'case':directory.name,**prefix_review(directory,previous_directory,previous_metrics['steps'])})
                files.append(previous_directory/'metrics.json')
                for role in ('full','cut'):
                    files.extend((previous_directory/role).glob('*.csv'))
        for control_directory,control_metrics in control_cases:
            if all(control_metrics[k] == m[k] for k in
                    ('steps','held_prewarm_steps','particles','full_nz','cut_nz','source_moment_sample_sha256')):
                matched=[]
                for csv_file in sorted((directory/'full').glob('*.csv')):
                    other=control_directory/'full'/csv_file.name
                    require(other.is_file() and csv_file.read_bytes() == other.read_bytes(),
                            'External full control differs: '+str(other))
                    matched.append(csv_file.name);files.append(other)
                files.append(control_directory/'metrics.json')
                external_control_replays.append({'case':directory.name,'control_case':str(control_directory),
                                                'all_full_control_csv_byte_identical':True,'files':matched})
        files += [directory/'metrics.json',directory/'overlap.csv']
        for role in ('full','cut'):
            files += list((directory/role).glob('*.csv'))
    require(not args.control_campaign or external_control_replays, 'No matching requested external full control')
    report = {
        'format':'aprl-held-cut-independent-review-v1','campaign':str(root),
        'structural_status':'PASS','physical_acceptance':'NOT_ASSIGNED_METRICS_ONLY',
        'numeric_precision':{
            'particle_momentum':'ux/uy/uz were emitted from C++ double with 17 significant digits; round-trip through float64, then promote to longdouble for stable gamma/energy calculations.',
            'boost_parameters':'boost_gamma/window_beta from each probe metrics.json were emitted from original double parameters. Their values are checked against production-run.yaml / production-geometry-preflight.json, recovered as float64, then promoted to longdouble without decimal-string reinterpretation.',
            'reported_energy_difference':'delta_energy_box_eV and delta_energy_lab_eV were computed as C++ long double; CSV columns remain directly parsed as numpy.longdouble, without float64 conversion.',
            'si_constants':'Electron mass and elementary charge retain independent decimal long-double parsing; their slight difference from C++ double constants and printed-value rounding can leave a small residual against producer delta-energy records.',
            'remaining_csv':'Other columns and field norms retain the recorded decimal values parsed as numpy.longdouble.',
            'source_types':['src/fieldvector.h: Double is double','src/eb_particles.h: properVelocity is FieldVector<Double>','verification/held_cut_bunch_probe.cpp: fresh uses setprecision(17); particleRow accumulates energy differences in long double']},
        'production_parameter_checks':checked,
        'reference_final_particles':reference_moments,
        'first_release_step':first_releases,
        'comparisons':comparisons,
        'matched_full_control_replays':full_replays,
        'external_full_control_replays':external_control_replays,
        'earlier_campaign_physical_prefix_replays':prefix_replays,
        'held_preparation_comparison':{'cases':preparation,
            'scope':'Each frozen candidate is compared with its own paired full PIC state. Changing prewarm changes the prepared reference and full initial fields; their equality across prewarms is not assumed. Raw prepared component/history norms are not an energy.'},
        'scope':['Generated shortened bunch and finite-time interface comparison only.',
                 'Nearby field plane is auxiliary; the production 4.5 m detector is not validated.',
                 'Reference donor/preparation costs must be counted; no production speedup or complexity claim.',
                 'Field norms use matched native plane events, not total energies of different domain sizes.',
                 'The auxiliary plane samples the transverse interior every four grid cells; it is not an integrated full-plane radiation-power estimate.',
                 'Per-particle energy differences use extended precision and rationalized gamma differences.',
                 'Native B is a half-step behind E at the same recorded z; these are not simultaneous lab E/B or a Poynting diagnostic.',
                 'Synthetic Gaussian source preserves production transverse RMS, scaled line density and sampled momentum/correlation moments; it is not the original full 6D ELEGANT distribution.',
                 'Net transverse size can focus because of the source correlations; agreement does not separately isolate space-charge expansion from zero-self-field dynamics.',
                 'Initial cases exercise incoming particles, not particle exits through the right sacrificial collar.',
                 'No numerical physical-pass threshold is introduced by this audit.'],
        'input_sha256':{str(p):digest(p) for p in sorted(set(files))},
    }
    save(root/'REVIEW.json',report)
    if not args.no_plot:
        plot(root,worlds,comparisons)
        if all('inlet' in c for c in comparisons):
            plot_inlets(root,comparisons,first_releases)
    print(json.dumps({'structural_status':'PASS','comparisons':len(comparisons),
                      'review':str(root/'REVIEW.json')},sort_keys=True))


if __name__ == '__main__':
    main()
