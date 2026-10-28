#!/usr/bin/env python3
"""Audit a matched CG/blocked Cowan-reference pair without running simulations.

Solver replacement and the existing frozen-inlet approximation are separate
comparisons. Only structural and producer residual contracts are pass/fail;
no field or weak-energy physical tolerance is invented by this reviewer.
"""
import argparse
import csv
import json
from pathlib import Path

import numpy as np

import review as audit
import convergence_review as convergence


EVENTS = ('initial-particles', 'release-particles', 'detector-crossings',
          'final-particles', 'history', 'field-plane')
MATCH_KEYS = (
    'steps', 'particles', 'nx', 'ny', 'full_nz', 'cut_nz',
    'dx_m', 'dy_m', 'dz_m', 'dt_s', 'reference_offset0_cells',
    'boost_gamma', 'window_beta', 'particle_beta_mean', 'particle_beta_min',
    'particle_beta_max', 'sigma_position_cells', 'charge_C',
    'cpml_cells', 'cpml_order', 'cpml_kappa', 'cpml_alpha', 'cpml_reflection',
    'lower_guard_cells', 'upper_guard_cells', 'particle_substeps',
    'initial_tolerance', 'initial_maximum_mean_transverse_beta',
    'undulator_K', 'undulator_period_m', 'undulator_periods',
    'undulator_entrance_lab_m', 'undulator_polarization_rad',
    'undulator_gaussian_fringe', 'undulator_fringe_cutoff',
    'activation_plane_lab_m', 'auxiliary_plane_lab_m',
    'production_detector_lab_m', 'source_moment_sample_sha256',
    'lab_u_mean_parameter', 'lab_u_sigma_parameter',
    'transverse_position_u_correlation_parameter',
)


def read_json(path):
    return json.loads(Path(path).read_text())


def provenance(directory, metrics, files):
    """Prefer campaign-frozen production inputs; never silently invent constants."""
    campaign = directory.parent
    if (campaign/'production-run.yaml').is_file():
        paths = [campaign/'production-run.yaml',
                 campaign/'production-geometry-preflight.json',
                 campaign/'production-aprl.log']
    else:
        paths = [audit.PRODUCTION/'run.yaml',
                 audit.PRODUCTION/'geometry-preflight.json',
                 audit.PRODUCTION/'aprl.log']
    audit.require(all(p.is_file() for p in paths), 'Missing production provenance')
    checks = audit.production_checks(metrics, audit.yaml.safe_load(paths[0].read_text()),
                                     read_json(paths[1]), paths[2].read_text(),
                                     metrics['initial_self_field_model'])
    sample_path = audit.WORK/metrics['source_moment_sample_path']
    audit.require(audit.digest(sample_path) == metrics['source_moment_sample_sha256'],
                  'Source-moment sample changed')
    sample = read_json(sample_path)
    audit.require(metrics['lab_u_mean_parameter'] == sample['lab_u_mean'] and
                  metrics['lab_u_sigma_parameter'] == sample['lab_u_sigma'] and
                  metrics['transverse_position_u_correlation_parameter'] ==
                  [sample['x_ux_correlation'], sample['y_uy_correlation']],
                  'Source-moment parameters differ from the input sample')
    files.extend(paths+[sample_path])
    return checks


def load_case(directory, label, expected_backend, files):
    directory = directory.resolve()
    m = read_json(directory/'metrics.json')
    backend = m.get('initial_backend', 'cg')  # Historical CG fixtures predate the key.
    audit.require(m.get('completed') is True, 'Incomplete probe: '+str(directory))
    audit.require(expected_backend is None or backend == expected_backend,
                  'Unexpected initial backend: '+str(directory))
    audit.require(m['mode'] == 'frozen', 'This audit expects the frozen-reference experiment')
    if expected_backend is not None:
        audit.require(m['initial_self_field_model'] == 'cowan-relativistic-poisson' and
                      m['held_prewarm_steps'] == 0,
                      'Solver pair must use the same Cowan0 reference')
    audit.require(m['released'] == m['auxiliary_detected'] == m['particles'] > 0,
                  'All particles must be released and observed at the auxiliary plane')
    audit.require(m['cut_discards'] == 0, 'Post-discard cohort analysis is not implemented')
    audit.require(m['minimum_cut_pic_clearance_m'] >= 2*m['dz_m'],
                  'Particle sampled the sacrificial boundary region')
    audit.require(m['cut_continuity_max'] < 1e-10 and m['full_continuity_max'] < 1e-10,
                  'Existing charge-continuity gate failed')
    residual = m['initial_relative_residual']
    audit.require(np.isfinite(residual) and 0 <= residual <= m['initial_tolerance'],
                  'Initial field failed its recorded relative-residual contract')
    worlds = {role: audit.load_world(directory/role, m)['tables'] for role in ('full','cut')}
    first = worlds['full']['initial-particles']
    ids = first['id']
    expected_steps = np.asarray(sorted(set(range(0,m['steps']+1,4))|{m['steps']}),
                                dtype=np.longdouble)
    pixel_i = list(range(m['cpml_cells'][0]+2,m['nx']-m['cpml_cells'][0]-2,4))
    pixel_j = list(range(m['cpml_cells'][1]+2,m['ny']-m['cpml_cells'][1]-2,4))
    for role, tables in worlds.items():
        for event in EVENTS[:4]:
            t = tables[event]
            audit.require(audit.size(t) == m['particles'] and np.array_equal(t['id'],ids),
                          role+'/'+event+' particle cohort differs')
            audit.require(np.array_equal(t['charge_C'],first['charge_C']),
                          role+'/'+event+' particle weights differ')
        audit.same_tables(first,tables['initial-particles'],'Matched initial source')
        audit.same_tables(worlds['full']['release-particles'],tables['release-particles'],
                          'Matched prescribed release events')
        audit.require(np.all(tables['final-particles']['step'] == m['steps']) and
                      np.all(tables['final-particles']['fraction'] == 0) and
                      np.all(tables['final-particles']['frozen'] == 0),
                      'Final particle time/state differs')
        audit.require(np.array_equal(tables['history']['step'],expected_steps),
                      'Incomplete history sampling')
        plane = tables['field-plane']
        for kind,half_step in (('E',0),('B',.5)):
            target=(plane['step']-half_step)*np.longdouble(m['dt_s'])
            bound=16*np.finfo(float).eps*np.maximum(np.abs(target),m['dt_s'])
            audit.require(np.all(np.abs(plane['time_'+kind+'_box_s']-target)<=bound),
                          'Plane E/B native time staggering differs')
        keys = set(zip(plane['step'],plane['pixel_i'],plane['pixel_j']))
        audit.require(len(keys) == audit.size(plane) == len(expected_steps)*len(pixel_i)*len(pixel_j),
                      'Duplicate or incomplete plane sampling')
        for step in expected_steps:
            selected = plane['step'] == step
            audit.require(np.count_nonzero(selected) == len(pixel_i)*len(pixel_j) and
                          np.array_equal(np.unique(plane['pixel_i'][selected]),pixel_i) and
                          np.array_equal(np.unique(plane['pixel_j'][selected]),pixel_j),
                          'Plane native pixel lattice differs')
        files.extend(directory/role/(event+'.csv') for event in EVENTS)
    releases = worlds['full']['release-particles']
    release_time = releases['step']+releases['fraction']
    prepared_audit = audit.prepared_review(directory/'prepared-reference.csv',m)
    inlet = audit.inlet_review(directory/'inlet.csv',m,expected_steps,np.min(release_time))
    files.extend([directory/'metrics.json',directory/'prepared-reference.csv',directory/'inlet.csv'])
    blocked_report = None
    if backend == 'blocked':
        report_name = m.get('blocked_reference_report_file')
        audit.require(isinstance(report_name,str) and report_name,
                      'Blocked backend did not record its independent report')
        report_path = (directory/report_name).resolve()
        audit.require(report_path.parent == directory, 'Blocked report is outside its case')
        blocked_report = read_json(report_path)
        nodes=[m[key]-2*m['cpml_cells'][axis]-1
               for axis,key in enumerate(('nx','ny','full_nz'))]
        audit.require(blocked_report.get('complete') is True and
                      blocked_report['particles']==m['particles'] and
                      blocked_report['interiorNodes']==nodes and
                      blocked_report['fieldTimes']=='E(0),B(-dt/2)' and
                      blocked_report['blockLayers']==min(m['reference_block_layers'],nodes[2]),
                      'Blocked report does not describe this completed reference')
        audit.require(blocked_report['beta']==m['reference_solver_beta'],
                      'Blocked report beta differs from the initializer')
        audit.require(m.get('reference_runtime_uses_full_donor') is False,
                      'Blocked candidate still uses the complete runtime reference donor')
        files.append(report_path)
    return {'name':label,'directory':str(directory),'metrics':m,'backend':backend,
            'worlds':worlds,'prepared':audit.table(directory/'prepared-reference.csv'),
            'first_release':float(np.min(release_time)),'last_release':float(np.max(release_time)),
            'production_checks':provenance(directory,m,files),
            'prepared_audit':prepared_audit,'inlet_audit':inlet,'blocked_report':blocked_report}


def matching_cases(a,b,allow_different_preparation=False):
    for key in MATCH_KEYS:
        audit.require(a['metrics'][key] == b['metrics'][key], 'Mismatched experiment setting: '+key)
    audit.require(a['metrics'].get('warm_horizon_steps',64) == b['metrics'].get('warm_horizon_steps',64),
                  'Reference geometry preparation horizon differs')
    if not allow_different_preparation:
        for key in ('initial_self_field_model','held_prewarm_steps','initial_maximum_iterations',
                    'compute_threads'):
            audit.require(a['metrics'][key] == b['metrics'][key], 'Mismatched solver-pair setting: '+key)
        for key in ('reference_solver_beta','reference_translation_beta'):
            if key in a['metrics'] or key in b['metrics']:
                audit.require(key in a['metrics'] and key in b['metrics'] and
                              a['metrics'][key]==b['metrics'][key],
                              'Mismatched solver-pair velocity: '+key)
    for event in ('initial-particles','release-particles'):
        audit.same_tables(a['worlds']['full'][event],b['worlds']['full'][event],
                          'Solver pair '+event)


def stable_energy_difference(reference,trial,metrics):
    """Per-ID energy difference without subtracting two ~600 MeV energies."""
    audit.require(np.array_equal(reference['id'],trial['id']),'Energy-difference IDs differ')
    boost,beta = convergence.precision_parameters(metrics)
    a,b = audit.particle_momentum(reference),audit.particle_momentum(trial)
    ga,gb = (np.sqrt(1+np.sum(v*v,axis=1,dtype=np.longdouble)) for v in (a,b))
    return boost*(np.sum((b-a)*(b+a),axis=1,dtype=np.longdouble)/(ga+gb)+
                  beta*(b[:,2]-a[:,2]))*audit.REST_EV


def distribution(values,charge):
    weights = np.abs(charge);weights /= np.sum(weights,dtype=np.longdouble)
    return {'count':len(values),'mean_eV':float(audit.weighted_mean(values,weights)),
            'rms_eV':float(np.sqrt(audit.weighted_mean(values*values,weights))),
            'sigma_eV':float(audit.weighted_sigma(values,weights)),
            'max_abs_eV':float(np.max(np.abs(values)))}


def comparisons(cg,blocked):
    inlet = [convergence.world_comparison(r,'full',r,'cut') for r in (cg,blocked)]
    replacement = [convergence.world_comparison(cg,role,blocked,role) for role in ('full','cut')]
    prepared = convergence.preparation_comparison(cg,blocked)
    energy_vectors = {}
    interaction = {'meaning':'Change of the matched cut-minus-full discrepancy after solver replacement.',
                   'particles':{},'native_field_plane':{}}
    for event in ('detector-crossings','final-particles'):
        a,b = cg['worlds'],blocked['worlds']
        vectors = {'cg_cut_minus_full_eV':stable_energy_difference(a['full'][event],a['cut'][event],cg['metrics']),
                   'blocked_cut_minus_full_eV':stable_energy_difference(b['full'][event],b['cut'][event],cg['metrics']),
                   'blocked_full_minus_cg_full_eV':stable_energy_difference(a['full'][event],b['full'][event],cg['metrics']),
                   'blocked_cut_minus_cg_cut_eV':stable_energy_difference(a['cut'][event],b['cut'][event],cg['metrics'])}
        vectors['change_of_inlet_discrepancy_eV'] = vectors['blocked_cut_minus_full_eV']-vectors['cg_cut_minus_full_eV']
        interaction['particles'][event] = distribution(vectors['change_of_inlet_discrepancy_eV'],a['full'][event]['charge_C'])
        energy_vectors[event] = {'id':a['full'][event]['id'],'charge_C':a['full'][event]['charge_C'],**vectors}
    for kind in ('E','B'):
        cfull,ccut = (audit.field_vector(cg['worlds'][r]['field-plane'],kind) for r in ('full','cut'))
        bfull,bcut = (audit.field_vector(blocked['worlds'][r]['field-plane'],kind) for r in ('full','cut'))
        difference = (bcut-bfull)-(ccut-cfull)
        # Do not add a tiny difference back onto a much larger reference merely
        # to subtract it again: that would obscure this interaction diagnostic.
        numerator=np.sum(difference*difference,dtype=np.longdouble)
        denominator=np.sum(cfull*cfull,dtype=np.longdouble)
        norm={'samples':len(cfull),
              'relative_l2':float(np.sqrt(numerator/denominator)) if denominator else None,
              'vector_rms_difference':float(np.sqrt(numerator/len(cfull))),
              'reference_vector_rms':float(np.sqrt(denominator/len(cfull))),
              'max_abs_component_difference':float(np.max(np.abs(difference))),
              'all_values_identical':bool(np.all(difference==0)),
              'unit':'V/m' if kind=='E' else 'T',
              'denominator':'CG full native field-plane norm over the same complete event/pixel series'}
        interaction['native_field_plane'][kind] = norm
    return inlet,replacement,prepared,interaction,energy_vectors


def ballistic_replacement(cg,blocked):
    result={'meaning':'Blocked minus CG in the response left after subtracting the identical analytic ballistic baseline.',
            'denominator':'Relative values use the CG self-field correction, not the much larger total beam width or focusing response.'}
    for role in ('full','cut'):
        a=cg[role]['self_field_correction_to_ballistic']
        b=blocked[role]['self_field_correction_to_ballistic']
        result[role]={key:{'difference':b[key]-a[key],
                          'relative_difference_to_cg_correction':(b[key]-a[key])/a[key] if a[key] else None}
                      for key in a}
    return result


def solver_cost(record):
    m = record['metrics']
    online = sum(m[k] for k in ('cut_preparation_wall_s','cut_advance_wall_s',
                                'cut_remap_and_feed_wall_s','live_reference_wall_s'))
    initialization_audits={k:m[k] for k in ('blocked_oracle_donor_fill_wall_s',
                                           'blocked_independent_gauss_wall_s') if k in m}
    return {'case':record['name'],'measurements':audit.cost(m),
            'candidate_online_s':online,
            'candidate_with_preparation_s':m['initialization_wall_s']+m['held_prewarm_wall_s']+online,
            'initialization_verification_subtimers_s':initialization_audits,
            'full_with_preparation_s':m['initialization_wall_s']+m['held_prewarm_wall_s']+
                m['full_preparation_wall_s']+m['full_advance_wall_s'],
            'reference_field_bytes':m['reference_field_bytes'],
            'initial_temporary_bytes':m['initial_temporary_bytes'],
            'reference_runtime_uses_full_donor':m.get('reference_runtime_uses_full_donor',True),
            'reference_block_layers':m.get('reference_block_layers'),
            'reference_solver_beta':m.get('reference_solver_beta'),
            'reference_translation_beta':m.get('reference_translation_beta'),
            'blocked_independent_gauss_bytes':m.get('blocked_independent_gauss_bytes'),
            'generator_preparation_s':record['blocked_report']['preparationSeconds']
                if record['blocked_report'] else None,
            'blocked_report':record['blocked_report'],
            'scope':'Named probe timers are preserved. Candidate total includes recorded initialization and feed; blocked initialization additionally includes the separately listed full-donor fill and independent Gauss verification. Nested generator timers are not added again. The comparison harness retains the complete donor/full-world oracle. Buffer capacity and I/O counters are not measured process peak RSS or production scalability.'}


def replay(current,old,files):
    result = {}
    for role in ('full','cut'):
        for event in EVENTS:
            a = Path(current['directory'])/role/(event+'.csv');b = old/role/(event+'.csv')
            audit.require(b.is_file() and audit.digest(a) == audit.digest(b),
                          'Default CG scientific replay differs: '+str(b))
            result[role+'/'+event+'.csv'] = audit.digest(a);files.append(b)
    return {'previous_case':str(old),'all_twelve_scientific_csv_byte_identical':True,'sha256':result}


def plot(root,inlet,replacement,ballistic,particles,steps):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    field_colors = ('#0072B2','#D55E00');energy_colors = ('#009E73','#CC79A7')
    fig,axes = plt.subplots(2,2,figsize=(12,8.5),layout='constrained')
    for column,pairs in enumerate((inlet,replacement)):
        labels = ['CG','Blocked'] if column == 0 else ['Full world','Cut world']
        x = np.arange(2);width=.34;row_values=[[],[]]
        for i,kind in enumerate(('E','B')):
            values = [p['field_plane']['all'][kind]['relative_l2']*100 for p in pairs]
            row_values[0].extend(values)
            bars = axes[0,column].bar(x+(i-.5)*width,values,width,label=kind,color=field_colors[i])
            axes[0,column].bar_label(bars,fmt='%.3g',padding=3,fontsize=9)
        for i,(key,label) in enumerate((('mean_energy_difference_eV','Mean'),('rms_energy_difference_eV','Per-particle RMS'))):
            values = [p['particles']['detector-crossings'][key] for p in pairs]
            row_values[1].extend(values)
            bars = axes[1,column].bar(x+(i-.5)*width,values,width,label=label,color=energy_colors[i])
            axes[1,column].bar_label(bars,fmt='%.3g',padding=3,fontsize=9)
        for row in range(2):
            ax=axes[row,column];ax.set_xticks(x,labels);ax.axhline(0,color='#666666',lw=.7)
            ax.grid(axis='y',alpha=.2);ax.set_axisbelow(True);ax.legend(fontsize=9)
            ax.margins(y=.25);ax.ticklabel_format(axis='y',style='sci',scilimits=(-3,4),useMathText=True)
            if all(value==0 for value in row_values[row]):
                ax.set_ylim(0,1);ax.set_yticks([0]);ax.ticklabel_format(axis='y',style='plain')
                ax.text(.5,.5,'All compared values are exactly zero',ha='center',
                        transform=ax.transAxes,color='#555555',fontsize=10)
            elif min(row_values[row])>=0:
                ax.set_ylim(bottom=0)
    axes[0,0].set_title('Existing inlet approximation: cut versus matched full')
    axes[0,1].set_title('Solver replacement: blocked versus CG')
    axes[1,0].set_title('Auxiliary-plane energy: cut minus matched full')
    axes[1,1].set_title('Auxiliary-plane energy: blocked minus CG')
    for ax in axes[0,:]:ax.set_ylabel('Native field relative L2 difference [%]')
    for ax in axes[1,:]:ax.set_ylabel('Lab energy difference [eV per electron]')
    fig.suptitle(f'Blocked reference comparison — {particles} particles, {steps} boost steps\n'
                 'Each field difference uses its own labeled reference; no total-energy normalization',fontsize=13)
    for suffix in ('png','pdf'):fig.savefig(root/('solver-and-inlet-comparison.'+suffix),dpi=170)
    plt.close(fig)
    fig,axes=plt.subplots(2,2,figsize=(12,8.5),layout='constrained')
    colors=('#0072B2','#56B4E9','#D55E00','#E69F00')
    names=('CG/full','CG/cut','Blocked/full','Blocked/cut')
    metrics=(('x_rms_m','Self-field correction to x RMS [pm]'),
             ('y_rms_m','Self-field correction to y RMS [pm]'),
             ('x_u_covariance_m',r'Self-field correction to cov(x,$u_x$) [$10^{-12}$ m]'),
             ('y_u_covariance_m',r'Self-field correction to cov(y,$u_y$) [$10^{-12}$ m]'))
    for ax,(key,title) in zip(axes.flat,metrics):
        values=[b[role]['self_field_correction_to_ballistic'][key]*1e12
                for b in ballistic for role in ('full','cut')]
        bars=ax.bar(np.arange(4),values,color=colors)
        ax.bar_label(bars,fmt='%.5g',padding=3,fontsize=9)
        ax.set_xticks(np.arange(4),names,rotation=15);ax.set_title(title,fontsize=11)
        ax.axhline(0,color='#666666',lw=.7);ax.grid(axis='y',alpha=.2);ax.set_axisbelow(True);ax.margins(y=.25)
    fig.suptitle('Matched auxiliary-plane response after subtracting analytic ballistic focusing\n'
                 'Same released particles and charge weights; original event interpolation remains',fontsize=13)
    for suffix in ('png','pdf'):fig.savefig(root/('ballistic-self-field-comparison.'+suffix),dpi=170)
    plt.close(fig)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cg-case',type=Path,required=True)
    parser.add_argument('--blocked-case',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--held-case',type=Path,help='Optional matched-geometry ordinary-held comparison, not a truth reference')
    parser.add_argument('--replay-case',type=Path,help='Previous CG case whose twelve scientific CSVs must replay byte-for-byte')
    parser.add_argument('--no-plot',action='store_true')
    args=parser.parse_args();root=args.output.resolve()
    audit.require(not root.exists(),'Use a fresh output directory: '+str(root))
    files=[Path(__file__).resolve(),Path(audit.__file__).resolve(),Path(convergence.__file__).resolve()]
    cg=load_case(args.cg_case,'CG','cg',files)
    blocked=load_case(args.blocked_case,'Blocked','blocked',files)
    matching_cases(cg,blocked)
    inlet,replacement,prepared,interaction,vectors=comparisons(cg,blocked)
    ballistic=[];baseline=None
    for record in (cg,blocked):
        result,current=convergence.ballistic_auxiliary_comparison(record)
        if baseline is not None:audit.same_tables(baseline,current,'Common analytic ballistic baseline')
        else:baseline=current
        ballistic.append(result)
    optional_held=None
    if args.held_case:
        held=load_case(args.held_case,'Ordinary held',None,files)
        audit.require(held['metrics']['initial_self_field_model']=='relativistic-poisson',
                      '--held-case requires the original held-prepared model')
        matching_cases(cg,held,allow_different_preparation=True)
        held_ballistic,held_baseline=convergence.ballistic_auxiliary_comparison(held)
        audit.same_tables(baseline,held_baseline,'Ordinary held common analytic ballistic baseline')
        optional_held={'held_case':str(args.held_case.resolve()),
            'prepared_against_blocked':convergence.preparation_comparison(held,blocked),
            'full_against_blocked':convergence.world_comparison(held,'full',blocked,'full'),
            'cut_against_blocked':convergence.world_comparison(held,'cut',blocked,'cut'),
            'analytic_ballistic_self_field_response':held_ballistic,
            'cost':solver_cost(held),
            'scope':'Preparation-history sensitivity; held is not declared an exact or converged physical reference.'}
    replay_result=replay(cg,args.replay_case.resolve(),files) if args.replay_case else None
    report={'format':'aprl-blocked-reference-review-v1','structural_status':'PASS',
            'producer_initial_residual_contract':'PASS','physical_acceptance':'NOT_ASSIGNED_METRICS_ONLY',
            'cg_case':str(args.cg_case.resolve()),'blocked_case':str(args.blocked_case.resolve()),
            'identical_initial_and_release_records':True,
            'production_parameter_checks':{r['name']:r['production_checks'] for r in (cg,blocked)},
            'initial_relative_residual':{r['name']:r['metrics']['initial_relative_residual'] for r in (cg,blocked)},
            'initial_residual_tolerance':cg['metrics']['initial_tolerance'],
            'prepared_solver_difference':prepared,
            'matched_inlet_approximation':inlet,'solver_replacement':replacement,
            'change_of_inlet_discrepancy':interaction,
            'analytic_ballistic_self_field_response':ballistic,
            'ballistic_solver_replacement':ballistic_replacement(*ballistic),
            'native_inlet_audit':{r['name']:r['inlet_audit'] for r in (cg,blocked)},
            'costs':[solver_cost(r) for r in (cg,blocked)],
            'optional_ordinary_held':optional_held,'default_cg_replay':replay_result,
            'precision':'Recover original binary64 momenta/boost parameters, then use long-double rationalized gamma differences; original extended-precision energy columns retain direct parsing.',
            'scope':['Different correct solvers need not be bitwise equal at finite CG tolerance. Field-difference/CG-tolerance ratios are not mathematical forward-error bounds.',
                     'The existing frozen-inlet approximation and the added solver replacement are reported separately; solver agreement does not remove the existing inlet discrepancy.',
                     'Native prepared and plane E/B are staggered in space/time and are not collocated laboratory fields or a Poynting diagnostic.',
                     'Energy differences are reported in eV per electron with same-ID weights; no weak net change is normalized to the total beam energy.',
                     'Ballistic subtraction verifies the original magnetic fringe is outside the whole-step interpolation interval; the remaining response includes existing mesh and event-interpolation errors.',
                     'All particles pass the local auxiliary plane. This is not the production4.5m plane, actual right-particle exit, MPI or long-time weak-radiation acceptance.',
                     'Generator buffer capacities and file counters are not actual peak RSS; the comparison probe retains complete donor/full-world audit arrays.'],
            'input_sha256':{str(p):audit.digest(p) for p in sorted(set(files))}}
    root.mkdir(parents=True)
    for event,columns in vectors.items():
        with (root/(event+'-energy-differences.csv')).open('x',newline='') as stream:
            writer=csv.writer(stream);names=list(columns);writer.writerow(names)
            for row in zip(*(columns[name] for name in names)):writer.writerow([str(v) for v in row])
    audit.save(root/'BLOCKED_REVIEW.json',report)
    if not args.no_plot:plot(root,inlet,replacement,ballistic,cg['metrics']['particles'],cg['metrics']['steps'])
    print(json.dumps({'structural_status':'PASS','review':str(root/'BLOCKED_REVIEW.json'),
                      'physical_acceptance':'NOT_ASSIGNED_METRICS_ONLY'},sort_keys=True))


if __name__ == '__main__':
    main()
