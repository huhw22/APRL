#!/usr/bin/env python3
"""Separate held preparation sensitivity from subsequent frozen-inlet error.

The longest held preparation is a comparison reference, not a physical truth.
All quantities are offline diagnostics; no propagation or source is changed.
"""
import argparse
import json
import re
from pathlib import Path

import numpy as np

import review as audit


def precision_parameters(metrics):
    return (np.longdouble(np.float64(metrics['boost_gamma'])),
            np.longdouble(np.float64(metrics['window_beta'])))


def prepared_vectors(record, scope):
    t=record['prepared']; m=record['metrics']; initial=record['worlds']['full']['initial-particles']
    weights=np.abs(initial['charge_C']); weights/=np.sum(weights,dtype=np.longdouble)
    centers=np.asarray([np.sum(initial[name]*weights,dtype=np.longdouble) for name in
                        ('x_m','y_m','z_box_m')],dtype=np.longdouble)
    spacing=np.asarray([m['dx_m'],m['dy_m'],m['dz_m']],dtype=np.longdouble)
    origin=np.asarray(m['prepared_reference_roi']['origin_m'],dtype=np.longdouble)
    centers=(centers-origin)/spacing
    sigmas=np.asarray(m['sigma_position_cells'],dtype=np.longdouble)
    selection=np.ones(audit.size(t),dtype=bool)
    stagger=np.asarray([[.5,0,0],[0,.5,0],[0,0,.5],[0,.5,.5],[.5,0,.5],[.5,.5,0]],dtype=np.longdouble)
    components=t['component'].astype(int)
    if scope != 'whole_saved_roi':
        axes=(2,) if scope == 'source_pm3sigma' else (0,1,2)
        multiplier=1 if scope == 'source_pm1sigma_xyz' else 3
        for axis in axes:
            coordinate=t[('i','j','global_k')[axis]]+stagger[components,axis]
            selection &= np.abs(coordinate-centers[axis]) <= multiplier*sigmas[axis]
    result={}
    for kind,components in [('E',(0,1,2)),('B',(3,4,5))]:
        parts=[]; keys=[]
        for component in components:
            mask=selection&(t['component']==component)
            parts.append(t['value'][mask])
            keys.extend(zip(t['component'][mask],t['i'][mask],t['j'][mask],t['global_k'][mask],
                            t['native_time_box_s'][mask]))
        # Staggered components can have different counts near the ROI edge.
        # Aggregate scalar native samples without pretending they are collocated.
        result[kind]=(np.concatenate(parts).reshape(-1,1),keys)
    return result


def preparation_comparison(base,trial):
    result={'reference_case':base['name'],'trial_case':trial['name'],
            'meaning':'Difference of prepared states at identical source state and native lattice/time; no claim that the reference is the exact steady state.'}
    for scope in ('whole_saved_roi','source_pm3sigma','source_pm3sigma_xyz','source_pm1sigma_xyz'):
        a,b=prepared_vectors(base,scope),prepared_vectors(trial,scope)
        result[scope]={}
        for kind in ('E','B'):
            audit.require(a[kind][1] == b[kind][1], 'Prepared native ROI coordinates differ')
            metrics=audit.norm_metrics(a[kind][0],b[kind][0])
            metrics['native_scalar_rms_difference']=metrics.pop('vector_rms_difference')
            metrics['reference_native_scalar_rms']=metrics.pop('reference_vector_rms')
            result[scope][kind]=metrics
    return result


def world_comparison(reference, reference_role, trial, trial_role):
    a=reference['worlds'][reference_role]; b=trial['worlds'][trial_role]
    audit.same_tables(a['initial-particles'],b['initial-particles'],'Cross-case initial cohort')
    audit.same_tables(a['release-particles'],b['release-particles'],'Cross-case prescribed releases')
    boost,beta=precision_parameters(reference['metrics'])
    initial=a['initial-particles']
    result={'reference_case':reference['name']+'/'+reference_role,
            'trial_case':trial['name']+'/'+trial_role,'initial_and_release_records_identical':True,
            'particles':{event:audit.particle_comparison(a[event],b[event],initial,boost,beta)
                         for event in ('detector-crossings','final-particles')},
            'field_plane':audit.fields_comparison(a['field-plane'],b['field-plane'],reference['first_release'])}
    for event in ('detector-crossings','final-particles'):
        for axis in ('x','y'):
            moments=result['particles'][event]
            key=axis+'_u_covariance_change_m'
            moments[key+'_difference']=moments['trial'][key]-moments['reference'][key]
    a_history,b_history=a['history'],b['history']
    audit.require(np.array_equal(a_history['step'],b_history['step']),'Cross-case history times differ')
    history=[]
    for i,step in enumerate(a_history['step']):
        history.append({'step':int(step),
            'mean_delta_energy_difference_eV':float(b_history['mean_delta_energy_lab_eV'][i]-a_history['mean_delta_energy_lab_eV'][i]),
            'x_rms_difference_m':float(b_history['rms_x_m'][i]-a_history['rms_x_m'][i]),
            'y_rms_difference_m':float(b_history['rms_y_m'][i]-a_history['rms_y_m'][i]),
            'energy_sigma_difference_eV':float(b_history['sigma_energy_lab_eV'][i]-a_history['sigma_energy_lab_eV'][i])})
    result['history_differences']=history
    last_release=max(reference['last_release'],trial['last_release'])
    tail=[x for x in history if x['step'] >= last_release]
    result['all_released_tail']={'release_completed_step':last_release,'samples':len(tail),
        'first_sample_step':tail[0]['step'] if tail else None,
        'last_sample_step':tail[-1]['step'] if tail else None,
        'mean_energy_difference_first_eV':tail[0]['mean_delta_energy_difference_eV'] if tail else None,
        'mean_energy_difference_last_eV':tail[-1]['mean_delta_energy_difference_eV'] if tail else None,
        'scope':'Short global-time tail only; the same particle ID has the same release age in paired final records. No dense per-particle age-resolved force or relaxation series was recorded.'}
    # Matched final times and identical release times establish same age per ID.
    ages=[]
    for role,t in [('reference',a),('trial',b)]:
        released=t['release-particles']; detected=t['detector-crossings']; final=t['final-particles']
        audit.require(np.array_equal(released['id'],detected['id']) and
                      np.array_equal(released['id'],final['id']),'Event cohort mismatch')
        release_time=released['step']+released['fraction']
        detected_age=detected['step']+detected['fraction']-release_time
        final_age=final['step']+final['fraction']-release_time
        ages.append((detected_age,final_age))
        result[role+'_release_age_steps']={'detector_min':float(np.min(detected_age)),
            'detector_max':float(np.max(detected_age)),
            'final_min':float(np.min(final_age)),'final_max':float(np.max(final_age))}
    audit.require(np.array_equal(ages[0][1],ages[1][1]),'Final per-ID release ages differ')
    result['detector_max_abs_age_difference_steps']=float(np.max(np.abs(ages[0][0]-ages[1][0])))
    result['final_per_id_release_age_identical']=True
    return result


def ballistic_auxiliary_comparison(record):
    """Exact field-free continuation from the prescribed release event.

    This uses actual released u/gamma, never the held longitudinal velocity.
    The same fixed lab plane is intersected analytically; no simulation changes.
    """
    m=record['metrics'];full=record['worlds']['full'];cut=record['worlds']['cut']
    release=full['release-particles'];initial=full['initial-particles']
    audit.same_tables(release,cut['release-particles'],'Ballistic release source')
    audit.require(np.array_equal(release['id'],initial['id']),'Ballistic initial/release IDs differ')
    boost,beta=precision_parameters(m);u=audit.particle_momentum(release)
    audit.require(np.array_equal(u,audit.particle_momentum(initial)),'Release momentum differs from source')
    gamma=np.sqrt(1+np.sum(u*u,axis=1,dtype=np.longdouble));velocity=u/gamma[:,None]
    gap=np.longdouble(m['auxiliary_plane_lab_m']-m['activation_plane_lab_m'])/boost
    audit.require(gap>0 and np.all(velocity[:,2]+beta>0),'Ballistic plane cannot be reached')
    release0=np.longdouble((m['reference_offset0_cells']+m['lower_guard_cells'])*m['dz_m'])
    release_residual=np.max(np.abs(release['z_box_m']+beta*audit.C*release['time_box_s']-release0))
    coordinate_bound=64*np.finfo(float).eps*max(abs(release0+gap),abs(release0))
    audit.require(release_residual<coordinate_bound,'Recorded release is not on the prescribed plane')
    plane_residuals={}
    for role,world in [('full',full),('cut',cut)]:
        event=world['detector-crossings']
        plane_residuals[role]=float(np.max(np.abs(event['z_box_m']+beta*audit.C*event['time_box_s']-(release0+gap))))
        audit.require(plane_residuals[role]<coordinate_bound,'Recorded auxiliary event is not on the target plane')
        audit.require(np.array_equal(event['id'],release['id']) and np.array_equal(event['charge_C'],release['charge_C']),
                      'Ballistic event cohort/weights differ')
    flight=gap/(audit.C*(velocity[:,2]+beta))
    dt=np.longdouble(m['dt_s']);release_time=release['time_box_s']
    baseline={key:values.copy() for key,values in release.items()}
    for axis,key in enumerate(('x_m','y_m','z_box_m')):
        baseline[key]=release[key]+velocity[:,axis]*audit.C*flight
    event_time=release_time+flight;event_step=event_time/dt
    baseline['step']=np.floor(event_step);baseline['fraction']=event_step-baseline['step']
    baseline['time_box_s']=event_time
    baseline['time_lab_s']=release['time_lab_s']+boost*flight*(1+beta*velocity[:,2])
    baseline['frozen']=np.zeros_like(release['frozen'])
    for key in ('delta_energy_box_eV','delta_energy_lab_eV'):
        baseline[key]=np.zeros_like(release[key])
    fringe_extent=(np.sqrt(-2*np.log(m['undulator_fringe_cutoff']))*m['undulator_period_m']/(2*np.pi)
                   if m['undulator_gaussian_fringe'] else 0.)
    fringe_start=m['undulator_entrance_lab_m']-fringe_extent
    whole_step_advance=m['boost_gamma']*m['dz_m']*(1+m['window_beta'])
    margin=fringe_start-m['auxiliary_plane_lab_m']-whole_step_advance
    audit.require(margin>0,'Auxiliary event interpolation could sample nonzero prescribed magnet')
    result={'case':record['name'],'cohort_count':audit.size(release),
        'formula':'dt_box = ((z_detector_lab-z_release_lab)/gamma_boost) / [c*(u_z/gamma_particle+beta_boost)]; r = r_release + c*(u/gamma_particle)*dt_box.',
        'release_momentum_preserved':True,'charge_weights_preserved':True,
        'maximum_release_plane_residual_m':float(release_residual),
        'maximum_auxiliary_plane_residual_m':plane_residuals,
        'ballistic_flight_steps_min':float(np.min(flight/dt)),
        'ballistic_flight_steps_max':float(np.max(flight/dt)),
        'prescribed_magnet_zero_audit':{'fringe_start_lab_m':fringe_start,
            'auxiliary_plane_lab_m':m['auxiliary_plane_lab_m'],
            'maximum_possible_lab_advance_in_one_whole_boost_step_m':whole_step_advance,
            'remaining_margin_after_one_whole_step_m':margin,
            'source':'src/eb_sources.cpp: Gaussian compact fringe returns zero outside sqrt(-2 log cutoff)/k; the probe interpolates plane events within at most one complete boost step.'},
        'scope':'Difference from field-free ballistic focusing isolates the self-field plus its existing discrete propagation and detector interpolation response; it is not an exact-continuum space-charge validation.'}
    charge=np.abs(release['charge_C']);weights=charge/np.sum(charge,dtype=np.longdouble)
    def moments(t):
        p=audit.particle_momentum(t);out={}
        for axis,key in enumerate(('x_m','y_m')):
            x=t[key];mean_x=audit.weighted_mean(x,weights);mean_u=audit.weighted_mean(p[:,axis],weights)
            out[('x','y')[axis]+'_rms_m']=audit.weighted_sigma(x,weights)
            out[('x','y')[axis]+'_u_covariance_m']=audit.weighted_mean((x-mean_x)*(p[:,axis]-mean_u),weights)
        out['transverse_rms_radius_m']=np.sqrt(out['x_rms_m']**2+out['y_rms_m']**2)
        return out
    reference_moments=moments(baseline)
    result['ballistic_moments']={key:float(value) for key,value in reference_moments.items()}
    corrections={}
    for role,t in [('full',full['detector-crossings']),('cut',cut['detector-crossings'])]:
        comparison=audit.particle_comparison(baseline,t,initial,boost,beta)
        current=moments(t);corrections[role]={key:current[key]-reference_moments[key] for key in current}
        result[role]={'moments':{key:float(value) for key,value in current.items()},
            'self_field_correction_to_ballistic':{key:float(value) for key,value in corrections[role].items()},
            'energy_change_mean_eV':comparison['mean_energy_difference_eV'],
            'energy_change_rms_eV':comparison['rms_energy_difference_eV'],
            'max_abs_position_difference_from_ballistic_m':comparison['max_abs_position_difference_m']}
    result['cut_vs_full_self_correction']={key:{'difference':float(corrections['cut'][key]-corrections['full'][key]),
        'relative_difference_to_full_correction':float((corrections['cut'][key]-corrections['full'][key])/corrections['full'][key])
            if corrections['full'][key] != 0 else None} for key in reference_moments}
    return result,baseline


def geometry_sensitivity(reference_directory, trial):
    """Compare the unchanged full control after a documented box-origin shift."""
    reference_directory=reference_directory.resolve()
    old=json.loads((reference_directory/'metrics.json').read_text());new=trial['metrics']
    for key in ('steps','particles','held_prewarm_steps','initial_self_field_model',
                'dx_m','dy_m','dz_m','dt_s','nx','ny','boost_gamma','window_beta','source_moment_sample_sha256'):
        audit.require(old[key]==new[key], 'Geometry control changed another setting: '+key)
    a=audit.load_world(reference_directory/'full',old)['tables'];b=trial['worlds']['full']
    shift=np.longdouble((new['reference_offset0_cells']-old['reference_offset0_cells'])*new['dz_m'])
    initial_a,initial_b=a['initial-particles'],b['initial-particles']
    for key in ('id','x_m','y_m','ux','uy','uz','charge_C','frozen'):
        audit.require(np.array_equal(initial_a[key],initial_b[key]),'Geometry source changed: '+key)
    source_alignment=np.max(np.abs(initial_b['z_box_m']-initial_a['z_box_m']-shift))
    bound=16*np.finfo(float).eps*max(np.max(np.abs(initial_b['z_box_m'])),abs(shift))
    audit.require(source_alignment<=bound,'Geometry source z changed beyond coordinate rounding')
    plane_a,plane_b=a['field-plane'],b['field-plane']
    for key in ('step','pixel_i','pixel_j','x_m','y_m','time_E_box_s','time_B_box_s'):
        audit.require(np.array_equal(plane_a[key],plane_b[key]),'Geometry plane events differ: '+key)
    alignment=np.max(np.abs(plane_b['z_box_m']-plane_a['z_box_m']-shift))
    audit.require(alignment<=16*np.finfo(float).eps*max(np.max(np.abs(plane_b['z_box_m'])),abs(shift)),
                  'Geometry plane z changed beyond coordinate rounding')
    fields={}
    for label,mask in [('all',np.ones(audit.size(plane_a),dtype=bool)),('initial',plane_a['step']==0)]:
        fields[label]={kind:audit.norm_metrics(audit.field_vector(plane_a,kind)[mask],
                                              audit.field_vector(plane_b,kind)[mask]) for kind in ('E','B')}
    boost,beta=precision_parameters(new); particles={}
    for event in ('detector-crossings','final-particles'):
        item=audit.particle_comparison(a[event],b[event],initial_a,boost,beta)
        differences=audit.vector(b[event],('x_m','y_m','z_box_m'))-audit.vector(a[event],('x_m','y_m','z_box_m'))
        differences[:,2]-=shift
        item['max_abs_position_difference_m']=float(np.max(np.abs(differences)))
        item['position_comparison']='Box z aligned by the documented constant origin shift; event times can respond to different fields.'
        particles[event]=item
    return {'reference_case':str(reference_directory),'trial_case':trial['name'],
            'reference_full_nz':old['full_nz'],'trial_full_nz':new['full_nz'],
            'box_z_origin_shift_m':float(shift),'initial_source_z_alignment_residual_m':float(source_alignment),
            'plane_z_alignment_residual_m':float(alignment),'fields':fields,'particles':particles,
            'scope':'Measured finite-domain preparation sensitivity of the full control. Dynamic endpoint guard proofs do not guarantee finite-domain elliptic initialization independence.'}


def plot(root, records, preparation_pairs, paired):
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(2,2,figsize=(12,8),layout='constrained')
    colors=['#0072B2','#E69F00','#009E73','#CC79A7','#56B4E9']
    positions=np.arange(len(preparation_pairs)); width=.36
    for index,kind in enumerate(('E','B')):
        axes[0,0].bar(positions+(index-.5)*width,
            [p['source_pm3sigma'][kind]['relative_l2']*100 for p in preparation_pairs],
            width=width,color=colors[index],label=kind)
    axes[0,0].set_xticks(positions,[p['trial_case']+' vs\n'+p['reference_case'] for p in preparation_pairs],fontsize=7)
    axes[0,0].set_ylabel('Prepared native field difference [%]')
    axes[0,0].set_title('Same source state: ±3σz\nAll sampled transverse interior')
    positions=np.arange(len(records))
    for index,kind in enumerate(('E','B')):
        axes[0,1].bar(positions+(index-.5)*width,
            [p['field_plane']['all'][kind]['relative_l2']*100 for p in paired],
            width=width,color=colors[index],label=kind)
    axes[0,1].set_xticks(positions,[r['name'] for r in records],rotation=12,fontsize=8)
    axes[0,1].set_ylabel('Cut vs paired full: plane relative L2 [%]')
    axes[0,1].set_title('Later frozen-inlet prediction and propagation')
    for index,(event,label) in enumerate([('detector-crossings','Auxiliary plane'),('final-particles','Final state')]):
        axes[1,0].bar(positions+(index-.5)*width,
            [p['particles'][event]['rms_energy_difference_eV'] for p in paired],
            width=width,color=colors[index],label=label)
    axes[1,0].set_xticks(positions,[r['name'] for r in records],rotation=12,fontsize=8)
    axes[1,0].set_ylabel('Per-particle energy difference RMS [eV]')
    axes[1,0].set_title('Candidate error on the matched particle cohort')
    for color,record,p in zip(colors,records,paired):
        h=p['history_differences']
        axes[1,1].plot([v['step'] for v in h],[v['mean_delta_energy_difference_eV'] for v in h],
                       color=color,label=record['name'])
    axes[1,1].axvline(records[0]['last_release'],color='.4',linestyle=':',label='All particles released')
    axes[1,1].set_xlabel('Boost-frame step')
    axes[1,1].set_ylabel('Cut − full mean energy change [eV]')
    axes[1,1].set_title('Differences may cancel in the mean; see RMS')
    maximum=records[0]['metrics']['steps']
    axes[1,1].set_xticks(sorted(set(range(0,maximum+1,20))|{maximum}))
    for axis in axes.flat:
        axis.grid(axis='y',alpha=.25);axis.legend(fontsize=8)
    fig.suptitle('Held preparation and quasi-steady inlet: separate tests\n'
                 'Native staggered fields; auxiliary plane; longest preparation is not ground truth')
    fig.savefig(root/'convergence-comparison.png',dpi=180)
    fig.savefig(root/'convergence-comparison.pdf')
    plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(12,5.3))
    fig.subplots_adjust(top=.75,bottom=.15,wspace=.22,left=.07,right=.985)
    scopes=[('whole_saved_roi','Whole saved region'),('source_pm3sigma','Within ±3σz; all sampled x/y'),
            ('source_pm3sigma_xyz','Within ±3σ in x/y/z'),('source_pm1sigma_xyz','Within ±1σ in x/y/z')]
    x=np.arange(len(preparation_pairs));width=.19
    for axis,kind in zip(axes,('E','B')):
        for index,((scope,label),color) in enumerate(zip(scopes,colors)):
            axis.bar(x+(index-1.5)*width,[p[scope][kind]['relative_l2']*100 for p in preparation_pairs],
                     width=width,color=color,label=label)
        axis.set_xticks(x,[p['trial_case']+' vs\n'+p['reference_case'] for p in preparation_pairs],fontsize=8)
        axis.set_ylabel(kind+' prepared native relative L2 difference [%]')
        axis.grid(axis='y',alpha=.25)
    handles,labels=axes[0].get_legend_handles_labels()
    fig.legend(handles,labels,loc='upper center',bbox_to_anchor=(.5,.895),ncol=2,fontsize=9,frameon=False)
    fig.suptitle('Prepared-state differences by fixed source-centered region\n'
                 'Native staggered samples at t=0; these are not particle-force measurements')
    fig.savefig(root/'preparation-regions.png',dpi=180)
    fig.savefig(root/'preparation-regions.pdf')
    plt.close(fig)


def plot_ballistic(root, results):
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(2,2,figsize=(12,8))
    fig.subplots_adjust(top=.82,bottom=.08,hspace=.38,wspace=.25)
    positions=np.arange(len(results));width=.36;colors=['#0072B2','#E69F00']
    measures=[('x_rms_m',r'$\Delta\sigma_x$ [pm]'),('y_rms_m',r'$\Delta\sigma_y$ [pm]'),
        ('x_u_covariance_m',r'$\Delta\mathrm{cov}(x,u_x)$ [$10^{-12}$ m]'),
        ('y_u_covariance_m',r'$\Delta\mathrm{cov}(y,u_y)$ [$10^{-12}$ m]')]
    for axis,(key,label) in zip(axes.flat,measures):
        maximum=0.
        for index,(role,name) in enumerate([('full','Full domain'),('cut','Cut domain')]):
            values=[r[role]['self_field_correction_to_ballistic'][key]*1e12 for r in results]
            maximum=max(maximum,max(values))
            axis.bar(positions+(index-.5)*width,values,width=width,color=colors[index],label=name)
        for position,result in zip(positions,results):
            value=max(result[role]['self_field_correction_to_ballistic'][key]*1e12 for role in ('full','cut'))
            relative=result['cut_vs_full_self_correction'][key]['relative_difference_to_full_correction']
            axis.annotate(f'{100*relative:+.3f}%',(position,value),xytext=(0,6),textcoords='offset points',
                          ha='center',fontsize=8,color='#333333')
        axis.set_ylim(0,maximum*1.2)
        axis.set_xticks(positions,[r['case'] for r in results]);axis.set_ylabel(label)
        axis.grid(axis='y',alpha=.25)
    handles,labels=axes[0,0].get_legend_handles_labels()
    fig.legend(handles,labels,loc='upper center',bbox_to_anchor=(.5,.895),ncol=2,frameon=False)
    fig.suptitle('Response to the grid self-field before the undulator\n'
                 'Subtract analytic ballistic focusing; percentages compare cut with full correction',y=.98)
    fig.savefig(root/'ballistic-self-field-response.png',dpi=180)
    fig.savefig(root/'ballistic-self-field-response.pdf')
    plt.close(fig)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,required=True)
    parser.add_argument('--campaign',type=Path,action='append',required=True)
    parser.add_argument('--horizon',type=int,default=256)
    parser.add_argument('--geometry-control',type=Path)
    parser.add_argument('--no-plot',action='store_true')
    args=parser.parse_args(); root=args.root.resolve()
    root.mkdir(parents=True,exist_ok=True)
    audit.require(not (root/'CONVERGENCE_REVIEW.json').exists(),'Use a fresh analysis directory')
    records=[]; sampling_records=[]; files=[Path(__file__).resolve(),Path(audit.__file__).resolve()]
    for campaign in args.campaign:
        campaign=campaign.resolve()
        audit.require((campaign/'COMPLETED.json').is_file(),'Unfinished campaign')
        review=json.loads((campaign/'REVIEW.json').read_text())
        audit.require(review['structural_status']=='PASS','Campaign audit did not pass')
        # Review must describe the exact existing raw data, not stale outputs.
        for path,sha in review['input_sha256'].items():
            audit.require(audit.digest(path)==sha,'Stale input to campaign review: '+path)
        files += [campaign/'plan.json',campaign/'COMPLETED.json',campaign/'REVIEW.json']
        plan=json.loads((campaign/'plan.json').read_text())
        production_log=campaign/'production-aprl.log'
        normalization=re.search(r'Macro-particle normalization: records=(\d+), represented electrons=([^,]+), electrons per macro range=\[([^,]+), ([^\]]+)\]',production_log.read_text())
        audit.require(normalization is not None,'Missing production particle normalization provenance')
        production_records=int(normalization.group(1));production_electrons=np.longdouble(normalization.group(2))
        production_weight=production_electrons/production_records
        production_range=[float(normalization.group(3)),float(normalization.group(4))]
        files.append(production_log)
        for specification in plan['cases']:
            directory=campaign/specification['name'];m=json.loads((directory/'metrics.json').read_text())
            if m.get('warm_horizon_steps',64) != args.horizon:
                continue
            audit.require(m['mode']=='frozen','Convergence aggregate expects frozen candidates')
            worlds={}
            for role in ('full','cut'):
                worlds[role]=audit.load_world(directory/role,m)['tables']
                files += list((directory/role).glob('*.csv'))
            release=worlds['full']['release-particles']
            times=release['step']+release['fraction']
            prepared=audit.table(directory/'prepared-reference.csv')
            audit.prepared_review(directory/'prepared-reference.csv',m)
            initial_cohort=worlds['full']['initial-particles']
            cohort_weights=np.abs(initial_cohort['charge_C']);cohort_weights/=np.sum(cohort_weights,dtype=np.longdouble)
            realized_rms_cells=[float(audit.weighted_sigma(initial_cohort[column],cohort_weights)/m[spacing])
                                for column,spacing in [('x_m','dx_m'),('y_m','dy_m'),('z_box_m','dz_m')]]
            represented=abs(np.longdouble(m['charge_C']))/np.longdouble('1.602176634e-19')
            weight=represented/m['particles']
            sampling_records.append({'case':directory.name,'probe_macroparticles':m['particles'],
                'configured_rms_cells':m['sigma_position_cells'],'realized_cohort_rms_cells':realized_rms_cells,
                'probe_represented_electrons':float(represented),'probe_electrons_per_macro':float(weight),
                'production_records':production_records,'production_represented_electrons':float(production_electrons),
                'production_mean_electrons_per_macro':float(production_weight),'production_electrons_per_macro_range':production_range,
                'macro_weight_ratio_probe_to_production':float(weight/production_weight),
                'macroparticles_at_matching_production_weight_approx':float(represented/production_weight),
                'scope':'The shortened source uses selected production moments as generator parameters and matching nominal line density, but its realized finite-cohort moments differ and its macro-particles are much heavier and sparser than production. CIC granularity is a candidate cause, not established without a particle-count study; no inverse-square-root extrapolation.'})
            records.append({'name':directory.name,'directory':str(directory),'metrics':m,
                'worlds':worlds,'prepared':prepared,'first_release':float(np.min(times)),
                'last_release':float(np.max(times))})
            files += [directory/'metrics.json',directory/'prepared-reference.csv',directory/'inlet.csv']
    audit.require(len(records)>=2,'Need at least two matched-horizon cases')
    audit.require(len({r['name'] for r in records})==len(records),'Duplicate case name')
    records.sort(key=lambda r:(r['metrics']['initial_self_field_model']=='cowan-relativistic-poisson',r['metrics']['held_prewarm_steps']))
    m0=records[0]['metrics']
    for record in records:
        m=record['metrics']
        for key in ('steps','particles','full_nz','cut_nz','reference_offset0_cells','nx','ny',
                    'dx_m','dy_m','dz_m','dt_s','boost_gamma','window_beta','source_moment_sample_sha256'):
            audit.require(m[key]==m0[key],'Cross-case geometry/source changed: '+key)
        for event in ('initial-particles','release-particles'):
            audit.same_tables(records[0]['worlds']['full'][event],record['worlds']['full'][event],
                              'Across preparations '+event)
        audit.require(m['released']==m['particles']==m['auxiliary_detected'],'Partial particle cohort')
    legacy=[r for r in records if r['metrics']['initial_self_field_model']=='relativistic-poisson']
    audit.require(legacy,'Need held-prepared original-model reference')
    longest=max(legacy,key=lambda r:r['metrics']['held_prewarm_steps'])
    adjacent=[preparation_comparison(a,b) for a,b in zip(legacy,legacy[1:])]
    cowan=[r for r in records if r['metrics']['initial_self_field_model']=='cowan-relativistic-poisson']
    prepared_pairs=adjacent+[preparation_comparison(longest,r) for r in cowan]
    paired=[world_comparison(r,'full',r,'cut') for r in records]
    ballistic=[];ballistic_reference=None
    for record in records:
        result,baseline=ballistic_auxiliary_comparison(record);ballistic.append(result)
        if ballistic_reference is not None:
            audit.same_tables(ballistic_reference,baseline,'Common analytic ballistic baseline')
        else:
            ballistic_reference=baseline
    import csv
    baseline_path=root/'analytic-ballistic-auxiliary-plane.csv'
    with baseline_path.open('x',newline='') as stream:
        writer=csv.writer(stream);names=list(ballistic_reference);writer.writerow(names)
        for i in range(audit.size(ballistic_reference)):
            writer.writerow([str(ballistic_reference[name][i]) for name in names])
    sensitivity=[world_comparison(longest,'full',r,'full') for r in records if r is not longest]
    total=[world_comparison(longest,'full',r,'cut') for r in cowan]
    held_family=[world_comparison(a,'full',b,'full')
                 for i,a in enumerate(legacy) for b in legacy[i+1:]]
    sensitivity_envelope={}
    for event in ('detector-crossings','final-particles'):
        sensitivity_envelope[event]={
            'maximum_pairwise_rms_energy_difference_eV':max((p['particles'][event]['rms_energy_difference_eV'] for p in held_family),default=None),
            'maximum_pairwise_abs_mean_energy_difference_eV':max((abs(p['particles'][event]['mean_energy_difference_eV']) for p in held_family),default=None),
            'meaning':'Observed sensitivity across the tested held preparations only; not a statistical confidence interval, guaranteed error bound or automatic acceptance tolerance.'}
    amplitudes=[]
    for record in records:
        m=record['metrics']; component_norms=m['prepared_reference_native_component_L2']
        item={'case':record['name'],'initial_model':m['initial_self_field_model'],
              'prewarm_steps':m['held_prewarm_steps'],'full_grid':{},'saved_native_regions':{}}
        for kind,indices in [('E',range(3)),('B',range(3,6))]:
            squared=sum(np.longdouble(component_norms[i])**2 for i in indices)
            item['full_grid'][kind]={'raw_stored_value_L2':float(np.sqrt(squared)),
                                   'raw_stored_value_squared_L2':float(squared)}
        for scope in ('whole_saved_roi','source_pm3sigma','source_pm3sigma_xyz','source_pm1sigma_xyz'):
            vectors=prepared_vectors(record,scope);item['saved_native_regions'][scope]={}
            for kind in ('E','B'):
                values=vectors[kind][0];squared=np.sum(values*values,dtype=np.longdouble)
                item['saved_native_regions'][scope][kind]={'native_scalar_samples':len(values),
                    'raw_native_scalar_L2':float(np.sqrt(squared)),
                    'raw_native_scalar_rms':float(np.sqrt(squared/len(values))),
                    'maximum_absolute_native_value':float(np.max(np.abs(values)))}
        amplitudes.append(item)
    legacy_amplitudes=[a for a in amplitudes if a['initial_model']=='relativistic-poisson']
    amplitude_increments=[]
    for a,b in zip(legacy_amplitudes,legacy_amplitudes[1:]):
        steps=b['prewarm_steps']-a['prewarm_steps']
        amplitude_increments.append({'from_case':a['case'],'to_case':b['case'],'additional_prewarm_steps':steps,
            'B_raw_squared_L2_difference':b['full_grid']['B']['raw_stored_value_squared_L2']-a['full_grid']['B']['raw_stored_value_squared_L2'],
            'B_raw_squared_L2_difference_per_step':(b['full_grid']['B']['raw_stored_value_squared_L2']-a['full_grid']['B']['raw_stored_value_squared_L2'])/steps})
    geometry_comparisons=[]
    if args.geometry_control:
        path=args.geometry_control.resolve(); old=json.loads((path/'metrics.json').read_text())
        matching=[r for r in records if r['metrics']['held_prewarm_steps']==old['held_prewarm_steps']
                  and r['metrics']['initial_self_field_model']==old['initial_self_field_model']]
        audit.require(matching,'No matching held preparation for requested geometry control')
        geometry_comparisons=[geometry_sensitivity(path,r) for r in matching]
        files += [path/'metrics.json']+list((path/'full').glob('*.csv'))
    data={'format':'aprl-held-quasisteady-convergence-review-v1','structural_status':'PASS',
        'physical_acceptance':'METRICS_ONLY_NO_AUTOMATIC_TOLERANCE',
        'comparison_reference':{'case':longest['name'],'warm_steps':longest['metrics']['held_prewarm_steps'],
            'meaning':'Longest tested preparation is an empirical comparison reference, not a converged or physical ground truth.'},
        'geometry':{key:m0[key] for key in ('steps','particles','full_nz','cut_nz','reference_offset0_cells','nx','ny','dx_m','dy_m','dz_m','dt_s')},
        'cases':[{'case':r['name'],'initial_model':r['metrics']['initial_self_field_model'],
                  'prewarm_steps':r['metrics']['held_prewarm_steps'],'directory':r['directory']} for r in records],
        'source_sampling_comparison':sampling_records,
        'prepared_native_amplitudes':amplitudes,
        'held_full_grid_B_squared_norm_increments':amplitude_increments,
        'amplitude_scope':'Raw native stored-value norms on fixed geometry. E/B occupy staggered positions and time layers; these are not radiation energy, a discrete conserved energy, or physical power. Each point is a separate preparation run returned to the same final source state, not consecutive samples of one trajectory. Three preparation times do not identify a unique cause or establish exponential instability.',
        'prepared_reference_differences':prepared_pairs,
        'analytic_ballistic_self_field_response':ballistic,
        'paired_cut_vs_full':paired,'full_preparation_sensitivity':sensitivity,
        'cowan_cut_vs_longest_held_full':total,
        'held_full_family_pairwise':held_family,'observed_held_preparation_sensitivity':sensitivity_envelope,
        'full_geometry_sensitivity':geometry_comparisons,
        'scope':['Preparation comparison holds source positions, momenta, native phase and domain geometry fixed.',
                 'Finite-time held reference convergence and frozen future prediction are separate questions.',
                 'Source ROIs use fixed configured sigmas about the common initial cohort charge centroid: z-only ±3sigma and Cartesian ±3sigma/±1sigma boxes; each sample uses its actual native stagger. No ROI threshold is selected from observed field differences.',
                 'These sampled source ROIs do not equal particle-position Lorentz force or an isolated space-charge effect.',
                 'No Lorentz-force series was recorded; particle energy/covariance are observed response, not an invented force diagnostic.',
                 'Same-ID final particles have the same release age; only a short all-released tail and an auxiliary-plane event are available.',
                 'Net size agreement alone does not prove retained space-charge expansion. The separate analytic ballistic comparison subtracts the strong original focusing at the same auxiliary plane and measures the grid self-field correction, including existing discretization and event interpolation.',
                 'Field L2 difference is a diagnostic, not the sole quasi-steady acceptance criterion; weak energy/radiation effects remain relevant.',
                 'The 64-particle surrogate is a sparse-source interface probe and does not directly accept or reject production quasi-steady behavior at the original macro-particle weight.',
                 'Older different-card eV differences provide context but are not transferred as a numerical tolerance for this case.',
                 'No claim of full-plane radiative power, large-run speedup or asymptotic convergence.'],
        'input_sha256':{str(p):audit.digest(p) for p in sorted(set(files))}}
    outputs=[baseline_path]
    if not args.no_plot:
        plot(root,records,prepared_pairs,paired)
        plot_ballistic(root,ballistic)
        outputs += sorted(root.glob('*.png'))+sorted(root.glob('*.pdf'))
    data['output_sha256']={str(p):audit.digest(p) for p in outputs}
    audit.save(root/'CONVERGENCE_REVIEW.json',data)
    print(json.dumps({'status':'PASS','cases':len(records),'reference':longest['name'],'result':str(root/'CONVERGENCE_REVIEW.json')}))


if __name__=='__main__':
    main()
