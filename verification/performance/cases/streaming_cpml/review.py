#!/usr/bin/env python3
"""Independent fixed-cohort review of the CPML particle-inflow/reference prototype."""
import argparse
import csv
import hashlib
import math
import json
from pathlib import Path
import numpy as np

C = np.longdouble('299792458')
ME = np.longdouble('9.1093837139e-31')  # Same SI electron mass as eb_field.cpp.
QE = np.longdouble('1.602176634e-19')
BOOST = np.longdouble(5)/3
BETA = np.longdouble('.8')


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):h.update(block)
    return h.hexdigest()


def audit_campaign(root):
    require(not (root/'FAILED.json').exists(),'Failed campaign: '+str(root))
    plan=json.loads((root/'plan.json').read_text())
    completed=json.loads((root/'COMPLETED.json').read_text())
    require(plan['format']=='aprl-streaming-cpml-campaign-v1','Unexpected campaign type')
    require(completed['cohort_and_charge_gates']=='COMPLETED','Incomplete campaign')
    require(plan['ranks']==plan['threads']==1 and plan['case_timeout_s']==180 and
            plan['total_timeout_s']==7200,'Unexpected resource budget')
    require('[binding] result=PASS failures=0' in (root/'binding.log').read_text(),'Binding gate failed')
    for path,expected in json.loads((root/'frozen.json').read_text()).items():
        require(digest(Path(path))==expected,'Changed frozen input: '+path)
    require(digest(Path(plan['binary']))==plan['binary_sha256'],'Binary changed')
    names=[item['name'] for item in plan['cases']]
    require(len(names)==len(set(names)),'Duplicate planned case')
    require([item['case'] for item in completed['cases']]==names,'Completion differs from plan')
    for case in completed['cases']:
        directory=root/case['case'];launch=json.loads((directory/'launch.json').read_text())
        require(launch['exit_code']==0 and not launch.get('timed_out') and not launch.get('error'),
                'Failed case process: '+str(directory))
        require(digest(directory/'metrics.json')==case['metrics_sha256'],'Changed completed metrics')
    return plan


def csv_records_by_id(path):
    with path.open(newline='') as stream:
        rows=list(csv.reader(stream))
    return rows[0],sorted(rows[1:],key=lambda row:int(row[0]))


def ratio(value,scale):
    return float(value/scale) if scale else None


def table(path):
    a = np.atleast_1d(np.genfromtxt(path, names=True, delimiter=','))
    require(a.size > 0, 'Empty table: '+str(path))
    for key in a.dtype.names:
        require(np.isfinite(a[key]).all(), 'Nonfinite column '+str(path)+'/'+key)
    return a


def sorted_particles(path):
    a = table(path)
    require(len(np.unique(a['id'])) == len(a), 'Repeated particle event: '+str(path))
    return np.sort(a, order='id')


def weights(a):
    w = np.abs(a['charge_C']).astype(np.longdouble)
    return w/w.sum()


def average(a, w):
    return np.sum(np.asarray(a, dtype=np.longdouble)*w, axis=0)


def rms(a, w):
    return float(np.sqrt(average((a-average(a,w))**2,w)))


def moments(a, initial):
    w = weights(a)
    sx, sy = rms(a['x_m'],w), rms(a['y_m'],w)
    radius = np.hypot(sx,sy)
    r0 = np.hypot(rms(initial['x_m'],w),rms(initial['y_m'],w))
    u = np.column_stack([a[k] for k in ('ux','uy','uz')]).astype(np.longdouble)
    g = np.sqrt(1+(u*u).sum(axis=1))
    lab_energy = BOOST*(g+BETA*u[:,2])*ME*C*C/QE
    bunching = {}
    centered_time = a['time_lab_s'].astype(np.longdouble)-average(a['time_lab_s'],w)
    for wavelength in (1e-6,5e-6,10e-6):
        value = np.sum(w*np.exp(-2j*np.pi*np.asarray(C*centered_time/wavelength,dtype=float)))
        bunching[str(wavelength)] = {'real':float(value.real),'imag':float(value.imag),
                                    'amplitude':float(abs(value))}
    return {'rms_x_m':sx,'rms_y_m':sy,'rms_radius_m':float(radius),
            'radius_growth_percent':float((radius/r0-1)*100),
            'mean_delta_energy_box_eV':float(average(a['delta_energy_box_eV'],w)),
            'mean_delta_energy_lab_eV':float(average(a['delta_energy_lab_eV'],w)),
            'mean_energy_loss_lab_eV':-float(average(a['delta_energy_lab_eV'],w)),
            'sigma_energy_lab_eV':rms(lab_energy,w),
            'mean_time_lab_s':float(average(a['time_lab_s'],w)),
            'bunching_centered_event_times':bunching}


def verify_energy(a, initial):
    u=np.column_stack([a[k] for k in ('ux','uy','uz')]).astype(np.longdouble)
    u0=np.column_stack([initial[k] for k in ('ux','uy','uz')]).astype(np.longdouble)
    g=np.sqrt(1+(u*u).sum(axis=1)); g0=np.sqrt(1+(u0*u0).sum(axis=1))
    dg=((u-u0)*(u+u0)).sum(axis=1)/(g+g0)
    box=dg*ME*C*C/QE
    lab=BOOST*(dg+BETA*(u[:,2]-u0[:,2]))*ME*C*C/QE
    require(np.max(np.abs(box-a['delta_energy_box_eV']))<1e-9,'Box energy column inconsistent with momentum')
    require(np.max(np.abs(lab-a['delta_energy_lab_eV']))<1e-9,'Lab energy column inconsistent with momentum')


def difference(a,b):
    a,b=np.asarray(a,dtype=float),np.asarray(b,dtype=float)
    require(a.shape==b.shape,'Field/cohort shape mismatch')
    d=b-a; scale=np.linalg.norm(a.ravel()); error=np.linalg.norm(d.ravel())
    return {'relative_l2':float(error/scale) if scale else None,
            'rms_abs':float(np.sqrt(np.mean(d*d))), 'max_abs':float(np.max(np.abs(d))),
            'reference_l2':float(scale), 'error_l2':float(error)}


def inspect(directory):
    m=json.loads((directory/'metrics.json').read_text())
    require(m['completed'],'Incomplete probe')
    require(m['mode'] in ('full','cpml','zero') and m['beta_box'] in (0.,.23),'Unexpected physical case')
    require(m['history_policy']=='world-position' and m['cpml_cells']==[0,0,20],
            'Changed CPML history or transverse boundaries')
    require(m['nx']==m['ny']==24 and m['full_nz']==424 and m['physical_full_nz']==384 and
            m['steps']==176 and m['particles']==1024 and m['seed']==710271,'Changed controlled fixture')
    require(abs(m['full_origin_z_m']-(-84e-6))<1e-18,'Changed full-domain physical placement')
    require(m['inner_guard_cells']==m['guard_cells'] and m['cpml_z_cells']==20 and
            m['physical_window_nz']==m['window_nz']-40,'Guard includes the absorbing layer')
    require(m['guard_cells'] in (16,24,40) and
            m['window_nz']==(424 if m['mode']=='full' else 56+2*m['guard_cells']),
            'Unexpected candidate geometry')
    require(m['mode'] not in ('full','zero') or m['guard_cells']==24,'Changed full/negative-control label')
    require(m['cpml_history_peak_abs']>=0 and
            bool(m['cpml_history_nonzero_seen'])==(m['cpml_history_peak_abs']>0),
            'Inconsistent CPML history audit')
    require(m['mode']=='full' or m['cpml_history_nonzero_seen'],'Candidate did not exercise CPML history')
    require(m['pic_pushes_in_cpml']==0 and m['min_pic_inner_clearance_m']>=0,'PIC push entered CPML')
    n=m['particles']; initial=sorted_particles(directory/'initial-particles.csv')
    released=sorted_particles(directory/'release-particles.csv')
    detected=sorted_particles(directory/'detector-crossings.csv')
    require(len(initial)==len(released)==len(detected)==n,'Incomplete cohort')
    require(np.array_equal(initial['id'],np.arange(1,n+1)),'Missing fixed particle IDs')
    require(abs(float(np.sum(initial['charge_C'].astype(np.longdouble)))-m['charge_C'])
            <=abs(m['charge_C'])*1e-14,'Initial total charge differs')
    for a in (released,detected):
        require(np.array_equal(a['id'],initial['id']),'Cohort IDs changed')
        require(np.array_equal(a['charge_C'],initial['charge_C']),'Particle charge changed')
        require(np.all(a['frozen']==0),'Frozen event at release/detector')
    for key in ('x_m','y_m','ux','uy','uz'):
        require(np.array_equal(released[key],initial[key]),'Release changed prescribed '+key)
    require(np.all(detected['time_box_s']>released['time_box_s']),'Detector preceded release')
    verify_energy(released,initial)
    verify_energy(detected,initial)
    for a,plane in ((released,m['release_plane_lab_m']),(detected,m['detector_plane_lab_m'])):
        zlab=BOOST*(a['z_box_m'].astype(np.longdouble)+BETA*C*a['time_box_s'])
        tlab=BOOST*(a['time_box_s'].astype(np.longdouble)+BETA*a['z_box_m']/C)
        require(np.max(np.abs(zlab-plane))<1e-16,'Event is not on requested fixed lab plane')
        require(np.max(np.abs(tlab-a['time_lab_s']))<1e-25,'Incorrect Lorentz event time')
    history=table(directory/'history.csv'); field=table(directory/'field-plane.csv')
    require(len(history)==m['steps']+1,'Incomplete step history')
    require(np.array_equal(history['step'],np.arange(m['steps']+1)),'Missing or repeated history step')
    require(len(field)==m['field_plane_samples'],'Incomplete field plane output')
    sample_keys=np.column_stack([field[k] for k in ('step','pixel_i','pixel_j')])
    require(len(np.unique(sample_keys,axis=0))==len(field),'Repeated field sample')
    require(len(field)==45*20*20 and np.array_equal(np.unique(field['step']),np.arange(0,177,4)),
            'Changed field sampling cadence or plane size')
    for step in np.arange(0,177,4):
        frame=field[field['step']==step]
        require(len(frame)==400 and set(zip(frame['pixel_i'],frame['pixel_j']))==
                {(i,j) for i in range(2,22) for j in range(2,22)},'Incomplete plane pixel grid')
    require(np.all(history['active']==history['admitted']-history['evicted']),
            'Active queue accounting failed')
    require(np.all(np.diff(history['admitted'])>=0) and np.all(np.diff(history['evicted'])>=0) and
            np.all(history['held']<=history['active']) and history['admitted'][-1]==n,
            'Invalid admission/retirement history')
    expected_origin=(np.full(177,-84e-6) if m['mode']=='full' else
                     (np.floor(218.-.8*np.arange(177))-m['guard_cells']-20)*1e-6)
    require(np.max(np.abs(history['origin_z_m']-expected_origin))<1e-18,'Incorrect integer moving origin')
    inner_lo=history['origin_z_m']+20e-6
    inner_hi=history['origin_z_m']+(m['window_nz']-20)*1e-6
    release_plane=m['release_plane_lab_m']/float(BOOST)-.8*1e-6*history['step']
    detector_plane=m['detector_plane_lab_m']/float(BOOST)-.8*1e-6*history['step']
    require(np.min(release_plane-inner_lo)>=0 and np.min(inner_hi-detector_plane)>=0,
            'Fixed-lab plane entered the absorbing layer')
    require(abs(float(np.min(release_plane-inner_lo))-m['min_frame_release_inner_m'])<1e-18 and
            abs(float(np.min(inner_hi-detector_plane))-m['min_frame_detector_inner_m'])<1e-18,
            'Frame/CPML clearance audit differs from independent geometry')
    for events,key,side in ((released,'min_event_release_inner_m','lower'),
                            (detected,'min_event_detector_inner_m','upper')):
        event_origin=history['origin_z_m'][events['step'].astype(int)]
        clearance=(events['z_box_m']-event_origin-20e-6 if side=='lower' else
                   event_origin+(m['window_nz']-20)*1e-6-events['z_box_m'])
        require(np.min(clearance)>=0 and abs(float(np.min(clearance))-m[key])<1e-18,
                'Event/CPML clearance audit differs: '+key)
    require(np.max(np.abs(field['time_E_box_s']-field['step']*m['dt_s']))<1e-25 and
            np.max(np.abs(field['time_B_box_s']-(field['step']-.5)*m['dt_s']))<1e-25,
            'Changed native field time staggering')
    require(m['admitted']==n,'Input queue lost particles')
    require(m['max_continuity_relative']<=1e-10,'Interior continuity gate failed')
    require(m['max_open_charge_balance_C']/abs(m['charge_C'])<=1e-10,'Open charge flux gate failed')
    require(m['max_remap_charge_balance_C']/abs(m['charge_C'])<=1e-10,'Grid-remap charge gate failed')
    # Gauss residuals, including initialization, are recorded without a threshold.
    return {'directory':directory,'metrics':m,'initial':initial,'released':released,
            'detected':detected,'history':history,'field':field,'moments':moments(detected,initial)}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--runs',type=Path,nargs='+',required=True)
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    plans=[audit_campaign(r) for r in args.runs]
    require(len({p['binary_sha256'] for p in plans})==1,'Compared executables differ')
    cases=[]
    for root,plan in zip(args.runs,plans):
        for item in plan['cases']:
            case=inspect(root/item['name'])
            require(case['metrics']['beta_box']==float(item['beta']) and
                    case['metrics']['mode']==item['mode'] and
                    case['metrics']['guard_cells']==item['guard'],'Executed case differs from plan')
            cases.append(case)
    identities=[(a['metrics']['beta_box'],a['metrics']['mode'],a['metrics']['guard_cells']) for a in cases]
    require(len(identities)==len(set(identities)),'Repeated physical case in reviewed campaigns')
    baselines={a['metrics']['beta_box']:a for a in cases if a['metrics']['mode']=='full'}
    require(set(baselines)=={a['metrics']['beta_box'] for a in cases},'Missing same-beta full control')
    source_maps=[json.loads((r/'source-files.json').read_text()) for r in args.runs]
    core=[{path:sha for path,sha in m.items() if '/src/' in path} for m in source_maps]
    require(core[0] and all(m==core[0] for m in core),'Core sources differ across campaigns')
    output={'format':'aprl-streaming-cpml-review-v1','cohort_charge_and_input_consistency_status':'PASS',
            'binary_sha256':plans[0]['binary_sha256'],'campaigns':[str(r.resolve()) for r in args.runs],
            'review_script_sha256':digest(Path(__file__).resolve()),
            'physical_status':'Measured CPML-window and supply differences; no physical-error threshold or production acceptance',
            'gauss_policy':'Recorded only, including initial residual; no Gauss threshold gate',
            'coverage':{'case_count':len(cases),'complete_planned_ten_case_set':
                set(identities)=={(b,mode,g) for b in (0.,.23)
                                 for mode,g in [('full',24),('cpml',16),('cpml',24),('cpml',40),('zero',24)]}},
            'cases':[]}
    for a in cases:
        m=a['metrics']; base=baselines[m['beta_box']]
        for key in a['initial'].dtype.names:
            require(a['initial'][key].tobytes()==base['initial'][key].tobytes(),'Input bits differ: '+key)
        for key in a['released'].dtype.names:
            require(a['released'][key].tobytes()==base['released'][key].tobytes(),'Prescribed release bits differ: '+key)
        for name in ('initial-particles.csv','release-particles.csv'):
            require(csv_records_by_id(a['directory']/name)==csv_records_by_id(base['directory']/name),
                    'Canonical prescribed input/release CSV text differs: '+name)
        for key in ('step','pixel_i','pixel_j','time_E_box_s','time_B_box_s','time_E_lab_s','x_m','y_m','z_box_m'):
            require(np.array_equal(a['field'][key],base['field'][key]),'Field sample coordinates differ: '+key)
        entry={'path':str(a['directory'].resolve()),'baseline_path':str(base['directory'].resolve()),
               'metrics':m,'moments':a['moments'],'same_initial_and_release_cohort':True,
               'initial_and_release_numeric_bits_and_canonical_csv_identical':True,
               'gauss_observations':{key:value for key,value in m.items() if 'gauss' in key},
               'cpml_and_clearance_observations':{key:value for key,value in m.items() if
                   'cpml' in key or 'inner' in key or key=='history_policy'},
               'initial_csv_raw_byte_identical':(a['directory']/'initial-particles.csv').read_bytes()==
                                                (base['directory']/'initial-particles.csv').read_bytes(),
               'release_csv_raw_byte_identical':(a['directory']/'release-particles.csv').read_bytes()==
                                                (base['directory']/'release-particles.csv').read_bytes(),
               'field_difference':{key:difference(base['field'][key],a['field'][key]) for key in
                                   ('Ex_V_m','Ey_V_m','Ez_V_m','Bx_T','By_T','Bz_T')},
               'particle_difference':{key:difference(base['detected'][key],a['detected'][key]) for key in
                                      ('x_m','y_m','ux','uy','uz','time_lab_s','delta_energy_box_eV','delta_energy_lab_eV')},
               'moment_difference':{key:a['moments'][key]-base['moments'][key] for key in a['moments']
                                    if key!='bunching_centered_event_times'}}
        for vector,keys in [('E',('Ex_V_m','Ey_V_m','Ez_V_m')),('B',('Bx_T','By_T','Bz_T'))]:
            entry['field_difference'][vector]=difference(
                np.column_stack([base['field'][k] for k in keys]),
                np.column_stack([a['field'][k] for k in keys]))
        entry['expansion_increment_relative_difference']=ratio(
            entry['moment_difference']['radius_growth_percent'],base['moments']['radius_growth_percent'])
        entry['sigma_energy_relative_difference']=ratio(
            entry['moment_difference']['sigma_energy_lab_eV'],base['moments']['sigma_energy_lab_eV'])
        w=weights(a['detected'])
        dt=a['detected']['time_lab_s']-base['detected']['time_lab_s']
        entry['arrival_time_difference']={'rms_s':float(np.sqrt(average(dt*dt,w))),
            'centered_rms_over_bunch_sigma':ratio(rms(dt,w),rms(base['detected']['time_lab_s'],w))}
        entry['arrival_current_bunching_amplitude_difference']={key:
            a['moments']['bunching_centered_event_times'][key]['amplitude']-
            base['moments']['bunching_centered_event_times'][key]['amplitude']
            for key in a['moments']['bunching_centered_event_times']}
        entry['field_by_time']=[]
        for step in np.unique(a['field']['step']):
            mask=a['field']['step']==step;item={'step':int(step)}
            for vector,keys,scale in [('E',('Ex_V_m','Ey_V_m','Ez_V_m'),1.),
                                      ('cB',('Bx_T','By_T','Bz_T'),float(C))]:
                ref=np.column_stack([base['field'][key][mask] for key in keys])*scale
                value=np.column_stack([a['field'][key][mask] for key in keys])*scale
                item[vector]={'reference_rms_V_m':float(np.sqrt(np.mean(np.sum(ref*ref,axis=1)))),
                              'candidate_rms_V_m':float(np.sqrt(np.mean(np.sum(value*value,axis=1)))),
                              'difference_rms_V_m':float(np.sqrt(np.mean(np.sum((value-ref)**2,axis=1))))}
            entry['field_by_time'].append(item)
        entry['workload_observation']={
            'field_cell_updates_over_full':ratio(m['field_cell_updates'],base['metrics']['field_cell_updates']),
            'overall_wall_s_over_full':ratio(m['overall_wall_s'],base['metrics']['overall_wall_s']),
            'evolution_wall_s_over_full':ratio(m['evolution_wall_s'],base['metrics']['evolution_wall_s']),
            'interpretation':'One run per different domain with research diagnostics; not same-physics production speedup'}
        output['cases'].append(entry)
    output['limitations']=['A full fixed reference cache and initial solve remain included.',
        'Particle detector moments use fixed IDs at their respective plane events, not a simultaneous bunch snapshot.',
        'Fields retain native E/B times; no mixed-time simultaneous lab E/B claim.',
        'The original CPML acts on total fields, including the reference-supplied near field; it is not residual-field CPML.',
        'Only longitudinal CPML is added; original transverse PEC and original numerical kernels remain.',
        'Outgoing particles are ballistic only after detection; later omitted exterior return fields remain an approximation.',
        'Gauss residuals are recorded, not used to change the original discretization or reject its established compromises.',
        'Positive particle delta energy means gain; reported mean energy loss has the opposite sign.',
        'One bunch length and one run per case cannot establish asymptotic or production speedup.']
    args.output.write_text(json.dumps(output,indent=2,allow_nan=False)+'\n')
    print(args.output)


if __name__=='__main__':
    main()
