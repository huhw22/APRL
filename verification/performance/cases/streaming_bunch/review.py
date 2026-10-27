#!/usr/bin/env python3
"""Independent fixed-cohort review of the charged-bunch window prototype."""
import argparse
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


def table(path):
    a = np.genfromtxt(path, names=True, delimiter=',')
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
    n=m['particles']; initial=sorted_particles(directory/'initial-particles.csv')
    released=sorted_particles(directory/'release-particles.csv')
    detected=sorted_particles(directory/'detector-crossings.csv')
    require(len(initial)==len(released)==len(detected)==n,'Incomplete cohort')
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
    require(m['admitted']==n,'Input queue lost particles')
    require(m['max_continuity_relative']<=1e-10,'Interior continuity gate failed')
    require(m['max_open_charge_balance_C']/abs(m['charge_C'])<=1e-10,'Open charge flux gate failed')
    require(m['max_remap_charge_balance_C']/abs(m['charge_C'])<=1e-10,'Grid-remap charge gate failed')
    require(m['initial_gauss_relative']<=4e-10,'Initial Cowan solve gate failed')
    return {'directory':directory,'metrics':m,'initial':initial,'released':released,
            'detected':detected,'history':history,'field':field,'moments':moments(detected,initial)}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--runs',type=Path,nargs='+',required=True)
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    plans=[json.loads((r/'plan.json').read_text()) for r in args.runs]
    require(len({p['binary_sha256'] for p in plans})==1,'Compared executables differ')
    cases=[]
    for root,plan in zip(args.runs,plans):
        require((root/'COMPLETED.json').is_file(),'Campaign did not complete: '+str(root))
        for item in plan['cases']:
            cases.append(inspect(root/item['name']))
    baselines={a['metrics']['beta_box']:a for a in cases if a['metrics']['mode']=='full'}
    output={'cohort_charge_and_input_consistency_status':'PASS','binary_sha256':plans[0]['binary_sha256'],
            'physical_status':'Measured guard/model differences; no production or pure-radiation acceptance',
            'cases':[]}
    for a in cases:
        m=a['metrics']; base=baselines[m['beta_box']]
        for key in a['initial'].dtype.names:
            require(np.array_equal(a['initial'][key],base['initial'][key]),'Input differs: '+key)
        for key in a['released'].dtype.names:
            require(np.array_equal(a['released'][key],base['released'][key]),'Prescribed release differs: '+key)
        for key in ('step','pixel_i','pixel_j','time_E_box_s','time_B_box_s','x_m','y_m','z_box_m'):
            require(np.array_equal(a['field'][key],base['field'][key]),'Field sample coordinates differ: '+key)
        entry={'path':str(a['directory'].resolve()),'metrics':m,'moments':a['moments'],
               'same_initial_and_release_cohort':True,
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
        entry['expansion_increment_relative_difference']=(
            entry['moment_difference']['radius_growth_percent']/base['moments']['radius_growth_percent'])
        entry['sigma_energy_relative_difference']=(
            entry['moment_difference']['sigma_energy_lab_eV']/base['moments']['sigma_energy_lab_eV'])
        w=weights(a['detected'])
        dt=a['detected']['time_lab_s']-base['detected']['time_lab_s']
        entry['arrival_time_difference']={'rms_s':float(np.sqrt(average(dt*dt,w))),
            'centered_rms_over_bunch_sigma':rms(dt,w)/rms(base['detected']['time_lab_s'],w)}
        entry['arrival_current_bunching_amplitude_difference']={key:
            a['moments']['bunching_centered_event_times'][key]['amplitude']-
            base['moments']['bunching_centered_event_times'][key]['amplitude']
            for key in a['moments']['bunching_centered_event_times']}
        output['cases'].append(entry)
    output['limitations']=['A full fixed reference cache and initial solve remain included.',
        'Particle detector moments use fixed IDs at their respective plane events, not a simultaneous bunch snapshot.',
        'Fields retain native E/B times; no mixed-time simultaneous lab E/B claim.',
        'Mur and missing outgoing-particle exterior evolution are approximations.',
        'One bunch length and one run per case cannot establish asymptotic or production speedup.']
    args.output.write_text(json.dumps(output,indent=2,allow_nan=False)+'\n')
    print(args.output)


if __name__=='__main__':
    main()
