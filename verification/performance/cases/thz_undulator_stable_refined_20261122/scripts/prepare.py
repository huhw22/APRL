"""Prepare fresh transverse-refined grids with a stable moving-inlet step ratio."""
from pathlib import Path
import hashlib
import json
import math
import shutil
import subprocess

import h5py
import numpy as np
import yaml
from scipy.constants import c, e, m_e, pi

ROOT = Path(__file__).resolve().parents[1]
WS = ROOT.parents[1]
SOURCE = WS/'THz/x7p8MV_phm180p0_a0p1050_BC.out'
BINARY = WS/'results/source_coupling_extended_20261122_r1/runtime/bin/aprl'
SHA = 'a04ebc4fb4f4ada60c1cf1fb667419a368bb232588c9cba98543fe7760d68c63'
NX, NY, H = 180, 120, 32e-6
CPML = 20
MAX_DZ_OVER_H = .6
GATE_STEPS = 256
F0, POINTS = 8.19e12, 40
OLD_DZ, OLD_LOWER, OLD_UPPER = 4e-6, 96, 192
NAMES = ['x', 'xp', 'y', 'yp', 't', 'p', 'particleID']


def save(path, data):
    path.write_text(json.dumps(data, ensure_ascii=False, indent=2)+'\n')


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(1<<20), b''):
            h.update(b)
    return h.hexdigest()


def main():
    assert sha(SOURCE) == SHA
    target = ROOT/SOURCE.name
    assert not target.exists()
    shutil.copy2(SOURCE, target)
    plain = ROOT/'beam.sdds-export.bin'
    dtype = np.dtype([(n,'<u8' if n=='particleID' else '<f8') for n in NAMES])
    assert int.from_bytes(plain.open('rb').read(4),'little') == 80000
    a = np.fromfile(plain,dtype=dtype,offset=4)
    assert len(a)==80000 and len(np.unique(a['particleID']))==80000
    assert all(np.isfinite(a[n]).all() for n in NAMES)
    records = np.zeros(len(a),dtype=[('plane_position_m','<f8',(2,)),('arrival_time_offset_s','<f8'),
                                    ('proper_velocity','<f8',(3,)),('source_id','<u8'),('macro_weight','<f8')])
    tref = np.mean(a['t'].astype(np.longdouble))
    records['plane_position_m'][:,0],records['plane_position_m'][:,1] = a['x'],a['y']
    records['arrival_time_offset_s'] = a['t'].astype(np.longdouble)-tref
    uz = a['p']/np.sqrt(1+a['xp']**2+a['yp']**2)
    records['proper_velocity'] = np.stack([uz*a['xp'],uz*a['yp'],uz],axis=1)
    records['source_id'],records['macro_weight'] = a['particleID'],1.
    with h5py.File(ROOT/'beam.h5','r') as f:
        saved=f['particles/records'][:]
        for name in records.dtype.names:
            assert np.array_equal(saved[name],records[name]),name
    gamma = np.sqrt(1+a['p']**2)
    gamma_mean = float(np.mean(gamma.astype(np.longdouble)))
    dt = records['arrival_time_offset_s']
    beta_z = uz/gamma
    centre = .001
    reference_time = float((np.longdouble(centre/c)+np.mean(beta_z.astype(np.longdouble)*dt))/np.mean(beta_z.astype(np.longdouble)))
    snapshot_z = beta_z*c*(reference_time-dt)
    input_times = -snapshot_z/(beta_z*c)
    rows=[]
    for divisor in (1,2,3):
        case=ROOT/f'f0-div-{divisor}';case.mkdir()
        f=8.19e12/divisor
        K=math.sqrt(2*(2*gamma_mean**2*c/(.05*f)-1))
        bg=gamma_mean/math.sqrt(1+K*K/2)
        bb=math.sqrt(1-1/bg**2)
        wavelength_box=(bg+math.sqrt(bg*bg-1))*c/F0
        requested_dz=wavelength_box/POINTS
        DZ=min(requested_dz,MAX_DZ_OVER_H*H)
        LOWER=max(4,math.ceil(OLD_LOWER*OLD_DZ/DZ))
        UPPER=max(8,math.ceil(OLD_UPPER*OLD_DZ/DZ))
        held_u=bg*(uz-bb*np.sqrt(1+uz**2))
        held_beta=held_u/np.sqrt(1+held_u**2)
        t_origin=float(np.min(input_times-2*DZ/(bg*c*(bb+held_beta))))
        raw_z=-bg*c*(bb+held_beta)*(input_times-t_origin)
        shift=-.5*(float(raw_z.min())+float(raw_z.max()))
        # Time guard uses the slowest input energy and measured initial timing,
        # the undulator's period-average longitudinal momentum, and 15% slack.
        gmin=float(gamma.min())
        perpendicular_bound=K*math.cosh(2*pi*float(np.max(abs(a['y'])))/.05)+float(np.max(abs(uz*a['xp'])))
        beta_plan=math.sqrt(1-(1+.5*perpendicular_bound**2)/gmin**2)
        entrance=.01;exit_z=entrance+1.
        duration=1.15*bg*(float(input_times.max())-t_origin+exit_z*(1/beta_plan-bb)/c)
        duration_steps=math.ceil(duration*c/DZ)
        duration=(duration_steps+.25)*DZ/c
        held_mean=float(np.mean(held_beta.astype(np.longdouble)))
        # Choose source/reference extent from the exact moving-reference query
        # bounds, while retaining the whole initial particle distribution.
        window=math.ceil(exit_z/(bg*DZ))+LOWER+UPPER+32
        window=math.ceil(window/32)*32
        nz=math.ceil(2*((1+held_mean)*(duration_steps+2)+LOWER+4-shift/DZ)/256)*256
        nz=max(nz,math.ceil(2*(window+shift/DZ+4)/256)*256,
               math.ceil((float(raw_z.max()-raw_z.min())/DZ+256)/256)*256)
        release_index=shift/DZ+nz/2
        initial_offset=math.floor(release_index)-LOWER
        assert initial_offset-(1+held_mean)*(duration_steps+2)>1
        assert initial_offset+window<nz-1
        # The entire finite-aperture arrival clock: bunch timing, maximum
        # sinusoidal magnetic slippage, and the transverse triangle bound.
        beta_min=math.sqrt(1-(1+perpendicular_bound**2)/gmin**2)
        transverse_bound=math.hypot(.002+.00224,.0015+.00128)
        start=exit_z/c+float(input_times.min())-4e-15
        end=exit_z/c+float(input_times.max())+exit_z*(1/beta_min-1)/c+transverse_bound/c+4e-15
        time_step=1/(POINTS*F0)
        bins=math.ceil((end-start)/time_step)
        bins=math.ceil(bins/256)*256
        config={
          'units':{'length':'m','time':'s'},
          'runtime':{'mode':'throughput','resource_monitor':{'enabled':True,'progress_interval_steps':100,
             'calibration_steps':1,'memory_safety_factor':1.25,'time_safety_factor':1.25},
             'shutdown_watchdog':{'enabled':True,'timeout_seconds':600,'finalize_timeout_seconds':60}},
          'output':{'overwrite':False,'manifest':'run-manifest.yaml'},
          'mesh':{'field_solver':'cowan-z','particle_coupling':'conservative-source',
             'longitudinal_shape':'linear-binomial3','transverse_x_shape':'linear-binomial3',
             'current_filter':'none','self_field_feedback':True,'one_way_source_diagnostic':False,
             'ballistic_outside_transverse':False,'cells':[NX,NY,nz],
             'cell_size':[H,H,DZ],'center':[0.,0.,0.],'duration':duration,'boost_gamma':bg,
             'particle_steps_per_undulator_period':512,'maximum_particle_substeps':4096,
             'moving_window':{'enabled':True,'cells':window,'lower_guard_cells':LOWER,'upper_guard_cells':UPPER,
                 'block_layers':32,'reference_directory':'reference','reference_surface_current':False,
                 'inlet_charge_match':True,'reference_through_cpml':True,'open_transverse_reference':True,
                 'cpml_startup':'osiris','cpml_startup_directory':'startup/arrays'}},
          'radiation_resolution':{'enabled':True,'maximum_photon_energy_eV':F0*4.135667696e-15,
             'warning_grid_points_per_wavelength':39.999,'warning_maxwell_samples_per_cycle':39.999},
          'boundary':{'type':'cpml','cells':[CPML,CPML,0],'polynomial_order':3.,'target_reflection':1e-8,
              'kappa_max':1.,'alpha_fraction':0.},
          'initial_self_field':{'enabled':True,'model':'cowan-relativistic-poisson',
              'relative_tolerance':1e-8,'maximum_iterations':10000,'maximum_mean_transverse_beta':.001},
          'beam':{'reference':{'input_plane_z':0.,'initial_center_z':centre,'frozen_lab_plane_injection':True,
                    'hold_transverse_until_activation':True,'activation_plane_z':0.},
              'input':{'type':'hdf5','file':str(ROOT/'beam.h5'),'electrons':.7e-9/e,'position_offset':[0.,0.,0.]}},
          'sources':{'incident_waves':[],'magnetic_elements':[{'type':'planar-undulator','characteristic':True,
              'entrance_z':entrance,'period':.05,'periods':20,'strength_parameter':K,
              'polarization_angle_rad':pi/2,'gaussian_fringe':False,'end_ramp_periods':.5,'phase_rad':pi/2}]},
          'detectors':{'enabled':True,'directory':'detectors','field_planes':[],
              'particle_planes':[{'name':name,'z':z,'buffer_records':8192,'compression':0}
                  for name,z in [('input',0.),('undulator-entry',entrance),('middle',.51),('exit',exit_z)]]},
          'trajectory':{'enabled':True,'directory':'events','basename':'particles','rhythm':1e-6,'buffer_records':16384},
          'energy_ledger':{'enabled':True,'directory':'ledger','filename':'energy-ledger.h5','sample_interval_steps':100,
              'buffer_records':128,'compression':0,'warning_relative_tolerance':.05,'discrete_operator_audit':True},
          'particle_retirement':{'enabled':False},'stop':{'mode':'after-last-element'},
          'radiation_screen':{'output_file':'radiation.h5','nx':25,'ny':19,'x0_m':-.002,'y0_m':-.0015,
              'dx_m':.004/24,'dy_m':.003/18,'z_m':exit_z,'time_origin_s':exit_z/c-reference_time,
              'time_start_s':start,'time_step_s':time_step,'time_bins':bins}}
        (case/'run.yaml').write_text(yaml.safe_dump(config,sort_keys=False))
        preflight=case/'preflight';preflight.mkdir()
        probe=json.loads(json.dumps(config))
        probe['mesh']['moving_window'].update(cpml_startup='standard')
        probe['mesh']['moving_window'].pop('cpml_startup_directory')
        probe['mesh']['duration']=8.25*DZ/c
        probe['detectors']={'enabled':False}
        probe['trajectory']={'enabled':False};probe['energy_ledger']={'enabled':False}
        probe.pop('radiation_screen')
        probe['stop']={'mode':'after-last-element'}
        (preflight/'run.yaml').write_text(yaml.safe_dump(probe,sort_keys=False))
        gate=case/'stability-gate';gate.mkdir()
        short=json.loads(json.dumps(config))
        short['mesh']['duration']=(GATE_STEPS-.75)*DZ/c
        short['mesh']['moving_window']['cpml_startup_directory']=str(case/'startup/arrays')
        short['detectors']={'enabled':False};short['trajectory']={'enabled':False}
        short.pop('radiation_screen')
        short['energy_ledger']['sample_interval_steps']=1
        (gate/'run.yaml').write_text(yaml.safe_dump(short,sort_keys=False))
        row=dict(divisor=divisor,target_frequency_Hz=f,K=K,B0_T=K*2*pi*m_e*c/(e*.05),boost_gamma=bg,
                 reference_cells=[NX,NY,nz],moving_window_cells=[NX,NY,window],initial_offset_estimate=initial_offset,
                 duration_s=duration,duration_steps=duration_steps,lab_reference_time_s=reference_time,
                 cell_size_m=[H,H,DZ],maxwell_dt_s=DZ/c,target_radiation_frequency_Hz=F0,
                 requested_radiation_points_per_cycle=POINTS,requested_dz_m=requested_dz,
                 actual_forward_points_per_cycle=wavelength_box/DZ,
                 dz_limited_by_transverse_spacing=bool(requested_dz>MAX_DZ_OVER_H*H),
                 dz_over_transverse_spacing=DZ/H,local_inlet_stability_ratio_limit=math.sqrt(math.sqrt(2)-1),
                 local_inlet_Q_max=2*(DZ/H)**2*(6-(DZ/H)**4)/((DZ/H)**2+2),
                 cpml_cells=[CPML,CPML,0],stability_gate_steps=GATE_STEPS,
                 lower_guard_cells=LOWER,upper_guard_cells=UPPER,radiation_dt_s=time_step,
                 lab_frame_time_origin_s=t_origin,held_beta=held_mean,
                 raw_initial_z_range_m=[float(raw_z.min()),float(raw_z.max())],
                 radiation_time_bins=bins,radiation_window_s=[start,start+bins*time_step],
                 virtual_screen_storage_GiB_per_rank=48*25*19*bins/2**30)
        save(case/'PARAMETERS.json',row);rows.append(row)
    save(ROOT/'PLAN.json',dict(source=str(SOURCE),source_sha256=SHA,particles=80000,charge_C=.7e-9,
         mean_gamma=gamma_mean,mean_total_energy_MeV=gamma_mean*.51099895,executable=str(BINARY),
         executable_sha256=sha(BINARY),source_revision='d279179a895a999ceb75412f98d6b86146e8c9b3',
         particle_input_sha256=sha(ROOT/'beam.h5'),rows=rows,
         virtual_screen='acceleration E/B only at exact undulator exit; 25x19 pixels, x +/-2 mm and y +/-1.5 mm',
         model='self-consistent APRL 2.1 main; F3 in x/z; OSIRIS four-side startup regenerated per case',
         end_fields='20 periods within 1 m; internal matched half-period ramps at both ends; no external fringe',
         target_radiation_frequency_Hz=F0,requested_radiation_points_per_cycle=POINTS,
         resolution_scope='Forward on-axis 8.19 THz for every resonance case, including harmonic cases; Cowan dz/dx=dz/dy=0.6 below charge-matched inlet local limit sqrt(sqrt(2)-1); time-screen samples 40 per 8.19 THz cycle',
         old_campaign=str(WS/'results/thz_undulator_aprl21_20261122_r1'),
         generation_note='Transverse spacing refined 40 to 32 um; dz=19.2 um; full extents and 640 um CPML thickness preserved. Original input particles, K and boost preserved. Recompute all native reference and OSIRIS startup, then require a closed 256-step self-consistent stability gate before the full radiation run. This local inlet condition is not a complete beam/CPML convergence proof.'))
    print(json.dumps(rows,indent=2))


if __name__=='__main__': main()
