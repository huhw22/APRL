#!/usr/bin/env python3
"""Prepare/run a distinct long-wait full-lifecycle queue comparison; never submit."""
import os
for name in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'): os.environ[name]='1'
import argparse
import copy
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import re
import shutil
import signal
import sys
import time
import h5py
import numpy as np
import yaml

HERE=Path(__file__).resolve().parent
C=299792458.;COUNT=31675;RANKS=16;THREADS=8


def require(test,message):
    if not test: raise RuntimeError(message)


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda:f.read(1024*1024),b''):h.update(chunk)
    return h.hexdigest()


def save(path,value):
    Path(path).write_text(json.dumps(value,indent=2,sort_keys=True,allow_nan=False)+'\n')


def module(path,name):
    spec=importlib.util.spec_from_file_location(name,path);result=importlib.util.module_from_spec(spec);spec.loader.exec_module(result);return result


def geometry(rows,card,reference_nz=4096):
    """Kinematic planning only. All-particle completion remains a runtime gate."""
    ld=np.longdouble;u=np.asarray(rows['proper_velocity'],dtype=ld)
    lab_gamma=np.sqrt(1+np.sum(u*u,axis=1));lab_beta=u[:,2]/lab_gamma
    gamma=float(card['mesh']['boost_gamma']);gb=math.sqrt(gamma-1)*math.sqrt(gamma+1);beta=gb/gamma
    dz=float(card['mesh']['cell_size'][2]);dt=dz/C
    release=float(card['beam']['reference']['activation_plane_z'])
    times=np.asarray(rows['arrival_time_offset_s'],dtype=ld);weights=np.asarray(rows['macro_weight'],dtype=ld)
    mean_beta=np.sum(weights*lab_beta)/np.sum(weights);mean_bt=np.sum(weights*lab_beta*times)/np.sum(weights)
    # v4 common snapshot must lie after every supplied plane arrival. This
    # coordinate-only change preserves all release events, momenta and weights.
    minimum_center=ld(release)+ld(C)*(mean_beta*np.max(times)-mean_bt)
    old_center=float(card['beam']['reference']['initial_center_z'])
    center=math.nextafter(max(old_center,float(minimum_center)),math.inf)
    common=((ld(center)-ld(release))/ld(C)+mean_bt)/mean_beta
    require(common >= np.max(times)-ld('1e-22'),'v4 would back-propagate an arrival')
    axial_plus=(np.sqrt(1+u[:,2]**2)+u[:,2])/(gamma+gb)
    held=(axial_plus*axial_plus-1)/(axial_plus*axial_plus+1)
    require(np.all((held>=0)&(held<1)),'Unsupported held velocity')
    input_times=times-common
    t0=np.min(input_times-2*dz/(gamma*C*(beta+held)))
    activation=gamma*(input_times-t0)
    initial=-C*(beta+held)*activation
    translation=-.5*(np.min(initial)+np.max(initial))
    field_z=float(card['detectors']['field_planes'][0]['z'])
    active=math.ceil((field_z-release)/(gamma*dz))+64+128+2
    require(active==289 and active>=8*RANKS,'Active field geometry changed')
    origin=-.5*reference_nz*dz
    offset=math.floor(float((translation-origin)/dz))-64
    active_origin=origin+offset*dz
    # This short horizon always has the same existing 4+ceil((1-beta)*N)=5 collar.
    safe_exit=release+gamma*(active_origin+(active-5)*dz-translation)
    rest2=1+u[:,0]**2+u[:,1]**2
    plus=(lab_gamma+u[:,2])/(gamma+gb);minus=rest2/plus
    beta_box=(plus-minus)/(plus+minus)
    require(np.all(beta+beta_box>0),'Ballistic release-to-exit event unreachable')
    ballistic_exit=activation+(safe_exit-release)/(gamma*C*(beta+beta_box))
    ballistic_steps=float(np.max(ballistic_exit)/dt)
    # Reuse the established1.1 ballistic scheduling allowance; this is not an
    # assertion about the nonlinear trajectory. Runtime requires every ID exit.
    planned_steps=math.ceil(1.1*ballistic_steps)+2
    require(planned_steps<=1000,'Long queue fixture exceeds agreed short horizon')
    duration=(planned_steps+1.5)*dt;query_steps=math.ceil(duration/dt)+1
    require(4+math.ceil((1-beta)*query_steps)==5,'Causal collar changed')
    bound_beta=float(np.max(held))+1e-12
    qmin=offset-(1+bound_beta)*query_steps-1;qmax=offset+active
    require(qmin>=0 and qmax<=reference_nz-1,'Reference scalar domain insufficient')
    centered=np.asarray(initial+translation,dtype=ld)
    require(np.min(centered)>origin+dz and np.max(centered)<origin+(reference_nz-1)*dz,'Source CIC touches longitudinal Dirichlet boundary')
    pending=(centered+np.asarray(held*dz,dtype=ld))<active_origin
    require(np.count_nonzero(pending)>=COUNT/2,'Long queue coverage below50%')
    reference_box=ld(translation)+ld(gamma)*(ld(center)-release)+ld(gb)*ld(C)*t0
    reference_lab0=ld(release)+ld(gamma)*(reference_box-translation)
    target=reference_lab0+ld(gb)*ld(C)*ld((planned_steps-.5)*dt)
    require(target>safe_exit and target>field_z,'Configured reference stop not downstream')
    return {'reference_cells_z':reference_nz,'active_cells_z':active,'initial_offset':offset,
        'lower_guard_cells':64,'upper_guard_cells':128,'source_records':COUNT,
        'old_initial_center_z_m':old_center,'minimum_v4_initial_center_z_m':float(minimum_center),
        'initial_center_z_m':center,'common_time_s':float(common),'frame_origin_lab_time_s':float(t0),
        'source_initial_untranslated_bounds_m':[float(np.min(initial)),float(np.max(initial))],
        'translation_z_m':float(translation),'reference_origin_z_m':origin,'active_origin_z_m':float(active_origin),
        'initial_pending_estimate':int(np.count_nonzero(pending)),'initial_pending_fraction_estimate':float(np.mean(pending)),
        'held_beta_range':[float(np.min(held)),float(np.max(held))],
        'last_release_step_estimate':float(np.max(activation)/dt),'ballistic_last_safe_exit_step_estimate':ballistic_steps,
        'safe_exit_lab_z_m':float(safe_exit),'initial_reference_lab_z_m':float(reference_lab0),
        'planning_factor':1.1,'complete_step_rounding_margin':2,'expected_stop_step':planned_steps,
        'stop_reference_lab_z_m':float(target),'duration_box_s':duration,'dt_box_s':dt,
        'reference_query_steps':query_steps,'reference_minimum_query_bound_cells':float(qmin),
        'reference_maximum_query_bound_cells':qmax,
        'scope':'Kinematic bounds and finite source-domain preparation; final all-particle release/three-plane/exit coverage must pass runtime validation, with no automatic horizon extension.'}


def prepare(args):
    parent=args.parent.resolve();root=args.root.resolve()
    require(not root.exists(),'Use a new independent long-queue directory')
    source_plan=json.loads((parent/'plan.json').read_text());science=json.loads((parent/'SCIENCE_PASS.json').read_text())
    require(science['pass'] and science['binary_sha256']==source_plan['binary_sha256'],'Same-binary short science not accepted')
    require(science['frozen_manifest_sha256']==sha(parent/'frozen.json'),'Parent science identity differs')
    for path,value in json.loads((parent/'frozen.json').read_text()).items():require(sha(path)==value,'Frozen parent changed: '+path)
    require(8<=args.stretch<=16,'Arrival stretch must be8–16')
    old_input=parent/'inputs/held.h5'
    with h5py.File(old_input) as f:
        rows=f['particles/records'][:];attrs=dict(f['particles'].attrs)
    require(len(rows)==COUNT and int(attrs['format_version'])==4,'Wrong source fixture')
    original=rows.copy();center=np.mean(rows['arrival_time_offset_s'].astype(np.longdouble))
    rows['arrival_time_offset_s']=np.asarray(center+np.longdouble(args.stretch)*(rows['arrival_time_offset_s'].astype(np.longdouble)-center),dtype=np.float64)
    for name in rows.dtype.names:
        if name!='arrival_time_offset_s':require(rows[name].tobytes()==original[name].tobytes(),'Changed non-arrival source field')
    old_card=parent/'science/undulator-16x8-baseline/run.yaml'
    base=yaml.safe_load(old_card.read_text());geo=geometry(rows,base)
    require(len(base['sources']['magnetic_elements'])==1 and base['sources']['magnetic_elements'][0]['periods']==1,'Expected original one-period undulator')
    root.mkdir(parents=True);scripts=root/'scripts';scripts.mkdir();inputs=root/'inputs';inputs.mkdir()
    # Use the already accepted and frozen comparison machinery unchanged.
    for name in ('prepare.py','run.py','ledger_review.py','smoke_review.py','duration_review.py'):
        shutil.copy2(parent/'scripts'/name,scripts/name)
    shutil.copy2(HERE/'long_queue.py',scripts/'long_queue.py');shutil.copy2(HERE/'long_queue.sbatch',scripts/'long_queue.sbatch')
    with h5py.File(inputs/'held-long.h5','x') as f:
        g=f.create_group('particles')
        for key,value in attrs.items():g.attrs[key]=value
        g.attrs['fixture_scope']='Controlled long-wait validation: arrival offsets stretched, same source IDs/xy/u/weights.'
        g.attrs['arrival_stretch_factor']=args.stretch;g.create_dataset('records',data=rows)
    cases={}
    for variant in ('baseline','candidate'):
        directory=root/variant;directory.mkdir();card=copy.deepcopy(base)
        card['beam']['input']['file']=str(inputs/'held-long.h5')
        card['beam']['reference']['initial_center_z']=geo['initial_center_z_m']
        card['mesh']['cells']=[160,160,geo['reference_cells_z']]
        card['mesh']['duration']=geo['duration_box_s'];card['mesh']['moving_window']['reference_directory']=str(directory/'reference')
        card['stop']={'mode':'reference-center-z','z':geo['stop_reference_lab_z_m']}
        card['output']['manifest']=str(directory/'run-manifest.yaml')
        for key,folder in (('detectors','detectors'),('trajectory','trajectory'),('energy_ledger','ledger')):card[key]['directory']=str(directory/folder)
        require(card['mesh']['moving_window']['lower_guard_cells']==64 and card['mesh']['moving_window']['upper_guard_cells']==128,'Original guards changed')
        (directory/'run.yaml').write_text(yaml.safe_dump(card,sort_keys=False));cases[variant]=str(directory)
    ledger=module(scripts/'ledger_review.py','long_queue_ledger_prepare')
    thresholds={'declared_before_execution':True,'field_relative_l2':2e-12,'particle_x_u_t_factor_eps64':128,
        'initial_pending_fraction_minimum':.5,'source_records':COUNT,'terminal_events':COUNT,
        'particle_planes':{'release':COUNT,'aux':COUNT,'exit':COUNT},'normal_physical_stop_required':True,
        'ledger':ledger.thresholds(COUNT),'case_timeout_s':300,'total_active_budget_s':840,'slurm_hard_limit_s':900}
    save(root/'THRESHOLDS.json',thresholds);save(root/'GEOMETRY.json',geo)
    plan={'format':'aprl-pending-long-queue-v1','parent':str(parent),'parent_science_sha256':sha(parent/'SCIENCE_PASS.json'),
        'parent_frozen_sha256':sha(parent/'frozen.json'),'binaries':source_plan['binaries'],'binary_sha256':source_plan['binary_sha256'],
        'binding_probe':source_plan['binding_probe'],'mpi_environment':source_plan['mpi_environment'],'cases':cases,
        'ranks':RANKS,'threads':THREADS,'arrival_stretch_factor':args.stretch,'arrival_center_s':float(center),
        'original_input_sha256':sha(old_input),'stretched_input_sha256':sha(inputs/'held-long.h5'),
        'source_card_sha256':sha(old_card),'geometry':geo,'thresholds':thresholds,
        'scope':'Distinct controlled long-wait fixture; never mix its runtime with original short fixture performance.',
        'source_changes':['Only arrival_time_offset_s stretched about its original mean; all source IDs, plane xy, momenta and weights identical.',
            'initial_center_z moved minimally to satisfy v4 forward common-time snapshot, with identical cards for both binaries.',
            'Explicit4096-cell disk scalar reference and complete-lifecycle horizon; active field289 cells and64/128 guards unchanged.']}
    save(root/'plan.json',plan)
    frozen=[root/'plan.json',root/'GEOMETRY.json',root/'THRESHOLDS.json',inputs/'held-long.h5',parent/'SCIENCE_PASS.json',
        parent/'frozen.json',Path(source_plan['binding_probe'])]+[p for p in scripts.iterdir() if p.is_file()]
    frozen += [Path(v)/'run.yaml' for v in cases.values()]
    for name,binary in source_plan['binaries'].items():
        binary=Path(binary);require(sha(binary)==source_plan['binary_sha256'][name],'Binary changed')
        frozen += [binary]+[p for p in (binary.parent/'runtime-libs').iterdir() if p.is_file()]
    save(root/'frozen.json',{str(p):sha(p) for p in frozen})
    print(json.dumps({'root':str(root),'geometry':geo,'submitted':False},indent=2))


def run(args):
    root=args.root.resolve();plan=json.loads((root/'plan.json').read_text());geo=plan['geometry']
    require(os.environ.get('SLURM_JOB_ID') and int(os.environ.get('SLURM_JOB_NUM_NODES','0'))==1,'Require a real one-node Slurm allocation')
    require(int(os.environ.get('SLURM_NTASKS','0'))==RANKS and int(os.environ.get('SLURM_CPUS_PER_TASK','0'))==THREADS,'Require16x8 allocation')
    execution=root/'execution';execution.mkdir(exist_ok=True);require(not (execution/'STARTED.json').exists(),'Already started')
    frozen=json.loads((root/'frozen.json').read_text());require(all(sha(p)==v for p,v in frozen.items()),'Frozen long-queue campaign changed')
    parent=Path(plan['parent']);proof=json.loads((parent/'SCIENCE_PASS.json').read_text())
    require(proof['pass'] and proof['binary_sha256']==plan['binary_sha256'],'Required same-binary science proof missing')
    sys.path.insert(0,str(root/'scripts'))
    driver=module(root/'scripts/run.py','long_queue_frozen_driver')
    helper=module(root/'scripts/smoke_review.py','long_queue_frozen_smoke');helper.RANKS=RANKS
    ledger=module(root/'scripts/ledger_review.py','long_queue_frozen_ledger')
    require(ledger.thresholds(COUNT)==plan['thresholds']['ledger'],'Ledger threshold/platform changed')
    start=time.monotonic();deadline=start+840
    save(execution/'STARTED.json',{'job_id':os.environ['SLURM_JOB_ID'],'time_unix':time.time(),'binary_sha256':plan['binary_sha256']})
    def timed_out(signum,frame):raise TimeoutError('Long-queue840s total budget exhausted')
    def cancelled(signum,frame):raise SystemExit(128+signum)
    signal.signal(signal.SIGALRM,timed_out);signal.signal(signal.SIGTERM,cancelled);signal.signal(signal.SIGINT,cancelled)
    signal.setitimer(signal.ITIMER_REAL,840)
    env=os.environ.copy();env.update(OMP_NUM_THREADS='8',OMP_DYNAMIC='FALSE',OMP_PLACES='cores',OMP_PROC_BIND='close',
        OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1',APRL_PROFILE_INTERVAL='0',HDF5_USE_FILE_LOCKING='FALSE',**plan['mpi_environment'])
    def command(directory):return ['srun','--mpi='+env.get('SRUN_MPI','pmix_v3'),'--exact','--nodes=1','--ntasks=16','--ntasks-per-node=16',
        '--cpus-per-task=8','--cpu-bind=cores','--kill-on-bad-exit=1','--wait=30','--chdir='+str(directory)]
    summaries={};planes={}
    try:
        record=driver.launch(command(execution)+[plan['binding_probe'],'16','8'],execution,env,45,'binding.log')
        require(record['exit_code']==0,'Binding launch failed');binding=driver.binding_result(execution/'binding.log',16,8)
        for variant in ('baseline','candidate'):
            require(deadline-time.monotonic()>=325,'Insufficient budget for next300s case')
            directory=Path(plan['cases'][variant]);binary=Path(plan['binaries'][variant]);case_env=env.copy()
            case_env['LD_LIBRARY_PATH']=str(binary.parent/'runtime-libs')+':'+env.get('LD_LIBRARY_PATH','')
            launched=driver.launch(command(directory)+[str(binary),str(directory/'run.yaml')],directory,case_env,300)
            require(launched['exit_code']==0,'Long-queue case interrupted or failed: '+variant)
            summaries[variant],planes[variant]=helper.inspect_case(directory)
            text=(directory/'aprl.log').read_text()
            m=re.search(r'Moving window: reference_nz=(\d+) field_nz=(\d+) initial_offset=(\d+)',text)
            require(m is not None and tuple(map(int,m.groups()))==(geo['reference_cells_z'],289,geo['initial_offset']),'Runtime geometry differs')
            require(summaries[variant]['steps']==geo['expected_stop_step'],'Configured stop full-step differs')
            if variant=='candidate':
                m=re.search(r'Pending admission queue: waiting=(\d+) resident=(\d+)',text)
                require(m is not None and int(m[1])>=COUNT/2,'Long pending coverage not exercised')
                summaries[variant]['actual_initial_pending']=int(m[1]);summaries[variant]['actual_initial_pending_fraction']=int(m[1])/COUNT
            save(directory/'LIFECYCLE_REVIEW.json',summaries[variant])
        comparison=helper.compare(Path(plan['cases']['baseline']),Path(plan['cases']['candidate']),planes['candidate'])
        comparison['ledger']=ledger.compare(Path(plan['cases']['baseline']),Path(plan['cases']['candidate']),COUNT)
        comparison['pass'] &= comparison['ledger']['pass'];save(root/'COMPARISON.json',comparison)
        require(comparison['pass'],'Long-queue strict field/particle/ledger comparison failed')
        require(all(sha(p)==v for p,v in frozen.items()),'Frozen files changed during execution')
        save(root/'LONG_QUEUE_PASS.json',{'pass':True,'binary_sha256':plan['binary_sha256'],'binding':binding,'cases':summaries,
            'comparison':comparison,'elapsed_s':time.monotonic()-start,'geometry':geo,'scope':plan['scope'],
            'limitations':['One31675-record controlled long-arrival fixture; not long-production radiation convergence or layout throughput comparison.']})
        print('LONG_QUEUE_PASS '+str(root/'LONG_QUEUE_PASS.json'),flush=True)
    except BaseException as error:
        save(execution/'FAILED_OR_INCOMPLETE.json',{'pass':False,'error':str(error),'completed_cases':list(summaries),'elapsed_s':time.monotonic()-start})
        raise
    finally:signal.setitimer(signal.ITIMER_REAL,0)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);sub=parser.add_subparsers(dest='action',required=True)
    p=sub.add_parser('prepare');p.add_argument('--parent',type=Path,required=True);p.add_argument('--root',type=Path,required=True);p.add_argument('--stretch',type=float,default=8.)
    p=sub.add_parser('run');p.add_argument('--root',type=Path,required=True)
    args=parser.parse_args();prepare(args) if args.action=='prepare' else run(args)
