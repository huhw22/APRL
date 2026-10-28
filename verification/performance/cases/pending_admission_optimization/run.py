#!/usr/bin/env python3
"""Execute a frozen science or bounded-full stage; never submits a job."""
import os
for name in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[name]='1'
import argparse
import importlib.util
import json
import math
from pathlib import Path
import re
import signal
import subprocess
import time
import h5py
import numpy as np
import yaml
from prepare import require, save, sha
from ledger_review import compare as compare_ledger, thresholds as ledger_thresholds


def load_module(path, name):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module


def launch(command,directory,env,limit,log_name='aprl.log'):
    full=['timeout','--signal=TERM','--kill-after=10',str(limit)]+command
    record={'command':full,'cwd':str(directory),'limit_s':limit,'started_unix':time.time(),
            'environment':{k:env.get(k) for k in ('OMP_NUM_THREADS','OMP_DYNAMIC','OMP_PLACES','OMP_PROC_BIND',
                'OMPI_MCA_btl','OMPI_MCA_btl_openib_warn_no_device_params_found','APRL_PROFILE_INTERVAL','LD_LIBRARY_PATH')}}
    start=time.monotonic(); process=None
    try:
        with (directory/log_name).open('x') as f:
            process=subprocess.Popen(full,cwd=directory,env=env,stdin=subprocess.DEVNULL,
                stdout=f,stderr=subprocess.STDOUT,start_new_session=True)
            record['exit_code']=process.wait(timeout=limit+15)
    finally:
        if process and process.poll() is None:
            os.killpg(process.pid,signal.SIGTERM)
            try: process.wait(timeout=10)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid,signal.SIGKILL); process.wait()
        record['seconds']=time.monotonic()-start
        save(directory/(log_name+'.launch.json'),record)
    return record


def binding_result(path,ranks,threads):
    text=path.read_text()
    rows=re.findall(r'\[binding\] rank=(\d+) thread=(\d+) host=(\S+) cpu=(\d+) package=(\d+) core=(\d+) affinity=(\S+) status=(\S+)',text)
    require('[binding] result=PASS failures=0' in text,'Binding probe failed')
    require(len(rows)==ranks*threads and all(r[-1]=='PASS' for r in rows),'Incomplete binding')
    require({(int(r[0]),int(r[1])) for r in rows}=={(i,j) for i in range(ranks) for j in range(threads)},'Rank/thread binding coverage')
    require(len({r[2] for r in rows})==1 and len({(r[2],r[4],r[5]) for r in rows})==ranks*threads,'Physical core overlap or multiple hosts')
    return {'pass':True,'ranks':ranks,'threads':threads,'host':rows[0][2],'physical_cores':len(rows),'sha256':sha(path)}


def scalar_fields(group):
    return {name:group[name][:] for name in ('records','time_s') if name in group}


def full_review(directory,item,duration_helper):
    record=json.loads((directory/'aprl.log.launch.json').read_text())
    text=(directory/'aprl.log').read_text(errors='replace')
    n=item['steps']
    status=duration_helper.check_duration_exit(text,record['exit_code'],steps=n)
    # The established parser recognizes expected MPI_Abort caused by duration
    # exhaustion after every output closes. Never accept scheduler/timeout kills.
    require(not re.search(r'CANCELLED|DUE TO TIME LIMIT|Timed out|Segmentation fault|Killed by signal|Exited with exit code (?:124|137|139|143)',text,re.I),
            'Unexpected timeout/cancellation/runtime failure')
    outputs=duration_helper.inspect_outputs(directory,False)
    require(re.search(r'Moving window completed: steps='+str(n)+r'\b',text),'Wrong full-step moving count')
    require('particles=41732775' in text,'Wrong production source count')
    require(re.search(r'Moving window: reference_nz=53504 field_nz=4099 .*ranks=16\b',text),'Wrong production mesh/layout')
    require('Compute threads per MPI rank: 8;' in text,'Wrong production threads')
    with h5py.File(directory/'detectors/exit-field.h5') as f:
        g=f['field_plane']; require(g['fields'].shape==(n+1,160,160),'Wrong field/time sample count')
    with h5py.File(directory/'ledger/energy-ledger.h5') as f:
        records=f['energy_ledger/records'][:]
        require(np.array_equal(records['step'],np.arange(0,n+1,10)),'Wrong energy-ledger cadence')
        require(math.isclose(float(records['time_box_s'][-1]),n*yaml.safe_load((directory/'run.yaml').read_text())['mesh']['cell_size'][2]/299792458.,rel_tol=1e-12),'Ledger final time differs')
    points={int(step):float(seconds) for step,seconds in re.findall(r'\[resource\] phase=(?:progress|final) step=(\d+) elapsed_s=([\d.eE+-]+)',text)}
    require(n in points and all(math.isfinite(x) and x>0 for x in points.values()),'Missing final timing')
    profiles=[]
    for line in text.splitlines():
        if line.startswith('[performance]'):
            profiles.append(dict(re.findall(r'(\w+)=([^ ]+)',line)))
    require(bool(profiles)==bool(item['profile_interval']),'Profile mode differs')
    preparation=json.loads((directory/'reference/preparation.json').read_text())
    result={'pass':True,'status':'TIMING_COMPLETE_PHYSICS_STOP_NOT_REACHED','simulation_complete':False,
        'steps':n,'exit_code':record['exit_code'],'duration_exit':status,'outputs':outputs,
        'progress_elapsed_s':points,'stepping_average_s':points[n]/n,
        'launcher_seconds_including_reference_io_teardown':record['seconds'],
        'reference_preparation_seconds':preparation['preparationSeconds'],'profile_enabled':bool(item['profile_interval']),
        'profile_records':profiles,'scope':'Known bounded duration exhaustion; not physical completion or production radiation acceptance.'}
    if n==100:
        require(50 in points and points[100]>points[50],'Missing50–100 timing interval')
        result['steps_50_100_average_s']=(points[100]-points[50])/50
    return result


def compare_full(reference,candidate,helper):
    require(helper.physical_card(reference/'run.yaml')==helper.physical_card(candidate/'run.yaml'),'Full pair cards differ')
    result={'fields':{},'ledger':{},'pass':True}
    with h5py.File(reference/'detectors/exit-field.h5') as fa,h5py.File(candidate/'detectors/exit-field.h5') as fb:
        a,b=fa['field_plane'],fb['field_plane']
        require(np.array_equal(a['time_s'][:],b['time_s'][:]),'Full-pair clocks differ')
        require(a['fields'].shape==b['fields'].shape,'Full-pair shapes differ')
        for field in ('electric_V_per_m','magnetic_T'):
            norm,error,maxdelta=np.longdouble(0),np.longdouble(0),np.longdouble(0)
            equal=True
            for i in range(len(a['fields'])):
                aa,bb=a['fields'][i][field],b['fields'][i][field]
                equal=equal and aa.tobytes()==bb.tobytes()
                x=np.asarray(aa,dtype=np.longdouble); d=np.asarray(bb,dtype=np.longdouble)-x
                norm+=np.sum(x*x); error+=np.sum(d*d); maxdelta=max(maxdelta,np.max(abs(d)))
            rel=float(np.sqrt(error/norm)) if norm else (0. if error==0 else None)
            passed=rel is not None and rel<=2e-12
            result['fields'][field]={'pass':passed,'relative_l2':rel,'max_abs':float(maxdelta),'bitwise_equal':equal}
            result['pass'] &= passed
    result['ledger']=compare_ledger(reference,candidate,41732775)
    result['pass'] &= result['ledger']['pass']
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--stage',choices=('science','full'),required=True)
    args=p.parse_args();root=args.root.resolve();stage=args.stage;execution=root/stage/'execution'
    execution.mkdir(exist_ok=True)
    require(os.environ.get('SLURM_JOB_ID') and int(os.environ.get('SLURM_JOB_NUM_NODES','0'))==1,'Require real one-node Slurm allocation')
    require(int(os.environ.get('SLURM_NTASKS','0'))==16 and int(os.environ.get('SLURM_CPUS_PER_TASK','0'))==8,'Require16x8 allocation')
    require(not (execution/'STARTED.json').exists(),'Stage already started; use a fresh campaign')
    plan=json.loads((root/'plan.json').read_text())
    declared=json.loads((root/'THRESHOLDS.json').read_text())
    require(declared['science_ledger']==ledger_thresholds(31675) and declared['full_ledger']==ledger_thresholds(41732775),'Ledger contract/platform precision differs from frozen thresholds')
    frozen=json.loads((root/'frozen.json').read_text())
    require(all(sha(path)==value for path,value in frozen.items()),'Frozen campaign changed')
    if stage=='full':
        science=json.loads((root/'SCIENCE_PASS.json').read_text())
        require(science['pass'] and science['binary_sha256']==plan['binary_sha256'],'Same-binary science gate has not passed')
        require(science['frozen_manifest_sha256']==sha(root/'frozen.json'),'Science used different frozen cards')
    helper=load_module(root/'scripts/smoke_review.py','pending_smoke_review')
    duration=load_module(root/'scripts/duration_review.py','pending_duration_review')
    save(execution/'STARTED.json',{'slurm_job_id':os.environ['SLURM_JOB_ID'],'stage':stage,'started_unix':time.time(),
        'binary_sha256':plan['binary_sha256'],'budget_s':plan['stage_budget_seconds']})
    start=time.monotonic();deadline=start+plan['stage_budget_seconds']
    def expired(signum,frame): raise TimeoutError('Stage1800s active budget exhausted')
    signal.signal(signal.SIGALRM,expired)
    def terminated(signum,frame): raise SystemExit(128+signum)
    signal.signal(signal.SIGTERM,terminated);signal.signal(signal.SIGINT,terminated)
    signal.setitimer(signal.ITIMER_REAL,plan['stage_budget_seconds'])
    env=os.environ.copy();env.update(OMP_NUM_THREADS='8',OMP_DYNAMIC='FALSE',OMP_PLACES='cores',OMP_PROC_BIND='close',
        OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1',HDF5_USE_FILE_LOCKING='FALSE',**plan['mpi_environment'])
    bindings={};summaries={};comparisons={};small_planes={}
    def launcher(ranks,directory):
        return ['srun','--mpi='+env.get('SRUN_MPI','pmix_v3'),'--exact','--nodes=1','--ntasks='+str(ranks),
            '--ntasks-per-node='+str(ranks),'--cpus-per-task=8','--cpu-bind=cores','--kill-on-bad-exit=1','--wait=30','--chdir='+str(directory)]
    try:
        runtime={}
        for variant,binary_name in plan['binaries'].items():
            binary=Path(binary_name); check_env=env.copy()
            check_env['LD_LIBRARY_PATH']=str(binary.parent/'runtime-libs')+':'+env.get('LD_LIBRARY_PATH','')
            command=['ldd','-r',str(binary)]
            checked=subprocess.run(command,env=check_env,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=30)
            runtime[variant]={'command':command,'exit_code':checked.returncode,'output':checked.stdout,
                              'LD_LIBRARY_PATH':check_env['LD_LIBRARY_PATH']}
            save(execution/'RUNTIME_LINKS.json',runtime)
            require(checked.returncode==0 and 'not found' not in checked.stdout and 'undefined symbol' not in checked.stdout,'Runtime dependency check failed: '+variant)
        if stage=='full':
            inp=Path(plan['full_input_file']);st=inp.stat()
            require((st.st_size,st.st_mtime_ns)==(plan['full_input_size_bytes'],plan['full_input_mtime_ns']),'Full input stat changed')
            require(sha(inp)==plan['full_input_sha256'],'Full input content changed')
            save(execution/'INPUT_VERIFIED.json',{'pass':True,'sha256':plan['full_input_sha256'],'records':41732775})
        for item in plan['cases'][stage]:
            directory=Path(item['directory']);ranks=item['ranks'];variant=item['variant']
            remaining=int(deadline-time.monotonic()-25)
            # Reserve the entire declared per-case allowance before starting a
            # costly full initialization. No automatic timeout extension.
            needed=item['timeout_s'] if stage=='full' else 60
            require(remaining>=needed,'Insufficient active budget to start '+item['name'])
            if ranks not in bindings:
                bindingdir=execution/f'binding-{ranks}x8';bindingdir.mkdir()
                r=launch(launcher(ranks,bindingdir)+[plan['binding_probe'],str(ranks),'8'],bindingdir,env,45,'binding.log')
                require(r['exit_code']==0,'Binding launcher failure')
                bindings[ranks]=binding_result(bindingdir/'binding.log',ranks,8)
            binary=Path(plan['binaries'][variant]);case_env=env.copy()
            case_env['LD_LIBRARY_PATH']=str(binary.parent/'runtime-libs')+':'+env.get('LD_LIBRARY_PATH','')
            case_env['APRL_PROFILE_INTERVAL']=str(item['profile_interval'])
            r=launch(launcher(ranks,directory)+[str(binary),str(directory/'run.yaml')],directory,case_env,
                min(item['timeout_s'],int(deadline-time.monotonic()-20)))
            if stage=='science':
                require(r['exit_code']==0,'Small simulation failed or was interrupted')
                helper.RANKS=ranks
                summaries[item['name']],small_planes[item['name']]=helper.inspect_case(directory,free=item['magnet']=='free-drift')
                if variant=='candidate':
                    base_name=item['name'].replace('-candidate','-baseline')
                    report=helper.compare(root/'science'/base_name,directory,small_planes[item['name']])
                    report['ledger']=compare_ledger(root/'science'/base_name,directory,31675)
                    report['pass'] &= report['ledger']['pass']
                    comparisons[item['name']]=report;save(directory/'PAIR_REVIEW.json',report)
                    require(report['pass'],'Same-layout scientific replay failed: '+item['name'])
            else:
                summaries[item['name']]=full_review(directory,item,duration)
                if variant=='candidate':
                    base_name=item['name'].replace('-candidate','-baseline')
                    report=compare_full(root/'full'/base_name,directory,helper)
                    comparisons[item['name']]=report;save(directory/'PAIR_REVIEW.json',report)
                    require(report['pass'],'Bounded scientific diagnostic comparison failed: '+item['name'])
            save(directory/'REVIEW.json',summaries[item['name']])
            save(execution/'PROGRESS.json',{'completed':list(summaries),'elapsed_s':time.monotonic()-start})
            print('REVIEWED '+item['name'],flush=True)
        if stage=='science':
            for magnet in ('undulator','free-drift'):
                name=f'{magnet}-16x8-candidate'
                report=helper.compare(root/'science'/f'{magnet}-1x8-candidate',root/'science'/name,small_planes[name])
                report['ledger']=compare_ledger(root/'science'/f'{magnet}-1x8-candidate',root/'science'/name,31675)
                report['pass'] &= report['ledger']['pass']
                comparisons[f'{magnet}-candidate-cross-layout']=report
                save(execution/f'{magnet}-CROSS_LAYOUT.json',report)
                require(report['pass'],'Candidate cross-layout science mismatch')
        result={'pass':True,'stage':stage,'binary_sha256':plan['binary_sha256'],'bindings':bindings,'cases':summaries,
            'comparisons':comparisons,'elapsed_s':time.monotonic()-start,'frozen_manifest_sha256':sha(root/'frozen.json'),
            'limits':plan['limits']}
        if stage=='full':
            a=summaries['throughput-baseline'];b=summaries['throughput-candidate']
            result['status']='TIMING_COMPLETE_PHYSICS_STOP_NOT_REACHED'
            result['simulation_complete']=False
            result['throughput']={key:{'baseline_s_per_step':a[key],'candidate_s_per_step':b[key],
                'reduction_fraction':1-b[key]/a[key]} for key in ('stepping_average_s','steps_50_100_average_s')}
            result['throughput_scope']='One same16x8 paired early100-step sample, profiling disabled; no repeat statistics or full-run extrapolation.'
        require(all(sha(path)==value for path,value in frozen.items()),'Frozen campaign changed during stage')
        save(root/('SCIENCE_PASS.json' if stage=='science' else 'TIMING_REVIEW.json'),result)
        print('STAGE_REVIEWED '+stage,flush=True)
    except BaseException as error:
        save(execution/'FAILED_OR_INCOMPLETE.json',{'pass':False,'error':str(error),'completed_cases':list(summaries),
            'elapsed_s':time.monotonic()-start,'scope':'No pending/partial case promoted to success; no auto-resume or extended allocation.'})
        raise
    finally:
        signal.setitimer(signal.ITIMER_REAL,0)


if __name__=='__main__': main()
