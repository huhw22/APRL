"""Run closed-file controls; no production result is inferred from this test."""
import json
import os
from pathlib import Path
import re
import sys

import h5py
import numpy as np

from prepare import require,save,sha
from launcher import launch
from binding import binding
from analyze_closed import verify,bitwise
import report_discrete


def compare_nodes(first,second):
    result={}
    for p in sorted((first/'detectors').glob('*.h5')):
        q=second/'detectors'/p.name
        with h5py.File(p,'r') as a,h5py.File(q,'r') as b:
            x=a[next(iter(a))];y=b[next(iter(b))]
            metrics={}
            for key in ['fields','time_s','records']:
                if key not in x:continue
                aa=x[key][:];bb=y[key][:]
                require(aa.shape==bb.shape and aa.dtype==bb.dtype,'Node output layout changed')
                if aa.dtype.names:
                    particle='source_id' in aa.dtype.names
                    if particle:
                        aa=np.sort(aa,order='source_id');bb=np.sort(bb,order='source_id')
                    for field in aa.dtype.names:
                        if field.endswith('_low') or (particle and field=='particle_id'):continue
                        if aa[field].dtype.kind not in 'fc':
                            require(np.array_equal(aa[field],bb[field]),'Node particle identity differs')
                        elif 'seconds' not in field:
                            av=aa[field].astype(np.longdouble);bv=bb[field].astype(np.longdouble)
                            if field+'_low' in aa.dtype.names:
                                av=av+aa[field+'_low'];bv=bv+bb[field+'_low']
                            difference=float(np.max(abs(av-bv),initial=0))
                            scale=float(np.max(abs(av),initial=0))
                            require(difference<=2e-11*max(scale,1e-30),'Node science differs: '+field)
                            metrics[key+'/'+field]=dict(max_abs=difference,reference_scale=scale)
                else:
                    difference=float(np.max(abs(aa-bb),initial=0));scale=float(np.max(abs(aa),initial=0))
                    require(difference<=2e-10*max(scale,1e-30),'Node field/time differs')
                    metrics[key]=dict(max_abs=difference,reference_scale=scale)
            result[p.name]=metrics
    return result


def main(root):
    root=Path(root).resolve();plan=json.loads((root/'plan.json').read_text())
    require(os.environ.get('SLURM_JOB_ID') and int(os.environ['SLURM_JOB_NUM_NODES'])==2,
        'Use the two-node bounded acceptance allocation')
    frozen=json.loads((root/'FROZEN_INPUTS.json').read_text())
    require(all(sha(p)==v for p,v in frozen.items()),'Frozen controls changed')
    identity=json.loads((Path(plan['build_record'])/'COMPLETE.json').read_text())
    require(identity.get('passed') is True and sha(plan['binary'])==identity['binary_sha256'],
        'Unaccepted/different binary')
    results={}
    for case in plan['cases']:
        d=Path(case['directory']);ranks=case['ranks'];nodes=case['nodes']
        command=['srun','--mpi=pmix_v3','--exclusive','--exact','--nodes='+str(nodes),
            '--ntasks='+str(ranks),'--ntasks-per-node=16','--cpus-per-task=8',
            '--cpu-bind=cores','--kill-on-bad-exit=1','--wait=30']
        env={k:v for k,v in os.environ.items() if not k.startswith('APRL_')};env['APRL_PROFILE_INTERVAL']='0'
        print('START',case['name'],flush=True)
        launch(command+[str(Path(plan['workspace'])/'build/cluster-tools/binding_probe'),str(ranks),'8'],
            d,env,120,'binding.log')
        save(d/'BINDING.json',binding(d/'binding.log',nodes,ranks,8))
        runtime=launch(command+[case['binary'],str(d/'run.yaml')],d,env,case['limit_seconds'])
        log=(d/'aprl.log').read_text()
        require(all(s in log for s in ['Configured stop reached:','phase=solve-complete','phase=mpi-finalize-complete']),
            'No complete normal stop')
        require(not re.search(r'Simulation failure|OpenFabrics device|Segmentation fault',log),'Simulation error')
        files=verify(d)
        rows,sums,actual_ranks,steps=report_discrete.load(d/'ledger',d/'ledger/energy-ledger.h5')
        require(steps==plan['expected_steps'] and actual_ranks==ranks,'Wrong geometry/layout')
        _,quantities,ledger=report_discrete.calculate(rows,sums,steps)
        if case['shape']=='linear-binomial3':
            denominator=float(np.max(abs(rows['delta_particle_kinetic_box_J'])))
            require(denominator>0,'No particle-change scale in drift control')
            ledger['max_total_over_max_particle_change']=float(np.max(abs(quantities['discrete_total_residual_J']))/denominator)
            ledger['max_coupling_over_max_particle_change']=float(np.max(abs(quantities['particle_grid_coupling_gap_J']))/denominator)
            require(ledger['max_total_over_max_particle_change']<1e-7 and
                ledger['max_coupling_over_max_particle_change']<1e-7,'Three-point particle/field work pairing failed')
        for p in (d/'detectors').glob('*-substeps.h5'):
            with h5py.File(p,'r') as f:
                g=f['particle_substep_plane'];a=g['records'][:]
                require(int(g['complete'][()])==1 and len(a)==int(g['committed_records'][()])==plan['expected_particles'],
                    'Missing closed particle crossing')
                require(len(np.unique(a['source_id']))==len(a),'Duplicate crossing')
        result=dict(program_integrity_pass=True,physical_acceptance=False,steps=steps,
            runtime_seconds=runtime['seconds'],closed_files=files,strict_ledger=ledger)
        save(d/'RESULT.json',result);results[case['name']]=result
        print('COMPLETE',case['name'],flush=True)
    old_linear=bitwise(root/'old-linear',root/'new-linear')
    old_notch=bitwise(root/'old-notch',root/'new-notch')
    dual=compare_nodes(root/'binomial3-single',root/'binomial3-dual')
    require(all(sha(p)==v for p,v in frozen.items()),'Frozen controls changed during execution')
    save(root/'COMPLETE.json',dict(program_integrity_pass=True,physical_acceptance=False,
        cases=results,linear_science_bitwise=old_linear,notch_science_bitwise=old_notch,
        three_point_node_consistency=dual,binary_sha256=identity['binary_sha256']))
    print('ALL COMPLETE',flush=True)


if __name__=='__main__':main(sys.argv[1])
