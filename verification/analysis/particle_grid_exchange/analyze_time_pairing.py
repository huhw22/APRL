#!/usr/bin/env python3
"""Decompose saved particle/mesh work without changing or replaying a PIC step."""
import argparse,csv,json,os
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import h5py,numpy as np
LD=np.longdouble


def read_rank(path):
    with path.open() as f:names=f.readline().strip().split(',')
    wanted=['step','applied_current_old_J','applied_current_mid_J','current_time_pairing_J',
            'unapplied_current_mid_J','field_operator_residual_J']
    a=np.loadtxt(path,delimiter=',',skiprows=1,usecols=[names.index(n) for n in wanted],dtype=LD)
    assert np.array_equal(a[:,0],np.arange(1,len(a)+1))
    assert path.with_suffix('.complete').read_text().split()==['1',str(len(a))]
    assert np.isfinite(a).all()
    bound=LD('1e-17')*np.maximum(abs(a[:,1])+abs(a[:,2]),LD('1e-100'))
    assert np.all(abs(a[:,3]-a[:,2]+a[:,1])<bound)
    return a[:,1:]


def case(run,out,workers):
    with h5py.File(run/'ledger/energy-ledger.h5') as f:
        g=f['energy_ledger'];assert int(g['complete'][()])==1
        rows=g['records'][:];assert len(rows)==int(g['committed_records'][()]);ranks=int(g.attrs['mpi_size'])
    files=[run/f'ledger/discrete-field-rank-{i:05d}.csv' for i in range(ranks)]
    before={str(p):[p.stat().st_size,p.stat().st_mtime_ns] for p in [run/'ledger/energy-ledger.h5',*files]}
    total=None
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for a in pool.map(read_rank,files):total=a if total is None else total+a
    index=rows['step'].astype(np.int64);assert index[0]==0 and index[-1]==len(total)
    cumulative=np.vstack([np.zeros((1,total.shape[1]),dtype=LD),np.cumsum(total,axis=0,dtype=LD)])[index]
    old,mid,pair,unapplied,field=cumulative.T
    get=lambda n:np.asarray(rows[n],dtype=LD)
    sources=get('boundary_source_work_J')+get('frozen_constraint_work_J')+get('carrier_constraint_work_J')
    particle=get('particle_grid_work_J')-sources
    residual_old=particle-old;temporal=old-mid;residual=particle-mid
    assert np.max(abs(residual_old+temporal-residual))<LD('1e-22')
    assert np.max(abs(temporal+pair))<LD('1e-22')
    columns={'step':index,'boost_time_s':rows['time_box_s'],'particle_grid_work_corrected_J':particle,
      'mesh_old_work_J':old,'mesh_mid_work_J':mid,'particle_minus_old_mesh_J':residual_old,
      'old_minus_mid_mesh_J':temporal,'particle_minus_mid_mesh_J':residual,'unapplied_mid_work_J':unapplied,
      'field_operator_residual_J':field,'particle_kinetic_change_J':get('delta_particle_kinetic_box_J')}
    out.mkdir()
    with (out/'decomposition.csv').open('w') as f:
        w=csv.writer(f);w.writerow(columns);w.writerows(zip(*columns.values()))
    final={n:float(v[-1]) for n,v in columns.items()};den=abs(residual[-1])
    result={'ranks':ranks,'steps':len(total),'final':final,
      'temporal_term_over_abs_original_gap':float(temporal[-1]/den),
      'remaining_old_field_term_over_abs_original_gap':float(residual_old[-1]/den),
      'scope':'Algebraic separation of actual saved mesh work on the unchanged trajectory. E_old is the field ledger begin state. This is not a midpoint-pusher trial or a claim of physical energy correction.',
      'source_stats_unchanged':before,'new_pic_simulations':0}
    for p,v in before.items():assert [Path(p).stat().st_size,Path(p).stat().st_mtime_ns]==v
    (out/'RESULT.json').write_text(json.dumps(result,indent=2)+'\n')
    print(run.name,json.dumps({k:v for k,v in result.items() if k!='source_stats_unchanged'}),flush=True)
    return result


def main():
    p=argparse.ArgumentParser();p.add_argument('--workspace',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--workers',type=int,default=4);a=p.parse_args()
    assert os.environ.get('SLURM_JOB_ID'),'Read large saved rank traces in a compute allocation'
    a.out.mkdir(exist_ok=False)
    base=a.workspace/'results/corrected_full_comparison_20261030'
    all_results={}
    for job,rel in [(19099,'r3_2nodes/undulator_4p5m'),(19098,'r2/free_drift_0p5m')]:
        all_results[str(job)]=case(base/rel,a.out/str(job),a.workers)
    (a.out/'COMPLETE.json').write_text(json.dumps({'analysis_integrity_passed':True,'analysis_job':os.environ['SLURM_JOB_ID'],'new_pic_simulations':0,'cases':all_results},indent=2)+'\n')


if __name__=='__main__':main()
