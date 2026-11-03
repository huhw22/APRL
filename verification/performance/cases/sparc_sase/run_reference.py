import os
from pathlib import Path
import h5py,numpy as np
from common import *
from launch import run
def main():
    root=Path(__file__).resolve().parents[1];plan=load(root/'plan.json')
    need(all(sha(p)==h for p,h in load(root/'FROZEN.json').items()),'Frozen reference changed')
    cmd=['srun','--mpi='+os.environ.get('SRUN_MPI','pmix_v3'),'--exact','--nodes=1',
         '--ntasks=64','--ntasks-per-node=64','--cpus-per-task=1','--cpu-bind=cores',
         '--kill-on-bad-exit=1','--wait=30',plan['genesis'],'reference.in']
    gen=root/'genesis';rec=run(cmd,gen,4*3600,'genesis.log')
    with h5py.File(gen/'reference.out.h5','r') as f:
        endpoint=float(f['Lattice/z'][-1]+f['Lattice/dz'][-1]);power=f['Field/power'][:]
        z=f['Lattice/zplot'][:]
        need(abs(endpoint-plan['geometry']['field_planes'][-1][1])<1e-8,'Reference lattice incomplete')
        need(np.isfinite(power).all(),'Nonfinite reference power')
        phases=f['Lattice/phaseshift'][:].ravel()
        observed=phases[np.abs(phases)>1e-12]
        expected=-np.array(load(gen/'lattice.json')['genesis_gap_phase_kicks_rad'])
        need(observed.shape==(5,) and np.allclose(observed,expected,atol=1e-12,rtol=0),
             'Reference segment phases were not applied once per gap')
        near=np.argmin(abs(z-(PREP+77*LU)))
        ratio=float(power[-1].sum()/max(power[near].sum(),1e-300))
        save(gen/'COMPLETE.json',dict(run=rec,lattice_end_m=endpoint,last_over_first_segment_power=ratio,
          gap_phase_kicks_verified=True,observed_autophase_rad=observed.tolist(),
          sufficient_gain_for_feedback_test=bool(ratio>10),physical_acceptance=False,
          scope='Independent reference integrity and gain report; no automatic physical precision certificate'))
if __name__=='__main__':main()
