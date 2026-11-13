"""Validate actual import, then run exactly one six-segment historical control."""
import math
import os
import sys
import traceback
from pathlib import Path
import h5py
import numpy as np
from common import W,read,save,need,frozen,srun,launch,dump_at,sha
from binding import binding


def native(root,plan):
    with h5py.File(root/'particles.h5') as f:a=f['electrons'][:]
    with h5py.File(root/'import-check/main_electrons_0.h5') as f:b=f['electrons'][:]
    need(a.shape==b.shape and np.isfinite(b).all(),'Native population/nonfinite changed')
    # The historical input intentionally has repeated x,y. Full-row sorting
    # checks the multiset exactly without requiring a unique transverse key.
    ia=np.lexsort(tuple(a[:,j] for j in reversed(range(7))))
    ib=np.lexsort(tuple(b[:,j] for j in reversed(range(7))))
    need(np.array_equal(a[ia],b[ib]),'Actual Puffin import changed one of seven input columns')
    first=dump_at(root/'import-check','electrons',1)
    with h5py.File(first) as f:v=f['electrons'][:]
    weighted=lambda x:np.sum(x[:,3].astype(np.longdouble)*x[:,6].astype(np.longdouble))/np.sum(x[:,6],dtype=np.longdouble)
    delta=float((weighted(v)-weighted(a))*plan['K'])
    expected=-plan['K']*.5*(1-math.cos(math.pi/80))*math.cos(2*math.pi/80)
    need(abs(delta-expected)<1e-5,'Real first-step magnetic trajectory does not match smooth entrance')
    result=dict(passed=True,records=len(a),all_seven_columns_bitwise=True,matching='Exact full-row multiset, including repeated transverse states and unequal weights',
        actual_first_step_mean_delta_ux=delta,on_axis_expected_delta_ux=expected,
        field_feedback_disabled_for_import_check=True,input_native_sha256=sha(root/'particles.h5'))
    save(root/'NATIVE_VALIDATED.json',result)


def main(root):
    root=Path(root);plan=read(root/'PLAN.json');need(bool(os.environ.get('SLURM_JOB_ID')),'Use Slurm')
    frozen(root/'FROZEN.json');ranks=plan['main_ranks'];nodes=plan['main_nodes']
    launch(srun(nodes,ranks)+[str(W/'build/cluster-tools/binding_probe'),str(ranks),'1'],root,root/'binding.log',120)
    save(root/'BINDING.json',binding(root/'binding.log',nodes,ranks,1))
    launch(srun(nodes,ranks)+[plan['binary'],'main.in'],root/'import-check',root/'preparation/import.log',plan['preflight_timeout_s'],True)
    need('Puffin simulation completed' in (root/'preparation/import.log').read_text(),'Import check did not close normally')
    native(root,plan);frozen(root/'FROZEN.json')
    save(root/'PREPARATION_COMPLETE.json',dict(complete=True,job=os.environ['SLURM_JOB_ID'],native_import_bitwise=True,
        matched_magnet_qualified=True,source_sample_unchanged=True,only_Puffin=True))
    save(root/'STARTED.json',dict(job=os.environ['SLURM_JOB_ID'],ranks=ranks,nodes=nodes,segments=6,
        historical_input_preserved=True,independent_random_benchmark=False,no_new_PIC=True))
    launch(srun(nodes,ranks)+[plan['binary'],'main.in'],root/'main',root/'main/puffin.log',plan['physics_timeout_s'],True)
    need('Puffin simulation completed' in (root/'main/puffin.log').read_text(),'Production did not close normally')
    save(root/'TRACKING_CLOSED.json',dict(complete=True,job=os.environ['SLURM_JOB_ID']))
    frozen(root/'FROZEN.json')
    launch([sys.executable,'-B',str(root/'scripts/analyze.py'),str(root)],root,root/'analysis.log',plan['analysis_timeout_s'])
    save(root/'COMPLETE.json',dict(complete=True,job=os.environ['SLURM_JOB_ID'],segments=6,
        only_Puffin=True,physical_acceptance=False))


if __name__=='__main__':
    try:main(sys.argv[1])
    except Exception:
        save(Path(sys.argv[1])/'FAILED.json',dict(error=traceback.format_exc(),job=os.environ.get('SLURM_JOB_ID')))
        raise
