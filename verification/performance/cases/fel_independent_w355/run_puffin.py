"""Run the validated w355 source with the unchanged qualified Puffin binary."""
import os
import sys
import traceback
from pathlib import Path
from common import W,read,save,need,frozen,srun,launch
from binding import binding


def main(root):
    root=Path(root).resolve();plan=read(root/'PLAN.json')
    need(bool(os.environ.get('SLURM_JOB_ID')),'Use Slurm')
    need(read(root.parent/'PREPARATION_COMPLETE.json')['complete'],'Campaign preparation did not close')
    need(read(root/'PREPARATION_COMPLETE.json')['complete'],'Source preparation did not close')
    frozen(root/'FROZEN_PHYSICS.json')
    ranks=plan['main_ranks'];nodes=plan['main_nodes']
    launch(srun(nodes,ranks)+[str(W/'build/cluster-tools/binding_probe'),str(ranks),'1'],root,root/'binding.log',120)
    save(root/'BINDING.json',binding(root/'binding.log',nodes,ranks,1))
    info=read(root/'INPUT_VALIDATED.json')
    save(root/'STARTED.json',dict(job=os.environ['SLURM_JOB_ID'],ranks=ranks,nodes=nodes,
        segments=plan['segments'],independent_physical_source=True,native_import_bitwise=True,
        actual_weight=info['actual_weight'],macro_records=info['records'],input_sha256=info['input_sha256'],
        binary_sha256=plan['original_binary_sha256'],field_route_validated=True,smooth_entrance_validated=True))
    launch(srun(nodes,ranks)+[plan['binary'],'main.in'],root/'main',root/'main/puffin.log',plan['physics_timeout_s'],True)
    need('Puffin simulation completed' in (root/'main/puffin.log').read_text(),'Physics did not normally close')
    save(root/'TRACKING_CLOSED.json',dict(complete=True,job=os.environ['SLURM_JOB_ID']))
    frozen(root/'FROZEN_PHYSICS.json')
    save(root/'COMPLETE.json',dict(complete=True,job=os.environ['SLURM_JOB_ID'],segments=plan['segments'],
        steps=plan['module_periods']*80,analysis_pending=True,physical_acceptance=False))


if __name__=='__main__':
    try:main(Path(sys.argv[1]))
    except Exception:
        save(Path(sys.argv[1])/'FAILED.json',dict(error=traceback.format_exc(),job=os.environ.get('SLURM_JOB_ID')))
        raise
