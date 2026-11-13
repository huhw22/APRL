"""One bounded eleven-segment Puffin run, followed by closed-data analysis."""
import os
import sys
import traceback
from pathlib import Path
from common import W,read,save,need,frozen,srun,launch
from binding import binding


def main(root):
    root=Path(root);plan=read(root/'PLAN.json')
    need(bool(os.environ.get('SLURM_JOB_ID')),'Use Slurm')
    frozen(root/'FROZEN_PHYSICS.json')
    need(read(root/'PREPARATION_COMPLETE.json')['complete'],'Preparation not closed')
    execution=plan['execution'];ranks=execution['main_ranks'];nodes=execution['main_nodes']
    launch(srun(nodes,ranks)+[str(W/'build/cluster-tools/binding_probe'),str(ranks),'1'],root,root/'binding.log',120)
    save(root/'BINDING.json',binding(root/'binding.log',nodes,ranks,1))
    save(root/'STARTED.json',dict(job=os.environ['SLURM_JOB_ID'],ranks=ranks,nodes=nodes,
        independent_physical_source=True,native_import_bitwise=True,segments=11,slippage_multiple=4,
        only_Puffin=True,APRL_input_saved=True,field_route_validated=True))
    launch(srun(nodes,ranks)+[plan['binary'],'main.in'],root/'main',root/'main/puffin.log',plan['physics_timeout_s'],True)
    need('Puffin simulation completed' in (root/'main/puffin.log').read_text(),'Physics did not normally close')
    save(root/'TRACKING_CLOSED.json',dict(complete=True,job=os.environ['SLURM_JOB_ID']))
    frozen(root/'FROZEN_PHYSICS.json')
    launch([sys.executable,'-B',str(root/'scripts/analyze.py'),str(root)],root,root/'analysis.log',plan['analysis_timeout_s'])
    save(root/'COMPLETE.json',dict(complete=True,job=os.environ['SLURM_JOB_ID'],segments=11,
        slippage_multiple=4,only_Puffin=True,APRL_input_saved=True,physical_acceptance=False))


if __name__=='__main__':
    try: main(sys.argv[1])
    except Exception:
        save(Path(sys.argv[1])/'FAILED.json',dict(error=traceback.format_exc(),job=os.environ.get('SLURM_JOB_ID')))
        raise
