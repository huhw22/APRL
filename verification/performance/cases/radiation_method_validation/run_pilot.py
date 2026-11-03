import json
import os
from pathlib import Path
from prepare import save,require,sha
from launcher import launch

def main():
    root=Path(__file__).resolve().parents[1]
    plan=json.loads((root/'pilot-plan.json').read_text())
    frozen=json.loads((root/'FROZEN.json').read_text())
    require(all(sha(p)==h for p,h in frozen.items()),'Frozen files changed')
    env=os.environ.copy()
    from check_retarded import verify
    save(root/'RETARDED_SOLVER_CHECK.json',verify())
    # One standalone MPI process avoids idle collective ranks for steady state.
    # Input/output remain inside the allocated compute node.
    results=[]
    for case in plan['genesis_design']:
        directory=Path(case['directory'])
        try:
            rec=launch([plan['genesis_binary'],'design.in'],directory,env,180,'genesis.log')
            results.append({**case,**rec,'complete':True})
        except Exception as exc:
            results.append({**case,'complete':False,'error':str(exc)})
        save(root/'GENESIS_DESIGN_STATUS.json',results)
    from analyze_dipole import analyze
    analyze(root)
    save(root/'PILOT_COMPLETE.json',{'complete':True,'three_method_acceptance':False})

if __name__=='__main__':main()
