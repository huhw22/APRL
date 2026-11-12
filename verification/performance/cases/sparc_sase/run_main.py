import os,sys,time,traceback
from pathlib import Path
import h5py,numpy as np,yaml
from common import *
from launch import run,srun
from loading_guard import validate_loading
def check_fields(directory,cfg):
    records=[]
    for plane in cfg['detectors']['field_planes']:
        p=directory/'detectors'/(plane['name']+'.h5')
        with h5py.File(p,'r') as f:
            g=f['field_plane'];count=int(g['committed_samples'][()])
            need(int(g['complete'][()])==1 and count==len(g['fields']),'Unclosed field output '+str(p))
        records.append(dict(path=str(p),samples=count,complete=True))
    return records
def main():
    root=Path(__file__).resolve().parents[1];plan=load(root/'plan.json')
    frozen=load(root/'FROZEN.json');need(all(sha(p)==h for p,h in frozen.items()),'Frozen identity changed')
    index=int(sys.argv[1]);need(0<=index<len(plan['cases']),'Invalid array index')
    case=plan['cases'][index];d=Path(case['directory'])
    cfg=yaml.safe_load((d/'run.yaml').read_text())
    validate_loading(cfg['beam']['input']['file'], report=d/'PIC_LOADING_AUDIT.json')
    ranks=plan['execution']['ranks']
    run(srun()+[str(WS/'build/cluster-tools/binding_probe'),str(ranks),'8'],d,90,'binding.log')
    marker=Path(plan['preflight'])/'PREFLIGHT_VALIDATED.json'
    if not marker.exists():marker=Path(plan['preflight'])/'PREFLIGHT_COMPLETE.json'
    need(load(marker)['pass_preflight'],'No passed preflight')
    # The independent GENESIS reference has its own bounded allocation.
    completed=[]
    for case in [case]:
        d=Path(case['directory']);print('PIC START',case['name'],time.time(),flush=True)
        rec=run(srun()+[plan['binary'],str(d/'run.yaml')],d,case['timeout_seconds'],'aprl.log')
        cfg=yaml.safe_load((d/'run.yaml').read_text());closed=check_fields(d,cfg)
        save(d/'PIC_COMPLETE.json',dict(run=rec,fields=closed,physical_acceptance=False))
        need(all(sha(p)==h for p,h in frozen.items()),'Frozen files changed')
    for case in plan['cases']:
        p=Path(case['directory'])/'PIC_COMPLETE.json'
        if p.exists():completed.append(dict(case=case['name'],**load(p)))
    save(root/'STATUS.json',dict(completed=completed,physical_acceptance=False))
    if len(completed)==len(plan['cases']):
        save(root/'COMPLETE.json',dict(cases=completed,physical_acceptance=False,analysis_pending=True))
if __name__=='__main__':main()
