"""Run independent native cases, then analyze only closed binary outputs."""
import json
import os
import subprocess
import sys
import time
from pathlib import Path
from prepare import sha,save

def main(root):
    root=Path(root);p=json.loads((root/'PLAN.json').read_text())
    assert os.environ.get('SLURM_JOB_ID'),'Compute-node Slurm allocation required'
    for path,h in json.loads((root/'FROZEN.json').read_text()).items():assert sha(path)==h
    for path,h in p['source_hashes'].items():assert sha(path)==h
    scale=p['scale'];closed=[]
    for c in p['cases']:
        d=root/'cases'/c['name'];t=time.monotonic()
        cmd=['mpirun','-np','4','--bind-to','core','--map-by','core',str(root/'bin/native_probe'),str(d),str(c['nx']),str(c['ny']),str(c['beta']),str(c['steps']),str(c['stride']),c['mode'],str(scale['gain_length_m']),str(scale['cooperation_length_m']),str(scale['transverse_scale_m'])]
        print('START',c['name'],flush=True)
        with (d/'run.log').open('w') as f:
            q=subprocess.run(['timeout','--signal=TERM','--kill-after=60','180']+cmd,stdout=f,stderr=subprocess.STDOUT)
        elapsed=time.monotonic()-t
        save(d/'LAUNCH.json',dict(command=cmd,returncode=q.returncode,elapsed_s=elapsed))
        assert q.returncode==0,(c['name'],q.returncode)
        complete=json.loads((d/'COMPLETE.json').read_text());assert complete['complete'] and complete['steps']==c['steps']
        nframes=c['steps']//c['stride']+1
        assert complete['frames']==nframes
        assert (d/'field_xy.bin').stat().st_size==nframes*81*81*16
        assert (d/'line_x.bin').stat().st_size==(c['steps']+1)*81*8
        assert len((d/'history.csv').read_text().splitlines())==c['steps']+2
        if c['mode']=='plane':assert json.loads((d/'PLANE_CHECK.json').read_text())['complete']
        closed.append(dict(name=c['name'],elapsed_s=elapsed,complete=complete));print('CLOSED',closed[-1],flush=True)
    for path,h in p['source_hashes'].items():assert sha(path)==h
    save(root/'PHYSICS_COMPLETE.json',dict(complete=True,job_id=os.environ['SLURM_JOB_ID'],cases=closed,source_hashes_unchanged=True))
    subprocess.run([sys.executable,'-B',str(root/'scripts/analyze.py'),str(root)],check=True)
    subprocess.run([sys.executable,'-B',str(root/'scripts/render.py'),str(root)],check=True)
    save(root/'ANALYSIS_COMPLETE.json',dict(complete=True,awaiting_manual_visual_review=True))

if __name__=='__main__':main(sys.argv[1])
