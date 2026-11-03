"""Resume only full-geometry timing after closed input/regression preparation."""
import os
from pathlib import Path
import h5py,numpy as np,yaml
from common import *
from launch import run,srun
def main():
    root=Path(__file__).resolve().parents[1];plan=load(root/'plan.json')
    need(load(root/'NOISE_VALIDATION.json')['pass_noise_variance'],'Noise has not passed')
    need(all(sha(p)==h for p,h in load(root/'FROZEN_CARDS.json').items()),'Prepared card changed')
    reused=Path(plan['reuse_input_and_tests'])
    need(load(reused/'required.log.launch.json')['exit_code']==0,'Required tests missing')
    env=os.environ.copy();env['APRL_TEARDOWN_TRACE']='1'
    run(srun()+[str(WS/'build/cluster-tools/binding_probe'),'32','8'],root,90,'binding.log',env)
    d=root/'case'
    rec=run(srun()+[plan['binary'],str(d/'run.yaml')],d,1800,'aprl.log',env,
            expected_duration_steps=120)
    cfg=yaml.safe_load((d/'run.yaml').read_text());closed=[]
    for p in cfg['detectors']['field_planes']:
        path=d/'detectors'/(p['name']+'.h5')
        with h5py.File(path,'r') as f:
            g=f['field_plane'];n=int(g['committed_samples'][()])
            need(int(g['complete'][()])==1 and n==len(g['fields']),'Unclosed bounded output')
            for start in range(0,n,16):
                a=g['fields'][start:start+16]
                for name in a.dtype.names:need(np.isfinite(a[name]).all(),'Nonfinite bounded field')
        closed.append(dict(path=str(path),samples=n))
    save(root/'PREFLIGHT_COMPLETE.json',dict(pass_preflight=True,binary_sha256=sha(plan['binary']),
      required=load(reused/'required.log.launch.json'),timing=rec,
      fields=closed,geometry=plan['geometry'],physical_acceptance=False,
      source_input_campaign=plan['root'],scope='Complete 120 field steps and closed finite outputs; full physical stop intentionally not reached'))
    print('PREFLIGHT_COMPLETE',flush=True)
if __name__=='__main__':main()
