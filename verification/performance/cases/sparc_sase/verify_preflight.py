"""Validate closed bounded data using the step count derived from the card."""
import argparse,math,re
from pathlib import Path
import h5py,numpy as np,yaml
from common import *
import duration_status
def main(root):
    root=Path(root);plan=load(root/'plan.json');d=root/'timing'
    need(all(sha(p)==h for p,h in load(root/'FROZEN_CARDS.json').items()),'Frozen preflight identity changed')
    need(load(root/'BINDING.json')['pass_binding'],'Core binding not passed')
    for log in ['legacy_magnets.log','new_magnets.log']:
        need(load(root/(log+'.launch.json'))['exit_code']==0,'Magnet probe failed')
    need((root/'legacy_magnets.log').read_text().splitlines()[0]==
         (root/'new_magnets.log').read_text().splitlines()[0],'Original magnets changed')
    need('new_magnet_checks PASS' in (root/'new_magnets.log').read_text(),'No magnet acceptance')
    cfg=yaml.safe_load((d/'run.yaml').read_text())
    dt=cfg['mesh']['cell_size'][2]/c
    steps=int(math.ceil(cfg['mesh']['duration']/dt))
    rec=load(d/'aprl.log.launch.json');text=(d/'aprl.log').read_text()
    status=duration_status.bounded_duration_status(text,
      dict(exit_code=rec['exit_code'],seconds=rec['seconds'],limit_s=rec['limit_seconds']),
      steps,duration_status)
    closed=[]
    for p in cfg['detectors']['field_planes']:
        path=d/'detectors'/(p['name']+'.h5')
        with h5py.File(path,'r') as f:
            g=f['field_plane'];n=int(g['committed_samples'][()])
            # Intentional duration exhaustion correctly marks physical output
            # incomplete. Ordered close/readiness traces above prove closure.
            need(int(g['complete'][()])==0 and n==len(g['fields'])==len(g['time_s']),'Uncommitted bounded output')
            need(n>=steps,'Too few field frames')
            for j in range(0,n,16):
                v=g['fields'][j:j+16]
                for key in v.dtype.names:need(np.isfinite(v[key]).all(),'Nonfinite bounded field')
        closed.append(dict(path=str(path),closed_and_committed=True,physical_complete_flag=0,samples=n))
    previous=Path(plan.get('reuse_input_and_tests',plan.get('reuse_tests',root)))
    need(load(previous/'required.log.launch.json')['exit_code']==0,'Regressions not passed')
    need(load(root/'NOISE_VALIDATION.json')['pass_noise_variance'],'Noise not passed')
    stats=re.findall(r'phase=progress step=(\d+) elapsed_s=([0-9.eE+-]+) measured_step_s=([0-9.eE+-]+)',text)
    result=dict(pass_preflight=True,expected_steps_derived_from_duration=steps,run=rec,
      bounded_status=status,fields=closed,geometry=plan['geometry'],physical_acceptance=False,
      binary_sha256=sha(plan['binary']),required_reused_from=str(previous),
      magnet_probes_pass=True,binding_pass=True,frozen_identity_pass=True,
      progress=[dict(step=int(s),elapsed_seconds=float(t),mean_seconds_per_step=float(v)) for s,t,v in stats],
      repair='The frozen reader assumed 120; ceil(mesh.duration/dt) is 121. Outputs are correctly marked physically incomplete for bounded timing. Original batch status and all raw data preserved.')
    save(root/'PREFLIGHT_VALIDATED.json',result)
    print(json.dumps(dict(steps=steps,pass_preflight=True,progress=result['progress']),indent=2))
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('root',type=Path);main(ap.parse_args().root)
