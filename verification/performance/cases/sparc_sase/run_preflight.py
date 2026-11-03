import json,math,os,re,subprocess,sys,traceback
from pathlib import Path
import h5py
from common import *
from launch import run,srun
def main():
    root=Path(__file__).resolve().parents[1];plan=load(root/'plan.json')
    need(all(sha(p)==x for p,x in load(root/'PREPARED.json').items()),'Prepared identity changed')
    env=os.environ.copy()
    ranks=plan['execution']['ranks'];cores=ranks*8
    rec=run(srun()+[str(WS/'build/cluster-tools/binding_probe'),str(ranks),'8'],root,90,'binding.log')
    rows=re.findall(r'\[binding\] rank=(\d+) thread=(\d+) host=(\S+) cpu=(\d+) package=(\d+) core=(\d+) affinity=(\S+) status=(\S+)',(root/'binding.log').read_text())
    need(len(rows)==cores and all(r[-1]=='PASS' for r in rows),'Binding failed')
    need(len({(r[2],r[4],r[5]) for r in rows})==cores,'Core overlap')
    save(root/'BINDING.json',dict(pass_binding=True,nodes=sorted({r[2] for r in rows}),cores=cores))
    old=run([str(WS/'build/sparc-sase-r1-20261103/legacy_magnet_checks')],root,60,'legacy_magnets.log')
    new=run([str(WS/'build/sparc-sase-r1-20261103/sparc_magnet_checks')],root,60,'new_magnets.log')
    need((root/'legacy_magnets.log').read_text().splitlines()[0]==(root/'new_magnets.log').read_text().splitlines()[0],
         'Legacy prescribed fields changed')
    need('new_magnet_checks PASS' in (root/'new_magnets.log').read_text(),'New magnet checks failed')
    if plan.get('reuse_input_and_tests'):
        previous=Path(plan['reuse_input_and_tests'])
        tests=dict(load(previous/'required.log.launch.json'),reused_from=str(previous))
        native=dict(load(previous/'native-noise/initial.log.launch.json'),reused_from=str(previous))
        need(tests['exit_code']==0 and native['exit_code']==0,'Cached preflight is incomplete')
        need(sha(root/'native-noise/initial.par.h5')==plan['reused_native_input_sha256'],'Native input changed')
        save(root/'REUSED_CHECKS.json',dict(required=tests,native=native))
    else:
        if plan.get('reuse_tests'):
            previous=Path(plan['reuse_tests'])
            tests=dict(load(previous/'required.log.launch.json'),reused_from=str(previous))
            need(tests['exit_code']==0,'Previous regressions failed')
            save(root/'REUSED_REQUIRED.json',tests)
        else:
            tests=run(['ctest','--test-dir',str(WS/'build/sparc-sase-r1-20261103'),
                   '-L','required','--output-on-failure','-j','1'],root,1200,'required.log')
        genenv=env.copy();genenv['OMP_NUM_THREADS']='1'
        native=run(srun()+[plan['genesis'],'initial.in'],root/'native-noise',900,'initial.log',genenv)
    need((root/'native-noise/initial.par.h5').exists(),'No native noise dump')
    run([sys.executable,str(root/'scripts/convert_input.py'),str(root)],root,600,'convert.log')
    need(load(root/'NOISE_VALIDATION.json')['pass_noise_variance'],'Noise failed')
    plan=load(root/'plan.json')
    import yaml,numpy as np
    cfg=yaml.safe_load((root/'timing/run.yaml').read_text())
    expected_steps=int(math.ceil(cfg['mesh']['duration']/(cfg['mesh']['cell_size'][2]/c)))
    timing_env=env.copy();timing_env['APRL_TEARDOWN_TRACE']='1'
    timing=run(srun()+[plan['binary'],str(root/'timing/run.yaml')],root/'timing',1800,'aprl.log',timing_env,
               expected_duration_steps=expected_steps)
    for plane in cfg['detectors']['field_planes']:
        with h5py.File(root/'timing/detectors'/(plane['name']+'.h5'),'r') as f:
            g=f['field_plane'];n=int(g['committed_samples'][()])
            need(int(g['complete'][()])==0 and n==len(g['fields'])==len(g['time_s']),
                 'Uncommitted bounded field (physical completion must remain false)')
            for i in range(0,n,16):
                data=g['fields'][i:i+16]
                for key in data.dtype.names:need(np.isfinite(data[key]).all(),'Nonfinite bounded field')
    save(root/'PREFLIGHT_COMPLETE.json',dict(pass_preflight=True,binary_sha256=sha(plan['binary']),
      required=tests,noise=native,timing=timing,geometry=plan['geometry'],
      expected_bounded_steps=expected_steps,physical_acceptance=False,
      scope='Full input and planned field geometry; intentionally partial physical propagation'))
    print('PREFLIGHT_COMPLETE',flush=True)
if __name__=='__main__':main()
