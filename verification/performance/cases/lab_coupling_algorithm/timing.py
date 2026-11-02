"""Same-allocation ABBA short-card timing; not a production speed estimate."""
import argparse,json,os,re,shutil,subprocess
from pathlib import Path
import numpy as np,yaml
from prepare import require,save,sha
def prepare(w,r):
    r.mkdir(parents=True,exist_ok=False);frozen=r/'frozen_source';frozen.mkdir();t=Path(__file__).resolve().parent;helper=t.parent/'coupled_lab_sbp_audit'
    for p in [t/'timing.py',t/'prepare.py',helper/'launcher.py']:shutil.copy2(p,frozen/p.name)
    science=w/'results/lab_coupling_algorithm_20261102/r3';require((science/'EXIT_CODE').read_text().strip()=='0','Candidate science incomplete');cases=[]
    for number,variant in enumerate(['baseline','candidate','candidate','baseline']):
        d=r/(str(number+1)+'-'+variant);d.mkdir();card=yaml.safe_load((science/'undulator-on/run.yaml').read_text().replace(str(science/'undulator-on'),str(d)));card['energy_ledger']['enabled']=False;card['energy_ledger']['discrete_operator_audit']=False;(d/'run.yaml').write_text(yaml.safe_dump(card,sort_keys=False));(d/'detectors').mkdir();cases.append(dict(name=d.name,variant=variant))
    p=dict(workspace=str(w),output=str(r),programs={'baseline':str(w/'build/coupled-lab-sbp-r2-20261102/aprl'),'candidate':str(w/'build/lab-coupling-algorithm-r3-20261102/aprl')},binding_probe=str(w/'build/cluster-tools/binding_probe'),cases=cases,physics_science_gate=str(science/'RESULT.json'),resources=dict(nodes=1,mpi=16,threads=8),heavy_ledgers_and_observers=False,detectors_identical=True);save(r/'PLAN.json',p)
    batch='''#!/usr/bin/env bash
#SBATCH --job-name=lab-gather-ABBA
#SBATCH --partition=ssmb
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=8
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=00:20:00
#SBATCH --no-requeue
set -Eeuo pipefail
source WORKSPACE/cluster.env
cd ROOTPATH
finish() { code=$?; trap - EXIT; printf '%s\\n' "$code" > EXIT_CODE; exit "$code"; }
trap finish EXIT
export OMP_NUM_THREADS=8 OMP_DYNAMIC=FALSE OMP_PROC_BIND=close OMP_PLACES=cores
export OMPI_MCA_btl=self,vader,tcp OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
timeout --signal=TERM --kill-after=20 1140 "$PYTHON_BIN" frozen_source/timing.py > driver.log 2>&1
'''.replace('WORKSPACE',str(w)).replace('ROOTPATH',str(r));(r/'run.sbatch').write_text(batch)
    files=list(frozen.iterdir())+[r/'PLAN.json',r/'run.sbatch']+[r/c['name']/'run.yaml' for c in cases]+[Path(x) for x in p['programs'].values()];save(r/'FROZEN.json',{str(p):sha(p) for p in files});print(json.dumps(p,indent=2))
def run():
    from launcher import launch,binding
    r=Path(__file__).resolve().parents[1];p=json.loads((r/'PLAN.json').read_text());env=os.environ.copy()
    for k in list(env):
        if k.startswith(('APRL_LAB_FIELD','APRL_WHOLE_COWAN','APRL_PARTICLE_MOMENTUM','APRL_SOURCE_EXIT','APRL_PROFILE')):env.pop(k)
    for f,h in json.loads((r/'FROZEN.json').read_text()).items():require(sha(f)==h,'Frozen changed '+f)
    srun=['srun','--mpi=pmix_v3','--exact','--nodes=1','--ntasks=16','--ntasks-per-node=16','--cpus-per-task=8','--cpu-bind=cores','--kill-on-bad-exit=1','--wait=30'];launch(srun+[p['binding_probe'],'16','8'],r,env,45,'binding.log');save(r/'BINDING.json',binding(r/'binding.log',1));reports=[]
    for c in p['cases']:
        d=r/c['name'];print('START',c['name'],flush=True);launch(srun+[p['programs'][c['variant']],str(d/'run.yaml')],d,env,270,'aprl.log');text=(d/'aprl.log').read_text();pairs=re.findall(r'phase=progress step=(\d+) elapsed_s=([\d.eE+\-]+)',text);points={int(k):float(v) for k,v in pairs};require(80 in points and 180 in points,'Missing stable timing interval');require('phase=final step=257' in text,'Unexpected endpoint');value=dict(name=c['name'],variant=c['variant'],steady_80_to_180_seconds_per_step=(points[180]-points[80])/100,launch_seconds=json.loads((d/'aprl.log.launch.json').read_text())['seconds']);reports.append(value);print(json.dumps(value),flush=True)
    A=[x['steady_80_to_180_seconds_per_step'] for x in reports if x['variant']=='baseline'];B=[x['steady_80_to_180_seconds_per_step'] for x in reports if x['variant']=='candidate'];save(r/'RESULT.json',dict(runs=reports,baseline_mean=float(np.mean(A)),candidate_mean=float(np.mean(B)),candidate_over_baseline=float(np.mean(B)/np.mean(A)),change_percent=float(100*(np.mean(B)/np.mean(A)-1)),production_speedup_claim=False,scope='one-node short undulator, same-allocation ABBA, heavy ledgers and observers disabled, detector settings identical'))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--workspace',type=Path);p.add_argument('--output',type=Path);a=p.parse_args();prepare(a.workspace.resolve(),a.output.resolve()) if a.output else run()
