"""Short undulator/drift replay across two nodes, original sixteen ranks."""
import argparse,json,shutil
from pathlib import Path
from prepare import sha,save,require
def main():
 p=argparse.ArgumentParser();p.add_argument('--workspace',required=True,type=Path);p.add_argument('--build-record',required=True,type=Path);p.add_argument('--output',required=True,type=Path);a=p.parse_args();w=a.workspace.resolve();b=a.build_record.resolve();r=a.output.resolve();require((b/'EXIT_CODE').read_text().strip()=='0','Closed unit job');program=Path(json.loads((b/'PLAN.json').read_text())['programs']['aprl']);require(sha(program)==json.loads((b/'SOURCE_BUILD.json').read_text())['program_sha256']['aprl'],'Program identity');r.mkdir(parents=True,exist_ok=False);t=Path(__file__).resolve().parent;baseline=w/'results/lab_energy_balance_20261102/coupled_sbp_r2';scripts=r/'frozen_source';shutil.copytree(w/'results/lab_coupling_algorithm_20261102/r3/frozen_source',scripts,ignore=shutil.ignore_patterns('__pycache__','*.pyc'));shutil.copy2(t/'run_flow.py',scripts/'run_flow.py')
 q=scripts/'post.py';text=q.read_text();require('candidate_magnetic_gather_changed=True' in text,'Expected prior declared algorithm-change reader');q.write_text(text.replace('candidate_magnetic_gather_changed=True','candidate_magnetic_gather_changed=True, candidate_electric_gather_changed=True, matched_Cowan_metric_ledger=False'))
 cases=[]
 for name in ['undulator-on','drift-on']:
  d=r/name;d.mkdir();(d/'run.yaml').write_text((baseline/name/'run.yaml').read_text().replace(str(baseline),str(r)))
  for sub in ['detectors','ledger','trajectory','lab-audit','lab-audit/synchronized','momentum']:(d/sub).mkdir()
  cases.append(dict(name=name,directory=str(d),timeout_seconds=600))
 (r/'bin').mkdir();(r/'bin/replay').symlink_to(w/'results/lab_energy_balance_20261102/r5/bin/replay');save(r/'PLAN.json',dict(workspace=str(w),baseline=str(baseline),output=str(r),programs={'wave':str(program)},binding_probe=str(w/'build/cluster-tools/binding_probe'),cases=cases,algorithm='Cowan-metric-E-and-momentum-adjoint-B',candidate_algorithm_changed=True,production_program_replaced=False,mpi_ranks=16,omp_threads=8,nodes=2,metric_field_diagnostic='Existing unweighted reference decomposition retained; not a new full moving-boundary matched-metric ledger'))
 script='''#!/usr/bin/env bash
#SBATCH --job-name=cowan-matched-flow
#SBATCH --partition=ssmb
#SBATCH --nodes=2
#SBATCH --ntasks=16
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=8
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=00:25:00
#SBATCH --no-requeue
set -Eeuo pipefail
source WORKSPACE/cluster.env
cd RESULTROOT
finish() { code=$?; trap - EXIT; printf '%s\\n' "$code" > EXIT_CODE; exit "$code"; }
trap finish EXIT
export OMP_NUM_THREADS=8 OMP_DYNAMIC=FALSE OMP_PROC_BIND=close OMP_PLACES=cores
export OMPI_MCA_pml=ob1 OMPI_MCA_btl=self,vader,tcp OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
timeout --signal=TERM --kill-after=20 1440 "$PYTHON_BIN" frozen_source/run_flow.py > driver.log 2>&1
'''.replace('WORKSPACE',str(w)).replace('RESULTROOT',str(r));(r/'run.sbatch').write_text(script);save(r/'FROZEN.json',{str(f):sha(f) for f in [*scripts.iterdir(),program,r/'PLAN.json',r/'run.sbatch',*[r/x['name']/'run.yaml' for x in cases]] if f.is_file()});print(r/'run.sbatch')
if __name__=='__main__':main()
