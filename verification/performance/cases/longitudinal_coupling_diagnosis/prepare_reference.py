import argparse,json,shutil
from pathlib import Path
from prepare import sha,save,require
def main():
 p=argparse.ArgumentParser();p.add_argument('--workspace',required=True,type=Path);p.add_argument('--input',required=True,type=Path);p.add_argument('--output',required=True,type=Path);a=p.parse_args();w=a.workspace.resolve();data=a.input.resolve();r=a.output.resolve();require((data/'EXIT_CODE').read_text().strip()=='0','Individual controls have not closed');r.mkdir(parents=True,exist_ok=False);f=r/'frozen_source';f.mkdir();t=Path(__file__).resolve().parent
 for name in ['fluid_reference.py','prepare.py']:shutil.copy2(t/name,f/name)
 script='''#!/usr/bin/env bash
#SBATCH --job-name=longitudinal-fluid-reference
#SBATCH --partition=ssmb
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:08:00
#SBATCH --no-requeue
set -Eeuo pipefail
source WORKSPACE/cluster.env
cd RESULTROOT
finish() { code=$?; trap - EXIT; printf '%s\\n' "$code" > EXIT_CODE; exit "$code"; }
trap finish EXIT
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
timeout --signal=TERM --kill-after=20 420 "$PYTHON_BIN" frozen_source/fluid_reference.py --input INPUTROOT --output RESULTROOT/analysis > driver.log 2>&1
'''.replace('WORKSPACE',str(w)).replace('RESULTROOT',str(r)).replace('INPUTROOT',str(data));(r/'run.sbatch').write_text(script);files=[*f.iterdir(),r/'run.sbatch',*[data/(n+suffix) for n in ['analytic_linear','analytic_cubic_continuum','analytic_quintic_continuum'] for suffix in ['.json','.csv.particles.csv']]];save(r/'FROZEN.json',{str(q):sha(q) for q in files});save(r/'PLAN.json',dict(input=str(data),scope='Read-only independent cold-fluid reference for individual longitudinal trajectories',production_modified=False));print(r/'run.sbatch')
if __name__=='__main__':main()
