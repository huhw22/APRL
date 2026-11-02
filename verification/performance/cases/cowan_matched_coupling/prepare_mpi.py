"""Freeze the real four-rank, two-node weighted gather check."""
import argparse,json,shutil
from pathlib import Path
from prepare import sha,save,require
def main():
 p=argparse.ArgumentParser();p.add_argument('--workspace',required=True,type=Path);p.add_argument('--build-record',required=True,type=Path);p.add_argument('--output',required=True,type=Path);a=p.parse_args();w=a.workspace.resolve();b=a.build_record.resolve();r=a.output.resolve();require((b/'EXIT_CODE').read_text().strip()=='0','Closed unit job');program=Path(json.loads((b/'PLAN.json').read_text())['programs']['checks']);require(sha(program)==json.loads((b/'SOURCE_BUILD.json').read_text())['program_sha256']['checks'],'Program identity');r.mkdir(parents=True,exist_ok=False)
 script='''#!/usr/bin/env bash
#SBATCH --job-name=cowan-matched-mpi
#SBATCH --partition=ssmb
#SBATCH --nodes=2
#SBATCH --ntasks=4
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=00:04:00
#SBATCH --no-requeue
set -Eeuo pipefail
source WORKSPACE/cluster.env
cd RESULTROOT
finish() { code=$?; trap - EXIT; printf '%s\\n' "$code" > EXIT_CODE; exit "$code"; }
trap finish EXIT
test "$(sha256sum PROGRAM | cut -d ' ' -f 1)" = PROGRAMSHA
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMPI_MCA_pml=ob1 OMPI_MCA_btl=self,vader,tcp
timeout --signal=TERM --kill-after=20 120 srun --mpi=pmix_v3 --exact --nodes=2 --ntasks=4 --ntasks-per-node=2 --cpus-per-task=1 --cpu-bind=cores --kill-on-bad-exit=1 --wait=30 PROGRAM > checks.log 2>&1
'''.replace('WORKSPACE',str(w)).replace('RESULTROOT',str(r)).replace('PROGRAMSHA',sha(program)).replace('PROGRAM',str(program));(r/'run.sbatch').write_text(script);save(r/'FROZEN.json',{str(program):sha(program),str(r/'run.sbatch'):sha(r/'run.sbatch')});save(r/'PLAN.json',dict(program=str(program),program_sha256=sha(program),scope='Cross-node actual MPI collar and weighted deposition-adjoint identity, not physical accuracy acceptance',ranks=4,nodes=2));print(r/'run.sbatch')
if __name__=='__main__':main()
