#!/usr/bin/env bash
#SBATCH --job-name=aprl-pending-mpi
#SBATCH --partition=ssmb
#SBATCH --nodes=2
#SBATCH --ntasks=16
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=8
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=00:15:00
#SBATCH --no-requeue
set -Eeuo pipefail
root=$(realpath "${1:?Frozen two-node validation directory}")
source /public/home/wfollowm/APRL_Workspace/cluster.env
# The one-node candidate science gate is mandatory before any MPI step.
"$PYTHON_BIN" - "$root" <<'PY'
import json,sys
from pathlib import Path
r=Path(sys.argv[1]);p=json.loads((r/'plan.json').read_text())
gate=json.loads(Path(p['parent_science_pass']).read_text())
assert gate['pass'] and gate['binary_sha256']['candidate']==p['binary_sha256']
PY
[[ ${SLURM_JOB_NUM_NODES:-0} == 2 && ${SLURM_NTASKS:-0} == 16 && ${SLURM_CPUS_PER_TASK:-0} == 8 ]]
mkdir "$root/execution"
cd "$root/execution"
finish() { code=$?; trap - EXIT; printf '%s\n' "$code" > EXIT_CODE; date --iso-8601=seconds > ENDED_AT; exit "$code"; }
trap finish EXIT
trap 'exit 143' TERM
trap 'exit 130' INT
date --iso-8601=seconds > STARTED_AT
scontrol show job "$SLURM_JOB_ID" > allocation.txt
export OMP_NUM_THREADS=8 OMP_DYNAMIC=FALSE OMP_PLACES=cores OMP_PROC_BIND=close
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export HDF5_USE_FILE_LOCKING=FALSE APRL_PROFILE_INTERVAL=0
ulimit -s unlimited
# Two sequential candidate-only cases, each <=240 s; 720 s active budget.
# The unused half of each exclusive node remains reserved, not oversubscribed.
timeout --signal=TERM --kill-after=20 760 "$PYTHON_BIN" "$root/scripts/mpi_review.py" --root "$root" > driver.log 2>&1
