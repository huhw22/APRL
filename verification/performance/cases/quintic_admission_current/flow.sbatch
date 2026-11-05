#!/usr/bin/env bash
#SBATCH --job-name=aprl-admission-flow
#SBATCH --partition=ssmb
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=8
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=00:25:00
#SBATCH --no-requeue
set -Eeuo pipefail
root=$(realpath "${1:?Prepared short-flow root}")
candidate_record=${2:-}
source /public/home/wfollowm/APRL_Workspace/cluster.env
finish() { code=$?; trap - EXIT; printf '%s\n' "$code" > "$root/EXIT_CODE"; exit "$code"; }
trap finish EXIT
export OMP_NUM_THREADS=8 OMP_DYNAMIC=FALSE OMP_PLACES=cores OMP_PROC_BIND=close
export OMPI_MCA_pml=ob1 OMPI_MCA_btl=self,vader,tcp HDF5_USE_FILE_LOCKING=FALSE
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 APRL_PROFILE_INTERVAL=0
ulimit -s unlimited
ulimit -c 0
args=("$root")
if [[ -n "$candidate_record" ]]; then args+=("$(realpath "$candidate_record")"); fi
timeout --signal=TERM --kill-after=20 1440 "$PYTHON_BIN" "$root/scripts/flow_run.py" \
  "${args[@]}" > "$root/driver.log" 2>&1
