#!/usr/bin/env bash
#SBATCH --job-name=aprl-pending-unit
#SBATCH --partition=ssmb
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=8
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=00:15:00
#SBATCH --no-requeue
set -Eeuo pipefail
build=$(realpath "${1:?candidate build}")
root=$(realpath "${2:?new result directory}")
source /public/home/wfollowm/APRL_Workspace/cluster.env
[[ ! -e "$root/STARTED" ]]
date --iso-8601=seconds > "$root/STARTED"
trap 'code=$?; printf "%s\n" "$code" > "$root/EXIT_CODE"' EXIT
trap 'exit 143' TERM
trap 'exit 130' INT
export LD_LIBRARY_PATH="$build/runtime-libs${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export OMPI_MCA_btl=self,vader,tcp OMPI_MCA_btl_openib_warn_no_device_params_found=0
export OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE OMP_PLACES=cores OMP_PROC_BIND=close
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
ulimit -s unlimited
for check in aprl_pending_particle_queue_checks aprl_pending_particle_statistics_checks; do
  timeout --signal=TERM --kill-after=10 120 srun --mpi=none --exact --nodes=1 --ntasks=1 --cpus-per-task=1 --cpu-bind=cores \
    "$build/verification/$check" > "$root/$check.log" 2>&1
done
cd "$build"
timeout --signal=TERM --kill-after=15 600 ctest -L required --output-on-failure -j 1 > "$root/required.log" 2>&1
printf 'PASS\n' > "$root/PASS"
