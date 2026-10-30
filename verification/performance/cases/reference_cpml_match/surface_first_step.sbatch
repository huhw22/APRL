#!/usr/bin/env bash
#SBATCH --job-name=aprl-surface-first
#SBATCH --partition=ssmb
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=00:03:00
#SBATCH --no-requeue
set -Eeuo pipefail
root=$(realpath "${1:?Result root}")
build=$(realpath "${2:?Independent build}")
source /public/home/wfollowm/APRL_Workspace/cluster.env
export OMP_NUM_THREADS=8 OMP_DYNAMIC=FALSE OMP_PLACES=cores OMP_PROC_BIND=close
export OMPI_MCA_btl=self,vader,tcp
export LD_LIBRARY_PATH="$build/runtime-libs:${LD_LIBRARY_PATH:-}"
timeout --signal=TERM --kill-after=5 90 srun --mpi=pmix_v3 --ntasks=1 --cpus-per-task=8 --cpu-bind=cores "$build/surface_first_step_probe" "$root/surface_500/reference" "$root/first_step/surface.csv" 160 10 0.22936572471392327 3328 -0.0020871988659157434 10 "$root/surface_500/reference/transverse-current.bin" 0.22914945799500638
printf '1\n' > "$root/first_step/COMPLETE"
