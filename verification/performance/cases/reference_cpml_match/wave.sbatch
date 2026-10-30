#!/usr/bin/env bash
#SBATCH --job-name=aprl-cpml-wave
#SBATCH --partition=ssmb
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=00:03:00
#SBATCH --no-requeue
set -Eeuo pipefail
root=$(realpath "${1:?Result root}")
build=$(realpath "${2:?Independent build}")
source /public/home/wfollowm/APRL_Workspace/cluster.env
export OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE OMPI_MCA_btl=self,vader,tcp
export LD_LIBRARY_PATH="$build/runtime-libs:${LD_LIBRARY_PATH:-}"
timeout --signal=TERM --kill-after=10 120 srun --mpi=pmix_v3 --nodes=1 --ntasks=1 --cpus-per-task=1 --cpu-bind=cores "$build/wave_absorption_probe" "$root/wave"
