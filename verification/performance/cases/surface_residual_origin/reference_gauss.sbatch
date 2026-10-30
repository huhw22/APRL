#!/usr/bin/env bash
#SBATCH --job-name=aprl-ref-gauss
#SBATCH --partition=ssmb
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=00:04:00
#SBATCH --no-requeue
set -Eeuo pipefail
root=$(realpath "${1:?Main result root}")
build=$(realpath "${2:?Independent diagnostic build}")
source /public/home/wfollowm/APRL_Workspace/cluster.env
export OMP_NUM_THREADS=8 OMP_DYNAMIC=FALSE OMP_PLACES=cores OMP_PROC_BIND=close OMPI_MCA_btl=self,vader,tcp
export LD_LIBRARY_PATH="$build/runtime-libs:${LD_LIBRARY_PATH:-}"
timeout --signal=TERM --kill-after=10 180 srun --mpi=pmix_v3 --ntasks=1 --cpus-per-task=8 --cpu-bind=cores "$build/reference_gauss" "$root/base/reference" "$root/reference_gauss.csv" 0.22949549820217161 0.2292034717621822 1677 150
