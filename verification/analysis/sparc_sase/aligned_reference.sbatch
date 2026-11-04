#!/usr/bin/env bash
#SBATCH --job-name=SPARC-radiation-reference
#SBATCH --partition=ssmb
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=00:25:00
#SBATCH --no-requeue
set -Eeuo pipefail
review_root="$1"
source /public/home/wfollowm/APRL_Workspace/cluster.env
cd "$review_root/genesis-aligned-r1"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 HDF5_USE_FILE_LOCKING=FALSE
export OMPI_MCA_btl=self,vader,tcp MPLCONFIGDIR="$review_root/mpl-cache-aligned"
timeout --signal=TERM --kill-after=20 700 srun --mpi=pmix_v3 --exact --nodes=1 --ntasks=64 \
  --ntasks-per-node=64 --cpus-per-task=1 --cpu-bind=cores --kill-on-bad-exit=1 --wait=30 \
  /public/software/Genesis-1.3-Version4-4.6.4/build/genesis4 reference.in > genesis.log 2>&1
printf '0\n' > TRACK_COMPLETE
timeout --signal=TERM --kill-after=20 700 srun --mpi=none --exact --nodes=1 --ntasks=1 \
  --cpus-per-task=8 --cpu-bind=cores "$PYTHON_BIN" "$review_root/aligned-tools/aligned_reference.py" analyze \
  /public/home/wfollowm/APRL_Workspace/results/sparc_sase_comparison_20261103/compact-production-r1 \
  "$review_root" "$review_root/genesis-aligned-r1" \
  /public/home/wfollowm/APRL_Workspace/source/aprl > comparison.log 2>&1
