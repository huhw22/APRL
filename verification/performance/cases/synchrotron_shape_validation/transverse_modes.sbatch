#!/bin/bash
#SBATCH --job-name=synch-8eV-modes
#SBATCH --partition=ssmb
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=01:00:00
#SBATCH --no-requeue

set -euo pipefail

output="$1"
shift

export OMP_NUM_THREADS=16
export OPENBLAS_NUM_THREADS=16
export MKL_NUM_THREADS=16
export NUMEXPR_NUM_THREADS=16
export OMP_DYNAMIC=FALSE
export OMP_PLACES=cores
export OMP_PROC_BIND=close
export MPLCONFIGDIR="${SLURM_TMPDIR:-/tmp}/aprl-8eV-modes-mpl"

python_bin=/public/home/wfollowm/micromamba/envs/physics/bin/python
script_root=/public/home/wfollowm/APRL_Workspace/source/aprl/verification/performance/cases/synchrotron_shape_validation

"$python_bin" "$script_root/analyze_100k_transverse_single_frequency.py" \
  "$output" "$@"
"$python_bin" "$script_root/review_100k_bend_modes.py" \
  "$output/transverse_maps_8eV.npz" "$output"
