#!/bin/bash
#SBATCH --job-name=sparc-bandwidth-scan
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --exclusive
#SBATCH --mem=24G
#SBATCH --time=00:40:00
#SBATCH --output=bandwidth-scan-%j.log

set -euo pipefail
campaign_root="${1:?pass completed eleven-segment campaign directory}"
analysis_output="${2:?pass a new output directory}"
script_directory="${3:?pass this versioned bandwidth_scan directory}"
export OMP_NUM_THREADS=8
export OPENBLAS_NUM_THREADS=1
export MPLCONFIGDIR="${SLURM_TMPDIR:-/tmp}/matplotlib-bandwidth-${SLURM_JOB_ID}"
mkdir -p "$MPLCONFIGDIR"
/public/home/wfollowm/micromamba/envs/physics/bin/python3 \
  "$script_directory/analyze_bandwidth_scan.py" "$campaign_root" "$analysis_output"
