#!/usr/bin/env bash
#SBATCH --job-name=lab-algorithm-screen
#SBATCH --partition=ssmb
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:03:00
#SBATCH --no-requeue
set -Eeuo pipefail
source /public/home/wfollowm/APRL_Workspace/cluster.env
root=$1
cd "$root"
finish() { code=$?; trap - EXIT; printf '%s\n' "$code" > EXIT_CODE; exit "$code"; }
trap finish EXIT
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
timeout 150 "$PYTHON_BIN" frozen_source/spectral_screen.py "$root/output" > screen.log 2>&1
