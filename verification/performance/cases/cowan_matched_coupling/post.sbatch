#!/usr/bin/env bash
#SBATCH --job-name=cowan-matched-readonly
#SBATCH --partition=ssmb
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:15:00
#SBATCH --no-requeue
set -Eeuo pipefail
source WORKSPACE/cluster.env
cd RESULTROOT
finish() { code=$?; trap - EXIT; printf '%s\n' "$code" > EXIT_CODE; exit "$code"; }
trap finish EXIT
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
timeout --signal=TERM --kill-after=20 840 "$PYTHON_BIN" frozen_source/post.py > driver.log 2>&1
timeout --signal=TERM --kill-after=20 30 "$PYTHON_BIN" frozen_source/review.py RESULTROOT > review.log 2>&1
