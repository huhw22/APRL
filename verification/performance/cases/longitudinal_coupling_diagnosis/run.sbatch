#!/usr/bin/env bash
#SBATCH --job-name=longitudinal-coupling-check
#SBATCH --partition=ssmb
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:15:00
#SBATCH --no-requeue
set -Eeuo pipefail
source WORKSPACE/cluster.env
cd RESULTROOT
finish() { code=$?; trap - EXIT; printf '%s\n' "$code" > EXIT_CODE; exit "$code"; }
trap finish EXIT
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export OMPI_MCA_pml=ob1 OMPI_MCA_btl=self,vader,tcp
timeout --signal=TERM --kill-after=20 840 "$PYTHON_BIN" frozen_source/run.py > driver.log 2>&1
