#!/bin/bash
#SBATCH --job-name=w355-independent-side-audit
#SBATCH --partition=ssmb
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=00:10:00
#SBATCH --no-requeue
set -Eeuo pipefail
export OMP_NUM_THREADS=8
audit_root="$1"
c++ -std=c++17 -O3 -fopenmp "$audit_root/scripts/side_check.cpp" -o "$audit_root/side_integrator"
timeout --signal=TERM --kill-after=30 500 "$audit_root/side_integrator" "$audit_root/../six-aprl/energy-momentum/worldtube" "$audit_root/../momentum-ledger-20261116-r1/lab/energy-momentum/tube-config.txt" "$audit_root/INDEPENDENT_SIDES.csv"
