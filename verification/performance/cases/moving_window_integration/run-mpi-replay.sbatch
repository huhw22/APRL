#!/bin/bash
#SBATCH --job-name=aprl-moving-mpi-replay
#SBATCH --partition=ssmb
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=128
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=00:05:00
#SBATCH --no-requeue
set -euo pipefail
[[ "$#" -gt 0 ]]
for campaign in "$@"; do
  root=$(realpath "$campaign")
  bash "$root/run-two-nodes.sbatch" "$root"
done
