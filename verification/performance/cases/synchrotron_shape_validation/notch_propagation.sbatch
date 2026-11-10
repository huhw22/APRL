#!/bin/bash
#SBATCH --job-name=synch-notch-z
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=8
#SBATCH --exclusive
#SBATCH --time=01:00:00

set -euo pipefail
source /public/home/wfollowm/APRL_Workspace/cluster.env
exec /public/home/wfollowm/micromamba/envs/physics/bin/python \
  "$1/scripts/run_notch_propagation.py" "$1"
