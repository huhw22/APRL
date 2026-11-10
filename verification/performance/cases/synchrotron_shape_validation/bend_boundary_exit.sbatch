#!/bin/bash
#SBATCH --job-name=synch-boundary
#SBATCH --partition=ssmb
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=8
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=01:00:00
#SBATCH --no-requeue

set -euo pipefail
source /public/home/wfollowm/APRL_Workspace/cluster.env
exec /public/home/wfollowm/micromamba/envs/physics/bin/python \
  "$1/scripts/run_100k_bend_boundary_exit.py" "$1"
