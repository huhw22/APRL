#!/bin/bash
#SBATCH --job-name=synch-nofb
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=8
#SBATCH --exclusive
#SBATCH --time=01:20:00

set -euo pipefail
source /public/home/wfollowm/APRL_Workspace/cluster.env
exec /public/home/wfollowm/micromamba/envs/physics/bin/python \
  "$1/scripts/run_feedback_control.py" "$1"
