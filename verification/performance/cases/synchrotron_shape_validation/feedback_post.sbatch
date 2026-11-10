#!/bin/bash
#SBATCH --job-name=synch-nofb-post
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=00:12:00

set -euo pipefail
export MPLCONFIGDIR="$1/mpl-cache-post"
exec /public/home/wfollowm/micromamba/envs/physics/bin/python \
  "$1/scripts/analyze_feedback_control.py" "$1"
