#!/bin/bash
#SBATCH --job-name=physics-labels
#SBATCH --partition=ssmb
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:20:00
#SBATCH --no-requeue
set -Eeuo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export MPLCONFIGDIR=/tmp/physics-labels-${SLURM_JOB_ID}
/public/home/wfollowm/micromamba/envs/physics/bin/python3 -B "$1/polish.py" "$2" /public/home/wfollowm/APRL_Workspace/results/physics_chinese_print_20261112_r1/DroidSansFallback.ttf /public/home/wfollowm/APRL_Workspace/results/physics_visual_validation_20261112_r1/vacuum/sphere_gaussian_r3/METADATA.json > /public/home/wfollowm/APRL_Workspace/results/display-polish-20261117.log 2>&1
