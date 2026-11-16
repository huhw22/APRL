#!/bin/bash
#SBATCH --job-name=physics-display
#SBATCH --partition=ssmb
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:40:00
#SBATCH --no-requeue
set -Eeuo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export MPLCONFIGDIR=/tmp/physics-display-${SLURM_JOB_ID}
template="$1"
root="$2"
font=/public/home/wfollowm/APRL_Workspace/results/physics_chinese_print_20261112_r1/DroidSansFallback.ttf
python=/public/home/wfollowm/micromamba/envs/physics/bin/python3
"$python" -B "$template/render.py" "$root" "$font" > "$root/../display-render-20261117.log" 2>&1
"$python" -B "$template/movies.py" "$root" "$font" /public/home/wfollowm/APRL_Workspace/results/physics_visual_validation_20261112_r1/vacuum/sphere_gaussian_r3/METADATA.json > "$root/../display-movies-20261117.log" 2>&1
