#!/bin/bash
#SBATCH --job-name=aprl-directory-check
#SBATCH --partition=ssmb
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=128
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=00:03:00
#SBATCH --no-requeue
set -euo pipefail
cd /public/home/wfollowm/APRL_Workspace
source ./cluster.env
export OMP_NUM_THREADS=1 OMP_PLACES=cores OMP_PROC_BIND=close
root=$(realpath "${1:?Pass a fresh prepared existing parent directory}")
[[ -f "$root/not-a-directory" ]]
exe=$(realpath "${2:?Pass aprl_output_directory_mpi_checks}")
scontrol show job "$SLURM_JOB_ID" > "$root/allocation.txt"
srun --mpi=pmix_v3 --exact --nodes=2 --ntasks=2 --ntasks-per-node=1 --cpus-per-task=1 --cpu-bind=cores --kill-on-bad-exit=1 --wait=15 --chdir="$root" "$PWD/build/cluster-tools/binding_probe" 2 1 > "$root/binding.log" 2>&1
for mode in new file parent-file; do
  if [[ "$mode" == new ]]; then args=(new "$root" fresh); else args=("$mode" "$root/not-a-directory"); fi
  timeout 45 srun --mpi=pmix_v3 --exact --nodes=2 --ntasks=2 --ntasks-per-node=1 --cpus-per-task=1 --cpu-bind=cores --kill-on-bad-exit=1 --wait=15 --chdir="$root" "$exe" "${args[@]}" > "$root/$mode.log" 2>&1
done
printf 'PASS\n' > "$root/PASS.txt"
