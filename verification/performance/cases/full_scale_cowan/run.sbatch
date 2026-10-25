#!/usr/bin/env bash
#SBATCH --job-name=aprl-cowan-full-pair
#SBATCH --partition=ssmb
#SBATCH --nodes=2
#SBATCH --ntasks=256
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=03:00:00
#SBATCH --no-requeue
set -Eeuo pipefail
root=$(realpath "${1:?Pass a new frozen full-scale comparison directory}")
source /public/home/wfollowm/APRL_Workspace/cluster.env
[[ ${SLURM_JOB_NUM_NODES:-0} == 2 && ${SLURM_NTASKS:-0} == 256 && ${SLURM_CPUS_PER_TASK:-0} == 1 ]] || exit 2
export OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE OMP_PLACES=cores OMP_PROC_BIND=close
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export APRL_PROFILE_INTERVAL=0 HDF5_USE_FILE_LOCKING=FALSE
# Keep the existing production transport setting; no system/network changes.
export OMPI_MCA_btl=self,vader,tcp
export OMPI_MCA_btl_openib_warn_no_device_params_found=0
ulimit -s unlimited
cd "$root"
trap 'rc=$?; echo "$rc" > "$root/JOB_EXIT_CODE"; date --iso-8601=seconds > "$root/JOB_ENDED_AT"' EXIT
date --iso-8601=seconds > JOB_STARTED_AT
sha256sum --check frozen.sha256 > frozen-verification.log
scontrol show job "$SLURM_JOB_ID" > allocation.txt
env | sort | sed -n '/^SLURM_/p;/^OMP_/p;/^APRL_PROFILE_INTERVAL=/p;/^OMPI_MCA_btl=/p' > execution-environment.txt
step=(srun --mpi="$SRUN_MPI" --exact --nodes=2 --ntasks=256 --ntasks-per-node=128 --cpus-per-task=1 --cpu-bind=cores --kill-on-bad-exit=1 --wait=60)
echo "Binding check $(date --iso-8601=seconds)"
timeout --signal=TERM --kill-after=30 120 "${step[@]}" \
  "$APRL_WORKSPACE/build/cluster-tools/binding_probe" 256 1 > binding.log 2>&1
timeout --signal=TERM --kill-after=30 300 "$PYTHON_BIN" scripts/run.py verify --root "$root" > verification.log 2>&1

for variant in control trial; do
  binary=$("$PYTHON_BIN" -c 'import json,sys; print(json.load(open(sys.argv[1]))["variants"][sys.argv[2]]["binary"])' plan.json "$variant")
  case_dir="$root/preflight/$variant"
  echo "Preflight $variant $(date --iso-8601=seconds)"
  set +e
  timeout --signal=TERM --kill-after=30 600 "${step[@]}" "$binary" "$case_dir/run.yaml" > "$case_dir/aprl.log" 2>&1
  rc=$?
  set -e
  echo "$rc" > "$case_dir/EXIT_CODE"
  [[ $rc == 0 ]] || exit "$rc"
done
timeout --signal=TERM --kill-after=30 300 "$PYTHON_BIN" scripts/run.py check-preflight --root "$root" > preflight-check.log 2>&1

for variant in control trial; do
  binary=$("$PYTHON_BIN" -c 'import json,sys; print(json.load(open(sys.argv[1]))["variants"][sys.argv[2]]["binary"])' plan.json "$variant")
  case_dir="$root/$variant"
  echo "Full-size $variant start $(date --iso-8601=seconds)"
  date --iso-8601=seconds > "$case_dir/STARTED_AT"
  set +e
  timeout --signal=TERM --kill-after=60 4800 "${step[@]}" "$binary" "$case_dir/run.yaml" > "$case_dir/aprl.log" 2>&1
  rc=$?
  set -e
  echo "$rc" > "$case_dir/EXIT_CODE"
  date --iso-8601=seconds > "$case_dir/ENDED_AT"
  # The 600-step duration stop is expected to return nonzero. Only the
  # checker may accept it; a wall-clock timeout must never count as success.
  [[ $rc != 124 && $rc != 137 ]] || exit "$rc"
  timeout --signal=TERM --kill-after=30 300 "$PYTHON_BIN" scripts/run.py check-run --root "$root" --variant "$variant" > "$case_dir/check.log" 2>&1
done
timeout --signal=TERM --kill-after=30 600 "$PYTHON_BIN" scripts/run.py compare --root "$root" > comparison.log 2>&1
echo "Bounded pair finished $(date --iso-8601=seconds)"
