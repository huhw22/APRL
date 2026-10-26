#!/usr/bin/env bash
#SBATCH --job-name=aprl-detector-sampling
#SBATCH --partition=ssmb
#SBATCH --nodes=2
#SBATCH --ntasks=256
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=00:10:00
#SBATCH --no-requeue
set -Eeuo pipefail
[[ $# == 2 && ${SLURM_JOB_NUM_NODES:-0} == 2 ]] || exit 2
root=$(realpath "$1")
build=$(realpath "$2")
[[ ! -e "$root/STARTED" ]] || exit 2
source /public/home/wfollowm/APRL_Workspace/cluster.env
# Same accepted two-node runtime environment; no network settings experiment.
export OMPI_MCA_btl=self,vader,tcp
export OMPI_MCA_btl_openib_warn_no_device_params_found=0
export OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE OMP_PLACES=cores OMP_PROC_BIND=close
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 APRL_PROFILE_INTERVAL=0
export LD_LIBRARY_PATH="$build/runtime-libs${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
finish() { code=$?; trap - EXIT; printf '%s\n' "$code" > "$root/EXIT_CODE"; date --iso-8601=seconds > "$root/ENDED"; exit "$code"; }
trap finish EXIT
date --iso-8601=seconds > "$root/STARTED"
scontrol show job "$SLURM_JOB_ID" > "$root/allocation.txt"
launch=(srun --mpi=pmix_v3 --exact --nodes=2 --ntasks=64 --ntasks-per-node=32
  --cpus-per-task=1 --cpu-bind=cores --kill-on-bad-exit=1 --wait=10)
timeout --signal=TERM --kill-after=10s 60s "${launch[@]}" \
  "$APRL_WORKSPACE/build/cluster-tools/binding_probe" 64 1 > "$root/binding.log" 2>&1
grep -q '\[binding\] result=PASS failures=0' "$root/binding.log"
timeout --signal=TERM --kill-after=15s 180s "${launch[@]}" \
  "$build/verification/aprl_detector_sampling_checks" "$root/fixed" > "$root/fixed.log" 2>&1
# The old library must fail this new regression; preserve its actual raw output.
set +e
timeout --signal=TERM --kill-after=15s 180s "${launch[@]}" \
  "$build/aprl_detector_sampling_legacy_checks" "$root/legacy" --sweep-only > "$root/legacy.log" 2>&1
legacy_exit=$?
set -e
printf '%s\n' "$legacy_exit" > "$root/legacy-exit-code"
[[ $legacy_exit == 1 ]] || { echo 'Old implementation did not fail as expected'; exit 1; }
"$PYTHON_BIN" "$APRL_SOURCE_DIR/verification/detector_sampling/report.py" \
  "$root/fixed" --control "$root/legacy" > "$root/spectral-review.log"
