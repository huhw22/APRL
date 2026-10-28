#!/usr/bin/env bash
#SBATCH --job-name=field-drift-spectra
#SBATCH --partition=ssmb
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --exclusive
#SBATCH --mem=128G
#SBATCH --time=00:35:00
#SBATCH --no-requeue
set -Eeuo pipefail
result_dir=$(realpath -- "${1:?Fresh prepared result directory}")
workspace=$(realpath -- "${2:?APRL workspace}")
[[ -n ${SLURM_JOB_ID:-} ]] || exit 2
python_bin=${PYTHON_BIN:-/public/home/wfollowm/micromamba/envs/physics/bin/python3}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
cd "$result_dir"
trap 's=$?; printf "%s\n" "$s" > EXIT_CODE' EXIT
sha256sum --check scripts/SHA256SUMS
scontrol show job "$SLURM_JOB_ID" > allocation.txt
step=(srun --mpi=none --exact --nodes=1 --ntasks=1 --cpus-per-task=8 --cpu-bind=verbose,cores --kill-on-bad-exit=1)
"${step[@]}" "$python_bin" - <<'PY' > binding.json
import json,os,platform
a=sorted(os.sched_getaffinity(0)); assert len(a)==8,a
print(json.dumps(dict(host=platform.node(),job_id=os.environ['SLURM_JOB_ID'],allocated_cpus=128,
    actual_cpus=8,fft_workers=8,blas_threads=1,affinity=a,purpose='Read existing fields only; no simulation')))
PY
"${step[@]}" timeout 120 "$python_bin" scripts/analyze_spectra.py --self-test > self-test.log 2>&1
"${step[@]}" timeout 900 "$python_bin" scripts/reduce.py --kernel scripts/analyze_spectra.py \
    --run "$workspace/results/elegant_field_only_20261027/production" \
    --output "$result_dir/undulator" --plane 4.5 --job 18932 --slab-cells 178 --workers 8 > undulator.log 2>&1
"${step[@]}" timeout 900 "$python_bin" scripts/reduce.py --kernel scripts/analyze_spectra.py \
    --run "$workspace/results/elegant_free_drift_20261027/production" \
    --output "$result_dir/free_drift" --plane 0.5 --job 18938 --slab-cells 175 --workers 8 > free_drift.log 2>&1
"$python_bin" - <<'PY'
import json,os
from datetime import datetime,timezone
json.dump(dict(complete=True,job_id=os.environ['SLURM_JOB_ID'],completed_utc=datetime.now(timezone.utc).isoformat(),
    scope='Both closed raw fields reduced; plot and review follow'),open('REDUCTION_COMPLETE.json','w'),indent=2)
PY
