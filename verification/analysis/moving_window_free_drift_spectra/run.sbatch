#!/usr/bin/env bash
#SBATCH --job-name=moving-drift-spectrum
#SBATCH --partition=ssmb
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --exclusive
#SBATCH --mem=128G
#SBATCH --time=00:20:00
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
    actual_cpus=8,fft_workers=8,blas_threads=1,affinity=a,purpose='Read completed 18993 fields only; no simulation')))
PY
"${step[@]}" timeout 120 "$python_bin" scripts/analyze_spectra.py --self-test > self-test.log 2>&1
"${step[@]}" timeout 1000 "$python_bin" scripts/reduce.py --kernel scripts/analyze_spectra.py \
    --run "$workspace/results/moving_window_comparison_20261029/free_drift_0p5m" \
    --output "$result_dir/new_free_drift" --plane 0.5 --job 18993 --workers 8 > reduction.log 2>&1
"$python_bin" - <<'PY'
import json,os
from datetime import datetime,timezone
json.dump(dict(complete=True,job_id=os.environ['SLURM_JOB_ID'],completed_utc=datetime.now(timezone.utc).isoformat(),
    scope='Closed 18993 field reduced; old 18938 reduced data reused; plot and review follow'),open('REDUCTION_COMPLETE.json','w'),indent=2)
PY
