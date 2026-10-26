#!/usr/bin/env bash
#SBATCH --job-name=raw-field-spectra
#SBATCH --partition=ssmb
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --exclusive
#SBATCH --mem=128G
#SBATCH --time=00:45:00
set -Eeuo pipefail
result_dir=$(realpath -- "${1:?Fresh result directory}")
workspace=$(realpath -- "${2:?APRL workspace}")
[[ -n ${SLURM_JOB_ID:-} ]] || { echo 'Submit through Slurm' >&2; exit 2; }
python_bin=${PYTHON_BIN:-/public/home/wfollowm/micromamba/envs/physics/bin/python3}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
cd "$result_dir"
sha256sum --check scripts/SPECTRA_SHA256SUMS
scontrol show job "$SLURM_JOB_ID" > slurm-allocation.txt
step=(srun --mpi=none --exact --nodes=1 --ntasks=1 --cpus-per-task=8 --cpu-bind=verbose,cores --kill-on-bad-exit=1)
"${step[@]}" "$python_bin" - <<'PY' > binding.json
import json,os,platform
a=sorted(os.sched_getaffinity(0)); assert len(a)==8,a
print(json.dumps(dict(host=platform.node(),job_id=os.environ['SLURM_JOB_ID'],allocated_cpus=128,
    actual_cpus=8,scipy_fft_workers=8,blas_threads=1,affinity=a,purpose='Existing raw field spectra; no simulation')))
PY
"${step[@]}" timeout 120 "$python_bin" scripts/analyze_spectra.py --self-test > self-test.log 2>&1
new_run="$workspace/results/elegant_retirement_20261024/compact-observed-window/production"
old_run="$workspace/results/elegant_4m_analysis_20261024/original/runs/20261021T070618Z-861230"
"${step[@]}" "$python_bin" scripts/analyze_spectra.py --run-dir "$new_run" \
    --field-file "$new_run/detectors/post-retirement-field.h5" --output-dir "$result_dir/retirement" \
    --slab-cells 182 --workers 8 > retirement.log 2>&1
"${step[@]}" "$python_bin" scripts/analyze_spectra.py --run-dir "$old_run" \
    --field-file "$old_run/detectors/exit-field.h5" --output-dir "$result_dir/original" \
    --slab-cells 178 --workers 8 > original.log 2>&1
"$python_bin" - <<'PY'
import json,os
from datetime import datetime,timezone
json.dump(dict(complete=True,job_id=os.environ['SLURM_JOB_ID'],
    completed_utc=datetime.now(timezone.utc).isoformat(),scope='Both committed raw fields reduced; plots separate'),
    open('REDUCTION_COMPLETE.json','w'),indent=2)
PY
