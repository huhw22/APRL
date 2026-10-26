#!/usr/bin/env bash
#SBATCH --job-name=retirement-field
#SBATCH --partition=ssmb
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --exclusive
#SBATCH --mem=16G
#SBATCH --time=00:20:00
set -Eeuo pipefail

run_dir=$(realpath -- "${1:?Production run directory}")
result_dir=$(realpath -- "${2:?Fresh result directory}")
field_file=$(realpath -- "${3:?Field-plane HDF5 file}")
[[ -n ${SLURM_JOB_ID:-} ]] || { echo 'Submit through Slurm' >&2; exit 2; }
python_bin=${PYTHON_BIN:-/public/home/wfollowm/micromamba/envs/physics/bin/python3}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export MPLBACKEND=Agg MPLCONFIGDIR="$result_dir/matplotlib-cache"
mkdir -p "$MPLCONFIGDIR"
cd "$result_dir"
sha256sum --check scripts/SHA256SUMS
scontrol show job "$SLURM_JOB_ID" > slurm-allocation.txt
step=(srun --mpi=none --exact --nodes=1 --ntasks=1 --cpus-per-task=1 --cpu-bind=verbose,cores --kill-on-bad-exit=1)
"${step[@]}" "$python_bin" - <<'PY' > binding.json
import json, os, platform
affinity = sorted(os.sched_getaffinity(0))
assert len(affinity) == 1, affinity
print(json.dumps(dict(host=platform.node(), job_id=os.environ['SLURM_JOB_ID'],
    allocated_cpus=128, actual_worker_cpus=1, affinity=affinity,
    purpose='Read existing detector fields; no simulation')))
PY
"${step[@]}" timeout 120 "$python_bin" scripts/analyze_fields.py --self-test > self-test.log 2>&1
"${step[@]}" "$python_bin" scripts/analyze_fields.py --run-dir "$run_dir" --field-file "$field_file" \
    --field-only --no-plots --block-frames 16 --output-dir "$result_dir/derived" > analysis.log 2>&1
"$python_bin" - <<'PY'
import json, os
from datetime import datetime, timezone
json.dump(dict(complete=True, job_id=os.environ['SLURM_JOB_ID'],
    completed_utc=datetime.now(timezone.utc).isoformat(),
    scope='Committed detector fields checked finite and reduced; plotting is a separate stage'),
    open('REDUCTION_COMPLETE.json', 'w'), indent=2)
PY
