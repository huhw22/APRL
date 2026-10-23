#!/usr/bin/env bash
#SBATCH --job-name=elegant-analysis
#SBATCH --partition=ssmb
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --exclusive
#SBATCH --mem=32G
#SBATCH --time=02:00:00
set -Eeuo pipefail

result_root=$(realpath -- "${1:?Provide the verified analysis result directory}")
script_root="$result_root/scripts"
[[ -n ${SLURM_JOB_ID:-} ]] || { echo 'Must run through Slurm' >&2; exit 2; }
python_bin=${PYTHON_BIN:-/public/home/wfollowm/micromamba/envs/physics/bin/python3}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export MPLBACKEND=Agg
mkdir -p "$result_root/derived"
analysis_dir=$(mktemp -d "$result_root/derived/job-${SLURM_JOB_ID}-XXXXXXXX")
export MPLCONFIGDIR="$analysis_dir/matplotlib-cache"
mkdir "$MPLCONFIGDIR"
cd "$script_root"
sha256sum --check SHA256SUMS
"$python_bin" - "$result_root/copy-manifest.json" <<'PY'
import json, sys
m = json.load(open(sys.argv[1]))
assert m['complete'] is True
assert all(f.get('verified', False) for f in m['files'] if 'size_bytes' in f)
print('Verified historical snapshot:', m['total_regular_bytes'], 'bytes')
PY

run_dir="$result_root/original/runs/20261021T070618Z-861230"
scontrol show job "$SLURM_JOB_ID" > "$analysis_dir/slurm-allocation.txt"
"$python_bin" - "$analysis_dir/run-metadata.json" <<'PY'
import json, os, platform, sys
from datetime import datetime, timezone
import h5py, numpy, matplotlib
json.dump(dict(started_utc=datetime.now(timezone.utc).isoformat(),
    host=platform.node(), python=sys.version, numpy=numpy.__version__,
    h5py=h5py.__version__, matplotlib=matplotlib.__version__,
    slurm_job_id=os.environ['SLURM_JOB_ID'], reserved_cpus=128,
    actual_sequential_worker_cpus=1, library_threads=1,
    purpose='Historical diagnostic postprocessing; no simulation'), open(sys.argv[1], 'w'), indent=2)
PY

# One CPU per sequential step. The allocation remains an exclusive 128-CPU node.
step=(srun --mpi=none --exact --nodes=1 --ntasks=1 --cpus-per-task=1 --cpu-bind=verbose,cores --kill-on-bad-exit=1)
"${step[@]}" "$python_bin" -c 'import os,platform; a=os.sched_getaffinity(0); print(platform.node(), "allowed CPUs:",sorted(a)); assert len(a)==1' > "$analysis_dir/binding.log" 2>&1
"${step[@]}" timeout 120 "$python_bin" "$script_root/analyze_particles.py" --self-test > "$analysis_dir/particle-self-test.log" 2>&1
"${step[@]}" timeout 120 "$python_bin" "$script_root/analyze_fields.py" --self-test > "$analysis_dir/field-self-test.log" 2>&1
"${step[@]}" "$python_bin" "$script_root/analyze_particles.py" --run-dir "$run_dir" --output-dir "$analysis_dir/particles" > "$analysis_dir/particles.log" 2>&1
"${step[@]}" "$python_bin" "$script_root/analyze_fields.py" --run-dir "$run_dir" --output-dir "$analysis_dir/fields" > "$analysis_dir/fields.log" 2>&1
"$python_bin" - "$analysis_dir" <<'PY'
import json, sys
from datetime import datetime, timezone
from pathlib import Path
p = Path(sys.argv[1])
json.dump(dict(complete=True, completed_utc=datetime.now(timezone.utc).isoformat(),
    scope='Five fixed lab particle planes, raw lab field plane and boost ledger',
    interpretation_review_required=True), (p/'COMPLETE.json').open('w'), indent=2)
print('Analysis products:', p)
PY
