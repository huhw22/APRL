#!/usr/bin/env bash
#SBATCH --job-name=und19099-review
#SBATCH --partition=ssmb
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --exclusive
#SBATCH --mem=128G
#SBATCH --time=00:25:00
#SBATCH --no-requeue
set -Eeuo pipefail
root=$(realpath "${1:?Fresh review directory}")
source /public/home/wfollowm/APRL_Workspace/cluster.env
cd "$root"
finish() { code=$?; trap - EXIT; printf '%s\n' "$code" > EXIT_CODE; date --iso-8601=seconds > ENDED_AT; exit "$code"; }
trap finish EXIT
date --iso-8601=seconds > STARTED_AT
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export MPLCONFIGDIR="$root/mpl-cache"
"$PYTHON_BIN" - <<'PY'
import hashlib,json
from pathlib import Path
for name,digest in json.load(open('FROZEN.json')).items():
 assert hashlib.sha256(Path(name).read_bytes()).hexdigest()==digest,name
PY
timeout --signal=TERM --kill-after=30 1350 srun --mpi=none --exact --nodes=1 --ntasks=1 --cpus-per-task=16 --cpu-bind=cores --kill-on-bad-exit=1 \
 "$PYTHON_BIN" scripts/analyze_undulator.py --run "$APRL_WORKSPACE/results/corrected_full_comparison_20261030/r3_2nodes/undulator_4p5m" --out "$root/science" > analysis.log 2>&1
