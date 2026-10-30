#!/usr/bin/env bash
#SBATCH --job-name=aprl-ID-recheck
#SBATCH --partition=ssmb
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:10:00
#SBATCH --no-requeue
set -Eeuo pipefail
repair=$(realpath "${1:?Prepared repair directory}")
case_dir=$(realpath "${2:?Closed simulation directory}")
source /public/home/wfollowm/APRL_Workspace/cluster.env
cd "$repair"
finish() { code=$?; trap - EXIT; printf '%s\n' "$code" > EXIT_CODE; date --iso-8601=seconds > ENDED_AT; exit "$code"; }
trap finish EXIT
date --iso-8601=seconds > STARTED_AT
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export MPLCONFIGDIR="$repair/mpl-cache"
"$PYTHON_BIN" - <<'PY'
import hashlib,json
from pathlib import Path
for name,digest in json.load(open('FROZEN.json')).items():
 assert hashlib.sha256(Path(name).read_bytes()).hexdigest()==digest,name
PY
timeout --signal=TERM --kill-after=30 540 srun --mpi=none --exact --nodes=1 --ntasks=1 --cpus-per-task=4 --cpu-bind=cores --kill-on-bad-exit=1 \
 "$PYTHON_BIN" scripts/check_completion.py "$case_dir" --output "$case_dir/completion-check.corrected.json" > recheck.log 2>&1
