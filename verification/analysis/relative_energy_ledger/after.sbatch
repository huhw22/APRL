#!/usr/bin/env bash
#SBATCH --job-name=aprl-relative-ledger
#SBATCH --partition=ssmb
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=00:05:00
#SBATCH --no-requeue
set -Eeuo pipefail
root=$(realpath "${1:?Relative-report root}")
source /public/home/wfollowm/APRL_Workspace/cluster.env
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 OMP_NUM_THREADS=1
"$PYTHON_BIN" - "$root" <<'PY'
import sys,json,hashlib
from pathlib import Path
root=Path(sys.argv[1])
for path,digest in json.loads((root/'SCRIPTS_SHA256.json').read_text()).items():
    assert hashlib.sha256(Path(path).read_bytes()).hexdigest()==digest,path
PY
timeout --signal=TERM --kill-after=10 240 "$PYTHON_BIN" "$root/scripts/render.py" --plan "$root/plan.json" --source-job-ended > "$root/after-job.log" 2>&1
