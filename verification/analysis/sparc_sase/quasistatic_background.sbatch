#!/usr/bin/env bash
#SBATCH --job-name=SPARC-background
#SBATCH --partition=ssmb
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=00:45:00
#SBATCH --no-requeue
set -Eeuo pipefail
analysis_root="$1"
source /public/home/wfollowm/APRL_Workspace/cluster.env
cd "$analysis_root"
finish() { result=$?; trap - EXIT; printf '%s\n' "$result" > EXIT_CODE; exit "$result"; }
trap finish EXIT
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=8 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export OMPI_MCA_btl=self,vader,tcp HDF5_USE_FILE_LOCKING=FALSE MPLCONFIGDIR="$analysis_root/mpl-cache"
scontrol show job "$SLURM_JOB_ID" > allocation.txt
"$PYTHON_BIN" - <<'PY'
from pathlib import Path
import hashlib,json
for manifest in ['INPUT_FROZEN.json','TOOLS_FROZEN.json']:
    for path,digest in json.loads(Path(manifest).read_text()).items():
        assert hashlib.sha256(Path(path).read_bytes()).hexdigest()==digest,path
PY
for control in base double; do
  timeout --signal=TERM --kill-after=20 120 srun --mpi=pmix_v3 --exact --nodes=1 --ntasks=1 \
    --cpus-per-task=8 --cpu-bind=cores --kill-on-bad-exit=1 --wait=30 \
    bin/quasistatic_probe controls.yaml "$control" > "control-$control.log" 2>&1
done
"$PYTHON_BIN" scripts/quasistatic_background.py controls "$analysis_root"
for method in boris energy; do
  timeout --signal=TERM --kill-after=20 1000 srun --mpi=pmix_v3 --exact --nodes=1 --ntasks=1 \
    --cpus-per-task=8 --cpu-bind=cores --kill-on-bad-exit=1 --wait=30 \
    bin/quasistatic_probe config.yaml "$method" > "$method-probe.log" 2>&1
  timeout --signal=TERM --kill-after=20 280 srun --mpi=none --exact --nodes=1 --ntasks=1 \
    --cpus-per-task=8 --cpu-bind=cores "$PYTHON_BIN" scripts/quasistatic_background.py reduce \
    /public/home/wfollowm/APRL_Workspace/results/sparc_sase_comparison_20261103 \
    "$analysis_root" "$method" > "$method-spectrum.log" 2>&1
done
"$PYTHON_BIN" scripts/quasistatic_background.py summarize "$analysis_root" > summary.log 2>&1
