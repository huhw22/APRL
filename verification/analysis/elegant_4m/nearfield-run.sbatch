#!/usr/bin/env bash
#SBATCH --job-name=elegant-nearfield
#SBATCH --partition=ssmb
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --exclusive
#SBATCH --mem=96G
#SBATCH --time=04:00:00
set -Eeuo pipefail
root=$(realpath "${1:?Result directory}")
source /public/home/wfollowm/APRL_Workspace/cluster.env
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export MPLCONFIGDIR="$root/original-card/mpl-cache" MPLBACKEND=Agg
mkdir -p "$MPLCONFIGDIR"
cd "$root"
sha256sum --check full-run.sha256
"$PYTHON_BIN" - "$root" <<'PY'
import json,sys,yaml
from pathlib import Path
p=Path(sys.argv[1]);r=json.loads((p/'projection-preflight.json').read_text())
assert r['unique_complete_ID_cohort']
assert r['relative_gamma_spread']<=.05 and r['rms_transverse_beta']<=.05
assert r['grids'][0]['padding']==[1,1,1] and r['grids'][0]['outside_charge_fraction']<=.01
c=yaml.safe_load((p/'reconstruction.yaml').read_text())
assert c['model']['padding_factor']==[1,1,1] and c['model']['fft_threads']==1
assert c['model']['maximum_outside_charge_fraction']==.01
assert c['output']['overwrite'] is False
print('Original model and thresholds preserved; complete-cohort projection preflight PASS')
PY
scontrol show job "$SLURM_JOB_ID" > "$root/original-card/allocation.txt"
step=(srun --mpi=none --exact --nodes=1 --ntasks=1 --cpus-per-task=1 --cpu-bind=verbose,cores --kill-on-bad-exit=1)
"${step[@]}" "$PYTHON_BIN" -c 'import os; a=os.sched_getaffinity(0); print("Allowed CPUs:",sorted(a)); assert len(a)==1' > "$root/original-card/binding.log" 2>&1
"${step[@]}" "$PYTHON_BIN" "$root/scripts/plot_reconstructed_flux.py" --self-test > "$root/original-card/plot-self-test.log" 2>&1
"${step[@]}" /usr/bin/time -v /public/home/wfollowm/APRL_Workspace/build/verification/postprocess/field_reconstruction/field_reconstruction "$root/reconstruction.yaml" > "$root/original-card/reconstruction.log" 2>&1
"${step[@]}" "$PYTHON_BIN" "$root/scripts/plot_reconstructed_flux.py" \
  --reconstructed-file "$root/original-card/reconstructed.h5" \
  --output-dir "$root/original-card/flux" > "$root/original-card/plot-flux.log" 2>&1
"$PYTHON_BIN" - "$root/original-card/COMPLETE.json" <<'PY'
from datetime import datetime,timezone
import json,os,sys
json.dump(dict(complete=True,slurm_job_id=os.environ['SLURM_JOB_ID'],
              completed_utc=datetime.now(timezone.utc).isoformat(),padding=[1,1,1],fft_threads=1,
              original_model_preserved=True,convergence_claim=False),open(sys.argv[1],'w'),indent=2)
PY
