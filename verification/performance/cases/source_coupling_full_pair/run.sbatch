#!/usr/bin/env bash
#SBATCH --job-name=aprl-source-full
#SBATCH --partition=ssmb
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=8
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=144:00:00
#SBATCH --no-requeue
set -Eeuo pipefail
root=$(realpath "${1:?Frozen full case directory}")
source /public/home/wfollowm/APRL_Workspace/cluster.env
[[ ! -e "$root/JOB_STARTED_AT" && ! -e "$root/run-manifest.yaml" ]]
cd "$root"
finish() { code=$?; trap - EXIT; printf '%s\n' "$code" > EXIT_CODE; date --iso-8601=seconds > JOB_ENDED_AT; exit "$code"; }
trap finish EXIT
trap 'exit 143' TERM
trap 'exit 130' INT
date --iso-8601=seconds > JOB_STARTED_AT
scontrol show job "$SLURM_JOB_ID" > allocation.txt
export OMP_NUM_THREADS=8 OMP_DYNAMIC=FALSE OMP_PLACES=cores OMP_PROC_BIND=close
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 APRL_PROFILE_INTERVAL=0
export HDF5_USE_FILE_LOCKING=FALSE OMPI_MCA_pml=ob1 OMPI_MCA_btl=self,vader,tcp
unset APRL_PARTICLE_MOMENTUM_PROBE APRL_PARTICLE_MOMENTUM_PROBE_FIRST APRL_PARTICLE_MOMENTUM_PROBE_STRIDE APRL_PARTICLE_MOMENTUM_PROBE_DT
unset APRL_INLET_CONSTRAINT_PROBE APRL_GRID_WORK_AUDIT APRL_FIELD_STAGE_AUDIT APRL_PARTICLE_WORK_GEOMETRY APRL_SOURCE_EXIT_AUDIT
unset APRL_DIAGNOSTIC_STEPS APRL_FORCE_SPLIT_DIRECTORY
ulimit -s unlimited
ulimit -c 0
read -r nodes ranks binary < <("$PYTHON_BIN" -c 'import json;p=json.load(open("run-plan.json"));print(p["nodes"],p["ranks"],p["binary"])')
[[ ${SLURM_JOB_NUM_NODES:-0} == "$nodes" && ${SLURM_NTASKS:-0} == "$ranks" && ${SLURM_CPUS_PER_TASK:-0} == 8 ]]
"$PYTHON_BIN" "$root/../scripts/startup_identity.py" "$root" > identity-check.json
step=(srun --mpi="$SRUN_MPI" --exact --nodes="$nodes" --ntasks="$ranks" --ntasks-per-node=16 --cpus-per-task=8 --cpu-bind=cores --kill-on-bad-exit=1)
timeout --signal=TERM --kill-after=30s 120s "${step[@]}" --wait=30 --chdir="$root" \
  "$APRL_WORKSPACE/build/cluster-tools/binding_probe" "$ranks" 8 > binding.log 2>&1
"$PYTHON_BIN" - "$root" <<'PY' > binding-check.json
import json,sys
from pathlib import Path
r=Path(sys.argv[1]);sys.path.insert(0,str(r.parent/'scripts'))
from completion_helpers import check_binding
print(json.dumps(check_binding(r/'binding.log',json.loads((r/'run-plan.json').read_text())),indent=2))
PY
date --iso-8601=seconds > PRODUCTION_STARTED_AT
if "${step[@]}" --wait=660 --chdir="$root" "$binary" "$root/run.yaml" > aprl.log 2>&1; then
  printf '0\n' > SIMULATION_EXIT_CODE
else
  code=$?; printf '%s\n' "$code" > SIMULATION_EXIT_CODE; exit "$code"
fi
"$PYTHON_BIN" "$root/../scripts/check_completion.py" "$root" > completion.log 2>&1
"$PYTHON_BIN" "$root/../scripts/check_source_coupling.py" "$root" > source-completion.log 2>&1
