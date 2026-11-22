#!/bin/bash
#SBATCH --job-name=thz-refined
#SBATCH --partition=ssmb
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=8
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --array=1-3%3
#SBATCH --time=24:00:00
#SBATCH --no-requeue
#SBATCH --output=/public/home/wfollowm/APRL_Workspace/results/thz_undulator_aprl21_stable_refined_20261122_r1/slurm-%A_%a.log
set -eo pipefail
source /public/home/wfollowm/APRL_Workspace/cluster.env
set -u
export OMP_NUM_THREADS=8 OMP_DYNAMIC=FALSE OMP_PLACES=cores OMP_PROC_BIND=close
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_MAX_THREADS=1
export OMPI_MCA_pml=ob1 OMPI_MCA_btl=self,vader,tcp
root="$APRL_WORKSPACE/results/thz_undulator_aprl21_stable_refined_20261122_r1"
export APRL_THZ_CASE="$root/f0-div-$SLURM_ARRAY_TASK_ID"
export MPLCONFIGDIR="$APRL_THZ_CASE/mpl-cache"
binary="$APRL_WORKSPACE/results/source_coupling_extended_20261122_r1/runtime/bin/aprl"
python=/public/home/wfollowm/micromamba/envs/physics/bin/python3
cd "$APRL_THZ_CASE"
"$python" -B "$root/scripts/verify_inputs.py"
test ! -f JOB_STARTED.json
"$python" -B -c 'import json,os,socket;from pathlib import Path;Path("JOB_STARTED.json").write_text(json.dumps({"job":os.environ["SLURM_JOB_ID"],"array":os.environ["SLURM_ARRAY_JOB_ID"],"case":os.environ["SLURM_ARRAY_TASK_ID"],"host":socket.gethostname(),"stage":"native-initial-reference"},indent=2)+"\n")'
trap 'echo "$?" > "$APRL_THZ_CASE/PIPELINE_EXIT_CODE"' EXIT
test "$(sha256sum "$binary" | cut -d ' ' -f1)" = 87ecf9e8ec701b3bb7ffee0ef9a562f1a267f06a7f8cbe5c3afb28f0eb196af0
test "$(sha256sum "$root/beam.h5" | cut -d ' ' -f1)" = 7263f0bdaca3ee9f7773c8e6e42a38f38034a287ab4e0cf2859c69a564cd731e
echo native-initial-reference > STAGE
cd preflight
set +e
srun --mpi="$SRUN_MPI" --exact --nodes=1 --ntasks=16 --cpus-per-task=8 --cpu-bind=cores --kill-on-bad-exit=1 --wait=60 "$binary" run.yaml > aprl.log 2>&1
code=$?
echo "$code" > EXIT_CODE
set -e
# This stage prepares the complete initial reference, then intentionally ends
# after nine steps. It is not a completed undulator simulation. Preserve its
# actual nonzero duration-limit exit, and reject every other failure.
test "$code" != 0
grep -q 'mesh.duration was exhausted before the configured physical stop; trajectory and detector outputs are readable but marked incomplete' aprl.log
test -f reference/metadata.txt
grep -q 'status complete' reference/metadata.txt
grep -q 'Moving window completed: steps=9 particle_exits=0 cropped_carriers=0' aprl.log
"$python" -B -c 'import json;from pathlib import Path;Path("REFERENCE_READY.json").write_text(json.dumps({"reference_complete":True,"physical_run_complete":False,"expected_duration_guard":True,"actual_exit_code":int(Path("EXIT_CODE").read_text()),"steps":9},indent=2)+"\n")'
cd "$APRL_THZ_CASE"
echo startup-fft > STAGE
srun --mpi=none --exact --ntasks=1 --cpus-per-task=8 --cpu-bind=cores "$python" -Bu "$root/scripts/startup.py" prepare > startup-prepare.log 2>&1
echo startup-modes > STAGE
OMP_NUM_THREADS=1 srun --mpi=none --exact --nodes=1 --ntasks=48 --cpus-per-task=1 --cpu-bind=cores --kill-on-bad-exit=1 "$python" -Bu "$root/scripts/startup.py" worker > startup-modes.log 2>&1
echo startup-export > STAGE
srun --mpi=none --exact --ntasks=1 --cpus-per-task=8 --cpu-bind=cores "$python" -Bu "$root/scripts/startup.py" finish 48 > startup-finish.log 2>&1
test -f startup/arrays/READY.json
echo self-consistent-stability-gate > STAGE
cd stability-gate
set +e
srun --mpi="$SRUN_MPI" --exact --nodes=1 --ntasks=16 --cpus-per-task=8 --cpu-bind=cores --kill-on-bad-exit=1 --wait=60 "$binary" run.yaml > aprl.log 2>&1
code=$?
echo "$code" > EXIT_CODE
set -e
test "$code" != 0
grep -q 'mesh.duration was exhausted before the configured physical stop; trajectory and detector outputs are readable but marked incomplete' aprl.log
grep -q 'Moving window completed: steps=256 particle_exits=0 cropped_carriers=0' aprl.log
"$python" -B "$root/scripts/check_stability.py" "$APRL_THZ_CASE"
cd "$APRL_THZ_CASE"
echo self-consistent-undulator > STAGE
set +e
srun --mpi="$SRUN_MPI" --exact --nodes=1 --ntasks=16 --cpus-per-task=8 --cpu-bind=cores --kill-on-bad-exit=1 --wait=60 "$binary" run.yaml > aprl.log 2>&1
code=$?
echo "$code" > EXIT_CODE
set -e
test "$code" = 0
grep -q 'Configured stop reached:' aprl.log
echo offline-spectrum > STAGE
srun --mpi=none --exact --ntasks=1 --cpus-per-task=8 --cpu-bind=cores "$python" -B "$root/scripts/spectrum.py" "$APRL_THZ_CASE" > spectrum.log 2>&1
echo validation-and-plots > STAGE
srun --mpi=none --exact --ntasks=1 --cpus-per-task=8 --cpu-bind=cores "$python" -B "$root/scripts/analyze.py" "$APRL_THZ_CASE" > analysis.log 2>&1
echo complete > STAGE
