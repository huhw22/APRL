Small edge acceptance, fixed 4 MPI x 4 OpenMP

This directory is a reviewable runner, not an executed campaign. The local
science_helpers.py and held.yaml.in are exact git-show copies from 856e172.
When versioning this template, omit science_helpers.py: the loader falls back to
../full_scale_cowan/run.py. Prepare freezes that helper in each result's scripts.
Only the scientific
inspection/comparison helpers are used; their historical campaign configuration
is not used. The output comparison checks every dataset and scientific attribute
bitwise, ignoring only the existing four run/source identity attributes.

Cases, independently identical across accepted control and every candidate:
  held-crossrank: existing 10410-particle hybrid held H5; ballistic mode; require
    actual cross-rank activation pushes and complete particle detector cohorts.
  held-retirement: existing 512-particle, 16x16x4096 non-ballistic fixture; require
    held startup, 512 retirement entries/completions and post-retirement fields.
  cpml-exit / outer-exit: eight generated particles near the +x surface; require
    the actual respective boundary-exit report. These new, tiny cards have not
    yet been executed. A failed coverage assertion stops the campaign.

72required already exercises held cross-rank release, but its dedicated hybrid
case is ballistic. It does not assert retirement lifecycle or simultaneous event
ordering. These extra cases test ordinary cached hits and held/received fallback.
They do NOT test exact or tolerance-near simultaneous CPML/outer/retirement/
MPI-interface events. The patch preserves both independent selectors unchanged;
do not claim that structural evidence is an executed tie-priority regression.

Prepare only after all candidate required checks pass. Create a JSON object with
2 to 5 entries: one named control and 1 to 4 named candidates. Every entry has
binary (absolute path), revision (full compiled Git hash),
and each non-control entry has required_directory (containing PASS, EXIT_CODE=0,
and build.sha256 for that exact executable). The control entry must identify:
  /public/home/wfollowm/APRL_Workspace/build/cowan-accepted-r2-20261025/aprl
  revision 77f1ef95267d037dcd629584e42aedc1b04668f3

Four cases produce 8 to 20 runs, depending on candidate count. For example,
control plus cowan, particle, diagnostic and combined produces 20 runs. All
candidates, including combined, require their own completed required evidence.

Example invocation after sourcing cluster.env (RESULT must be new):
  "$PYTHON_BIN" /public/home/wfollowm/APRL_Workspace/source/aprl/verification/performance/cases/optimization_edges/run.py prepare --root "$RESULT" --variants "$VARIANTS_JSON"

Run the frozen copy on authorized gnode1 with APRL_ALLOW_UNSCHEDULED_CPU_TESTS=1:
  "$PYTHON_BIN" "$RESULT/scripts/run.py" run --root "$RESULT"
If a resource-checking wrapper already owns descriptor 9 for the canonical CPU
calibration lock, append --lock-fd 9. The runner verifies the inherited descriptor
and reacquires the same open-file lock. Otherwise it acquires that lock itself.
After acquiring the lock, the runner also requires no existing APRL/MPI process
and at least 95% aggregate CPU idle during a bounded two-second observation.

Every binary, CMakeCache, actual core/executable flags.make, local runtime library,
input, card and runner is hashed. Core flags must match each declared revision;
the actual binary revision is checked in its manifest after launch. Actual
compiler flags must match except revision.
Each run rechecks frozen files, has a 180 s timeout plus 10 s kill grace; the whole
campaign has a 1800 s budget. Binding and runtime dependency checks occur first.
LAUNCHED.json records process launch only. STARTED_AND_COMPLETED additionally
requires real APRL startup/physical-stop messages, compiled revision, 4 MPI ranks,
4 threads, scientific steps, committed H5 outputs and case-specific coverage.
There is no automatic rerun, resume, submission, system setting or production
change. These are correctness fixtures, not performance or radiation convergence
measurements. Use the unchanged three generated workload cards for speed claims.
