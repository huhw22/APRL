# Low-cost lab energy/momentum algorithm exploration

Experimental tools only; production src/defaults remain unchanged.
See METHOD.md for the derivation, the unchanged accounting definitions,
the actual magnetic-gather change, and the unaccepted field candidates.

prepare.py takes --workspace, --output, --build and requires fresh directories.
It materializes the verified coupled_sbp_r2 source/archive, changes the two
real magnetic gather calls, and rebuilds only source_coupling.cpp.o. Bx/By
use two MPI collar planes with the original message count/tags. The original
constraint observer and real global-face clipping remain; full-record physical
particle conservation is not claimed. run.py checks scalar identities/binding,
runs the two short cards, and uses the same numerical budget formulas with
honest metadata indicating an intentional candidate algorithm change.

timing.py --workspace ... --output ... creates same-allocation ABBA cards.
It requires completed science, disables heavy ledgers/observers equally, and
keeps detectors and physics settings identical between the four runs.
Two samples per method and noisy baseline cannot establish production speedup.

spectral_screen.py and screen.sbatch run separate mathematical/Fourier and
local-stencil field prototypes. They contain no PIC particles, CPML, moving
window or reference injector; energy identities do not establish physics.
compare.py reads closed baseline/candidate results and uses one fixed BASELINE
particle-change denominator per case for all displayed comparison curves.

Artifacts: /public/home/wfollowm/APRL_Workspace/results/lab_coupling_algorithm_20261102/
Actual complete science: r3 (19230); theory_r1 (19231); timing_r1 (19232).
Failed r1 preparation and r2/19229 collar-range attempt remain as evidence.
No full simulation, production replacement, remote push or double-node claim.
