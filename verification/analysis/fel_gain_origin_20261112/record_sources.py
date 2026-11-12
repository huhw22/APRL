"""Record the local evidence actually inspected, without modifying any source."""
import json
import subprocess
import sys
from pathlib import Path

from prepare import save, sha

R=Path(sys.argv[1]).resolve()
W=Path('/public/home/wfollowm/APRL_Workspace')
G=Path('/public/software/Genesis-1.3-Version4-4.6.4')
F=Path('/public/home/wfollowm/Puffin')
paths=[G/name for name in ['src/Loading/QuietLoading.cpp','src/Loading/ShotNoise.cpp',
    'src/Loading/ImportBeam.cpp','src/IO/readBeamHDF5.cpp','src/Core/BeamSolver.cpp',
    'src/Core/EFieldSolver.cpp','src/Core/FieldSolver.cpp']]
paths += [F/name for name in ['README.md','BUILD.md','CMakeLists.txt',
    'puffin/lib/common/puffin_fftwInfo.f90']]
paths += [W/name for name in [
    'results/sparc_sase_11segment_2slip_20261111_r1/NOISE_VALIDATION.json',
    'results/sparc_sase_11segment_2slip_20261111_r1/scripts/prepare.py',
    'results/sparc_sase_11segment_2slip_20261111_r1/analysis-energy-binomial3/theory-literature-comparison/RESULT.json',
    'results/energy_critical_control_20261105/bf_origin_review/loading_correlation_test_20261105/completed-review/README.md',
    'results/conversation_handoff_20261109/README.md',
    'results/sparc_sase_comparison_20261103/energy_review_20261104/THEORY_GAP.md']]
paths += [R/name for name in ['README.md','STATIC_SAMPLE.json','input_family_harmonics.png']]
paths += [p for p in Path(__file__).parent.iterdir() if p.is_file()]
save(R/'SOURCE_IDENTITIES.json',dict(
    source_hashes={str(p):sha(p) for p in paths},
    puffin_git_revision=subprocess.check_output(['git','-C',str(F),'rev-parse','HEAD'],text=True).strip(),
    puffin_git_status=subprocess.check_output(['git','-C',str(F),'status','--short'],text=True),
    installed_and_ADI_fixed_BeamSolver_identical=(sha(G/'src/Core/BeamSolver.cpp')==sha(W/'results/binomial3_fel_20261106/resonance_origin_20261108/genesis_reassessment_r1/adi_only_control_r3/source/fixed/src/Core/BeamSolver.cpp')),
    source_changes_this_task='Only new verification/analysis/fel_gain_origin_20261112 templates; no physical kernels, default, inputs or active-job cards changed.',
    chart_inspected=True,
    validation_scope='Central-cell static audit and source/caller inspection; full-input and new GENESIS controls pending.'))
print('Inspected-source identities saved')
