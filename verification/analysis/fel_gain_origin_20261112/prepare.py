"""Freeze same-input GENESIS controls; never edit the active APRL campaigns."""
import hashlib
import json
import re
import shutil
import sys
from pathlib import Path

W = Path('/public/home/wfollowm/APRL_Workspace')
PIC = W/'results/sparc_sase_11segment_2slip_20261111_r1'
ADI = W/'results/binomial3_fel_20261106/resonance_origin_20261108/genesis_reassessment_r1/adi_only_control_r3'
INSTALLED = Path('/public/software/Genesis-1.3-Version4-4.6.4/build/genesis4')


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda: f.read(4*1024*1024), b''):
            h.update(b)
    return h.hexdigest()


def save(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')


def main(root):
    root = Path(root).resolve()
    root.mkdir(exist_ok=False)
    shutil.copytree(Path(__file__).parent, root/'scripts',
                    ignore=shutil.ignore_patterns('__pycache__'))
    p = json.loads((PIC/'plan.json').read_text())
    fixed = ADI/'build/fixed/genesis4'
    assert sha(fixed) == json.loads((ADI/'build/fixed/BUILD_COMPLETE.json').read_text())['binary_sha256']
    assert sha(PIC/'particles.h5') == json.loads((PIC/'INPUT_READY.json').read_text())['input_sha256']
    cases = []
    for name, sc, binary, dz in [
        ('sc-off-fixed', False, fixed, .007),
        ('longrange-fixed', True, fixed, .007),
        ('longrange-installed', True, INSTALLED, .007),
        ('sc-off-halfstep', False, fixed, .0035),
    ]:
        d = root/name
        d.mkdir()
        (d/'common.par.h5').symlink_to(PIC/'genesis/common.par.h5')
        original_lattice = (PIC/'genesis/reference.lat').read_text()
        seen = [0]
        def mark(match):
            seen[0] += 1
            return match.group(0)+(',AUDIT' if seen[0] in [1,6] else '')
        lattice = re.sub('R0,R1,R2,R3,U,R3,R2,R1,R0',mark,original_lattice)
        assert seen[0] == 11
        assert lattice.replace(',AUDIT','') == original_lattice
        (d/'reference.lat').write_text('AUDIT: MARKER={dumpfield=1};\n'+lattice)
        ef = ('&efield\nlongrange=true\nrmax=2\nnz=1\nnphi=0\nngrid=100\n&end\n'
              if sc else '')
        deck = f'''&setup
rootname=reference
lattice=reference.lat
beamline=FEL
lambda0={p['physics']['wavelength_m']:.17g}
gamma0={p['physics']['gamma']:.17g}
delz={dz:.17g}
shotnoise=false
one4one=false
npart=1024
nbins=16
beam_global_stat=true
field_global_stat=true
&end
&importbeam
file=common.par.h5
time=true
&end
&write
beam=imported
&end
&field
power=0
dgrid=.0012
ngrid=161
&end
{ef}&track
zstop={p['geometry']['field_planes'][5][1]:.17g}
output_step=1
&end
&write
beam=exit
&end
'''
        (d/'reference.in').write_text(deck)
        cases.append(dict(name=name, directory=str(d), binary=str(binary),
                          binary_sha256=sha(binary), longrange=sc,
                          short_range_effective=False, delz_m=dz, mpi=32,
                          timeout_seconds=1800))
    plan = dict(root=str(root), pic_campaign=str(PIC), cases=cases,
                exact_input=str(PIC/'particles.h5'),
                exact_input_sha256=sha(PIC/'particles.h5'),
                genesis_input=str(PIC/'genesis/common.par.h5'),
                genesis_input_sha256=sha(PIC/'genesis/common.par.h5'),
                active_slices=p['active_slices'], npart=1024, nbins=16,
                original_plan=p,
                research_question='Separate start-up source strength, linear growth, loading correlations and reference-code implementation effects.',
                current_APRL_jobs_untouched=[19743,19744],
                scope='Same frozen 6D input, seed, current and six-segment phase-matched lattice. No new APRL propagation. Four GENESIS controls use one exclusive 128-core allocation, 32 MPI ranks per control.',
                comparison_limits=['GENESIS integrated radiation power versus cached PIC full-record band power are separately labelled.',
                    'Actual input-particle or imported-cell grouping is retained for BF; physical transport can change Eulerian cell membership.',
                    'One seed only; no stochastic precision claim.',
                    'The installed shortRange() call is commented out; nz=1 is not an enabled microbunch space-charge model.'])
    save(root/'PLAN.json', plan)
    files = [root/'PLAN.json', PIC/'particles.h5', PIC/'genesis/common.par.h5', fixed, INSTALLED]
    files += [q for q in (root/'scripts').iterdir() if q.is_file()]
    files += [Path(c['directory'])/name for c in cases for name in ['reference.in','reference.lat']]
    save(root/'FROZEN.json', {str(q):sha(q) for q in files})
    print(root)


if __name__ == '__main__':
    main(sys.argv[1])
