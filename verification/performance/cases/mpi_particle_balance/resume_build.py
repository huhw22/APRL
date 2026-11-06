"""Resume acceptance after a test-only correction; compiled core stays immutable."""
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys


def main(root):
    root=Path(root).resolve();plan=json.loads((root/'plan.json').read_text())
    parent=Path(plan['compiled_record']);source=Path(plan['source']);build=Path(plan['build'])
    frozen=json.loads((root/'FROZEN.json').read_text())
    sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
    def verify():
        assert all(sha(p)==v for p,v in frozen.items()),'Frozen resumption changed'
        compiled=json.loads((parent/'SOURCE_BUILD.json').read_text())['source_hashes']
        assert all(sha(source/p)==v for p,v in compiled.items()),'Core changed: rebuild instead of resuming'
        assert sha(build/'aprl')==plan['binary_sha256']
    def run(command,name,limit):
        print('RUN',name,flush=True)
        with (root/name).open('x') as stream:subprocess.run(command,stdout=stream,stderr=subprocess.STDOUT,timeout=limit,check=True)
    verify();runtime=Path('/public/home/wfollowm/aprl_dev/aprl-runtime-lab')
    compiled_modules={}
    for test in ['balance_checks','compact_checks']:
        filename='checks.cpp' if test=='balance_checks' else 'compact_checks.cpp'
        if not (source/'verification/performance/cases/mpi_particle_balance'/filename).exists():continue
        module=build/(test+'-'+root.name);compiled_modules[test]=module
        run(['/usr/bin/mpic++','-std=c++11','-O3','-fopenmp','-DAPRL_HAS_OPENMP=1',
        '-DAPRL_HAS_MOVING_WINDOW=1','-DOMPI_SKIP_MPICXX','-DMPICH_SKIP_MPICXX',
        '-I'+str(source/'src'),'-I'+str(runtime/'include'),
        str(source/'verification/performance/cases/mpi_particle_balance'/filename),
        str(build/'libdirect_eb_core.a'),str(runtime/'lib/libhdf5.so'),
        str(runtime/'lib/libyaml-cpp.so.0.8.0'),'-lfftw3','-Wl,-rpath,'+str(runtime/'lib'),
        '-o',str(module)],test+'-compile.log',120)
    modules=['binomial3_checks','binomial3_admission_checks','notch_checks','substep_checks','admission_checks','quintic_checks']
    for ranks in [1,2,4]:
        for name in modules:
            run(['mpiexec','-n',str(ranks),str(build/name),str(root/(name+'-'+str(ranks)))],name+'-'+str(ranks)+'.json',180)
            assert json.loads((root/(name+'-'+str(ranks)+'.json')).read_text())['passed']
    for test,module in compiled_modules.items():
        for ranks in ([1,2,3,4,8] if test=='balance_checks' else [1,2,4,8]):
            name=test+'-'+str(ranks)+'.json'
            run(['mpiexec','-n',str(ranks),str(module),str(root/(test+'-'+str(ranks)))],name,180)
            assert json.loads((root/name).read_text())['passed']
    run(['ctest','--test-dir',str(build),'-L','required','--output-on-failure','-j','1'],'required.log',600)
    verify()
    (root/'COMPLETE.json').write_text(json.dumps(dict(passed=True,job=os.environ['SLURM_JOB_ID'],
        binary=str(build/'aprl'),binary_sha256=sha(build/'aprl'),compiled_record=str(parent),
        snapshot=str(source),production_acceptance=False,
        test_only_correction=plan['test_only_correction'],core_recompiled=False),indent=2)+'\n')
    print('COMPLETE',flush=True)


if __name__=='__main__':main(sys.argv[1])
