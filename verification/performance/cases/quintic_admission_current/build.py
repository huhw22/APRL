"""Build an isolated admission repair and run meaningful source regressions."""
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys


def main(root):
    root = Path(root).resolve()
    assert os.environ.get('SLURM_JOB_ID')
    plan = json.loads((root / 'plan.json').read_text())
    source, build = Path(plan['source']), Path(plan['build'])
    runtime = Path('/public/home/wfollowm/aprl_dev/aprl-runtime-lab')
    cmake = '/public/home/wfollowm/software/cmake-3.31.12/bin/cmake'
    frozen = json.loads((root / 'FROZEN.json').read_text())

    def verify():
        for path, sha in frozen.items():
            assert hashlib.sha256(Path(path).read_bytes()).hexdigest() == sha, path

    def run(cmd, name, limit):
        print('RUN', name, flush=True)
        with (root / name).open('x') as log:
            subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT,
                           timeout=limit, check=True)

    verify()
    assert not build.exists()
    run([cmake, '-S', str(source), '-B', str(build),
         '-DCMAKE_BUILD_TYPE=Release', '-DCMAKE_CXX_COMPILER=/usr/bin/mpic++',
         '-DCMAKE_PREFIX_PATH=' + str(runtime), '-DHDF5_ROOT=' + str(runtime),
         '-DBUILD_TESTING=ON', '-DAPRL_ENABLE_OPENMP=ON',
         '-DAPRL_ENABLE_MOVING_WINDOW=ON', '-DAPRL_REQUIRE_PARALLEL_HDF5=ON'],
        'configure.log', 120)
    run([cmake, '--build', str(build), '--parallel', '12'], 'compile.log', 900)
    common = ['/usr/bin/mpic++', '-std=c++11', '-O3', '-fopenmp',
              '-DAPRL_HAS_OPENMP=1', '-DAPRL_HAS_MOVING_WINDOW=1',
              '-DOMPI_SKIP_MPICXX', '-DMPICH_SKIP_MPICXX',
              '-I' + str(source / 'src'), '-I' + str(runtime / 'include')]
    libraries = [str(build / 'libdirect_eb_core.a'),
                 str(runtime / 'lib/libhdf5.so'),
                 str(runtime / 'lib/libyaml-cpp.so.0.8.0'), '-lfftw3',
                 '-Wl,-rpath,' + str(runtime / 'lib')]
    cases = source / 'verification/performance/cases'
    for name, path in [
            ('admission_checks', cases / 'quintic_admission_current/checks.cpp'),
            ('quintic_checks', cases / 'longitudinal_main_flow/checks.cpp'),
            ('exit_checks', cases / 'source_coupling_integration/exit_checks.cpp')]:
        run(common + [str(path)] + libraries + ['-o', str(build / name)],
            name + '-compile.log', 120)
    run([str(build / 'exit_checks')], 'exit-checks.json', 60)
    for ranks in [1, 2, 4]:
        for name in ['admission_checks', 'quintic_checks']:
            run(['mpiexec', '-n', str(ranks), str(build / name),
                 str(root / (name + '-' + str(ranks)))],
                name + '-' + str(ranks) + '.json', 120)
    run(['ctest', '--test-dir', str(build), '-L', 'required',
         '--output-on-failure', '-j', '1'], 'required.log', 600)
    verify()
    result = dict(passed=True, job=os.environ['SLURM_JOB_ID'],
                  binary=str(build / 'aprl'),
                  binary_sha256=hashlib.sha256((build / 'aprl').read_bytes()).hexdigest(),
                  snapshot=str(source), production_acceptance=False)
    (root / 'COMPLETE.json').write_text(json.dumps(result, indent=2) + '\n')
    print('COMPLETE', flush=True)


if __name__ == '__main__':
    main(sys.argv[1])
