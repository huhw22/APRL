"""Build the isolated MPI particle-balance candidate and old-mode regressions."""
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys


def main(root):
    root = Path(root).resolve()
    assert os.environ.get('SLURM_JOB_ID'), 'Build inside the bounded Slurm allocation'
    plan = json.loads((root / 'plan.json').read_text())
    source, build = Path(plan['source']), Path(plan['build'])
    runtime = Path('/public/home/wfollowm/aprl_dev/aprl-runtime-lab')
    cmake = '/public/home/wfollowm/software/cmake-3.31.12/bin/cmake'
    frozen = json.loads((root / 'FROZEN.json').read_text())

    def verify():
        for path, digest in frozen.items():
            assert hashlib.sha256(Path(path).read_bytes()).hexdigest() == digest, path

    def run(command, name, limit):
        print('RUN', name, flush=True)
        with (root / name).open('x') as stream:
            subprocess.run(command, stdout=stream, stderr=subprocess.STDOUT,
                           timeout=limit, check=True)

    verify()
    assert not build.exists(), 'Refusing to reuse or replace a build directory'
    run([cmake, '-S', str(source), '-B', str(build), '-DCMAKE_BUILD_TYPE=Release',
         '-DCMAKE_CXX_COMPILER=/usr/bin/mpic++', '-DCMAKE_PREFIX_PATH=' + str(runtime),
         '-DHDF5_ROOT=' + str(runtime), '-DBUILD_TESTING=ON', '-DAPRL_ENABLE_OPENMP=ON',
         '-DAPRL_ENABLE_MOVING_WINDOW=ON', '-DAPRL_REQUIRE_PARALLEL_HDF5=ON'],
        'configure.log', 120)
    run([cmake, '--build', str(build), '--parallel', '12'], 'compile.log', 900)
    common = ['/usr/bin/mpic++', '-std=c++11', '-O3', '-fopenmp', '-DAPRL_HAS_OPENMP=1',
              '-DAPRL_HAS_MOVING_WINDOW=1', '-DOMPI_SKIP_MPICXX', '-DMPICH_SKIP_MPICXX',
              '-I' + str(source / 'src'), '-I' + str(runtime / 'include')]
    libraries = [str(build / 'libdirect_eb_core.a'), str(runtime / 'lib/libhdf5.so'),
                 str(runtime / 'lib/libyaml-cpp.so.0.8.0'), '-lfftw3',
                 '-Wl,-rpath,' + str(runtime / 'lib')]
    cases = source / 'verification/performance/cases'
    modules = {
        'binomial3_checks': cases / 'binomial3_control/checks.cpp',
        'binomial3_admission_checks': cases / 'binomial3_control/admission_checks.cpp',
        'notch_checks': cases / 'energy_critical_control/checks.cpp',
        'substep_checks': cases / 'energy_critical_control/substep_checks.cpp',
        'admission_checks': cases / 'quintic_admission_current/checks.cpp',
        'quintic_checks': cases / 'longitudinal_main_flow/checks.cpp',
        'balance_checks': cases / 'mpi_particle_balance/checks.cpp',
        'compact_checks': cases / 'mpi_particle_balance/compact_checks.cpp',
    }
    for name, path in modules.items():
        assert path.is_file(), str(path)
        run(common + [str(path)] + libraries + ['-o', str(build / name)],
            name + '-compile.log', 120)
    for ranks in [1, 2, 4]:
        for name in modules:
            run(['mpiexec', '-n', str(ranks), str(build / name),
                 str(root / (name + '-' + str(ranks)))],
                name + '-' + str(ranks) + '.json', 180)
            assert json.loads((root / (name + '-' + str(ranks) + '.json')).read_text())['passed']
    for ranks in [3, 8]:
        name = 'balance_checks'
        run(['mpiexec', '-n', str(ranks), str(build / name),
             str(root / (name + '-' + str(ranks)))],
            name + '-' + str(ranks) + '.json', 180)
        assert json.loads((root / (name + '-' + str(ranks) + '.json')).read_text())['passed']
    name='compact_checks'
    run(['mpiexec','-n','8',str(build/name),str(root/(name+'-8'))],name+'-8.json',180)
    assert json.loads((root/(name+'-8.json')).read_text())['passed']
    run(['ctest', '--test-dir', str(build), '-L', 'required',
         '--output-on-failure', '-j', '1'], 'required.log', 600)
    verify()
    complete = dict(passed=True, job=os.environ['SLURM_JOB_ID'], binary=str(build / 'aprl'),
        binary_sha256=hashlib.sha256((build / 'aprl').read_bytes()).hexdigest(),
        snapshot=str(source), production_acceptance=False,
        requested_optimization='MPI particle-weighted nonuniform integer z slabs',
        modules=list(modules), module_ranks={name: ([1,2,3,4,8] if name=='balance_checks' else
            ([1,2,4,8] if name=='compact_checks' else [1,2,4])) for name in modules},
        required_regressions=True,
        production_default_replaced=False)
    (root / 'COMPLETE.json').write_text(json.dumps(complete, indent=2) + '\n')
    print('COMPLETE', flush=True)


if __name__ == '__main__':
    main(sys.argv[1])
