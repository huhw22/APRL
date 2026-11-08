"""Build an immutable one-way-source validation snapshot inside Slurm."""
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time


def main(root):
    root = Path(root).resolve()
    assert os.environ.get('SLURM_JOB_ID'), 'Use the bounded compute allocation'
    plan = json.loads((root / 'plan.json').read_text())
    source, build = Path(plan['source']), Path(plan['build'])
    frozen = json.loads((root / 'FROZEN.json').read_text())
    runtime = Path('/public/home/wfollowm/aprl_dev/aprl-runtime-lab')
    cmake = '/public/home/wfollowm/software/cmake-3.31.12/bin/cmake'

    def verify():
        for path, digest in frozen.items():
            assert hashlib.sha256(Path(path).read_bytes()).hexdigest() == digest, path

    def run(command, name, limit):
        print('RUN', name, flush=True)
        with (root / name).open('x') as stream:
            subprocess.run(command, stdout=stream, stderr=subprocess.STDOUT,
                           timeout=limit, check=True)

    started = time.time()
    verify()
    assert not build.exists(), 'Refusing to overwrite an existing binary/build'
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
    for name, relative in plan['modules'].items():
        run(common + [str(source / relative)] + libraries + ['-o', str(build / name)],
            name + '-compile.log', 120)
        for ranks in [1, 2, 4]:
            target = root / (name + '-' + str(ranks))
            run(['mpiexec', '-n', str(ranks), str(build / name), str(target)],
                target.name + '.json', 180)
            data = json.loads((root / (target.name + '.json')).read_text())
            assert data['passed'], (name, ranks, data)
    run(['ctest', '--test-dir', str(build), '-L', 'required',
         '--output-on-failure', '-j', '1'], 'required.log', 600)
    verify()
    result = dict(passed=True, job=os.environ['SLURM_JOB_ID'],
                  elapsed_s=time.time() - started, binary=str(build / 'aprl'),
                  binary_sha256=hashlib.sha256((build / 'aprl').read_bytes()).hexdigest(),
                  source_snapshot=str(source), parent_commit=plan['parent_commit'],
                  modules=plan['modules'], module_ranks=[1, 2, 4],
                  required_regressions=True, production_default_replaced=False,
                  physical_acceptance=False)
    (root / 'COMPLETE.json').write_text(json.dumps(result, indent=2) + '\n')
    print('COMPLETE', flush=True)


if __name__ == '__main__':
    main(sys.argv[1])
