"""Freeze code and build helpers; never modify the accepted R5 installation."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--source', required=True, type=Path)
    parser.add_argument('--root', required=True, type=Path)
    parser.add_argument('--build', required=True, type=Path)
    args = parser.parse_args()
    source, root, build = args.source.resolve(), args.root.resolve(), args.build.resolve()
    root.mkdir(parents=True, exist_ok=False)
    snapshot = root / 'source'
    shutil.copytree(source, snapshot, ignore=shutil.ignore_patterns('.git', '__pycache__', '*.pyc'))
    scripts = root / 'scripts'
    scripts.mkdir()
    for name in ['build.py', 'build.sbatch']:
        shutil.copy2(source / 'verification/performance/cases/oneway_emission_build' / name,
                     scripts / name)
    plan = dict(source=str(snapshot), build=str(build), source_worktree=str(source),
                parent_commit='8760cfb',
                source_commit=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=source, text=True).strip(),
                source_status=subprocess.check_output(['git', 'status', '--short'], cwd=source, text=True),
                accepted_binary='/public/home/wfollowm/APRL_Workspace/build/mpi-particle-balance-r5-20261106/aprl',
                modules={
                    'oneway_checks': 'verification/performance/cases/oneway_emission/checks.cpp',
                    'binomial3_checks': 'verification/performance/cases/binomial3_control/checks.cpp',
                    'binomial3_admission_checks': 'verification/performance/cases/binomial3_control/admission_checks.cpp',
                    'compact_checks': 'verification/performance/cases/mpi_particle_balance/compact_checks.cpp',
                })
    assert all((snapshot / item).is_file() for item in plan['modules'].values())
    (root / 'plan.json').write_text(json.dumps(plan, indent=2) + '\n')
    frozen = {str(p): hashlib.sha256(p.read_bytes()).hexdigest()
              for p in sorted(root.rglob('*')) if p.is_file()}
    (root / 'FROZEN.json').write_text(json.dumps(frozen, indent=2) + '\n')
    print(json.dumps({'root': str(root), 'frozen_files': len(frozen), 'source_commit': plan['source_commit']}))


if __name__ == '__main__':
    main()
