"""Prepare four original 257-step drift cards for program/diagnostic regression."""
import argparse
import copy
import json
from pathlib import Path
import shutil

import yaml
from prepare import require, save, sha


def prepare_smoke(workspace, output, binary, build_record):
    workspace, output, binary, build_record = (Path(x).resolve() for x in (workspace, output, binary, build_record))
    output.mkdir(parents=True, exist_ok=False)
    scripts = output/'scripts'
    scripts.mkdir()
    template = Path(__file__).resolve().parent
    for name in ['prepare.py', 'prepare_smoke.py', 'run_smoke.py', 'smoke.sbatch']:
        shutil.copy2(template/name, scripts/name)
    for folder, name in [('space_charge_physics', 'launcher.py'),
                         ('source_coupling_integration', 'binding.py'),
                         ('source_coupling_integration', 'analyze_closed.py')]:
        shutil.copy2(template.parent/folder/name, scripts/name)
    for name in ['report.py', 'review_discrete.py']:
        shutil.copy2(template.parent/'discrete_energy_ledger'/name,
                     scripts/('report_discrete.py' if name=='report.py' else name))
    original = workspace/'results/longitudinal_main_flow_20261103/flow-r2/drift-linear'
    require(json.loads((original/'RESULT.json').read_text())['program_integrity_pass'], 'Baseline small card is not closed')
    old = workspace/'build/sparc-sase-r1-20261103/aprl'
    cases = []
    for name, program, diagnostic, shape in [
            ('old', old, False, 'linear'), ('default', binary, False, 'linear'),
            ('substeps', binary, True, 'linear'), ('notch', binary, True, 'linear-notch')]:
        directory = output/name
        directory.mkdir()
        card = yaml.safe_load((original/'run.yaml').read_text())
        beam = copy.deepcopy(card['beam'])
        def rewrite(value):
            if isinstance(value, dict):
                return {k:rewrite(v) for k,v in value.items()}
            if isinstance(value, list):
                return [rewrite(v) for v in value]
            return value.replace(str(original), str(directory)) if isinstance(value, str) else value
        card = rewrite(card)
        card['beam'] = beam
        card['mesh']['longitudinal_shape'] = shape
        if name != 'old':
            card['detectors']['substep_particle_planes'] = diagnostic
        (directory/'run.yaml').write_text(yaml.safe_dump(card, sort_keys=False))
        for folder in ['detectors', 'trajectory', 'ledger']:
            (directory/folder).mkdir()
        cases.append(dict(name=name, directory=str(directory), binary=str(program),
            substeps=diagnostic, shape=shape, limit_seconds=300))
    plan = dict(workspace=str(workspace), directory=str(output), binary=str(binary),
        build_record=str(build_record), expected_particles=31675, expected_steps=257,
        cases=cases, physical_acceptance=False,
        scope='Original drift-card software/diagnostic regression, not SPARC BF physics.')
    save(output/'plan.json', plan)
    paths = [old, output/'plan.json', *scripts.iterdir()]
    for case in cases:
        path = Path(case['directory'])/'run.yaml'
        paths.extend([path, Path(yaml.safe_load(path.read_text())['beam']['input']['file'])])
    save(output/'FROZEN_INPUTS.json', {str(p):sha(p) for p in paths if p.is_file()})
    print(output)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('workspace', type=Path)
    p.add_argument('output', type=Path)
    p.add_argument('--binary', type=Path, required=True)
    p.add_argument('--build-record', type=Path, required=True)
    a = p.parse_args()
    prepare_smoke(a.workspace, a.output, a.binary, a.build_record)
