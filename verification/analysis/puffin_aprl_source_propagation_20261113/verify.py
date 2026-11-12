"""Independent delivery checks on closed caches and frozen identities; no PIC run."""
import ast
import csv
import hashlib
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import simpson

W = Path('/public/home/wfollowm/APRL_Workspace')
REPO = W / 'source/aprl'
HERE = Path(__file__).resolve().parent


def read(path):
    return json.loads(Path(path).read_text())


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def need(value, message):
    if not value:
        raise RuntimeError(message)


def main(root):
    root = Path(root).resolve()
    control = root / 'physical-route-control-r1'
    need(read(control / 'COMPLETE.json')['complete'], 'Control not closed')
    frozen = {}
    for manifest in [root / 'FROZEN.json', control / 'FROZEN.json',
                     W / 'results/puffin_same_sample_20261113_r3/FROZEN.json',
                     W / 'results/sparc_random_startup_20261112_r2/w64-a/FROZEN.json']:
        for path, expected in read(manifest).items():
            need(path not in frozen or frozen[path] == expected, 'Manifest conflict: ' + path)
            frozen[path] = expected
    plan = read(W / 'results/sparc_random_startup_20261112_r2/PLAN.json')
    frozen[plan['pic_binary']] = plan['pic_binary_sha256']
    for path, expected in frozen.items():
        need(sha(path) == expected, 'Frozen identity mismatch: ' + path)

    # All local headers used by the real-kernel probe must match the frozen
    # source from which the observed production PIC binary was built.
    dep = subprocess.run(['c++', '-std=c++11', '-I' + str(REPO / 'src'), '-MM',
                          str(HERE / 'cowan_wave_probe.cpp'), str(REPO / 'src/eb_field.cpp'),
                          str(REPO / 'src/boostframe.cpp')], check=True, capture_output=True,
                         text=True).stdout.replace('\\\n', ' ')
    dependencies = sorted({Path(token) for line in dep.splitlines()
                           for token in line.split(':', 1)[1].split()
                           if token.startswith(str(REPO / 'src') + '/')})
    production = W / 'results/longitudinal_shape_lab_balance_20261110_r4'
    production_hashes = read(production / 'FROZEN.json')
    source_identity = {}
    for path in dependencies:
        original = production / 'source' / path.relative_to(REPO)
        expected = production_hashes[str(original)]
        need(sha(original) == expected and sha(path) == expected,
             'Probe source differs from observed production build: ' + str(path))
        source_identity[str(path)] = expected

    syntax = []
    for path in sorted(HERE.glob('*.py')):
        ast.parse(path.read_text(), filename=str(path))
        syntax.append(path.name)
    for path in sorted(HERE.glob('*.sbatch')):
        subprocess.run(['bash', '-n', str(path)], check=True)
        syntax.append(path.name)

    summary = read(root / 'RESULT.json')
    stages = {row['stage']: row['band_energy_J'] for row in summary['field_measurement_stages']}
    with np.load(control / 'analysis/fields.npz') as cache:
        freq = cache['frequency_Hz']
        need(np.all(np.diff(freq) > 0), 'Frequency order')
        for key in ['Puffin', 'PIC', 'Puffin_same_pixels', 'Puffin_transverse_detector']:
            need(np.isfinite(cache[key]).all(), 'Nonfinite field cache: ' + key)
        xp, yp = cache['PIC_x_m'], cache['PIC_y_m']
        dx, dy = xp[1] - xp[0], yp[1] - yp[0]
        ix = np.flatnonzero(abs(xp) < .00056 - 1e-12)
        iy = np.flatnonzero(abs(yp) < .00056 - 1e-12)
        with np.load(control / 'phase-analysis/bilinear_field.npz') as other:
            need(np.array_equal(freq, other['frequency_Hz']), 'Field frequency mismatch')
            p = cache['Puffin']
            x, y = cache['x_m'], cache['y_m']
            # Independent literal half-grid four-node evaluation.
            i = np.searchsorted(x, xp) - 1
            j = np.searchsorted(y, yp) - 1
            literal = sum(p[:, j[:, None] + b, i[None, :] + a, :] / 4
                          for a in (0, 1) for b in (0, 1))
            v = other['Puffin_bilinear']
            relative = float(np.max(abs(literal - v)) / np.max(abs(v)))
            need(relative < 3e-14, 'Literal interpolation disagrees')
            need(np.isfinite(v).all(), 'Nonfinite interpolated field')

        # Trapz on a twice refined frequency grid would require new data;
        # instead use the actual shared Simpson frequency quadrature, but
        # independent explicit full-pixel masks rather than partial weights.
        def energy(a):
            values = abs(a[:, iy[:, None], ix[None, :], :]) ** 2
            return float(simpson(values.sum(axis=(1, 2, 3)) * dx * dy, x=freq))

        pic_energy = energy(cache['PIC'])
        puf_energy = energy(literal)
        need(abs(pic_energy / stages['PIC stored detector'] - 1) < 2e-13, 'PIC integral mismatch')
        need(abs(puf_energy / stages['Both-axis linear sampling'] - 1) < 2e-13, 'Puffin integral mismatch')

    need(read(control / 'PHYSICAL_ROUTE_VALIDATED.json')['native_all_seven_columns_bitwise'],
         'Native import not verified')
    need('ROUTE_PASS' in (control / 'route.log').read_text(), 'Runtime route not verified')
    with (root / 'COWAN_ACTUAL_KERNEL_PROBE.csv').open() as stream:
        modes = list(csv.DictReader(stream))
    need(len(modes) == 6 and all(float(row['E_max_error']) < 3e-14 and
                                float(row['cB_max_error']) < 4e-15 for row in modes),
         'Actual Cowan kernel probe failed')
    need(summary['new_seeds'] == 0 and summary['new_PIC_runs'] == 0 and
         summary['one_corrected_Puffin_run'] == 19770 and not summary['physical_acceptance'],
         'Scope or acceptance changed')
    result = dict(passed=True, frozen_identities=len(frozen), actual_production_source_matches=source_identity,
                  syntax_checked=syntax, actual_Cowan_modes=len(modes),
                  literal_interpolation_relative_error=relative,
                  independent_pixel_integrals_J=dict(PIC=pic_energy, Puffin_bilinear=puf_energy),
                  new_seeds=0, new_PIC_runs=0, physical_acceptance=False,
                  scope='Hashes of frozen inputs/programs/scripts; closed small field caches. Raw large output H5 not rehashed.')
    (root / 'VERIFICATION.json').write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == '__main__':
    main(sys.argv[1])
