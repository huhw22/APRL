"""Read only closed short-flow results; linear science must remain bitwise equal."""
import json
from pathlib import Path

import h5py
import numpy as np
from scipy.constants import c

from flow_prepare import require, save
from analyze_closed import verify, bitwise
import report_discrete


def analyze_case(directory, particles=31675, expected_steps=257):
    directory = Path(directory)
    closed = verify(directory)
    rows, sums, ranks, steps = report_discrete.load(directory / 'ledger',
                                                  directory / 'ledger/energy-ledger.h5')
    require(steps == expected_steps, 'Unexpected short-flow duration')
    require(np.all(rows['active_macroparticles'] + rows['removed_macroparticles'] == particles),
            'Cohort accounting differs')
    _, _, strict = report_discrete.calculate(rows, sums, steps)
    for p in (directory / 'detectors').glob('*.h5'):
        with h5py.File(p) as f:
            if 'particle_plane' not in f:
                continue
            a = f['particle_plane/records'][:]
            require(len(a) == len(np.unique(a['particle_id'])) == particles,
                    'Missing/duplicate particle-plane event: ' + str(p))
    for p in (directory / 'reference').glob('source-coupling-rank-*.csv'):
        require(p.with_suffix('.complete').read_text().split() == ['1', str(steps)],
                'Unclosed coupled steps')
        a = np.genfromtxt(p, delimiter=',', names=True)
        require(np.all(a['self_force_defect'] < 1e-13) and
                np.all(a['relative_fixed_point_work_defect'] < 1e-8), 'Coupling did not converge')
    return dict(program_integrity_pass=True, physical_acceptance=False,
                steps=steps, ranks=ranks, closed_scientific_files=closed,
                strict_ledger=strict,
                final_field_collocated_J=float(rows['field_energy_interior_J'][-1]),
                final_frozen_constraint_work_J=float(rows['frozen_constraint_work_J'][-1]))


def compare(root):
    root = Path(root)
    linear = bitwise(root / 'linear-old', root / 'linear-new')
    a, b = root / 'quintic-old', root / 'quintic-new'
    with h5py.File(a / 'detectors/exit-field.h5') as fa, h5py.File(b / 'detectors/exit-field.h5') as fb:
        ga, gb = fa['field_plane'], fb['field_plane']
        require(ga['time_s'][:].tobytes() == gb['time_s'][:].tobytes(), 'Changed output timestamps')
        require(ga['fields'].shape == gb['fields'].shape, 'Changed field output geometry')
        error = {'E': np.longdouble(0), 'cB': np.longdouble(0)}
        scale = {'E': np.longdouble(0), 'cB': np.longdouble(0)}
        for start in range(0, len(ga['fields']), 8):
            x, y = ga['fields'][start:start+8], gb['fields'][start:start+8]
            # Field records store separate vector members, not six scalars.
            # Use E and cB for a common unit in the combined relative norm.
            for key, member, factor in [('E', 'electric_V_per_m', 1.),
                                        ('cB', 'magnetic_T', c)]:
                xx, yy = x[member].astype(np.longdouble), y[member].astype(np.longdouble)
                error[key] += np.sum(((xx-yy)*np.longdouble(factor))**2)
                scale[key] += np.sum((xx*np.longdouble(factor))**2)
    result = dict(linear_scientific_bitwise_pass=True, linear_files=linear,
                  quintic_raw_field_relative_L2=float(np.sqrt(sum(error.values()) / max(sum(scale.values()), np.longdouble('1e-300')))),
                  quintic_raw_E_relative_L2=float(np.sqrt(error['E'] / max(scale['E'], np.longdouble('1e-300')))),
                  quintic_raw_B_relative_L2=float(np.sqrt(error['cB'] / max(scale['cB'], np.longdouble('1e-300')))),
                  field_norm_definition='Combined E,cB Euclidean norm over all stored raw field samples; no near-field subtraction.',
                  quintic_difference_expected=True,
                  physical_acceptance=False,
                  scope='Short complete program regression; does not establish the production mJ anomaly is removed.')
    save(root / 'COMPARISON.json', result)
    return result


if __name__ == '__main__':
    import sys
    print(json.dumps(compare(Path(sys.argv[1])), indent=2))
