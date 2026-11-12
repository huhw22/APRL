"""Reject phase-clone columns before using SPARC data as a random 3D PIC beam.

This detects exact repeated transverse positions, including shuffled copies.
It is a necessary loading check, not a statistical-independence certificate.
Historical frozen diagnostic scripts are deliberately not changed.
"""
import argparse
import json
from pathlib import Path

import h5py
import numpy as np


def inspect_cell(xy, state, longitudinal):
    xy = np.asarray(xy)
    state = np.asarray(state)
    longitudinal = np.asarray(longitudinal)
    if (xy.ndim != 2 or xy.shape[1] != 2 or state.ndim != 2
            or len(xy) != len(state) or longitudinal.shape != (len(xy),)
            or not len(xy) or not all(np.isfinite(v).all()
                                      for v in (xy, state, longitudinal))):
        raise ValueError('Invalid/nonfinite SPARC loading arrays')
    _, inverse, counts = np.unique(xy, axis=0, return_inverse=True, return_counts=True)
    low = np.full(len(counts), np.inf)
    high = np.full(len(counts), -np.inf)
    np.minimum.at(low, inverse, longitudinal)
    np.maximum.at(high, inverse, longitudinal)
    column = (counts > 1) & (high > low)
    return dict(records=len(xy), unique_xy=len(counts),
                unique_transverse_energy_states=len(np.unique(state, axis=0)),
                column_particles=int(counts[column].sum()),
                column_groups=int(column.sum()))


def validate_loading(path, kind='crossings', particles_per_cell=1024, report=None):
    if particles_per_cell <= 0:
        raise ValueError('particles_per_cell must be positive')
    result = dict(input=str(Path(path).resolve()), kind=kind, checked_cells=0,
                  checked_records=0, column_particles=0, column_groups=0,
                  min_unique_xy=None, min_unique_transverse_energy_states=None,
                  first_bad_cell=None, passed=False,
                  scope='Exact within-cell transverse-column rejection only; '
                        'does not certify independence, physical shot noise or FEL accuracy.')
    with h5py.File(path, 'r') as source:
        if kind == 'native':
            count = int(source['slicecount'][0])
            def cells():
                for i in range(count):
                    group = source[f'slice{i+1:06d}']
                    xy = np.column_stack([group[k][:] for k in ('x', 'y')])
                    state = np.column_stack([group[k][:] for k in ('x', 'y', 'px', 'py', 'gamma')])
                    yield xy, state, group['theta'][:]
        elif kind == 'crossings':
            records = source['particles/records']
            if not len(records) or len(records) % particles_per_cell:
                raise ValueError('Incomplete SPARC input cell')
            def cells():
                for start in range(0, len(records), particles_per_cell):
                    a = records[start:start+particles_per_cell]
                    xy = a['plane_position_m']
                    yield xy, np.column_stack([xy, a['proper_velocity']]), a['arrival_time_offset_s']
        else:
            raise ValueError('Unknown loading format')
        for i, (xy, state, longitudinal) in enumerate(cells()):
            cell = inspect_cell(xy, state, longitudinal)
            result['checked_cells'] += 1
            result['checked_records'] += cell['records']
            for key in ('column_particles', 'column_groups'):
                result[key] += cell[key]
            for key in ('unique_xy', 'unique_transverse_energy_states'):
                name = 'min_' + key
                result[name] = cell[key] if result[name] is None else min(result[name], cell[key])
            if cell['column_particles'] and result['first_bad_cell'] is None:
                result['first_bad_cell'] = i
    if not result['checked_records']:
        raise ValueError('Empty SPARC input')
    result['passed'] = result['column_particles'] == 0
    if report:
        with Path(report).open('x') as stream:
            stream.write(json.dumps(result, indent=2) + '\n')
    if not result['passed']:
        raise RuntimeError(
            'Rejected SPARC phase-clone input: '
            f"{result['column_particles']}/{result['checked_records']} particles "
            'share exact (x,y) at different longitudinal positions/phases. '
            'Use the independent 6D physical_electron_source workflow and its '
            'explicit macro-particle sampling. Low-order BF and same-input '
            'GENESIS checks do not authorize this loading as a random 3D PIC beam.')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('kind', choices=['native', 'crossings'])
    parser.add_argument('input', type=Path)
    parser.add_argument('--report', type=Path)
    parser.add_argument('--particles-per-cell', type=int, default=1024)
    args = parser.parse_args()
    validate_loading(args.input, args.kind, args.particles_per_cell, args.report)
