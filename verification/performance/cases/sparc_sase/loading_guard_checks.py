"""Small offline rejection controls; does not launch GENESIS or PIC."""
import tempfile
from pathlib import Path

import h5py
import numpy as np

from loading_guard import inspect_cell, validate_loading


def main():
    rng = np.random.default_rng(20261112)
    base = rng.normal(size=(64, 5))
    state = np.repeat(base, 16, axis=0)
    phase = np.tile(np.arange(16, dtype=float), 64)
    original = inspect_cell(state[:, :2], state, phase)
    assert original['column_particles'] == 1024 and original['column_groups'] == 64
    order = rng.permutation(1024)
    assert inspect_cell(state[order, :2], state[order], phase[order]) == original
    # Matching global b1 says nothing about the repeated transverse columns.
    assert abs(np.exp(2j*np.pi*phase/16).mean()) < 1e-15
    independent = rng.normal(size=(1024, 5))
    assert inspect_cell(independent[:, :2], independent, phase)['column_particles'] == 0
    dtype = [('plane_position_m', '<f8', (2,)),
             ('proper_velocity', '<f8', (3,)), ('arrival_time_offset_s', '<f8')]
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        for name, values, reject in [('copies', state, True), ('independent', independent, False)]:
            path = tmp/(name+'.h5')
            records = np.zeros(1024, dtype=dtype)
            records['plane_position_m'] = values[:, :2]
            records['proper_velocity'] = values[:, 2:]
            records['arrival_time_offset_s'] = phase
            with h5py.File(path, 'x') as f:
                f.create_dataset('particles/records', data=records)
            try:
                audit = validate_loading(path, report=tmp/(name+'.json'))
            except RuntimeError:
                assert reject and (tmp/(name+'.json')).exists()
            else:
                assert not reject and audit['checked_records'] == 1024
        native = tmp/'native.h5'
        with h5py.File(native, 'x') as f:
            f['slicecount'] = [1]
            group = f.create_group('slice000001')
            for j, key in enumerate(('x', 'y', 'px', 'py', 'gamma')):
                group[key] = state[:, j]
            group['theta'] = phase
        try:
            validate_loading(native, 'native')
        except RuntimeError:
            pass
        else:
            raise AssertionError('Native cloned dump was accepted')
    print('PASS: cloned, shuffled-clone, zero-b1, independent and native/reuse guards')


if __name__ == '__main__':
    main()
