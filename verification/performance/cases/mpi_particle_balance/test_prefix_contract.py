#!/usr/bin/env python3
"""Offline prefix transaction controls; no MPI, submission, or live HDF5 reads.

The legacy 17 bounded-stop controls exercise the NEW classifier. Tiny temporary
detectors are synthetic transaction fixtures; their ledger is the unchanged
first 20 steps of an already closed small control. These are reader tests, never
new physical results. Run with the physics Python environment.
"""
import copy
import csv
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import types
import unittest

import h5py
import numpy as np
import yaml


HERE = Path(__file__).resolve().parent
CASES = HERE.parent
WORKSPACE = Path('/public/home/wfollowm/APRL_Workspace')
SUPPORT = WORKSPACE / 'results/mpi_particle_balance_20261106/full-prefix-r1/scripts'
CLOSED = WORKSPACE / 'results/mpi_particle_balance_20261106/smoke-r1/drift-newoff'
OBSERVED = WORKSPACE / 'results/pending_admission_20261028/pair-r1/full/profile-baseline'
STEPS = 20
sys.path.insert(0, str(SUPPORT))
sys.path.insert(0, str(HERE))


def load_module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


DURATION = load_module(CASES / 'full_scale_followup/run.py', 'duration_review')
PREFIX = load_module(HERE / 'prefix_closed.py', 'prefix_closed')
FULL = load_module(HERE / 'analyze_closed.py', 'prefix_contract_full_reader')

# Reuse the exact legacy negative-control definitions without importing its
# execution runner or replacing the full-run reader's default behavior.
shim = types.ModuleType('run')
shim.bounded_duration_status = PREFIX.bounded_duration_status
shim.load_module = load_module
previous_run = sys.modules.get('run')
sys.modules['run'] = shim
try:
    LEGACY = load_module(CASES / 'pending_admission_optimization/test_bounded_contract.py',
                         'prefix_contract_legacy_controls')
finally:
    if previous_run is None:
        sys.modules.pop('run', None)
    else:
        sys.modules['run'] = previous_run

DurationContract = LEGACY.DurationContract
OutputTransactions = LEGACY.OutputTransactions


class PrefixTransactions(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # Refuse to read the source fixture until its whole PIC result is closed.
        assert json.loads((CLOSED / 'PIC_COMPLETE.json').read_text()).get('complete') is True
        cls.log = (OBSERVED / 'aprl.log').read_text()
        cls.runtime = json.loads((OBSERVED / 'aprl.log.launch.json').read_text())
        source_card = yaml.safe_load((CLOSED / 'run.yaml').read_text())
        cls.cadence = int(source_card['energy_ledger']['sample_interval_steps'])
        with h5py.File(CLOSED / 'ledger/energy-ledger.h5', 'r') as handle:
            group = handle['energy_ledger']
            assert int(group['complete'][()]) == 1
            assert int(group['committed_records'][()]) == len(group['records'])
            cls.attrs = dict(group.attrs)
            cls.ranks = int(group.attrs['mpi_size'])
            # Read a bounded prefix only; step20 is on the accepted cadence.
            raw = group['records'][:STEPS + 1]
            cls.rows = raw[raw['step'] <= STEPS]
        assert len(cls.rows) and int(cls.rows['step'][-1]) == STEPS
        cls.particles = int(np.max(cls.rows['active_macroparticles'] +
                                   cls.rows['removed_macroparticles']))
        assert cls.particles >= 2
        cls.discrete, cls.coupling = {}, {}
        for rank in range(cls.ranks):
            for folder, prefix, destination in [('ledger', 'discrete-field', cls.discrete),
                                                ('reference', 'source-coupling', cls.coupling)]:
                path = CLOSED / folder / f'{prefix}-rank-{rank:05d}.csv'
                marker = path.with_suffix('.complete').read_text().split()
                assert marker[0] == '1' and int(marker[1]) >= STEPS
                with path.open() as stream:
                    lines = [next(stream) for _ in range(STEPS + 1)]
                assert list(csv.DictReader(lines))[-1]['step'] == str(STEPS)
                destination[path.name] = ''.join(lines)

    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='aprl-prefix-contract-')
        self.root = Path(self.temp.name)
        for name in ['detectors', 'ledger', 'reference', 'trajectory']:
            (self.root / name).mkdir()
        self.case = dict(directory=str(self.root), ranks=self.ranks, expected_steps=STEPS,
                         bounded_duration_prefix=True, mode='off', balance_enabled=False,
                         expected_particles=self.particles, require_all_particles_through_planes=False,
                         input_source_ids_contiguous=True, expected_source_id_first=1,
                         expected_source_id_last=self.particles)
        self.write_json('run-plan.json', self.case)
        self.write_json('aprl.log.launch.json', self.runtime)
        (self.root / 'aprl.log').write_text(self.log)
        card = {'detectors': {'field_planes': [{'name': 'prefix-field'}],
                             'particle_planes': [{'name': 'prefix-particles'}],
                             'energy_planes': [], 'substep_particle_planes': True},
                'energy_ledger': {'filename': 'energy-ledger.h5',
                                  'sample_interval_steps': self.cadence}}
        (self.root / 'run.yaml').write_text(yaml.safe_dump(card))
        with h5py.File(self.root / 'ledger/energy-ledger.h5', 'w') as handle:
            group = handle.create_group('energy_ledger')
            for key, value in self.attrs.items():
                group.attrs[key] = value
            group['complete'] = np.uint8(0)
            group['committed_records'] = np.uint64(len(self.rows))
            group['records'] = self.rows
        for name, text in self.discrete.items():
            path = self.root / 'ledger' / name
            path.write_text(text)
            path.with_suffix('.complete').write_text(f'0 {STEPS}\n')
        for name, text in self.coupling.items():
            path = self.root / 'reference' / name
            path.write_text(text)
            path.with_suffix('.complete').write_text(f'1 {STEPS}\n')
        with h5py.File(self.root / 'detectors/prefix-field.h5', 'w') as handle:
            group = handle.create_group('field_plane')
            group['complete'] = np.uint8(0)
            group['committed_samples'] = np.uint64(2)
            group['time_s'] = np.array([0., 1.])
            group['fields'] = np.zeros((2, 1, 1))
        dtype = [('particle_id', '<u8'), ('source_id', '<u8'), ('weight', '<f8'), ('time_s', '<f8')]
        records = np.array([(1, 1, 1., 0.), (self.particles, self.particles, 1., 1.)], dtype=dtype)
        for name, groupname in [('prefix-particles.h5', 'particle_plane'),
                                ('prefix-particles-substeps.h5', 'particle_substep_plane')]:
            with h5py.File(self.root / 'detectors' / name, 'w') as handle:
                group = handle.create_group(groupname)
                group['complete'] = np.uint8(0)
                group['committed_records'] = np.uint64(len(records))
                group['records'] = records

    def tearDown(self):
        self.temp.cleanup()

    def write_json(self, relative, value):
        (self.root / relative).write_text(json.dumps(value) + '\n')

    def rejected(self, callback):
        with self.assertRaises((RuntimeError, ValueError, AssertionError, KeyError, OSError)):
            callback()

    def test_closed_prefix_complete0_accepted_without_physical_promotion(self):
        result = PREFIX.verify_prefix_case(self.case)
        self.assertTrue(result['prefix_complete'])
        self.assertFalse(result['physical_complete'])
        self.assertFalse(result['complete_source_through_all_planes'])
        self.assertFalse((self.root / 'PIC_COMPLETE.json').exists())

    def test_prefix_ledger_physical_complete1_rejected(self):
        with h5py.File(self.root / 'ledger/energy-ledger.h5', 'r+') as handle:
            handle['energy_ledger/complete'][()] = 1
        self.rejected(lambda: PREFIX.load_prefix_ledger(self.root, self.case))

    def test_prefix_detector_physical_complete1_rejected(self):
        with h5py.File(self.root / 'detectors/prefix-field.h5', 'r+') as handle:
            handle['field_plane/complete'][()] = 1
        self.rejected(lambda: PREFIX.verify_prefix_case(self.case))

    def test_truncated_discrete_final_step_rejected(self):
        path = self.root / 'ledger/discrete-field-rank-00000.csv'
        path.write_text('\n'.join(path.read_text().splitlines()[:-1]) + '\n')
        self.rejected(lambda: PREFIX.load_prefix_ledger(self.root, self.case))

    def test_discrete_marker_wrong_flag_or_step_rejected(self):
        path = self.root / 'ledger/discrete-field-rank-00000.complete'
        for text in [f'1 {STEPS}\n', f'0 {STEPS - 1}\n']:
            with self.subTest(marker=text):
                path.write_text(text)
                self.rejected(lambda: PREFIX.load_prefix_ledger(self.root, self.case))

    def test_source_marker_wrong_flag_or_step_rejected(self):
        path = self.root / 'reference/source-coupling-rank-00000.complete'
        for text in [f'0 {STEPS}\n', f'1 {STEPS - 1}\n']:
            with self.subTest(marker=text):
                path.write_text(text)
                self.rejected(lambda: PREFIX.verify_prefix_case(self.case))

    def test_missing_close_proof_rejected_before_hdf5(self):
        text = '\n'.join(line for line in self.log.splitlines()
                         if 'phase=detector-close-complete' not in line)
        (self.root / 'aprl.log').write_text(text)
        self.rejected(lambda: PREFIX.verify_prefix_case(self.case))

    def test_other_failure_and_timeout_cannot_be_expected_prefix(self):
        for change in ['other-failure', 'timeout']:
            with self.subTest(change=change):
                (self.root / 'aprl.log').write_text(self.log +
                    ('\nSimulation failure on MPI rank 3: invalid particle owner\n'
                     if change == 'other-failure' else ''))
                runtime = copy.deepcopy(self.runtime)
                if change == 'timeout':
                    runtime['seconds'] = runtime['limit_s']
                self.rejected(lambda: PREFIX.duration_status(self.root, self.case, runtime))

    def test_fake_pic_completion_rejected(self):
        self.write_json('PIC_COMPLETE.json', {'complete': True})
        self.rejected(lambda: PREFIX.verify_prefix_case(self.case))

    def test_default_full_reader_still_rejects_missing_pic_and_complete0(self):
        case = dict(self.case)
        case.pop('bounded_duration_prefix')
        self.rejected(lambda: FULL.verify_case(case))
        self.write_json('PIC_COMPLETE.json', {'complete': True})
        self.rejected(lambda: FULL.verify_case(case))
        self.rejected(lambda: FULL.load_ledger(self.root))

    def test_numeric_prefix_route_needs_marker_and_explicit_plan(self):
        self.rejected(lambda: FULL.load_ledger(self.root))
        self.write_json('PREFIX_COMPLETE.json', {'prefix_complete': True, 'physical_complete': False})
        self.assertEqual(FULL.load_ledger(self.root)[-1], STEPS)
        case = dict(self.case)
        case.pop('bounded_duration_prefix')
        self.write_json('run-plan.json', case)
        self.rejected(lambda: FULL.load_ledger(self.root))

    def test_numeric_prefix_route_rejects_physical_complete_marker(self):
        self.write_json('PREFIX_COMPLETE.json', {'prefix_complete': True, 'physical_complete': True})
        self.rejected(lambda: FULL.load_ledger(self.root))

    def test_full_ledger_loader_keeps_its_complete1_contract(self):
        self.rejected(lambda: PREFIX.report_discrete.load(self.root / 'ledger',
                                                         self.root / 'ledger/energy-ledger.h5'))
        with h5py.File(self.root / 'ledger/energy-ledger.h5', 'r+') as handle:
            handle['energy_ledger/complete'][()] = 1
        for path in (self.root / 'ledger').glob('discrete-field-rank-*.complete'):
            path.write_text(f'1 {STEPS}\n')
        self.assertEqual(PREFIX.report_discrete.load(self.root / 'ledger',
                         self.root / 'ledger/energy-ledger.h5')[-1], STEPS)

    def test_uncommitted_or_nonfinite_detector_rejected(self):
        path = self.root / 'detectors/prefix-field.h5'
        with h5py.File(path, 'r+') as handle:
            handle['field_plane/committed_samples'][()] = 3
        self.rejected(lambda: PREFIX.verify_prefix_case(self.case))
        with h5py.File(path, 'r+') as handle:
            handle['field_plane/committed_samples'][()] = 2
            handle['field_plane/fields'][0, 0, 0] = np.nan
        self.rejected(lambda: PREFIX.verify_prefix_case(self.case))

    def test_missing_detector_and_duplicate_or_unknown_id_rejected(self):
        path = self.root / 'detectors/prefix-particles.h5'
        with h5py.File(path, 'r+') as handle:
            records = handle['particle_plane/records'][:]
            original = records.copy()
            records['source_id'][1] = records['source_id'][0]
            handle['particle_plane/records'][:] = records
        self.rejected(lambda: PREFIX.verify_prefix_case(self.case))
        with h5py.File(path, 'r+') as handle:
            original['source_id'][1] = self.particles + 1
            handle['particle_plane/records'][:] = original
        self.rejected(lambda: PREFIX.verify_prefix_case(self.case))
        path.unlink()
        self.rejected(lambda: PREFIX.verify_prefix_case(self.case))


if __name__ == '__main__':
    unittest.main()
