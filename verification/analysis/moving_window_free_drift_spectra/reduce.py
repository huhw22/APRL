#!/usr/bin/env python3
"""Reduce completed field-only / free-drift runs with the accepted FFT kernel.

Read-only inputs. Run the large scan on a Slurm compute node. This wrapper
checks transaction closure and records the corrected MPI sampling provenance;
it does not change the native stored time layers or subtract a background.
"""
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path

import h5py
import numpy as np


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--kernel', type=Path, required=True)
    p.add_argument('--run', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--plane', type=float, required=True)
    p.add_argument('--job', type=int, required=True)
    p.add_argument('--slab-cells', type=int, default=None)
    p.add_argument('--workers', type=int, default=8)
    args = p.parse_args()
    if not os.environ.get('SLURM_JOB_ID'):
        p.error('Large field scans require a Slurm allocation')
    if args.output.exists():
        p.error('Use a new output directory')
    spec = importlib.util.spec_from_file_location('accepted_spectra', args.kernel)
    kernel = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(kernel)
    field = args.run / 'detectors/exit-field.h5'
    before = field.stat()
    with h5py.File(field, 'r') as f:
        g = f['field_plane']
        n = int(g['committed_samples'][()])
        assert int(g['complete'][()]) == 1
        assert g['fields'].shape == (n, 160, 160)
        assert g['time_s'].shape == (n,)
        assert np.isclose(float(g.attrs['plane_z_m']), args.plane, rtol=0, atol=1e-14)
        sampling = g.attrs['spatial_sampling']
        if isinstance(sampling, bytes): sampling = sampling.decode()
        assert sampling == 'global Yee trilinear stencil; adjacent MPI half-cell planes included'
        transaction = {'complete': 1, 'committed_samples': n,
                       'field_shape': list(g['fields'].shape),
                       'time_shape': list(g['time_s'].shape),
                       'spatial_sampling': sampling}
    log = (args.run / 'aprl.log').read_text()
    assert 'Configured stop reached:' in log
    assert 'phase=detector-close-complete' in log
    config = kernel.yaml.safe_load((args.run / 'run.yaml').read_text())
    moving = bool(config.get('mesh', {}).get('moving_window', {}).get('enabled', False))
    if not moving and args.slab_cells is None:
        p.error('Fixed-window run requires --slab-cells')
    # Legacy slab argument affects only obsolete comb metadata, never the FFT.
    kernel.analyze(args.run, field, args.output, args.slab_cells or 1,
                   args.workers, 256)
    after = field.stat()
    assert (before.st_size, before.st_mtime_ns) == (after.st_size, after.st_mtime_ns)
    summary_path = args.output / 'summary.json'
    summary = json.loads(summary_path.read_text())
    summary['notes'] = [n.replace(
        'Native spectral flux inherits Yee time staggering and the known MPI-interface diagnostic defect; retain its sign.',
        'This run used the corrected MPI detector stencil. Native spectral flux still inherits mixed Yee time layers in the Lorentz-transformed lab E/B; retain its sign and do not label it pure radiation.')
        for n in summary['notes']]
    if moving:
        summary['slab_cells'] = None
        summary['slab_cells_note'] = 'Not applicable to moving-window layout; no fixed-slab comb model used.'
        completed = json.loads((args.run / 'completion-check.json').read_text())
        assert completed['pass'] and completed['output_transactions_complete']
        summary['moving_window_completion_sha256'] = sha(args.run / 'completion-check.json')
    summary['moving_window'] = moving
    summary['sampling_comb_spacing_eV'] = None
    summary['sampling_comb_note'] = 'Legacy comb defect repaired in these runs; no comb removal or notch filtering.'
    summary['simulation_job'] = args.job
    summary['input_transaction'] = transaction
    summary['input_unchanged_during_read'] = True
    summary['source_mtime_ns'] = before.st_mtime_ns
    summary['source_config_sha256'] = sha(args.run / 'run.yaml')
    summary['fft_kernel_sha256'] = sha(args.kernel)
    summary['batch_exit_code'] = int((args.run / 'EXIT_CODE').read_text())
    summary['field_use_scope'] = 'Committed closed field with every scalar checked finite; not certification of the whole simulation batch or a quiet tail.'
    summary['production_physical_stop'] = [x for x in log.splitlines() if x.startswith('Configured stop reached:')]
    summary['production_particle_boundary_reports'] = [x for x in log.splitlines() if x.startswith(('Frozen lab-plane injection:', 'CPML physical-domain exits:'))]
    kernel.dump(summary_path, summary)
    print('PASS: closed full field, fixed MPI sampling, finite scan, both-window Parseval, unchanged input', flush=True)


if __name__ == '__main__':
    main()
