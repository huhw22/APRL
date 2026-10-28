#!/usr/bin/env python3
"""Independent long-double checks of saved old/new field spectra."""
import argparse
import importlib.util
import json
from pathlib import Path
import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--include-undulators', action='store_true')
    args = parser.parse_args()
    root = args.root.resolve()
    spec = importlib.util.spec_from_file_location('independent', root/'scripts/independent_review.py')
    checker = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(checker)
    names = ('old_free_drift', 'new_free_drift')
    expected_jobs = [18938, 18993]
    pairs = [('old_free_drift', 'new_free_drift')]
    if args.include_undulators:
        names += ('old_undulator', 'new_undulator')
        expected_jobs += [18932, 18992]
        pairs.append(('old_undulator', 'new_undulator'))
    reviewed = {name: checker.review_case(root/name) for name in names}
    summaries = [json.loads((root/name/'summary.json').read_text()) for name in names]
    by_name = dict(zip(names, summaries))
    for old, new in pairs:
        for key in ('samples', 'time_start_s', 'time_end_s', 'sample_interval_s',
                    'physical_shape_yx', 'physical_area_m2'):
            checker.require(by_name[old][key] == by_name[new][key], f'{old}/{new}: mismatched {key}')
        with np.load(root/old/'spectra.npz') as oa, np.load(root/new/'spectra.npz') as na:
            checker.require(np.array_equal(oa['frequency_Hz'], na['frequency_Hz']),
                            f'{old}/{new}: mismatched frequency grid')
    checker.require([s['simulation_job'] for s in summaries] == expected_jobs, 'Wrong source jobs')
    checker.require(all(s['source_attributes']['frame']=='laboratory' for s in summaries), 'Wrong frame')
    checker.require(all(s['physical_area_m2']==summaries[0]['physical_area_m2'] for s in summaries), 'Different apertures')
    for name, summary in zip(names, summaries):
        checker.require(summary['input_transaction']['complete']==1 and summary['input_unchanged_during_read'],
                        'Incomplete or modified source')
        arrays = np.load(root/name/'spectra.npz')
        # Check half-open band integration independently of the plotting code.
        energy = arrays['photon_energy_eV'].astype(np.longdouble)
        de = checker.H_EV_S/(int(summary['samples'])*np.longdouble(summary['sample_interval_s']))
        bands = {}
        for window in ('rect','hann'):
            flux = arrays[window+'_native_flux_J_per_eV'].astype(np.longdouble)
            bands[window] = {}
            for label, lo, hi in (('all',0,np.inf),('below_1_eV',0,1),
                                  ('90_to_93_eV',90,93),('91p3_to_91p7_eV',91.3,91.7)):
                mask = (energy>=lo)&(energy<hi)
                value = flux[mask].sum()*de
                bands[window][label] = dict(bins=int(mask.sum()), native_energy_J=float(value),
                    record_average_power_W=float(value/(summary['samples']*np.longdouble(summary['sample_interval_s']))))
        reviewed[name]['half_open_bands'] = bands
    record = dict(status='PASS', cases=reviewed, matching_full_record=True,
        matching_record_scope='Within each old/new pair, not between different detector planes',
        matching_pairs=pairs,
        scope='Independent long-double integration of saved spectra against time sums, axis/area/Cauchy checks; no re-FFT or reread of large H5.',
        limitations='Total native stored lab E/B; no proof of pure radiation, simultaneous time layers, quiet-tail coverage or grid convergence.')
    (root/'INDEPENDENT_REVIEW.json').write_text(json.dumps(record, indent=2, allow_nan=False)+'\n')
    print('PASS: saved spectra independently checked; full records and frequency grids match within old/new pairs.')


if __name__=='__main__':
    main()
