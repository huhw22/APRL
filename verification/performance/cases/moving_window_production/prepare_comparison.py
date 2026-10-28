#!/usr/bin/env python3
"""Record old/new full-field comparison scope; no simulation or FFT is run."""
import argparse
import hashlib
import json
from pathlib import Path

import h5py
import yaml


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def require(ok, message):
    if not ok:
        raise ValueError(message)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    args = parser.parse_args()
    root = args.root.resolve()
    plan = json.loads((root / 'plan.json').read_text())
    old_analysis = Path(plan['workspace']) / 'results/field_free_drift_spectra_20261028'
    old_final = json.loads((old_analysis / 'FINAL.json').read_text())
    for name, expected in old_final['sha256'].items():
        require(sha(old_analysis / name) == expected, 'Old analysis changed: ' + name)
    report = {
        'format': 'aprl-moving-full-comparison-plan-v1',
        'status': 'PREPARED_NOT_PHYSICALLY_ACCEPTED',
        'binary_sha256': plan['binary_sha256'],
        'old_analysis_final_sha256': sha(old_analysis / 'FINAL.json'),
        'cases': {},
        'aperture': {'physical_pixels_xy': [140, 140], 'excluded_edge_cells': 10,
                     'pixel_pitch_m': 5e-6, 'area_m2': 4.9e-7},
        'comparisons': [
            'Match each new case to its own old case, retaining all six raw laboratory E/B components.',
            'Check complete/committed/finite output and exact physical sampling coordinates before spectral analysis.',
            'Use absolute laboratory timestamps; no fitted time shift. If lengths differ, report both common-time-span and full-record spectra; no invented quiet tail.',
            'FFT each physical pixel before area integration; preserve DC and signed native E/B cross-spectrum. Rectangular window primary, normalized Hann sensitivity secondary.',
            'Compare 91.3–91.7eV band energy and average power with their explicit Ndt denominator; compare peak position/width and neighboring spectral background separately.',
            'Compare full/low-frequency spectra, time-resolved integrated flux, and transverse peak-band intensity/energy maps within the same aperture.',
            'Report pulse energy, elapsed time and node-hours separately; early100-step speedup is not a full-run prediction.'
        ],
        'interpretation': [
            'Fixed-window versus moving-window is a physical-observable comparison, not a bitwise implementation replay.',
            'Initial reference, longitudinal boundary handling and grid phase differ as documented in CARD_AUDIT; large low-frequency self-field differences alone are not radiation-peak failure.',
            'Do not subtract the no-undulator case and label the result pure radiation; do not remove MPI-like spectral lines or DC before comparison.',
            'Native Lorentz-transformed E/B retain their stored Yee time staggering; no artificial phase alignment.',
            'No new percentage acceptance threshold is invented. Report differences and judge feasibility against the requested source observables before further method changes.',
            'Particle streaming, memory reclamation, new particle planes and snapshots are outside this round.'
        ]}
    for name, old_name, plane, job in (
            ('undulator_4p5m', 'undulator', 4.5, 18932),
            ('free_drift_0p5m', 'free_drift', .5, 18938)):
        current = Path(plan['cases'][name])
        card = yaml.safe_load((current / 'run.yaml').read_text())
        summary_path = old_analysis / old_name / 'summary.json'
        summary = json.loads(summary_path.read_text())
        field = Path(summary['source_file'])
        stat = field.stat()
        require((stat.st_size, stat.st_mtime_ns) ==
                (summary['source_size_bytes'], summary['source_mtime_ns']), 'Old field changed')
        require(card['detectors']['field_planes'][0]['z'] == plane, 'Wrong new plane')
        require(not card['detectors']['particle_planes'], 'Unexpected particle-plane output')
        with h5py.File(field, 'r') as handle:
            group = handle['field_plane']
            count = len(group['fields'])
            require(int(group['complete'][()]) == 1 and
                    int(group['committed_samples'][()]) == count, 'Old field not committed')
            require(group['fields'].shape[1:] == (160, 160), 'Wrong old aperture')
            require(group.attrs['plane_z_m'] == plane, 'Wrong old plane')
            require(count == summary['samples'], 'Old record length changed')
            metadata = {'samples': count, 'time_start_s': float(group['time_s'][0]),
                        'time_end_s': float(group['time_s'][-1]),
                        'nominal_sample_interval_s': summary['sample_interval_s']}
        report['cases'][name] = {
            'old_job': job, 'old_raw_field': str(field),
            'old_field_size_bytes': stat.st_size, 'old_field_mtime_ns': stat.st_mtime_ns,
            'old_summary_sha256': sha(summary_path), 'old_metadata': metadata,
            'new_directory': str(current), 'new_card_sha256': sha(current / 'run.yaml'),
            'new_raw_field': str(current / 'detectors/exit-field.h5'),
            'old_integrity_scope': 'Prior full finite/Parseval audit plus current unchanged stat and committed metadata; no new full-array scan.'}
    target = root / 'COMPARISON_PLAN.json'
    require(not target.exists(), 'Comparison protocol already recorded')
    target.write_text(json.dumps(report, indent=2, ensure_ascii=False, allow_nan=False) + '\n')
    print(target)


if __name__ == '__main__':
    main()
