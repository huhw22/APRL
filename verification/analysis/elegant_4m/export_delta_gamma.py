#!/usr/bin/env python3
"""Export ID-matched delta-gamma moments from an existing full particle scan.

Usage: export_delta_gamma.py --summary particle-summary.json
       --inventory INVENTORY.json --output-dir NEW_RESULT_DIR

No HDF5 data are rescanned. Positive delta gamma means entry minus downstream.
The inventory must refer to unchanged files and a complete, equal-weight cohort.
"""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--summary', type=Path, required=True)
    parser.add_argument('--inventory', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    summary = json.loads(args.summary.read_text())
    inventory = json.loads(args.inventory.read_text())
    entry = next(p for p in summary['planes'] if p['label'] == 'entry')
    assert len(summary['planes']) == inventory['usable_particle_crossing_file_count']
    indexed = {p['source_file']: p for p in inventory['planes']}
    rows = []
    for plane in summary['planes']:
        meta = indexed[plane['source_file']]
        quality = plane['quality']
        for key in ('complete', 'commit_equals_shape', 'unique_ids_1_through_N',
                    'all_downstream', 'finite_all_required_numeric_fields',
                    'same_ID_cohort_as_entry',
                    'same_source_ids_and_represented_weights_as_entry'):
            assert quality[key] is True, (plane['label'], key)
        assert quality['raw_relative_weight_min'] == quality['raw_relative_weight_max'] == 1
        assert quality['charge_C_min'] == quality['charge_C_max'] == -1.602176634e-19
        assert meta['complete'] == 1
        assert plane['records'] == entry['records'] == meta['records'] == meta['committed_records']
        assert plane['z_m'] == meta['plane_z_m']
        assert meta['size_and_mtime_equal_copy_manifest']
        stat = Path(plane['source_file']).stat()
        assert (stat.st_size, stat.st_mtime_ns) == (meta['size_bytes'], meta['mtime_ns'])
        rest_keV = plane['electron_rest_mass_kg'] * 299792458**2 / 1.602176634e-19 / 1000
        mean = plane['mean_ID_matched_loss_keV'] / rest_keV
        std = plane['rms_ID_matched_loss_about_mean_keV'] / rest_keV
        independent_mean = entry['mean_gamma'] - plane['mean_gamma']
        # Comparing two rounded absolute means requires their own ULP scale.
        bound = 8 * (math.ulp(entry['mean_gamma']) + math.ulp(plane['mean_gamma']))
        assert abs(mean - independent_mean) <= bound
        assert math.isfinite(mean) and math.isfinite(std) and std >= 0
        ratio = std / mean if mean != 0 else None
        rows.append(dict(
            label=plane['label'], z_m=plane['z_m'], particle_count=plane['records'],
            avg_delta_gamma=mean, std_delta_gamma=std, std_over_avg=ratio,
            avg_final_minus_entry_gamma=-mean,
            std_over_avg_final_minus_entry=-ratio if ratio is not None else None,
            avg_energy_loss_keV=plane['mean_ID_matched_loss_keV'],
            std_energy_loss_keV=plane['rms_ID_matched_loss_about_mean_keV'],
            mean_gamma_difference_check=independent_mean,
            mean_check_absolute_difference=abs(mean - independent_mean),
            mean_check_rounding_bound=bound,
            electron_rest_energy_keV=rest_keV, source_file=plane['source_file']))
    result = dict(
        simulation_job_id=inventory['simulation_job_id'],
        definition='delta_gamma_i(z) = gamma_i(entry) - gamma_i(z); positive is loss',
        frame='laboratory', reference_plane_z_m=entry['z_m'],
        std_definition='sqrt(mean((delta_gamma - mean(delta_gamma))**2)); ddof=0',
        ratio_definition='std(delta_gamma)/avg(delta_gamma); undefined at zero average',
        selection='All ID-matched records, no slice/chirp/aperture selection; unit weights',
        provenance={str(p.resolve()): digest(p) for p in (args.summary, args.inventory, Path(__file__))},
        reuse='Prior full-scan ID-matched moments, current metadata/stat check; no new full HDF5 scan',
        limitations=[
            'Reference is the -0.03 m entry detector, not the -0.05 m original input.',
            'Moments describe stored records, not isolated radiation losses.',
            'Historical full-step boost-momentum interpolation introduces a negative lab-energy bias, '
            'notably at 1 m and 3 m; all of the observed differences are not proven to have that origin.',
            'Recent field-only runs 18932/18938/18992/18993 have no particle-plane outputs.',
        ], rows=rows)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / 'delta_gamma.json').write_text(json.dumps(result, indent=2) + '\n')
    with (args.output_dir / 'delta_gamma.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    for row in rows:
        print(row['z_m'], row['avg_delta_gamma'], row['std_delta_gamma'], row['std_over_avg'])


if __name__ == '__main__':
    main()
