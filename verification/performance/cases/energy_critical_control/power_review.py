"""Closed-only field-power comparison using the unchanged SPARC FFT reducer."""
import hashlib
import json
import os
from pathlib import Path
import re
import sys

import h5py
import numpy as np
from scipy.constants import c, h, e, mu_0

sys.path.insert(0, str(Path(__file__).resolve().parent))
from analyze_explicit import pic_plane, save, table, need


def sha(path):
    value = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(4*1024*1024), b''):
            value.update(block)
    return value.hexdigest()


def run(out):
    out = Path(out).resolve()
    plan = json.loads((out/'plan.json').read_text())
    frozen = json.loads((out/'FROZEN.json').read_text())
    need(all(sha(p) == digest for p,digest in frozen.items()), 'Frozen analysis identity changed')
    root = Path(plan['campaign'])
    geometry = json.loads((root/'GEOMETRY.json').read_text())
    prior = json.loads(Path(plan['prior_plan']).read_text())
    prior['geometry'].update(dt_s=geometry['dt_box_s'], dx_m=geometry['dx_m'],
        physical_half_width_m=(geometry['nx']/2-10)*geometry['dx_m'])
    prior['physics']['boost_gamma'] = geometry['boost_gamma']
    photon = h*c/e/prior['physics']['wavelength_m']
    constants_text = Path(plan['compiled_constants_source']).read_text()
    mu0_match = re.search(r'const Double mu0\s*=\s*([\d.eE+-]+)', constants_text)
    need(mu0_match is not None, 'Cannot identify compiled SI permeability')
    compiled_mu0 = float(mu0_match.group(1))
    results, rows, source_identity = {}, [], {}
    times = []
    for mode in ['off','on']:
        marker = json.loads((root/mode/'PIC_COMPLETE.json').read_text())
        need(marker.get('analysis_complete') is True and marker.get('stage') == 'full', 'Case not closed')
        raw = root/mode/'detectors/p1.h5'
        stat = raw.stat()
        source_identity[mode] = dict(path=str(raw), size=stat.st_size, mtime_ns=stat.st_mtime_ns)
        with h5py.File(raw,'r') as f:
            need(int(f['field_plane/complete'][()]) == 1, 'Raw plane not closed')
            times.append(f['field_plane/time_s'][:])
        case_out = out/mode
        case_out.mkdir(exist_ok=False)
        data = pic_plane(root, case_out, prior, 'p1', case=mode)
        need(not data['missing_ticks'], 'Unexpected missing field sample')
        duration = data['fft_samples']*data['dt_s']
        spectrum = np.load(case_out/'p1-spectra.npz')
        energy = spectrum['energy_eV']
        de = energy[1]-energy[0]
        for window in ['rect','hann']:
            q = data['results'][window]
            flux = spectrum[window+'_native_flux_J_per_eV']
            below = float(flux[energy<.5*photon].sum()*de)
            row = dict(mode=mode,window=window,record_duration_s=duration,
                native_full_frequency_energy_J=q['raw_full_frequency_Poynting_J'],
                native_full_frequency_record_mean_power_W=q['raw_full_frequency_Poynting_J']/duration,
                native_below_half_resonance_energy_J=below,
                native_resonant_band_energy_J=q['native_Poynting_band_J'],
                native_resonant_band_record_mean_power_W=q['native_Poynting_band_J']/duration,
                transverse_E_proxy_resonant_band_energy_J=q['E_transverse_band_proxy_J'],
                transverse_E_proxy_resonant_band_record_mean_power_W=q['E_transverse_band_proxy_J']/duration,
                parseval_flux_relative=q['parseval_flux_relative'],parseval_E_relative=q['parseval_E_relative'])
            rows.append(row)
        # Independent online plane integral uses a trapezoid, whereas the
        # established FFT convention uses all N sample rectangles.
        check = np.trapz(spectrum['native_total_power_W'], dx=data['dt_s'])
        # Preserve the earlier FFT convention, which imports the installed
        # SciPy CODATA constants. The compiled PIC uses an older explicit
        # mu0, so convert only this independent diagnostic comparison.
        check_pic_constants = check*mu_0/compiled_mu0
        ledger = json.loads(Path(plan['ledger_review']).read_text())
        expected = ledger['cases'][mode]['p1']['physical_field_energy_lab_J']
        need(abs(check_pic_constants-expected) <= 1e-10*max(abs(check_pic_constants),abs(expected)), 'FFT/raw power disagrees with online energy plane after SI constant conversion')
        data['raw_power_trapezoid_J'] = float(check)
        data['raw_power_trapezoid_with_compiled_mu0_J'] = float(check_pic_constants)
        data['online_energy_comparison_relative'] = float((check_pic_constants-expected)/expected)
        results[mode] = data
        after = raw.stat()
        need((after.st_size,after.st_mtime_ns)==(stat.st_size,stat.st_mtime_ns),'Closed raw source changed')
    need(np.array_equal(*times), 'Paired raw times differ')
    table(out/'power_comparison.csv',rows)
    save(out/'RESULT.json',dict(complete=True,physical_acceptance=False,
        target_photon_eV=photon,band_eV=[.95*photon,1.05*photon],
        aperture_side_m=2*prior['geometry']['physical_half_width_m'],
        aperture_area_m2=(2*prior['geometry']['physical_half_width_m'])**2,
        unchanged_prior_fft_reducer=True,rows=rows,cases=results,sources=source_identity,
        scipy_mu0=mu_0,compiled_mu0=compiled_mu0,
        cautions=['Raw total lab fields; no near-field subtraction.',
                  'Original stored Yee half-time magnetic field; no phase/time correction.',
                  'Signed Poynting spectrum is not a uniquely separated radiation spectrum.',
                  'Windowed values are energy-normalized finite-record diagnostics; powers use N*dt.',
                  'A source shape response is not a prediction of final radiated power.']))
    need(all(sha(p)==digest for p,digest in frozen.items()),'Frozen analysis identity changed')
    print(json.dumps(rows,indent=2),flush=True)


if __name__ == '__main__':
    run(sys.argv[1])
