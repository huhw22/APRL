"""Compare a closed first-segment GENESIS reference with the two PIC caches.

The SPARC GENESIS optical-field reducer retains its transformations, windows and
normalization; its dump-coordinate lookup also accepts a terminal boundary. A
names-only adapter accepts another prefix without rewriting its HDF5 files.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import sys

import h5py
import numpy as np


def require(value, message):
    if not value:
        raise RuntimeError(message)


def save(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False,
        default=lambda x: x.item() if isinstance(x, np.generic) else str(x))+'\n')


def digest(path):
    value = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(4*1024*1024), b''):
            value.update(block)
    return value.hexdigest()


def identity(path):
    path = Path(path)
    stat = path.stat()
    return dict(path=str(path.resolve()), bytes=stat.st_size, mtime_ns=stat.st_mtime_ns)


def run(workspace, reference, output, closure_record, reference_root='reference',
        input_validation=None, spectra_root=None, reference_time_shift_s=0.):
    workspace, reference, output, closure_record = map(
        lambda p: Path(p).resolve(), [workspace, reference, output, closure_record])
    # Both markers are required before any GENESIS HDF5 file is opened.
    require((reference/'TRACK_COMPLETE').read_text().strip() == '0',
            'GENESIS has not published successful tracking completion')
    closure = json.loads(closure_record.read_text())
    require(closure.get('run', {}).get('exit_code') == 0,
            'GENESIS launcher did not exit normally')
    if input_validation is not None:
        input_validation = Path(input_validation).resolve()
        require(input_validation.is_file(), 'Missing input-validation record')
        input_audit = json.loads(input_validation.read_text())
    else:
        input_audit = None
    campaign = workspace/'results/energy_critical_control_20261105/campaign-r4'
    if spectra_root is None:
        spectra_root = campaign.parent/'power-review/spectra-r2'
    spectra_root = Path(spectra_root).resolve()
    spectral_review = json.loads((spectra_root/'RESULT.json').read_text())
    require(spectral_review.get('complete') is True, 'PIC spectral review not complete')
    summaries = {mode: json.loads((spectra_root/mode/'p1-summary.json').read_text())
                 for mode in ['off', 'on']}
    for mode in ['off', 'on']:
        marker = json.loads((campaign/mode/'PIC_COMPLETE.json').read_text())
        require(marker.get('analysis_complete') is True, 'PIC final closure missing')
        require(summaries[mode]['missing_ticks'] == [], 'Unexpected missing PIC field samples')
    first = summaries['on']
    for key in ['fft_samples', 'dt_s', 'first_retarded_s', 'last_retarded_s', 'z_m']:
        require(summaries['off'][key] == first[key], 'PIC support differs: '+key)
    geometry = json.loads((campaign/'GEOMETRY.json').read_text())
    prior_plan_path = workspace/'results/sparc_sase_comparison_20261103/compact-production-r1/plan.json'
    plan = json.loads(prior_plan_path.read_text())
    # The fresh run retains the original GENESIS slice grid and padding.
    # Only the PIC public clock changes when selecting central slices.
    # Positive shift delays the optical packet in the common physical clock.
    # A terminal native integer reindex can be undone here without changing
    # any optical samples, particle states, or PIC observation times.
    require(np.isfinite(reference_time_shift_s), 'Nonfinite reference time shift')
    plan['geometry']['reference_time_offset_s'] = geometry['reference_time_offset_lab_s']-reference_time_shift_s
    local_helper = Path(__file__).resolve().parent/'analyze_explicit.py'
    helper = local_helper if local_helper.exists() else (
        workspace/'source/aprl/verification/analysis/sparc_sase/analyze_explicit.py')
    sys.path.insert(0, str(helper.parent))
    from analyze_explicit import genesis, table, lattice_dump_positions
    source = reference/(reference_root+'.out.h5')
    source_before = identity(source)
    with h5py.File(source, 'r') as handle:
        z = handle['Lattice/z'][:]
        dz = handle['Lattice/dz'][:]
        dumps = handle['Meta/Fielddumps/intstep'][:].astype(int)
    require(len(dumps) > 0, 'No GENESIS field dump')
    positions = lattice_dump_positions(z,dz,dumps)
    dump_index = int(np.argmin(abs(positions-first['z_m'])))
    step = int(dumps[dump_index])
    require(abs(positions[dump_index]-first['z_m']) < 1e-8, 'No dump at current first-segment exit')
    dump = reference/(reference_root+'.'+str(step)+'.fld.h5')
    dump_before = identity(dump)
    output.mkdir(exist_ok=False)
    adapter = output/'reference-adapter'
    adapter.mkdir()
    (adapter/'reference.out.h5').symlink_to(source)
    (adapter/('reference.'+str(step)+'.fld.h5')).symlink_to(dump)
    # Only cache inputs are linked. The established reducer reads them but
    # writes all new comparison products in this independent output folder.
    (output/'p1-spectra.npz').hardlink_to(spectra_root/'on/p1-spectra.npz')
    shutil.copy2(spectra_root/'on/p1-summary.json', output/'p1-summary.json')
    cache_identity = {str(spectra_root/mode/name): digest(spectra_root/mode/name)
        for mode in ['off','on'] for name in ['p1-summary.json','p1-spectra.npz']}
    reference_rows = genesis(campaign, output, plan, [first], reference_directory=adapter)
    require(len(reference_rows) == 1, 'Unexpected number of reduced reference planes')
    gen = reference_rows[0]
    duration = first['fft_samples']*first['dt_s']
    rows = []
    for window in ['rect','hann']:
        gen_energy = gen['results'][window]['band_J']
        require(gen_energy > 0, 'Nonpositive GENESIS target-band energy')
        for mode in ['off','on']:
            stats = summaries[mode]['results'][window]
            rows.append(dict(window=window, case=mode,
                native_band_energy_J=stats['native_Poynting_band_J'],
                native_record_mean_band_power_W=stats['native_Poynting_band_J']/duration,
                transverse_E_proxy_record_mean_band_power_W=stats['E_transverse_band_proxy_J']/duration,
                GENESIS_record_mean_band_power_W=gen_energy/duration,
                native_over_GENESIS=stats['native_Poynting_band_J']/gen_energy,
                E_proxy_over_GENESIS=stats['E_transverse_band_proxy_J']/gen_energy))
    table(output/'genesis_power_comparison.csv', rows)
    save(output/'PIC_SUMMARIES.json', summaries)
    save(output/'ADAPTED_REFERENCE_PLAN.json', plan)
    save(output/'INPUT_VALIDATION_COPY.json', input_audit)
    require(identity(source) == source_before and identity(dump) == dump_before,
            'Closed GENESIS source changed while reducing')
    require(all(digest(path) == value for path,value in cache_identity.items()),
            'Original PIC spectral cache changed')
    result = dict(complete=True, physical_acceptance=False, new_pic_runs=0,
        target_photon_eV=spectral_review['target_photon_eV'],
        band_eV=spectral_review['band_eV'], record_duration_s=duration,
        aperture_side_m=spectral_review['aperture_side_m'],
        rows=rows, genesis=gen, field_dump_step=step,
        reference_output=source_before, reference_dump=dump_before,
        reference_run_closure=dict(path=str(closure_record), sha256=digest(closure_record)),
        input_validation_record=None if input_validation is None else dict(
            path=str(input_validation), sha256=digest(input_validation)),
        reducer=dict(path=str(helper), sha256=digest(helper),
            change='Terminal-boundary dump-coordinate lookup only; existing interior coordinates, transforms, reference_lag, aperture, windows and normalization preserved'),
        PIC_cache_identities=cache_identity,
        current_PIC_time_offset_s=geometry['reference_time_offset_lab_s'],
        reference_time_shift_s=reference_time_shift_s,
        adapted_GENESIS_time_offset_s=plan['geometry']['reference_time_offset_s'],
        retained_original_GENESIS_grid_bunch_length_m=plan['physics']['bunch_length_m'],
        cautions=['GENESIS is a forward near-resonance radiation-envelope model; PIC data remain raw total E/B without near-field subtraction.',
                  'Native PIC B retains its original stored Yee half-time layer.',
                  'All band powers use the same PIC N*dt denominator and matched physical aperture/absolute-time window.',
                  'Input consistency is documented by the supplied launch validation; reduction alone cannot certify a source model.',
                  'One physical-noise realization is not a 5% SASE statistical acceptance.'])
    save(output/'RESULT.json', result)
    plot(output, spectra_root, result)
    print(json.dumps(rows, indent=2), flush=True)


def plot(output, spectra_root, result):
    os.environ.setdefault('MPLCONFIGDIR','/tmp/aprl-genesis-estimate-mpl')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False})
    arrays = {mode: np.load(spectra_root/mode/'p1-spectra.npz') for mode in ['off','on']}
    gen = np.load(output/'p1-genesis.npz')
    duration = result['record_duration_s']
    photon = result['target_photon_eV']
    styles = {'off':('#24648c','PIC unfiltered'), 'on':('#b18428','PIC seven-point')}
    fig, axes = plt.subplots(2,2,figsize=(12.2,8),layout='constrained')
    for j, window in enumerate(['rect','hann']):
        ax = axes[0,j]
        for mode,(color,label) in styles.items():
            energy = arrays[mode]['energy_eV']
            density = arrays[mode][window+'_native_flux_J_per_eV']/duration
            mask = (energy >= .9*photon)&(energy <= 1.1*photon)
            ax.semilogy(energy[mask],np.where(density[mask]>0,density[mask],np.nan),
                        color=color,label=label+' native E×B')
        eg = gen['energy_eV']
        dg = gen[window+'_J_per_eV']/duration
        maskg = (eg >= .9*photon)&(eg <= 1.1*photon)
        ax.semilogy(eg[maskg],np.where(dg[maskg]>0,dg[maskg],np.nan),
                    color='#668d6b',label='GENESIS radiation')
        ax.axvspan(.95*photon,1.05*photon,color='#777777',alpha=.12)
        ax.set(title='Rectangular window' if window=='rect' else 'Energy-normalized Hann',
               xlabel='Photon energy (eV)',ylabel='Record-mean band density (W/eV)')
        ax.legend(fontsize=8);ax.grid(alpha=.2)
    x = np.arange(2)
    for k,(case,color,label) in enumerate([('off','#24648c','PIC unfiltered'),
            ('on','#b18428','PIC seven-point'),('genesis','#668d6b','GENESIS')]):
        values=[]
        for window in ['rect','hann']:
            row=next(q for q in result['rows'] if q['window']==window and q['case']=='on')
            if case=='genesis':value=row['GENESIS_record_mean_band_power_W']
            else:value=next(q['native_record_mean_band_power_W'] for q in result['rows']
                            if q['window']==window and q['case']==case)
            values.append(value)
        bars=axes[1,0].bar(x+(k-1)*.25,values,width=.23,color=color,label=label)
        axes[1,0].bar_label(bars,labels=[f'{v:.3g}' for v in values],padding=3,fontsize=8)
    axes[1,0].set(yscale='log',xticks=x,xticklabels=['Rectangular','Normalized Hann'],
                  ylabel='Same-record band power (W)',title='Integrated resonance band (±5%)')
    for k,(case,(color,label)) in enumerate(styles.items()):
        values=[next(q['native_over_GENESIS'] for q in result['rows']
                     if q['window']==window and q['case']==case) for window in ['rect','hann']]
        bars=axes[1,1].bar(x+(k-.5)*.3,values,width=.28,color=color,label=label)
        axes[1,1].bar_label(bars,labels=[f'{v:.3g}' for v in values],padding=3,fontsize=8)
    axes[1,1].axhline(1,color='#777777',lw=1)
    axes[1,1].set(yscale='log',xticks=x,xticklabels=['Rectangular','Normalized Hann'],
                  ylabel='PIC band / GENESIS',title='Comparison to the radiation reference')
    for ax in axes[1]:ax.legend(fontsize=8);ax.grid(alpha=.2)
    fig.suptitle('First-segment power estimate — common aperture, band and record\nPIC raw total fields; GENESIS radiation envelope; no PIC near-field subtraction',fontsize=12)
    fig.savefig(output/'genesis_power_comparison.png',dpi=200)
    fig.savefig(output/'genesis_power_comparison.pdf')
    plt.close(fig)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('workspace',type=Path)
    parser.add_argument('reference',type=Path)
    parser.add_argument('output',type=Path)
    parser.add_argument('--closure-record',type=Path,required=True)
    parser.add_argument('--reference-root',default='reference')
    parser.add_argument('--input-validation',type=Path)
    parser.add_argument('--spectra-root',type=Path)
    parser.add_argument('--reference-time-shift-s',type=float,default=0.)
    args=parser.parse_args()
    run(args.workspace,args.reference,args.output,args.closure_record,
        args.reference_root,args.input_validation,args.spectra_root,args.reference_time_shift_s)
