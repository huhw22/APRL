"""Read completed energy-coupled fields; reuse the verified radiation reference."""
import csv
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import yaml

from analyze_explicit import pic_plane, save


def read(p):
    return json.loads(Path(p).read_text())


def run(root, old, out):
    root, old, out = map(Path, [root, old, out])
    out.mkdir(exist_ok=False)
    complete = read(root / 'energy/PIC_COMPLETE.json')
    assert complete['run']['exit_code'] == 0
    assert len(complete['fields']) == 8 and all(x['complete'] for x in complete['fields'])
    assert read(old / 'VALIDATION.json')['reference_phase_error_corrected']
    for p, digest in read(root / 'FROZEN.json').items():
        assert hashlib.sha256(Path(p).read_bytes()).hexdigest() == digest, p
    cards = [yaml.safe_load((root / case / 'run.yaml').read_text()) for case in ['boris', 'energy']]
    def changes(a, b, path=''):
        if isinstance(a, dict) and isinstance(b, dict):
            return sum((changes(a.get(k), b.get(k), path + '/' + str(k)) for k in a.keys() | b.keys()), [])
        return [] if a == b else [path]
    differences = changes(*cards)
    assert set(differences) == {'/mesh/particle_coupling', '/energy_ledger/directory',
      '/detectors/directory', '/mesh/moving_window/reference_directory', '/output/manifest'}, differences
    save(out / 'CARD_DIFFERENCES.json', differences)
    plan = read(root / 'plan.json')
    energy = [pic_plane(root, out, plan, name, case='energy') for name, z in plan['geometry']['field_planes']]
    save(out / 'pic.json', energy)
    explicit, gen = read(old / 'pic.json'), read(old / 'genesis.json')
    rows = []
    for a, b, g in zip(explicit, energy, gen):
        assert a['plane'] == b['plane'] == g['plane']
        for k in ['fft_samples', 'dt_s', 'first_retarded_s', 'last_retarded_s', 'missing_ticks']:
            assert a[k] == b[k], (a['plane'], k)
        name = b['plane'];T = b['fft_samples'] * b['dt_s']
        with np.load(old / (name + '-spectra.npz')) as ca, np.load(out / (name + '-spectra.npz')) as cb:
            for k in ['energy_eV', 'retarded_time_s', 'x_m', 'y_m']:
                assert np.array_equal(ca[k], cb[k]), (name, k)
            x = cb['energy_eV'];dx = x[1] - x[0]
            mask = abs(x / (sum(read(old / 'RADIATION_RESULT.json')['band_eV']) / 2) - 1) <= .05
            for window in ['rect', 'hann']:
                av = a['results'][window]['native_Poynting_band_J']
                bv = b['results'][window]['native_Poynting_band_J']
                gv = g['results'][window]['band_J']
                assert abs(cb[window + '_native_flux_J_per_eV'][mask].sum() * dx / bv - 1) < 1e-12
                l1 = float(abs(cb[window + '_native_flux_J_per_eV'][mask] - ca[window + '_native_flux_J_per_eV'][mask]).sum() * dx / abs(av))
                rows.append(dict(plane=name,z_m=b['z_m'],window=window,record_ps=T*1e12,
                  explicit_W=av/T,energy_coupled_W=bv/T,GENESIS_W=gv/T,
                  energy_over_explicit=bv/av,energy_over_GENESIS=bv/gv,explicit_over_GENESIS=av/gv,
                  energy_change_pct=100*(bv/av-1),band_spectrum_L1_over_explicit=l1,
                  explicit_peak_eV=a['results'][window]['peak_energy_eV'],
                  energy_peak_eV=b['results'][window]['peak_energy_eV'],
                  GENESIS_peak_eV=g['results'][window]['peak_energy_eV']))
    with (out / 'comparison.csv').open('w') as f:
        writer = csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    save(out / 'COMPARISON.json', dict(rows=rows,physical_acceptance=False,
      same_time_grid=True,same_spatial_grid=True,same_reference=True,
      near_field_subtraction=False,reference_directory=str(old),
      scope='Radiation power only; full raw PIC E/B; energy-coupling name is not a full physics certificate'))
    plot(out, old, rows, explicit, energy, gen)
    save(out / 'COMPLETE.json',dict(complete=True,physical_acceptance=False,
      max_parseval_relative=max(p['results'][k]['parseval_flux_relative'] for p in energy for k in ['rect','hann']),
      scope='Energy-coupled versus explicit radiation; original GENESIS reference not rerun'))
    print(json.dumps(rows[-2:],indent=2),flush=True)


def plot(out, old, rows, explicit, energy, gen):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size':12,'axes.titlesize':13,'axes.labelsize':12,
      'legend.fontsize':10,'xtick.labelsize':11,'ytick.labelsize':11,
      'axes.spines.top':False,'axes.spines.right':False,'figure.facecolor':'white',
      'savefig.facecolor':'white','lines.linewidth':1.8})
    colors=['#24648c','#8e4b75','#d47928']
    fig,axes=plt.subplots(2,2,figsize=(12.4,8.2),layout='constrained')
    for i,window in enumerate(['rect','hann']):
        rr=[r for r in rows if r['window']==window];z=[r['z_m'] for r in rr]
        for key,label,color,marker in zip(['explicit_W','energy_coupled_W','GENESIS_W'],
          ['Explicit','Energy coupled','GENESIS'],colors,['o','^','s']):
            axes[i,0].semilogy(z,[r[key] for r in rr],marker+'-',color=color,label=label)
        axes[i,1].plot(z,[r['energy_over_explicit'] for r in rr],'^-',color=colors[1],label='Energy / explicit')
        axes[i,1].axhspan(.95,1.05,color='#cccccc',alpha=.45,label='5% screening band')
        axes[i,1].axhline(1,color='#777777',lw=.9)
        axes[i,0].set_ylabel('Record-mean band power (W)');axes[i,1].set_ylabel('Power ratio')
        for j in range(2):
            axes[i,j].set_xlabel('Fixed laboratory plane z (m)');axes[i,j].grid(alpha=.2)
            axes[i,j].set_title(f'({chr(97+2*i+j)}) '+('Rectangular' if i==0 else 'Energy-normalized Hann'))
            axes[i,j].legend()
    fig.suptitle('SPARC short-bunch SASE: two completed pushers',fontsize=15)
    fig.savefig(out/'radiation_power.png',dpi=210);plt.close(fig)
    photon=sum(read(old/'RADIATION_RESULT.json')['band_eV'])/2
    with np.load(old/'far-spectra.npz') as a,np.load(out/'far-spectra.npz') as b,np.load(old/'far-genesis.npz') as g:
        T=energy[-1]['fft_samples']*energy[-1]['dt_s']
        fig,axes=plt.subplots(2,2,figsize=(12.4,8.2),layout='constrained')
        for i,window in enumerate(['rect','hann']):
            series=[(a['energy_eV'],a[window+'_native_flux_J_per_eV']/T),
                    (b['energy_eV'],b[window+'_native_flux_J_per_eV']/T),
                    (g['energy_eV'],g[window+'_J_per_eV']/T)]
            for j in range(2):
                ax=axes[i,j];lim=(0,5) if j==0 else (photon*.95,photon*1.05)
                hi=max(float(np.max(y[(x>=lim[0])&(x<=lim[1])])) for x,y in series)
                for (x,y),label,color,style in zip(series,['Explicit: native E x B','Energy coupled: native E x B','GENESIS radiation'],colors,['-','--','-']):
                    ax.semilogy(x,np.where(y>0,y,np.nan),style,color=color,label=label)
                ax.axvline(photon,color='#777777',ls=':',lw=1);ax.set_xlim(lim)
                ax.set_ylim(max(hi*(1e-10 if j==0 else 1e-5),1e-20),hi*2)
                ax.set_xlabel('Photon energy (eV)');ax.set_ylabel('Record-mean spectrum (W/eV)')
                ax.set_title(f'({chr(97+2*i+j)}) '+('Rectangular' if i==0 else 'Energy-normalized Hann')+(': wide view' if j==0 else ': resonant band'))
                ax.grid(alpha=.18)
                if i==0 and j==0:ax.legend(loc='upper right')
        fig.suptitle('Outlet at 15.246 m: full physical aperture',fontsize=15)
        fig.savefig(out/'spectra_comparison.png',dpi=210);plt.close(fig)


if __name__=='__main__':
    run(*sys.argv[1:4])
