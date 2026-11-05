"""Closed energy/momentum CSV waveform comparison; no HDF5 or dynamics."""
import csv,json,hashlib
from pathlib import Path
import numpy as np

TOP=Path('/public/home/wfollowm/APRL_Workspace/results/sparc_sase_comparison_20261103')
OUT=TOP/'momentum_startup_assessment_20261105/completed_waveform'
BASE={'energy':TOP/'energy_startup_origin_20261105/r2/analysis',
      'momentum':TOP/'momentum_startup_assessment_20261105/deferred/analysis/momentum'}

def read(p):return list(csv.DictReader(p.open()))
def fit(X,y,w):
    b=np.linalg.solve(X.T@(w[:,None]*X),X.T@(w*y));res=y-X@b
    total=w@((y-w@y)**2)
    return dict(coefficients=b.tolist(),profile_R2=float(1-w@(res*res)/total),
                residual_rms=float(np.sqrt(w@(res*res)))),X@b

report={'case_scope':['energy','momentum'],'new_simulation':False,
        'HDF5_opened':False,'source_manifest':[],'cases':{}}
phase_profiles={};phase_fits={};harmonic_arrays={}
for case,base in BASE.items():
    assert json.loads((base/'COMPLETE.json').read_text())['complete']
    assert json.loads((base/'RESULT.json').read_text())['case']==case
    for name in ['phase_bins.csv','harmonics.csv','RESULT.json']:
        p=base/name;report['source_manifest'].append(dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
    rows=read(base/'phase_bins.csv');x=np.array([float(r['phase_rad']) for r in rows]);w=np.array([float(r['charge_fraction']) for r in rows]);w/=w.sum()
    case_report={}
    for field in ['phase_within_mean_rad','energy_within_mean_eV']:
        y=np.array([float(r[field]) for r in rows])
        models={'periodic_ramp':np.column_stack([np.ones(len(x)),x-np.pi]),
            'harmonic_1':np.column_stack([np.ones(len(x)),np.cos(x),np.sin(x)]),
            'harmonics_1_to_3':np.column_stack([np.ones(len(x))]+[fn(n*x) for n in range(1,4) for fn in [np.cos,np.sin]]),
            'harmonics_1_to_8':np.column_stack([np.ones(len(x))]+[fn(n*x) for n in range(1,9) for fn in [np.cos,np.sin]]),
            'harmonics_10_20':np.column_stack([np.ones(len(x))]+[fn(n*x) for n in [10,20] for fn in [np.cos,np.sin]]),
            'harmonics_1_to_8_plus_10_20':np.column_stack([np.ones(len(x))]+[fn(n*x) for n in [*range(1,9),10,20] for fn in [np.cos,np.sin]])}
        case_report[field]={}
        for name,X in models.items():
            ff,curve=fit(X,y,w);case_report[field][name]=ff
            if field=='phase_within_mean_rad':phase_fits[(case,name)]=curve
        if field=='phase_within_mean_rad':phase_profiles[case]=(x,y,w)
    harmonics=read(base/'harmonics.csv');harmonic_arrays[case]={}
    for quantity in ['delta_gamma_within','phase_residual_within']:
        rows=[r for r in harmonics if r['plane']=='p1' and r['quantity']==quantity]
        nall=np.array([int(r['harmonic']) for r in rows]);zall=np.array([complex(float(r['global_real']),float(r['global_imag'])) for r in rows])
        harmonic_arrays[case][quantity]=(nall,zall)
        use=nall<=8;n=nall[use];z=zall[use]
        jumps=np.linspace(0,2*np.pi,20001)
        basis=-1j*np.exp(1j*jumps[:,None]*n[None,:])/n[None,:]
        slopes=np.real(np.sum(basis.conj()*z,axis=1))/np.sum(abs(basis)**2,axis=1)
        err=np.sum(abs(z[None,:]-slopes[:,None]*basis)**2,axis=1);best=int(np.argmin(err))
        case_report[quantity+'_n1_to_n8']={'periodic_ramp_spectral_fit_fraction':float(1-err[best]/np.sum(abs(z)**2)),
            'slope':float(slopes[best]),'jump_rad':float(jumps[best]),
            'amplitude_ratios_to_n1':(abs(z)/abs(z[0])).tolist()}
    report['cases'][case]=case_report

def corr(x,y,w):
    x=x-w@x;y=y-w@y
    return float((w@(x*y))/np.sqrt((w@(x*x))*(w@(y*y))))
xe,ye,we=phase_profiles['energy'];xm,ym,wm=phase_profiles['momentum']
assert np.array_equal(xe,xm)
comparison={'raw_profile_correlation':corr(ye,ym,.5*(we+wm))}
for kind in ['periodic_ramp','harmonic_1','harmonics_1_to_3','harmonics_1_to_8']:
    comparison[kind+'_profile_correlation']=corr(phase_fits[('energy',kind)],phase_fits[('momentum',kind)],.5*(we+wm))
for quantity in ['delta_gamma_within','phase_residual_within']:
    ze=harmonic_arrays['energy'][quantity][1][:8];zm=harmonic_arrays['momentum'][quantity][1][:8]
    ratio=np.vdot(ze,zm)/np.vdot(ze,ze)
    comparison[quantity+'_complex_n1_to_n8']=dict(best_complex_scale_real=float(ratio.real),best_complex_scale_imag=float(ratio.imag),
        coefficient_residual_fraction=float(np.linalg.norm(zm-ratio*ze)/np.linalg.norm(zm)),
        normalized_complex_overlap=float(abs(np.vdot(ze,zm))/np.linalg.norm(ze)/np.linalg.norm(zm)))
    all_e=harmonic_arrays['energy'][quantity][1];all_m=harmonic_arrays['momentum'][quantity][1]
    comparison[quantity+'_selected_harmonics']=[dict(n=n,
        energy_abs=float(abs(all_e[n-1])),momentum_abs=float(abs(all_m[n-1])),
        momentum_over_energy=float(abs(all_m[n-1])/abs(all_e[n-1]))) for n in [1,2,3,10,20,30]]
    comparison[quantity+'_n10_20_30_spectral_norm']=dict(
        energy=float(np.linalg.norm(all_e[np.array([10,20,30])-1])),
        momentum=float(np.linalg.norm(all_m[np.array([10,20,30])-1])),
        momentum_over_energy=float(np.linalg.norm(all_m[np.array([10,20,30])-1])/np.linalg.norm(all_e[np.array([10,20,30])-1])))
report['comparison']=comparison
report['limits']='Profile fits describe 48 conditional means; spectral fits describe n=1..8 only. Neither is a particle-variance fraction, force-source share, or independent growth prediction. Raw reconstructed energy has known diagnostic limitations. Shared startup shape does not exclude shared numerical error.'
(OUT/'RESULT.json').write_text(json.dumps(report,indent=2)+'\n')

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False})
fig,axes=plt.subplots(1,3,figsize=(15,4.6),layout='constrained')
style={'energy':('#24648c','Energy coupled'),'momentum':('#8e4b75','Momentum matched')}
for case,(color,label) in style.items():
    x,y,w=phase_profiles[case]
    axes[0].plot(x,y,'o-',ms=3,color=color,label=label)
    axes[1].plot(x,phase_fits[(case,'harmonics_1_to_3')],'-',color=color,label=label)
    n,z=harmonic_arrays[case]['phase_residual_within']
    axes[2].semilogy(n,abs(z),'o-',ms=3,color=color,label=label)
axes[0].set_title('(a) Full conditional phase mean')
axes[1].set_title('(b) Harmonics 1-3 of that mean')
for ax in axes[:2]:
    ax.set_xlabel('Entrance optical phase (rad)');ax.set_ylabel('Mean within-family phase residual (rad)')
    ax.set_xticks([0,np.pi,2*np.pi],['0',r'$\pi$',r'$2\pi$'])
axes[2].set_title('(c) Exact-particle phase harmonics')
axes[2].set_xlabel('Harmonic n');axes[2].set_ylabel('Complex harmonic magnitude (rad)')
for ax in axes:ax.grid(alpha=.18);ax.legend(fontsize=9)
fig.suptitle('First-segment startup waveform: shared low harmonics, different fine structure',fontsize=13)
fig.savefig(OUT/'waveform_comparison.png',dpi=220,bbox_inches='tight');plt.close(fig)
print(json.dumps(report,indent=2))
