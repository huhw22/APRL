"""Closed-substep phase attribution and clone-head/tail diagnostics.

This is a descriptive endpoint analysis. It does not identify a force source.
The exact clone-family decomposition does not require an orbit model; the
separately reported no-feedback reference uses the earlier averaged/ramped
undulator approximation.
"""
import json
import hashlib
from pathlib import Path
import sys

import h5py
import numpy as np
from scipy.constants import c, e, m_e
from scipy.integrate import cumulative_simpson

LD=np.longdouble


def save(path,value):
    Path(path).write_text(json.dumps(value,indent=2,allow_nan=False,
        default=lambda x:x.item() if isinstance(x,np.generic) else x.tolist())+'\n')


def identity(path):
    path=Path(path);s=path.stat();h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
    return dict(path=str(path.resolve()),size=s.st_size,mtime_ns=s.st_mtime_ns,sha256=h.hexdigest())


def mean(x,w):return np.sum(x*w)/np.sum(w)
def rms(x,w):return float(np.sqrt(mean(abs(x)**2,w)))
def corr(x,y,w):
    x=x-mean(x,w);y=y-mean(y,w)
    d=np.sqrt(mean(x*x,w)*mean(y*y,w))
    return float(mean(x*y,w)/d) if d>0 else None
def familymean(x,w):
    return np.sum(x.reshape(-1,16)*w.reshape(-1,16),axis=1)/np.sum(w.reshape(-1,16),axis=1)
def within(x,w):return x-familymean(x,w).repeat(16)
def phase(t,wl):
    angle=np.asarray(np.remainder(np.asarray(t,dtype=LD)*LD(c/wl),1)*LD(2*np.pi),float)
    return np.exp(1j*angle)
def cs(x,w):
    z=mean(x,w);return dict(real=float(z.real),imag=float(z.imag),abs=float(abs(z)),angle_rad=float(np.angle(z)))
def summarize(p,w,group):
    group=np.asarray(group-group.min(),dtype=np.int64);W=np.bincount(group,weights=np.asarray(w,float))
    C=(np.bincount(group,weights=np.asarray(w*p.real,float))+
       1j*np.bincount(group,weights=np.asarray(w*p.imag,float)))
    good=W>0;local=C[good]/W[good];glob=mean(p,w)
    return dict(global_BF=float(abs(glob)),local_RMS_BF=rms(local,W[good]),
        phase_alignment=float(abs(glob)/(np.sum(W[good]*abs(local))/W.sum())))
def harmonic_fit(x,w,phi):
    X=np.c_[within(phi.real,w),within(phi.imag,w)]
    # Float64 solve of a 2x2 system; moments/centering retain long double.
    A=np.asarray(X.T@(w[:,None]*X),float);b=np.asarray(X.T@(w*x),float)
    coeff=np.linalg.solve(A,b);fit=X@coeff
    variance=mean(x*x,w)
    return dict(cos_coefficient=float(coeff[0]),sin_coefficient=float(coeff[1]),
        fitted_peak=float(np.hypot(*coeff)),variance_fraction=float(1-mean((x-fit)**2,w)/variance) if variance>0 else 0.)
def ramp_fit(x,w,r):
    a=mean(x*r,w)/mean(r*r,w);fit=a*r;v=mean(x*x,w)
    return dict(slope_per_optical_cycle=float(a),variance_fraction=float(1-mean((x-fit)**2,w)/v) if v>0 else 0.)


def read(path):
    with h5py.File(path,'r') as f:
        g=f['particle_substep_plane'];assert int(g['complete'][()])==1
        a=g['records'][:];assert len(a)==int(g['committed_records'][()])
    a=a[np.argsort(a['source_id'])];assert np.all(np.diff(a['source_id'])>0)
    assert np.all(a['crossing_direction']==1)
    for n in a.dtype.names:
        if a[n].dtype.kind=='f':assert np.isfinite(a[n]).all()
    return a


def loading_geometry(a,geometry):
    gamma=LD(geometry['boost_gamma']);beta=np.sqrt(1-1/gamma**2)
    uz=a['proper_velocity'][:,2].astype(LD)
    light=(np.hypot(1,uz)+uz)/(gamma+np.sqrt(gamma-1)*np.sqrt(gamma+1))
    held_uz=.5*(light-1/light);held_beta=held_uz/np.hypot(1,held_uz)
    t=a['arrival_time_offset_s'].astype(LD).reshape(-1,16)
    factor=gamma*(beta+held_beta.reshape(-1,16)[:,0])*LD(c)
    ts=np.sort(t,axis=1);length=factor*(ts[:,-1]-ts[:,0]);spacing=factor[:,None]*np.diff(ts,axis=1)
    weights=a['macro_weight'].astype(LD).reshape(-1,16).sum(axis=1)
    def stats(x):return dict(min=float(np.min(x)),median=float(np.median(x)),max=float(np.max(x)))
    return dict(clones_per_family=16,macro_electrons=stats(a['macro_weight']),
        family_electrons=stats(weights),held_boost_rod_length_m=stats(length),
        within_rod_boost_spacing_over_dz=stats(spacing/LD(geometry['dz_m'])),
        transverse_CIC_spacing_m=geometry['dx_m'],
        interpretation='Input clones have identical x,y,u and correlated longitudinal positions over one optical period. Equal-boost-time held-source lengths use the exact existing held z kinematics gamma_b*(beta_b+beta_held)*c*delta_arrival. This is a loading-correlation diagnostic, not proof of self-force or a physical electron rod.')


def run(workspace):
    workspace=Path(workspace).resolve();root=workspace/'results/energy_critical_control_20261105'
    out=root/'bf_origin_review';out.mkdir(exist_ok=True)
    campaign=root/'campaign-r4';src=campaign/'inputs/particles.h5'
    with h5py.File(src,'r') as f:a=f['particles/records'][:]
    ids=a['source_id'];w=a['macro_weight'].astype(LD);N=len(a);assert N%16==0
    assert np.all(w.reshape(-1,16)==w.reshape(-1,16)[:,:1])
    assert np.all(a['proper_velocity'].reshape(-1,16,3)==a['proper_velocity'].reshape(-1,16,3)[:,:1])
    assert np.all(a['plane_position_m'].reshape(-1,16,2)==a['plane_position_m'].reshape(-1,16,2)[:,:1])
    selection=json.loads((campaign/'INPUT_SELECTION.json').read_text())
    plan=json.loads((workspace/'results/sparc_sase_comparison_20261103/compact-production-r1/plan.json').read_text())
    wl=plan['physics']['wavelength_m'];omega=LD(2*np.pi*c/wl);k=2*np.pi/wl
    group=(ids-1)//10240;optical=(ids-1)//1024
    masks=dict(full=np.ones(N,bool),core=(optical>=selection['core_first_slice_inclusive'])&(optical<selection['core_end_slice_exclusive']))
    sources=[src,*[campaign/mode/'detectors'/(p+'-particles-substeps.h5') for mode in ['off','on'] for p in ['entry','und-entry','mid','p1']]]
    before=[identity(p) for p in sources]
    u0=a['proper_velocity'].astype(LD);g0=np.sqrt(1+np.sum(u0*u0,axis=1));q0=1+u0[:,0]**2+u0[:,1]**2
    prep=.01;length=77*.028;K=plan['physics']['K'];ku=2*np.pi/.028
    zz=np.linspace(0,length,77*512+1);dist=np.minimum(zz,length-zz)
    env=np.where(dist<.028,.5*(1-np.cos(np.pi*dist/.028)),1.)
    integ=cumulative_simpson((env*np.cos(ku*zz))**2,x=zz,initial=0.)
    yentry=a['plane_position_m'][:,1]+prep*u0[:,1]/u0[:,2]
    def lag(g,q):
        uz=np.sqrt(g*g-q);return q/(uz*(g+uz))
    def reference_flight(z):
        d,l=min(z,prep),max(0.,z-prep)
        extra=K*K*float(np.interp(l,zz,integ))/l*(1+ku*ku*yentry*yentry) if l else 0.
        return d*lag(g0,q0)+l*lag(g0,q0+extra)
    result=dict(scope='Closed actual-substep endpoint phase analysis; not force-source attribution',
        scalar_energy='gamma_lab and gamma_box are affine endpoint scalars; never gamma reconstructed from interpolated u',
        exact_clone_statement='In continuous dynamics identical initial position, proper velocity, and charge within each 16-clone family imply identical spatial orbit/flight time in static prescribed magnets alone. This symmetry does not by itself bound activation/substep-clock effects of the discrete prescribed-field pusher; a numerical no-feedback control is still needed for that.',
        nofeedback_reference='Earlier averaged-undulator/ramped-end and paraxial transverse-focus approximation; not an exact 3D orbit integration.',
        modes={},new_simulations=0)
    result['input_loading_geometry']=loading_geometry(a,json.loads((campaign/'GEOMETRY.json').read_text()))
    profile_rows=[];harmonics=[]
    for mode in ['off','on']:
        assert json.loads((campaign/mode/'PIC_COMPLETE.json').read_text())['analysis_complete']
        particles={p:read(campaign/mode/'detectors'/(p+'-particles-substeps.h5')) for p in ['entry','und-entry','mid','p1']}
        for p,v in particles.items():assert np.array_equal(v['source_id'],ids) and np.array_equal(v['weight'],w)
        entry,und=particles['entry'],particles['und-entry']
        tentry=entry['time_s'].astype(LD);tund=und['time_s'].astype(LD)
        phaseentry=phase(tentry,wl);phiund=phase(tund-LD(prep/c),wl)
        r=within((tund-tund.min())/LD(wl/c),w)
        tfmin=tund.reshape(-1,16).min(axis=1).repeat(16)
        phirel=phase(tund-tfmin,wl)
        order=np.argsort(tund.reshape(-1,16),axis=1);rank=np.empty((N//16,16),int)
        np.put_along_axis(rank,order,np.tile(np.arange(16),(N//16,1)),axis=1);rank=rank.ravel()
        # Transverse density is only a descriptive covariate, not a self-force.
        xy=a['plane_position_m'].reshape(-1,16,2)[:,0];fw=w.reshape(-1,16).sum(1)
        radius2=np.sum(xy*xy,axis=1);dx=plan['geometry']['dx_m']
        frac=np.remainder(xy/dx,1);self_shape=np.prod(frac**2+(1-frac)**2,axis=1)
        neighbors=np.empty(len(fw),dtype=float)
        for sl in range(len(fw)//64):
            q=xy[sl*64:(sl+1)*64];dd=q[:,None,:]-q[None,:,:]
            matrix=np.exp(-np.sum(dd*dd,axis=2)/(2*dx*dx))
            neighbors[sl*64:(sl+1)*64]=matrix@np.asarray(fw[sl*64:(sl+1)*64],float)-fw[sl*64:(sl+1)*64]
        mode_rows={}
        for pname,z in [('und-entry',.01),('mid',1.088),('p1',2.166)]:
            p=particles[pname];tz=p['time_s'].astype(LD)
            total=omega*(tz-tentry-LD(z/c));total-=mean(total,w)
            fam=familymean(total,w);yp=total-fam.repeat(16)
            start=entry if pname=='und-entry' else und
            dz=z-(0 if pname=='und-entry' else prep)
            delta=omega*(tz-start['time_s'].astype(LD)-LD(dz/c));delta-=mean(delta,w)
            y=within(delta,w)
            dgl=p['gamma_lab'].astype(LD)-start['gamma_lab'].astype(LD)
            dgb=p['gamma_box'].astype(LD)-start['gamma_box'].astype(LD)
            duz=p['proper_velocity'][:,2].astype(LD)-start['proper_velocity'][:,2].astype(LD)
            duzb=p['proper_velocity_box'][:,2].astype(LD)-start['proper_velocity_box'][:,2].astype(LD)
            u=p['proper_velocity'].astype(LD);us=start['proper_velocity'].astype(LD)
            angle=np.sum(u[:,:2]**2,axis=1)/u[:,2]**2-np.sum(us[:,:2]**2,axis=1)/us[:,2]**2
            values=dict(phase_rad=y,delta_energy_lab_eV=within(dgl,w)*LD(m_e*c*c/e),
                delta_gamma_box=within(dgb,w),delta_uz_lab=within(duz,w),delta_uz_box=within(duzb,w),
                delta_angle_squared=within(angle,w))
            ref_flight=reference_flight(z);ref_phi=phaseentry*np.exp(1j*np.asarray(k*(ref_flight-mean(ref_flight,w)),float))
            direct=phase(tz-LD(z/c),wl)
            actual=phaseentry*np.exp(1j*np.asarray(total,float))
            mode_rows[pname]={}
            for scope,mask in masks.items():
                wm=w[mask];gm=group[mask];tr=total[mask];fm=fam.repeat(16)[mask];wi=yp[mask]
                variants=dict(input_phase=phaseentry[mask],no_feedback_approx=ref_phi[mask],
                    family_common_only=phaseentry[mask]*np.exp(1j*np.asarray(fm,float)),
                    within_family_only=phaseentry[mask]*np.exp(1j*np.asarray(wi,float)),actual=actual[mask])
                summaries={n:summarize(v,wm,gm) for n,v in variants.items()}
                direct_b=summarize(direct[mask],wm,gm)
                assert abs(summaries['actual']['global_BF']-direct_b['global_BF'])<2e-8
                assert abs(summaries['actual']['local_RMS_BF']-direct_b['local_RMS_BF'])<2e-8
                row=dict(z_m=z,records=int(mask.sum()),bunching=summaries,
                    segment_start='entry' if pname=='und-entry' else 'und-entry',
                    segment_within_phase_rms_rad=rms(y[mask],wm),
                    segment_within_energy_lab_rms_eV=rms(values['delta_energy_lab_eV'][mask],wm),
                    phase_energy_correlation=corr(y[mask],values['delta_energy_lab_eV'][mask],wm),
                    phase_uz_lab_correlation=corr(y[mask],values['delta_uz_lab'][mask],wm),
                    phase_angle_squared_correlation=corr(y[mask],values['delta_angle_squared'][mask],wm),
                    fits={})
                for label,x in values.items():
                    xm=x[mask];rm=r[mask];pm=phiund[mask];pr=phirel[mask]
                    row['fits'][label]=dict(ramp=ramp_fit(xm,wm,rm),
                        absolute_optical_first_harmonic=harmonic_fit(xm,wm,pm),
                        family_relative_first_harmonic=harmonic_fit(xm,wm,pr))
                    for harmonic in range(1,17):
                        carrier=pm**harmonic;moment=mean(xm*carrier,wm)
                        harmonics.append(dict(mode=mode,plane=pname,scope=scope,quantity=label,harmonic=harmonic,
                            real=float(moment.real),imag=float(moment.imag),abs=float(abs(moment))))
                # Temporal-rank profile directly distinguishes rod head/tail
                # from clone storage order and from a fitted sinusoid.
                for j in range(16):
                    mm=mask&(rank==j); ww=w[mm]
                    profile_rows.append(dict(mode=mode,plane=pname,scope=scope,temporal_rank=j,
                        relative_time_cycles=float(mean(r[mm],ww)),
                        relative_cos=float(mean(within(phirel.real,w)[mm],ww)),
                        relative_sin=float(mean(within(phirel.imag,w)[mm],ww)),
                        absolute_cos=float(mean(within(phiund.real,w)[mm],ww)),
                        absolute_sin=float(mean(within(phiund.imag,w)[mm],ww)),
                        phase_rad=float(mean(y[mm],ww)),
                        delta_energy_lab_eV=float(mean(values['delta_energy_lab_eV'][mm],ww)),
                        delta_uz_box=float(mean(values['delta_uz_box'][mm],ww))))
                # Family slope covariates: no particle or family selection
                # depends on final BF, energy change, or survival.
                af=np.sum((r*y*w).reshape(-1,16),axis=1)/np.sum((r*r*w).reshape(-1,16),axis=1)
                ef=np.sum((r*values['delta_energy_lab_eV']*w).reshape(-1,16),axis=1)/np.sum((r*r*w).reshape(-1,16),axis=1)
                fmask=mask.reshape(-1,16)[:,0]; fww=fw[fmask]
                row['family_ramp_covariates']=dict(
                    phase_slope_mean=float(mean(af[fmask],fww)),phase_slope_rms=rms(af[fmask],fww),
                    positive_phase_slope_charge_fraction=float(np.sum(fww[af[fmask]>0])/fww.sum()),
                    corr_phase_slope_energy_slope=corr(af[fmask],ef[fmask],fww),
                    corr_phase_slope_radius_squared=corr(af[fmask],radius2[fmask],fww),
                    corr_phase_slope_charge=corr(af[fmask],fw[fmask],fww),
                    corr_phase_slope_transverse_CIC_self_shape=corr(af[fmask],self_shape[fmask],fww),
                    corr_phase_slope_same_slice_neighbor_density=corr(af[fmask],neighbors[fmask],fww),
                    family_electrons_range=[float(fww.min()),float(fww.max())],
                    limitation='Endpoint associations; macrocharge range is narrow and varies with longitudinal current. Gaussian 40um same-slice neighbor score is a descriptive proxy, not a resolved Green-function self-field.')
                pr=profile_rows[-16:]
                row['temporal_rank_conditional_profile_fits']={}
                def profile_fit(yv,X):
                    coef=np.linalg.lstsq(X,yv,rcond=None)[0];fit=X@coef
                    return dict(coefficients=coef.tolist(),variance_fraction=float(1-np.sum((yv-fit)**2)/np.sum((yv-yv.mean())**2)))
                for key in ['phase_rad','delta_energy_lab_eV','delta_uz_box']:
                    yv=np.array([q[key] for q in pr])
                    row['temporal_rank_conditional_profile_fits'][key]=dict(
                        ramp=profile_fit(yv,np.c_[np.ones(16),[q['relative_time_cycles'] for q in pr]]),
                        family_relative_first_harmonic=profile_fit(yv,np.c_[np.ones(16),[q['relative_cos'] for q in pr],[q['relative_sin'] for q in pr]]),
                        absolute_optical_first_harmonic=profile_fit(yv,np.c_[np.ones(16),[q['absolute_cos'] for q in pr],[q['absolute_sin'] for q in pr]]))
                segcommon=familymean(delta,w).repeat(16)
                startphase=phase(start['time_s'].astype(LD),wl)
                base=startphase*np.exp(1j*np.asarray(segcommon,float))
                model= row['fits']['phase_rad']
                ramp=model['ramp']['slope_per_optical_cycle']*r
                co=model['absolute_optical_first_harmonic']
                firstfit=within(co['cos_coefficient']*phiund.real+co['sin_coefficient']*phiund.imag,w)
                models=dict(common_only=base,shared_fitted_ramp=base*np.exp(1j*np.asarray(ramp,float)),
                    shared_fitted_first_harmonic=base*np.exp(1j*np.asarray(firstfit,float)),
                    per_family_fitted_ramp=base*np.exp(1j*np.asarray(af.repeat(16)*r,float)))
                row['descriptive_phase_reconstructions']={label:summarize(v[mask],wm,gm) for label,v in models.items()}
                row['descriptive_phase_reconstruction_limit']='Fitted from the observed endpoint phase; not an independent prediction or force-cause fraction.'
                mode_rows[pname][scope]=row
        result['modes'][mode]=mode_rows
    result['profiles']=profile_rows;result['harmonics']=harmonics
    after=[identity(p) for p in sources];assert before==after
    result['sources']=before;result['originals_unchanged']=True
    save(out/'PARTICLE_PHASE_REVIEW.json',result)
    render(result,out)
    print(json.dumps({m:{p:{sc:dict(BF=d['bunching']['actual']['local_RMS_BF'],within_phase_rms=d['segment_within_phase_rms_rad'],phase_energy_corr=d['phase_energy_correlation'],phase_ramp=d['fits']['phase_rad']['ramp'],phase_h1=d['fits']['phase_rad']['absolute_optical_first_harmonic'],ramp_cov=d['family_ramp_covariates']) for sc,d in rr.items()} for p,rr in result['modes'][m].items()} for m in result['modes']},indent=2))


def render(result,out):
    import os
    os.environ.setdefault('MPLCONFIGDIR','/tmp/aprl-bf-particle-phase-mpl')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size':10.5,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(2,2,figsize=(11.8,7.6),layout='constrained')
    for j,pname in enumerate(['mid','p1']):
        for i,key in enumerate(['phase_rad','delta_energy_lab_eV']):
            scale=1. if i==0 else 1e-3
            for mode,color,label in [('off','#24648c','Unfiltered'),('on','#b18428','Seven-point')]:
                rows=[p for p in result['profiles'] if p['mode']==mode and p['plane']==pname and p['scope']=='core']
                xx=np.array([p['relative_time_cycles'] for p in rows]);yy=np.array([p[key] for p in rows])*scale
                axes[i,j].plot(xx,yy,'o-',ms=3.5,color=color,label=label)
            fits=result['modes']['on'][pname]['core']['temporal_rank_conditional_profile_fits'][key]
            ramp=fits['ramp']['coefficients'];sine=fits['family_relative_first_harmonic']['coefficients']
            axes[i,j].plot(xx,(ramp[0]+ramp[1]*xx)*scale,'--',color='#555555',label='Ramp fit (seven-point)')
            pred=sine[0]+sine[1]*np.array([p['relative_cos'] for p in rows])+sine[2]*np.array([p['relative_sin'] for p in rows])
            axes[i,j].plot(xx,pred*scale,':',color='#668d6b',lw=2,label='First-harmonic fit')
            axes[i,j].set_title(('z=1.088 m' if pname=='mid' else 'z=2.166 m')+
                f"; profile R² ramp={fits['ramp']['variance_fraction']:.3f}, sine={fits['family_relative_first_harmonic']['variance_fraction']:.3f}")
            axes[i,j].set_ylabel('Within-family phase change (rad)' if i==0 else 'Within-family energy change (keV)')
            axes[i,j].set_xlabel('Entrance head → tail position (optical cycles)')
            axes[i,j].axhline(0,color='#999999',lw=.7);axes[i,j].grid(alpha=.2)
            axes[i,j].legend(fontsize=7.5)
    fig.suptitle('Observed clone-family head/tail modulation — central 77 slices\nActual-substep arrival times and affine scalar energy; endpoint associations, not force attribution',fontsize=12)
    fig.savefig(out/'particle_phase_profiles.png',dpi=200)
    fig.savefig(out/'particle_phase_profiles.pdf')
    plt.close(fig)


if __name__=='__main__':run(sys.argv[1])
