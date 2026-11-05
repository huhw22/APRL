"""Bounded raw-plane E/B fingerprint of an axial Nyquist Cowan mode."""
from pathlib import Path
import json,hashlib,csv
import h5py,numpy as np
D=Path('/public/home/wfollowm/APRL_Workspace/results/sparc_sase_comparison_20261103')
OUT=D/'energy_critical_mode_20261105/diagnostic'
C=np.longdouble(299792458.);G=np.longdouble('168.80411067346705');BETA=np.longdouble('.9999824528019369')
DZ=np.longdouble('8.293564128124092e-6');DT=DZ/C
origin=json.loads((D/'momentum_field_origin_20261105/independent_source/ACTUAL_ADMISSION_EVENTS.json').read_text())
BOX0=np.longdouble(str(origin['box_translation_m']));WIN0=np.longdouble(str(origin['initial_inlet_m']))
with h5py.File(D/'compact-production-r1/energy/detectors/entry-transport-energy.h5') as f:
    T0=np.longdouble(f['energy_plane'].attrs['time_lab_origin_s'])
targets=[0,10,20,50,100,200,400,800,1000,1500,2000,4000,6000,8000,10000,12000,16000,20000,23000]
rows=[];sources=[]
for plane in ['p1','p2','p3','p4','p5','p6','free-015','far']:
    path=D/'compact-production-r1/energy/detectors'/(plane+'.h5');before=path.stat()
    with h5py.File(path,'r') as f:
        g=f['field_plane'];assert int(g['complete'][()])==1
        t=g['time_s'][:].astype(np.longdouble);assert len(t)==len(g['fields'])==int(g['committed_samples'][()]);z=np.longdouble(g.attrs['plane_z_m'])
        # Fixed-lab records at the same Maxwell step differ by beta*z/c.
        n_float=(t-T0-BETA*z/C)*G/DT;n=np.rint(n_float).astype(int)
        assert np.max(np.abs(n_float-n))<1e-3
        for wanted in targets:
            hit=np.flatnonzero(n==wanted)
            if len(hit)!=1:continue
            index=int(hit[0]);v=g['fields'][index];E=v['electric_V_per_m'].astype(np.longdouble);H=v['magnetic_T'].astype(np.longdouble)*C
            ex=G*(E[...,0]-BETA*H[...,1]);ey=G*(E[...,1]+BETA*H[...,0]);hx=G*(H[...,0]+BETA*E[...,1]);hy=G*(H[...,1]-BETA*E[...,0])
            xi=(BOX0-WIN0)/DZ+z/(G*DZ)+(1-BETA)*wanted;frac=float(xi%1)
            ee=1-2*frac;hh=-2*min(frac,1-frac)
            for aperture,sl in [('physical',(slice(10,-10),slice(10,-10))),('beam_core_0p8mm',(slice(30,50),slice(30,50)))]:
                for component,ev,hv in [('x',ex[sl],hy[sl]),('y',ey[sl],-hx[sl])]:
                    for centered in [False,True]:
                        a=ev-ev.mean() if centered else ev;b=hv-hv.mean() if centered else hv
                        aa=np.sum(a*a);bb=np.sum(b*b);ab=np.sum(a*b);norm=aa+bb
                        if norm==0:continue
                        residual=np.sum((ee*b-hh*a)**2)/((ee*ee+hh*hh)*norm)
                        physical_back=np.sum((b+a)**2)/(2*norm);physical_forward=np.sum((b-a)**2)/(2*norm)
                        rows.append(dict(plane=plane,z_lab_m=float(z),maxwell_step=wanted,index=index,region=aperture,spatial_mean_removed=centered,component=component,grid_fraction=frac,nyquist_predicted_H_over_E=hh/ee if abs(ee)>1e-12 else None,observed_regression_H_over_E=float(ab/aa) if aa else None,pair_correlation=float(ab/np.sqrt(aa*bb)) if aa*bb else None,nyquist_template_orthogonal_energy_fraction=float(residual),nyquist_template_residual_amplitude=float(np.sqrt(residual)),nyquist_template_projection_score=float(1-residual),quasistatic_reference_beta_orthogonal_energy_fraction=float((bb+np.longdouble('.5147976239380371')**2*aa-2*np.longdouble('.5147976239380371')*ab)/((1+np.longdouble('.5147976239380371')**2)*norm)),backward_wave_orthogonal_energy_fraction=float(physical_back),forward_wave_orthogonal_energy_fraction=float(physical_forward),E_boost_rms_V_per_m=float(np.sqrt(aa/a.size)),cB_boost_rms_V_per_m=float(np.sqrt(bb/b.size)),pair_mean_square=float(norm/a.size)))
    after=path.stat();assert before.st_size==after.st_size and before.st_mtime_ns==after.st_mtime_ns
    sources.append(dict(path=str(path),bytes=after.st_size,mtime_ns=after.st_mtime_ns,selected_frames=sum(1 for q in targets if np.any(n==q))))
    print(plane,'read',sources[-1]['selected_frames'],'frames',flush=True)
result=dict(scope='Energy-conserving linear SPARC only, bounded closed raw-plane sampling; no new PIC',theory='Axial critical native H_y=-E_x; H_x=E_y. After native node/half-node z interpolation: e=1-2f,h=-2min(f,1-f), so e*H-h*E=0 for H=(cBy,-cBx).',limitation='Projection onto a two-component signature is compatibility, not a unique modal energy fraction. Finite transverse structure, envelopes, physical fields and source response can contribute; sparse z samples cannot establish full kz spectrum.',frame='Inverse Lorentz to native-staggered boost samples; no time centering correction',sources=sources,rows=rows)
(OUT/'NYQUIST_FINGERPRINT.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
with (OUT/'NYQUIST_FINGERPRINT.csv').open('w') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
print('DONE',len(rows),flush=True)
