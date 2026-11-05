from pathlib import Path
import json,h5py,numpy as np
root=Path('/public/home/wfollowm/APRL_Workspace/results/sparc_sase_comparison_20261103')
out=root/'momentum_field_origin_20261105/independent_diagnostic';rows=[]
c=299792458.;G=168.80411067346705;beta=.9999824528019369
names=['Ex_boost','Ey_boost','Ez_boost','cBx_boost','cBy_boost','cBz_boost','forward_x','forward_y','backward_x','backward_y']
for case in ['energy','momentum']:
 for plane in ['p1','far']:
  with h5py.File(root/'compact-production-r1'/case/'detectors'/(plane+'.h5'),'r') as f:
   g=f['field_plane'];assert int(g['complete'][()])==1
   z=float(g.attrs['plane_z_m'])
   for start in [10000,20000]:
    v=g['fields'][start:start+256];t=g['time_s'][start:start+256];E=v['electric_V_per_m'];C=v['magnetic_T']*c
    vboost=np.stack([G*(E[...,0]-beta*C[...,1]),G*(E[...,1]+beta*C[...,0]),E[...,2],G*(C[...,0]+beta*E[...,1]),G*(C[...,1]-beta*E[...,0]),C[...,2]],axis=-1)
    phys=vboost[:,10:-10,10:-10];ch=np.stack([phys[...,0]+phys[...,4],phys[...,1]-phys[...,3],phys[...,0]-phys[...,4],phys[...,1]+phys[...,3]],axis=-1)
    a=np.concatenate([phys,ch],axis=-1)
    mean2=np.mean(a*a,axis=(0,1,2));mean=np.mean(a,axis=0);var=np.mean((a-mean[None])**2,axis=(0,1,2));diff=np.mean(np.diff(a,axis=0)**2,axis=(0,1,2));numcorr=np.sum(a[1:]*a[:-1],axis=(0,1,2));den=np.sqrt(np.sum(a[1:]**2,axis=(0,1,2))*np.sum(a[:-1]**2,axis=(0,1,2)))
    tf=np.fft.rfft(a,axis=0);tw=np.abs(tf)**2;tw[1:-1]*=2;power=np.sum(tw,axis=(1,2));fr=np.fft.rfftfreq(256);tweight=np.sum(power,axis=0)
    # Full plane spatial transform including absorbing surround; preserves exact grid Nyquist bins.
    sf=np.fft.fftn(vboost,axes=(1,2));sp=np.mean(np.abs(sf)**2,axis=0);freq=np.fft.fftfreq(80);near=(np.abs(freq[:,None])>=.4)|(np.abs(freq[None,:])>=.4);high=(np.abs(freq[:,None])>=.25)|(np.abs(freq[None,:])>=.25);spall=np.sum(sp,axis=(0,1))
    component=[]
    for j,name in enumerate(names):
     component.append(dict(component=name,rms=float(np.sqrt(mean2[j])),temporal_mean_fraction=float(1-var[j]/mean2[j]),adjacent_difference_rms_over_rms=float(np.sqrt(diff[j]/mean2[j])),adjacent_correlation=float(numcorr[j]/den[j]),temporal_high_above_quarter_sampling_fraction=float(np.sum(power[fr>=.25,j])/tweight[j]),temporal_near_nyquist_fraction=float(np.sum(power[fr>=.4,j])/tweight[j]),spatial_high_fraction=float(np.sum(sp[high,j])/spall[j]) if j<6 else None,spatial_near_nyquist_fraction=float(np.sum(sp[near,j])/spall[j]) if j<6 else None))
    rows.append(dict(case=case,plane=plane,z_m=z,start_index=start,end_index=start+255,time_begin_s=float(t[0]),retarded_begin_s=float(t[0]-z/c),dt_min_s=float(np.min(np.diff(t))),dt_max_s=float(np.max(np.diff(t))),temporal_nyquist_Hz=float(1/(2*np.mean(np.diff(t)))),component_mean_square={n:float(q) for n,q in zip(names,mean2)},components=component))
    print(case,plane,start,'RMS',dict(zip(names,np.sqrt(mean2))),flush=True)
(out/'MODE_PROBE.json').write_text(json.dumps(rows,indent=2)+'\n')
