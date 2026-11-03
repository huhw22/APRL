"""Theory-only audit of the exported mixed-time Lorentz observable."""
import argparse,json
from pathlib import Path
import numpy as np
from scipy.constants import c,e,epsilon_0,h,mu_0
from theory import arrival_length,velocity_fields,energy_density_spectrum
from prepare import save

def full_fields(points,u,p):
    E,B=velocity_fields(points,u,p)
    initial=(np.sum(points[:,:2]**2,axis=1)/
             (np.linalg.norm(points,axis=1)+p['screen']))/c
    finish=arrival_length(p['theta'],points,p)[0]/c
    it,ip=np.nonzero((u[:,None]>initial)&(u[:,None]<finish))
    if len(it):
        pts=points[ip];target=c*u[it]
        lo=np.zeros(len(it));hi=np.full(len(it),p['theta'])
        for _ in range(60):
            theta=(lo+hi)/2
            val,R,n,beta,normal=arrival_length(theta,pts,p)
            hi=np.where(val>target,theta,hi);lo=np.where(val<=target,theta,lo)
        theta=(lo+hi)/2
        val,R,n,beta,normal=arrival_length(theta,pts,p)
        kappa=1-np.sum(n*beta,axis=-1)
        acc=(-e*p['beta']**2/(4*np.pi*epsilon_0*p['rho']))*np.cross(n,np.cross(n-beta,normal))/(kappa**3*R)[:,None]
        E[it,ip]+=acc;B[it,ip]+=np.cross(n,acc)/c
    return E,B

def to_boost(E,B,g,b):
    outE=E.copy();outB=B.copy()
    outE[...,0]=g*(E[...,0]-b*c*B[...,1])
    outE[...,1]=g*(E[...,1]+b*c*B[...,0])
    outB[...,0]=g*(B[...,0]+b*E[...,1]/c)
    outB[...,1]=g*(B[...,1]-b*E[...,0]/c)
    return outE,outB

def to_lab(E,B,g,b):
    return to_boost(E,B,g,-b)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('case',type=Path);a=ap.parse_args()
    plan=json.loads((a.case/'PLAN.json').read_text())
    src=np.load(a.case/'analysis/field_spectra.npz')
    result=json.loads((a.case/'analysis/RESULT.json').read_text())
    times=src['time_s'];n=len(times);dt=(times[-1]-times[0])/(n-1)
    first=result['time_first_relative_to_magnet_light_s']
    u=first+np.arange(n)*dt
    xx,yy=np.meshgrid(src['x_m'],src['y_m'])
    roi=src['roi'];p=plan['parameters']
    pts=np.c_[xx[roi],yy[roi],np.full(int(roi.sum()),p['screen'])]
    E,B=full_fields(pts,u,p)
    g=929.1305492858543;b=np.sqrt((g-1)*(g+1))/g
    Eb,Bb=to_boost(E,B,g,b)
    photon=h*np.fft.rfftfreq(n,dt)/e;de=photon[1]-photon[0]
    area=result['registered_aperture_area_m2']/roi.sum()
    def spectrum(v,w):
        F=np.conj(np.fft.rfft(v,axis=0))*dt
        G=np.conj(np.fft.rfft(w,axis=0))*dt
        d=energy_density_spectrum(F,G).sum(axis=1)*area
        d[0]*=.5
        if n%2==0:d[-1]*=.5
        lo,hi=plan['registered_band_eV'];band=(photon>=lo)&(photon<=hi)
        return float(d[band].sum()*de)
    rows={'synchronous_point_samples_J':spectrum(E,B)}
    for sign in [-1,1]:
        shift=sign*1.2543262415359034e-6/c/2
        p2=dict(p,screen=p['screen']+g*b*c*shift)
        pts2=pts.copy();pts2[:,2]=p2['screen']
        E2,B2=full_fields(pts2,u+shift/(g*(1+b)),p2)
        E2b,B2b=to_boost(E2,B2,g,b)
        rawE,rawB=to_lab(Eb,B2b,g,b)
        rows['B_boost_half_time_'+str(sign)+'_J']=spectrum(rawE,rawB)
    rows['continuous_fourier_LW_J']=result['complete_LW_reference_band_J']
    rows['PIC_J']=result['PIC_band_J']
    rows['scope']='Theory sampling audit; no C++ changes and no claim that this alone explains PIC disagreement'
    save(a.case/'analysis/SAMPLED_LW_AUDIT.json',rows)
    print(json.dumps(rows,indent=2))

if __name__=='__main__':main()
