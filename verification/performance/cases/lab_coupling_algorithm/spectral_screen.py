"""Analytical algorithm screening, not a PIC result or a physical acceptance."""
import argparse,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
def coefficients(d):
    r=(d[2]/d)**2;den=r[0]*r[1]+r[1]*r[2]+r[2]*r[0];b=r*(1-np.prod(r)/den)/8
    delta=np.array([r[1]*r[2]*(1/16-r[1]*r[2]/(8*den)),r[2]*r[0]*(1/16-r[2]*r[0]/(8*den)),r[0]*r[1]*(1/16-r[0]*r[1]/(8*den))]);a=np.array([1-2*b[(i+1)%3]-2*b[(i+2)%3]-4*delta[i] for i in range(3)])
    return a,b,delta
def symbols(theta,coeff):
    a,b,delta=coeff;v=[]
    for i in range(3):
        j,k=(i+1)%3,(i+2)%3;v.append(a[i]+2*b[j]*np.cos(theta[...,j])+2*b[k]*np.cos(theta[...,k])+4*delta[i]*np.cos(theta[...,j])*np.cos(theta[...,k]))
    return np.stack(v,axis=-1)
def cross(k):
    a=np.zeros((3,3));a[0,1]=-k[2];a[0,2]=k[1];a[1,0]=k[2];a[1,2]=-k[0];a[2,0]=-k[1];a[2,1]=k[0];return a
def scalar_periodic_check(d):
    """Independent local-stencil prototype; no FFT or global solve in update."""
    r=(d[2]/d)**2;h=d[2];rng=np.random.default_rng(20261102);shape=(3,12,10,24)
    def forward(x,a):return (np.roll(x,-1,axis=a)-x)/d[a]
    def backward(x,a):return (x-np.roll(x,1,axis=a))/d[a]
    def curl(e,diff):return np.stack([diff(e[2],1)-diff(e[1],2),diff(e[0],2)-diff(e[2],0),diff(e[1],0)-diff(e[0],1)])
    def lap(e,a):return 2*e-np.roll(e,1,axis=a+1)-np.roll(e,-1,axis=a+1)
    def S(e):
        z=lap(e,2);return e-(r[0]*lap(z,0)+r[1]*lap(z,1))/16
    def inner(a,b):return np.sum(a.astype(np.longdouble)*b,dtype=np.longdouble)
    def faraday(e):return -h*curl(S(e),forward)
    def stock(e,b):
        k=faraday(e);bm=b+.5*k
        P=inner(.5*(e[0]+np.roll(e[0],-1,axis=2)),bm[1])-inner(.5*(e[1]+np.roll(e[1],-1,axis=2)),bm[0])
        return .5*(inner(e,S(e))+inner(b,b+k)),P
    e=curl(rng.normal(size=shape),backward);e/=np.sqrt(np.mean(e*e));b=curl(rng.normal(size=shape),forward);b/=np.sqrt(np.mean(b*b));H0,P0=stock(e,b);maxH=maxP=0.
    for n in range(160):
        b+=faraday(e);e+=h*curl(b,backward);H,P=stock(e,b);maxH=max(maxH,float(abs(H-H0)/H0));maxP=max(maxP,float(abs(P-P0)/H0))
    # With arbitrary current, verify its independently calculated weighted work.
    maxWork=0.
    for n in range(32):
        old=e.copy();H0,_=stock(e,b);j=.01*rng.normal(size=shape);b+=faraday(e);e+=h*curl(b,backward)-j;H,_=stock(e,b);work=inner(j,S((old+e)*.5));maxWork=max(maxWork,float(abs(H-H0+work)/H0))
    assert maxH<1e-12 and maxP<1e-12 and maxWork<1e-12
    return dict(vacuum_steps=160,current_steps=32,relative_energy_max=maxH,momentum_over_initial_energy_max=maxP,weighted_current_identity_max=maxWork,field_update='local finite differences and a scalar 13-site spatial filter; no FFT/global iteration',PIC_particles=False,CPML_or_injection_tested=False)
def main():
    p=argparse.ArgumentParser();p.add_argument('output',type=Path);x=p.parse_args();x.output.mkdir(parents=True,exist_ok=False);d=np.array([5e-6,5e-6,1.2543262415359034e-6]);c=299792458.;G=929.1305492858543;beta=np.sqrt(1-1/G**2);dt=d[2]/c;yee_dt=1/(c*np.sqrt(np.sum(1/d**2)));coeff=coefficients(d)
    rng=np.random.default_rng(20261102);metric_error=0.;wrong_metric_error=0.;momentum_metric_error=0.
    for theta in rng.uniform(-np.pi,np.pi,(512,3)):
        k=2*np.sin(theta/2)/d;s=symbols(theta,coeff);C=cross(k);Cs=cross(k*s);ME=np.diag([s[1]*s[2],s[2]*s[0],s[0]*s[1]]);MB=np.diag(s);scale=np.max(abs(C));metric_error=max(metric_error,float(np.max(abs(ME@C-Cs@MB))/scale));wrong_metric_error=max(wrong_metric_error,float(np.max(abs(C-Cs))/scale))
        E=rng.normal(size=3)+1j*rng.normal(size=3);B=rng.normal(size=3)+1j*rng.normal(size=3);E-=k*(k@E)/(k@k);ks=k*s;B-=ks*(ks@B)/(ks@ks);N=np.array([[0,s[1],0],[-s[0],0,0],[0,0,0]])
        defect=np.real(np.vdot(1j*C@B,N@B)-np.vdot(E,N@(1j*Cs@E)));momentum_metric_error=max(momentum_metric_error,float(abs(defect)/(scale*(np.vdot(E,E).real+np.vdot(B,B).real))))
    assert metric_error<1e-14 and wrong_metric_error>1e-3 and momentum_metric_error<1e-14
    q=np.stack(np.meshgrid(*([np.linspace(-np.pi,np.pi,65)]*3),indexing='ij'),axis=-1);S=symbols(q,coeff);M=np.stack([S[...,1]*S[...,2],S[...,2]*S[...,0],S[...,0]*S[...,1]],axis=-1);K=2*np.sin(q/2)/d;stability=float(np.max(dt**2*c*c*np.sum(K*K*S,axis=-1)/4));assert stability<=1+1e-12
    theta=np.linspace(0,1.2e-3,601);hc=1.2398419843320026e-6;lab_energy=91.5;lab_lambda=hc/lab_energy;k0=2*np.pi/lab_lambda;klab=np.stack([k0*np.sin(theta),np.zeros_like(theta),k0*np.cos(theta)],axis=-1);kb=klab.copy();kb[:,2]=G*(klab[:,2]-beta*k0);wprime=G*c*(k0-beta*klab[:,2]);phase=kb*d;K=2*np.sin(phase/2)/d;s=symbols(phase,coeff);assert np.max(abs(phase))<np.pi
    r=(d[2]/d)**2;Q=np.sin(phase/2)**2;scalar=1-Q[:,2]*(r[0]*Q[:,0]+r[1]*Q[:,1]);omega_s=2/dt*np.arcsin(.5*dt*c*np.sqrt(scalar*np.sum(K*K,axis=-1)));err_s=(omega_s/wprime-1)*100
    omega_c=2/dt*np.arcsin(.5*dt*c*np.sqrt(np.sum(K*K*s,axis=-1)));omega_y=2/yee_dt*np.arcsin(.5*yee_dt*c*np.sqrt(np.sum(K*K,axis=-1)));err_c=(omega_c/wprime-1)*100;err_y=(omega_y/wprime-1)*100
    fig,ax=plt.subplots(figsize=(8,4.6),layout='constrained');ax.plot(theta*1e3,err_c,label='Current Cowan');ax.plot(theta*1e3,err_y,label='Yee at its CFL limit');ax.plot(theta*1e3,err_s,label='Experimental scalar-filter pair');ax.set(xlabel='Lab photon angle from z (mrad)',ylabel='Boost-frame frequency error (%)',title='Analytical dispersion at lab 91.5 eV, azimuth 0');ax.axhline(0,color='black',lw=.6);ax.legend();ax.grid(alpha=.2);fig.savefig(x.output/'dispersion.png',dpi=180);fig.savefig(x.output/'dispersion.pdf');plt.close(fig)
    np.savetxt(x.output/'dispersion.csv',np.column_stack([theta,err_c,err_y,(omega_c-wprime)*G/(c*k0)*100,(omega_y-wprime)*G/(c*k0)*100]),delimiter=',',header='lab_angle_rad,Cowan_boost_frequency_error_percent,Yee_boost_frequency_error_percent,Cowan_lab_frequency_projection_error_percent,Yee_lab_frequency_projection_error_percent',comments='')
    result=dict(kind='analytical Fourier screening, not a PIC run',weighted_vacuum_adjoint_relative_error=metric_error,native_metric_adjoint_mismatch=wrong_metric_error,periodic_metric_M_E='diag(Sy*Sz,Sz*Sx,Sx*Sy)',periodic_metric_M_B='diag(Sx,Sy,Sz)',weighted_current_work='J dot (M_E E), not J dot E',metric_S_min=float(np.min(S)),metric_M_E_min=float(np.min(M)),metric_M_E_max=float(np.max(M)),Cowan_CFL_symbol_max=stability,Yee_dt_over_Cowan=yee_dt/dt,Yee_step_count_factor=dt/yee_dt,axial_Cowan_frequency_error_percent=float(err_c[0]),axial_Yee_frequency_error_percent=float(err_y[0]),Cowan_error_at_one_over_beam_gamma_percent=float(np.interp(1/1174,theta,err_c)),Yee_error_at_one_over_beam_gamma_percent=float(np.interp(1/1174,theta,err_y)),physics_changed=False,weighted_metric_is_not_a_physical_accuracy_acceptance=True)
    T=r[0]*np.sin(q[...,0]/2)**2+r[1]*np.sin(q[...,1]/2)**2;Z=np.sin(q[...,2]/2)**2;scalar_stability=float(np.max((Z+T)*(1-T*Z)));assert r[0]+r[1]<np.sqrt(2)-1 and scalar_stability<=1+1e-12
    result.update(Cowan_weighted_momentum_vacuum_error=momentum_metric_error,scalar_pair_CFL_symbol_max=scalar_stability,scalar_pair_production_aspect_sum=float(r[0]+r[1]),scalar_pair_stability_sufficient_bound=float(np.sqrt(2)-1),scalar_pair_filter='S=1-sin(kz*dz/2)^2 * [(dz/dx)^2*sin(kx*dx/2)^2+(dz/dy)^2*sin(ky*dy/2)^2]',scalar_pair_E_metric='S',scalar_pair_B_metric='identity',scalar_pair_current_work='J dot S(Emean)',scalar_pair_axial_frequency_error_percent=float(err_s[0]),scalar_pair_error_at_one_over_beam_gamma_percent=float(np.interp(1/1174,theta,err_s)),scalar_pair_periodic_stencil_test=scalar_periodic_check(d))
    np.savetxt(x.output/'dispersion.csv',np.column_stack([theta,err_c,err_y,err_s]),delimiter=',',header='lab_angle_rad,Cowan_boost_frequency_error_percent,Yee_boost_frequency_error_percent,scalar_pair_boost_frequency_error_percent',comments='')
    (x.output/'RESULT.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
