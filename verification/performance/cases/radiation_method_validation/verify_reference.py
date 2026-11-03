"""Small independent normalization and quadrature checks, no PIC execution."""
import numpy as np
from scipy.constants import c,e,h,mu_0,epsilon_0
from theory import constants,arc_spectrum,energy_density_spectrum
from prepare import require
def verify():
    p=constants(screen=.309);en=np.array([8.,10.,12.]);pts=np.array([[0,0,.309],[7.5e-6,7.5e-6,.309]])
    a,b=arc_spectrum(pts,en,p,256);aa,bb=arc_spectrum(pts,en,p,512)
    error=float(np.linalg.norm(a-aa)/np.linalg.norm(aa));require(error<1e-8,'Arc quadrature not converged')
    n=4096;dt=1e-17;t=(np.arange(n)-n//2)*dt
    wave=np.exp(-.5*(t/1e-15)**2)*np.cos(2*np.pi*c/150e-9*t)
    E=np.zeros((n,3));B=np.zeros_like(E);E[:,0]=wave;B[:,1]=wave/c
    ef=np.fft.rfft(E,axis=0)*dt;bf=np.fft.rfft(B,axis=0)*dt
    photon=h*np.fft.rfftfreq(n,dt)/e;v=energy_density_spectrum(ef,bf);v[[0,-1]]*=.5
    spectral=float(v.sum()*(photon[1]-photon[0]));direct=float(np.sum(wave**2)*dt/(mu_0*c))
    parseval=abs(spectral/direct-1);require(parseval<1e-12,'Independent plane wave Parseval failed')
    # Integrate x'(z) independently, rather than reusing arc_geometry.
    from scipy.integrate import quad
    x=quad(lambda z:z/np.sqrt(p['rho']**2-z*z),0,p['length'],epsabs=1e-20)[0]
    expected=2*p['rho']*np.sin(p['theta']/2)**2
    require(abs(x-expected)<1e-15,'Independent orbit integral differs')
    require(abs(mu_0*epsilon_0*c*c-1)<1e-10,'SI constants inconsistent')
    return {'arc_256_vs_512_relative_L2':error,'synthetic_plane_wave_spectral_vs_time_relative_error':parseval,
      'independent_orbit_x_difference_m':float(x-expected),'scope':'Reference implementation checks, not PIC physical acceptance'}
if __name__=='__main__':
    import json
    print(json.dumps(verify(),indent=2))
