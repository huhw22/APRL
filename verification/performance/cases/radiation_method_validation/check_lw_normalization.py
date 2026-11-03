"""Independent Larmor angular integral and zero-frequency endpoint identity."""
import numpy as np
from scipy.constants import c,e,epsilon_0
from scipy.special import roots_legendre
from theory import constants,arc_spectrum

def verify():
    p=constants(length=.0024,screen=1e6)
    g,b,R=p['gamma'],p['beta'],p['rho']
    u,w=roots_legendre(80)
    phi=(np.arange(64)+.5)*2*np.pi/64
    uu,pp=np.meshgrid(u,phi,indexing='ij')
    den=1+b*uu
    n=np.stack([np.sqrt(1-uu**2)*np.cos(pp)/(g*den),
                np.sqrt(1-uu**2)*np.sin(pp)/(g*den),
                (uu+b)/den],-1)
    velocity=np.array([0.,0.,b])
    acceleration=np.array([b*b*c/R,0.,0.])
    kappa=1/(g*g*den)
    cross=np.cross(n,np.cross(n-velocity,acceleration))
    dp=e*e/(16*np.pi**2*epsilon_0*c)*np.sum(cross*cross,axis=-1)/kappa**5
    jacobian=1/(g*g*den*den)
    integrated=float(np.sum(dp*jacobian*w[:,None])*(2*np.pi/64))
    larmor=e*e*c*b**4*g**4/(6*np.pi*epsilon_0*R*R)
    field,_=arc_spectrum(np.array([[0.,0.,p['screen']]]),np.array([0.]),p,512)
    actual=field[0,0,0]*p['screen']
    expected=(-e/(4*np.pi*epsilon_0*c))*(-b*np.sin(p['theta'])/(1-b*np.cos(p['theta'])))
    result={'Larmor_power_W':larmor,'angular_integral_power_W':integrated,
       'Larmor_relative_error':abs(integrated/larmor-1),
       'zero_frequency_far_field_endpoint_relative_error':float(abs(actual/expected-1)),
       'scope':'Independent theory normalization checks; not PIC radiation acceptance'}
    assert result['Larmor_relative_error']<1e-8
    assert result['zero_frequency_far_field_endpoint_relative_error']<1e-7
    return result

if __name__=='__main__':
    import json
    print(json.dumps(verify(),indent=2))
