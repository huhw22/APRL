"""Independent inverse checks using forward targets and scalar Brent roots."""
import numpy as np
from scipy.optimize import brentq
from theory import constants, arrival_length, velocity_fields

def verify():
    p=constants(screen=.309)
    theta=np.linspace(p['theta']*.00001,p['theta']*.99999,31)
    points=np.c_[np.linspace(-100e-6,100e-6,31),
                 3e-6+17e-6*np.sin(np.arange(31)),np.full(31,p['screen'])]
    targets=arrival_length(theta,points,p)[0]
    roots=np.array([brentq(lambda q:arrival_length(q,point,p)[0]-target,
                          0,p['theta'],xtol=1e-17,rtol=1e-14)
                    for point,target in zip(points,targets)])
    assert max(abs(roots-theta))<1e-13
    # Dense observer times exercise incoming, arc, and outgoing branches.
    us=np.linspace(-.7e-15,9e-15,1700)
    E,B=velocity_fields(points,us,p)
    assert np.isfinite(E).all() and np.isfinite(B).all()
    # Very high gamma and points close to the arc test poor time conditioning.
    for gamma in [1174.171072271572,10000.]:
        q=constants(gamma=gamma,screen=.309)
        velocity_fields(np.array([[0.,2.5e-6,.309],[7.5e-6,-7.5e-6,.309]]),us,q)
    return {'forward_vs_independent_brent_max_angle_error_rad':float(max(abs(roots-theta))),
            'finite_velocity_samples':int(E.size),'pass':True,
            'scope':'Theory inverse and finite-value checks, not PIC accuracy acceptance'}

if __name__=='__main__':
    import json
    print(json.dumps(verify(),indent=2))
