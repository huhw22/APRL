"""Independent spherical Maxwell/Lorentz mean-field reference, before crossing.

In spherical symmetry B=0 in the bunch rest frame. Each shell sees constant
enclosed charge before crossing. Relativistic radial momentum is integrated;
this avoids imposing the nonrelativistic approximation on the PIC comparison.
No APRL stencil, field initialization or particle pusher is called here.
"""
import numpy as np
from scipy.constants import c, e, m_e, epsilon_0
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from scipy.special import gammaincinv

G_BOOST=929.1305492858543
G_BEAM=float(np.hypot(1.,1174.1852350615736))
B_BEAM=float(np.sqrt((G_BEAM-1)*(G_BEAM+1))/G_BEAM)
DZ=1.2543262415359034e-6
DT=DZ/c


class Sphere:
    def __init__(self, kind, electrons=1e8, sigma=1e-4/np.sqrt(5), nr=64, nmu=16, nphi=32):
        self.kind=kind;self.electrons=electrons;self.sigma=sigma
        self.nr=nr;self.nmu=nmu;self.nphi=nphi
        self.eta=electrons*e*e/(4*np.pi*epsilon_0*m_e*c*c*sigma)
        self.t0=sigma/(c*np.sqrt(self.eta))
        q,w=np.polynomial.legendre.leggauss(nr);self.q=(q+1)/2;wr=w/2
        mu,wm=np.polynomial.legendre.leggauss(nmu)
        phi=(np.arange(nphi)+.5)*2*np.pi/nphi
        self.x0=np.sqrt(2*gammaincinv(1.5,self.q)) if kind=='gaussian' else np.sqrt(5)*self.q**(1/3)
        self.shell=np.repeat(np.arange(nr),nmu*nphi)
        mm=np.tile(np.repeat(mu,nphi),nr);pp=np.tile(phi,nr*nmu)
        self.directions=np.c_[np.sqrt(1-mm*mm)*np.cos(pp),np.sqrt(1-mm*mm)*np.sin(pp),mm]
        self.weights=np.repeat(wr,nmu*nphi)*np.tile(np.repeat(wm/2,nphi),nr)/nphi
        self.mu=mu
        self.chirp=-.2 if kind=='converging' else 0.
        initial_beta=self.chirp*self.x0*np.sqrt(self.eta)
        initial_p=initial_beta/np.sqrt(self.eta)/np.sqrt(1-initial_beta**2)
        def rhs(t,s):
            x,p=s[:nr],s[nr:]
            return np.r_[p/np.sqrt(1+self.eta*p*p),self.q/(x*x)]
        self.solution=solve_ivp(rhs,[0,6],np.r_[self.x0,initial_p],method='DOP853',
            rtol=2e-12,atol=2e-13,dense_output=True)
        assert self.solution.success
        self.entry_center_time=1.1*np.max(self.x0)*sigma/c
        # All entry events must occur after the common cold/chirped origin.
        assert B_BEAM*c*self.entry_center_time>np.max(self.x0)*sigma
        # Positive shell ordering is essential for constant enclosed charge.
        self.min_shell_spacing=min(float(np.min(np.diff(self.solution.sol(t)[:nr]))) for t in np.linspace(0,4,81))
        assert self.min_shell_spacing>0
        inv0=np.sqrt(1+self.eta*initial_p**2)+self.eta*self.q/self.x0
        self.energy_invariant_error=max(float(np.max(np.abs(np.sqrt(1+self.eta*self.solution.sol(t)[nr:]**2)
            +self.eta*self.q/self.solution.sol(t)[:nr]-inv0))) for t in np.linspace(0,4,41))
        assert self.energy_invariant_error<1e-10

    def shell_state(self, shell, tau):
        tau=np.asarray(tau);state=self.solution.sol(tau)
        if tau.ndim==0:return state[shell],state[self.nr+shell]
        return state[shell,np.arange(len(tau))],state[self.nr+shell,np.arange(len(tau))]

    def event(self,zplane):
        times=np.empty((self.nr,self.nmu))
        for i in range(self.nr):
            for j,mu in enumerate(self.mu):
                def distance(tau):
                    radius=self.solution.sol(tau)[i]*self.sigma
                    return mu*radius+B_BEAM*c*(tau*self.t0-self.entry_center_time)-zplane/G_BEAM
                assert distance(0)<0 and distance(6)>0
                times[i,j]=brentq(distance,0,6,xtol=5e-15)
        tau=np.repeat(times.ravel(),self.nphi)
        radius,p=self.shell_state(self.shell,tau)
        rest_u=(np.sqrt(self.eta)*p)[:,None]*self.directions
        gamma=np.sqrt(1+np.sum(rest_u*rest_u,axis=1))
        lab_u=rest_u.copy();lab_u[:,2]=G_BEAM*(rest_u[:,2]+B_BEAM*gamma)
        lab_t=B_BEAM*zplane/c+(tau*self.t0-self.entry_center_time)/G_BEAM
        xy=(radius*self.sigma)[:,None]*self.directions[:,:2]
        return {'tau':tau,'lab_time_s':lab_t,'xy_m':xy,'lab_u':lab_u,'rest_u':rest_u,
                'rest_position_m':(radius*self.sigma)[:,None]*self.directions}

    def manifest(self):
        return {'kind':self.kind,'electrons':self.electrons,'sigma_m':self.sigma,
                'nr':self.nr,'nmu':self.nmu,'nphi':self.nphi,'eta':self.eta,'t0_s':self.t0,
                'entry_center_rest_time_s':self.entry_center_time,'beam_gamma_lab':G_BEAM,
                'boost_gamma':G_BOOST,'chirp_times_t0':self.chirp,
                'minimum_shell_spacing_through_tau4':self.min_shell_spacing,
                'maximum_single_shell_energy_invariant_error_mc2':self.energy_invariant_error,
                'reference_scope':'Relativistic spherical continuum, constant enclosed charge before shell crossing; discrete quiet-start quadrature in radius and angle'}
