"""Independent dense spherical force lookup for uniform/chirped spheres.

Interpolate enclosed charge against radius cubed. This avoids the O(dr^2)
piecewise-linear q(r) error of the earlier sparse shell diagnostic. The PIC
algorithm and the independent same-cohort momentum comparison are untouched.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.constants import c,e,epsilon_0

class DenseSphereForce:
    def __init__(self,model,nq=512,nt=1025):
        assert model.kind in ('uniform','converging')
        self.model=model;self.q=np.linspace(0,1,nq+1);self.t=np.linspace(0,3,nt)
        q=self.q[1:];x=np.sqrt(5)*q**(1/3);beta=model.chirp*x*np.sqrt(model.eta)
        p=beta/np.sqrt(model.eta)/np.sqrt(1-beta**2)
        def rhs(t,s):
            x,p=s[:nq],s[nq:]
            return np.r_[p/np.sqrt(1+model.eta*p*p),q/(x*x)]
        sol=solve_ivp(rhs,[0,3],np.r_[x,p],method='DOP853',rtol=2e-12,atol=2e-13,dense_output=True)
        assert sol.success
        self.r3=np.vstack([np.zeros(nt),(model.sigma*sol.sol(self.t)[:nq])**3]);self.nq=nq
        assert np.all(np.diff(self.r3,axis=0)>0)
    def field(self,position,tau):
        position=np.asarray(position,dtype=float);tau=np.asarray(tau,dtype=float)
        assert np.min(tau)>=0 and np.max(tau)<=3
        t=tau/(self.t[1]-self.t[0]);k=np.minimum(t.astype(int),len(self.t)-2);f=t-k
        r2=np.sum(position*position,axis=1);r3=r2*np.sqrt(r2)
        def radial(i):return self.r3[i,k]*(1-f)+self.r3[i,k+1]*f
        lo=np.zeros(len(t),dtype=int);hi=np.full(len(t),self.nq,dtype=int)
        for _ in range(int(np.ceil(np.log2(self.nq)))+1):
            mid=(lo+hi)//2;above=radial(mid)<r3;lo=np.where(above,mid,lo);hi=np.where(above,hi,mid)
        hi=np.maximum(hi,lo+1);hi=np.minimum(hi,self.nq);lo=np.minimum(lo,self.nq-1)
        a,b=radial(lo),radial(hi);q=np.clip(self.q[lo]+(self.q[hi]-self.q[lo])*(r3-a)/(b-a),0,1)
        coefficient=-self.model.electrons*e/(4*np.pi*epsilon_0)*np.divide(q,r3,out=np.zeros_like(q),where=r3>0)
        return coefficient[:,None]*position
