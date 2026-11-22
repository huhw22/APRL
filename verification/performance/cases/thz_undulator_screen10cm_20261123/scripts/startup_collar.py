"""One-time discrete vacuum-trace continuation, retaining the actual four faces and corners."""
from pathlib import Path
import sys
import numpy as np
from scipy import sparse as sp
from scipy.sparse.linalg import splu
from startup_operator import Grid, BETA, R, E_OFF, B_OFF
NAMES=('ex','ey','ez','bx','by','bz')
def in_collar(sh,ofs,nx,ny,L):
    x=np.arange(sh[0])[:,None]+ofs[0]
    y=np.arange(sh[1])[None,:]+ofs[1]
    return (x<L)|(y<L)|(x>nx-L)|(y>ny-L)

def masks(g):
    me=[in_collar(s,o,g.nx,g.ny,g.layer).ravel() for s,o in zip(g.eshapes,E_OFF)]
    mb=[in_collar(s,o,g.nx,g.ny,g.layer).ravel() for s,o in zip(g.bshapes,B_OFF)]
    return me,mb

def solve(g,vac):
    from time import perf_counter
    tic=perf_counter()
    t=g.t
    f,a=g.harmonic_matrices()
    me,mb=masks(g)
    ue=np.concatenate(me)&np.concatenate(g.free_e)
    ub=np.concatenate(mb)
    vE=np.concatenate([v.ravel() for v in vac[:3]])
    vB=np.concatenate([v.ravel() for v in vac[3:]])
    vE[np.concatenate(me)]=0.;vB[ub]=0.
    ff=f[ub,:][:,ue].tocsr();aa=a[ue,:][:,ub].tocsr()
    ext_e=f[ub,:]@vE;ext_b=a[ue,:]@vB
    matrix=(t*t*sp.eye(ue.sum(),format='csc')-aa@ff).tocsc()
    matrix.eliminate_zeros()
    rhs=t*ext_b+aa@ext_e
    lu=splu(matrix,permc_spec='COLAMD',diag_pivot_thresh=.02)
    e=lu.solve(rhs)
    for _ in range(2):
        b=(ff@e+ext_e)/t
        e+=lu.solve(t*(aa@b+ext_b-t*e))
    b=(ff@e+ext_e)/t
    defect=np.max(abs(aa@b+ext_b-t*e))
    vE[ue]=e;vB[ub]=b
    fields=g.reshape(vE,g.eshapes)+g.reshape(vB,g.bshapes)
    mem={}
    for k,u in enumerate(g.F+g.A):
        if u.axis==2:continue
        family=fields[:3] if k<6 else fields[3:]
        g0=u.D@family[u.src].ravel()
        gain=u.drive/(1-u.decay*np.exp(1j*g.omega))
        omega_phase=g.omega*(1 if k<6 else .5)
        arr=(np.exp(1j*omega_phase)*gain*g0).reshape((g.bshapes if k<6 else g.eshapes)[u.target])
        mem['psi'+str(k)]=arr
    return fields,mem,dict(seconds=perf_counter()-tic,ampere_residual=float(defect),unknowns=int(ue.sum()),nonzeros=int(matrix.nnz),factor_nnz=int(lu.nnz))

