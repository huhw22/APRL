"""Fourier-z reduction of APRL staggered-time Cowan update, transverse PEC and corner CPML."""
from dataclasses import dataclass
from pathlib import Path
import struct
import numpy as np
from scipy import sparse as sp

import os
import yaml
CASE=Path(os.environ['APRL_THZ_CASE']).resolve()
CFG=yaml.safe_load((CASE/'run.yaml').read_text())
NX,NY,NZ=CFG['mesh']['cells']
L=CFG['boundary']['cells'][0]
META=CASE/'preflight/reference/metadata.txt'
META_BITS={r[0]:r[1:] for line in META.read_text().splitlines() if (r:=line.split())}
def values(key):
    return np.array([struct.unpack('>d',bytes.fromhex(x))[0] for x in META_BITS[key+'_bits']])
H=values('spacing')
DT,=values('dt')
C,EPS,MU=values('si')
BETA,=values('beta')
ALPHA=values('cowan_alpha')
B=values('cowan_beta')
D=values('cowan_delta')[[1,2,0]]
R=C*DT/H
V=1/(C*C*EPS*MU)
E_OFF = [(0.5,0,0),(0,0.5,0),(0,0,0.5)]
B_OFF = [(0,0.5,0.5),(0.5,0,0.5),(0.5,0.5,0)]

@dataclass
class Term:
    target: int
    src: int
    axis: int
    sign: int
    scale: float
    D: object
    decay: np.ndarray
    drive: np.ndarray

def add_neighbors(n, normal):
    # PEC halo for E normal (clamp) or tangential (odd), as in electricValue().
    i=np.arange(n)
    minus=i-1; plus=i+1
    minus[0]=0 if normal else 1; plus[-1]=n-1 if normal else n-2
    v=np.ones(2*n)
    if not normal:
        v[0]=v[-1]=-1
    return sp.csr_matrix((v,(np.r_[i,i],np.r_[minus,plus])),shape=(n,n))

class Grid:
    def __init__(self,zeta,nx=NX,ny=NY,layer=L):
        self.zeta,self.nx,self.ny,self.layer=zeta,nx,ny,layer
        self.omega=BETA*R[2]*zeta
        self.t=-2j*np.sin(self.omega/2)
        self.eshapes=[(nx,ny+1),(nx+1,ny),(nx+1,ny+1)]
        self.bshapes=[(nx+1,ny),(nx,ny+1),(nx,ny)]
        self.free_e=[]
        for comp,shape in enumerate(self.eshapes):
            mask=np.ones(shape,dtype=bool)
            if comp!=0: mask[[0,-1],:]=False
            if comp!=1: mask[:,[0,-1]]=False
            self.free_e.append(mask.ravel())
        self.F=[]; self.A=[]
        pairs=[(0,2,1,-1),(0,1,2,1),(1,0,2,-1),(1,2,0,1),(2,1,0,-1),(2,0,1,1)]
        for family,terms,shapes_in,shapes_out in [('F',self.F,self.eshapes,self.bshapes),
                                                ('A',self.A,self.bshapes,self.eshapes)]:
            faraday=(family=='F')
            for target,src,axis,sign in pairs:
                shape=shapes_in[src]
                result_shape=shapes_out[target]
                I=[sp.eye(n,format='csr') for n in shape]
                if axis==2:
                    diff=(2j*np.sin(zeta/2))*sp.eye(np.prod(shape),format='csr')
                else:
                    n=shape[axis]
                    if faraday:
                        difference=sp.diags([-np.ones(n-1),np.ones(n-1)],[0,1],shape=(n-1,n),format='csr')
                    else:
                        difference=sp.diags([np.ones(n),-np.ones(n)],[0,-1],shape=(n+1,n),format='csr')
                    diff=sp.kron(difference,I[1],format='csr') if axis==0 else sp.kron(I[0],difference,format='csr')
                if faraday:
                    neighbors=[add_neighbors(shape[a],src==a) for a in (0,1)]
                    fullI=sp.eye(np.prod(shape),format='csr')
                    if axis==2:
                        smooth=(ALPHA[2]*fullI+B[0]*sp.kron(neighbors[0],I[1],format='csr')+
                                B[1]*sp.kron(I[0],neighbors[1],format='csr')+
                                D[2]*sp.kron(neighbors[0],neighbors[1],format='csr'))
                    else:
                        other=1-axis
                        transverse=(sp.kron(neighbors[0],I[1],format='csr') if other==0
                                    else sp.kron(I[0],neighbors[1],format='csr'))
                        cosz=2*np.cos(zeta)
                        smooth=(ALPHA[axis]+B[2]*cosz)*fullI+(B[other]+D[axis]*cosz)*transverse
                    diff=diff@smooth
                else:
                    diff=diff.multiply(self.free_e[target][:,None]).tocsr()
                decay=np.ones(result_shape)
                if layer and axis<2:
                    nc=[nx,ny][axis]
                    coord=np.arange(result_shape[axis])+(0.5 if faraday else 0.0)
                    depth=np.maximum.reduce([1-coord/layer,1-(nc-coord)/layer,np.zeros_like(coord)])
                    maxrate=-4*R[axis]*np.log(1.e-8)/(2*layer)
                    raw=np.exp(-maxrate*depth**3)
                    decay[:]=raw[:,None] if axis==0 else raw[None,:]
                decay=decay.ravel()
                terms.append(Term(target,src,axis,sign if faraday else -sign,
                                  R[axis]*(1 if faraday else V),diff.tocsr(),decay,decay-1.0))

    def harmonic_matrices(self):
        out=[]
        for terms in (self.F,self.A):
            blocks=[[None]*3 for _ in range(3)]
            for u in terms:
                # Exact harmonic multiplier at the phase of the current raw difference.
                h=u.drive/(1-u.decay*np.exp(1j*self.omega))
                v=u.sign*u.scale*u.D.multiply((1+h)[:,None])
                blocks[u.target][u.src]=v
            out.append(sp.bmat(blocks,format='csr'))
        return out

    def enforce(self,E):
        for x,mask in zip(E,self.free_e): x.reshape(-1)[~mask]=0

    @staticmethod
    def reshape(vector,shapes):
        arrays=[]
        start=0
        for shape in shapes:
            end=start+np.prod(shape)
            arrays.append(vector[start:end].reshape(shape))
            start=end
        assert start==vector.size
        return arrays

    def initial_memory(self,E,Bfield):
        """Memories before the first update, for co-temporal E,cB complex phasors."""
        result=[]
        for terms,data,phase in ((self.F,E,np.exp(1j*self.omega)),
                                 (self.A,Bfield,np.exp(.5j*self.omega))):
            result.append([(u.drive/(1-u.decay*np.exp(1j*self.omega)))*phase*
                           (u.D@data[u.src].ravel()) for u in terms])
        return result
