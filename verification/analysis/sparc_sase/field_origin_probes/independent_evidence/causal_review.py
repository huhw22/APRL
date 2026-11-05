"""Independent bulk-stencil and Jordan-limit review; no PIC simulation."""
from pathlib import Path
import hashlib,json
import numpy as np
base=Path('/public/home/wfollowm/APRL_Workspace')
out=base/'results/sparc_sase_comparison_20261103/momentum_field_origin_20261105/independent_evidence'
source=base/'source/aprl/src/eb_field.cpp'
# Axial Fourier phase theta_z=pi. S_z=1 follows from the exact coefficient
# normalization; all transverse offsets give the same value in this mode.
# Use real grid samples, not the source agent's complex Fourier block.
pattern=np.array([1,-1,1,-1],dtype=np.int64)
def step(x):
    e=x[0]*pattern;h=x[1]*pattern
    hnew=h-(np.roll(e,-1)-e)
    enew=e-(hnew-np.roll(hnew,1))
    assert np.all(enew==enew[0]*pattern) and np.all(hnew==hnew[0]*pattern)
    return np.array([enew[0],hnew[0]],np.int64)
R=np.stack([step(np.array([1,0])),step(np.array([0,1]))],axis=1)
assert np.array_equal(R,np.array([[-3,-2],[2,1]]))
N=R+np.eye(2,dtype=np.int64)
assert np.array_equal(N@N,np.zeros((2,2),np.int64)) and np.linalg.matrix_rank(N)==1
Mr=np.ones((2,2),dtype=np.int64)
assert np.array_equal(R.T@Mr@R,Mr)
# Hphysical=i Hnative_Fourier, for the native z+1/2 By phase.
T=np.diag([1,1j]);A=np.array([[-3,-2j],[-2j,1]],complex)
assert np.max(abs(T@A@np.linalg.inv(T)-R))==0
M=np.array([[1,1j],[-1j,1]],complex)
assert np.max(abs(T.conj().T@Mr@T-M))==0
rows=[]
for n in [0,1,10,100,1000]:
    v=np.array([1,0],np.int64)
    for _ in range(n):v=step(v)
    native=int(v@v);Q=int(v@Mr@v)
    assert native==1+4*n+8*n*n and Q==1
    rows.append(dict(steps=n,E_grid=int(v[0]),H_grid=int(v[1]),norm_squared=native,semidefinite_Q=Q))
# Reconstruct smoothing from its 9 explicit offsets (not cosine shortcut).
dx=dy=4e-5;dz=8.293564128124092e-6
r=np.array([(dz/dx)**2,(dz/dy)**2,1.]);den=r[0]*r[1]+r[1]*r[2]+r[2]*r[0]
beta=r*(1-np.prod(r)/den)/8
D=np.array([r[1]*r[2]*(1/16-r[1]*r[2]/(8*den)),r[2]*r[0]*(1/16-r[2]*r[0]/(8*den)),r[0]*r[1]*(1/16-r[0]*r[1]/(8*den))])
alpha=np.array([1-2*beta[(i+1)%3]-2*beta[(i+2)%3]-4*D[i] for i in range(3)])
def stencil(theta,i):
    j,k=(i+1)%3,(i+2)%3
    v=complex(alpha[i])
    for a in [-1,1]:v+=beta[j]*np.exp(1j*a*theta[j])+beta[k]*np.exp(1j*a*theta[k])
    for a in [-1,1]:
      for b in [-1,1]:v+=D[i]*np.exp(1j*(a*theta[j]+b*theta[k]))
    return v
rng=np.random.default_rng(77311);maximum=0.;nyquist_max=0.
for th in rng.uniform(-np.pi,np.pi,(1000,3)):
    a=2*np.sqrt(r)*np.sin(th/2);b=a*np.array([stencil(th,i) for i in range(3)])
    q=a@b;exact=4*(1-np.prod(1-r*np.sin(th/2)**2));maximum=max(maximum,float(abs(q-exact)))
    th[2]=np.pi;a=2*np.sqrt(r)*np.sin(th/2);b=a*np.array([stencil(th,i) for i in range(3)])
    nyquist_max=max(nyquist_max,float(abs(a@b-4)))
assert maximum<5e-15 and nyquist_max<5e-15
result=dict(new_pic_simulation=False,reviewed_source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),real_sample_update_matrix=R.tolist(),native_phase_matrix_matches=True,jordan_nilpotent_square_zero=True,jordan_rank=int(np.linalg.matrix_rank(N)),semidefinite_Q_eigenvalues=[0,2],invariant_sign_checked=True,norm_growth_formula_checked=True,explicit_sample_results=rows,nine_offset_symbol_factorization_max_error=maximum,nyquist_plane_q_minus_4_max_error=nyquist_max,
 causal_boundary=dict(proved='Unchanged translation-invariant bulk Cowan stencil at cdt=dz contains a defective -1 root at exact longitudinal Nyquist; a generalized initial component grows linearly in amplitude and quadratically in native norm while a semidefinite modified invariant remains constant.',not_proved=['An exponentially unstable vacuum root','Dominance of production k_z=pi energy','A unique seed from high-order shapes, inlet, pending particles or feedback','An eigenmode of the finite CPML and moving-boundary global operator','A source-free late-production growth test']),
 production_scalar_evidence='TE/U_native approaches 2 in both versions, consistent with strongly under-resolved time oscillations but insufficient to identify a unique spatial Fourier mode. Production still has active particles and nonzero J through the final step.')
(out/'CAUSAL_REVIEW.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
