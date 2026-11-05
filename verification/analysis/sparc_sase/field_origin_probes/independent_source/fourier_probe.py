"""Small Fourier-symbol probe of unchanged Cowan vacuum update; no PIC run."""
from pathlib import Path
import json
import numpy as np
OUT=Path(__file__).resolve().parent
DX=DY=4e-5;DZ=8.293564128124092e-6
r=np.array([(DZ/DX)**2,(DZ/DY)**2,1.])
den=r[0]*r[1]+r[1]*r[2]+r[2]*r[0]
beta=r*(1-np.prod(r)/den)/8
# delta[i] is smoothing term in the two axes transverse to derivative i.
delta=np.array([r[1]*r[2]*(1/16-r[1]*r[2]/(8*den)),r[2]*r[0]*(1/16-r[2]*r[0]/(8*den)),r[0]*r[1]*(1/16-r[0]*r[1]/(8*den))])
alpha=np.array([1-2*beta[(i+1)%3]-2*beta[(i+2)%3]-4*delta[i] for i in range(3)])
def symbol(theta):
 theta=np.array(theta);s=np.sin(theta/2);a=2*np.sqrt(r)*s
 smooth=np.array([alpha[i]+2*beta[(i+1)%3]*np.cos(theta[(i+1)%3])+2*beta[(i+2)%3]*np.cos(theta[(i+2)%3])+4*delta[i]*np.cos(theta[(i+1)%3])*np.cos(theta[(i+2)%3]) for i in range(3)])
 return a,a*smooth
rng=np.random.default_rng(20261105);factor_error=0
for theta in rng.uniform(-np.pi,np.pi,(1000,3)):
 a,b=symbol(theta);q=np.dot(a,b);factor_error=max(factor_error,abs(q-4*(1-np.prod(1-r*np.sin(theta/2)**2))))
assert factor_error<5e-15
plane=[]
for tx,ty in [(0,0),(np.pi/60,np.pi/60),(.7,.8),(np.pi,np.pi)]:
 a,b=symbol([tx,ty,np.pi]);plane.append(dict(theta_x=tx,theta_y=ty,q=float(np.dot(a,b))))
# Exact axial 2x2 block in native Fourier phases for E_x and H_y=c*B_y.
kappa=2.;A=np.array([[1-kappa*kappa,-1j*kappa],[-1j*kappa,1]],complex)
M=np.array([[1,1j*kappa/2],[-1j*kappa/2,1]],complex)
assert np.max(np.abs(A.conj().T@M@A-M))==0
assert np.max(np.abs((A+np.eye(2))@(A+np.eye(2))))==0
rows=[]
for n in [0,1,10,100,1000]:
 v=np.linalg.matrix_power(A,n)@np.array([1.,0.],complex)
 native=float(np.vdot(v,v).real);invariant=float(np.vdot(v,M@v).real)
 assert abs(native-(1+4*n+8*n*n))<1e-6
 rows.append(dict(steps=n,native_norm_squared=native,discrete_semidefinite_energy=invariant))
p=5;a=p/6;b=3*p/40+p*(p-1)/72;d=5*p/112+p*(p-1)/80+p*(p-1)*(p-2)/1296
comp=1+a+b+d;ratio=(2/np.pi)**4*comp
result=dict(new_pic_simulation=False,production_aspect_ratios=r.tolist(),coefficient_alpha=alpha.tolist(),coefficient_beta=beta.tolist(),coefficient_delta=delta.tolist(),factorization_max_error=factor_error,nyquist_plane=plane,axial_generalized_mode=rows,compensation_nyquist_gain=comp,continuum_quintic_current_over_CIC_at_nyquist=ratio,limitation='Demonstrates a marginal vacuum mode of the unchanged stencil. Does not prove this mode dominates the actual production fields or identify which new source operation seeds it.')
(OUT/'FOURIER_PROBE.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
