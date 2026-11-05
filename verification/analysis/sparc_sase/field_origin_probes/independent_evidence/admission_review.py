"""Exact rational independent B-spline continuity/admission check."""
from pathlib import Path
from fractions import Fraction as F
from math import comb,factorial
import json
out=Path(__file__).resolve().parent
# Closed cardinal splines, independent of the C++ recurrence and quadrature.
def B(x,m):
    x=F(x)
    if x<=-F(m+1,2) or x>=F(m+1,2):return F(0)
    return sum((F((-1)**j*comb(m+1,j))*max(F(0),x+F(m+1,2)-j)**m for j in range(m+2)),F(0))/factorial(m)
def CDF(x,m):
    x=F(x)
    if x<=-F(m+1,2):return F(0)
    if x>=F(m+1,2):return F(1)
    return sum((F((-1)**j*comb(m+1,j))*max(F(0),x+F(m+1,2)-j)**(m+1) for j in range(m+2)),F(0))/factorial(m+1)
a=F(5,6);b=F(47,72);d=F(1571,3024)
c=[-d/64,b/16+3*d/32,-a/4-b/4-15*d/64,1+a/2+3*b/8+5*d/16,-a/4-b/4-15*d/64,b/16+3*d/32,-d/64]
def rho(z,k,high):
    if not high:return max(F(0),1-abs(z-k))
    return sum((w*B(z-(k+j-3),5) for j,w in enumerate(c)),F(0))
def current(lo,hi,k,high):
    if not high:return max(F(0),min(hi,F(k+1))-max(lo,F(k)))
    return sum((w*(CDF(hi-(k+j-3)-F(1,2),4)-CDF(lo-(k+j-3)-F(1,2),4)) for j,w in enumerate(c)),F(0))
beta=F('0.51479762393803552');rows=[];details=[]
for phase in ['0.01','0.1','0.2','0.37','0.5']:
    lo=-F(phase);hi=lo+beta
    assert lo<0<hi
    for high in [False,True]:
      rr=[]
      for k in range(-9,11):
        change=rho(hi,k,high)-rho(lo,k,high)
        full=change+current(lo,hi,k,high)-current(lo,hi,k-1,high)
        clip=change+current(F(0),hi,k,high)-current(F(0),hi,k-1,high)
        exact=rho(F(0),k,high)-rho(lo,k,high)
        assert full==0 and clip==exact
        rr.append(dict(k=k,full_residual=full,clipped_residual=clip))
        details.append(dict(phase=phase,high=high,node=k,rho_begin=float(rho(lo,k,high)),rho_end=float(rho(hi,k,high)),full_residual=float(full),clipped_residual=float(clip),exact_expected=float(exact)))
      retained=sum((abs(z['clipped_residual']) for z in rr if z['k']>=1),F(0))
      if not high:assert retained==0
      else:assert retained>0
      rows.append(dict(begin=float(lo),end=float(hi),method='quintic-compensated' if high else 'linear',full_L1_exact=0,clipped_L1=float(sum((abs(z['clipped_residual']) for z in rr),F(0))),retained_old_k_ge_1_L1=float(retained),retained_after_shift_new_k_ge_2_L1=float(retained),max_retained=float(max(abs(z['clipped_residual']) for z in rr if z['k']>=1))))
result=dict(new_pic_simulation=False,arithmetic='Exact rational cardinal B-spline closed forms and exact antiderivatives, no numerical quadrature or source recurrence.',beta=float(beta),full_continuity_exact_all_cases=True,clipped_residual_identity_exact_all_cases=True,identity='R_k = rho_k(z=0) - rho_k(z_begin)',linear_retained_defect_exact_zero=True,quintic_retained_defect_nonzero_all_phases=True,shift='Old k>=1 maps to new k>=2; completion changes Ez[0] only, affecting Gauss at new node1, not new k>=2.',rows=rows,details=details,
 conclusion='A definite per-admission high-order continuity defect survives the moving-inlet single-layer completion. Its quantitative share of the production field anomaly is not yet established without a corrected short replay.')
(out/'ADMISSION_REVIEW.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k!='details'},indent=2))
