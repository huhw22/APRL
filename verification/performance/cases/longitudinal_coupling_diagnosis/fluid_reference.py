"""Independent periodic cold two-fluid reference and individual characteristics.

The reference has no PIC deposition/gather and no Cowan operator.  It integrates
the 1-D relativistic cold-fluid equations with a spectral Poisson solve and RK4.
Only the weak, stationary control is supported; no production beam is inferred.
"""
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
import argparse,json,time
from pathlib import Path
import numpy as np
from prepare import save,require
C=299792458.;Q=1.602176634e-19;M=9.1093837015e-31;EPS=8.8541878128e-12;L=24*1.2543262415359034e-6;AREA=60e-6*50e-6
def integrate(probes,amplitude,grid,substeps,modes):
 k=2*np.pi/L;omega=np.sqrt(192e5*Q*Q/(AREA*L*M*EPS));ratio=omega/(k*C);a=k*amplitude;xi=np.arange(grid)*2*np.pi/grid;dens=[]
 for sign in [1,-1]:
  origin=xi.copy()
  for _ in range(8):origin-=(origin+sign*a*np.sin(origin)-xi)/(1+sign*a*np.cos(origin))
  dens.append(1/(1+sign*a*np.cos(origin)))
 signs=np.sign(probes['charge_C']);y=np.r_[dens[0],dens[1],np.zeros(2*grid),probes['initial_z_m']*k,np.zeros(len(probes))];freq=np.fft.fftfreq(grid)*grid;inv=np.zeros(grid,dtype=complex);inv[1:]=1/(1j*freq[1:]);select=np.abs(freq)<=modes;selected=freq[select];dt=1.2543262415359034e-6/C;h=omega*dt/substeps
 def rhs(s):
  fluid=s[:4*grid].reshape(4,grid);r=fluid[:2];u=fluid[2:];velocity=u/np.sqrt(1+(ratio*u)**2);e=.5*np.fft.fft(r[0]-r[1])*inv;field=np.fft.ifft(e).real;density_rhs=-np.fft.ifft(1j*freq*np.fft.fft(r*velocity,axis=1),axis=1).real;u_rhs=-velocity*np.fft.ifft(1j*freq*np.fft.fft(u,axis=1),axis=1).real+np.array([1.,-1.])[:,None]*field
  positions=s[4*grid:4*grid+len(probes)];proper=s[4*grid+len(probes):];ep=(np.exp(1j*positions[:,None]*selected[None,:])@(e[select]/grid)).real
  return np.r_[density_rhs.ravel(),u_rhs.ravel(),proper/np.sqrt(1+(ratio*proper)**2),signs*ep]
 for _ in range(160*substeps):
  p=rhs(y);q=rhs(y+.5*h*p);v=rhs(y+.5*h*q);s=rhs(y+h*v);y+=h*(p+2*q+2*v+s)/6.
 return dict(z=y[4*grid:4*grid+len(probes)]/k,u=y[4*grid+len(probes):]*ratio,omega=omega,grid=grid,substeps=substeps,modes=modes)
def main():
 p=argparse.ArgumentParser();p.add_argument('--input',required=True,type=Path);p.add_argument('--output',required=True,type=Path);a=p.parse_args();root=a.input.resolve();out=a.output.resolve();out.mkdir(parents=True,exist_ok=False);base='analytic_linear';points=np.genfromtxt(root/(base+'.csv.particles.csv'),delimiter=',',names=True);plan=json.loads((root/(base+'.json')).read_text());amplitude=plan['result']['displacement_amplitude_m'];require(np.all(points['initial_uz']==0),'Reference is stationary only');records=[];solutions=[]
 for grid,substeps,modes in [(128,16,16),(256,32,32),(256,64,32)]:
  print('START_REFERENCE',grid,substeps,modes,flush=True);t=time.time();s=integrate(points,amplitude,grid,substeps,modes);solutions.append(s);records.append(dict(grid=grid,substeps=substeps,modes=modes,elapsed_s=time.time()-t));np.savez(out/f'reference_{grid}_{substeps}.npz',z=s['z'],u=s['u']);print('COMPLETE_REFERENCE',grid,substeps,flush=True)
 fine=solutions[-1];velocity_scale=np.max(abs(fine['u']));displacement_scale=np.max(abs(fine['z']-points['initial_z_m']));convergence=[dict(configuration=records[i],u_relative_Linf=float(np.max(abs(s['u']-fine['u']))/velocity_scale),z_change_relative_Linf=float(np.max(abs(s['z']-fine['z']))/displacement_scale)) for i,s in enumerate(solutions[:-1])];require(max(x['u_relative_Linf'] for x in convergence)<1e-7,'Fluid reference velocity not converged');require(max(x['z_change_relative_Linf'] for x in convergence)<1e-6,'Fluid reference position not converged');comparisons={}
 for name in ['analytic_linear','analytic_cubic_continuum','analytic_quintic_continuum']:
  data=np.genfromtxt(root/(name+'.csv.particles.csv'),delimiter=',',names=True);require(np.array_equal(data[['id','charge_C','mass_kg','initial_z_m','initial_uz']],points[['id','charge_C','mass_kg','initial_z_m','initial_uz']]),'Input phase space changed');du=data['final_uz']-fine['u'];dz=data['final_z_m']-fine['z'];comparisons[name]=dict(velocity_relative_L2=float(np.linalg.norm(du)/np.linalg.norm(fine['u'])),velocity_relative_Linf=float(np.max(abs(du))/velocity_scale),position_change_relative_L2=float(np.linalg.norm(dz)/np.linalg.norm(fine['z']-points['initial_z_m'])),max_absolute_uz_error=float(np.max(abs(du))),max_absolute_z_error_m=float(np.max(abs(dz))));np.savetxt(out/(name+'.csv'),np.column_stack([data['id'],data['initial_z_m'],fine['z'],data['final_z_m'],fine['u'],data['final_uz']]),delimiter=',',header='id,initial_z_m,reference_final_z_m,PIC_final_z_m,reference_uz,PIC_uz',comments='')
 save(out/'RESULT.json',dict(reference_converged=True,reference='Independent spectral cold-fluid + RK4 characteristics; relativistic 1D equations',configuration=records,convergence=convergence,comparisons=comparisons,velocity_scale=float(velocity_scale),position_change_scale_m=float(displacement_scale),production_modified=False,full_3D_physical_acceptance=False));print(json.dumps(comparisons,indent=2),flush=True)
if __name__=='__main__':main()
