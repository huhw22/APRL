"""Independent prescribed-orbit Lienard-Wiechert reference, SI throughout.

No APRL pusher, field solver, deposition, or initial-field routine is imported.
The arc is bounded by planar hard-edge z faces, exactly as the test card.
"""
import numpy as np
from scipy.constants import c,e,m_e,epsilon_0,h,hbar,mu_0
from scipy.special import roots_hermitenorm

def constants(gamma=1174.171072271572,rho=3.,length=.009,screen=.059):
    beta=np.sqrt((gamma-1)*(gamma+1))/gamma
    theta=np.arcsin(length/rho)
    result=dict(gamma=gamma,beta=beta,rho=rho,length=length,screen=screen,
      theta=theta,field_T=gamma*beta*m_e*c/(e*rho),
      critical_eV=1.5*hbar*c*gamma**3/(rho*e))
    return {k:float(v) for k,v in result.items()}

def arc_geometry(theta,p):
    theta=np.asarray(theta);r=p['rho'];b=p['beta']
    pos=np.stack([2*r*np.sin(theta/2)**2,np.zeros_like(theta),r*np.sin(theta)],axis=-1)
    vel=np.stack([b*np.sin(theta),np.zeros_like(theta),b*np.cos(theta)],axis=-1)
    normal=np.stack([np.cos(theta),np.zeros_like(theta),-np.sin(theta)],axis=-1)
    return pos,vel,normal

def arrival_length(theta,points,p):
    """c*(arrival_time - screen_z/c), avoiding ultrarelativistic cancellation."""
    pos,vel,normal=arc_geometry(theta,p)
    dr=points-pos;R=np.linalg.norm(dr,axis=-1)
    th=np.asarray(theta);minus_sin=th**3*(1/6-th**2/120+th**4/5040-th**6/362880)
    invbeta_minus1=1/(p['gamma']**2*p['beta']*(1+p['beta']))
    arrival=p['rho']*(th*invbeta_minus1+minus_sin)+(dr[...,0]**2+dr[...,1]**2)/(R+dr[...,2])
    return arrival,R,dr/R[...,None],vel,normal

def arc_spectrum(points,energy_eV,p,ntheta=512):
    """Exact finite-distance acceleration-field Fourier integral over the arc.

    Fourier convention integral E(t) exp(+i omega*(t-z_screen/c)) dt.
    Jacobian dt_observer=(1-n.beta) dt_emitter is included explicitly.
    """
    x,w=np.polynomial.legendre.leggauss(ntheta);theta=(x+1)*p['theta']/2;w=w*p['theta']/2
    outE=np.zeros((len(energy_eV),len(points),3),complex);outB=np.zeros_like(outE)
    omega=np.asarray(energy_eV)*e/hbar
    for start in range(0,len(points),128):
        pts=points[start:start+128,None,:]
        arrival,R,n,beta,normal=arrival_length(theta,pts,p)
        kappa=1-np.sum(n*beta,axis=-1)
        kernel=(-e*p['beta']/(4*np.pi*epsilon_0*c))*np.cross(n,np.cross(n-beta,normal))/(kappa**2*R)[...,None]
        bk=np.cross(n,kernel)/c
        phase=np.exp(1j*omega[:,None,None]*arrival[None,:,:]/c)*w[None,None,:]
        outE[:,start:start+len(pts)]=np.einsum('fpt,ptd->fpd',phase,kernel,optimize=True)
        outB[:,start:start+len(pts)]=np.einsum('fpt,ptd->fpd',phase,bk,optimize=True)
    return outE,outB

def uniform_velocity_fields(points,u,p,outgoing):
    """Exact Heaviside field of the appropriate infinite straight worldline."""
    pts=points[None,:,:];tt=np.asarray(u)[:,None];b=p['beta'];g=p['gamma'];Z=p['screen']
    if outgoing:
        end,vel,_=arc_geometry(np.array(p['theta']),p);axis=vel/b
        dr=pts-end;para=np.sum(dr*axis,axis=-1)
        perp=dr-para[...,None]*axis
        along=para-b*(Z+c*tt-p['rho']*p['theta']/b)
    else:
        axis=np.array([0.,0.,1.]);perp=pts.copy();perp[...,2]=0.
        along=Z/(g*g*(1+b))-b*c*tt
    R=perp+along[...,None]*axis
    denominator=along**2+np.sum(perp**2,axis=-1)/(g*g)
    E=(-e/(4*np.pi*epsilon_0*g*g))*R/denominator[...,None]**1.5
    B=np.cross(b*axis,E)/c
    return E,B

def velocity_fields(points,u,p):
    """Piecewise exact velocity term, with a safeguarded retarded arc solve."""
    u=np.asarray(u);points=np.asarray(points)
    start=np.sum(points[:,:2]**2,axis=1)/(np.linalg.norm(points,axis=1)+p['screen'])/c
    finish=arrival_length(np.array(p['theta']),points,p)[0]/c
    incoming=u[:,None]<=start[None,:];arc=(u[:,None]>start[None,:])&(u[:,None]<finish[None,:])
    E,B=uniform_velocity_fields(points,u,p,True)
    Ei,Bi=uniform_velocity_fields(points,u,p,False)
    E[incoming]=Ei[incoming];B[incoming]=Bi[incoming]
    it,ip=np.nonzero(arc)
    if len(it):
        pts=points[ip];target=c*u[it]
        lo=np.zeros(len(it));hi=np.full(len(it),p['theta'])
        # A monotone bracket avoids moving an already converged Newton root
        # back to a distant midpoint when strict bracket membership rejects it.
        # Sixty bisections bound angular error independently of a small residual
        # or the ultrarelativistic conditioning of observer time.
        theta=(lo+hi)/2
        for _ in range(60):
            val,R,n,beta,_=arrival_length(theta,pts,p);delta=val-target
            hi=np.where(delta>0,theta,hi);lo=np.where(delta<=0,theta,lo)
            theta=(lo+hi)/2
        val,R,n,beta,_=arrival_length(theta,pts,p)
        if np.max(abs(val-target))>2e-16:raise ValueError('Retarded arc inversion failed')
        kappa=1-np.sum(n*beta,axis=-1)
        Ev=(-e/(4*np.pi*epsilon_0*p['gamma']**2))*(n-beta)/(kappa**3*R**2)[:,None]
        E[it,ip]=Ev;B[it,ip]=np.cross(n,Ev)/c
    return E,B

def velocity_spectrum(points,energy_eV,p,time_first=-.9e-15,time_span=20e-15,ntime=8192,sigma_time=0.):
    """Same finite record as PIC; uniform direct quadrature, not a fitted pulse."""
    from scipy.ndimage import gaussian_filter1d
    dt=time_span/ntime;u=time_first+np.arange(ntime)*dt
    margin=int(np.ceil(8*sigma_time/dt)) if sigma_time else 0
    ue=time_first+np.arange(-margin,ntime+margin)*dt
    omega=np.asarray(energy_eV)*e/hbar
    phase=np.exp(1j*omega[:,None]*u[None,:])*dt
    Eout=np.zeros((len(omega),len(points),3),complex);Bout=np.zeros_like(Eout)
    for start in range(0,len(points),64):
        E,B=velocity_fields(points[start:start+64],ue,p)
        if margin:
            E=gaussian_filter1d(E,sigma_time/dt,axis=0,truncate=8,mode='nearest')[margin:margin+ntime]
            B=gaussian_filter1d(B,sigma_time/dt,axis=0,truncate=8,mode='nearest')[margin:margin+ntime]
        Eout[:,start:start+len(E[0])]=np.einsum('ft,tpd->fpd',phase,E,optimize=True)
        Bout[:,start:start+len(B[0])]=np.einsum('ft,tpd->fpd',phase,B,optimize=True)
    return Eout,Bout

def gaussian_point_reference(x,y,energies,p,sigma_xy=20e-6,sigma_z=20e-9,order=9,**window):
    nodes,w=roots_hermitenorm(order);w=w/np.sqrt(2*np.pi)
    xx,yy=np.meshgrid(nodes*sigma_xy,nodes*sigma_xy,indexing='xy');ww=np.outer(w,w).ravel()
    pts=np.c_[x-xx.ravel(),y-yy.ravel(),np.full(xx.size,p['screen'])]
    a,ab=arc_spectrum(pts,energies,p)
    v,vb=velocity_spectrum(pts,energies,p,sigma_time=sigma_z/(p['beta']*c),**window)
    form=np.exp(-.5*(np.asarray(energies)*e/hbar*sigma_z/(p['beta']*c))**2)
    return (np.einsum('fpd,p->fd',a,ww)*form[:,None],np.einsum('fpd,p->fd',ab,ww)*form[:,None],
            np.einsum('fpd,p->fd',v,ww),np.einsum('fpd,p->fd',vb,ww))

def energy_density_spectrum(E,B):
    """One-sided dW/(dE_eV dA), SI Fourier fields integrated in seconds."""
    return 2*e/(h*mu_0)*np.real(E[...,0]*np.conj(B[...,1])-E[...,1]*np.conj(B[...,0]))

if __name__=='__main__':
    import json
    from pathlib import Path
    p=constants();energies=np.array([4.,6.,8.,10.,12.]);rows=[]
    for x,y in [(-100e-6,0.),(-50e-6,0.),(0.,0.),(-50e-6,50e-6)]:
        a,ab,v,vb=gaussian_point_reference(x,y,energies,p)
        row={'x_m':x,'y_m':y,'energy_eV':energies.tolist(),
         'velocity_over_acceleration_E_norm':(np.linalg.norm(v,axis=-1)/np.linalg.norm(a,axis=-1)).tolist(),
         'radiation_J_per_eV_per_m2_per_electron_squared':energy_density_spectrum(a,ab).tolist(),
         'total_over_radiation_flux':(energy_density_spectrum(a+v,ab+vb)/energy_density_spectrum(a,ab)).tolist()}
        rows.append(row)
    out={'parameters':p,'gaussian_sigma_xy_m':20e-6,'gaussian_sigma_z_m':20e-9,'pilot':rows,
      'scope':'Independent theory planning only; no PIC acceptance'}
    target=Path(__file__).resolve().parents[1]/'THEORY_PILOT.json';target.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
