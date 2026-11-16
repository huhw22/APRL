"""Exact integration of bilinear measured currents on clipped spacetime cells.

No conservation residual is used to infer any surface term. Laboratory
coordinates are obtained with stable light-front Lorentz expressions.
"""
import json
from pathlib import Path
import numpy as np
LD=np.longdouble

def read_series(directory):
    directory=Path(directory)
    assert (directory/'COMPLETE').read_text().split()==['1','257']
    m=json.loads((directory/'METADATA.json').read_text(),parse_float=LD)
    st=np.genfromtxt(directory/'states.csv',delimiter=',',names=True,dtype=LD)
    q=np.fromfile(directory/'field_lines.bin',dtype='<f8').reshape(len(st),int(m['nz']),6).astype(LD)
    assert np.isfinite(q).all() and np.all(np.diff(st['time_box_s'])>0)
    return m,st,q

def clip(poly,a,b,c):
    """Keep a+b*x+c*y >= 0; polygon coordinates within a unit cell."""
    if not poly:return []
    out=[];p=poly[-1];vp=a+b*p[0]+c*p[1]
    for n in poly:
        vn=a+b*n[0]+c*n[1]
        if (vp>=0)!=(vn>=0):
            t=vp/(vp-vn);out.append((p[0]+t*(n[0]-p[0]),p[1]+t*(n[1]-p[1])))
        if vn>=0:out.append(n)
        p,vp=n,vn
    return out

def moments(poly):
    """Area, integral x, integral y, integral xy, via exact triangle moments."""
    if len(poly)<3:return (LD(0),)*4
    a=poly[0];A=Mx=My=Mxy=LD(0)
    for i in range(1,len(poly)-1):
        b,c=poly[i],poly[i+1]
        area=abs((b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0]))/2
        sx=a[0]+b[0]+c[0];sy=a[1]+b[1]+c[1]
        A+=area;Mx+=area*sx/3;My+=area*sy/3
        Mxy+=area*(a[0]*a[1]+b[0]*b[1]+c[0]*c[1]+sx*sy)/12
    return A,Mx,My,Mxy

def interval(bounds,a,b,low,high):
    """Clip x interval by low <= a+b*x <= high."""
    l,h=bounds
    if b==0:return (l,h) if low<=a<=high else (LD(1),LD(0))
    p,q=(low-a)/b,(high-a)/b
    return max(l,min(p,q)),min(h,max(p,q))

def polynomial_trace(q00,q10,q01,q11,u,v,l,h):
    """Integral bilinear q(x,u+v*x) dx."""
    qx=q10-q00;qy=q01-q00;qxy=q11-q10-q01+q00
    a=q00+qy*u;b=qx+qy*v+qxy*u;c=qxy*v
    return a*(h-l)+b*(h*h-l*l)/2+c*(h*h*h-l*l*l)/3

class Tube:
    def __init__(self,m,states,values):
        self.m=m;self.t=np.array(states['time_box_s'],dtype=LD);self.o=np.array(states['origin_z_box_m'],dtype=LD);self.q=values
        self.g=LD(m['gamma']);self.b=np.sqrt((self.g-1)*(self.g+1))/self.g
        self.eta=1/(self.g*(1+self.b));self.c=LD(m['c_m_per_s']);self.dz=LD(m['dz_box_m'])
        self.zorig=LD(m['frame_lab_z_origin_m']);self.boxorig=LD(m['frame_box_z_origin_m'])
        self.tau_ref=LD(m['frame_lab_t_origin_s'])-self.zorig/self.c
        self.nz=values.shape[1];self.bz=self.g*self.dz;self.bt=-self.eta*self.dz/self.c
        self.dt=np.diff(self.t);do=np.diff(self.o)
        self.az=self.g*(do+self.b*self.c*self.dt)
        self.at=self.eta*(self.dt-do/self.c)
        self.zbase=self.zorig+self.g*(self.o[:-1]+self.dz/2-self.boxorig+self.b*self.c*self.t[:-1])
        self.tbase=self.eta*(self.t[:-1]-(self.o[:-1]+self.dz/2-self.boxorig)/self.c)

    def gate(self,z0,z1):
        # Relativity of simultaneity: t_lab at fixed z varies as t_box/gamma.
        # Express 1-beta stably as eta/gamma.
        ts=[self.t[0]/self.g-(self.eta/self.g)*(z-self.zorig)/self.c for z in [z0,z1]]
        te=[self.t[-1]/self.g-(self.eta/self.g)*(z-self.zorig)/self.c for z in [z0,z1]]
        return max(ts),min(te)

    def plane(self,z,lo,hi):
        total=coverage=LD(0)
        for n in range(len(self.dt)):
            k0=(z-self.zbase[n])/self.bz;k1=(z-self.zbase[n]-self.az[n])/self.bz
            for k in range(max(0,int(np.floor(min(k0,k1)))),min(self.nz-2,int(np.floor(max(k0,k1))))+1):
                u=(z-self.zbase[n]-self.bz*k)/self.bz;v=-self.az[n]/self.bz
                bounds=interval((LD(0),LD(1)),u,v,LD(0),LD(1))
                a=self.tbase[n]+self.bt*k+self.bt*u;b=self.at[n]+self.bt*v
                l,h=interval(bounds,a,b,lo,hi)
                if h<=l:continue
                q=self.q[:,:,1]
                total+=polynomial_trace(q[n,k],q[n+1,k],q[n,k+1],q[n+1,k+1],u,v,l,h)*self.dt[n]/self.g
                coverage+=(h-l)*self.dt[n]/self.g
        return total,coverage

    def cap(self,tau,z0,z1):
        total=coverage=LD(0)
        if z1<=z0:return total,coverage
        for n in range(len(self.dt)):
            k0=(tau-self.tbase[n])/self.bt;k1=(tau-self.tbase[n]-self.at[n])/self.bt
            for k in range(max(0,int(np.floor(min(k0,k1)))),min(self.nz-2,int(np.floor(max(k0,k1))))+1):
                u=(tau-self.tbase[n]-self.bt*k)/self.bt;v=-self.at[n]/self.bt
                bounds=interval((LD(0),LD(1)),u,v,LD(0),LD(1))
                a=self.zbase[n]+self.bz*k+self.bz*u;b=self.az[n]+self.bz*v
                l,h=interval(bounds,a,b,z0,z1)
                if h<=l:continue
                q=self.q[:,:,0]
                total+=polynomial_trace(q[n,k],q[n+1,k],q[n,k+1],q[n+1,k+1],u,v,l,h)*abs(b)
                coverage+=(h-l)*abs(b)
        return total,coverage

    def sides(self,z0,z1,lo,hi):
        total=np.zeros(4,dtype=LD);coverage=LD(0)
        if z1<=z0:return total,coverage
        for n in range(len(self.dt)):
            # Spatial clipping limits the k loop; time-gate clipping is exact.
            kmin=(z0-self.zbase[n]-max(self.az[n],LD(0)))/self.bz
            kmax=(z1-self.zbase[n]-min(self.az[n],LD(0)))/self.bz
            for k in range(max(0,int(np.floor(kmin))),min(self.nz-2,int(np.floor(kmax)))+1):
                t0=self.tbase[n]+self.bt*k;z=self.zbase[n]+self.bz*k
                ts=[t0,t0+self.at[n],t0+self.bt,t0+self.at[n]+self.bt]
                if max(ts)<lo or min(ts)>hi:continue
                poly=[(LD(0),LD(0)),(LD(1),LD(0)),(LD(1),LD(1)),(LD(0),LD(1))]
                for a,b,c in [(t0-lo,self.at[n],self.bt),(hi-t0,-self.at[n],-self.bt),(z-z0,self.az[n],self.bz),(z1-z,-self.az[n],-self.bz)]:
                    poly=clip(poly,a,b,c)
                    if not poly:break
                A,Mx,My,Mxy=moments(poly)
                if A==0:continue
                q00=self.q[n,k,2:];q10=self.q[n+1,k,2:];q01=self.q[n,k+1,2:];q11=self.q[n+1,k+1,2:]
                weight=self.dt[n]*self.dz
                total+=(q00*A+(q10-q00)*Mx+(q01-q00)*My+(q11-q10-q01+q00)*Mxy)*weight
                coverage+=A*weight
        return total,coverage

    def balance(self,z0,z1,lo,hi):
        entry,ce=self.plane(z0,lo,hi);exit,cx=self.plane(z1,lo,hi)
        early,ca=self.cap(lo,z0,z1);late,cb=self.cap(hi,z0,z1)
        sides,cs=self.sides(z0,z1,lo,hi)
        length=z1-z0;duration=hi-lo
        if length>0:
            errors=[abs(ca-length)/length,abs(cb-length)/length,abs(cs-length*duration)/(length*duration)]
        else:errors=[LD(0)]
        errors += [abs(ce-duration)/duration,abs(cx-duration)/duration]
        assert max(errors)<LD('2e-10'),('Uncovered boundary surface',errors)
        return dict(field_in=-entry,field_out=exit,cap_early=-early,cap_late=late,
            xminus=sides[0],xplus=sides[1],yminus=sides[2],yplus=sides[3],
            field_net=exit-entry,cap_change=late-early,side_sum=sum(sides),
            field_gain=exit-entry+late-early+sum(sides),max_coverage_error=max(errors))

def synthetic_test():
    # A nontrivial manufactured conserved laboratory current. Each of the
    # end-port difference, side flux and cap change is nonzero. It is a
    # quadrature/sign test, not a PIC or Maxwell solution.
    gamma=LD('929.1305492858543');beta=np.sqrt((gamma-1)*(gamma+1))/gamma;eta=1/(gamma*(1+beta));c=LD(299792458)
    dt=LD('4.1839819784122236e-15');dz=c*dt;nt,nz=97,129
    meta=dict(gamma=gamma,c_m_per_s=c,dz_box_m=dz,frame_lab_t_origin_s=LD(0),frame_lab_z_origin_m=LD(0),frame_box_z_origin_m=LD(0))
    states=np.zeros(nt,dtype=[('time_box_s',LD),('origin_z_box_m',LD)])
    states['time_box_s']=np.arange(nt,dtype=LD)*dt;states['origin_z_box_m']=-64*dz-np.arange(nt,dtype=LD)*dz
    q=np.zeros((nt,nz,6),dtype=LD)
    for n in range(nt):
        zbox=states['origin_z_box_m'][n]+(np.arange(nz,dtype=LD)+LD('.5'))*dz
        tau=eta*(states['time_box_s'][n]-zbox/c);z=gamma*(zbox+beta*c*states['time_box_s'][n]);T=LD('1e-15')
        q[n,:,0]=2*tau/T # u-Sz/c, integrated across unit area
        q[n,:,1]=c*LD('.3')*tau/T-z/T
        q[n,:,2]=LD('.2')/T;q[n,:,3]=LD('.2')/T
        q[n,:,4]=-LD('.7')/T;q[n,:,5]=-LD('.7')/T
    tube=Tube(meta,states,q);z0,z1=LD('-.02'),LD('.02');lo,hi=tube.gate(z0,z1)
    v=tube.balance(z0,z1,lo,hi);expected=(z1-z0)*(hi-lo)/LD('1e-15')
    scale=max(abs(v['field_net']),abs(v['cap_change']),abs(v['side_sum']))
    residual=v['field_gain'];relative=abs(residual)/scale
    assert relative<LD('1e-10'),(relative,v)
    assert abs(v['field_net']+expected)/expected<LD('1e-10')
    assert abs(v['cap_change']-2*expected)/expected<LD('1e-10')
    assert abs(v['side_sum']+expected)/expected<LD('1e-10')
    wrong=v['field_net']+v['side_sum']-v['cap_change']
    assert abs(wrong)/scale>1,'Wrong cap sign not detected'
    # Explicit Lorentz-coordinate forward/inverse and the Jacobian.
    assert max(abs((tube.at*tube.bz-tube.bt*tube.az)/(tube.dt*tube.dz)-1))<LD('1e-15')
    return dict(pass_=True,manufactured_complete_sum_relative_error=float(relative),
        end_difference=float(v['field_net']),side_flux=float(v['side_sum']),cap_difference=float(v['cap_change']),
        wrong_cap_sign_negative_control_detected=True,max_surface_coverage_relative_error=float(v['max_coverage_error']),
        scope='Independent exact spacetime geometry and outward boundary signs, with all three nonzero channels; not a Maxwell/PIC physics test.')

if __name__=='__main__':print(json.dumps(synthetic_test(),indent=2))
