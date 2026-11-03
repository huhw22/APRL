#include "eb_sources.h"
#include "eb_particles.h"
#include <cmath>
#include <cstdint>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <stdexcept>
using namespace aprl;
typedef FieldVector<Double> V;
static V vec(double x,double y,double z){V a(0.);a[0]=x;a[1]=y;a[2]=z;return a;}
static void need(bool x,const char* msg){if(!x)throw std::runtime_error(msg);}
static void close(double x,double y,double tol,const char* msg){
  need(std::isfinite(x)&&std::abs(x-y)<=tol,msg);
}
int main(){
  std::uint64_t hash=UINT64_C(1469598103934665603);
  for(int m=0;m<4;m++){
    SIMagneticElement a=SIMagneticElement::planarUndulatorFromK(1.092714657005,.0234,.05,10,.37);
    if(m<2){a.gaussianFringe=(m==0);a.prepare();}
    else {a.type=SIMagnetType::UniformDipole;a.length=.04;a.peakMagneticField=.3;a.prepare();}
    for(int j=0;j<600;j++){
      V p=vec(.0001*std::sin(j),.00012*std::cos(j),.015+j*.0005);
      SIFieldValue f;a.fieldsLab(p,0.,f);
      for(int k=0;k<3;k++){
        std::uint64_t bits;std::memcpy(&bits,&f.magnetic[k],sizeof bits);
        hash=(hash^bits)*UINT64_C(1099511628211);
      }
    }
  }
  std::cout<<"legacy_field_hash "<<hash<<"\n";
#ifdef SPARC_NEW
  SIMagneticElement q;q.type=SIMagnetType::Quadrupole;q.center[2]=.1;
  q.length=.053;q.quadrupoleGradient=-9.;q.gaussianFringe=false;q.prepare();
  SIFieldValue f;q.fieldsLab(vec(2e-4,-3e-4,.12),0,f);
  close(f.magnetic[0],.0027,1e-17,"Quadrupole Bx");
  close(f.magnetic[1],-.0018,1e-17,"Quadrupole By");
  need(f.magnetic[1]<0.,"Electron x force must focus");
  q.fieldsLab(vec(2e-4,-3e-4,.099),0,f);need(f.magnetic.norm2()==0.,"Quad upstream support");
  q.fieldsLab(vec(2e-4,-3e-4,.154),0,f);need(f.magnetic.norm2()==0.,"Quad downstream support");
  // Actual unchanged Boris pusher against paraxial thick-lens transfer.
  for(int axis=0;axis<2;axis++){
    q.polarizationAngle=0.;q.prepare();
    RelativisticParticleSI p;
    double gamma=298.26088438022987,uz=std::sqrt(gamma*gamma-1),beta=uz/gamma;
    p.position=vec(0.,0.,.1);p.position[axis]=1e-9;
    p.properVelocity=vec(0.,0.,uz);p.charge=-SI::elementaryCharge;p.mass=SI::electronMass;
    double dt=.053/(beta*SI::c*50000);
    for(int i=0;i<50000;i++){
      q.fieldsLab(p.position,0.,f);
      RelativisticBorisPusher::push(p,f.electric,f.magnetic,dt);
    }
    double k=9*SI::elementaryCharge/(SI::electronMass*SI::c*uz),a=std::sqrt(k),v=a*.053;
    double expected=axis==0?std::cos(v):std::cosh(v);
    double slope=axis==0?-a*std::sin(v):a*std::sinh(v);
    close(p.position[axis]/1e-9,expected,2e-6,"Boris quadrupole matrix position");
    close((p.properVelocity[axis]/p.properVelocity[2])/1e-9,slope,2e-6,"Boris quadrupole matrix angle");
    close(std::sqrt(1+p.properVelocity.norm2())/gamma,1.,5e-11,"Magnet changed lab gamma");
  }
  q.polarizationAngle=std::acos(-1.)/2.;q.prepare();
  q.fieldsLab(vec(2e-4,-3e-4,.12),0,f);
  close(f.magnetic[0],-.0027,1e-17,"Rotated quadrupole Bx");
  close(f.magnetic[1],.0018,1e-17,"Rotated quadrupole By");
  double divmax=0.,integral=0.;const double lu=.028,L=77*lu,dy=1e-7;
  for(double phase: {0.,.71,2.3}){
    SIMagneticElement a=SIMagneticElement::planarUndulatorFromK(2.060074987322757,lu,0.,77,std::acos(-1.)/2.);
    a.gaussianFringe=false;a.endRampPeriods=1.;a.undulatorPhase=phase;a.prepare();
    for(int i=1;i<800;i++){
      double z=L*(i+.113)/800.;
      SIFieldValue p,m,r,s;
      a.fieldsLab(vec(0.,30e-6+dy,z),0,p);a.fieldsLab(vec(0.,30e-6-dy,z),0,m);
      a.fieldsLab(vec(0.,30e-6,z+dy),0,r);a.fieldsLab(vec(0.,30e-6,z-dy),0,s);
      double div=(p.magnetic[1]-m.magnetic[1]+r.magnetic[2]-s.magnetic[2])/(2*dy);
      divmax=std::max(divmax,std::abs(div));
    }
    integral=0.;
    for(int i=0;i<77*512;i++){a.fieldsLab(vec(0.,0.,(i+.5)*lu/512.),0,f);integral+=f.magnetic[1]*lu/512.;}
    need(std::abs(integral)<1e-12,"Ramped magnet nonzero transverse field integral");
    a.fieldsLab(vec(0.,0.,0.),0,f);need(f.magnetic.norm2()==0.,"Nonzero entrance end field");
    a.fieldsLab(vec(0.,0.,L),0,f);need(f.magnetic.norm2()==0.,"Nonzero exit end field");
  }
  need(divmax<1e-4,"Internal end fields violate divergence");
  std::cout<<std::setprecision(16)<<"new_magnet_checks PASS divmax_T_per_m "<<divmax<<"\n";
#endif
}
