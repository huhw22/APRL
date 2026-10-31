#include "source_coupling.h"
#include "particle_work_geometry.h"
#include <iomanip>
#include <iostream>
#include <random>
using namespace source_trial;
long double dotMesh(const EBFieldGrid& f,const Electric& old) {
 const YeeComponent* e[3]={&f.ex(),&f.ey(),&f.ez()};const YeeComponent* j[3]={&f.jx(),&f.jy(),&f.jz()};long double w=0;
 for(unsigned a=0;a<3;++a)for(std::size_t i=0;i<e[a]->size();++i)w+=(long double)j[a]->data()[i]*((long double)old[a].data()[i]+e[a]->data()[i])/2;
 const auto& g=f.geometry();return w*(long double)g.dt*g.dx*g.dy*g.dz;
}
int main() {
 try {
  std::mt19937 random(20261031);std::uniform_real_distribution<Double> unit(-1.,1.);const Double dt=1.2543262415359034e-6/SI::c;
  long double worstMomentum=0,worstEnergy=0,worstEnergyRoundoff=0,worstLab=0;unsigned maxRoots=0;
  BoostFrameTransform frame;frame.setFromGamma(929.1305492858543,SI::c,0.);
  for(unsigned n=0;n<512;++n) {
   Vec u(0.),e(0.),b(0.);const Double scale=n%3==0?.001:(n%3==1?.25:1.2);
   for(unsigned a=0;a<3;++a){u[a]=scale*unit(random);e[a]=1e9*unit(random);b[a]=1e3*unit(random);}
   const Vec before(u);maxRoots=std::max(maxRoots,momentum(u,e,b,-SI::elementaryCharge,SI::electronMass,dt/5));
   Vec v(before);v+=u;v*=SI::c/(BoostFrameTransform::gammaFromProperVelocity(before)+BoostFrameTransform::gammaFromProperVelocity(u));
   Vec expected(cross(v,b));expected+=e;expected*=-SI::elementaryCharge*dt/(5*SI::electronMass*SI::c);
   for(unsigned a=0;a<3;++a)worstMomentum=std::max(worstMomentum,(long double)std::abs(u[a]-before[a]-expected[a])/(1+std::abs(expected[a])));
   long double work=0;for(unsigned a=0;a<3;++a)work+=(long double)(-SI::elementaryCharge)*(dt/5)*e[a]*v[a];
   const auto delta=(long double)SI::electronMass*SI::c*SI::c*deltaGamma(before,u);
   worstEnergy=std::max(worstEnergy,std::abs(delta-work)/std::max(std::abs(work),1e-18L));
   const long double roundingScale=std::numeric_limits<Double>::epsilon()*(long double)SI::electronMass*SI::c*SI::c*(1+before.norm2()+u.norm2());
   worstEnergyRoundoff=std::max(worstEnergyRoundoff,std::abs(delta-work)/roundingScale);
   Vec labE(0.),labB(0.);for(unsigned a=0;a<3;++a)labB[a]=unit(random);
   frame.fieldsLabToBox(labE,labB,e,b);u=before;momentum(u,e,b,-SI::elementaryCharge,SI::electronMass,dt/5);
   const long double labDefect=deltaGamma(before,u)+(long double)frame.beta()*((long double)u[2]-before[2]);
   worstLab=std::max(worstLab,std::abs(labDefect)/gammaLD(before));
  }
  // This is a local algebra implementation check. Small electric work can
  // cancel despite finite particle momentum, so use its actual floating
  // arithmetic scale here; retain the work-relative number and independently
  // test the small physical exchange in the dynamical PIC cases.
  if(worstMomentum>5e-14L||worstEnergyRoundoff>128||worstLab>5e-14L){std::cerr<<std::setprecision(21)<<"local checks: momentum="<<worstMomentum<<" energy/work="<<worstEnergy<<" energy/roundoff="<<worstEnergyRoundoff<<" pure-lab-magnet="<<worstLab<<'\n';throw std::runtime_error("Mean-velocity Lorentz/work invariant failed");}
  EBGridGeometry g(24,24,24,5e-6,5e-6,1.2543262415359034e-6,dt,EBMaxwellSolver::CowanZ);EBFieldGrid f(g),previous(g);Vec origin(0.);
  YeeComponent* e[3]={&f.ex(),&f.ey(),&f.ez()};YeeComponent* old[3]={&previous.ex(),&previous.ey(),&previous.ez()};
  for(unsigned a=0;a<3;++a)for(std::size_t k=0;k<e[a]->nz();++k)for(std::size_t j=0;j<e[a]->ny();++j)for(std::size_t i=0;i<e[a]->nx();++i) {
   (*e[a])(i,j,k)=1e5*(std::sin(.3*i+.5*j+.2*k+a)+.02*i*j);
   (*old[a])(i,j,k)=1e5*(std::cos(.2*i-.1*j+.4*k+a)-.01*j*k);
  }
  const auto saved=copyElectric(previous);ChargeConservingCurrentDepositor deposit(f,origin);long double worstOrbit=0,worstGauss=0;
  YeeComponent* bs[3]={&f.bx(),&f.by(),&f.bz()};
  for(unsigned a=0;a<3;++a)for(std::size_t k=0;k<bs[a]->nz();++k)for(std::size_t j=0;j<bs[a]->ny();++j)for(std::size_t i=0;i<bs[a]->nx();++i)(*bs[a])(i,j,k)=1e-3*std::sin(.3*i+.5*j-.2*k+a);
  Double worstMagnetic=0.;
  for(unsigned n=0;n<192;++n) {
   Vec x(0.),y(0.);const Double ds[3]={g.dx,g.dy,g.dz};
   for(unsigned a=0;a<3;++a){x[a]=(10.3+2*unit(random))*ds[a];y[a]=x[a]+(n%4==a?0.:2.9*unit(random))*ds[a];}
   if(n%3==0)for(unsigned a=0;a<3;++a)y[a]=x[a]+1e-5*(y[a]-x[a]);
   const Double q=-SI::elementaryCharge*(1+n%9);f.clearCurrent();deposit.depositSegment(x,y,q);
   const auto magnetic=magneticSample(f,origin,x);Vec dummy(0.),standard(0.),after(0.);f.sampleFieldsPosition(x,origin,dummy,standard);f.sampleFieldsPosition(y,origin,dummy,after);
   for(unsigned a=0;a<3;++a)worstMagnetic=std::max(worstMagnetic,std::abs(magnetic.value[a]-standard[a]));
   after-=magnetic.value;const auto bound=magneticDifferenceBound(f,origin,y,magnetic);
   if(std::sqrt(after.norm2())>bound+1e-18)throw std::runtime_error("Magnetic interpolation Lipschitz bound failed");
   const auto average=orbitElectric(f,&saved,origin,x,y);const auto measured=chargeWork(q,average,x,y);const auto actual=dotMesh(f,saved);
   const auto independent=.5L*(trajectoryElectricWork(f,origin,x,y,q)+trajectoryElectricWork(previous,origin,x,y,q));
   const auto scale=std::max(std::abs(actual),1e-20L);worstOrbit=std::max(worstOrbit,std::abs(measured-actual)/scale);worstGauss=std::max(worstGauss,std::abs(measured-independent)/scale);
  }
  if(worstOrbit>5e-12L||worstGauss>5e-12L||worstMagnetic>1e-18)throw std::runtime_error("Orbit/magnetic interpolation differs from original interpolation/deposition");
  std::cout<<std::setprecision(21)<<"{\"checks_pass\":true,\"lorentz_samples\":512,\"orbit_samples\":192,\"max_scaled_momentum_equation_error\":"<<worstMomentum<<",\"max_scaled_electric_energy_error\":"<<worstEnergy<<",\"max_energy_error_in_arithmetic_roundoff_scales\":"<<worstEnergyRoundoff<<",\"max_scaled_pure_lab_magnetic_energy_error\":"<<worstLab<<",\"max_local_scalar_iterations\":"<<maxRoots<<",\"max_scaled_orbit_deposition_error\":"<<worstOrbit<<",\"max_scaled_orbit_independent_gauss_error\":"<<worstGauss<<"}\n";return 0;
 }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
