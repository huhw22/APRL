#include "particle_work_geometry.h"
#include "eb_deposition.h"
#include <iostream>
#include <random>
#include <stdexcept>
using namespace aprl;
using Vec=FieldVector<Double>;
long double dot(const YeeComponent& a,const YeeComponent& b) {
  long double s=0;for(std::size_t k=0;k<a.nz();++k)for(std::size_t j=0;j<a.ny();++j)for(std::size_t i=0;i<a.nx();++i)s+=(long double)a(i,j,k)*b(i,j,k);return s;
}
long double mesh(const EBFieldGrid& f) {
  const auto& g=f.geometry();return (long double)g.dt*g.dx*g.dy*g.dz*(dot(f.jx(),f.ex())+dot(f.jy(),f.ey())+dot(f.jz(),f.ez()));
}
int main() {
 try {
  const Double dz=1.2543262415359034e-6;
  EBGridGeometry g(24,24,24,5e-6,5e-6,dz,dz/SI::c,EBMaxwellSolver::CowanZ);
  EBFieldGrid f(g);Vec origin(0.);ChargeConservingCurrentDepositor deposit(f,origin);
  std::mt19937 random(20261031);std::uniform_real_distribution<Double> loc(7.1,15.9),shift(-2.8,2.8);
  unsigned checks=0;long double maximum=0;
  for(bool uniform:{true,false}) {
   YeeComponent* es[3]={&f.ex(),&f.ey(),&f.ez()};
   for(unsigned a=0;a<3;++a)for(std::size_t k=0;k<es[a]->nz();++k)for(std::size_t j=0;j<es[a]->ny();++j)for(std::size_t i=0;i<es[a]->nx();++i)
    (*es[a])(i,j,k)=1e5*(uniform ? (a+1.) : std::sin(.37*i+.21*j+.43*k+a)+.02*i*j-.003*i*k);
   for(unsigned n=0;n<96;++n) {
    Vec x(0.),y(0.);const Double spacing[3]={g.dx,g.dy,g.dz};
    for(unsigned a=0;a<3;++a){x[a]=loc(random)*spacing[a];y[a]=x[a]+shift(random)*spacing[a];}
    const Double q=(n%2 ? -1. : 1.)*SI::elementaryCharge*(1+n%7);
    f.clearCurrent();deposit.depositSegment(x,y,q);
    long double independent=trajectoryElectricWork(f,origin,x,y,q),actual=mesh(f);
    const long double scale=std::max(std::abs(actual),1e-20L),error=std::abs(independent-actual)/scale;
    maximum=std::max(maximum,error);if(error>2e-12L)throw std::runtime_error("Independent quadrature disagrees with actual deposition");
    if(uniform) {
      long double exact=0;for(unsigned a=0;a<3;++a)exact+=(long double)q*1e5*(a+1)*((long double)y[a]-x[a]);
      if(std::abs(exact-actual)>2e-12L*std::max(std::abs(exact),1e-20L))throw std::runtime_error("Uniform-field analytic work differs");
    }
    ++checks;
   }
  }
  Vec x(0.),m(0.),y(0.);const Double spacing[3]={g.dx,g.dy,g.dz};
  for(unsigned a=0;a<3;++a){x[a]=(8.1+.3*a)*spacing[a];m[a]=(11.4-.2*a)*spacing[a];y[a]=(9.7+.5*a)*spacing[a];}
  f.clearCurrent();deposit.depositSegment(x,m,-SI::elementaryCharge);deposit.depositSegment(m,y,-SI::elementaryCharge);
  long double curved=trajectoryElectricWork(f,origin,x,m,-SI::elementaryCharge)+trajectoryElectricWork(f,origin,m,y,-SI::elementaryCharge);
  if(std::abs(curved-mesh(f))>2e-12L*std::max(std::abs(curved),1e-20L))throw std::runtime_error("Subpath sum differs");
  auto chord=trajectoryElectricWork(f,origin,x,y,-SI::elementaryCharge);
  if(std::abs(curved-chord)<1e-23L)throw std::runtime_error("Curved-path negative control is degenerate");
  std::cout<<std::setprecision(21)<<"{\"passed\":true,\"random_3D_checks\":"<<checks<<",\"max_scaled_work_difference\":"<<maximum<<",\"curved_minus_chord_J\":"<<curved-chord<<",\"subpath_sum_verified\":true,\"new_PIC_simulations\":0}\n";
  return 0;
 }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
