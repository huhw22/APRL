#include "yee_sbp_geometry.h"
#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>
#include <limits>
#include <stdexcept>
using namespace aprl;
struct CachedField {
 std::size_t nx,ny,nz;std::vector<long double> values;
 CachedField(std::size_t x,std::size_t y,std::size_t z):nx(x+8),ny(y+8),nz(z+8),values(3*nx*ny*nz,std::numeric_limits<long double>::quiet_NaN()){}
 std::size_t index(unsigned a,std::ptrdiff_t i,std::ptrdiff_t j,std::ptrdiff_t k)const{if(i < -4||j < -4||k < -4||i>=(std::ptrdiff_t)nx-4||j>=(std::ptrdiff_t)ny-4||k>=(std::ptrdiff_t)nz-4)throw std::runtime_error("Test collar");return ((a*nz+k+4)*ny+j+4)*nx+i+4;}
 long double operator()(unsigned a,std::ptrdiff_t i,std::ptrdiff_t j,std::ptrdiff_t k)const{const auto v=values[index(a,i,j,k)];if(!std::isfinite(v))throw std::runtime_error("Uninitialized test field");return v;}
 template<class F>void fill(unsigned margin,F f){for(std::ptrdiff_t k=-(int)margin;k<(std::ptrdiff_t)nz-8+margin;++k)for(std::ptrdiff_t j=-(int)margin;j<(std::ptrdiff_t)ny-8+margin;++j)for(std::ptrdiff_t i=-(int)margin;i<(std::ptrdiff_t)nx-8+margin;++i)for(unsigned a=0;a<3;++a)values[index(a,i,j,k)]=f(a,i,j,k);}
};
int main(){const double dx=5e-6,dy=5e-6,dz=1.2543262415359034e-6;EBGridGeometry g(12,11,9,dx,dy,dz,.8/(SI::c*std::sqrt(1/(dx*dx)+1/(dy*dy)+1/(dz*dz))),EBMaxwellSolver::Yee);YeeSBPGeometry m(g,2,2,0,g.nz);long double maxSBP=0.,maxEnergy=0.,maxMomentum=0.,maxTime=0.;
 for(unsigned seed=0;seed<64;++seed){
  auto analyticE=[&](unsigned a,std::ptrdiff_t i,std::ptrdiff_t j,std::ptrdiff_t k){return std::sin(.13L*(seed+1)*(a+1)+.31L*i+.22L*j-.19L*k)+.37L*std::cos(.17L*i-.47L*j+.23L*k+.11L*seed*a);};
  auto analyticB=[&](unsigned a,std::ptrdiff_t i,std::ptrdiff_t j,std::ptrdiff_t k){return (std::cos(.19L*(seed+2)*(a+1)-.28L*i+.33L*j+.21L*k)+.29L*std::sin(.13L*i+.37L*j-.41L*k+.17L*seed*a))/SI::c;};
  CachedField E(g.nx,g.ny,g.nz),B(g.nx,g.ny,g.nz),BP(g.nx,g.ny,g.nz),EN(g.nx,g.ny,g.nz),BN(g.nx,g.ny,g.nz);E.fill(4,analyticE);B.fill(4,analyticB);
  auto e=[&](unsigned a,std::ptrdiff_t i,std::ptrdiff_t j,std::ptrdiff_t k){return E(a,i,j,k);};auto b=[&](unsigned a,std::ptrdiff_t i,std::ptrdiff_t j,std::ptrdiff_t k){return B(a,i,j,k);};
  BP.fill(3,[&](unsigned a,std::ptrdiff_t i,std::ptrdiff_t j,std::ptrdiff_t k){return b(a,i,j,k)+m.faraday(e,a,i,j,k);});auto bp=[&](unsigned a,std::ptrdiff_t i,std::ptrdiff_t j,std::ptrdiff_t k){return BP(a,i,j,k);};
  EN.fill(2,[&](unsigned a,std::ptrdiff_t i,std::ptrdiff_t j,std::ptrdiff_t k){return e(a,i,j,k)+m.ampere(bp,a,i,j,k);});auto en=[&](unsigned a,std::ptrdiff_t i,std::ptrdiff_t j,std::ptrdiff_t k){return EN(a,i,j,k);};
  BN.fill(1,[&](unsigned a,std::ptrdiff_t i,std::ptrdiff_t j,std::ptrdiff_t k){return bp(a,i,j,k)+m.faraday(en,a,i,j,k);});auto bn=[&](unsigned a,std::ptrdiff_t i,std::ptrdiff_t j,std::ptrdiff_t k){return BN(a,i,j,k);};
  auto em=[&](unsigned a,std::ptrdiff_t i,std::ptrdiff_t j,std::ptrdiff_t k){return .5L*(e(a,i,j,k)+en(a,i,j,k));};
  auto de=[&](unsigned a,std::ptrdiff_t i,std::ptrdiff_t j,std::ptrdiff_t k){return en(a,i,j,k)-e(a,i,j,k);};
  auto bm=[&](unsigned a,std::ptrdiff_t i,std::ptrdiff_t j,std::ptrdiff_t k){return .5L*(b(a,i,j,k)+bp(a,i,j,k));};
  auto bmn=[&](unsigned a,std::ptrdiff_t i,std::ptrdiff_t j,std::ptrdiff_t k){return .5L*(bp(a,i,j,k)+bn(a,i,j,k));};
  const auto work=m.curlEnergyWork(e,b),surface=YeeSBPGeometry::sum(m.energyPort(e,b));const auto h0=.5L*SI::epsilon0*m.electricInner(e,e)+.5L/SI::mu0*m.magneticInner(b,bp),h1=.5L*SI::epsilon0*m.electricInner(en,en)+.5L/SI::mu0*m.magneticInner(bp,bn);const auto qh=YeeSBPGeometry::sum(m.energyPort(em,bp));
  const auto pe=YeeSBPGeometry::sum(m.electricStress(e)),pn=YeeSBPGeometry::sum(m.electricStress(en)),pb=YeeSBPGeometry::sum(m.magneticStress(bp));const auto force=.5L*(m.gaussElectricImpulse(e)+m.gaussElectricImpulse(en))+m.gaussMagneticImpulse(bp),p0=m.momentum(e,bm),p1=m.momentum(en,bmn);
  maxSBP=std::max(maxSBP,std::abs(work+surface)/(std::abs(work)+std::abs(surface)+1e-100L));maxEnergy=std::max(maxEnergy,std::abs(h1-h0+qh)/(std::abs(h0)+std::abs(h1)+std::abs(qh)+1e-100L));maxMomentum=std::max(maxMomentum,std::abs(p1-p0+.5L*(pe+pn)+pb+force)/(std::abs(p0)+std::abs(p1)+std::abs(pe)+std::abs(pn)+std::abs(pb)+std::abs(force)+1e-100L));const auto pm=YeeSBPGeometry::sum(m.electricStress(em)),pd=YeeSBPGeometry::sum(m.electricStress(de));maxTime=std::max(maxTime,std::abs(pm-.5L*(pe+pn)+.25L*pd)/(std::abs(pm)+std::abs(pe)+std::abs(pn)+std::abs(pd)+1e-100L));
 }
 std::cout<<std::setprecision(21)<<"{\"cases\":64,\"energy_surface_SBP_relative\":"<<maxSBP<<",\"finite_step_H_relative\":"<<maxEnergy<<",\"finite_step_P_relative\":"<<maxMomentum<<",\"stress_time_pairing_relative\":"<<maxTime<<"}\n";
 return maxSBP<1e-12L&&maxEnergy<1e-12L&&maxMomentum<1e-12L&&maxTime<1e-12L?0:1;
}
