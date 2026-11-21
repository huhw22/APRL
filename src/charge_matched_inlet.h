#ifndef APRL_CHARGE_MATCHED_INLET_H
#define APRL_CHARGE_MATCHED_INLET_H
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>
#include "eb_field.h"
#include "eb_particles.h"
#include "longitudinal_shape.h"
#include "transverse_source_shape.h"
namespace aprl {
// A generator-independent completion of the NEW inlet layer only. The
// supplied background may be approximate; the incoming particle states and
// the retained fields determine its normal electric component. No particle,
// retained E/B, current, CPML history, or Maxwell/Boris stencil is modified.
class ChargeMatchedInlet {
public:
  struct Report {long double nativeWorkJ=0.;double maxBefore=0.,maxAfter=0.,maxDeltaEz=0.;std::size_t particles=0;};
  ChargeMatchedInlet(const EBGridGeometry& g,std::size_t gx,std::size_t gy,bool quintic=false,bool linearNotch=false,bool linearBinomial3=false,bool transverseXBinomial3=false)
    :g_(g),gx_(gx),gy_(gy),quintic_(quintic||linearNotch||linearBinomial3),linearNotch_(linearNotch),linearBinomial3_(linearBinomial3),xBinomial_(transverseXBinomial3),rho_((g.nx+1)*(g.ny+1),0.) {
    if(static_cast<unsigned>(quintic)+static_cast<unsigned>(linearNotch)+static_cast<unsigned>(linearBinomial3)>1)
      throw std::invalid_argument("Ambiguous charge-matched inlet source shape");
    if(xBinomial_ && !linearBinomial3_)
      throw std::invalid_argument("Transverse three-point inlet requires longitudinal linear-binomial3");
    if(g.nz<2 || gx+1>=g.nx-gx || gy+1>=g.ny-gy)
      throw std::invalid_argument("Charge-matched inlet needs an interior physical aperture");
  }
  void begin(const FieldVector<Double>& origin) {
    origin_=origin;std::fill(rho_.begin(),rho_.end(),0.);particles_=0;
  }
  void particle(const RelativisticParticleSI& p) {
    const double z=(p.position[2]-origin_[2])/g_.dz,wz=quintic_?longitudinal_shape::nodeWeight(z,1,linearBinomial3_?longitudinal_shape::Shape::LinearBinomial3:(linearNotch_?longitudinal_shape::Shape::LinearNotch:longitudinal_shape::Shape::QuinticCompensated)):std::max(0.,1.-std::abs(z-1.));
    if(wz==0.)return;
    const double x=(p.position[0]-origin_[0])/g_.dx,y=(p.position[1]-origin_[1])/g_.dy;
    if(!std::isfinite(x)||!std::isfinite(y)||!std::isfinite(z)||!std::isfinite(p.charge))
      throw std::runtime_error("Nonfinite charge-matched inlet source");
    if(y<0. || y>g_.ny)return;
    // The wider x shape can overlap a native node before its centre enters.
    if(xBinomial_ ? (x<=-2. || x>=g_.nx+2.) : (x<0. || x>g_.nx))return;
    const long long ix=static_cast<long long>(std::floor(x));
    const auto j=static_cast<std::size_t>(std::floor(y));
    const double fx=x-ix,fy=y-j,scale=p.charge*(1./(g_.dx*g_.dy*g_.dz));
    const auto wx=xBinomial_ ? transverse_source_shape::nodeWeights(fx) :
                               std::array<double,4>{{1.-fx,fx,0.,0.}};
    const unsigned nX=xBinomial_?4:2;
    const auto first=ix-(xBinomial_?1:0);
    for(unsigned b=0;b<2;++b)for(unsigned a=0;a<nX;++a){
      const auto i=first+a;
      if(i>=0 && i<=static_cast<long long>(g_.nx) && j+b<=g_.ny)
        rho_[(j+b)*(g_.nx+1)+static_cast<std::size_t>(i)]+=scale*wx[a]*(b?fy:1.-fy)*wz;
    }
    ++particles_;
  }
  Report complete(EBFieldGrid& f,bool audit) {
    Report report;report.particles=particles_;
    const long double volume=static_cast<long double>(g_.dx)*g_.dy*g_.dz;
    for(std::size_t j=gy_+1;j<g_.ny-gy_;++j)for(std::size_t i=gx_+1;i<g_.nx-gx_;++i) {
      const long double divxy=(static_cast<long double>(f.ex()(i,j,1))-f.ex()(i-1,j,1))/g_.dx+
        (static_cast<long double>(f.ey()(i,j,1))-f.ey()(i,j-1,1))/g_.dy;
      const long double rho=rho_[j*(g_.nx+1)+i],old=f.ez()(i,j,0),retained=f.ez()(i,j,1);
      const double value=static_cast<double>(retained+g_.dz*(divxy-rho/SI::epsilon0));
      if(!std::isfinite(value))throw std::runtime_error("Nonfinite charge-matched inlet field");
      if(audit) {
        report.maxBefore=std::max(report.maxBefore,static_cast<double>(std::abs(SI::epsilon0*(divxy+(retained-old)/g_.dz)-rho)));
        report.maxAfter=std::max(report.maxAfter,static_cast<double>(std::abs(SI::epsilon0*(divxy+(retained-value)/g_.dz)-rho)));
        report.maxDeltaEz=std::max(report.maxDeltaEz,std::abs(value-static_cast<double>(old)));
        report.nativeWorkJ+=.5L*SI::epsilon0*volume*(value-old)*(value+old);
      }
      f.ez()(i,j,0)=value;
    }
    return report;
  }
private:
  EBGridGeometry g_;std::size_t gx_,gy_,particles_=0;FieldVector<Double> origin_;
  bool quintic_,linearNotch_,linearBinomial3_;
  const bool xBinomial_;
  std::vector<double> rho_;
};
}
#endif
