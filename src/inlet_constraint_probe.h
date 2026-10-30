#ifndef APRL_INLET_CONSTRAINT_PROBE_H
#define APRL_INLET_CONSTRAINT_PROBE_H
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <vector>
#include "eb_deposition.h"
#include "eb_particles.h"
namespace aprl {
// Read-only diagnostic: use the original CIC charge diagnostic, then observe
// eps*div(E)-rho and its evolution on the first retained inlet nodes. No
// correction, field/particle write, MPI call, or convolution-history reset.
class InletConstraintProbe {
public:
  InletConstraintProbe(const std::string& directory,const EBGridGeometry& g)
    :g_(g),rho_(g.nx+1,g.ny+1,g.nz+1),layers_(std::min<std::size_t>(8,g.nz-1)),
      i0_(std::max<std::size_t>(1,g.nx/2>24?g.nx/2-24:1)),i1_(std::min(g.nx-1,g.nx/2+24)),
      j0_(std::max<std::size_t>(1,g.ny/2>24?g.ny/2-24:1)),j1_(std::min(g.ny-1,g.ny/2+24)) {
    const auto file=directory+"/inlet-constraint.csv";std::ifstream exists(file);
    if(exists.good())throw std::runtime_error("Inlet probe file exists");
    out_.open(file);if(!out_)throw std::runtime_error("Cannot create inlet probe");
    out_<<std::setprecision(21)<<"step,stage,k,charge_C,gauss_signed_C,gauss_L1_C,gauss_rms_C_per_m3,continuity_max_C_per_m3,ampere_gauss_max_C_per_m3,axis_Ez_V_per_m\n";
  }
  void begin(const EBFieldGrid& f,const std::vector<RelativisticParticleSI>& particles,const FieldVector<Double>& origin,std::size_t step) {
    step_=step;charge(f,particles,origin);oldRho_.clear();oldDiv_.clear();
    for(std::size_t k=1;k<=layers_;++k)for(std::size_t j=j0_;j<=j1_;++j)for(std::size_t i=i0_;i<=i1_;++i) {
      oldRho_.push_back(rho_(i,j,k));oldDiv_.push_back(divergence(f.ex(),f.ey(),f.ez(),i,j,k));
    }
    record(f,"begin");
  }
  void finish(const EBFieldGrid& f,const std::vector<RelativisticParticleSI>& particles,const FieldVector<Double>& origin) {
    charge(f,particles,origin);record(f,"electric_complete");
  }
  void close(){out_.flush();if(!out_)throw std::runtime_error("Cannot flush inlet probe");out_.close();}
private:
  EBGridGeometry g_;YeeComponent rho_;std::size_t layers_,i0_,i1_,j0_,j1_,step_=0;
  std::vector<long double> oldRho_,oldDiv_;std::ofstream out_;
  long double divergence(const YeeComponent& x,const YeeComponent& y,const YeeComponent& z,std::size_t i,std::size_t j,std::size_t k) const {
    return (static_cast<long double>(x(i,j,k))-x(i-1,j,k))/g_.dx+
      (static_cast<long double>(y(i,j,k))-y(i,j-1,k))/g_.dy+
      (static_cast<long double>(z(i,j,k))-z(i,j,k-1))/g_.dz;
  }
  void charge(const EBFieldGrid& f,const std::vector<RelativisticParticleSI>& particles,const FieldVector<Double>& origin) {
    rho_.fill(0.);ChargeConservingCurrentDepositor depositor(const_cast<EBFieldGrid&>(f),origin);
    for(const auto& p:particles) {
      // Outside particles cannot touch retained charge nodes k>=1. Exclude
      // the distant part of this rank because only the inlet band is read.
      if(p.position[2]<origin[2]||p.position[2]>origin[2]+(layers_+1)*g_.dz)continue;
      depositor.depositCharge(p.position,p.charge,rho_);
    }
  }
  void record(const EBFieldGrid& f,const char* stage) {
    const bool end=std::string(stage)=="electric_complete";std::size_t n=0;
    const long double volume=static_cast<long double>(g_.dx)*g_.dy*g_.dz;
    for(std::size_t k=1;k<=layers_;++k) {
      long double q=0.,signedG=0.,l1=0.,square=0.,continuity=0.,ampere=0.;std::size_t count=0;
      for(std::size_t j=j0_;j<=j1_;++j)for(std::size_t i=i0_;i<=i1_;++i,++n) {
        const long double div=divergence(f.ex(),f.ey(),f.ez(),i,j,k),rho=rho_(i,j,k),defect=SI::epsilon0*div-rho;
        q+=rho*volume;signedG+=defect*volume;l1+=std::abs(defect)*volume;square+=defect*defect;++count;
        if(end) {
          const long double dj=divergence(f.jx(),f.jy(),f.jz(),i,j,k);
          continuity=std::max(continuity,std::abs(rho-oldRho_[n]+g_.dt*dj));
          ampere=std::max(ampere,std::abs(SI::epsilon0*(div-oldDiv_[n])+g_.dt*dj));
        }
      }
      out_<<step_<<','<<stage<<','<<k<<','<<q<<','<<signedG<<','<<l1<<','<<std::sqrt(square/count)<<','<<continuity<<','<<ampere<<','<<f.ez()(g_.nx/2,g_.ny/2,k)<<'\n';
    }
    if(!out_)throw std::runtime_error("Cannot write inlet probe");
  }
};
}
#endif
