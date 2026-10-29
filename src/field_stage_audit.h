#ifndef APRL_FIELD_STAGE_AUDIT_H
#define APRL_FIELD_STAGE_AUDIT_H

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include "eb_field.h"

namespace aprl {
// A diagnostic observer of completed operations. Fields, sources, convolution
// histories and the MPI schedule are never changed. Energies use exactly the
// existing physical-cell diagnostic, so operation deltas can be checked
// against its ledger. This is not a new Maxwell energy invariant.
class FieldStageAudit {
public:
  FieldStageAudit(const std::string& directory,int rank,std::size_t globalNz,
      std::size_t localOffset,std::size_t cx,std::size_t cy)
    :rank_(rank),globalNz_(globalNz),offset_(localOffset),cx_(cx),cy_(cy) {
    std::ostringstream name;
    name<<directory<<"/field-stages-rank-"<<std::setw(5)<<std::setfill('0')<<rank<<".csv";
    std::ifstream old(name.str());
    if(old.good())throw std::runtime_error("Field-stage audit output exists");
    file_.open(name.str());
    if(!file_)throw std::runtime_error("Cannot create field-stage audit output");
    file_<<std::setprecision(21)
      <<"step,stage,particles,carriers,center_energy_J,left2_energy_J,right2_energy_J,right8_energy_J,"
        "native_electric_energy_J,native_magnetic_energy_J,current_sampled,integral_native_J2_A2_per_m,max_abs_J_any_A_per_m2,"
        "upper_tangential_max_delta_V_per_m,power_xminus_W,power_xplus_W,power_yminus_W,power_yplus_W,power_zminus_W,power_zplus_W\n";
  }
  void begin(const EBFieldGrid& f,std::size_t step,std::size_t particles,std::size_t carriers) {
    step_=step+1;populations(particles,carriers);
    upper_.clear();
    if(offset_+f.geometry().nz==globalNz_)
      for(std::size_t j=cy_;j<f.geometry().ny-cy_;++j)
        for(std::size_t i=cx_;i<f.geometry().nx-cx_;++i) {
          upper_.push_back(f.ex()(i,j,f.geometry().nz));
          upper_.push_back(f.ey()(i,j,f.geometry().nz));
        }
    sample("begin",f);
  }
  void populations(std::size_t particles,std::size_t carriers) {
    particles_=particles;carriers_=carriers;
  }
  static void observe(void* context,const char* stage,const EBFieldGrid& f) {
    static_cast<FieldStageAudit*>(context)->sample(stage,f);
  }
  void sample(const char* stage,const EBFieldGrid& f) {
    const auto& g=f.geometry();
    const bool sampleCurrent=std::string(stage)=="electric_mpi";
    const YeeComponent* e[3]={&f.ex(),&f.ey(),&f.ez()};
    const YeeComponent* b[3]={&f.bx(),&f.by(),&f.bz()};
    const YeeComponent* current[3]={&f.jx(),&f.jy(),&f.jz()};
    const long double volume=static_cast<long double>(g.dx)*g.dy*g.dz;
    long double u=0.,left=0.,right=0.,right8=0.,ue=0.,ub=0.,j2=0.;
#ifdef APRL_HAS_OPENMP
#pragma omp parallel for reduction(+:u,left,right,right8,ue,ub,j2) schedule(static)
#endif
    for(std::ptrdiff_t kk=0;kk<static_cast<std::ptrdiff_t>(g.nz);++kk) {
      const auto k=static_cast<std::size_t>(kk),globalK=offset_+k;
      for(std::size_t j=cy_;j<g.ny-cy_;++j)
        for(std::size_t i=cx_;i<g.nx-cx_;++i) {
          const long double cell=static_cast<long double>(f.fieldEnergyDensityCell(i,j,k))*volume;
          u+=cell;if(globalK<2)left+=cell;
          if(globalK+2>=globalNz_)right+=cell;
          if(globalK+8>=globalNz_)right8+=cell;
          for(unsigned a=0;a<3;++a) {
            long double e2=0.,b2=0.,js=0.;
            for(unsigned p=0;p<2;++p)for(unsigned q=0;q<2;++q) {
              const std::size_t x=i+(a==0?0:p),y=j+(a==0?p:(a==1?0:q)),z=k+(a==2?0:q);
              const long double ev=(*e[a])(x,y,z);e2+=ev*ev;
              if(sampleCurrent) {
                const long double jv=(*current[a])(x,y,z);js+=jv*jv;
              }
            }
            for(unsigned p=0;p<2;++p) {
              const long double bv=(*b[a])(i+(a==0?p:0),j+(a==1?p:0),k+(a==2?p:0));
              b2+=bv*bv;
            }
            ue+=.125L*SI::epsilon0*e2*volume;
            ub+=.25L/SI::mu0*b2*volume;
            if(sampleCurrent)j2+=js*.25L*volume;
          }
        }
    }
    long double power[6]={};
    const auto i1=g.nx-cx_,j1=g.ny-cy_;
    const long double yz=static_cast<long double>(g.dy)*g.dz,xz=static_cast<long double>(g.dx)*g.dz,xy=static_cast<long double>(g.dx)*g.dy;
    for(std::size_t k=0;k<g.nz;++k) {
      for(std::size_t j=cy_;j<j1;++j) {
        power[0]-=static_cast<long double>(f.poyntingComponentCell(cx_,j,k,0))*yz;
        power[1]+=static_cast<long double>(f.poyntingComponentCell(i1-1,j,k,0))*yz;
      }
      for(std::size_t i=cx_;i<i1;++i) {
        power[2]-=static_cast<long double>(f.poyntingComponentCell(i,cy_,k,1))*xz;
        power[3]+=static_cast<long double>(f.poyntingComponentCell(i,j1-1,k,1))*xz;
      }
    }
    for(std::size_t j=cy_;j<j1;++j)for(std::size_t i=cx_;i<i1;++i) {
      if(offset_==0)power[4]-=static_cast<long double>(f.poyntingComponentCell(i,j,0,2))*xy;
      if(offset_+g.nz==globalNz_)power[5]+=static_cast<long double>(f.poyntingComponentCell(i,j,g.nz-1,2))*xy;
    }
    long double upperDelta=0.;std::size_t n=0;
    if(!upper_.empty())for(std::size_t j=cy_;j<j1;++j)for(std::size_t i=cx_;i<i1;++i) {
      upperDelta=std::max(upperDelta,std::abs(static_cast<long double>(f.ex()(i,j,g.nz))-upper_[n++]));
      upperDelta=std::max(upperDelta,std::abs(static_cast<long double>(f.ey()(i,j,g.nz))-upper_[n++]));
    }
    Double maxJ=0.;
    if(sampleCurrent)for(unsigned a=0;a<3;++a)
      for(std::size_t index=0;index<current[a]->size();++index)
        maxJ=std::max(maxJ,std::abs(current[a]->data()[index]));
    file_<<step_<<','<<stage<<','<<particles_<<','<<carriers_<<','<<u<<','<<left<<','<<right<<','<<right8
      <<','<<ue<<','<<ub<<','<<sampleCurrent<<','<<j2<<','<<maxJ<<','<<upperDelta;
    for(const auto value:power)file_<<','<<value;
    file_<<'\n';if(!file_)throw std::runtime_error("Cannot write field-stage audit");
  }
  void close() {file_.flush();if(!file_)throw std::runtime_error("Cannot flush field-stage audit");file_.close();}
private:
  int rank_;
  std::size_t globalNz_,offset_,cx_,cy_,step_=0,particles_=0,carriers_=0;
  std::vector<Double> upper_;
  std::ofstream file_;
};
}
#endif
